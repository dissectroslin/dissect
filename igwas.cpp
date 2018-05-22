/*****************************************************************************
 *  DISSECT: a tool for performing genomic analysis with large sample sizes
 *  Copyright (C) 2014-2015 Oriol Canela-Xandri and Albert Tenesa
 *                          The Roslin Institute (University of Edinburgh)
 *
 *  This file is part of DISSECT.
 *
 *  DISSECT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DISSECT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DISSECT.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include "igwas.h"
#include "analysis.h"
#include "reml.h"
#include "auxiliar.h"
#include "labeledmatrix.h"

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream> 
#include <cmath>
#include <cstdlib>
#include <cstdio>

IGWAS::IGWAS()
{
  this->currentFile = "";
  this->currentGRMFile = "";
  this->currentGRMBase = NULL;
  
  this->useCovariateMatrix = false;
  this->currentCovariance = NULL;
  this->currentCovarianceMatrix = NULL;
  int dfReducedModel = -1;
  
  this->testingCovariate = NULL;
  
  this->accumulatedVariance.clear();
  this->nTests = 0;
  
  if(options.genotypeFile != "")
  {
    this->genotypeFiles.push_back(options.genotypeFile);
  }
  else if(options.genotypeListFile != "")
  {
    getListFromFile(options.genotypeListFile, this->genotypeFiles);
  }
  else if(options.genotypeAndGRMsListFile != "")
  {
    std::vector< std::vector<std::string> > files;
    getTableFromFile(options.genotypeAndGRMsListFile, files, 2);
    std::map<std::string, std::vector<std::string> > grms2genotypes;
    for(int i = 0; i<files.size(); i++) //We load genotype and grm files in this way, because we want all genotype with same grms together in the list.
    {
      misc.checkFileExists(files[i][0] + ".bed");
      misc.checkFileExists(files[i][0] + ".bim");
      misc.checkFileExists(files[i][0] + ".fam");
      misc.checkFileExists(files[i][1] + ".grm.ids");
      misc.checkFileExists(files[i][1] + ".grm.dat");
      misc.checkFileExists(files[i][1] + ".grm.snps");
      grms2genotypes[ files[i][1] ].push_back( files[i][0] );
    }
    for( std::map<std::string, std::vector<std::string> >::iterator it = grms2genotypes.begin(); it != grms2genotypes.end(); ++it )
    {
      std::vector<std::string> genoFiles = it->second;
      std::string grmFile = it->first;
      for(int i = 0; i<genoFiles.size(); i++)
      {
        if( this->grmFiles.count(genoFiles[i]) != 0 )
        {
          misc.error("Error: There is at least one genotype file repeated in [ " + options.genotypeAndGRMsListFile + " ] file. Genotype files have to be unique.", 0);
        }
        this->genotypeFiles.push_back(genoFiles[i]);
        this->grmFiles[ genoFiles[i] ] = grmFile;
      }
    }
  }
  else
  {
    misc.error("Error: An internal error was happened when loading the genotype for a GWAS analysis. No files specified.", 0);
  }
  
}

IGWAS::~IGWAS()
{
}

void IGWAS::computeGWAS()
{
  std::string backupOutFile = options.outFile;
  if(this->genotypeFiles.size() > 1)
  {
    misc.message.tab = "  ";
  }
  
  for(int i = 0; i<this->genotypeFiles.size(); i++)
  {
    //Load needed data: genotypes, covars, phenotypes and grms
    this->currentFile = this->genotypeFiles[i];
    Genotype * genotype = new Genotype(this->genotypeFiles[i]);
    //Phenotype * phenotype = new Phenotype(cyclicDistribution, options.phenotypesFile, options.phenotypeColumn);
    //Phenotype * phenotype = NULL;
    std::set<std::string> setPrefilterCovInds(genotype->individualIds.begin(), genotype->individualIds.end());
    Covariate * covariate = new Covariate(options.covarsFile, options.qCovarsFile, genotype->individualIds, true, setPrefilterCovInds);
  
    genotype->normalizeGenotypes();
    
    std::vector<std::string> commonIndividualsInGenotypeOrder;
    if(this->grmFiles.count( this->currentFile ) == 0 && options.discreteRandomEffectsFileColumns.size() == 0 && options.GCTAGRMsFile == "")
    {
      //Get shared individuals and SNPs in group and filter
      std::vector<std::string> commonIndividuals = intersectionStringVectors(2, &genotype->individualIds, &covariate->individualIds);
      commonIndividualsInGenotypeOrder = orderVectorAsTemplate(genotype->individualIds, commonIndividuals);
      
      setPrefilterCovInds = std::set<std::string>(commonIndividualsInGenotypeOrder.begin(), commonIndividualsInGenotypeOrder.end());
      if(options.igwasTestCovariates == "")
      {
        misc.error("Error: A file defining the covariates to test is expected.", 0);
      }
      this->testingCovariate = new Covariate(options.igwasTestCovariates, options.igwasTestQCovariates, std::vector<std::string>(), false, setPrefilterCovInds);
      if (this->testingCovariate->nIndividuals == 0)
      {
        misc.error("Error: There are not overlapping individuals between covariates to test, other covariates and the genotype files.", 0);
      }
      this->testingCovariate->parseRawCovariates(std::vector<std::string>(), 0, 0);
      commonIndividualsInGenotypeOrder = orderVectorAsTemplate(commonIndividualsInGenotypeOrder, this->testingCovariate->individualIds);
      
      this->useCovariateMatrix = false;
      
      misc.message << "Fitting SNPs in file " << this->currentFile << "..." << std::endl;
    }
    else
    {
      this->accumulatedVariance.clear();
      this->nTests = 0;
      
      bool success = computeCovariance(genotype, NULL, covariate);
      if(success == false)
      {
        misc.message << "Error: The covariance file specified for the genotypes in " << this->genotypeFiles[i] << " cannot be computed. Skipping this file.";
        continue;
      }
      
      commonIndividualsInGenotypeOrder = orderVectorAsTemplate(genotype->individualIds, currentCovariance->individualIds);
      if(commonIndividualsInGenotypeOrder != currentCovariance->individualIds)
      {
        misc.error("Error: An internal error was happened. Unexpected order of individuals when preparing the covariance matrix.", 0);
      }
      
      this->useCovariateMatrix = true;
      
      misc.message << "Fitting SNPs in file " << this->currentFile << " corrected by a covariance..." << std::endl;
    }
    
    //phenotype->filterIndividuals(commonIndividualsInGenotypeOrder);
    covariate->filterIndividuals(commonIndividualsInGenotypeOrder);
    genotype->filterSNPsAndIndividuals(genotype->SNPIds, commonIndividualsInGenotypeOrder, true);
    
    if( this->useCovariateMatrix == true )
    {
      if(this->testingCovariate != NULL)
      {
        misc.error("Error: An internal error has happened. Unexpected testing covariate matrices pointer on igwas analysis preparation.", 0);
      }
      //If GRM is diagonal, correct phenotypes, genotypes and covariates
      if(this->currentCovariance->diagonalized == true)
      {
        /*Matrix * tempPhenos = new Matrix();
        tempPhenos->multiply(this->currentCovariance->eigenVectors, 'T', phenotype->phenotypes, 'N');
        delete phenotype->phenotypes;
        phenotype->phenotypes = tempPhenos;*/
        
        Matrix * tempCovars = new Matrix();
        tempCovars->multiply(this->currentCovariance->eigenVectors, 'T', covariate->covariates, 'N');
        delete covariate->covariates;
        covariate->covariates = tempCovars;
        
        Matrix * tempGenos = new Matrix();
        tempGenos->multiply(genotype->genotypes, 'N', this->currentCovariance->eigenVectors, 'N'); //Geno*EigV = (EigVt*Genot)t
        delete genotype->genotypes;
        genotype->genotypes = tempGenos;
      }
    }
    
    if( this->genotypeFiles.size() > 1 )
    {
      options.outFile += "." + getFileName(this->genotypeFiles[i]);
    }
    
    //Start analysis
    if(options.analysis != recursiveGWASAnalysis) //Perform grouped GWAS or standard GWAS
    {
      misc.setGetElapsedTime("GWAS");
      misc.message << "Starting GWAS analysis..." << std::endl;
      if(options.regionBy != ungrouped)
      {
        /*genotype->groupSNPs(options.regionBy);
        computeGroupedGWAS(genotype, phenotype, covariate);*/
        misc.error("Not valid analysis.", 0);
      }
      else
      {
        if(options.parallelGWAS == false)
        {
          //computeIndividualGWAS(genotype, NULL, covariate);
          misc.error("Not valid analysis.", 0);
        }
        else
        {
          computeIndividualGWASGroupedCommunicator(genotype, NULL, covariate);
        }
      }
      misc.message << "GWAS analysis finished after " << misc.setGetElapsedTime("GWAS", true) << "." << std::endl;
    }
    else //Perform recursive GWAS
    {
      /*misc.setGetElapsedTime("recursiveGWAS");
      this->significantSNPIds.clear();
      std::vector<std::string> subsetToAnalyzeSNPIds = genotype->SNPIds;
      int iteration = 1;
      std::string backupOutFile2 = options.outFile;
      options.fixedGroupSize = floor(double(genotype->nIndividuals)*options.relationFitSNPsIndividuals);
      if(options.fixedGroupSize < 1)
      {
        misc.error("Error: The resultant SNPs group size for a recursive GWAS is less than 1. This could happen because sample size is too small or because the parameter passed by --rgwas-ratio is too small.", 0);
      }
      misc.message << "Starting recursive GWAS analysis using groups of " << options.fixedGroupSize << " SNPs..." << std::endl;
      misc.message.tab = "  ";
      while( misc.gt( subsetToAnalyzeSNPIds.size() != this->significantSNPIds.size() && (options.recursiveGWASMaxIterations < 1 || options.recursiveGWASMaxIterations >= iteration) ) )
      {
        if(iteration > options.significanceThresholdsFilterSNPs.size())
        {
          this->significanceThreshold = options.significanceThresholdsFilterSNPs[ options.significanceThresholdsFilterSNPs.size() - 1 ];
        }
        else
        {
          this->significanceThreshold = options.significanceThresholdsFilterSNPs[ iteration - 1 ];
        }
        
        options.outFile = backupOutFile2 + ".iter" + getString(iteration);
        subsetToAnalyzeSNPIds = this->significantSNPIds;
        misc.message << "Performing iteration " << iteration << "." << std::endl;
        
        genotype->groupSNPs(byOrderedFixedSize, subsetToAnalyzeSNPIds);
        computeGroupedGWAS(genotype, phenotype, covariate);
        
        misc.message << "Analysis finished. " << ((iteration==1)?genotype->nSNPs:subsetToAnalyzeSNPIds.size()) << " SNPs analyzed." << std::endl;
        
        iteration++;
        
        if( misc.gt(this->significantSNPIds.size() == 0) )
        {
          misc.message << "No significant SNPs found. Stopping iterations." << std::endl;
          break;
        }
      }
      options.outFile = backupOutFile2;
      misc.message.tab = "";
      misc.message << "recursive GWAS analysis finished after " << iteration - 1 << " iterations. It needed " << misc.setGetElapsedTime("recursiveGWAS", true) << "." << std::endl;
      */
      misc.error("Not valid analysis.", 0);
    }
    
    options.outFile = backupOutFile;
    
    delete genotype;
    //delete phenotype;
    delete covariate;
  }
  
  this->currentFile = "";
  
  if(this->currentGRMBase != NULL)
  {
    delete this->currentGRMBase;
    this->currentGRMBase = NULL;
  }
  if( this->currentCovariance != NULL )
  {
    delete this->currentCovariance;
    this->currentCovariance = NULL;
    this->currentCovarianceMatrix = NULL;
  }
  
  if( this->testingCovariate != NULL )
  {
    delete this->testingCovariate;
    this->testingCovariate = NULL;
  }
  
  if(this->genotypeFiles.size() > 1)
  {
    misc.message.tab = "";
    misc.message << "GWAS analysis finished on all files." << std::endl;
  }
}



void IGWAS::computeIndividualGWASGroupedCommunicator(Genotype * gcgenotype, Phenotype * gcphenotype, Covariate * gccovariate)
{

  std::map<std::string, GLMResults > results;
  std::map<int, GLMResults > groupedResults;
  std::map<std::string, std::vector<SNP> > resultSNPInfo;
  std::vector<std::string> unfittedSNPs;

  Communicator * globalCommunicator = communicator;  
  Communicator * groupedCommunicator = new Communicator(globalCommunicator, basicGroupedCommunicator);
 
  
  std::map< int, std::vector<int> > SNPidxs;
  
  Matrix * genotype = gcgenotype->genotypesRedistributionToGroupedCommunicatorMatrices(groupedCommunicator, SNPidxs);
  //Matrix * phenotype = gcphenotype->phenotypes->copyToGroupedCommunicator(groupedCommunicator);
  Matrix * covariate = gccovariate->covariates->copyToGroupedCommunicator(groupedCommunicator);
  
  Matrix * globalCovarianceMatrix = NULL;
  Matrix * distributedTestingCovariate = NULL;
  if( this->useCovariateMatrix == true )
  {
    globalCovarianceMatrix = this->currentCovarianceMatrix;
    this->currentCovarianceMatrix = this->currentCovarianceMatrix->copyToGroupedCommunicator(groupedCommunicator);
  }
  else
  {
    distributedTestingCovariate = this->testingCovariate->covariates->copyToGroupedCommunicator(groupedCommunicator);
  }
  
  //std::stringstream logFile;
  //std::stringstream logOut;
  //misc.message.output = &logOut;
  //misc.message.log = &logFile;
  
  communicator = groupedCommunicator;
  
  /*GLMResults reducedResults;
  Matrix * br = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  computeGLM(phenotype, covariate, br, reducedResults);
  delete br;
  */
  
  Matrix * gTransposed = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  gTransposed->transpose(genotype);
  
  if( SNPidxs[ communicator->group ].size() != genotype->nGlobRows )
  {
    misc.error("Error: An internal error was happened. When performing individual SNP test in GWAS, discordance in grouped communicator size with local genotype size.", 0);
  }
  
  int stepsForPercentageOutput = (genotype->nGlobRows/10);
  stepsForPercentageOutput = (stepsForPercentageOutput == 0?1:stepsForPercentageOutput);
  
  for(int isnp = 0; isnp < genotype->nGlobRows; isnp++ )
  {
    if( communicator->group == 0 && isnp % stepsForPercentageOutput == 1 )
    {
      misc.message << (100*isnp)/genotype->nGlobRows << "% completed." << std::endl;
    }
    
    //Get data for one SNP
    Matrix * snpData = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, genotype->nGlobCols, 1, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    snpData->fillWithConstant(0.);
    snpData->add(gTransposed, 0., 1., subMatrix(0, 0, snpData->nGlobRows, 1), subMatrix(0, isnp, gTransposed->nGlobRows, 1));
    
    Matrix * incidence = NULL;
    if( distributedTestingCovariate != NULL )
    {
      incidence = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
      incidence->joinMatricesHorizontally(covariate, distributedTestingCovariate);
    }
    else
    {
      incidence = new Matrix(covariate);
    }
    
    GLMResults SNPResults;
    Matrix * fixedEffects = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    bool success = computeGLM(snpData, incidence, fixedEffects, SNPResults);

    //If success, collect results
    GLMResults reducedResults = GLMResults();;
    if(success == false)
    {
      SNPResults = GLMResults();
      SNPResults.success = false;
    }
    else{
      Matrix * br = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
      if(distributedTestingCovariate != NULL)
      {
        computeGLM(snpData, covariate, br, reducedResults);
      }
      else
      {
        computeGLM(snpData, covariate, br, reducedResults, true);
      }
      delete br;
    }
    if(communicator->mpiRoot)
    {
      bool succesGroupSignificance = computeGroupSignificance(snpData, reducedResults, SNPResults);
      
      computeGroupVariance(SNPResults, snpData, NULL, -1);
      
      int globalISNPdx = SNPidxs[ communicator->group ][isnp];
      groupedResults[ globalISNPdx ] = SNPResults;
    }
        
    delete snpData;    
    delete incidence;
    delete fixedEffects;
  }
  delete gTransposed;
  
  delete genotype;
  //delete phenotype;
  delete covariate;
  if( this->useCovariateMatrix == true )
  {
    delete this->currentCovarianceMatrix;
  }
  if( distributedTestingCovariate != NULL )
  {
    delete distributedTestingCovariate;
  }
  
  if(communicator->mpiRoot == false)
  {
    groupedResults.clear();
  }
  
  communicator = globalCommunicator;
  delete groupedCommunicator;
  
  if( this->useCovariateMatrix == true )
  {
    this->currentCovarianceMatrix = globalCovarianceMatrix;
  }
  
  misc.message << "100%" << std::endl;
  gatherResults(results, groupedResults, resultSNPInfo, gcgenotype, unfittedSNPs);
  
  storeResults(results, gccovariate, resultSNPInfo);
  
  //misc.message.output = &std::cout;
  //misc.message.log = misc.logFile;
  //misc.message <<
  //misc.message.output = &logOut;
  //misc.message.log = &logFile;
  
  if(communicator->mpiRoot && unfittedSNPs.size() != 0)
  {
    misc.message << "Warning: There are " << unfittedSNPs.size() << " SNPs which cannot be fitted. They are stored in file [ " << (options.outFile + ".gwas.unfitted") << " ]." << std::endl;
    Message message(options.outFile + ".gwas.unfitted");
    for(int i = 0; i < unfittedSNPs.size(); i++)
    {
      message << unfittedSNPs[i] << std::endl;
    }
  }
}


bool IGWAS::computeGLM(Matrix * y, Matrix * X, Matrix * b, GLMResults & results, bool onlyE)
{
  if( this->useCovariateMatrix == false )
  {
    return computeGLMWithoutCovariance(y, X, b, results);
  }
  else
  {
    return computeGLMWithCovariance(y, X, b, results, onlyE);
  }
}

bool IGWAS::computeGLMWithoutCovariance(Matrix * y, Matrix * X, Matrix * b, GLMResults & results)
{
  if(X->nGlobCols > X->nGlobRows)
  {
    misc.message << "Warning: There are more fixed effects than data points. I am unable to fit the model." << std::endl;
    results.success = false;
    return false;
  }
  
  //Adjust linear model
  Matrix * XtX_i = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  XtX_i->multiply(X, 'T', X, 'N');
  
  bool inverted = XtX_i->symmetricInvert();
  if(inverted == false)
  {
    //XtX_i->multiply(X, 'T', X, 'N');
    //inverted = XtX_i->invert();
    //if(inverted == false)
    //{
    //  return false;
    //}
    results = GLMResults();
    results.success = false;
    return false;
  }
  
  Matrix * Xty = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Xty->multiply(X, 'T', y, 'N');
  
  b->multiply(XtX_i, 'N', Xty, 'N');
  
  //Compute some matrices for SSE and t-statistics estimation
  Matrix * btXty = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  btXty->multiply(b, 'T', Xty, 'N');
  Matrix * yty = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  yty->multiply(y, 'T', y, 'N');
  
  if(btXty->nGlobRows != 1 || btXty->nGlobCols != 1 || yty->nGlobRows != 1 || yty->nGlobCols != 1)
  {
    misc.error("Error: An internal error was happened when computing GLM.", 0);
  }
  
  //Collect results
  std::vector<double> globalbtXty;
  std::vector<double> globalyty;
  std::vector<double> globalXtX_iDiagonal;
  results.SE.clear();
  results.tStatistics.clear();
  results.tStatisticPValues.clear();
  
  b->matrixToStandardVector(results.b);
  btXty->matrixToStandardVector(globalbtXty);
  yty->matrixToStandardVector(globalyty);
  globalXtX_iDiagonal = XtX_i->diagonal();
  
  if(communicator->mpiRoot)
  {
    double n_q_1 = double(y->nGlobRows) - double(b->nGlobRows);
    
    results.btXty = globalbtXty[0];
    
    results.SSE = globalyty[0] - results.btXty;
    results.MSE = results.SSE/double(n_q_1);

    for(int i = 0; i < results.b.size(); i++)
    {
      double temp = sqrt(results.MSE*globalXtX_iDiagonal[i]);
      results.SE.push_back(temp);
      results.tStatistics.push_back(results.b[i]/temp);
      results.tStatisticPValues.push_back( 2*tStatCDF(n_q_1, fabs(results.b[i]/temp) ) );
    }
  }
  
  delete Xty;
  delete XtX_i;
  delete btXty;
  delete yty;
  
  results.type = OLSModelType;
  results.success = true;
  
  return true;
}

bool IGWAS::computeGLMWithCovariance(Matrix * y, Matrix * X, Matrix * b, GLMResults & results, bool onlyE)
{
  if(X->nGlobCols > X->nGlobRows)
  {
    misc.message << "Warning: There are more fixed effects than data points. I am unable to fit the model." << std::endl;
    results.success = false;
    return false;
  }
  if( this->currentCovarianceMatrix == NULL )
  {
    misc.error("Error: An internal error was happened. Invalid GRM pointer when fitting the liner model in a GWAS.", 0);
  }
  
  REML reml;

  Matrix * remly = new Matrix(y);
  Matrix * remlX = new Matrix(X);
  Matrix * remlKernel = new Matrix(this->currentCovarianceMatrix);
  
  std::vector<Matrix*> kernels;
  kernels.push_back(remlKernel);
  
  std::vector<KernelType> kernelTypes;
  kernelTypes.push_back(kernelCovarianceMatrix);
  
  std::vector<double> weights;
  weights.push_back(1.);
  
  
  bool prepared = reml.prepare(remly, remlX, kernels, kernelTypes, options.initialh2Trait,  weights);
  if(onlyE == true)
  {
    int temp = reml.V->variances.size();
    std::set<std::string> tempIds;
    for(int j = 0; j < reml.V->elements.size(); j++)
    {
      tempIds.insert(reml.V->elements[j].subCovarianceId);
    }
    for(std::set<std::string>::iterator it = tempIds.begin(); it != tempIds.end(); ++it)
    {
      if( *it != "E" )
      {
        reml.V->deleteElementsWithSubCovarianceId( *it );
      }
    }
    reml.V->setInitialVariancesWithCurrentValues();
    if((temp - reml.V->variances.size()) != 1)
    {
      misc.error("Error: An internal error was happened. Invalid model and unexpected number of degrees of freedom.", 0);
    }
  }
  reml.useMLinsteadOfREML = true;
  /*if( this->nTests != 0 )
  {
    if( this->accumulatedVariance.size() != reml.V->variances.size() )
    {
      misc.error("Error: An internal error was happened. Unexpected variance structure discordance in GWAS REML fits.", 0);
    }
    for(int i = 0; i<reml.V->variances.size(); i++)
    {
      reml.V->variances[i].initialVariance = this->accumulatedVariance[ i ]/double(this->nTests);
    }
  }*/
  
  double logLikelihood = 0.;
  if( prepared == true )
  {
    double temp = options.changeAIStepThreshold;
    bool temp2 = options.firstStepEM;
    /*if( this->nTests > 10 )
    {
      options.changeAIStepThreshold = 1.;
      options.firstStepEM = false;
    }*/
    misc.message.active = options.verboseLog;
    
    logLikelihood = reml.computeREML();
    
    options.changeAIStepThreshold = temp;
    options.firstStepEM = temp2;
    misc.message.active = true;
  }
  else
  {
    misc.message << "Sorry, a problem was happened while preparing data for performing REML. The MLM cannot be fitted. Please, check the logs." << std::endl;
    results.success = false;
    return false;
  }
  if(reml.success == false)
  {
    results.success = false;
    return false;
  }
  
  Matrix * temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  temp->multiply(reml.ViX, 'T', reml.y, 'N');
  b->multiply(reml.XtViX_i, 'N', temp, 'N');
  delete temp;
  
  /*if(this->nTests == 0)
  {
    for(int i = 0; i<reml.V->variances.size(); i++)
    {
      this->accumulatedVariance.push_back(reml.V->variances[i].variance);
    }
  }
  else
  {
    if( this->accumulatedVariance.size() != reml.V->variances.size() )
    {
      misc.error("Error: An internal error was happened. Unexpected variance structure discordance in GWAS REML fits.", 0);
    }
    for(int i = 0; i<reml.V->variances.size(); i++)
    {
      this->accumulatedVariance[ i ] += reml.V->variances[i].variance;
    }
  }*/
  this->nTests++;
    
  results = GLMResults();
  
  b->matrixToStandardVector(results.b);
  std::vector<double> globalXtViX_iDiagonal = reml.XtViX_i->diagonal();
  
  if(communicator->mpiRoot)
  {
    for(int i = 0; i < results.b.size(); i++)
    {
      double temp = sqrt(globalXtViX_iDiagonal[i]);
      results.SE.push_back(temp);
      double chi2Test = (results.b[i]*results.b[i])/globalXtViX_iDiagonal[i];
      results.chi2Statistics.push_back(chi2Test);
      results.chi2StatisticsPValues.push_back( chi1_CDF(1, chi2Test) );
    }
    results.logLikelihood = logLikelihood;
  }
  
  if( reml.useMLinsteadOfREML == false )
  {
    results.type = REMLModelType;
  }
  else
  {
    results.type = MLModelType;
  }
  results.success = true;
  
  return true;
}

bool IGWAS::computeGroupSignificance(Matrix * y, GLMResults & reducedResults, GLMResults & results)
{
  if( results.type == OLSModelType )
  {
    if(reducedResults.success == false || results.success == false)
    {
      return false;
    }
    
    if( misc.gt( (results.btXty - reducedResults.btXty) < 0. ) )
    {
      return false;
    }
    
    if(communicator->mpiRoot)
    {
      double h = double(results.b.size()) - double(reducedResults.b.size());
      double n_q_1 = double(y->nGlobRows) - double(results.b.size());
      results.SSR = results.btXty - reducedResults.btXty;
      results.MSR = results.SSR/double(h);
      results.FStatistic = results.MSR / results.MSE;
      results.FStatisticPValue = FStatCDF(double(h), double(n_q_1), results.FStatistic);
    }
  }
  else if( results.type == MLModelType && reducedResults.type == MLModelType )
  {
    if(reducedResults.success == false || results.success == false)
    {
      return false;
    }
    
    if(communicator->mpiRoot)
    {
      double LogRatio = 2.0*(results.logLikelihood - reducedResults.logLikelihood);
      if(LogRatio < 0.)
      {
        LogRatio = 0.;
        results.FStatisticPValue = -1.;
      }
      else
      {
        //int df = int(results.b.size()) - int(reducedResults.b.size());
        int df;
        if( this->testingCovariate == NULL)
        {
          df = 1;
        }
        else
        {
          misc.error("Error: An internal error was happened.", 0);
        }
        results.FStatisticPValue = chi1_CDF(df, LogRatio);
      }
    }
  }
  else
  {
    return true;
  }
  
  return true;
}

void IGWAS::computeGroupVariance(GLMResults & results, Matrix * genotypes, LabeledMatrix * groupEffects, int idx)
{
  if( options.computeGroupVariance == false && options.saveGroupEffects == false )
  {
    return;
  }
  if( results.success == false )
  {
    if( groupEffects != NULL && options.saveGroupEffects == true )
    {
      Matrix * groupEffectMatrix = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, genotypes->nGlobRows, 1);
      groupEffectMatrix->fillWithConstant(0.);
      groupEffects->insertCol(groupEffectMatrix, idx);
    }
    
    return;
  }
  
  int nSNPs = genotypes->nGlobCols;
  Matrix * SNPEffectsMatrix = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nSNPs, 1);
  SNPEffectsMatrix->fillWithConstant(0.);
  std::vector<double> SNPEffects;
  if(communicator->mpiRoot == true)
  {
    if(nSNPs >= results.b.size())
    {
      misc.error("Error: An internal error was happened. Unexpected number of SNPs for computing the group variance.", 0);
    }
    int shift = results.b.size();
    shift -= nSNPs;
    SNPEffects = std::vector<double>(results.b.begin() + shift, results.b.end());
  }
  SNPEffectsMatrix->scatterMatrix(&(SNPEffects[0]));
  
  Matrix * groupEffectMatrix = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  if( this->useCovariateMatrix == true )
  {
    //If GRM is diagonal, genotypes are corrected
    if(this->currentCovariance->diagonalized == true)
    {
      Matrix * tempGenos = new Matrix();
      
      //tempGenos->multiply(genotypes, 'T', this->currentCovariance->eigenVectors, 'T');
      //groupEffectMatrix->multiply(tempGenos, 'T', SNPEffectsMatrix, 'N');
      
      //The previous multiplications can be trannslated to:
      tempGenos->multiply(genotypes, 'N', SNPEffectsMatrix, 'N');
      groupEffectMatrix->multiply(this->currentCovariance->eigenVectors, 'N', tempGenos, 'N');
      
      delete tempGenos;
    } //This correction can be moved to correcting the full matrix at the end of the analysis.
  }
  else
  {
    groupEffectMatrix->multiply(genotypes, 'N', SNPEffectsMatrix, 'N');
  }
  results.groupVariance = computeVariance(groupEffectMatrix);
  
  if( groupEffects != NULL && options.saveGroupEffects == true )
  {
    groupEffects->insertCol(groupEffectMatrix, idx);
  }
  
  delete SNPEffectsMatrix;
}

void IGWAS::storeResults(std::map<std::string, GLMResults > & effects, Covariate * covariate, std::map<std::string, std::vector<SNP> > & effectsSNPs)
{
  this->significantSNPIds.clear();
  
  if(communicator->mpiRoot)
  {
    Message messageMean(options.outFile + ".gwas.mean");
    Message messageDiscrete(options.outFile + ".gwas.discrete");
    Message messageQuantitative(options.outFile + ".gwas.quantitative");
    Message messageSNPs(options.outFile + ".gwas.snps");
    messageMean << "GROUP NAME BETA SE PV" << std::endl;
    messageDiscrete << "GROUP NAME BETA SE PV" << std::endl;
    messageQuantitative << "GROUP NAME BETA SE PV" << std::endl;
    messageSNPs << "GROUP SNP ALLELE MEAN STDEV BETA NBETA SE PV GROUPPV" << ((options.computeGroupVariance==false)?"":" GROUPVAR") << std::endl;
    
    for(std::map<std::string, GLMResults >::iterator it = effects.begin(); it != effects.end(); ++it)
    {
      std::string group = it->first;
      GLMResults groupResults = it->second;
      std::vector<double> globalEffects = groupResults.b;
      std::vector<double> globalPValues;
      if( groupResults.type == OLSModelType )
      {
        globalPValues = groupResults.tStatisticPValues;
      }
      else if( groupResults.type == REMLModelType || groupResults.type == MLModelType )
      {
        globalPValues = groupResults.chi2StatisticsPValues;
      }
      else
      {
        misc.error("Error: An internal error was happened. Invalid model type.", 0);
      }
      std::vector<SNP> SNPs = effectsSNPs[group];
      
      int nTestingCategories = 0;
      if(this->testingCovariate != NULL)
      {
        nTestingCategories = this->testingCovariate->meanNames.size() + this->testingCovariate->discreteCovarNames.size() + this->testingCovariate->quantitativeCovarNames.size();
      }
      if( globalEffects.size() != (covariate->meanNames.size() + covariate->discreteCovarNames.size() + covariate->quantitativeCovarNames.size() + nTestingCategories) )
      {
        misc.error("Error: An internal error was happened. The size of the covariate names does not have the same dimension of fitted model fixed effects.", 0);
      }
      int shift = 0;
      
      for(int i = 0; i < covariate->meanNames.size(); i++)
      {
        messageMean << group << " " << covariate->meanNames[i];
        messageMean << " " << globalEffects[ i + shift ];
        messageMean << " " << groupResults.SE[ i + shift ];
        messageMean << " " << globalPValues[ i + shift ];
        messageMean << std::endl;
      }
      shift += covariate->meanNames.size();
      
      for(int i = 0; i < covariate->discreteCovarNames.size(); i++)
      {
        messageDiscrete << group << " " << covariate->discreteCovarNames[i];
        messageDiscrete << " " << globalEffects[ i + shift ];
        messageDiscrete << " " << groupResults.SE[ i + shift ];
        messageDiscrete << " " << globalPValues[ i + shift ];
        messageDiscrete << std::endl;
      }
      shift += covariate->discreteCovarNames.size();
      
      for(int i = 0; i < covariate->quantitativeCovarNames.size(); i++)
      {
        messageQuantitative << group << " " << covariate->quantitativeCovarNames[i];
        messageQuantitative << " " << globalEffects[ i + shift ];
        messageQuantitative << " " << groupResults.SE[ i + shift ];
        messageQuantitative << " " << globalPValues[ i + shift ];
        messageQuantitative << std::endl;
      }
      shift += covariate->quantitativeCovarNames.size();
      
      for(int i = 0; i < SNPs.size(); i++)
      {
        messageSNPs << group;
        messageSNPs << " " << SNPs[i].name;
        messageSNPs << " " << SNPs[i].allele2;
        messageSNPs << " " << std::setprecision(3) << 2.*SNPs[i].p2;
        messageSNPs << " " << std::setprecision(3) << SNPs[i].standardDev;
        messageSNPs << " " << "NA"; //globalEffects[ i + shift ];
        messageSNPs << " " << "NA"; //std::setprecision(5) << globalEffects[ i + shift ]/SNPs[i].standardDev;
        messageSNPs << " " << "NA"; //groupResults.SE[ i + shift ];
        messageSNPs << " " << "NA"; //globalPValues[ i + shift ];
        messageSNPs << " " << groupResults.FStatisticPValue;
        if( options.computeGroupVariance == true )
        {
          messageSNPs << " " << groupResults.groupVariance;
        }
        messageSNPs << std::endl;
        
        /*if(globalPValues[ i + shift ] < this->significanceThreshold)
        {
          this->significantSNPIds.push_back(SNPs[i].name);
        }*/
      }
    }
  }

  if(communicator->mpiRoot)
  {
    if(this->correlatedSNPIds.size() != 0)
    {
      std::vector<std::string> temp(this->correlatedSNPIds.begin(), this->correlatedSNPIds.end());
      std::vector<std::string> intersection = intersectionStringVectors(2, &(this->significantSNPIds), &temp);
      if(intersection.size() != 0)
      {
        std::vector<std::string> uncorrelatedSignifcantSNPIds = differenceBetweenTwoVectors(this->significantSNPIds, intersection);
        misc.message << this->significantSNPIds.size() - uncorrelatedSignifcantSNPIds.size() << " correlated SNPs removed." << std::endl;
        this->significantSNPIds = uncorrelatedSignifcantSNPIds;
        
        Message messageCorrelated(options.outFile + ".gwas.correlatedSNPs");
        for(int i = 0; i < intersection.size(); i++)
        {
          messageCorrelated << intersection[i] << std::endl;
        }
      }
      this->correlatedSNPIds.clear();
    }
  }
}

void IGWAS::getLessSignificantCorrelatedSNPs(double threshold, Genotype* genotypes, GLMResults results, int shift)
{
  if(threshold <= 0.)
  {
    return;
  }
  
  Matrix * correlations = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Matrix * normalization = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  correlations->multiply(genotypes->genotypes, 'N', genotypes->genotypes, 'T');
  normalization->multiply(genotypes->missings, 'N', genotypes->missings, 'T');
  correlations->elementWiseDivision(normalization);
  delete normalization;
  
  std::vector<int> idxSNPs1;
  std::vector<int> idxSNPs2;
  correlations->getGlobalIndexOutsideRange(-threshold, threshold, idxSNPs1, idxSNPs2);

  std::vector<std::string> correlatedSNPs;
  if(communicator->mpiRoot == true)
  {
    if( (correlations->nGlobRows + shift) != results.b.size())
    {
      misc.error("Error: An internal error was happened. Discordant number of elements when searching for uncorrelated SNPs.", 0);
    }
    std::set<int> idxsToRemove;
    for(int i = 0; i<idxSNPs1.size(); i++)
    {
      int idxS1 = idxSNPs1[i];
      int idxS2 = idxSNPs2[i];
      
      if( (idxS1 + shift) >= results.tStatisticPValues.size() || (idxS2 + shift) >= results.tStatisticPValues.size() )
      {
        misc.error("Error: An internal error was happened. Incorrect indices.", 0);
      }
      
      if( idxS1 > idxS2 ) //Because the matrix is symmetric.
      {
        if( results.tStatisticPValues[ idxS1 + shift ] < results.tStatisticPValues[ idxS2 + shift ] )
        {
          idxsToRemove.insert(idxS2);
        }
        else
        {
          idxsToRemove.insert(idxS1);
        }
      }
    }
    
    for (std::set<int>::iterator it = idxsToRemove.begin(); it != idxsToRemove.end(); ++it)
    {
      this->correlatedSNPIds.insert( genotypes->SNPIds[*it] );
    }
  }
  
  delete correlations;
}

void IGWAS::gatherResults(std::map<std::string, GLMResults > & results, std::map<int, GLMResults > & groupedResults, std::map<std::string, std::vector<SNP> > & resultSNPInfo, Genotype * genotype, std::vector<std::string> & unfittedSNPs)
{
  results.clear();
  resultSNPInfo.clear();
  
  std::vector<int> SNPIdxs;
  std::vector<double> b;
  std::vector<double> SE;
  std::vector<double> PValues;
  std::vector<int> modelType; // 0: OLSModelType, 1: REMLModelType, 2: MLModelType
  std::vector<double> FStatisticPValue;
  std::vector<double> groupVariance;
  
  for(std::map<int, GLMResults >::iterator it = groupedResults.begin(); it != groupedResults.end(); ++it)
  {
    int gsi = it->first; //global SNP index
    if(it->second.success == false)
    {
      continue;
    }
    for(int i = 0; i<it->second.b.size(); i++)
    {
      SNPIdxs.push_back(gsi);
      b.push_back(it->second.b[i]);
      SE.push_back(it->second.SE[i]);
      if( it->second.type == OLSModelType )
      {
        PValues.push_back(it->second.tStatisticPValues[i]);
      }
      else if( it->second.type == REMLModelType || it->second.type == MLModelType )
      {
        PValues.push_back(it->second.chi2StatisticsPValues[i]);
      }
      else
      {
        misc.error("Error: An internal error was happened. Invalid model type when gathering results (1).", 0);
      }
    }
    if( it->second.type == OLSModelType )
    {
      modelType.push_back(0);
    }
    else if( it->second.type == REMLModelType )
    {
      modelType.push_back(1);
    }
    else if( it->second.type == MLModelType )
    {
      modelType.push_back(2);
    }
    else
    {
      misc.error("Error: An internal error was happened. Invalid model type when gathering results (2).", 0);
    }
    FStatisticPValue.push_back(it->second.FStatisticPValue);
    groupVariance.push_back(it->second.groupVariance);
  }
  
  int * globSNPIdxs;
  double * globb;
  double * globSE;
  double * globPValues;
  int * globmodelType;
  double * globFStatisticPValue;
  double * globGroupVariance;
  
  int vectorsSize;
  int scalarsSize;
  globSNPIdxs = communicator->asymmetricGather(&(SNPIdxs[0]), SNPIdxs.size(), &vectorsSize);
  globb = communicator->asymmetricGather(&(b[0]), b.size(), &vectorsSize);
  globSE = communicator->asymmetricGather(&(SE[0]), SE.size(), &vectorsSize);
  globPValues = communicator->asymmetricGather(&(PValues[0]), PValues.size(), &vectorsSize);
  globmodelType = communicator->asymmetricGather(&(modelType[0]), modelType.size(), &scalarsSize);
  globFStatisticPValue = communicator->asymmetricGather(&(FStatisticPValue[0]), FStatisticPValue.size(), &scalarsSize);
  globGroupVariance = communicator->asymmetricGather(&(groupVariance[0]), groupVariance.size(), &scalarsSize);
  
  if(communicator->mpiRoot && vectorsSize != 0)
  {
    for(int isnp = 0; isnp<genotype->nSNPs; isnp++)
    {
      results[ genotype->SNPs[isnp].name ] = GLMResults();
    }
    
    int iPrevSNP = -1;
    int prevNVars = -1;
    int nVars = -1;
    int globSNPVarGroupIdx = 0; //The index of the number of SNPs already gathered. Used for indexing scalar variables.
    for(int i = 0; i < vectorsSize; i++)
    {
      int isnp = globSNPIdxs[i];
      std::string snpname = genotype->SNPs[isnp].name;
      
      if( iPrevSNP != isnp )
      {
        if( results[snpname].b.size() != 0 && iPrevSNP != -1 )
        {
          misc.error("Error: An internal error was happened while gathering SNP information from GWAS done in a grouped communicator. Repeated SNP.", 0);
        }
        if( prevNVars != nVars && prevNVars != -1)
        {
          misc.error("Error: An internal error was happened while gathering SNP information from GWAS done in a grouped communicator. Inconsistent number of vars (1): " + i2s(prevNVars) + " " + i2s(nVars) + " " + i2s(i), 0);
        }
        if( globSNPVarGroupIdx >= scalarsSize )
        {
          misc.error("Error: An internal error was happened while gathering SNP information from GWAS done in a grouped communicator. Inconsistent number of vars (2).", 0);
        }
        
        results[ snpname ].FStatisticPValue = globFStatisticPValue[globSNPVarGroupIdx];
        results[ snpname ].groupVariance = globGroupVariance[globSNPVarGroupIdx];
        if( globmodelType[globSNPVarGroupIdx] == 0 )
        {
          results[ snpname ].type = OLSModelType;
        }
        else if( globmodelType[globSNPVarGroupIdx] == 1 )
        {
          results[ snpname ].type = REMLModelType;
        }
        else if( globmodelType[globSNPVarGroupIdx] == 2 )
        {
          results[ snpname ].type = MLModelType;
        }
        else
        {
          misc.error("Error: An internal error was happened. Invalid model type when gathering results (3).", 0);
        }
        results[ snpname ].success = true;
        
        globSNPVarGroupIdx++;
        prevNVars = nVars;
        nVars = 0;
      }
      
      results[ snpname ].b.push_back(globb[i]);
      results[ snpname ].SE.push_back(globSE[i]);
      if( results[ snpname ].type == OLSModelType )
      {
        results[ snpname ].tStatisticPValues.push_back(globPValues[i]);
      }
      else if( results[ snpname ].type == REMLModelType || results[ snpname ].type == MLModelType )
      {
        results[ snpname ].chi2StatisticsPValues.push_back(globPValues[i]);
      }
      else
      {
        misc.error("Error: An internal error was happened. Invalid model type when gathering results (4).", 0);
      }
      
      iPrevSNP = isnp;
      nVars++;
    }
    
    if( nVars*scalarsSize != vectorsSize )
    {
      misc.error("Error: An internal error was happened while gathering SNP information from GWAS done in a grouped communicator. Inconsistent number of vars (3).", 0);
    }
      
    delete [] globSNPIdxs;
    delete [] globb;
    delete [] globSE;
    delete [] globPValues;
    delete [] globmodelType;
    delete [] globFStatisticPValue;
    delete [] globGroupVariance;
    
    
    
    unfittedSNPs.clear();
    for(int isnp = 0; isnp<genotype->nSNPs; isnp++)
    {
      std::string snpname = genotype->SNPs[isnp].name;
      if( results[ snpname ].success == true )
      {
        std::vector<SNP> temp;
        temp.push_back(genotype->SNPs[isnp]);
        resultSNPInfo[ snpname ] = temp;
      }
      else
      {
        results.erase(snpname);
        unfittedSNPs.push_back(snpname);
      }
    }
  }

}

bool IGWAS::computeCovariance(Genotype * genotype, Phenotype * phenotype, Covariate * covariate)
{
  std::vector<Kernel*> kernels;
  
  bool newGRM = true; //This has to be improved, if there are no GRMs, then it can be optimized to not recompute the covariance matrix each time.
  if(this->grmFiles.count( this->currentFile ) != 0)
  {
    if( this->grmFiles[this->currentFile] != this->currentGRMFile ) //Only load a new GRM if file changes.
    {
      if(this->currentGRMBase != NULL)
      {
        delete this->currentGRMBase;
      }
      this->currentGRMFile = this->grmFiles[this->currentFile];
      this->currentGRMBase = new Kernel(this->currentGRMFile);
      this->currentGRMBase->normalize();
      this->currentGRMBase->name = "GRM";
      
      if( this->currentCovariance != NULL )
      {
        delete this->currentCovariance;
        this->currentCovariance = NULL;
        this->currentCovarianceMatrix = NULL;
      }
    }
    else
    {
      newGRM = false;
    }
    
    Kernel *temp = new Kernel(this->currentGRMBase);
    kernels.push_back(temp);
  }
  else
  {
    if( this->currentGRMFile == "" )
    {
      newGRM = false;
    }
    else
    {
      newGRM == true;
      this->currentGRMFile = "";
      delete this->currentGRMBase;
      this->currentGRMBase = NULL;
    }
  }
  
  //If the GRM have not changed, the other random effects are the same, and the individuals are the same, skip recomputing the covariance.
  if( this->currentCovariance != NULL && newGRM == false )
  {
    std::vector<std::string> testIndividualsInGenotypeOrder = orderVectorAsTemplate(genotype->individualIds, this->currentCovariance->individualIds);
    if( testIndividualsInGenotypeOrder == this->currentCovariance->individualIds )
    {
      misc.message << "Skipping covariance computation and decomposition. Using the previous one..." << std::endl;
      for(int iker = 0; iker<kernels.size(); iker++)
      {
        delete kernels[iker];
      }
      return true;
    }
  }
  
  //Add other random effects if specified.
  std::vector<std::string> dummyReducedModels;
  std::vector<std::string> dummyIndividualBLUPNames;
  addKernelsUsingOptions(kernels, dummyReducedModels, dummyIndividualBLUPNames);
  
  //Get shared individuals and filter
  std::vector<std::string> commonIndividuals = intersectionStringVectors(2, &genotype->individualIds, &covariate->individualIds);
  for(int iker = 0; iker<kernels.size(); iker++)
  {
    commonIndividuals = intersectionStringVectors(2, &commonIndividuals, &(kernels[iker]->individualIds));
  }
  std::vector<std::string> commonIndividualsInGenotypeOrder = orderVectorAsTemplate(genotype->individualIds, commonIndividuals);
  for(int iker = 0; iker<kernels.size(); iker++)
  {
    if( kernels[iker]->diagonalized == true && commonIndividualsInGenotypeOrder != kernels[iker]->individualIds )
    {
      misc.error("Error: When using diagonal GRMs, the individuals in phenotype, genotype, covars and GRM files must be all the same without missing phenotypes. Aborting.", 0);
    }
    kernels[iker]->normalize();
    if( kernels[iker]->diagonalized == true )
    {
      kernels[iker]->filterIndividuals(commonIndividualsInGenotypeOrder, false);
    }
    else
    {
      kernels[iker]->filterIndividuals(commonIndividualsInGenotypeOrder); //if it is not diagonal, N has to be kept for Kernel method sanitizeKernel().
    }
  }
  
  if(kernels.size() == 0)
  {
    misc.error("Error: An internal error has happened. There are no kernels for computing the covariance matrix.", 0);
  }
  
  //Clear previous covariance matrix and create a copy of the parameters of the first kernel.
  if( this->currentCovariance != NULL )
  {
    delete this->currentCovariance;
    this->currentCovariance = NULL;
    this->currentCovarianceMatrix = NULL;
  }
  
  //Avoid REML if there is only one matrix
  if(kernels.size() == 1)
  {
    this->currentCovariance = kernels[0];
    if( options.gwasUseAlwaysDiagonalCovariances == true )
    {
      this->currentCovariance->diagonalizeKernel();
    }
    this->currentCovarianceMatrix = this->currentCovariance->getNormalizedKernel();
    
    kernels.clear();
    
    return true;
  }
  
  misc.error("Error: Internal error. This part is currently not valid in a igwas. It has to be worked more on the options and this part of code.", 0);
  
  //Prepare, run REML, and compute covariance matrix
  this->currentCovariance = new Kernel();
  this->currentCovariance->copyParameters(kernels[0]);
  
  std::vector<double> weights;
  weights.assign( kernels.size() , 1./double(kernels.size()) );
  
  std::vector<int> phenotypeColumns;
  phenotypeColumns.push_back(options.phenotypeColumn);
  if(phenotypeColumns.size() != 1)
  {
    misc.error("Error: An internal error was happened. Unexpected number of phenotypes.", 0);
  }

  std::vector<double> heritabilities;
  heritabilities.push_back(options.initialh2Trait);
  
  std::vector<std::pair<std::string, std::string> > covariateFiles;
  covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFile, options.qCovarsFile));
  
  REML * reml = new REML(false);

  bool prepared = reml->prepare(singleREMLType, kernels, weights, phenotypeColumns, heritabilities, covariateFiles);
  kernels.clear();
  if( prepared == true )
  {
    reml->computeREML();
  }
  else
  {
    delete reml;
    delete this->currentCovariance;
    this->currentCovariance = NULL;
    misc.message << "An error happened preparing the REML analysis, variances cannot be computed." << std::endl;
    return false;
  }
  
  if( reml->success == false )
  {
    delete reml;
    delete this->currentCovariance;
    this->currentCovariance = NULL;
    misc.message << "REML did not converge, variances cannot be computed." << std::endl;
    return false;
  }
  
  reml->V->computeCovariance();
  
  if( reml->V->m->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Unexpected diagonal covariance matrix.", 0);
  }
  if( misc.gt(reml->V->m->nGlobRows != this->currentCovariance->nIndividuals || reml->V->m->nGlobCols != this->currentCovariance->nIndividuals) )
  {
    misc.error("Error: An internal error was happened. Unexpected covariance matrix dimensions.", 0);
  }
  this->currentCovariance->kernel = new Matrix(reml->V->m);
  this->currentCovariance->type = kernelCovarianceMatrix;
  this->currentCovariance->normalized = true;
  this->currentCovariance->diagonalized = false;
  this->currentCovariance->asymmetric = false;

  double variancesTotal = 0.;
  for(int i=0; i<reml->V->variances.size(); i++)
  {
    if( reml->V->variances[i].typeEffect != ParameterAttributes::environment && reml->V->variances[i].typeEffect != ParameterAttributes::other )
    {
      variancesTotal += reml->V->variances[i].variance;
    }
  }
  communicator->broadcast(&variancesTotal);
  
  this->currentCovariance->kernel->scaleBy(1./variancesTotal);
  if( options.gwasUseAlwaysDiagonalCovariances == true )
  {
    this->currentCovariance->diagonalizeKernel();
  }
  this->currentCovarianceMatrix = this->currentCovariance->getNormalizedKernel();
  
  delete reml;
  
  return true;
}

void IGWAS::debugWrite(Matrix * y, Matrix * X, Matrix * b)
{
  if(options.debug)
  {
    y->debugWrite(options.outFile + ".y");
    X->debugWrite(options.outFile + ".X");
    b->debugWrite(options.outFile + ".b");
  }
}