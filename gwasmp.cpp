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

#include "gwas.h"
#include "analysis.h"
#include "reml.h"
#include "auxiliar.h"
#include "labeledmatrix.h"
#include "mpresiduals.h"

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream> 
#include <cmath>
#include <cstdlib>
#include <cstdio>

void GWAS::loadGenotypeResidualFiles()
{
  std::map<std::string, std::vector<std::string> > residuals2genotypes;
  std::string usedFile;
  if( options.genotypeAndResidualsListFile != "" )
  {
    usedFile = options.genotypeAndResidualsListFile;
    std::vector< std::vector<std::string> > files;
    getTableFromFile(options.genotypeAndResidualsListFile, files, 2);
    for(int i = 0; i<files.size(); i++) //We load genotype and grm files in this way, because we want all genotype with same grms together in the list.
    {
      misc.checkFileExists(files[i][0] + ".bed");
      misc.checkFileExists(files[i][0] + ".bim");
      misc.checkFileExists(files[i][0] + ".fam");
      misc.checkFileExists(files[i][1] + ".rowids");
      misc.checkFileExists(files[i][1] + ".colids");
      misc.checkFileExists(files[i][1] + ".dat");
      residuals2genotypes[ files[i][1] ].push_back( files[i][0] );
    }
  }
  else if(options.genotypeBGENAndResidualsListFile != "")
  {
    usedFile = options.genotypeBGENAndResidualsListFile;
    this->genotypeFilesType = GenotypeAttributes::probabilities;
    std::vector< std::vector<std::string> > files;
    getTableFromFile(options.genotypeBGENAndResidualsListFile, files, 2);
    for(int i = 0; i<files.size(); i++) //We load genotype and grm files in this way, because we want all genotype with same grms together in the list.
    {
      misc.checkFileExists(files[i][0] + ".bgen");
      misc.checkFileExists(files[i][1] + ".rowids");
      misc.checkFileExists(files[i][1] + ".colids");
      misc.checkFileExists(files[i][1] + ".dat");
      residuals2genotypes[ files[i][1] ].push_back( files[i][0] );
    }
  }
  else
  {
    misc.error("Error: A file containing a list of genotype and residuals files have to be specified.", 0);
  }
  
  
  for( std::map<std::string, std::vector<std::string> >::iterator it = residuals2genotypes.begin(); it != residuals2genotypes.end(); ++it )
  {
    std::vector<std::string> genoFiles = it->second;
    std::string grmFile = it->first;
    for(int i = 0; i<genoFiles.size(); i++)
    {
      if( this->residualFiles.count(genoFiles[i]) != 0 )
      {
        misc.error("Error: There is at least one genotype file repeated in [ " + usedFile + " ] file. Genotype files have to be unique.", 0);
      }
      this->genotypeFiles.push_back(genoFiles[i]);
      this->residualFiles[ genoFiles[i] ] = grmFile;
    }
  }

}

void GWAS::computeMultiplePhenotypeGWAS()
{
  std::string backupOutFile = options.outFile;
  if(this->genotypeFiles.size() > 1)
  {
    misc.message.tab = "  ";
  }
  
  LabeledMatrix * residuals = NULL;
  Matrix * groupedCommunicatorDistributedResiduals = NULL;
  for(int i = 0; i<this->genotypeFiles.size(); i++)
  {
    if( this->genotypeFiles.size() > 1 )
    {
      options.outFile += "." + getFileName(this->genotypeFiles[i]);
    }
    
    //Load genotypes
    this->currentFile = this->genotypeFiles[i];
    Genotype * genotype = NULL;
    if(this->genotypeFilesType == GenotypeAttributes::calls)
    {
      genotype = new Genotype(this->genotypeFiles[i]);
      genotype->normalizeGenotypes();
    }
    else
    {
      std::set<std::string> keepSNPs;
      std::set<std::string> keepIndividualIds = getIndividualIdsSpecifiedByOptionKeep();
      genotype = new Genotype(this->genotypeFiles[i], keepSNPs, keepIndividualIds, GenotypeAttributes::probabilities);
      genotype->normalizeGenotypes();
    }
    
    if(this->residualFiles.count( this->currentFile ) == 0)
    {
      misc.error("Error: Currently only allowed with covariance matrix. Work in progress. If you are interested on this type of analysis. Please, contact us.", 0);
    }
    else
    {
      if( this->residualFiles[this->currentFile] != this->currentResidualFile ) //Only load a new residual if file changes.
      {
        this->currentResidualFile = this->residualFiles[this->currentFile];
        if(residuals != NULL)
        {
          delete residuals;
          residuals = NULL;
        }
        residuals = new LabeledMatrix(this->currentResidualFile);
        residuals->matrix->centerMatrixRowsColumns(column);
        
        if(misc.gt(residuals->rowLabels.size() < 2))
        {
          misc.error("Error: There are less than 2 individuals for running the analysis. Please, check the overlap between covariates, phenotypes, and GRMs.", 0);
        }
      }
      else
      {
        if(residuals == NULL)
        {
          misc.error("Error: An internal error has happened. The residuals are missing.", 0);
        }
      }
    }
    std::vector<std::string> commonIndividualsInGenotypeOrder = orderVectorAsTemplate(genotype->individualIds, residuals->rowLabels);
    if( commonIndividualsInGenotypeOrder != residuals->rowLabels )
    {
      misc.error("Error: When performing a multiple phenotype analysis, all the genotype files, corrected by the same GRM, have to include the same individuals in the same order.", 0);
    }
    genotype->filterSNPsAndIndividuals(genotype->SNPIds, commonIndividualsInGenotypeOrder, false);
    if( genotype->individualIds != residuals->rowLabels )
    {
      misc.error("Error: When performing a multiple phenotype analysis, all the genotype files, corrected by the same GRM have to include the same individuals in the same order.", 0);
    }
    
    //Start analysis
    if(options.analysis != recursiveGWASAnalysis) //Perform grouped GWAS or standard GWAS
    {
      misc.setGetElapsedTime("GWAS");
      misc.message << "Starting analysis..." << std::endl;

      computeIndividualGWASMultiplePhenotypesGroupedCommunicator(genotype, residuals); //, Matrix ** groupedCommunicatorDistributedResiduals

      misc.message << "Analysis finished after " << misc.setGetElapsedTime("GWAS", true) << "." << std::endl;
    }
    
    options.outFile = backupOutFile;
    
    delete genotype;
  }
  if(residuals != NULL)
  {
    delete residuals;
    residuals = NULL;
  }
  
  this->currentFile = "";
  
  if(this->currentGRMBase != NULL)
  {
    delete this->currentGRMBase;
    this->currentGRMBase = NULL;
  }
  
  if(this->genotypeFiles.size() > 1)
  {
    misc.message.tab = "";
    misc.message << "GWAS analysis finished on all files." << std::endl;
  }
}

void GWAS::computeMultiplePhenotypeGWASGroupedCommunicator()
{
  std::string backupOutFile = options.outFile;
  if(this->genotypeFiles.size() > 1)
  {
    misc.message.tab = "  ";
  }
  
  Communicator * globalCommunicator = communicator;  
  Communicator * groupedCommunicator = new Communicator(globalCommunicator, basicGroupedCommunicator);
  communicator = groupedCommunicator;
  
  LabeledMatrix * residuals = NULL;
  Matrix * groupedCommunicatorDistributedResiduals = NULL;
  for(int i = 0; i<this->genotypeFiles.size(); i++)
  {
    if( this->genotypeFiles.size() > 1 )
    {
      options.outFile += "." + getFileName(this->genotypeFiles[i]);
    }
    
    //Get SNPs to keep for this group.
    this->currentFile = this->genotypeFiles[i];
    std::vector<std::string> allSNPIds;
    if(this->genotypeFilesType == GenotypeAttributes::calls)
    {
      Genotype loadSNPs;
      std::set<std::string> keepSNPs;
      std::vector<int> SNPsBEDIdxsToKeep;
      loadSNPs.readBIMFile(this->genotypeFiles[i], keepSNPs, SNPsBEDIdxsToKeep);
      allSNPIds = loadSNPs.SNPIds;
    }
    else
    {
#if defined(BGEN) && defined(ZLIB)
      Genotype loadSNPs;
      allSNPIds = loadSNPs.getBGENFileSNPs(this->genotypeFiles[i]);
#else
      misc.error("Error: An internal error has happened. The current version of DISSECT is compiled without zlib and bgen support. These libraries are required for using this analysis with the current inputs.", 0);
#endif
    }
    int allNSNPs = allSNPIds.size();
    int nSNPsBaseInGroup = allNSNPs/groupedCommunicator->nGroups;
    int remainingSNPs = allNSNPs - (nSNPsBaseInGroup*groupedCommunicator->nGroups);
    if(remainingSNPs < 0 || remainingSNPs >= groupedCommunicator->nGroups)
    {
      misc.error("Error: An internal error was happened. Unexpected remainder when distributing SNPs between grouped processes.", 0);
    }
    std::set<std::string> currentGroupSNPIds;
    int isShift = 0;
    for(int ig = 0; ig < groupedCommunicator->nGroups; ig++)
    {
      int nSNPsInGroup = nSNPsBaseInGroup;
      if( ig < remainingSNPs )
      {
        nSNPsInGroup++;
      }
      if(ig == communicator->group)
      {
        for(int is = 0; is < nSNPsInGroup; is++)
        {
          int gsi = is + isShift; //Global SNP index.
          currentGroupSNPIds.insert(allSNPIds[gsi]);
        }
      }
      isShift += nSNPsInGroup;
    }
    if(isShift != allNSNPs)
    {
      misc.error("Error: An internal error has happened. Not all SNPs have been redistributed between grouped processes.", 0);
    }
    
    //Load genotypes
    Genotype * genotype = NULL;
    std::set<std::string> keepIndividualIds = getIndividualIdsSpecifiedByOptionKeep();
    if(this->genotypeFilesType == GenotypeAttributes::calls)
    {
      genotype = new Genotype(this->genotypeFiles[i], currentGroupSNPIds, keepIndividualIds);
      genotype->normalizeGenotypes();
    }
    else
    {
      genotype = new Genotype(this->genotypeFiles[i], currentGroupSNPIds, keepIndividualIds, GenotypeAttributes::probabilities);
      genotype->normalizeGenotypes();
    }
    
    if(this->residualFiles.count( this->currentFile ) == 0)
    {
      misc.error("Error: Currently only allowed with covariance matrix. Work in progress. If you are interested on this type of analysis. Please, contact us.", 0);
    }
    else
    {
      if( this->residualFiles[this->currentFile] != this->currentResidualFile ) //Only load a new residual if file changes.
      {
        this->currentResidualFile = this->residualFiles[this->currentFile];
        if(residuals != NULL)
        {
          delete residuals;
          residuals = NULL;
        }
        residuals = new LabeledMatrix(this->currentResidualFile);
        residuals->matrix->centerMatrixRowsColumns(column);
        
        if(misc.gt(residuals->rowLabels.size() < 2))
        {
          misc.error("Error: There are less than 2 individuals for running the analysis. Please, check the overlap between covariates, phenotypes, and GRMs.", 0);
        }
      }
      else
      {
        if(residuals == NULL)
        {
          misc.error("Error: An internal error has happened. The residuals are missing.", 0);
        }
      }
    }
    std::vector<std::string> commonIndividualsInGenotypeOrder = orderVectorAsTemplate(genotype->individualIds, residuals->rowLabels);
    if( commonIndividualsInGenotypeOrder != residuals->rowLabels )
    {
      misc.error("Error: When performing a multiple phenotype analysis, all the genotype files, corrected by the same GRM, have to include the same individuals in the same order.", 0);
    }
    genotype->filterSNPsAndIndividuals(genotype->SNPIds, commonIndividualsInGenotypeOrder, false);
    if( genotype->individualIds != residuals->rowLabels )
    {
      misc.error("Error: When performing a multiple phenotype analysis, all the genotype files, corrected by the same GRM have to include the same individuals in the same order.", 0);
    }
    
    //Start analysis
    if(options.analysis != recursiveGWASAnalysis) //Perform grouped GWAS or standard GWAS
    {
      misc.setGetElapsedTime("GWAS");
      misc.message << "Starting analysis..." << std::endl;

      computeIndividualGWASMultiplePhenotypes(genotype, residuals, globalCommunicator); //, Matrix ** groupedCommunicatorDistributedResiduals

      misc.message << "Analysis finished after " << misc.setGetElapsedTime("GWAS", true) << "." << std::endl;
    }
    
    options.outFile = backupOutFile;
    
    delete genotype;
  }
  if(residuals != NULL)
  {
    delete residuals;
    residuals = NULL;
  }
  
  this->currentFile = "";
  
  communicator = globalCommunicator;
  delete groupedCommunicator;
  
  if(this->genotypeFiles.size() > 1)
  {
    misc.message.tab = "";
    misc.message << "GWAS analysis finished on all files." << std::endl;
  }
}

void GWAS::computeIndividualGWASMultiplePhenotypesGroupedCommunicator(Genotype * gcgenotype, LabeledMatrix * gcphenotypes)
{
  std::map<std::string, GLMResults > results;
  std::map<int, GLMResults > groupedResults;
  std::map<std::string, std::vector<SNP> > resultSNPInfo;
  std::vector<std::string> unfittedSNPs;
  
  int nSNPs = gcgenotype->nSNPs;

  Communicator * globalCommunicator = communicator;  
  Communicator * groupedCommunicator = new Communicator(globalCommunicator, basicGroupedCommunicator);
 
  
  std::map< int, std::vector<int> > SNPidxs;

  //Compute some matrices which whill be conserved.
  misc.setGetElapsedTime("ytys");
  misc.message << "Computing yty elements..." << std::endl;
  Matrix * temp = new Matrix(gcphenotypes->matrix);
  temp->elementWiseMultiplication(gcphenotypes->matrix);
  Matrix * mOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, temp->nGlobRows, 1);
  mOnes->fillWithConstant(1.);
  Matrix * gcytys = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  gcytys->multiply(temp, 'T', mOnes, 'N'); //1 column matrix, each row should contain the yty for each phenotype.
  delete temp;
  delete mOnes;
  misc.message << "Computation finished after " << misc.setGetElapsedTime("ytys", true) << "." << std::endl;
  
  
  //Redistribute data between different communicators.
  misc.setGetElapsedTime("redistribute");
  misc.message << "Redistributing data between processes..." << std::endl;
  
  Matrix * genotype = gcgenotype->genotypesRedistributionToGroupedCommunicatorMatrices(groupedCommunicator, SNPidxs);
  Matrix * phenotypes = gcphenotypes->matrix->copyToGroupedCommunicator(groupedCommunicator);
  Matrix * ytys = gcytys->copyToGroupedCommunicator(groupedCommunicator);
  
  Matrix * globalCovarianceMatrix = NULL;
  if( this->useCovariateMatrix == true )
  {
    globalCovarianceMatrix = this->currentCovarianceMatrix;
    this->currentCovarianceMatrix = this->currentCovarianceMatrix->copyToGroupedCommunicator(groupedCommunicator);
  }
  
  misc.message << "Data redistributed after " << misc.setGetElapsedTime("redistribute", true) << "." << std::endl;
  
  communicator = groupedCommunicator;
  
  
  //Compute matrices needed after:
  if( communicator->group == 0 )
  {
    misc.setGetElapsedTime("XtXs");
    misc.message << "Computing XtX elements..." << std::endl;
  }
  temp = new Matrix(genotype);
  temp->elementWiseMultiplication(genotype);
  mOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, temp->nGlobCols, 1);
  mOnes->fillWithConstant(1.);
  Matrix * XtXs = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  XtXs->multiply(temp, 'N', mOnes, 'N'); //1 column matrix, each row should contain the XtX for each SNP.
  delete mOnes;
  delete temp;
  if( communicator->group == 0 )
  {
    misc.message << "Computation finished after " << misc.setGetElapsedTime("XtXs", true) << "." << std::endl;
  }
  
  if( communicator->group == 0 )
  {
    misc.setGetElapsedTime("Xtys");
    misc.message << "Computing Xty matrices..." << std::endl;
  }
  Matrix * Xtys = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Xtys->multiply(genotype, 'N', phenotypes, 'N'); //MxN matrix where column n contain the Xty for the phenotype n and the row m the snp m.
  if( communicator->group == 0 )
  {
    misc.message << "Computation finished after " << misc.setGetElapsedTime("Xtys", true) << "." << std::endl;
  }
  
  std::vector< std::vector<double> > gXtys;
  Xtys->matrixToStandardVector(gXtys);
  std::vector<double> gytys;
  ytys->matrixToStandardVector(gytys);
  std::vector<double> gXtXs;
  XtXs->matrixToStandardVector(gXtXs);
  
  if( misc.gt(gXtXs.size() != genotype->nGlobRows || gytys.size() != phenotypes->nGlobCols || gXtys.size() != genotype->nGlobRows || gXtys.size() == 0 ) )
  {
    misc.error("Error: An internal error was happened. Unexpected matrix dimensions when performing GWAS on multiple phenotypes.", 0);
  }
  if( misc.gt(gXtys[0].size() != phenotypes->nGlobCols) )
  {
    misc.error("Error: An internal error was happened. Unexpected matrix dimensions when performing GWAS on multiple phenotypes.", 0);
  }
  
  
  if( SNPidxs[ communicator->group ].size() != genotype->nGlobRows )
  {
    misc.error("Error: An internal error was happened. When performing individual SNP test in GWAS, discordance in grouped communicator size with local genotype size.", 0);
  }
  
  if( communicator->group == 0 )
  {
    misc.setGetElapsedTime("mpGWAS");
    misc.message << "Starting GWAS for " + i2s(Xtys->nGlobCols) + " phenotypes and " + i2s(nSNPs) + " SNPs..." << std::endl;
  }
  
  int stepsForPercentageOutput = (genotype->nGlobRows/10);
  stepsForPercentageOutput = (stepsForPercentageOutput == 0?1:stepsForPercentageOutput);
  
  for(int isnp = 0; isnp < genotype->nGlobRows; isnp++ )
  {
    if( communicator->group == 0 && isnp % stepsForPercentageOutput == 1 )
    {
      misc.message << (100*isnp)/genotype->nGlobRows << "% completed." << std::endl;
    }
    
    int nIndividuals = genotype->nGlobCols;
    
    GLMResults SNPResults;
    bool success = computeGLMWithoutCovarianceMultiplePhenos(gXtXs[isnp], gXtys[isnp], gytys, nIndividuals, SNPResults);

    if(communicator->mpiRoot)
    {
      int globalISNPdx = SNPidxs[ communicator->group ][isnp];
      groupedResults[ globalISNPdx ] = SNPResults;
    }
  }
  
  delete genotype;
  delete phenotypes;
  delete ytys;
  delete XtXs;
  
  if(communicator->mpiRoot == false)
  {
    groupedResults.clear();
  }
  
  communicator = globalCommunicator;
  delete groupedCommunicator;
  
  misc.message << "100%" << std::endl;
  misc.message << "GWAS finished after " << misc.setGetElapsedTime("mpGWAS", true) << std::endl;
  gatherResults(results, groupedResults, resultSNPInfo, gcgenotype, unfittedSNPs);
  
  storeResultsMultiplePhenotype(results, resultSNPInfo, gcphenotypes->colLabels);
  
  //delete intermediate matrices
  delete gcytys;
  
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

void GWAS::computeIndividualGWASMultiplePhenotypes(Genotype * gcgenotype, LabeledMatrix * gcphenotypes, Communicator * globalCommunicator)
{
  std::map<std::string, GLMResults > results;
  std::map<std::string, std::vector<SNP> > resultSNPInfo;
  std::vector<std::string> unfittedSNPs;
  
  int nSNPs = gcgenotype->nSNPs;

  //Compute some matrices which whill be conserved.
  misc.setGetElapsedTime("ytys");
  misc.message << "Computing yty elements..." << std::endl;
  Matrix * temp = new Matrix(gcphenotypes->matrix);
  temp->elementWiseMultiplication(gcphenotypes->matrix);
  Matrix * mOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, temp->nGlobRows, 1);
  mOnes->fillWithConstant(1.);
  Matrix * gcytys = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  gcytys->multiply(temp, 'T', mOnes, 'N'); //1 column matrix, each row should contain the yty for each phenotype.
  delete temp;
  delete mOnes;
  misc.message << "Computation finished after " << misc.setGetElapsedTime("ytys", true) << "." << std::endl;
  
  
  //Redistribute data between different communicators.
  Matrix * genotype = gcgenotype->genotypes;
  Matrix * phenotypes = gcphenotypes->matrix;
  Matrix * ytys = gcytys;
  
  Matrix * globalCovarianceMatrix = NULL;
  if( this->useCovariateMatrix == true )
  {
    misc.error("Error: This analysis is not implemented.", 0);
  }
  
  //Compute matrices needed after:
  if( communicator->group == 0 )
  {
    misc.setGetElapsedTime("XtXs");
    misc.message << "Computing XtX elements..." << std::endl;
  }
  temp = new Matrix(genotype);
  temp->elementWiseMultiplication(genotype);
  mOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, temp->nGlobCols, 1);
  mOnes->fillWithConstant(1.);
  Matrix * XtXs = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  XtXs->multiply(temp, 'N', mOnes, 'N'); //1 column matrix, each row should contain the XtX for each SNP.
  delete mOnes;
  delete temp;
  if( communicator->group == 0 )
  {
    misc.message << "Computation finished after " << misc.setGetElapsedTime("XtXs", true) << "." << std::endl;
  }
  
  if( communicator->group == 0 )
  {
    misc.setGetElapsedTime("Xtys");
    misc.message << "Computing Xty matrices..." << std::endl;
  }
  Matrix * Xtys = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Xtys->multiply(genotype, 'N', phenotypes, 'N'); //MxN matrix where column n contain the Xty for the phenotype n and the row m the snp m.
  if( communicator->group == 0 )
  {
    misc.message << "Computation finished after " << misc.setGetElapsedTime("Xtys", true) << "." << std::endl;
  }
  
  std::vector< std::vector<double> > gXtys;
  Xtys->matrixToStandardVector(gXtys);
  std::vector<double> gytys;
  ytys->matrixToStandardVector(gytys);
  std::vector<double> gXtXs;
  XtXs->matrixToStandardVector(gXtXs);
  
  if( misc.gt(gXtXs.size() != genotype->nGlobRows || gytys.size() != phenotypes->nGlobCols || gXtys.size() != genotype->nGlobRows || gXtys.size() == 0 ) )
  {
    misc.error("Error: An internal error was happened. Unexpected matrix dimensions when performing GWAS on multiple phenotypes.", 0);
  }
  if( misc.gt(gXtys[0].size() != phenotypes->nGlobCols) )
  {
    misc.error("Error: An internal error was happened. Unexpected matrix dimensions when performing GWAS on multiple phenotypes.", 0);
  }
  
  if( communicator->group == 0 )
  {
    misc.setGetElapsedTime("mpGWAS");
    misc.message << "Starting GWAS for " + i2s(Xtys->nGlobCols) + " phenotypes and " + i2s(nSNPs) + " SNPs..." << std::endl;
  }
  
  int stepsForPercentageOutput = (genotype->nGlobRows/10);
  stepsForPercentageOutput = (stepsForPercentageOutput == 0?1:stepsForPercentageOutput);
  
  for(int isnp = 0; isnp < genotype->nGlobRows; isnp++ )
  {
    if( communicator->group == 0 && isnp % stepsForPercentageOutput == 1 )
    {
      misc.message << (100*isnp)/genotype->nGlobRows << "% completed." << std::endl;
    }
    
    int nIndividuals = genotype->nGlobCols;
    
    GLMResults SNPResults;
    bool success = computeGLMWithoutCovarianceMultiplePhenos(gXtXs[isnp], gXtys[isnp], gytys, nIndividuals, SNPResults);

    if(communicator->mpiRoot)
    {
      results[ gcgenotype->SNPIds[isnp] ] = SNPResults;
    }
  }
  
  genotype = NULL;
  phenotypes = NULL;
  delete ytys;
  delete XtXs;
  
  if( communicator->group == 0 )
  {
    misc.message << "100%" << std::endl;
    misc.message << "GWAS finished after " << misc.setGetElapsedTime("mpGWAS", true) << std::endl;
  }
  
  if(communicator->mpiRoot)
  {
    for(int isnp = 0; isnp<gcgenotype->nSNPs; isnp++)
    {
      std::string snpname = gcgenotype->SNPs[isnp].name;
      if( results[ snpname ].success == true )
      {
        std::vector<SNP> temp;
        temp.push_back(gcgenotype->SNPs[isnp]);
        resultSNPInfo[ snpname ] = temp;
      }
      else
      {
        results.erase(snpname);
        unfittedSNPs.push_back(snpname);
      }
    }
  }
  
  //Save results
  std::stringstream partialResultsSS;
  storeResultsMultiplePhenotype(results, resultSNPInfo, gcphenotypes->colLabels, partialResultsSS, communicator->group == 0);
  
  Communicator * tempCommunicator = communicator;
  communicator = globalCommunicator;
#if defined(BOOSTLIB) && defined(ZLIB)
  std::string partialResults = compressData(partialResultsSS.str());
  communicator->storeArraysMPI(options.outFile + ".multipheno.gwas.snps.gz", partialResults);
#else
  std::string partialResults = partialResultsSS.str();
  communicator->storeArraysMPI(options.outFile + ".multipheno.gwas.snps", partialResults);
#endif
  communicator = tempCommunicator;
  
  //Save unfitted SNPs
  std::stringstream ssUnfittedSNPs;
  if(communicator->mpiRoot && unfittedSNPs.size() != 0)
  {
    misc.message << "Warning: There are " << unfittedSNPs.size() << " SNPs which cannot be fitted. They are stored in file [ " << (options.outFile + ".gwas.unfitted") << " ]." << std::endl;
    //Message message(options.outFile + ".gwas.unfitted");
    for(int i = 0; i < unfittedSNPs.size(); i++)
    {
      ssUnfittedSNPs << unfittedSNPs[i] << std::endl;
    }
  }
  tempCommunicator = communicator;
  communicator = globalCommunicator;
  communicator->storeArraysMPI(options.outFile + ".gwas.unfitted", ssUnfittedSNPs.str());
  communicator = tempCommunicator;
}


bool GWAS::computeGLMWithoutCovarianceMultiplePhenos(double sXtX,  std::vector<double> & phenosXtys, std::vector<double> & phenoytys, int nIndividuals, GLMResults & results)
{
  if(sXtX <= 0.)
  {
    results.success = false;
    return false;
  }
  double sXtX_i = 1./sXtX;
  
  if(communicator->mpiRoot)
  {
    double n_q_1 = double(double(nIndividuals) - 1.);
    
    results.btXty = -1;
    
    results.SSE = -1;
    results.MSE = -1;
    
    results.b = std::vector<double>(phenosXtys.size());
    results.SE = std::vector<double>(phenosXtys.size());
    results.tStatistics = std::vector<double>(phenosXtys.size());
    results.tStatisticPValues = std::vector<double>(phenosXtys.size());
    
    for(int ipheno = 0; ipheno < phenosXtys.size(); ipheno++)
    {
      results.b[ipheno] = sXtX_i*phenosXtys[ipheno];
      
      double btXty = results.b[ipheno]*phenosXtys[ipheno];
      
      double SSE = phenoytys[ipheno] - btXty;
      double MSE = SSE/double(n_q_1);
      
      double temp = sqrt(MSE*sXtX_i);
      results.SE[ipheno] = temp;
      results.tStatistics[ipheno] = results.b[ipheno]/temp;
      results.tStatisticPValues[ipheno] = 2*tStatCDF(n_q_1, fabs(results.b[ipheno]/temp) ) ;
    }
  }
  results.type = OLSModelType;
  results.success = true;
  
  return true;
}

void GWAS::storeResultsMultiplePhenotype(std::map<std::string, GLMResults > & effects, std::map<std::string, std::vector<SNP> > & effectsSNPs, std::vector<std::string> & phenoLabels)
{
  if(communicator->mpiRoot)
  {
    Message messageSNPs(options.outFile + ".multipheno.gwas.snps");
    storeResultsMultiplePhenotype(effects, effectsSNPs, phenoLabels, *messageSNPs.output, true);
  }
}

void GWAS::storeResultsMultiplePhenotype(std::map<std::string, GLMResults > & effects, std::map<std::string, std::vector<SNP> > & effectsSNPs, std::vector<std::string> & phenoLabels, std::ostream & messageSNPs, bool outputHeader)
{
  if(communicator->mpiRoot)
  {
    //Message messageSNPs(options.outFile + ".multipheno.gwas.snps");
    if(outputHeader == true)
    {
      messageSNPs << "SNP ALLELE MEAN STDEV";
      for(int ipheno = 0; ipheno < phenoLabels.size(); ipheno++)
      {
        std::string phenoLabel = phenoLabels[ipheno];
        messageSNPs << /*" BETA-" + phenoLabel <<*/ " NBETA-" + phenoLabel << /*" SE-" + phenoLabel <<*/ " NSE-" + phenoLabel << " PV-" + phenoLabel;
      }
      messageSNPs << std::endl;
    }
    
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
      
      if( globalEffects.size() != phenoLabels.size() || globalPValues.size() != phenoLabels.size() || groupResults.SE.size() != phenoLabels.size())
      {
        misc.error("Error: An internal error was happened. Unexpected number of SNP results.", 0);
      }
      
      if(effectsSNPs[group].size() != 1)
      {
        misc.error("Error: An internal error was happened. Unexpected number of SNPs.", 0);
      }
      SNP snp = effectsSNPs[group][0];
      
      messageSNPs << snp.name;
      messageSNPs << " " << snp.allele2;
      messageSNPs << " " << std::setprecision(3) << 2.*snp.p2;
      messageSNPs << " " << std::setprecision(3) << snp.standardDev;
      for(int ipheno = 0; ipheno < phenoLabels.size(); ipheno++)
      {
        //messageSNPs << " " << globalEffects[ ipheno ];
        messageSNPs << " " << std::setprecision(5) << globalEffects[ ipheno ]/snp.standardDev;
        //messageSNPs << " " << groupResults.SE[ ipheno ];
        messageSNPs << " " << std::setprecision(5) << groupResults.SE[ ipheno ]/snp.standardDev;
        messageSNPs << " " << globalPValues[ ipheno ];
      }
      messageSNPs << std::endl;
    }
  }
}
