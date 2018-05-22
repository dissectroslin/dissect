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

#include "matrix.h"
#include "options.h"
#include "communicator.h"
#include "misc.h"
#include "global.h"
#include "reml.h"
#include "covariate.h"
#include "genotype.h"
#include "phenotype.h"
#include "covariancematrix.h"
#include "auxiliar.h"
#include "message.h"

#include <fstream>
#include <iomanip>
#include <cmath>

REML::REML(bool argWriteResults)
{
  this->useMLinsteadOfREML = options.useMLinsteadOfREML;
  
  this->y = NULL;
  this->X = NULL;
  
  this->V = NULL;
  
  this->ViX = NULL;
  this->XtViX_i = NULL;
  this->P = NULL;
  this->PDiagonal = NULL;
  this->Py = NULL;
  this->subVPy = NULL;
  this->AI = NULL;
  this->yPsubVPy_trPsubV = NULL;
  this->yPsubVPy_trVisubV = NULL;
  this->VisubVTraces = NULL;
  
  this->secondDerivativesMatrixType = REMLParameters::undefined;
  
  this->success = false;
  this->logLikelihoodConverged = false;
  this->variancesConverged = false;
  this->gradientConverged = false;
  
  this->writeResults = argWriteResults;
  this->writeBLUE = options.writeBLUEInReducedModels;
  this->writeSummary = false;
  
  this->warnings.clear();
  
  this->currentStep = 0;
  
  this->singlePrecisionInversion = options.allowSinglePrecisionInversion;
  
  this->nPhenotypes = 0;
  
  this->usingDiagonalKernels = false;
  
  this->sREMLML = "REML";
}

REML::~REML()
{
  deleteyList();
  if(this->y != NULL)
  {
    delete this->y;
  }
  if(this->X != NULL)
  {
    delete this->X;
  }
  if(this->V != NULL)
  {
    delete this->V;
  }
  for(std::map<std::string, Genotype *>::iterator it = this->SNPsBLUPGenotypes.begin(); it != this->SNPsBLUPGenotypes.end(); ++it)
  {
    delete it->second;
  }
  this->SNPsBLUPGenotypes.clear();
  
  Matrix * deleted = NULL;
  for(std::map<std::string, Matrix *>::iterator it = this->GRMEigenVectors.begin(); it != this->GRMEigenVectors.end(); ++it)
  {
    if(deleted == NULL)
    {
      deleted = it->second;
      delete it->second;
    }
    else
    {
      if(deleted != it->second)
      {
        misc.error("Error: An internal error has happened. During a check of the eigenvectors, detected an unexpected eigenvector. This is an unexpected behaviour that could lead to unexpected results. Please, contact us, and do not use any of the outputs from this analysis.", 0);
      }
    }
  }
  deleted = NULL;
  this->GRMEigenVectors.clear();
  
  deleteIntermediateMatrices();
}

void REML::deleteIntermediateMatrices()
{
  if( this->ViX != NULL )
  {
    delete this->ViX;
    this->ViX = NULL;
  }
  if( this->XtViX_i != NULL )
  {
    delete this->XtViX_i;
    this->XtViX_i = NULL;
  }
  if( this->P != NULL )
  {
    delete this->P;
    this->P = NULL;
  }
  if( this->PDiagonal != NULL )
  {
    delete this->PDiagonal;
    this->PDiagonal = NULL;
  }
  if( this->Py != NULL )
  {
    delete this->Py;
    this->Py = NULL;
  }
  if( this->subVPy != NULL )
  {
    delete this->subVPy;
    this->subVPy = NULL;
  }
  if( this->AI != NULL )
  {
    delete this->AI;
    this->AI = NULL;
  }
  if( this->yPsubVPy_trPsubV != NULL )
  {
    delete this->yPsubVPy_trPsubV;
    this->yPsubVPy_trPsubV = NULL;
  }
  if( this->yPsubVPy_trVisubV != NULL )
  {
    delete this->yPsubVPy_trVisubV;
    this->yPsubVPy_trVisubV = NULL;
  }
  if( this->VisubVTraces != NULL )
  {
    delete this->VisubVTraces;
    this->VisubVTraces = NULL;
  }
}

bool REML::prepare(REMLType type, std::vector<Kernel*> & kernels, std::vector<double> weights, std::vector<int> phenotypeColumns, std::vector<double> heritabilities, std::vector<std::pair<std::string, std::string> > covariateFiles)
{
  //Perform some checks and parameter initialization
  
  if(kernels.size() < 1)
  {
    misc.error("Error: An internal error was happened. No GRMs specified for REML analysis.", 0);
  }
  
  if(phenotypeColumns.size() < 1)
  {
    misc.error("Error: An internal error was happened. No Phenotypes specified for REML analysis.", 0);
  }
  
  if(covariateFiles.size() != phenotypeColumns.size())
  {
    misc.error("Error: An internal error was happened. Different number of covariate files than phenotypes columns to analyze specified for REML analysis.", 0);
  }  

  if(weights.size() == 0)
  {
    double equallyDistributedWeight = 1./double(kernels.size());
    for(int i = 0; i<kernels.size(); i++)
    {
      weights.push_back(equallyDistributedWeight);
    }
  }
  else
  {
    if( kernels.size() != weights.size() )
    {
      misc.error("Error: An internal error was happened. Kernel weights not properly defined.", 0);
    }
  }
  
  if(heritabilities.size() == 0)
  {
    for(int i = 0; i<phenotypeColumns.size(); i++)
    {
      heritabilities.push_back(0.5);
    }
  }
  else
  {
    if( heritabilities.size() != phenotypeColumns.size() )
    {
      misc.error("Error: An internal error was happened. Initial heritabilities not properly defined.", 0);
    }
  }
  
  //Set reml type
  
  this->type = type;
  
  //Set the number of phenotypes
  
  this->nPhenotypes = phenotypeColumns.size();
  
  //Sanitize kernels
  
  bool sanitized = true;
  for(int i = 0; i<kernels.size(); i++)
  {
    bool temp = kernels[i]->sanitizeKernel();
    if( temp == false )
    {
      sanitized = false;
    }
  }
  if( sanitized == false )
  {
    for(int i = 0; i<kernels.size(); i++)
    {
      delete kernels[i];
      kernels[i] = NULL;
    }
    kernels.clear();
    return false;
  }
  
  //Load the data for the phenotype (y).
  
  std::vector<Phenotype*> phenotypes;
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    Phenotype * temp = new Phenotype(cyclicDistribution, options.phenotypesFile, phenotypeColumns[i]);
    phenotypes.push_back(temp);
  }
  filterPhenotypesUsingCouples(phenotypes);
  
  //Load the data for the covariance (X).
  
  std::vector< std::set<std::string> > tempSharedIndividualsSets;
  if( communicator->mpiRoot == true ) //Precompute shared individuals with kernels and phenotypes to avoid X matrices with dependent columns after filtering.
  {
    std::vector<std::string> tempSharedKernelIndividualsVector;
    tempSharedKernelIndividualsVector = kernels[0]->individualIds;
    for(int i = 1; i<kernels.size(); i++)
    {
      tempSharedKernelIndividualsVector = intersectionStringVectors(2, &tempSharedKernelIndividualsVector, &(kernels[i]->individualIds));
    }
    for(int ip = 0; ip<this->nPhenotypes; ip++)
    {
      std::vector<std::string> temp = intersectionStringVectors(2, &tempSharedKernelIndividualsVector, &(phenotypes[ ip ]->individualIds));
      if(temp.size() == 0)
      {
        misc.error("Error: An error has happened. The overlap between the covariance matrices and the phenotypes is empty for the phenotype column " + i2s(phenotypeColumns[ip]) + ".", 0);
      }
      tempSharedIndividualsSets.push_back( std::set<std::string>(temp.begin(), temp.end()) );
    }
  }
  else
  {
    for(int ip = 0; ip<this->nPhenotypes; ip++)
    {
      tempSharedIndividualsSets.push_back( std::set<std::string>() );
    }
  }
  
  this->covariateNames.clear();
  
  std::vector<Covariate*> covariates;
  int previousIndex = 0;
  for(int i = 0; i<covariateFiles.size(); i++)
  {
    if(options.joinCovariatesVertically == true) //This option is untested
    {
      Covariate * tempCovariate = new Covariate(covariateFiles[i].first, covariateFiles[i].second, phenotypes[i]->individualIds, false, tempSharedIndividualsSets[ i ]);
      for(int j = 0; j<covariates.size(); j++)
      {
        covariates[j]->syncronizeDiscreteCovariateCategoriesWith(tempCovariate);
      }
      covariates.push_back(tempCovariate);
    }
    else
    {
      Covariate * tempCovariate = new Covariate(covariateFiles[i].first, covariateFiles[i].second, phenotypes[i]->individualIds, true, tempSharedIndividualsSets[ i ]);
      covariates.push_back(tempCovariate);
      this->addCovariatesNames(tempCovariate, addSuffix("p" + i2s(i + 1) + "-"), previousIndex);
      previousIndex += tempCovariate->covariates->nGlobCols;
    }
  }
  
  if(options.joinCovariatesVertically == true) //This option is untested
  {
    for(int i = 0; i<covariates.size(); i++)
    {
      covariates[i]->parseRawCovariates(phenotypes[i]->individualIds, covariates.size(), i);
    }
    this->addCovariatesNames(covariates[0]);
  }
  
  //Load the data for the environmental weights.
  
  std::vector<Phenotype *> environmentalWeights;
  if( options.environmentalWeightsFile != "" )
  {
    for(int i = 0; i<this->nPhenotypes; i++) //This can be optimized for not reading several times the same. Although currently this option is only abailable for single phenotypes.
    {
      Phenotype * temp = new Phenotype(cyclicDistribution, options.environmentalWeightsFile, options.environmentalWeightsCol);
      environmentalWeights.push_back(temp);
    }
  }
  
  //Search for common individuals
  
  std::vector< std::vector<std::string> > commonIndividualsInGRMOrder(this->nPhenotypes, std::vector<std::string>() );
  std::vector< std::vector<int> > idxsKeptKernel0(this->nPhenotypes, std::vector<int>() );
  int totalIndividuals;
  if(communicator->mpiRoot)
  {
    for(int i = 0; i<this->nPhenotypes; i++)
    {
      std::vector<std::string> tempIndividuals = intersectionStringVectors(3, &(kernels[0]->individualIds), &phenotypes[i]->individualIds, &covariates[i]->individualIds);
      if( environmentalWeights.size() != 0 )
      {
        tempIndividuals = intersectionStringVectors(2, &tempIndividuals, &(environmentalWeights[i]->individualIds));
      }
      for(int j = 1; j<kernels.size(); j++)
      {
        tempIndividuals = intersectionStringVectors(2, &tempIndividuals, &(kernels[j]->individualIds));
      }
      
      std::vector<std::string> tempOrderedIndividuals = orderVectorAsTemplate(kernels[0]->individualIds, tempIndividuals);
      commonIndividualsInGRMOrder[i] = tempOrderedIndividuals;
      
      if( tempOrderedIndividuals.size() < 1)
      {
        misc.error("Error: There are not enough individuals for one of the traits for performing the analysis.", 0);
      }
      misc.message << tempOrderedIndividuals.size() << " individuals are kept for the analysis for phenotype " << i + 1 << "." << std::endl;
      
      std::vector<int> tempIdxsKept = extractMapValues(tempOrderedIndividuals, kernels[0]->individualIdsIdx);
      idxsKeptKernel0[i] = tempIdxsKept;
    }
  }

  int nIndividualsBeforeFilterKernel0 = kernels[0]->nIndividuals;
  std::vector<int> nIndividualsTraits;
  int nTotalIndividuals = 0;
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    int nIndividualsTrait = commonIndividualsInGRMOrder[i].size();
    communicator->broadcast(&nIndividualsTrait, 1);
    nIndividualsTraits.push_back(nIndividualsTrait);
    
    nTotalIndividuals += nIndividualsTrait;
  }
  
  ////////////////////////////////////////////////////
  // Init phenotype matrix
  
  std::vector<Matrix*> tempyMatrices;
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    phenotypes[i]->filterIndividuals(commonIndividualsInGRMOrder[i]);
    tempyMatrices.push_back(new Matrix(phenotypes[i]->phenotypes));
    if( i == 0)
    {
      this->y = new Matrix(phenotypes[i]->phenotypes);
    }
    else
    {
      Matrix * temp = new Matrix(this->y);
      this->y->joinMatricesVertically(temp, phenotypes[i]->phenotypes);
      delete temp;
    }
    delete phenotypes[i];
  }
  phenotypes.clear();
  
  ////////////////////////////////////////////////////
  // Filter environmental weights (if any) and store the matrix as a global vector
  std::vector<double*> globalEnvironmentalWeights;
  if( environmentalWeights.size() != 0 )
  {
    for(int i = 0; i<this->nPhenotypes; i++)
    {
      environmentalWeights[i]->filterIndividuals(commonIndividualsInGRMOrder[i]);
      int nWeights = environmentalWeights[i]->phenotypes->nGlobRows;
      if( options.scaleEnvironmentalWeightTrace == true )
      {
        Matrix * rowOfOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, 1, nWeights);
        rowOfOnes->fillWithConstant(1.);
        Matrix * tempTotal = new Matrix();
        tempTotal->multiply(rowOfOnes, 'N', environmentalWeights[i]->phenotypes, 'N');
        double elementsSum = 0.;
        tempTotal->gatherMatrix(&elementsSum);
        communicator->broadcast(&elementsSum);
        delete rowOfOnes;
        delete tempTotal;
        environmentalWeights[i]->phenotypes->scaleBy(double(nWeights)/elementsSum);
      }
      if(communicator->mpiRoot == true)
      {
        double * gWeights = new double [nWeights];
        environmentalWeights[i]->phenotypes->gatherMatrix(gWeights);
        globalEnvironmentalWeights.push_back(gWeights);
      }
      else
      {
        globalEnvironmentalWeights.push_back(NULL);
      }
      delete environmentalWeights[i];
    }
  }
  environmentalWeights.clear();
  
  ////////////////////////////////////////////////////
  // Check whether Kernel diagonalization is coherent. If all kernels are diagonal, they must have the same individuals and this must be shared with phenotypes and covariates.
  // If only some kernels are diagonal it is an error if options.forceUseDiagonalizedKernels == false (default option). Kernels are de-diagonalized otherwise.
  this->usingDiagonalKernels = false;
  bool allKernelsDiagonal = true;
  bool someKernelsDiagonal = false;
  for(int i = 0; i<kernels.size(); i++)
  {
    if(kernels[i]->diagonalized == true)
    {
      someKernelsDiagonal = true;
    }
    else
    {
      allKernelsDiagonal = false;
    }
  }
  
  if( allKernelsDiagonal == false && someKernelsDiagonal == true )
  {
    if(options.forceUseDiagonalizedKernels == true)
    {
      for(int i = 0; i<kernels.size(); i++)
      {
        kernels[i]->recoverKernelFromEigenDecomposition();
      }
    }
    else
    {
      misc.error("Error: Sorry, this analysis cannot be performed with diagonal GRMs/Kernels. You can use the option --force-use-diag-kernels to force converting diagonal GRMs/Kernels to their non-diagonalized form before starting the analysis. However, this have some limitations: any analysis that depends on the normalization matrix cannot be performed from diagonal Kernels/GRMs (e.g. regional analysis)", 0);
    }
  }
  else if( allKernelsDiagonal == true )
  {
    for(int i = 0; i<kernels.size(); i++)
    {
      if(kernels[i]->individualIds != commonIndividualsInGRMOrder[i] || kernels[i]->individualIds != kernels[0]->individualIds)
      {
        misc.error("Error: Sorry, this analysis cannot be performed with diagonal GRMs/Kernels. All kernels have to share same individuals and these also have to be present on phenotype and covariate data.", 0);
      }
    }
    
    if(kernels.size() == 1 && this->nPhenotypes == 1)
    {
      this->GRMEigenVectors[kernels[0]->name] = new Matrix(kernels[0]->eigenVectors);
      this->usingDiagonalKernels = true;
      
      transformUsingEigenVectors(kernels[0]->name, 'T', &(this->y), this->nPhenotypes);
      
      if(this->SNPsBLUPGenotypes.count(kernels[0]->name) != 0)
      {
        if(this->SNPsBLUPGenotypes[ kernels[0]->name ]->spaceChangedKernelName == "")
        {
          this->SNPsBLUPGenotypes[ kernels[0]->name ]->normalizeGenotypes();
          this->SNPsBLUPGenotypes[ kernels[0]->name ]->filterSNPsAndIndividuals(this->SNPsBLUPGenotypes[ kernels[0]->name ]->SNPIds, kernels[0]->individualIds);
          if( this->SNPsBLUPGenotypes[ kernels[0]->name ]->individualIds != kernels[0]->individualIds )
          {
            misc.error("Error: The order of individuals in the Kernel/GRM is different that the order in the genotypes file. Sorry, at this version DISSECT needs individuals in both files must be in the same order in both files.", 0);
          }
          Matrix *temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
          temp->multiply(this->SNPsBLUPGenotypes[ kernels[0]->name ]->genotypes, 'N', kernels[0]->eigenVectors, 'N');
          delete this->SNPsBLUPGenotypes[ kernels[0]->name ]->genotypes;
          this->SNPsBLUPGenotypes[ kernels[0]->name ]->genotypes = temp;
          
          this->SNPsBLUPGenotypes[ kernels[0]->name ]->changeSpace(kernels[0]->name);
        }
        else
        {
          if(this->SNPsBLUPGenotypes[ kernels[0]->name ]->spaceChangedKernelName != kernels[0]->name)
          {
            misc.error("Error: An internal error has happened. Trying to change the space of a genotypes matrix which has already changed using a different kernel.", 0);
          }
          if( this->SNPsBLUPGenotypes[ kernels[0]->name ]->individualIds != kernels[0]->individualIds )
          {
            misc.error("Error: An internal error has happened. Unexpected individuals order between genotypes and diagonalized kernels.", 0);
          }
        }
      }
    }
    else if( kernels.size() == 1  && this->nPhenotypes > 1 )
    {
      this->GRMEigenVectors[kernels[0]->name] = new Matrix(kernels[0]->eigenVectors);
      this->usingDiagonalKernels = true;
      
      transformUsingEigenVectors(kernels[0]->name, 'T', &(this->y), this->nPhenotypes);
      
      if(this->SNPsBLUPGenotypes.count(kernels[0]->name) != 0)
      {
        misc.error("Error: An internal error was happened. Unexpected genotype data to be corrected when preparing REML for performing the analysis with diagonal Kernels.", 0);
      }
    }
    else
    {
      misc.error("Error: Sorry, this type of analysis is not implemented, yet.", 0);
    }
    
    this->singlePrecisionInversion = false;
  }
  
  ////////////////////////////////////////////////////
  // Init covariate matrix and compute initial variances.
  std::vector<double> phenotypeVariances;
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    covariates[i]->filterIndividuals(commonIndividualsInGRMOrder[i]);

    if( kernels.size() == 1 && this->GRMEigenVectors.count(kernels[0]->name) != 0 && this->usingDiagonalKernels == true ) //Correct in case we are using diagonal kernels.
    {
      transformUsingEigenVectors(kernels[0]->name, 'T', &(covariates[i]->covariates), 1);
    }
    
    double phenotypeVariance = computeInitialVariance(tempyMatrices[i], covariates[i]->covariates);
    phenotypeVariances.push_back(phenotypeVariance);
    delete tempyMatrices[i];
    
    if( i == 0 )
    {
      this->X = new Matrix(covariates[i]->covariates);
    }
    else
    {
      Matrix * temp = new Matrix(this->X);
      if(options.joinCovariatesVertically == true) //Untested option
      {
        if(temp->nGlobCols != covariates[i]->covariates->nGlobCols)
        {
          misc.error("Error: The covariate or quantitative covariate files of at least two traits do not have the same number of columns.", 0); //.
        }
        this->X->joinMatricesVertically(temp, covariates[i]->covariates);
      }
      else
      {
        subMatrix sm1 = subMatrix(0, 0, temp->nGlobRows, temp->nGlobCols);
        subMatrix sm2 = subMatrix(temp->nGlobRows, temp->nGlobCols, covariates[i]->covariates->nGlobRows, covariates[i]->covariates->nGlobCols);
        this->X->joinMatrices(temp, sm1, covariates[i]->covariates, sm2, 0.);
      }
      delete temp;
    }
    delete covariates[i];
  }
  covariates.clear();
  tempyMatrices.clear();
  

  ////////////////////////////////////////////////////
  // Create the covariance matrix
  communicator->broadcast(&nTotalIndividuals, 1);
  
  this->dimension = nTotalIndividuals;
  
  this->V = new CovarianceMatrix(this->dimension, allKernelsDiagonal, this->nPhenotypes);
  
  ////////////////////////////////////////////////////
  // Insert the covariance matrices
  std::vector<std::string> kernelNames; //For storing kernels names before deleting them.
  std::vector<KernelType> kernelTypes;
  std::map<std::string, double> kernelAverages;
  this->vIndividuals.clear();
  for(int i = 0; i<kernels.size(); i++)
  {
    for(int j = 0; j<this->nPhenotypes; j++)
    {
      kernels[i]->normalize();
      Kernel * filteredKernel = new Kernel(kernels[i]); //Kernel for removing trait individuals
      filteredKernel->filterIndividuals(commonIndividualsInGRMOrder[j], false);
      
      if(i == 0)
      {
        this->vIndividuals.push_back( filteredKernel->individuals );
      }
      
      this->V->insertCovarianceMatrix(kernels[i]->name + "_" + i2s(j+1), filteredKernel);
      
      for(int k = j + 1; k<this->nPhenotypes; k++)
      {
        Kernel * kernelCov = new Kernel(kernels[i]); //Kernel for intersection trait1-trait2
        kernelCov->filterIndividualsAsymmetric(commonIndividualsInGRMOrder[j], commonIndividualsInGRMOrder[k], false);
        
        this->V->insertCovarianceMatrix(kernels[i]->name + "_" + i2s(j+1) + "_" + i2s(k+1), kernelCov);
      }
    }
    
    kernelNames.push_back(kernels[i]->name);
    kernelTypes.push_back(kernels[i]->type);
    if(kernels[i]->type == kernelSquaredExponential)
    {
      kernelAverages[kernels[i]->name] = options.expKernelParameterInitialFactor/kernels[i]->getNormalizedKernel()->elementsAverage();
    }
    delete kernels[i];
    kernels[i] = NULL;
  }
  
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    Matrix * temp;
    if( allKernelsDiagonal == false )
    {
      temp = new Matrix(cyclicDistribution, nIndividualsTraits[i], nIndividualsTraits[i], communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    }
    else
    {
      temp = new Matrix(diagonalDistribution, nIndividualsTraits[i], nIndividualsTraits[i], communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    }
    
    if( globalEnvironmentalWeights.size() == 0 )
    {
      temp->fillDiagonal(1.);
    }
    else
    {
      temp->fillWithConstant(0.);
      temp->setDiagonal(globalEnvironmentalWeights[i], temp->nGlobCols);
    }
    
    this->V->insertCovarianceMatrix("E_" + i2s(i+1), temp);
  }
  
  std::map<std::pair<int, int>, bool> computeEnvironmentalCovariances;
  if(options.environmentalCovariance)
  {
    for(int i = 0; i<this->nPhenotypes; i++)
    {
      for(int j = i + 1; j<this->nPhenotypes; j++)
      {
        if(globalEnvironmentalWeights.size() != 0)
        {
          misc.error("Error: An internal error was happened. Environmental weights are not currently allowed on multireml analysis. If you are interested on this, please, contact us.", 0);
        }
        
        int commonIndividualsBetweenTraits = ( intersectionStringVectors(2, &(commonIndividualsInGRMOrder[i]), &(commonIndividualsInGRMOrder[j])) ).size();
        int tempBiggestSample = (commonIndividualsInGRMOrder[i].size()>commonIndividualsInGRMOrder[j].size()?commonIndividualsInGRMOrder[i].size():commonIndividualsInGRMOrder[j].size());
        double factor = double(commonIndividualsBetweenTraits)/double(tempBiggestSample);

        if( misc.gt(factor < 0.1) )
        {
          misc.message << "Less than 10% of individuals were measured for both traits. The residual covariance component for traits " << i + 1 << " and " << j + 1 << " is discarded." << std::endl;
          computeEnvironmentalCovariances[std::pair<int, int>(i, j)] = false;
          continue;
        }
        
        Matrix * covar;
        if( allKernelsDiagonal == false )
        {
          Matrix * temp = new Matrix(cyclicDistribution, nIndividualsBeforeFilterKernel0, nIndividualsBeforeFilterKernel0, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
          temp->fillDiagonal(1.);
          covar = new Matrix(cyclicDistribution);
          temp->filterRowsAndColumns(covar, idxsKeptKernel0[i], idxsKeptKernel0[j]);
          delete temp;
          
          covar->symmetric = false;
          covar->uplo = 'B';
        }
        else
        {
          if( nIndividualsTraits[i] != nIndividualsTraits[j] )
          {
            misc.error("Error: An internal error was happened. traits with different numbers of individuals using diagonal Kernels is not allowed.", 0);
          }
          covar = new Matrix(diagonalDistribution, nIndividualsTraits[i], nIndividualsTraits[j], communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
          covar->fillDiagonal(1.);
        }
  
        this->V->insertCovarianceMatrix("E_" + i2s(i+1) + "_" + i2s(j+1), covar);
        
        computeEnvironmentalCovariances[std::pair<int, int>(i, j)] = true;
      }
    }
  }
  else
  {
    for(int i = 0; i<this->nPhenotypes; i++)
    {
      for(int j = i + 1; j<this->nPhenotypes; j++)
      {
        computeEnvironmentalCovariances[std::pair<int, int>(i, j)] = false;
      }
    }
  }
  
  
  for(int i = 0; i < globalEnvironmentalWeights.size(); i++)
  {
    if(globalEnvironmentalWeights[i] != NULL)
    {
      delete [] globalEnvironmentalWeights[i];
    }
  }
  globalEnvironmentalWeights.clear();
  
  ////////////////////////////////////////////////////
  // Insert variance groups
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    this->V->insertVarianceGroup("Phenotype_" + i2s(i+1), phenotypeVariances[i]);
    for(int j = i + 1; j<this->nPhenotypes; j++)
    {
      this->V->insertVarianceGroup("Phenotype_" + i2s(i+1) + "_" + i2s(j+1), 0.5*sqrt(phenotypeVariances[i]*phenotypeVariances[j]));
    }
  }
  
  ////////////////////////////////////////////////////
  // Insert variances
  
  for(int i = 0; i<kernels.size(); i++)
  {
    for(int j = 0; j<this->nPhenotypes; j++)
    {
      this->V->insertVariance("Var(" + kernelNames[i] + "_p" + i2s(j+1) + ")", "Phenotype_" + i2s(j+1), ParameterAttributes::variance, ParameterAttributes::genetic, phenotypeVariances[j]*heritabilities[j]*weights[i], std::set<std::string>());
      if( kernelTypes[i] == kernelSquaredExponential )
      {
        this->V->insertVariance("alpha0(" + kernelNames[i] + "_p" + i2s(j+1) + ")", "Phenotype_" + i2s(j+1), ParameterAttributes::parameter, ParameterAttributes::other, kernelAverages[kernelNames[i]], std::set<std::string>() );
      }
      
      for(int k = j + 1; k<this->nPhenotypes; k++)
      {
        std::set<std::string> constrainedDependingOnProductOfName;
	if( options.useCorrelations == false )
	{
	  constrainedDependingOnProductOfName.insert("Var(" + kernelNames[i] + "_p" + i2s(j+1) + ")");
	  constrainedDependingOnProductOfName.insert("Var(" + kernelNames[i] + "_p" + i2s(k+1) + ")");
	  this->V->insertVariance("Covar(" + kernelNames[i] + "_p" + i2s(j+1) + "-" + i2s(k+1) + ")", "Phenotype_" + i2s(j+1) + "_" + i2s(k+1), ParameterAttributes::covariance, ParameterAttributes::genetic, 0.5*sqrt( phenotypeVariances[j]*heritabilities[j]*weights[i]*phenotypeVariances[k]*heritabilities[k]*weights[i] ), constrainedDependingOnProductOfName );
	}
	else
	{
	  this->V->insertVariance("Cor(" + kernelNames[i] + "_p" + i2s(j+1) + "-" + i2s(k+1) + ")", "Phenotype_" + i2s(j+1) + "_" + i2s(k+1), ParameterAttributes::correlation, ParameterAttributes::genetic, 0.5, constrainedDependingOnProductOfName );
	}
	if( kernelTypes[i] == kernelSquaredExponential )
        {
          this->V->insertVariance("alpha0(" + kernelNames[i] + "_p" + i2s(j+1) + "-" + i2s(k+1) + ")", "Phenotype_" + i2s(j+1) + "_" + i2s(k+1), ParameterAttributes::parameter, ParameterAttributes::other, kernelAverages[kernelNames[i]], std::set<std::string>() );
        }
      }
    }
  }

  for(int i = 0; i<this->nPhenotypes; i++)
  {
    this->V->insertVariance("Var(E_p" + i2s(i+1) + ")", "Phenotype_" + i2s(i+1), ParameterAttributes::variance, ParameterAttributes::environment, phenotypeVariances[i]*(1.-heritabilities[i]), std::set<std::string>());
    
    for(int j = i + 1; j<this->nPhenotypes; j++)
    {
      if( computeEnvironmentalCovariances[std::pair<int, int>(i, j)] == true )
      {
        std::set<std::string> constrainedDependingOnProductOfName;
	if( options.useCorrelations == false )
	{
	  constrainedDependingOnProductOfName.insert("Var(E_p" + i2s(i+1) + ")");
	  constrainedDependingOnProductOfName.insert("Var(E_p" + i2s(j+1) + ")");
	  this->V->insertVariance("Covar(E_p" + i2s(i+1) + "-" + i2s(j+1) + ")", "Phenotype_" + i2s(i+1) + "_" + i2s(j+1), ParameterAttributes::covariance, ParameterAttributes::environment, 0.5*sqrt(phenotypeVariances[i]*(1.-heritabilities[i])*phenotypeVariances[j]*(1.-heritabilities[j])), constrainedDependingOnProductOfName );
	}
	else
	{
	  this->V->insertVariance("Cor(E_p" + i2s(i+1) + "-" + i2s(j+1) + ")", "Phenotype_" + i2s(i+1) + "_" + i2s(j+1), ParameterAttributes::correlation, ParameterAttributes::environment, 0.5, constrainedDependingOnProductOfName );
	}
      }
    }
  }
  
  this->V->updateVariancesDependenceIndexs();
  
  ////////////////////////////////////////////////////
  // Insert Elements
  
  for(int i = 0; i<kernels.size(); i++)
  {
    int shiftRows = 0;
    int shiftColumns = 0;
    for(int j = 0; j<this->nPhenotypes; j++)
    {
      this->V->insertElement(kernelNames[i], kernelNames[i] + "_" + i2s(j + 1), kernelTypes[i], kernelNames[i] + "_" + i2s(j + 1), 1., std::pair<int, int>(j, j), subMatrix(shiftRows, shiftColumns, nIndividualsTraits[j], nIndividualsTraits[j]), subMatrix(0, 0, nIndividualsTraits[j], nIndividualsTraits[j]));
      this->V->appendVarianceToElement(kernelNames[i] + "_" + i2s(j + 1), "Var(" + kernelNames[i] + "_p" + i2s(j + 1) + ")", ParameterAttributes::nochange);
      if( kernelTypes[i] == kernelSquaredExponential )
      {
        this->V->appendVarianceToElement(kernelNames[i] + "_" + i2s(j + 1), "alpha0(" + kernelNames[i] + "_p" + i2s(j+1) + ")", ParameterAttributes::nochange, ParameterAttributes::insideMatrix);
      }
      
      for(int k = j + 1; k<this->nPhenotypes; k++)
      {
        shiftColumns += nIndividualsTraits[ k - 1 ];
        this->V->insertElement(kernelNames[i], kernelNames[i] + "_" + i2s(j + 1) + "_" + i2s(k + 1), kernelTypes[i], kernelNames[i] + "_" + i2s(j + 1) + "_" + i2s(k + 1), 1., std::pair<int, int>(j, k), subMatrix(shiftRows, shiftColumns, nIndividualsTraits[j], nIndividualsTraits[k]), subMatrix(0, 0, nIndividualsTraits[j], nIndividualsTraits[k]));
	if( options.useCorrelations == false )
	{
	  this->V->appendVarianceToElement(kernelNames[i] + "_" + i2s(j + 1) + "_" + i2s(k + 1), "Covar(" + kernelNames[i] + "_p" + i2s(j + 1) + "-" + i2s(k + 1) + ")", ParameterAttributes::nochange);
	}
	else
	{
	  this->V->appendVarianceToElement(kernelNames[i] + "_" + i2s(j + 1) + "_" + i2s(k + 1), "Cor(" + kernelNames[i] + "_p" + i2s(j + 1) + "-" + i2s(k + 1) + ")", ParameterAttributes::nochange);
	  this->V->appendVarianceToElement(kernelNames[i] + "_" + i2s(j + 1) + "_" + i2s(k + 1), "Var(" + kernelNames[i] + "_p" + i2s(j + 1) + ")", ParameterAttributes::squareRoot);
	  this->V->appendVarianceToElement(kernelNames[i] + "_" + i2s(j + 1) + "_" + i2s(k + 1), "Var(" + kernelNames[i] + "_p" + i2s(k + 1) + ")", ParameterAttributes::squareRoot);
	}
	if( kernelTypes[i] == kernelSquaredExponential )
        {
          this->V->appendVarianceToElement(kernelNames[i] + "_" + i2s(j + 1) + "_" + i2s(k + 1), "alpha0(" + kernelNames[i] + "_p" + i2s(j+1) + "-" + i2s(k+1) + ")", ParameterAttributes::nochange, ParameterAttributes::insideMatrix);
        }
      }
      
      shiftRows += nIndividualsTraits[j];
      shiftColumns = shiftRows;
    }
  }

  int shiftRows = 0;
  int shiftColumns = 0;
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    this->V->insertElement("E", "E_" + i2s(i + 1), kernelEnvirontmental, "E_" + i2s(i + 1), 1., std::pair<int, int>(i, i), subMatrix(shiftRows, shiftColumns, nIndividualsTraits[i], nIndividualsTraits[i]), subMatrix(0, 0, nIndividualsTraits[i], nIndividualsTraits[i]));
    this->V->appendVarianceToElement("E_" + i2s(i + 1), "Var(E_p" + i2s(i + 1) + ")", ParameterAttributes::nochange);
    
    for(int j = i + 1; j<this->nPhenotypes; j++)
    {
      shiftColumns += nIndividualsTraits[ j - 1 ];
      if( computeEnvironmentalCovariances[std::pair<int, int>(i, j)] == true )
      {
        this->V->insertElement("E", "E_" + i2s(i + 1) + "_" + i2s(j + 1), kernelEnvirontmental, "E_" + i2s(i + 1) + "_" + i2s(j + 1), 1., std::pair<int, int>(i, j), subMatrix(shiftRows, shiftColumns, nIndividualsTraits[i], nIndividualsTraits[j]), subMatrix(0, 0, nIndividualsTraits[i], nIndividualsTraits[j]));
	if( options.useCorrelations == false )
	{
	  this->V->appendVarianceToElement("E_" + i2s(i + 1) + "_" + i2s(j + 1), "Covar(E_p" + i2s(i + 1) + "-" + i2s(j + 1) + ")", ParameterAttributes::nochange);
	}
	else
	{
	  this->V->appendVarianceToElement("E_" + i2s(i + 1) + "_" + i2s(j + 1), "Cor(E_p" + i2s(i + 1) + "-" + i2s(j + 1) + ")", ParameterAttributes::nochange);
	  this->V->appendVarianceToElement("E_" + i2s(i + 1) + "_" + i2s(j + 1), "Var(E_p" + i2s(i + 1) + ")", ParameterAttributes::squareRoot);
	  this->V->appendVarianceToElement("E_" + i2s(i + 1) + "_" + i2s(j + 1), "Var(E_p" + i2s(j + 1) + ")", ParameterAttributes::squareRoot);
	}
      }
    }
    
    shiftRows += nIndividualsTraits[i];
    shiftColumns = shiftRows;
  }
  
  
  //////////////////////////
  // Final adjustments
  if( this->GRMEigenVectors.size() != 0 && this->usingDiagonalKernels == false )
  {
    misc.error("Error: An internal error has happened. Inconsistent variables related to using diagonalized kernels.", 0);
  }
  if(this->usingDiagonalKernels == true)
  {
    misc.message << "Fitting REML with diagonalized matrices." << std::endl;
    if( this->GRMEigenVectors.size() != 1 )
    {
      misc.error("Error: An internal error has happened. Inconsistent variables related to using diagonalized kernels.", 0);
    }
    
    Matrix * evec = NULL;
    for(std::map<std::string, Matrix *>::iterator it = this->GRMEigenVectors.begin(); it != this->GRMEigenVectors.end(); ++it)
    {
      if(evec == NULL)
      {
        evec = it->second;
      }
      else
      {
        if(evec != it->second)
        {
          misc.error("Error: An internal error has happened. Unexpected different eigenvector matrices.", 0);
        }
      }
    }
    this->GRMEigenVectors[ "E" ] = evec;
  }
  
  this->individualBLUPNames.push_back("E");
  
  kernels.clear();
  
  this->V->setVarianceInitialValuesFromFile(options.initialVariancesFile);
  
  return true;
}

bool REML::prepare(Matrix* yparam, Matrix* Xparam, std::vector<Matrix*> & kernels, std::vector<KernelType> kernelTypes, double heritability, std::vector<double> weights)
{
  //Perform some checks and parameter initialization
  
  if(kernels.size() < 1)
  {
    misc.error("Error: An internal error was happened. No GRMs specified for REML analysis.", 0);
  }
  
  if(weights.size() == 0)
  {
    double equallyDistributedWeight = 1./double(kernels.size());
    for(int i = 0; i<kernels.size(); i++)
    {
      weights.push_back(equallyDistributedWeight);
    }
  }
  else
  {
    if( kernels.size() != weights.size() )
    {
      misc.error("Error: An internal error was happened. Kernel weights not properly defined.", 0);
    }
  }
  
  if(heritability < 0)
  {
    heritability = 0.5;
  }
  
  //Set reml type
  
  this->type = rawREMLType;
  
  //Set the number of phenotypes
  
  this->nPhenotypes = 1;
  
  //Check dimensions
  
  for(int i = 0; i<kernels.size(); i++)
  {
    if( kernels[i]->nGlobRows != yparam->nGlobRows || kernels[i]->nGlobCols != yparam->nGlobRows )
    {
      misc.error("Error: An internal error was happened. Different number of individuals in kernels and phenotypes when computing raw REML.", 0);
    }
  }
  if( yparam->nGlobRows != Xparam->nGlobRows )
  {
    misc.error("Error: An internal error was happened. Different number of individuals in covariates and phenotypes when computing raw REML.", 0);
  }
  
  int nTotalIndividuals = yparam->nGlobRows;
  
  ////////////////////////////////////////////////////
  // Set phenotype and covariate matrices
  
  this->y = yparam;
  this->X = Xparam;
  
  ////////////////////////////////////////////////////
  // Diagonal kernels?
  this->usingDiagonalKernels = false;
  bool allKernelsDiagonal = true;
  bool someKernelsDiagonal = false;
  for(int i = 0; i<kernels.size(); i++)
  {
    if(kernels[i]->distribution == diagonalDistribution)
    {
      someKernelsDiagonal = true;
    }
    else
    {
      allKernelsDiagonal = false;
    }
  }
  
  if( allKernelsDiagonal == false && someKernelsDiagonal == true )
  {
    misc.error("Error: An internal error was happened. Diagonal kernels mixed with nondiagonal ones have been found.", 0);
  }
  
  if( allKernelsDiagonal == true )
  {
    this->singlePrecisionInversion = false;
  }
  
  ////////////////////////////////////////////////////
  // Create the covariance matrix
  communicator->broadcast(&nTotalIndividuals, 1);
  
  this->dimension = nTotalIndividuals;
  
  this->V = new CovarianceMatrix(this->dimension, allKernelsDiagonal, this->nPhenotypes);
  
  ////////////////////////////////////////////////////
  // Insert the covariance matrices
  std::vector<std::string> kernelNames; //For storing kernels names before deleting them.
  std::map<std::string, double> kernelAverages;
  this->vIndividuals.clear();
  for(int i = 0; i<kernels.size(); i++)
  {
    this->V->insertCovarianceMatrix("K_" + i2s(i + 1), kernels[i]);
    kernelNames.push_back("K_" + i2s(i + 1));
    if(kernelTypes[i] == kernelSquaredExponential)
    {
      kernelAverages["K_" + i2s(i + 1)] = options.expKernelParameterInitialFactor/kernels[i]->elementsAverage();
    }
    kernels[i] = NULL;
  }
  this->vIndividuals.push_back( std::vector<Individual>() );
  
  Matrix * temp;
  if( allKernelsDiagonal == false )
  {
    temp = new Matrix(cyclicDistribution, nTotalIndividuals, nTotalIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  }
  else
  {
    temp = new Matrix(diagonalDistribution, nTotalIndividuals, nTotalIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  }
  temp->fillDiagonal(1.);
  
  this->V->insertCovarianceMatrix("E", temp);

  
 ////////////////////////////////////////////////////
  // Insert variance groups
  double phenotypeVariance = computeInitialVariance(this->y, this->X);
  this->V->insertVarianceGroup("Phenotype_1", phenotypeVariance);
  
  ////////////////////////////////////////////////////
  // Insert variances
  
  for(int i = 0; i<kernels.size(); i++)
  {
    this->V->insertVariance("Var(" + kernelNames[i] + ")", "Phenotype_1", ParameterAttributes::variance, ParameterAttributes::genetic, phenotypeVariance*heritability*weights[i], std::set<std::string>());
    if( kernelTypes[i] == kernelSquaredExponential )
    {
      this->V->insertVariance("alpha0(" + kernelNames[i] + ")", "Phenotype_1", ParameterAttributes::parameter, ParameterAttributes::other, kernelAverages[kernelNames[i]], std::set<std::string>() );
    }
  }
  this->V->insertVariance("Var(E)", "Phenotype_1", ParameterAttributes::variance, ParameterAttributes::environment, phenotypeVariance*(1.-heritability), std::set<std::string>());
  
  this->V->updateVariancesDependenceIndexs();
  
  ////////////////////////////////////////////////////
  // Insert Elements
  
  for(int i = 0; i<kernels.size(); i++)
  {
    this->V->insertElement(kernelNames[i], kernelNames[i], kernelTypes[i], kernelNames[i], 1., std::pair<int, int>(0, 0), subMatrix(0, 0, nTotalIndividuals, nTotalIndividuals), subMatrix(0, 0, nTotalIndividuals, nTotalIndividuals));
    this->V->appendVarianceToElement(kernelNames[i], "Var(" + kernelNames[i] + ")", ParameterAttributes::nochange);
    if( kernelTypes[i] == kernelSquaredExponential )
    {
      this->V->appendVarianceToElement(kernelNames[i], "alpha0(" + kernelNames[i] + ")", ParameterAttributes::nochange, ParameterAttributes::insideMatrix);
    }
  }

  int shiftRows = 0;
  int shiftColumns = 0;

  this->V->insertElement("E", "E", kernelEnvirontmental, "E", 1., std::pair<int, int>(0, 0), subMatrix(0, 0, nTotalIndividuals, nTotalIndividuals), subMatrix(0, 0, nTotalIndividuals, nTotalIndividuals));
  this->V->appendVarianceToElement("E", "Var(E)", ParameterAttributes::nochange);
 
  kernels.clear();
  
  return true;
}

void REML::createCovarMatrix(Matrix * srcCovarMatrix, Matrix * resultCovarMatrix)
{
  int dimension1 = srcCovarMatrix->nGlobRows;
  int dimension2 = srcCovarMatrix->nGlobCols;
  Matrix * srcCovarMatrixTransp = new Matrix(cyclicDistribution);
  srcCovarMatrixTransp->transpose(srcCovarMatrix);
  resultCovarMatrix->joinMatrices(srcCovarMatrix, subMatrix(0, dimension1, dimension1, dimension2), srcCovarMatrixTransp, subMatrix(dimension1, 0, dimension2, dimension1));
  delete srcCovarMatrixTransp;
}

double REML::computeInitialVariance(Matrix* yMatrix, Matrix* XMatrix)
{
  Matrix * tempXtXi = new Matrix();
  tempXtXi->multiply(XMatrix, 'T', XMatrix, 'N');
  bool success = tempXtXi->symmetricInvert();
  if(success == false)
  {
    delete tempXtXi;
    
    misc.message << "WARNING: Xt*X is not invertible. It is possible that X has linearly dependent columns." << std::endl;
    
    return computeVariance(yMatrix);
  }
  Matrix *tempXty = new Matrix();
  tempXty->multiply(XMatrix, 'T', yMatrix, 'N');
  
  Matrix *tempBetas = new Matrix();
  tempBetas->multiply(tempXtXi, 'N', tempXty, 'N');
  
  Matrix *tempResid = new Matrix();
  tempResid->multiply(XMatrix, 'N', tempBetas, 'N');
  tempResid->add(yMatrix, -1., 1.);
  
  double variance = computeVariance(tempResid);
  
  delete tempXtXi;
  delete tempXty;
  delete tempBetas;
  delete tempResid;
  
  return variance;
}

void REML::transformUsingEigenVectors(std::string name, char t, Matrix ** m, int numberBlocks)
{
  if( this->usingDiagonalKernels == false || this->GRMEigenVectors.count(name) == 0 )
  {
    misc.error("Error: An internal error was happened. The matrix cannot be transformed using eigenvectors. This is not a diagonal analysis or eigenvectors not present.", 0);
  }
  
  Matrix * eigenVectors = this->GRMEigenVectors[name];
  
  if( numberBlocks == 1)
  {
    Matrix * temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    temp->multiply(eigenVectors, t, *m, 'N');
    delete *m;
    *m = temp;
  }
  else
  {
    if( numberBlocks*eigenVectors->nGlobCols != (*m)->nGlobRows )
    {
      misc.error("Error: An internal error was happened. The matrix cannot be transformed using eigenvectors. Dimensions do not match.", 0);
    }
    
    Matrix * temp = new Matrix((*m)->distribution, (*m)->nGlobRows, (*m)->nGlobCols);
    temp->fillWithConstant(0.);
    
    int rowShift = 0;
    for(int ib = 0; ib < numberBlocks; ib++)
    {
      Matrix * mBlock = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, eigenVectors->nGlobRows, (*m)->nGlobCols);
      mBlock->fillWithConstant(0.);
      subMatrix smDest(mBlock);
      subMatrix smSrc(rowShift, 0, eigenVectors->nGlobRows, (*m)->nGlobCols);
      mBlock->add(*m, 1., 1., smDest, smSrc);
      
      Matrix * transformedmBlock = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
      transformedmBlock->multiply(eigenVectors, t, mBlock, 'N');
      
      temp->add(transformedmBlock, 1., 1., smSrc, smDest);
      
      rowShift += eigenVectors->nGlobRows;
    }
    
    delete *m;
    *m = temp;
  }
}

void REML::addCovariatesNames(Covariate * cm, std::string prefix, int shift)
{
  if( communicator->mpiRoot == true )
  {
    CovariateNames tempCovariateNames;
    int idx = shift;
    for( int i = 0; i < cm->meanNames.size(); i++ )
    {
      std::string tempName = prefix + cm->meanNames[i];
      std::pair<std::string, int> tempPair(tempName, idx);
      tempCovariateNames.meanNames.push_back(tempPair);
      idx++;
    }
    for( int i = 0; i < cm->discreteCovarNames.size(); i++ )
    {
      std::string tempName = prefix + cm->discreteCovarNames[i];
      std::pair<std::string, int> tempPair(tempName, idx);
      tempCovariateNames.discreteCovarNames.push_back(tempPair);
      idx++;
    }
    for( int i = 0; i < cm->quantitativeCovarNames.size(); i++ )
    {
      std::string tempName = prefix + cm->quantitativeCovarNames[i];
      std::pair<std::string, int> tempPair(tempName, idx);
      tempCovariateNames.quantitativeCovarNames.push_back(tempPair);
      idx++;
    }
    this->covariateNames.push_back(tempCovariateNames);

    //Check there are not repeated indices.    
    std::set<int> checkIndices;
    for(int idxPheno = 0; idxPheno<this->covariateNames.size(); idxPheno++)
    {
      for( int i = 0; i < this->covariateNames[idxPheno].meanNames.size(); i++ )
      {
        int testIdx = this->covariateNames[idxPheno].meanNames[i].second;
        if( checkIndices.find(testIdx) != checkIndices.end() )
        {
          misc.error("Error: An internal error was happened. When creating covariates names indexes. Repeated index.", 0);
        }
        checkIndices.insert(testIdx);
      }
      for( int i = 0; i < this->covariateNames[idxPheno].discreteCovarNames.size(); i++ )
      {
        int testIdx = this->covariateNames[idxPheno].discreteCovarNames[i].second;
        if( checkIndices.find(testIdx) != checkIndices.end() )
        {
          misc.error("Error: An internal error was happened. When creating covariates names indexes. Repeated index.", 0);
        }
        checkIndices.insert(testIdx);
      }
      for( int i = 0; i < this->covariateNames[idxPheno].quantitativeCovarNames.size(); i++ )
      {
        int testIdx = this->covariateNames[idxPheno].quantitativeCovarNames[i].second;
        if( checkIndices.find(testIdx) != checkIndices.end() )
        {
          misc.error("Error: An internal error was happened. When creating covariates names indexes. Repeated index.", 0);
        }
        checkIndices.insert(testIdx);
      }
    }//End for checking
  } //End mpi root
}

std::string REML::addSuffix(std::string suffix)
{
  if(this->nPhenotypes > 1)
  {
    return suffix;
  }
  return "";
}

void REML::setyList(std::vector<Matrix *> & ys)
{
  if(ys.size() == 0)
  {
    misc.error("Error: An internal error has happened. The list of phenotypes is empty.", 0);
  }
  deleteyList();
  this->yList = ys;
  
  if( this->usingDiagonalKernels == true )
  {
    misc.message << "Transforming phenotypes..." << std::endl;
    if(this->GRMEigenVectors.size() != 1)
    {
      misc.error("Error: An internal error has happened. Unexpected number of eigenvectors when preparing sampled phenotypes.", 0);
    }
    std::string name = this->GRMEigenVectors.begin()->first;
    
    for(int i = 0; i < this->yList.size(); i++)
    {
      transformUsingEigenVectors(name, 'T', &(this->yList[i]), this->nPhenotypes);
    }
  }
  
  this->y = this->yList[0];
  ys.clear();
}

void REML::deleteyList()
{
  bool deleteY = true;
  for(int i = 0; i < this->yList.size(); i++)
  {
    if(this->yList[i] == this->y)
    {
      deleteY = false;
    }
    delete this->yList[i];
  }
  if(deleteY == true && this->y != NULL)
  {
    delete this->y;
  }
  this->y = NULL;
  this->yList.clear();
}

void REML::computeREMLWithReducedModels()
{
  double fullModelLogLikelihood = computeREML();
  
  if(this->success == false)
  {
    return;
  }
  
  if( options.skipComputeReducedModels == false )
  {
    misc.message << "Computing reduced models for testing..." << std::endl;
    
    int remlStepsToUnfixExpKernelParameterTemp = options.remlStepsToUnfixExpKernelParameter;
    options.remlStepsToUnfixExpKernelParameter = -1;
    
    std::set<std::string> fullModelWarnings = this->warnings;
    
    this->V->setInitialVariancesWithCurrentValues();
    
    std::vector< std::string > reducedModelsName;
    std::vector< int > reducedModelsNumberOfVariances;
    std::vector< double > reducedModelsLogLikelihoods;
    std::vector< bool > reducedModelsConverged;
    std::vector< std::set<std::string> > reducedModelsWarnings;
    
    int fullModelNumberOfVariances = this->V->variances.size();
    this->writeResults = false;
    bool writeSummaryTemp = this->writeSummary;
    this->writeSummary = true;
    bool tempOptionsFirstStepEM = options.firstStepEM;
    options.firstStepEM = false;            //If expected start near convergence, avoid EM step can improve performance.
    this->singlePrecisionInversion = false; //If expected low number of iterations, maybe it is better not use single precission.
    
    for( int i = 0; i < this->elementsToTest.size(); i++)
    {
      this->V->backupFullModel();
      this->V->deleteElementsWithSubCovarianceId( this->elementsToTest[i] );
      this->V->setInitialVariancesWithCurrentValues();
      
      std::string baseOutput = options.outFile;      
      options.outFile += "." + spacetab2underscore(this->elementsToTest[i]) + "_removed";
      
      double reducedModelLogLikelihood = computeREML();
      
      options.outFile = baseOutput;
      
      if(this->success == true)
      {
        reducedModelsName.push_back( this->elementsToTest[i]  + " removed" );
        reducedModelsNumberOfVariances.push_back( this->V->variances.size() );
        reducedModelsLogLikelihoods.push_back( reducedModelLogLikelihood );
        reducedModelsConverged.push_back(true);
        reducedModelsWarnings.push_back( this->warnings );
      }
      else
      {
        reducedModelsName.push_back( "Warning: The model with " + this->elementsToTest[i]  + " removed has not converged." );
        reducedModelsNumberOfVariances.push_back( -1 );
        reducedModelsLogLikelihoods.push_back( 0. );
        reducedModelsConverged.push_back(false);
        reducedModelsWarnings.push_back( std::set<std::string>() );
        misc.message << "Sorry, the REML for the reduced model has not converged. The LRT for " << this->elementsToTest[i] << " cannot be computed." << std::endl;
      }
      
      this->V->restoreBackupFullModel();  
    }
    
    for( int i = 0; i < options.restrictedCovariances.size(); i++ )
    {
      std::string baseId = options.restrictedCovariances[i].baseId;
      std::string p1 = i2s(options.restrictedCovariances[i].p1);
      std::string p2 = i2s(options.restrictedCovariances[i].p2);
      
      std::string element = baseId + "_" + p1 + "_" + p2;
      std::string covarName;
      std::string outCovarName;
      if(options.useCorrelations == false)
      {
	covarName = "Covar(" + baseId + "_p" + p1 + "-" + p2 + ")";
        outCovarName = "Correlation from Covar(" + baseId + "_p" + p1 + "-" + p2 + ")";
      }
      else
      {
	covarName = "Cor(" + baseId + "_p" + p1 + "-" + p2 + ")";
        outCovarName = "Cor(" + baseId + "_p" + p1 + "-" + p2 + ")";
      }
      std::string varName1 = "Var(" + baseId + "_p" + p1 + ")";
      std::string varName2 = "Var(" + baseId + "_p" + p2 + ")";
      
      if( this->V->checkVarianceExists(covarName) == false )
      {
        misc.message << "The variance " << covarName << " does not exist in the current model. Skipping..." << std::endl;
        continue;
      }
      
      this->V->backupFullModel();
      
      this->V->clearElementVariances(element);
      this->V->deleteVariance(covarName);
      this->V->changeElementConstantFactor(element, options.restrictedCovariances[i].correlation);
      this->V->appendVarianceToElement(element, varName1, ParameterAttributes::squareRoot);
      this->V->appendVarianceToElement(element, varName2, ParameterAttributes::squareRoot);
      
      std::string baseOutput = options.outFile;      
      options.outFile += "." + spacetab2underscore(outCovarName) + "_fixedto_" + getString(options.restrictedCovariances[i].correlation);
      
      double reducedModelLogLikelihood = computeREML();
      
      options.outFile = baseOutput;
      
      if(this->success == true)
      {
        reducedModelsName.push_back( outCovarName  + " fixed to "  + getString(options.restrictedCovariances[i].correlation) );
        reducedModelsNumberOfVariances.push_back( this->V->variances.size() );
        reducedModelsLogLikelihoods.push_back( reducedModelLogLikelihood );
        reducedModelsConverged.push_back(true);
        reducedModelsWarnings.push_back( this->warnings );
      }
      else
      {
        reducedModelsName.push_back( "Warning: The model with " + outCovarName  + " fixed to "  + getString(options.restrictedCovariances[i].correlation) + " has not converged." );
        reducedModelsNumberOfVariances.push_back( -1 );
        reducedModelsLogLikelihoods.push_back( 0. );
        reducedModelsConverged.push_back(false);
        reducedModelsWarnings.push_back( std::set<std::string>() );
        misc.message << "Sorry, the REML for the reduced model has not converged. The LRT for fixing " << outCovarName << " to " << options.restrictedCovariances[i].correlation << " cannot be computed." << std::endl;
      }
      
      this->V->restoreBackupFullModel();
    }
    
    for( int i = 0; i < options.reducedModelsOnlyCovariances.size(); i++)
    {
      this->V->backupFullModel();
      bool elementExists = false;
      std::set<std::string> tempIds;
      for(int j = 0; j < this->V->elements.size(); j++)
      {
        tempIds.insert(this->V->elements[j].subCovarianceId);
      }
      for(std::set<std::string>::iterator it = tempIds.begin(); it != tempIds.end(); ++it)
      {
        if( *it != "E" && *it != options.reducedModelsOnlyCovariances[i] )
        {
          this->V->deleteElementsWithSubCovarianceId( *it );
        }
        if( *it == options.reducedModelsOnlyCovariances[i] )
        {
          elementExists = true;
        }
      }
      if(elementExists == false)
      {
        misc.message << "WARNING: The covariance element " << options.reducedModelsOnlyCovariances[i] << " does not exist. The reduced model with only this element will not be computed. Please, check the command line options." << std::endl; 
        this->V->restoreBackupFullModel();  
        continue;
      }
      this->V->setInitialVariancesWithCurrentValues();
      
      std::string baseOutput = options.outFile;      
      options.outFile += ".only_" + spacetab2underscore(options.reducedModelsOnlyCovariances[i]);
      
      double reducedModelLogLikelihood = computeREML();
      
      options.outFile = baseOutput;
      
      if(this->success == true)
      {
        reducedModelsName.push_back( "Only " + options.reducedModelsOnlyCovariances[i] );
        reducedModelsNumberOfVariances.push_back( this->V->variances.size() );
        reducedModelsLogLikelihoods.push_back( reducedModelLogLikelihood );
        reducedModelsConverged.push_back(true);
        reducedModelsWarnings.push_back( this->warnings );
      }
      else
      {
        reducedModelsName.push_back( "Warning: The model with only " + options.reducedModelsOnlyCovariances[i]  + " has not converged." );
        reducedModelsNumberOfVariances.push_back( -1 );
        reducedModelsLogLikelihoods.push_back( 0. );
        reducedModelsConverged.push_back(false);
        reducedModelsWarnings.push_back( std::set<std::string>() );
        misc.message << "Sorry, the REML for the reduced model has not converged. The LRT for the model with only " << options.reducedModelsOnlyCovariances[i] << " cannot be computed." << std::endl;
      }
      
      this->V->restoreBackupFullModel();  
    }
    
    if(communicator->mpiRoot)
    {
      Message message(options.outFile + ((this->useMLinsteadOfREML==false)?".reml.tests":".ml.tests"));
      for(int i = 0; i < reducedModelsName.size(); i++)
      {
        double LogRatio = 2.0*(fullModelLogLikelihood - reducedModelsLogLikelihoods[i]);
        
        int df = fullModelNumberOfVariances - reducedModelsNumberOfVariances[i];
        
        if( i != 0)
        {
          message << std::endl;
        }
        message << "#  " << reducedModelsName[i] << std::endl;
        for(std::set<std::string>::iterator it = fullModelWarnings.begin(); it != fullModelWarnings.end(); ++it)
        {
          message << "#  WARNING (full model): " << *it << std::endl;
        }
        for(std::set<std::string>::iterator it = reducedModelsWarnings[i].begin(); it != reducedModelsWarnings[i].end(); ++it)
        {
          message << "#  WARNING (reduced model): " << *it << std::endl;
        }
        message << "#------------------------------------------------------------" << std::endl;
        message << "Full_Model_logL    " << std::setprecision(10) << fullModelLogLikelihood << std::endl;
        if(reducedModelsConverged[i] == true && LogRatio >= 0.)
        {
          message << "Reduced_Model_logL " << std::setprecision(10) << reducedModelsLogLikelihoods[i] << std::endl;
          message << "LRT                " << std::setprecision(5) << LogRatio << std::endl;
          message << "df                 " << std::setprecision(1) << df << std::endl;
          message << "(0.5)*Pval         " << std::setprecision(4) << 0.5*chi1_CDF(df, LogRatio) << std::endl;
        }
        else if(LogRatio < 0.)
        {
          message << "Reduced_Model_logL " << std::setprecision(10) << reducedModelsLogLikelihoods[i] << std::endl;
          message << "LRT                " << std::setprecision(5) << LogRatio << std::endl;
          message << "df                 " << std::setprecision(1) << df << std::endl;
          message << "(0.5)*Pval         " << "NA" << std::endl;
        }
        else
        {
          message << "Reduced_Model_logL " << "NA" << std::endl;
          message << "LRT                " << "NA" << std::endl;
          message << "df                 " << std::setprecision(1) << df << std::endl;
          message << "(0.5)*Pval         " << "NA" << std::endl;
        }
      }
    }
    
    this->writeSummary = writeSummaryTemp;
    options.firstStepEM = tempOptionsFirstStepEM;
    options.remlStepsToUnfixExpKernelParameter = remlStepsToUnfixExpKernelParameterTemp;
  }
}

double REML::computeREML()
{
  bool iterate = true;
  int nIterations = 0;
  double logLikelihoodDifference;
  
  this->V->reinitializeVariances();
  
  this->logLikelihood = -1.e50;
  
  this->currentStep = 0;
  
  deleteIntermediateMatrices();
  this->ViX = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->XtViX_i = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->P = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->PDiagonal = new Matrix(diagonalDistribution);
  this->Py = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->subVPy = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->dimension, this->V->variances.size());
  this->AI = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->V->variances.size(), this->V->variances.size());
  this->yPsubVPy_trPsubV = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->V->variances.size(), 1);
  this->yPsubVPy_trVisubV = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->V->variances.size(), 1);
  this->VisubVTraces = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->V->variances.size(), 1);
  this->VisubVTraces->fillWithConstant(0.);
  
  this->secondDerivativesMatrixType = REMLParameters::undefined;
  
  this->delta = std::vector<double>(this->V->variances.size(), 0.);
  
  this->constrainedVariances.clear();
  
  if( this->useMLinsteadOfREML == false )
  {
    this->sREMLML = "REML";
  }
  else
  {
    this->sREMLML = "ML";
  }
  
  misc.message << "Starting " + this->sREMLML + " iterations" << ((this->usingDiagonalKernels == true)?" (using diagonal covariance matrix)":"") << "..." << std::endl;
  misc.message << "Initial variance values:            " << this->V->getStringVarianceValues() << std::endl;
  misc.message << std::setw(14) << "Step   M " << std::setw(3+4) << "Time" << std::setw(options.logFieldWidth+1) << "Prev. LgL " << this->V->getStringVarianceNames() << std::endl;
  
  this->success = true;
  
  this->warnings.clear();
  
  misc.setGetElapsedTime("REMLIteration");
  misc.setGetElapsedTime("REMLAnalysis");
  
  int variancesConstrainedMethod = 0;
  
  //Fix variances or parameters?
  fixVariancesAndParameters();
  
  //Start REML iterations
  while(iterate)
  {
    this->logLikelihoodConverged = false;
    this->variancesConverged = false;
    this->gradientConverged = false;
    
    std::stringstream stepResults;
    this->stepModifications = "";
    //REML step
    stepResults << std::setw(3) << nIterations + 1;
    if( ( ( nIterations == 0 && options.firstStepEM == true ) || options.REMLMethod == 1) && this->useMLinsteadOfREML == false )
    {
      stepResults << std::setw(5) << "EM";
      emREMLStep();
    }
    else
    {
      stepResults << std::setw(5) << ((this->useMLinsteadOfREML == false)?"AI":"F");
      aiREMLStep();
    }
    if(this->success == false)
    {
      iterate = false;
      break;
    }
    
    //delta loglikelihood
    this->logLikelihoodDifference = computeLogLikelihood();
    
    //Constrain variances
    int nConstrained = this->V->constrainVariancesM1(this->constrainedVariances);
    variancesConstrainedMethod = 0;
    std::stringstream constrainedMessage;
    if(nConstrained != 0)
    {
      constrainedMessage << " (" << nConstrained << " of " << this->V->variances.size() << " parameters constrained)";
      variancesConstrainedMethod = 1;
    }
    double test = double(nConstrained)/double(this->V->variances.size());
    if(test > 0.5 && nIterations == 0)
    {
      misc.message << "Error: More than half of the parameters are constrained in the first step. " + this->sREMLML + " stopped.";
      this->success = false;
    }
    else if( test > 0.5 && nIterations != 0 )
    {
      if(options.remlGCTAMode == false)
      {
        double scalingFactor = this->V->constrainVariancesM3(this->oldVariances, this->delta);
        nConstrained = 0;
        constrainedMessage.str( std::string() );
        constrainedMessage.clear();
        constrainedMessage << " (parameters constrained by a scaling factor: " << scalingFactor << ")";
        variancesConstrainedMethod = 2;
        
        if( options.allowFixingVariancesToZero == true )
        {
          int nFixedToZero = this->V->fixVariancesToZero();
          if( nFixedToZero != 0 )
          {
            constrainedMessage.str( std::string() );
            constrainedMessage.clear();
            constrainedMessage << " (parameters constrained by a scaling factor, " << nFixedToZero << " variances fixed to 0.)";
          }
        }
      }
      else
      {
        misc.error("Error: More than half of the parameters are constrained in the first step. " + this->sREMLML + " stopped.", 0);
      }
    }
    
    //Write REML step results
    if(this->singlePrecisionInversion == true)
    {
      this->stepModifications += "s";
    }
    stepResults << std::setw(5) << this->stepModifications;
    stepResults << std::setprecision(1) << std::setw(8) << misc.setGetElapsedTime("REMLIteration");
    stepResults << std::setprecision(options.logOutputPrecision + 2) << std::setw(options.logFieldWidth) << this->logLikelihood;
    misc.message << stepResults.str() << this->V->getStringVarianceValues() << constrainedMessage.str() << std::endl;
    
    
    //UnfixVariances?
    unfixVariancesAndParameters();
    
    //REML converged?
    this->logLikelihoodConverged = ((this->logLikelihoodDifference < 1e-4) && (this->logLikelihoodDifference > -1e-2));
    this->variancesConverged = (allVariancesRelativeDifferencesLowerThan(options.varianceConvergenceThreshold) == true);
    communicator->broadcast(&this->logLikelihoodConverged);
    communicator->broadcast(&this->variancesConverged);
    this->gradientConverged = gradientLowerThanThreshold(options.gradientConvergenceThreshold);
    //if(varianceDifferenceNorm() < options.varianceConvergenceThreshold && this->logLikelihoodDifference < 1e-4 && this->logLikelihoodDifference > -1e-2)
    if( this->logLikelihoodConverged && this->variancesConverged && variancesConstrainedMethod != 2 && ( variancesConstrainedMethod != 1 || options.allowConvergenceWithConstrainedVars == true ) && getNVariancesFixed() == 0)
    {
      if(this->singlePrecisionInversion == false)
      {
        iterate = false;
        if(variancesConstrainedMethod == 1)
        {
          std::string temp = "";
          std::string tempsep = "";
          for (std::set<std::string>::iterator it = this->constrainedVariances.begin(); it != this->constrainedVariances.end(); ++it)
          {
            temp += tempsep + *it;
            tempsep = ", ";
          }
          this->warnings.insert("There are " + i2s(nConstrained) + " parameters with constrained values < " + temp + " >. Please, check the logs.");
        }
        if( this->gradientConverged == false )
	{
	  this->warnings.insert("The gradient did not converge to a threshold of " + getString(options.gradientConvergenceThreshold));
	}
      }
      else
      {
        this->singlePrecisionInversion = false;
      }
    }
    
    //Switch to double precision?
    if( this->singlePrecisionInversion == true && allVariancesRelativeDifferencesLowerThan(options.varianceConvergenceThreshold/10.) == true )
    {
      this->singlePrecisionInversion = false;
    }
    
    //Other iteration?
    nIterations++;
    if( (nIterations>=options.maxREMLIterations) && (iterate == true) )
    {
      this->success = false;
    }
    
    //If something failed, stop iterations.
    if(this->success == false)
    {
      iterate = false;
    }
    
    this->currentStep++;
  }
  
  //REML converged?
  if(this->success == true && this->type != rawREMLType)
  {
    misc.message << this->sREMLML + " analysis has been finished with success after " << misc.setGetElapsedTime("REMLAnalysis", true) << "! (logL: " << std::setprecision(10) << this->logLikelihood << ")." << std::endl;
    
    if(this->warnings.size() != 0)
    {
      for(std::set<std::string>::iterator it = this->warnings.begin(); it != this->warnings.end(); ++it)
      {
        misc.message << "WARNING: " << *it << std::endl;
      }
    }
    
    if(this->writeResults == true)
    {
      computeSummary();
    
      if(options.computeIndividualsBLUP)
      {
        computeIndividualsBLUP();
      }
      if(options.computeSNPsBLUP)
      {
        computeSNPsBLUP();
      }
      if(options.computeBLUE)
      {
        computeBLUE();
      }
      experimentalEpistasisPrediction();
    }
    else
    {
      if(this->writeSummary == true)
      {
        computeSummary();
      }
      if(this->writeBLUE == true)
      {
        computeBLUE();
      }
    }
  }
  
  if(this->success == false)
  {
    if( variancesConstrainedMethod != 2 )
    {
      if( this->logLikelihoodConverged == false )
      {
        misc.message << "Sorry, " + this->sREMLML + " failed to converge..." << std::endl;
      }
      else
      {
        misc.message << "Sorry, despite the logL seems to converge to: " << std::setprecision(10) << this->logLikelihood << ", " + this->sREMLML + " failed to converge..."  << std::endl;
      }
      misc.message << "You can try to set a less stringent variance and gradient convergence thresholds by using the --variance-threshold and --gradient-threshold options. You can also increase the maximum number of iterations with the --reml-maxit option." << std::endl;
    }
    else
    {
      misc.message << "Sorry, all parameters have been constrained by a scaling factor. " + this->sREMLML + " failed to converge..." << std::endl;
    }
    
    if( options.useLogLogisticScale == false)
    {
      if( this->nPhenotypes == 1 )
      {
        misc.message << "You can also try with --use-log-logistic option to improve convergence." << std::endl;
      }
      else
      {
        misc.message << "You can also try with --use-correlations and --use-log-logistic options to improve convergence." << std::endl;
      }
    }
    
    if( variancesConstrainedMethod != 2 )
    {
      showNonConvergenceCauses();
    }
  }
  
  if( options.verboseLog == true )
  {
    verboseLogOutput1();
  }
  
  if(this->type != rawREMLType)
  {
    deleteIntermediateMatrices();
  }
  
  return this->logLikelihood;
}

void REML::computePMatrix()
{
  if( this->V->mStatus != inverseCovarianceMatrix && this->V->mStatus != inverseCovarianceMatrixByBlocks )
  {
    misc.error("Error: An internal error was happened. P matrix cannot be computed if covariance matrix is not inverted.", 0);
  }
  
  if( this->V->mStatus == inverseCovarianceMatrix )
  {
    this->ViX->multiply(this->V->m, 'N', this->X, 'N');
  }
  else
  {
    BlockMatrix ViXtemp;
    ViXtemp.multiply(this->V->mInBlocks, this->X);
    delete this->ViX;
    this->ViX = ViXtemp.block2distributed();
  }
  
  this->XtViX_i->multiply(this->X, 'T', this->ViX, 'N');
  this->XtViX_i->symmetric = true;
  this->XtViX_i->uplo = 'B';
  
  bool inverted = this->XtViX_i->symmetricInvert(&this->logDetXtViX);
  if(inverted == false)
  {
    this->XtViX_i->multiply(this->X, 'T', this->ViX, 'N');
    this->XtViX_i->symmetric = true;
    this->XtViX_i->uplo = 'B';
    inverted = this->XtViX_i->invert(&this->logDetXtViX);
    if(inverted == false)
    {
      misc.message << "Error: XtViX matrix can not be inverted. " + this->sREMLML + " iterations can not continue. Please check that there are not covariates linearly dependent. It may be also possible that after filtering missing individuals which do not overlap between input files, there is a covariate category not represented. You can solve this situation by prefiltering the non-overlapping individuals in the covariates file." << std::endl;
      this->success = false;
      return;
    }
  }
  
  if(this->V->m->distribution != diagonalDistribution || this->V->mStatus != inverseCovarianceMatrix)
  {
    Matrix * temp = new Matrix(cyclicDistribution);
    temp->multiply(this->ViX, 'N', this->XtViX_i, 'N');
    this->P->multiply(temp, 'N', this->ViX, 'T');
    this->P->symmetric = true;
    this->P->uplo = 'B';
    if( this->V->mStatus == inverseCovarianceMatrix )
    {
      this->P->add(this->V->m, -1., 1.); //These two steps maybe can be optimized to a single one addapting the multiply function.
    }
    else
    {
      Matrix * Vitemp = this->V->mInBlocks.block2distributed();
      this->P->add(Vitemp, -1., 1.); //These two steps maybe can be optimized to a single one addapting the multiply function.
      delete Vitemp;
    }
    
    delete temp;
    
    this->Py->multiply(this->P, 'N', this->y, 'N');
  }
  else
  {
    if(this->P != NULL)
    {
      delete this->P;
      this->P = NULL;
    }
    
    computePxMatrix(this->Py, this->y);
    
    this->PDiagonal->diagonalOfABAt(this->ViX, this->XtViX_i);
    this->PDiagonal->add(this->V->m, -1., 1.);
  }
}

void REML::computePxMatrix(Matrix * Px, Matrix * x)
{
  if( this->V->mStatus != inverseCovarianceMatrix )
  {
    misc.error("Error: An internal error was happened. Px matrix cannot be computed if covariance matrix is not inverted.", 0);
  }
  
  Matrix * ViXtx = new Matrix();
  ViXtx->multiply(this->ViX, 'T', x, 'N');
  Matrix * temp = new Matrix();
  temp->multiply(this->XtViX_i, 'N', ViXtx, 'N');
  delete ViXtx;

  Px->multiply(this->ViX, 'N', temp, 'N');
  temp->multiply(this->V->m, 'N', x, 'N');
  
  Px->add(temp, -1., 1.);
  
  delete temp;
}

/*void REML::computePyMatrix()
{
  this->Py->multiply(this->P, 'N', this->y, 'N');
}*/

void REML::computeytPy()
{
  Matrix * temp = new Matrix(cyclicDistribution);
  
  temp->multiply(this->y, 'T', this->Py, 'N');
  temp->gatherMatrix(&(this->ytPy));
  communicator->broadcast(&(this->ytPy), 1);
  
  delete temp;
}

void REML::computeSubVPyMatrix()
{
  this->subVPy->fillWithConstant(0.);
  Matrix * temp = new Matrix();
  for(int i = 0; i<this->V->variances.size(); i++)
  {
    Matrix * subV = this->V->computeDerivateCovariance(i);
    temp->multiply(subV, 'N', this->Py, 'N');
    this->subVPy->add(temp, 1., 1., subMatrix(0, i, this->dimension, 1), subMatrix(temp));
    //this->subVPy->multiply(subV, 'N', this->Py, 'N', 1., subMatrix(0, i, this->dimension, 1));
  }
  delete temp;
}


void REML::computeAIMatrix()
{
  if(this->useMLinsteadOfREML == false)
  {
    if( this->V->checkMoreThanOneVariance() == false || options.forceUseREMLAIWhenNoLinearCovariance == true )
    {
      Matrix * PsubVPy = new Matrix(cyclicDistribution);

      if(this->P != NULL)
      {
        PsubVPy->multiply(P, 'N', this->subVPy, 'N');
      }
      else
      {
        computePxMatrix(PsubVPy, this->subVPy);
      }
      this->AI->multiply(this->subVPy, 'T', PsubVPy, 'N', 0.5);
      
      computeAIMatrixCrossedDerivatesCorrection();
      
      this->secondDerivativesMatrixType = REMLParameters::AI;
      
      delete PsubVPy;
    }
    else
    {
      computeREMLFMatrix();
    }
  }
  else
  {
    computeMLFMatrix();
  }
  
  std::vector<int> fixedVariancesIdxs;
  std::vector<double> operatorDiagonal;
  for(int i = 0; i<this->V->variances.size(); i++)
  {
    if(this->V->variances[i].fixed == true)
    {
      fixedVariancesIdxs.push_back(i);
      operatorDiagonal.push_back(0.);
    }
    else
    {
      operatorDiagonal.push_back(1.);
    }
  }
  
  if(fixedVariancesIdxs.size() != 0) //If true, compute the AI inverse assuming that the fixed variances are not present, then invert, and then add 0's rows/columns for the fixed variances indices.
  {
    Matrix *tempAI = new Matrix();
    this->AI->filterRowsAndColumns(tempAI, fixedVariancesIdxs, fixedVariancesIdxs, false);
    
    tempAI->symmetric = true;
    bool inverted = tempAI->invert();
    if(inverted == false)
    {
      this->success = false;
      delete tempAI;
      return;
    }
    
    Matrix * mOperatorTemp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->AI->nGlobRows, this->AI->nGlobCols);
    Matrix * mOperator = new Matrix();
    mOperatorTemp->fillWithConstant(0.);
    mOperatorTemp->setDiagonal(&(operatorDiagonal[0]), operatorDiagonal.size());
    std::vector<int> idxTemp;
    mOperatorTemp->filterRowsAndColumns(mOperator, fixedVariancesIdxs, idxTemp, false);
    
    mOperatorTemp->multiply(tempAI, 'N', mOperator, 'N');
    this->AI->multiply(mOperator, 'T', mOperatorTemp, 'N');
    this->AI->symmetric = true;
    
    delete mOperatorTemp;
    delete mOperator;
    delete tempAI;
  }
  else
  {
    this->AI->symmetric = true;
    bool inverted = this->AI->invert();
    if(inverted == false)
    {
      this->success = false;
      return;
    }
  }
}

void REML::computeMLFMatrix()
{
  std::vector<Matrix*> VisubVs;
  this->VisubVTraces->fillWithConstant(0.);
  std::vector<double> VisubVTracesVector;
  if( this->V->mStatus == inverseCovarianceMatrix )
  {
    for(int i = 0; i<this->V->variances.size(); i++)
    {
      Matrix * VisubV = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
      VisubV->multiply(this->V->m, 'N', this->V->computeDerivateCovariance(i), 'N');
      VisubVs.push_back(VisubV);
      VisubVTracesVector.push_back(VisubV->trace());
    }
  }
  else if(this->V->mStatus == inverseCovarianceMatrixByBlocks)
  {
    for(int i = 0; i<this->V->variances.size(); i++)
    {
      BlockMatrix VisubVtemp;
      VisubVtemp.multiply(this->V->mInBlocks, this->V->computeDerivateCovariance(i));
      Matrix * VisubV = VisubVtemp.block2distributed();
      VisubVs.push_back(VisubV);
      VisubVTracesVector.push_back(VisubV->trace());
    }
  }
  else
  {
    misc.error("Error: An internal error was happened. VinvVi matrix cannot be computed if covariance matrix is not inverted.", 0);
  }
  this->VisubVTraces->scatterMatrix(&(VisubVTracesVector[0]));
  
  this->AI->fillWithConstant(0.);
  double * gAI = new double [this->AI->nGlobRows*this->AI->nGlobCols];

  for(int i = 0; i<this->V->variances.size(); i++)
  {
    for(int j = i; j<this->V->variances.size(); j++)
    {
      double VisubViVisubVjTrace;

      VisubViVisubVjTrace = VisubVs[i]->traceOfMatrixProduct(VisubVs[j]);

      gAI[j*this->AI->nGlobRows + i] = 0.5*VisubViVisubVjTrace;
      gAI[i*this->AI->nGlobRows + j] = 0.5*VisubViVisubVjTrace;
    }
  }
  this->AI->scatterMatrix(gAI);
  
  for(int i = 0; i<this->V->variances.size(); i++)
  {
    delete VisubVs[i];
  }
  
  this->secondDerivativesMatrixType = REMLParameters::MLF;
  
  delete [] gAI;
}

void REML::computeREMLFMatrix()
{
  std::vector<Matrix*> subViPs;
  
  Matrix *Ptemp = NULL;
  if( this->P != NULL )
  {
    Ptemp = this->P;
  }
  else
  {
    Ptemp = this->PDiagonal;
  }
  for(int i = 0; i<this->V->variances.size(); i++)
  {
    Matrix * subViP = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    subViP->multiply(Ptemp, 'N', this->V->computeDerivateCovariance(i), 'N');
    subViPs.push_back(subViP);
  }
  
  this->AI->fillWithConstant(0.);
  double * gAI = new double [this->AI->nGlobRows*this->AI->nGlobCols];

  for(int i = 0; i<this->V->variances.size(); i++)
  {
    for(int j = i; j<this->V->variances.size(); j++)
    {
      double PsubViPsubVjTrace;

      PsubViPsubVjTrace = subViPs[i]->traceOfMatrixProduct(subViPs[j]);

      gAI[j*this->AI->nGlobRows + i] = 0.5*PsubViPsubVjTrace;
      gAI[i*this->AI->nGlobRows + j] = 0.5*PsubViPsubVjTrace;
    }
  }
  this->AI->scatterMatrix(gAI);
  
  for(int i = 0; i<this->V->variances.size(); i++)
  {
    delete subViPs[i];
  }
  
  this->secondDerivativesMatrixType = REMLParameters::REMLF;
  
  delete [] gAI;
}

void REML::computeAIMatrixCrossedDerivatesCorrection()
{
  if( this->V->checkMoreThanOneVariance() == false )
  {
    return;
  }
  
  Matrix * AICrossedDerivates = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->AI->nGlobRows, this->AI->nGlobCols );
  AICrossedDerivates->fillWithConstant(0.);
  
  double * gAICrossedDerivates = new double [AICrossedDerivates->nGlobRows*AICrossedDerivates->nGlobCols];
  
  for(int i = 0; i<this->V->variances.size(); i++)
  {
    for(int j = i; j<this->V->variances.size(); j++)
    {
      Matrix * Vij = this->V->computeDerivateCovariance(i, j);
      double AICrossedDerivatesij = 0;
      if(Vij != NULL)
      {
        Matrix * VijPy = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
        VijPy->multiply(Vij, 'N', this->Py, 'N');
        Matrix * ytPVijPyMatrix = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
        ytPVijPyMatrix->multiply(this->Py, 'T', VijPy, 'N');
        delete VijPy;
        double ytPVijPy = 0.;
        ytPVijPyMatrix->gatherMatrix(&ytPVijPy);
        delete ytPVijPyMatrix;
        
        double VijPTrace;
        if(this->P != NULL)
        {
          VijPTrace = this->P->traceOfMatrixProduct(Vij);
        }
        else
        {
          if(Vij->distribution != diagonalDistribution)
          {
            misc.error("Error: An internal error was happened. Expected diagonal matrix when computing AIMatrixCrossedDerivatesCorrection()", 0);
          }
          VijPTrace = this->PDiagonal->traceOfMatrixProduct(Vij);
        }
        delete Vij;
        
        AICrossedDerivatesij = 0.25*(VijPTrace - ytPVijPy);
      }

      gAICrossedDerivates[j*AICrossedDerivates->nGlobRows + i] = AICrossedDerivatesij;
      gAICrossedDerivates[i*AICrossedDerivates->nGlobRows + j] = AICrossedDerivatesij;
    }
  }
  
  AICrossedDerivates->scatterMatrix(gAICrossedDerivates);
  
  this->AI->add(AICrossedDerivates);
  
  delete [] gAICrossedDerivates;
  delete AICrossedDerivates;
  
}

void REML::computeyPsubVPy_trPsubVVector(double scale)
{
  Matrix * yPsubVPy = new Matrix(cyclicDistribution);
  yPsubVPy->multiply(this->subVPy, 'T', this->Py, 'N');
  
  double * vectorTraces = new double [this->V->variances.size()];
  Matrix * temp = new Matrix(cyclicDistribution);
  if(this->P != NULL)
  {
    for(int i=0; i<this->V->variances.size(); i++)
    {
      vectorTraces[i] = this->P->traceOfMatrixProduct(this->V->computeDerivateCovariance(i));
    }
  }
  else
  {
    for(int i=0; i<this->V->variances.size(); i++)
    {
      Matrix * VDerivate = this->V->computeDerivateCovariance(i);
      if(VDerivate->distribution != diagonalDistribution)
      {
        misc.error("Error: An internal error was happened. Expected diagonal matrix when computing computeyPsubVPy_trPsubVVector()", 0);
      }
      vectorTraces[i] = this->PDiagonal->traceOfMatrixProduct(VDerivate);
    }
  }
  this->yPsubVPy_trPsubV->scatterMatrix(vectorTraces);
  
  this->yPsubVPy_trPsubV->add(yPsubVPy, -scale, scale);

  delete [] vectorTraces;
  delete temp;  
  delete yPsubVPy;
}

void REML::computeyPsubVPy_trVisubVVector(double scale)
{
  this->yPsubVPy_trVisubV->fillWithConstant(0.);
  this->yPsubVPy_trVisubV->add(this->VisubVTraces);
  
  Matrix * yPsubVPy = new Matrix(cyclicDistribution);
  yPsubVPy->multiply(this->subVPy, 'T', this->Py, 'N');
  this->yPsubVPy_trVisubV->add(yPsubVPy, -scale, scale);
  
  delete yPsubVPy;
}

double REML::computeLogLikelihood()
{
  double previuosLogLikelihood = this->logLikelihood;
  double logLikelihoodDifference = this->logLikelihood;
  if( this->useMLinsteadOfREML == false )
  {
    this->logLikelihood = -0.5*(this->logDetV + this->logDetXtViX + this->ytPy);
  }
  else
  {
    this->logLikelihood = -0.5*(this->logDetV + this->ytPy);
  }
  logLikelihoodDifference = this->logLikelihood - logLikelihoodDifference;
  
  this->logLikelihoodRelativeDifference = fabs(logLikelihoodDifference/previuosLogLikelihood);
  
  return logLikelihoodDifference;
}

void REML::aiREMLStep()
{
  //Compute the gradients and AI or equivalent matrices.
  
  this->V->computeCovariance();
  this->success = this->V->invertCovariance(&this->logDetV, this->singlePrecisionInversion);
  if(this->success == false)
  {
    return;
  }
  
  bool clearyList = false;
  if(this->yList.size() == 0)
  {
    this->yList.push_back(this->y);
    clearyList = true;
  }
  
  Matrix * meanAI = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->V->variances.size(), this->V->variances.size());
  meanAI->fillWithConstant(0.);
  Matrix * meanGradient = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->V->variances.size(), 1);
  meanGradient->fillWithConstant(0.);
  double meanytPy = 0.;
  for(int i = 0; i < this->yList.size(); i++)
  {
    this->y = yList[i];
    
    computePMatrix();
    if(this->success == false)
    {
      return;
    }
    
    //computePyMatrix();
    computeytPy();
    meanytPy += this->ytPy/double(this->yList.size());
    
    computeSubVPyMatrix();
    
    if( this->useMLinsteadOfREML == false || i == 0 ) //Compute this only it it is the first iteration or not computing the ML F, which does not depend on this->y.
    {
      computeAIMatrix();
      if(this->success == false)
      {
        return;
      }
    }
    
    meanAI->add(this->AI, 1., 1./double(this->yList.size()));
    if( this->useMLinsteadOfREML == false )
    {
      computeyPsubVPy_trPsubVVector(0.5);
      meanGradient->add(this->yPsubVPy_trPsubV, 1., 1./double(this->yList.size()));
    }
    else
    {
      computeyPsubVPy_trVisubVVector(0.5);
      meanGradient->add(this->yPsubVPy_trVisubV, 1., 1./double(this->yList.size()));
    }
  }
  
  if(clearyList == true)
  {
    this->yList.clear();
  }
  
  delete this->AI;
  this->AI = meanAI;
  this->AI->symmetric = true;
  meanAI = NULL;
  
  this->ytPy = meanytPy;
  
  Matrix * delta = new Matrix(cyclicDistribution);
  Matrix * deltaUntransformed = NULL;
  Matrix * newVariances = new Matrix(cyclicDistribution, this->V->variances.size(), 1);
  if( this->useMLinsteadOfREML == false )
  {
    delete this->yPsubVPy_trPsubV;
    this->yPsubVPy_trPsubV = meanGradient;
    delta->multiply(this->AI, 'N', this->yPsubVPy_trPsubV, 'N');
  }
  else
  {
    delete this->yPsubVPy_trVisubV;
    this->yPsubVPy_trVisubV = meanGradient;
    delta->multiply(this->AI, 'N', this->yPsubVPy_trVisubV, 'N');
  }
  meanGradient = NULL;
  
  //Change of scale depending on options
  
  //delta->showGlobal("delta");
  std::vector<double> vv = this->V->getVectorVariances();
  this->oldVariances = vv;
  
  if(options.useLogLogisticScale == true)
  {
    Matrix * jacobianI = new Matrix(cyclicDistribution, this->V->variances.size(), this->V->variances.size());
    jacobianI->fillWithConstant(0.);
    
    std::vector<double> jacobianIv;
    for(int i = 0; i<vv.size(); i++)
    {
      if(this->V->variances[i].type == ParameterAttributes::correlation)
      {
	vv[i] = logisticInv(vv[i]);
	double temp = (1.+exp(-vv[i]));
	temp = 2.*exp(-vv[i])/(temp*temp);
	jacobianIv.push_back( 1./temp );
      }
      else
      {
	jacobianIv.push_back( 1./vv[i] );
	vv[i] = log(vv[i]);
      }
    }
    jacobianI->setDiagonal(&(jacobianIv[0]), jacobianIv.size());
    
    deltaUntransformed = new Matrix(delta);
    delta->multiply(jacobianI, 'N', deltaUntransformed, 'N');
    delete jacobianI;
  }
  
  
  //Update variances
  
  newVariances->scatterMatrix(&(vv[0]));
  bool EMPartialStep = false;
  if(this->logLikelihoodRelativeDifference > options.changeAIStepThreshold)
  {
    if( options.allowSwitchFromAItoEM == true && this->useMLinsteadOfREML == false ) //Switch to a EM step.
    {
      emPartialREMLStep(newVariances, delta, 0.5);
      this->stepModifications += "e";
      EMPartialStep = true;
    }
    else
    {
      newVariances->add(delta, 1., options.stepWeightingConstant);
      this->stepModifications += "q";
    }
  }
  else
  {
    newVariances->add(delta, 1., 1.);
  }
  
  //Change of scale depending on options
  
  newVariances->gatherMatrix(&(vv[0]));
  if( options.useLogLogisticScale == true && EMPartialStep == false )
  {
    for(int i = 0; i<vv.size(); i++)
    {
      if(this->V->variances[i].type == ParameterAttributes::correlation)
      {
	vv[i] = logistic(vv[i]);
      }
      else
      {
	vv[i] = exp(vv[i]);
      }
    }
    deltaUntransformed->gatherMatrix(&this->delta[0]);
    this->stepModifications += "l";
  }
  else
  {
    delta->gatherMatrix(&this->delta[0]);
  }
  communicator->broadcast(&(vv[0]), vv.size());
  
  //If there is a squared exponential kernel and a parameter changes over a threshold, override the step for an EM step.
  if( options.switchToEMLargeParameterChange == true && EMPartialStep == false )
  {
    for(int i = 0; i < this->V->variances.size(); i++)
    {
      if(this->V->variances[i].type == ParameterAttributes::parameter)
      {
        double temp = fabs((vv[i] - this->oldVariances[i])/this->oldVariances[i]);
        if( misc.gt( temp > options.switchToEMParameterChangeThreshold ) )
        {
          EMPartialStep = true;
        }
      }
    }
    if( EMPartialStep == true )
    {
      emPartialREMLStep(newVariances, delta, 0.5);
      
      removeCharacters(this->stepModifications, "l");
      this->stepModifications += "e";
      
      newVariances->gatherMatrix(&(vv[0]));
      communicator->broadcast(&(vv[0]), vv.size());
      
      delta->gatherMatrix(&this->delta[0]);
    }
  }
  
  //Store variances in the V matrix
  this->V->storeVectorVariances(vv);
  
  communicator->broadcast(&this->delta[0], vv.size()); 
  
  delete newVariances;
  delete delta;
  if(deltaUntransformed != NULL)
  {
    delete deltaUntransformed;
  }
}

void REML::emREMLStep()
{
  this->V->computeCovariance();
  this->success = this->V->invertCovariance(&this->logDetV, this->singlePrecisionInversion);
  if(this->success == false)
  {
    return;
  }
  
  computePMatrix();
  if(this->success == false)
  {
    return;
  }
  
  computeytPy();
  //computePyMatrix();
  
  computeSubVPyMatrix();
  
  computeyPsubVPy_trPsubVVector();
  
  double * globyPsubVPy_trPsubV = new double [this->V->variances.size()];
  yPsubVPy_trPsubV->gatherMatrix(globyPsubVPy_trPsubV);
  
  communicator->broadcast(globyPsubVPy_trPsubV, this->V->variances.size());
  
  this->oldVariances = this->V->getVectorVariances();
  
  for(int i=0; i<this->V->variances.size(); i++)
  {
    this->V->variances[i].variance = double(this->dimension)*this->V->variances[i].variance + this->V->variances[i].variance*this->V->variances[i].variance*globyPsubVPy_trPsubV[i];
    this->V->variances[i].variance = this->V->variances[i].variance/double(this->dimension);
  }
  
  std::vector<double> vv = this->V->getVectorVariances();
  communicator->broadcast(&(vv[0]), vv.size());
  this->V->storeVectorVariances(vv);
  
  
  delete [] globyPsubVPy_trPsubV;
}

void REML::emPartialREMLStep(Matrix * newVariances, Matrix * delta, double scale)
{
  double * globyPsubVPy_trPsubV = new double [this->V->variances.size()];
  yPsubVPy_trPsubV->gatherMatrix(globyPsubVPy_trPsubV);
  
  communicator->broadcast(globyPsubVPy_trPsubV, this->V->variances.size());
  
  this->oldVariances = this->V->getVectorVariances();
  
  std::vector<double> globalDelta(this->V->variances.size());
  std::vector<double> globalNewVariances(this->V->variances.size());
  for(int i=0; i<this->V->variances.size(); i++)
  {
    globalDelta[i] = this->V->variances[i].variance;
    this->V->variances[i].variance = double(this->dimension)*this->V->variances[i].variance + this->V->variances[i].variance*this->V->variances[i].variance*(globyPsubVPy_trPsubV[i]/scale);
    this->V->variances[i].variance = this->V->variances[i].variance/double(this->dimension);
    globalDelta[i] = this->V->variances[i].variance - globalDelta[i];
    globalNewVariances[i] = this->V->variances[i].variance;
  }
  
  std::vector<double> vv = this->V->getVectorVariances();
  communicator->broadcast(&(vv[0]), vv.size());
  this->V->storeVectorVariances(vv);
  
  if( communicator->mpiRoot && (globalDelta.size() != delta->nGlobRows || this->V->variances.size() != newVariances->nGlobRows) )
  {
    misc.error("Error: An internal error was happened when performing a partial EM step.", 0);
  }
  delta->scatterMatrix(&(globalDelta[0]));
  newVariances->scatterMatrix(&(globalNewVariances[0]));
  
  delete [] globyPsubVPy_trPsubV;
}

double REML::varianceDifferenceNorm()
{
  if(this->V->variances.size() != this->delta.size())
  {
    misc.error("Error: An internal error was happened in varianceDifferenceNorm().", 0);
  }
  
  double result = 0.;
  double norm = 0.;
  for(int i=0; i<this->V->variances.size(); i++)
  {
    double temp = this->V->variances[i].variance - this->oldVariances[i];
    result += temp * temp;
    norm += this->V->variances[i].variance * this->V->variances[i].variance;
  }
  
  result = result/norm;
  communicator->broadcast(&result);
  
  return result;
}

bool REML::allVariancesRelativeDifferencesLowerThan(double threshold)
{
  bool result = true;
  
  for(int i=0; i<this->V->variances.size(); i++)
  {
    double temp = (this->V->variances[i].variance - this->oldVariances[i]);
    temp /= this->oldVariances[i];
    if( fabs(temp) > threshold )
    {
      result = false;
    }
  }
  
  communicator->broadcast(&result);
  
  return result;
}

bool REML::gradientLowerThanThreshold(double threshold)
{
  bool result = true;
  
  std::vector<double> gradient;
  if( this->useMLinsteadOfREML == false )
  {
    this->yPsubVPy_trPsubV->matrixToStandardVector(gradient);
  }
  else
  {
    this->yPsubVPy_trVisubV->matrixToStandardVector(gradient);
  }
  
  if(communicator->mpiRoot == true)
  {
    for(int i = 0; i<gradient.size(); i++)
    {
      if( (fabs(gradient[i]) > threshold) && this->constrainedVariances.find(this->V->variances[i].name) == this->constrainedVariances.end() )
      {
        result = false;
      }
    }
  }
  
  communicator->broadcast(&result);
  
  return result;
}

void REML::showNonConvergenceCauses()
{
  misc.message << "The convergence status of different parameters is:" << std::endl;
  if( this->logLikelihoodConverged == false )
  {
    misc.message << "  The logL has not converged. The last difference has been: " << this->logLikelihoodDifference << std::endl;
  }
  else
  {
    misc.message << "  The logL has converged." << std::endl;
  }
  
  if( this->oldVariances.size() != 0 ) //Otherwise, REML/ML failed before performing the first step.
  {
    if( this->variancesConverged == false )
    {
      misc.message << "  The parameters have not converged. There is at least one parameter with an relative difference larger than: " << options.varianceConvergenceThreshold << std::endl;
      for(int i=0; i<this->V->variances.size(); i++)
      {
        double temp = (this->V->variances[i].variance - this->oldVariances[i]);
        temp /= this->oldVariances[i];
        misc.message << "    The " <<  this->V->variances[i].name << " relative difference is: " << fabs(temp) << std::endl;
      }
    }
    else
    {
      misc.message << "  The parameters have converged. All their relative differences are smaller than: " << options.varianceConvergenceThreshold << std::endl;
    }
  }
  
  if( this->gradientConverged == false )
  {
    misc.message << "  Gradients have not converged. There is at least one derivate with a difference larger than: " << options.gradientConvergenceThreshold << std::endl;
    std::vector<double> gradient;
    if( this->useMLinsteadOfREML == false )
    {
      this->yPsubVPy_trPsubV->matrixToStandardVector(gradient);
    }
    else
    {
      this->yPsubVPy_trVisubV->matrixToStandardVector(gradient);
    }
    
    if(communicator->mpiRoot == true)
    {
      for(int i = 0; i<gradient.size(); i++)
      {
        if( this->constrainedVariances.find(this->V->variances[i].name) == this->constrainedVariances.end() )
        {
          misc.message << "    The absolute value of the derivate by " <<  this->V->variances[i].name << " is: " << fabs(gradient[i]) << std::endl;
        }
      }
    }
  }
  else
  {
    misc.message << "  Gradients have converged. All their absolute differences are smaller than: " << options.gradientConvergenceThreshold << std::endl;
  }
  if( this->constrainedVariances.size() != 0 )
  {
    std::string temp = "";
    std::string tempsep = "";
    for (std::set<std::string>::iterator it = this->constrainedVariances.begin(); it != this->constrainedVariances.end(); ++it)
    {
      temp += tempsep + *it;
      tempsep = ", ";
    }
    misc.message << "  The parameters < " + temp + " > have been constrained and have been ignored for testing the gradient convergence." << std::endl;
  }
}

void REML::fixVariancesAndParameters()
{
  for(int i = 0; i < this->V->variances.size(); i++)
  {
    if(this->V->variances[i].type == ParameterAttributes::parameter)
    {
      this->V->variances[i].fixed = true;
    }
  }
}

void REML::unfixVariancesAndParameters()
{
  if( allVariancesRelativeDifferencesLowerThan(options.varianceConvergenceThreshold*10) == true || this->currentStep > (options.remlStepsToUnfixExpKernelParameter - 2) )
  {
    for(int i = 0; i < this->V->variances.size(); i++)
    {
      if(this->V->variances[i].type == ParameterAttributes::parameter)
      {
        if( this->V->variances[i].fixed == true )
        {
          this->logLikelihoodRelativeDifference = 1.; //Force a short AI step or EM step on the next REML step.
        }
        this->V->variances[i].fixed = false;
      }
    }
  }
}

int REML::getNVariancesFixed()
{
  int nVariancesFixed = 0;
  for(int i = 0; i < this->V->variances.size(); i++)
  {
    if(this->V->variances[i].fixed == true)
    {
      nVariancesFixed++;
    }
  }
  return nVariancesFixed;
}

void REML::computeSummary()
{
  Message message(options.outFile + ((this->useMLinsteadOfREML==false)?".reml":".ml"));
  std::vector< std::vector<double> > gAI;
  this->AI->symmetrizeTriangularMatrix();
  this->AI->matrixToStandardVector(gAI);
  
  if(options.REMLMethod == 1)
  {
    misc.error("Error: summary for EM method still not implemented.", 0);
  }
  
  if(communicator->mpiRoot)
  {
    if( this->warnings.size() != 0 )
    {
      for(std::set<std::string>::iterator it = this->warnings.begin(); it != this->warnings.end(); ++it)
      {
        message << "#  WARNING: " << *it << std::endl;
      }
      message << std::endl;
    }
    
    message << "#  Summary results:" << std::endl;
    message << "#-----------------------------\n" << std::endl;
    
    //Variances and their errors
    for(int i = 0; i< this->V->variances.size(); i++)
    {
      message << std::setw(0)<< spacetab2underscore(this->V->variances[i].name) << std::setw(options.logFieldWidth) << this->V->variances[i].variance << std::setw(options.logFieldWidth) << sqrt(gAI[i][i]) << std::endl;
    }
    
    //Group variances
    std::map<std::string, std::vector<int> > groupedGeneticVariances;
    std::map<std::string, int > groupedEnvironmentVariances;
    int nonEnvironmentVariancesExistence = false;
    for(int i = 0; i< this->V->variances.size(); i++)
    {
      if(this->V->variances[i].type != ParameterAttributes::variance )
      {
        continue;
      }
      ParameterAttributes::VarianceTypeEffect typeEffect = this->V->variances[i].typeEffect;
      if(typeEffect == ParameterAttributes::environment)
      {
        if(groupedEnvironmentVariances.count(this->V->variances[i].group) == 0)
        {
          groupedEnvironmentVariances[ this->V->variances[i].group ] = i;
        }
        else
        {
          misc.error("Error: An internal error was happened when computing the summary statistics. The structure of variances is unexpected. Two environment variances in the same group.", 0);
        }
      }
      else if(typeEffect == ParameterAttributes::genetic)
      {
        groupedGeneticVariances[ this->V->variances[i].group ].push_back(i);
        nonEnvironmentVariancesExistence = true;
      }
    }
    
    if( nonEnvironmentVariancesExistence == true )
    {
      if( groupedGeneticVariances.size() != groupedEnvironmentVariances.size() )
      {
        misc.error("Error: An internal error was happened when computing the summary statistics. The structure of variances is unexpected. Not all groups have an environment and a genetic variance.", 0);
      }
      
      //Compute heritabilities
      for(std::map<std::string, int >::iterator itEnv = groupedEnvironmentVariances.begin(); itEnv != groupedEnvironmentVariances.end(); ++itEnv)
      {
        if( groupedGeneticVariances.count(itEnv->first) != 1 )
        {
          misc.error("Error: An internal error was happened when computing the summary statistics. The structure of variances is unexpected. The group names differ between environment and genetic variances.", 0);
        }
        if( groupedGeneticVariances[itEnv->first].size() < 1 )
        {
          misc.error("Error: An internal error was happened when computing the summary statistics. The structure of variances is unexpected. The number of genetic variances in the group is 0.", 0);
        }
        message << "\n# " << itEnv->first << ":\n" << std::endl;
        
        int environmentVarianceIdx = itEnv->second;
        double environmentVariance = this->V->variances[environmentVarianceIdx].variance;
        double totalVariance = environmentVariance;
        
        double varianceTotalVariance = gAI[environmentVarianceIdx][environmentVarianceIdx];
        
        for(int idx = 0; idx < groupedGeneticVariances[itEnv->first].size(); idx++)
        {
          int geneticVarianceIdx = groupedGeneticVariances[itEnv->first][idx];
          double geneticVariance = this->V->variances[geneticVarianceIdx].variance;
          
          totalVariance += geneticVariance;
          varianceTotalVariance += gAI[geneticVarianceIdx][environmentVarianceIdx] + gAI[environmentVarianceIdx][geneticVarianceIdx];
          for(int idx2 = 0; idx2 < groupedGeneticVariances[itEnv->first].size(); idx2++)
          {
            int geneticVarianceIdx2 = groupedGeneticVariances[itEnv->first][idx2];
            varianceTotalVariance += gAI[geneticVarianceIdx][geneticVarianceIdx2];
          }
        }
        message << std::setw(0) << "Var(" << itEnv->first << ")" << std::setw(options.logFieldWidth) << totalVariance << std::setw(options.logFieldWidth) << sqrt(varianceTotalVariance) << std::endl;
        
        for(int idx = 0; idx < groupedGeneticVariances[itEnv->first].size(); idx++)
        {
          int geneticVarianceIdx = groupedGeneticVariances[itEnv->first][idx];
          double geneticVariance = this->V->variances[geneticVarianceIdx].variance;
          
          double varianceGeneticVariance = gAI[geneticVarianceIdx][geneticVarianceIdx];
          double cov = gAI[geneticVarianceIdx][environmentVarianceIdx];
          for(int idx2 = 0; idx2 < groupedGeneticVariances[itEnv->first].size(); idx2++)
          {
            int geneticVarianceIdx2 = groupedGeneticVariances[itEnv->first][idx2];
            cov += gAI[geneticVarianceIdx][geneticVarianceIdx2];
          }
          
          double h2 = geneticVariance/totalVariance;
          double varh2;
          varh2 = varianceGeneticVariance/(geneticVariance*geneticVariance);
          varh2 += varianceTotalVariance/(totalVariance*totalVariance);
          varh2 -= 2.*cov/(geneticVariance*totalVariance);
          varh2 *= h2*h2;
          
          message << std::setw(0) << spacetab2underscore(this->V->variances[geneticVarianceIdx].name) << "/Var(" << itEnv->first << ")" << std::setw(options.logFieldWidth) << h2 << std::setw(options.logFieldWidth) << sqrt(varh2) << std::endl;
        }
      }
    }
  }
  
  message << std::endl;
  
  if(communicator->mpiRoot)
  {
    int fieldWidth = 15;
    if( this->secondDerivativesMatrixType == REMLParameters::AI )
    {
      message << "#  AI Matrix inverse:" << std::endl;
    }
    else
    {
      message << "#  F Matrix inverse:" << std::endl;
    }
    message << "#-----------------------------\n" << std::endl;
    
    message << std::setw(fieldWidth*2);
    for(int i = 0; i< this->V->variances.size(); i++)
    {
      message << " " << this->V->variances[i].name << std::setw(fieldWidth);
    }
    message << std::endl;
    for(int i = 0; i< this->V->variances.size(); i++)
    {
      message << std::setw(fieldWidth) << this->V->variances[i].name;
      for(int j = 0; j< i + 1; j++)
      {
        message << std::setw(fieldWidth) << gAI[i][j];
      }
      message << std::endl;
    }
    message << std::endl;
  }
  
}

void REML::computeBLUE()
{
  misc.message << "Computing BLUEs..." << std::endl;
  
  Matrix * blue = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Matrix * temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  temp->multiply(this->ViX, 'T', this->y, 'N');
  blue->multiply(this->XtViX_i, 'N', temp, 'N');
  
  std::vector<double> globalBLUE;
  blue->matrixToStandardVector(globalBLUE);
  
  std::vector< std::vector<double> > errorBLUE;
  this->XtViX_i->matrixToStandardVector(errorBLUE);
  
  delete blue;
  delete temp;
  
  if(communicator->mpiRoot)
  {
    int nTest = 0;
    for(int idxPheno = 0; idxPheno<this->covariateNames.size(); idxPheno++)
    {
      nTest += (this->covariateNames[idxPheno].meanNames.size() + this->covariateNames[idxPheno].discreteCovarNames.size() + this->covariateNames[idxPheno].quantitativeCovarNames.size());
    }
    if( globalBLUE.size() !=  nTest)
    {
      misc.error("Error: An internal error was happened. The size of covariate names is not of same dimension of BLUE.", 0);
    }

    for(int idxPheno = 0; idxPheno<this->covariateNames.size(); idxPheno++)
    {
      Message message(options.outFile + addSuffix(".pheno" + i2s(idxPheno + 1)) + ".blue.mean");
      message << "NAME BETA STD" << std::endl;
      for(int i = 0; i < this->covariateNames[idxPheno].meanNames.size(); i++)
      {
        int idx = this->covariateNames[idxPheno].meanNames[i].second;
        message << this->covariateNames[idxPheno].meanNames[i].first << " " << globalBLUE[ idx ] << " " << sqrt(errorBLUE[ idx ][ idx ]) << std::endl;
      }
      
      message.redirect(options.outFile + addSuffix(".pheno" + i2s(idxPheno + 1)) + ".blue.discrete");
      message << "NAME BETA STD" << std::endl;
      for(int i = 0; i < this->covariateNames[idxPheno].discreteCovarNames.size(); i++)
      {
        int idx = this->covariateNames[idxPheno].discreteCovarNames[i].second;
        message << this->covariateNames[idxPheno].discreteCovarNames[i].first << " " << globalBLUE[ idx ] << " " << sqrt(errorBLUE[ idx ][ idx ]) << std::endl;
      }
      
      message.redirect(options.outFile + addSuffix(".pheno" + i2s(idxPheno + 1)) + ".blue.quantitative");
      message << "NAME BETA STD" << std::endl;
      for(int i = 0; i < this->covariateNames[idxPheno].quantitativeCovarNames.size(); i++)
      {
        int idx = this->covariateNames[idxPheno].quantitativeCovarNames[i].second;
        message << this->covariateNames[idxPheno].quantitativeCovarNames[i].first << " " << globalBLUE[ idx ] << " " << sqrt(errorBLUE[ idx ][ idx ]) << std::endl;
      }
    } //End for over different phenotypes
  } //End mpiroot if
}

void REML::computeIndividualsBLUP()
{
  for(int idxBLUP = 0; idxBLUP < individualBLUPNames.size(); idxBLUP++)
  {
    std::string name = individualBLUPNames[idxBLUP];
    
    Matrix * temp;
    
    misc.message << "Computing Individual BLUPs (" << name << ")..." << std::endl;
    
    temp = this->V->multiply(name, this->Py);
    if(this->GRMEigenVectors.count(name) != 0)
    {
      transformUsingEigenVectors(name, 'N', &(temp), this->nPhenotypes);
    }
    else
    {
      if( this->usingDiagonalKernels == true )
      {
        misc.error("Error: An internal error has happened. The analysis uses diagonal kernels, but the eigenvectors were not defined.", 0);
      }
    }
    std::vector<double> blup;
    temp->matrixToStandardVector(blup);
    
    delete temp;
    
    std::vector<double> blupErrors = computeBLUPErrors(name);
    
    if(communicator->mpiRoot)
    {
      int nTest = 0;
      for(int indPheno = 0; indPheno<this->vIndividuals.size(); indPheno++)
      {
        nTest += this->vIndividuals[indPheno].size();
      }
      
      if(blup.size() == 0)
      {
        misc.error("Error: An internal error was happened. BLUPs array is empty.", 0);
      }
      if(blup.size() != nTest)
      {
        misc.error("Error: An internal error was happened. Dimensions differ when computing individual BLUPs.", 0);
      }
      if(blupErrors.size() != 0 && blupErrors.size() != blup.size())
      {
        misc.error("Error: An internal error was happened. Dimensions differ between BLUP errors and BLUPs.", 0);
      }
      
      int shift = 0;
      for(int indPheno = 0; indPheno<this->vIndividuals.size(); indPheno++)
      {
        Message message(options.outFile + "." + spacetab2underscore(name) + addSuffix(".pheno" + i2s(indPheno + 1)) + ".blup.indiv");
        message << "FID" << " " << "IID" << " " << "BLUP";
        if( blupErrors.size() != 0 )
        {
          message << " " << "STD";
        }
        message << std::endl;
        for(int i = 0; i < vIndividuals[indPheno].size(); i++)
        {
          message << this->vIndividuals[indPheno][i].familyID << " " << this->vIndividuals[indPheno][i].individualID << " " << blup[ shift + i ];
          if( blupErrors.size() != 0 )
          {
            message << " " << blupErrors[ shift + i ];
          }
          message << std::endl;
        }
        shift += vIndividuals[indPheno].size();
      } //End for storing different traits
    } //End mpiRoot if
  } //End for over different subCovariances
}

std::vector<double> REML::computeBLUPErrors(std::string name)
{
  if( options.computeIndividualsBLUPErrors == false )
  {
    return std::vector<double>();
  }
  if(this->P == NULL)
  {
    return std::vector<double>();
  }
  if(this->GRMEigenVectors.count(name) != 0)
  {
    return std::vector<double>(); //To implement
  }
  
  Matrix * cov = this->V->getSubCovariance(name);
  
  Matrix * temp = new Matrix();
  temp->diagonalOfABAt(cov, this->P);
  
  std::vector<double> covPcovDiag = temp->diagonal();
  std::vector<double> covDiag = cov->diagonal();
  
  delete cov;
  delete temp;
  
  std::vector<double> errors;
  if( communicator->mpiRoot == true )
  {
    errors = std::vector<double>(covDiag.size(), 0);
    for(int i = 0; i<covDiag.size(); i++)
    {
      //errors[i] = sqrt(covDiag[i] - covPcovDiag[i]);
      //misc.message << "Testing blup error..." << std::endl;
      errors[i] = sqrt(covPcovDiag[i]);
    }
  }
  return errors;
}

void REML::computeSNPsBLUP()
{
  for(std::map<std::string, std::vector<std::string> >::iterator it = this->SNPsBLUPGenotypeFiles.begin(); it != this->SNPsBLUPGenotypeFiles.end(); ++it)
  {
    std::string name = it->first;
    std::vector<std::string> files = it->second;
    
    for(int i = 0; i < files.size(); i++)
    {
      if( this->usingDiagonalKernels == true )
      {
        misc.error("Error: An internal error has happened. Computing SNP BLUPs from genotype files when using diagonal kernels, is not implemented, yet.", 0);
      }
      
      misc.message << "Using genotypes specified in [ " + files[i] + " ] for computing the SNP BLUPs from the GRM < " + name + " >" << std::endl;
      
      Genotype* genotypes = new Genotype(files[i]);
      
      std::string baseOutput = options.outFile;
      options.outFile += "." + getFileName(files[i]);
      
      computeSNPsBLUP(name, genotypes);
      
      options.outFile = baseOutput;
      
      delete genotypes;
    }
  }
  
  for(std::map<std::string, Genotype *>::iterator it = this->SNPsBLUPGenotypes.begin(); it != this->SNPsBLUPGenotypes.end(); ++it)
  {
    std::string name = it->first;
    Genotype* genotypes = it->second;

    computeSNPsBLUP(name, genotypes);
    
  } //End different subcovariances iteration.
}

void REML::computeSNPsBLUP(std::string name, Genotype* genotypes)
{
  misc.message << "Computing SNP BLUPs (" << name << ")..." << std::endl;

  if(genotypes == NULL)
  {
    misc.error("Error: An internal error was happened when computing SNP BLUPs. No genotypes are loaded.", 0);
  }
  
  if( this->mSNPIds.count(name) == 0 )
  {
    misc.error("Error: An internal error was happened. " + name + " is not defined in the mSNPIds on " + this->sREMLML + ".", 0);
  }
  std::multiset<std::string> test1(this->mSNPIds[name].begin(), this->mSNPIds[name].end());
  std::multiset<std::string> test2(genotypes->SNPIds.begin(), genotypes->SNPIds.end());
  if( test1 != test2 )
  {
    misc.message << "WARNING: SNPs in genotype data for computing SNP BLUPs differ from SNPs used for computing GRM." << std::endl;
  }
  std::vector<std::string> SNPIdsToKeep = orderVectorAsTemplate(genotypes->SNPIds, this->mSNPIds[name]);
  if( SNPIdsToKeep != genotypes->SNPIds )
  {
    if( options.removeNonOvarlapingBLUPSNPs == true )
    {
      misc.message << "WARNING: SNPs in genotype data for computing SNP BLUPs include SNPs not present in the GRM < " + name + " >. These will be filtered. You can use the option --blup-no-filter-snps to avoid this."  << std::endl;
    }
    else
    {
      misc.message << "WARNING: SNPs in genotype data for computing SNP BLUPs include SNPs not present in the GRM < " + name + " >. You are using the option --blup-no-filter-snps and those will not be filtered."  << std::endl;
      SNPIdsToKeep = genotypes->SNPIds;
    }
  }
  
  if( misc.gt( SNPIdsToKeep.size() == 0 ) )
  {
    misc.message << "WARNING: The overlapping between SNPs in genotype data for computing SNP BLUPs and the GRM < " + name + " > is empty. Skipping..."  << std::endl;
    return;
  }
  int nSNPIdsToKeep = SNPIdsToKeep.size();
  communicator->broadcast( &nSNPIdsToKeep );

  if(communicator->mpiRoot == true )
  {
    if( this->nPhenotypes != this->vIndividuals.size() )
    {
      misc.error("Error: An internal error was happened. The number of phenotypes disagrees with the stores samples size.", 0);
    }
  }
  
  //Compute the product of the blocks genotype*PyBlock. Where PyBlock es the Py values for each phenotype
  int PyRowShift = 0;
  std::vector< Matrix* > snpBLUPsSubMatrices;
  Matrix * snpBLUPsTemp;
  Matrix *nNonMissings = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nSNPIdsToKeep, 1);
  nNonMissings->fillWithConstant(0.);
  int nTotalIndividuals = 0;
  std::vector<SNP> filteredGenotypesSNPs;
  std::vector<double> BLUPErrors;
  for(int phenoIdx = 0; phenoIdx < this->nPhenotypes; phenoIdx++)
  {
    
    std::vector<std::string> tempIndividualIds;
    for(int i = 0; i < this->vIndividuals[phenoIdx].size(); i++)
    {
      tempIndividualIds.push_back(this->vIndividuals[phenoIdx][i].familyID + "@" + this->vIndividuals[phenoIdx][i].individualID);
    }
    
    genotypes->normalizeGenotypes();
    Genotype * filteredGenotypes = NULL;
    bool newFilteredGenotypes = false;
    if( misc.gt(SNPIdsToKeep != genotypes->SNPIds || tempIndividualIds != genotypes->individualIds) ) //To save memory, create new file only if filtering is required.
    {
      filteredGenotypes =  new Genotype();
      genotypes->filterSNPsAndIndividuals(SNPIdsToKeep, tempIndividualIds, true, filteredGenotypes);
      newFilteredGenotypes = true;
    }
    else
    {
      filteredGenotypes = genotypes;
    }
    if(filteredGenotypes->individualIds != tempIndividualIds)
    {
      misc.error("Error: An error was happened when computing SNP BLUPs. The individuals in the genotype data are not the same or in the same order than in GRMs used for " + this->sREMLML + " analysis.", 0);
    }
    
    nTotalIndividuals += filteredGenotypes->nIndividuals;
    
    Matrix *tempColumnWithOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, filteredGenotypes->nIndividuals, 1);
    Matrix *tempNNonMissings = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    tempColumnWithOnes->fillWithConstant(1.);
    tempNNonMissings->multiply(filteredGenotypes->missings, 'N', tempColumnWithOnes, 'N');
    nNonMissings->add(tempNNonMissings);
    delete tempColumnWithOnes;
    delete tempNNonMissings;
    
    //
    Matrix * PyBlock = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, filteredGenotypes->nIndividuals, 1);
    PyBlock->fillWithConstant(0.);
    subMatrix smDest(PyBlock);
    subMatrix smSrc(PyRowShift, 0, filteredGenotypes->nIndividuals, 1);
    PyBlock->add(this->Py, 1., 1., smDest, smSrc);
    PyRowShift += filteredGenotypes->nIndividuals;
    //
    
    snpBLUPsTemp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    snpBLUPsTemp->multiply(filteredGenotypes->genotypes, 'N', PyBlock, 'N');
    snpBLUPsSubMatrices.push_back(snpBLUPsTemp);

    if( phenoIdx == 0 )
    {
      filteredGenotypesSNPs = filteredGenotypes->SNPs;
    }
    
    if( options.computeIndividualsBLUPErrors == true && this->nPhenotypes == 1 && this->P != NULL && this->GRMEigenVectors.size() == 0) //Currently, the errors can only be computed in this restricted conditions, for testing pourposes.
    {
      Matrix * temp = new Matrix();
      temp->diagonalOfABAt(filteredGenotypes->genotypes, this->P);
      
      std::vector<double> missingsVecTemp;
      nNonMissings->matrixToStandardVector(missingsVecTemp); //This is assuming that there is only one phenotype (this->nPhenotypes == 1), otherwise this would be wrong.
      BLUPErrors = temp->diagonal();
      if(missingsVecTemp.size() != BLUPErrors.size())
      {
        misc.error("Error: An internal error was happened. Unexpected dimensions for SNP BLUP errors.", 0);
      }
      for(int ibe = 0; ibe < BLUPErrors.size(); ibe++)
      {
        BLUPErrors[ibe] *= (nTotalIndividuals*nTotalIndividuals)/(missingsVecTemp[ibe]*missingsVecTemp[ibe]);
      }
      delete temp;
    }
    
    if(newFilteredGenotypes == true)
    {
      delete filteredGenotypes;
    }
    delete PyBlock;
  }
  
  if( nTotalIndividuals != this->dimension || nTotalIndividuals != this->V->dimension )
  {
    misc.error("Error: An internal error was happened when computing SNP BLUPs. The total number of individuals is different than the covariance dimensions.", 0);
  }
  
  //Compute the covariance matrix for the current GRM.
  std::vector< std::vector<double> > covarianceMatrix( this->nPhenotypes, std::vector<double>(this->nPhenotypes, 0.) );
  for(int phenoIdx1 = 0; phenoIdx1 < this->nPhenotypes; phenoIdx1++)
  {
    for(int phenoIdx2 = phenoIdx1; phenoIdx2 < this->nPhenotypes; phenoIdx2++)
    {
      int elementIdx = -1;
      if( phenoIdx1 == phenoIdx2)
      {
        elementIdx = this->V->getElement( name + "_" + i2s(phenoIdx1 + 1) );
      }
      else
      {
        elementIdx = this->V->getElement( name + "_" + i2s(phenoIdx1 + 1) + "_" + i2s(phenoIdx2 + 1) );
      }
      if( elementIdx == -1 )
      {
        misc.error("Error: An internal error was happened when computing SNP BLUPs. The element does not exist.", 0);
      }
      double elementVariance = this->V->computeElementMultiplyingVariance( this->V->elements[ elementIdx ] );
      communicator->broadcast(&elementVariance);
      
      covarianceMatrix[phenoIdx1][phenoIdx2] = elementVariance;
      covarianceMatrix[phenoIdx2][phenoIdx1] = elementVariance;
    }
  }
  
  //Compute the total SNPBLUPs combining the products genotype*PyBlock with the computed variances and store the results.
  for(int phenoIdx1 = 0; phenoIdx1 < this->nPhenotypes; phenoIdx1++)
  {
    Matrix * snpBLUPsTotal = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nSNPIdsToKeep, 1);
    snpBLUPsTotal->fillWithConstant(0.);
    for(int phenoIdx2 = 0; phenoIdx2 < this->nPhenotypes; phenoIdx2++)
    {
      snpBLUPsTotal->add(snpBLUPsSubMatrices[phenoIdx2], 1., covarianceMatrix[phenoIdx1][phenoIdx2]);
    }
    int tempNGRMSNPs = this->mSNPIds[name].size();
    communicator->broadcast(&tempNGRMSNPs);
    snpBLUPsTotal->elementWiseDivision(nNonMissings, double(nTotalIndividuals)/double(tempNGRMSNPs));
    std::vector<double> snpBLUPs;
    snpBLUPsTotal->matrixToStandardVector(snpBLUPs);
    delete snpBLUPsTotal;
    
    if(communicator->mpiRoot)
    {
      if( snpBLUPs.size() != filteredGenotypesSNPs.size() )
      {
        misc.error("Error: An internal error was happened when storing SNP BLUPs.", 0);
      }
      
      Message message(options.outFile + "." + spacetab2underscore(name) + addSuffix(".pheno" + i2s(phenoIdx1 + 1)) + ".blup.snps");
      message << "SNP ALLELE BLUP STDEV MEAN NBLUP" << ((BLUPErrors.size() != 0)?" BLUPE":"") << std::endl;
      for(int i = 0; i < snpBLUPs.size(); i++)
      {
        message << filteredGenotypesSNPs[i].name << " " << filteredGenotypesSNPs[i].allele2;
        message << " " << std::setprecision(14) << snpBLUPs[i];
        message << " " << std::setprecision(14) << filteredGenotypesSNPs[i].standardDev;
        message << " " << std::setprecision(14) << 2.*filteredGenotypesSNPs[i].p2;
        message << " " << std::setprecision(14) << snpBLUPs[i]/filteredGenotypesSNPs[i].standardDev;
        if(BLUPErrors.size() != 0)
        {
          message << " " << std::setprecision(14) << sqrt(BLUPErrors[i]*covarianceMatrix[phenoIdx1][phenoIdx1]*covarianceMatrix[phenoIdx1][phenoIdx1]/(tempNGRMSNPs*tempNGRMSNPs));
        }
        message << std::endl;
      }
    }
  }
    
  delete nNonMissings;
  for(int i =0; i<snpBLUPsSubMatrices.size(); i++)
  {
    delete snpBLUPsSubMatrices[i];
  }
  snpBLUPsSubMatrices.clear();
}

void REML::verboseLogOutput1()
{
  misc.message << "Variances on final step (logLikelihood: " << std::setprecision(10) << this->logLikelihood << "):" << std::endl;
  for(int i = 0; i< this->V->variances.size(); i++)
  {
    misc.message << std::setw(0)<< spacetab2underscore(this->V->variances[i].name) << std::setw(options.logFieldWidth) << std::setprecision(10) << this->V->variances[i].variance << std::endl;
    
  }
  
  std::vector< std::vector<double> > gAI;
  int fieldWidth = 15;
  this->AI->symmetrizeTriangularMatrix();
  this->AI->matrixToStandardVector(gAI);
  
  if(communicator->mpiRoot)
  {
    misc.message << std::setw(fieldWidth*2);
    for(int i = 0; i< this->V->variances.size(); i++)
    {
      misc.message << " " << this->V->variances[i].name << std::setw(fieldWidth);
    }
    misc.message << std::endl;
    for(int i = 0; i< this->V->variances.size(); i++)
    {
      misc.message << std::setw(fieldWidth) << this->V->variances[i].name;
      for(int j = 0; j< i + 1; j++)
      {
        misc.message << std::setw(fieldWidth) << gAI[i][j];
      }
      misc.message << std::endl;
    }
    misc.message << std::endl;
  }
}

void REML::experimentalEpistasisPrediction()
{
  if( options.computeEpistasisVariance == false || options.epistatitcPredictionGRMFile == "" )
  {
    return;
  }
  
  for(int idxBLUP = 0; idxBLUP < individualBLUPNames.size(); idxBLUP++)
  {
    std::string name = individualBLUPNames[idxBLUP];
    
    misc.message << "Computing Epistasis Predictions..." << std::endl;
    
    Kernel * grm = new Kernel(options.epistatitcPredictionGRMFile);
    Kernel * epistaticGRMBase = new Kernel(grm, kernelEpistaticGRM, NULL);
    delete grm;
    
    std::vector<std::string> tempIndividualIds;
    for(int i = 0; i < this->vIndividuals[0].size(); i++)
    {
      tempIndividualIds.push_back(this->vIndividuals[0][i].familyID + "@" + this->vIndividuals[0][i].individualID);
    }
    epistaticGRMBase->filterIndividualsAsymmetric(epistaticGRMBase->individualIds, tempIndividualIds, false);
    
    Matrix * predictions = new Matrix();
    predictions->multiply(epistaticGRMBase->kernel, 'N', this->Py, 'N');
    
    std::vector<double> blup;
    predictions->matrixToStandardVector(blup);
    
    delete predictions;
    
    if(communicator->mpiRoot)
    {
      
      if(blup.size() == 0)
      {
        misc.error("Error: An internal error was happened. BLUPs array is empty.", 0);
      }
      
      Message message(options.outFile + "." + spacetab2underscore(name) + addSuffix(".pheno" + i2s(0 + 1)) + ".epi.blup.indiv");
      message << "FID" << " " << "IID" << " " << "BLUP";

      message << std::endl;
      for(int i = 0; i < epistaticGRMBase->individualIdsRows.size(); i++)
      {
        message << epistaticGRMBase->individualIdsRows[i] << " " << blup[ i ] << std::endl;
      }
    } //End mpiRoot if
    
    delete epistaticGRMBase;
  } //End for over different subCovariances
}

double REML::logistic(double x)
{
  return ( 2./(1.+exp(-x)) ) - 1;
}

double REML::logisticInv(double x)
{
  return log( (x+1.) /( 2.-(x+1.) ));
}