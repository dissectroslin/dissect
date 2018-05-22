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

MPResiduals::MPResiduals()
{
}

MPResiduals::~MPResiduals()
{
}

LabeledMatrix * MPResiduals::computeResiduals(std::string grmFile)
{
  misc.setGetElapsedTime("residualsComputation");
  misc.message << "Computing the residuals..." << std::endl;
  
  //Load grm
  Kernel * currentGRMBase = new Kernel(grmFile);
  currentGRMBase->normalize();
  currentGRMBase->name = "GRM";

  
  //Load phenotypes and covariates
  std::vector<std::string> th;
  th.push_back("FID");
  th.push_back("IID");
  LabeledMatrix * phenotypes = new LabeledMatrix(options.phenotypesFile, th);
  std::vector<std::string> prefilterCovInds = intersectionStringVectors(2, &currentGRMBase->individualIds, &phenotypes->rowLabels);
  std::set<std::string> setPrefilterCovInds(prefilterCovInds.begin(), prefilterCovInds.end());
  Covariate * covariate = new Covariate(options.covarsFile, options.qCovarsFile, phenotypes->rowLabels, true, setPrefilterCovInds);
  
  
  //Get shared individuals and filter
  std::vector<std::string> commonIndividuals = intersectionStringVectors(3, &phenotypes->rowLabels, &covariate->individualIds, &currentGRMBase->individualIds);
  std::vector<std::string> commonIndividualsInGRMOrder = orderVectorAsTemplate(currentGRMBase->individualIds, commonIndividuals);
  if( currentGRMBase->diagonalized == true )
  {
    if( commonIndividualsInGRMOrder != currentGRMBase->individualIds )
    {
      misc.error("Error: When using diagonal GRMs, the individuals in phenotype, covars and GRM files must be all the same without missing phenotypes. Aborting.", 0);
    }
  }
  else
  {
    currentGRMBase->filterIndividuals(commonIndividualsInGRMOrder, false);
    currentGRMBase->diagonalizeKernel();
  }
  phenotypes->filterRowsAndCols(commonIndividualsInGRMOrder, phenotypes->colLabels);
  covariate->filterIndividuals(commonIndividualsInGRMOrder);

  //Correct phenotypes, and covariates using GRM eigenvectors
  Matrix * tempPhenos = new Matrix();
  tempPhenos->multiply(currentGRMBase->eigenVectors, 'T', phenotypes->matrix, 'N');
  delete phenotypes->matrix;
  phenotypes->matrix = tempPhenos;
  
  Matrix * tempCovars = new Matrix();
  tempCovars->multiply(currentGRMBase->eigenVectors, 'T', covariate->covariates, 'N');
  delete covariate->covariates;
  covariate->covariates = tempCovars;
  
  //Get residuals
  LabeledMatrix * residuals = new LabeledMatrix(phenotypes);
  residuals->matrix->fillWithConstant(0.);
  
  int nPhenotypes = phenotypes->matrix->nGlobCols;
  std::vector<std::string> convergedResiduals;
  std::vector<std::string> nonConvergedResiduals;
  for(int ipheno = 0; ipheno<nPhenotypes; ipheno++)
  {
    //Prepare, run REML, and compute covariance matrix
    Matrix * remly = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, phenotypes->matrix->nGlobRows, 1);
    remly->fillWithConstant(0.);
    remly->add(phenotypes->matrix, 0., 1., subMatrix(0, 0, remly->nGlobRows, 1), subMatrix(0, ipheno, phenotypes->matrix->nGlobRows, 1));
    
    Matrix * remlX = new Matrix(covariate->covariates);
    Matrix * remlKernel = new Matrix(currentGRMBase->getNormalizedKernel());
    
    std::vector<Matrix*> kernels;
    kernels.push_back(remlKernel);
    
    std::vector<KernelType> kernelTypes;
    kernelTypes.push_back(kernelCovarianceMatrix);
    
    std::vector<double> weights;
    weights.push_back(1.);
    
    REML reml;
    bool prepared = reml.prepare(remly, remlX, kernels, kernelTypes, options.initialh2Trait,  weights);
    
    kernels.clear();
    if( prepared == true )
    {
      reml.computeREML();
    }
    else
    {
      if(communicator->mpiRoot == true)
      {
        misc.message << "An error happened preparing the REML analysis, residuals for " + phenotypes->colLabels[ipheno] + " cannot be computed." << std::endl;
        nonConvergedResiduals.push_back(phenotypes->colLabels[ipheno]);
      }
      continue;
    }
    if( reml.success == false )
    {
      if(communicator->mpiRoot == true)
      {
        misc.message << "REML did not converge, residuals for " + phenotypes->colLabels[ipheno] + " cannot be computed.." << std::endl;
        nonConvergedResiduals.push_back(phenotypes->colLabels[ipheno]);
      }
      continue;
    }
    
    Matrix * residual = reml.V->multiply("E", reml.Py);
    if(communicator->mpiRoot == true)
    {
      convergedResiduals.push_back( phenotypes->colLabels[ipheno] );
    }
    residuals->insertCol(residual, ipheno);
    delete residual;
  }
  if(misc.gt(residuals->colLabels != convergedResiduals))
  {
    residuals->filterRowsAndCols(residuals->rowLabels, convergedResiduals);
  }
  
  delete phenotypes;
  delete covariate;

  if(communicator->mpiRoot && nonConvergedResiduals.size() != 0)
  {
    misc.message << "There are " << nonConvergedResiduals.size() << " phenotypes for which the residuals cannot be computed. They are stored in file [ " << (options.outFile + ".gwas.phenos.unfitted") << " ]." << std::endl;
    Message message(options.outFile + ".gwas.phenos.unfitted");
    for(int i = 0; i < nonConvergedResiduals.size(); i++)
    {
      message << nonConvergedResiduals[i] << std::endl;
    }
  }

  Matrix * tempPhenos2 = new Matrix();
  tempPhenos2->multiply(currentGRMBase->eigenVectors, 'N', residuals->matrix, 'N');
  delete residuals->matrix;
  residuals->matrix = tempPhenos2;
  

  delete currentGRMBase;
  currentGRMBase = NULL;
  
  misc.message << "Residuals computed after " << misc.setGetElapsedTime("residualsComputation", true) << "." << std::endl;
  
  if(misc.gt(residuals->colLabels.size() == 0))
  {
    misc.error("Error: The residual estimation failed for all phenotypes.", 0);
  }
  
  return residuals;
}