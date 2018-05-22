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

#include "glmm.h"
#include "glm.h"
#include "matrix.h"
#include "options.h"
#include "communicator.h"
#include "misc.h"
#include "global.h"
#include "reml.h"
#include "covariate.h"
#include "genotype.h"
#include "phenotype.h"
#include "auxiliar.h"
#include "message.h"

#include <fstream>
#include <iomanip>
#include <cmath>

GLMM::GLMM(REML * ireml, GLMParameters::LinkFunction iLinkFunction)
{
  this->linkFunction = iLinkFunction;
  
  this->reml = ireml;
  this->reml->useMLinsteadOfREML = true;
  this->reml->V->computeCovariance();
  bool fisuccess = this->reml->V->invertCovariance(NULL, true);
  if(fisuccess == false)
  {
    misc.error("Error: An error has happened before starting iterations. Sorry, the covariance matrix is not invertible.", 0);
  }
  
  this->y = this->reml->y;
  this->reml->y = NULL;
  this->X = new Matrix(this->reml->X);
  this->betas = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->X->nGlobCols, 1);
  this->betas->fillWithConstant(0.);
  //this->randomEffects = NULL;
  
  this->previousSample = NULL;
  
  this->success = false;
  
  this->seed = options.randomSeed;
  
  this->nIterations = -1;
}

GLMM::~GLMM()
{
  if(this->y != NULL)
  {
    delete this->y;
  }
  if(this->X != NULL)
  {
    delete this->X;
  }
  if(this->betas != NULL)
  {
    delete this->betas;
  }
  if(this->reml != NULL)
  {
    delete this->reml;
  }
  /*if(this->randomEffects != NULL)
  {
    delete this->randomEffects;
  }*/
  
  deleteIntermediateMatrices();
}

void GLMM::deleteIntermediateMatrices()
{
  if(this->previousSample != NULL)
  {
    delete this->previousSample;
    this->previousSample = NULL;
  }
}

Matrix * GLMM::MHSampling(Matrix * u, bool getStartingPoint)
{
  if(this->reml->V->mStatus != inverseCovarianceMatrix)
  {
    misc.error("Error: An internal error has happened. The covariance matrix is not inverted when performing MH sampling.", 0);
  }
  std::vector<double> gDiagonalInv = this->reml->V->m->diagonal();
  std::vector<double> gSample(gDiagonalInv.size(), 0.);
  std::vector<double> gUniformSample(gDiagonalInv.size(), 0.);
  if( communicator->mpiRoot == true )
  {
    for(int i = 0; i < gDiagonalInv.size(); i++)
    {
      gDiagonalInv[i] = 1./gDiagonalInv[i];
      gSample[i] = box_muller(0., sqrt(fabs(gDiagonalInv[i])), &seed);
      gUniformSample[i]  = unif_rand_dbl(&seed); //ATTENTION: Check that this is generating the expected distribution between 0 and 1. I checked: test100().
    }
  }
  Matrix * sample = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->reml->V->m->nGlobRows, 1);
  sample->scatterMatrix(&(gSample[0]));
  if( getStartingPoint == true )
  {
    return sample;
  }
  Matrix * diagonalInv = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->reml->V->m->nGlobRows, 1);
  diagonalInv->scatterMatrix(&(gDiagonalInv[0]));
  Matrix * uniformSample = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->reml->V->m->nGlobRows, 1);
  uniformSample->scatterMatrix(&(gUniformSample[0]));
  
  if(this->reml->V->m->distribution != diagonalDistribution)
  {
    Matrix * covarianceWithoutDiagonal = new Matrix(this->reml->V->m);
    double * zeros;
    if( communicator->mpiRoot == true )
    {
      zeros = new double[ covarianceWithoutDiagonal->nGlobRows ];
      for(int i = 0; i < covarianceWithoutDiagonal->nGlobRows; i++)
      {
        zeros[i] = 0.;
      }
    }
    covarianceWithoutDiagonal->setDiagonal(zeros, covarianceWithoutDiagonal->nGlobRows);
    if( communicator->mpiRoot == true )
    {
      delete [] zeros;
    }
    Matrix * means = new Matrix();
    means->multiply(covarianceWithoutDiagonal, 'N', u, 'N', -1.);
    means->elementWiseMultiplication(diagonalInv);
    sample->add(means, 1., 1.);
    delete covarianceWithoutDiagonal;
    delete means;
  }
  
  //the following code depends on distributions and a lot of things. it has to be passed to matrix.cpp
  Matrix * newB = computeProbabilities(sample);
  Matrix * oldB = computeProbabilities(u);
  this->y->checkMatrixStructure(newB);
  this->y->checkMatrixStructure(oldB);
  #pragma omp parallel for
  for(int c = 0; c<this->y->nCols; c++)
  {
    for(int r = 0; r<this->y->nRows; r++)
    {
      if(this->y->m[c*this->y->nRows + r] < 0.5)
      {
        newB->m[c*newB->nRows + r] = 1. - newB->m[c*newB->nRows + r];
        oldB->m[c*oldB->nRows + r] = 1. - oldB->m[c*oldB->nRows + r];
      }
    }
  }
  Matrix * ratio = new Matrix(newB);
  ratio->elementWiseDivision(oldB);
  
  
  ratio->checkMatrixStructure(uniformSample);
  ratio->checkMatrixStructure(sample);
  ratio->checkMatrixStructure(u);
  #pragma omp parallel for
  for(int c = 0; c<ratio->nCols; c++)
  {
    for(int r = 0; r<ratio->nRows; r++)
    {
      if(ratio->m[c*ratio->nRows + r] < uniformSample->m[c*ratio->nRows + r])
      {
        sample->m[c*ratio->nRows + r] = u->m[c*ratio->nRows + r];
      }
    }
  }
  
  delete diagonalInv;
  delete ratio;
  delete oldB;
  delete newB;
  delete uniformSample;
  
  return sample;
}

Matrix * GLMM::computeProbabilities(Matrix * randomEffects)
{
  Matrix * probabilities;
  
  if( randomEffects != NULL )
  {
    if( randomEffects->nGlobRows != this->X->nGlobRows )
    {
      misc.error("Error: An internal error was happened. Invalid dimensions for incidence matrix and/or random effects matrix.", 0);
    }
    
    probabilities = new Matrix(randomEffects);
    
    Matrix * XBetas = new Matrix();
    XBetas->multiply(this->X, 'N', this->betas, 'N');
    
    for( int c = 0; c < probabilities->nGlobCols; c++ )
    {
      probabilities->add(XBetas, 1., 1., subMatrix(0, c, probabilities->nGlobRows, 1), subMatrix(XBetas));
    }
    
    delete XBetas;
  }
  else
  {
    probabilities = new Matrix();
    probabilities->multiply(this->X, 'N', this->betas, 'N');
  }
  
  if( this->linkFunction == GLMParameters::logit )
  {
    probabilities->applyInverseLogistic();
  }
  else
  {
    misc.error("Error: An internal error was happened. Unsupported link function.", 0);
  }
  
  return probabilities;
}

bool GLMM::fit()
{
  for(int i = 0; i < 20; i++)
  {
    iteration(300);
  }
}

bool GLMM::iteration(int nSamples)
{
  misc.setGetElapsedTime("samples");
  misc.message << "Sampling..." << std::endl;
  
  Matrix * randomEffects = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->y->nGlobRows, nSamples);
  randomEffects->fillWithConstant(0.);
  
  std::vector<Matrix *> samples;
  
  if(this->previousSample == NULL)
  {
    this->previousSample = MHSampling(NULL, true);
  }
  
  for(int i = 0; i<nSamples; i++)
  {
    samples.push_back(this->previousSample);
    randomEffects->add(this->previousSample, 1., 1., subMatrix(0, i, this->previousSample->nGlobRows, 1), subMatrix(this->previousSample));
    this->previousSample = MHSampling(this->previousSample);
  }
  delete this->previousSample;
  this->previousSample = NULL;
  
  misc.message << nSamples << " samples generated after " << misc.setGetElapsedTime("samples", true) << "." << std::endl;
  
  //Fit GLM
  misc.message << "Fitting GLM..." << std::endl;
  GLM * glm = new GLM(&y, &X, &this->betas, &randomEffects, this->linkFunction, false);
  bool success = glm->fit();
  if(success == false)
  {
    misc.message << "Sorry, GLM fitting failed. Stoping GLMM iterations." << std::endl;
    return false;
  }
  misc.message << "GLM have been fitted after " << glm->nIterations << " iterations." << std::endl;
  this->y = glm->y;
  glm->y = NULL;
  this->X = glm->X;
  glm->X = NULL;
  this->betas = glm->betas;
  glm->betas = NULL;
  
  delete glm;
  
  //Fit ML
  this->reml->setyList(samples);
  this->reml->y = samples[0];
  this->reml->singlePrecisionInversion = true;
  this->reml->computeREML();
  this->reml->deleteyList();
  
  this->reml->V->setInitialVariancesWithCurrentValues();
  
  this->betas->showGlobal("Betas");
  
  return true;
}