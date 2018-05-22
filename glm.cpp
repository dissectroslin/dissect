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

GLM::GLM(Matrix **iy, Matrix **iX, Matrix **iBetas, Matrix **iRandomEffects, GLMParameters::LinkFunction iLinkFunction, bool copyMatrices)
{
  this->linkFunction = iLinkFunction;
  
  this->y = NULL;
  this->X = NULL;
  this->betas = NULL;
  this->randomEffects = NULL;
  if( copyMatrices == true )
  {
    this->y = new Matrix(*iy);
    this->X = new Matrix(*iX);
    if( iBetas != NULL )
    {
      this->betas = new Matrix(*iBetas);
    }
    else
    {
      this->betas = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->X->nGlobCols, 1);
      this->betas->fillWithConstant(0.);
    }
    if( iRandomEffects != NULL )
    {
      this->randomEffects = new Matrix(*iRandomEffects);
    }
  }
  else
  {
    this->y = *iy;
    this->X = *iX;
    if( iBetas != NULL )
    {
      this->betas = *iBetas;
    }
    else
    {
      this->betas = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->X->nGlobCols, 1);
      this->betas->fillWithConstant(0.);
    }
    if( iRandomEffects != NULL )
    {
      this->randomEffects = *iRandomEffects;
    }
    *iy = NULL;
    *iX = NULL;
    *iBetas = NULL;
    *iRandomEffects = NULL;
  }
  
  this->probabilities = NULL;
  this->meanProbabilities = NULL;
  this->gradient = NULL;
  this->secondDerivativesInv = NULL;
  
  this->success = false;
  
  this->nIterations = -1;
}

GLM::~GLM()
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
  if(this->randomEffects != NULL)
  {
    delete this->randomEffects;
  }
  
  deleteIntermediateMatrices();
}

void GLM::deleteIntermediateMatrices()
{
  if(this->probabilities != NULL)
  {
    delete this->probabilities;
    this->probabilities = NULL;
  }
  if(this->meanProbabilities != NULL)
  {
    delete this->meanProbabilities;
    this->meanProbabilities = NULL;
  }
  if(this->gradient != NULL)
  {
    delete this->gradient;
    this->gradient = NULL;
  }
  if(this->secondDerivativesInv != NULL)
  {
    delete this->secondDerivativesInv;
    this->secondDerivativesInv = NULL;
  }
}


void GLM::computeProbabilities()
{
  if(this->probabilities != NULL)
  {
    delete this->probabilities;
    this->probabilities = NULL;
  }
  if(this->meanProbabilities != NULL)
  {
    delete this->meanProbabilities;
    this->meanProbabilities = NULL;
  }
  
  if( this->randomEffects != NULL )
  {
    if( this->randomEffects->nGlobRows != this->X->nGlobRows )
    {
      misc.error("Error: An internal error was happened. Invalid dimensions for incidence matrix and/or random effects matrix.", 0);
    }
    
    this->probabilities = new Matrix(this->randomEffects);
    
    Matrix * XBetas = new Matrix();
    XBetas->multiply(this->X, 'N', this->betas, 'N');
    
    for( int c = 0; c < this->probabilities->nGlobCols; c++ )
    {
      this->probabilities->add(XBetas, 1., 1., subMatrix(0, c, this->probabilities->nGlobRows, 1), subMatrix(XBetas));
    }
    
    delete XBetas;
  }
  else
  {
    this->probabilities = new Matrix();
    this->probabilities->multiply(this->X, 'N', this->betas, 'N');
  }
  
  if( this->linkFunction == GLMParameters::logit )
  {
    this->probabilities->applyInverseLogistic();
    if( this->probabilities->nGlobCols != 1 )
    {
      this->meanProbabilities = new Matrix();
      Matrix * normVector = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->probabilities->nGlobCols, 1);
      normVector->fillWithConstant(1./double(this->probabilities->nGlobCols));
      this->meanProbabilities->multiply(this->probabilities, 'N', normVector, 'N');
      
      delete normVector;
    }
    else
    {
      this->meanProbabilities = new Matrix(this->probabilities);
    }
  }
  else
  {
    misc.error("Error: An internal error was happened. Unsupported link function.", 0);
  }
}

void GLM::computeLogLikelihoodGradient()
{
  if( this->linkFunction == GLMParameters::logit )
  {
    Matrix * temp = new Matrix(this->y);
    temp->add(this->meanProbabilities, 1., -1.);
    
    this->gradient->multiply(this->X, 'T', temp, 'N');
    
    delete temp;
  }
  else
  {
    misc.error("Error: An internal error was happened. Insuported link function when computing the gradient.", 0);
  }
}

bool GLM::computeLogLikelihoodSecondDerivatives()
{
  bool sdsuccess = false;
  if( this->linkFunction == GLMParameters::logit )
  {
    Matrix * temp = new Matrix(this->probabilities); //Here it will be computed pi*(1-pi);
    Matrix * tempOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->probabilities->nGlobRows, this->probabilities->nGlobCols);
    tempOnes->fillWithConstant(1.);
    temp->add(tempOnes, -1., 1.);
    delete tempOnes;
    
    temp->elementWiseMultiplication(this->probabilities);
    
    double * globalDiagonal;
    if(communicator->mpiRoot)
    {
      globalDiagonal = new double [temp->nGlobRows];
    }
    
    if( temp->nGlobCols != 1 )
    {
      Matrix * mean = new Matrix();
      Matrix * normVector = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, temp->nGlobCols, 1);
      normVector->fillWithConstant(1./double(this->probabilities->nGlobCols));
      mean->multiply(temp, 'N', normVector, 'N');
      
      mean->gatherMatrix(globalDiagonal);
      
      delete normVector;
      delete mean;
    }
    else
    {
      temp->gatherMatrix(globalDiagonal);
    }
    
    Matrix * diagonal = new Matrix(diagonalDistribution, this->probabilities->nGlobRows, this->probabilities->nGlobRows);
    diagonal->setDiagonal(globalDiagonal, temp->nGlobRows);
    
    delete temp;
    if(communicator->mpiRoot)
    {
      delete [] globalDiagonal;
    }
    
    temp = new Matrix();
    temp->multiply(this->X, 'T', diagonal, 'N', -1.);
    this->secondDerivativesInv->multiply(temp, 'N', this->X, 'N');
    this->secondDerivativesInv->symmetric = true;
    sdsuccess = this->secondDerivativesInv->invert();
  }
  else
  {
    misc.error("Error: An internal error was happened. Unsupported link function when computing the gradient.", 0);
  }
  
  return sdsuccess;
}

bool GLM::newtonIteration()
{
  computeProbabilities();
  computeLogLikelihoodGradient();
  bool sdsuccess = computeLogLikelihoodSecondDerivatives();
  if(sdsuccess == false)
  {
    return false;
  }
  
  Matrix * delta = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  delta->multiply(this->secondDerivativesInv, 'N', this->gradient, 'N');
  
  this->betas->add(delta, 1., -1.);
  
  return true;
}

bool GLM::fit()
{
  this->probabilities = new Matrix();
  this->meanProbabilities = new Matrix();
  this->gradient = new Matrix();
  this->secondDerivativesInv = new Matrix();

  bool iterate = true;
  this->nIterations = 0;
  this->success = true;
  
  std::vector<double> oldBetas;
  while(iterate)
  {
    this->betas->matrixToStandardVector(oldBetas);
    this->success = newtonIteration();
    
    if( this->nIterations != 0 && allParametersRelativeDifferencesLowerThan(oldBetas, 0.001) == true )
    {
      iterate = false;
    }
    
    if( iterate == true && this->nIterations > 30 )
    {
      this->success = false;
      iterate = false;
    }
    
    this->nIterations++;
  }
  
  deleteIntermediateMatrices();
  
  return this->success;
}

bool GLM::allParametersRelativeDifferencesLowerThan(std::vector<double> & oldBetas, double threshold)
{
  bool result = true;
  
  std::vector<double> betas;
  this->betas->matrixToStandardVector(betas);
  
  if(communicator->mpiRoot == true)
  {
    for(int i=0; i<betas.size(); i++)
    {
      double temp = (betas[i] - oldBetas[i]);
      temp /= oldBetas[i];
      if( fabs(temp) > threshold )
      {
        result = false;
      }
    }
  }
  
  communicator->broadcast(&result);
  
  return result;
}