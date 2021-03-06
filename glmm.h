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

#ifndef GLMM_H
#define GLMM_H

#include "matrix.h"
#include "glm.h"
#include "reml.h"

#include <vector>
#include <set>
#include <utility>
#include <fstream>

class Matrix;


class GLMM
{
public:
  GLMParameters::LinkFunction linkFunction;
  
  REML * reml;
  
  bool success;
  int nIterations;
  
  Matrix * y;					///<Phenotypes matrix.
  Matrix * X;					///<Covariates matrix.
  Matrix * betas;                               ///<Parameters.
  //Matrix * randomEffects;                       ///<random effects.
  
  Matrix * previousSample;                      ///<
  
  long seed;
  
  /**
   * 
   * \param irelm Expects a prepared REML.
   */
  GLMM(REML * ireml, GLMParameters::LinkFunction iLinkFunction = GLMParameters::logit);
  ~GLMM();
  
  void deleteIntermediateMatrices();
  
  Matrix * MHSampling(Matrix * u, bool getStartingPoint = false);
  Matrix * computeProbabilities(Matrix * randomEffects);
  
  bool fit();
  bool iteration(int nSamples);
  /*
  
  
  
  void computeLogLikelihoodGradient();
  bool computeLogLikelihoodSecondDerivatives();

  bool newtonIteration();
  bool fit();
  bool allParametersRelativeDifferencesLowerThan(std::vector<double> & oldBetas, double threshold);
  */
};

#endif
