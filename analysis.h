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

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "options.h"
#include "misc.h"
#include "genotype.h"
#include "kernel.h"
#include "phenotype.h"
#include "reml.h"
#include "simulatephenotype.h"
#include "predictphenotype.h"
#include "pca.h"
#include "gwas.h"
#include "igwas.h"
#include "results.h"
#include "covariancematrix.h"
#include "auxiliar.h"

/*struct REMLResults
{
  std::vector<Variance> variances;
  double logLikelihood;
};*/

class Analysis
{
public:
  //REMLResults previousREMLResults;
  
  Analysis();
  ~Analysis();
  
  void makeGRM();
  void makeGRMAndStoreMostRelated();
  void makeREML();
  void makeMultivarREML();
  void makeSimulatePhenotype();
  void makePredictPhenotype();
  void makePCA();
  void makeGWAS();
  void makeRecursiveGWAS();
  void makeComputeAccuracyBySNP();
  void makeEffectsAnalysis();
  void makeGLMMAnalysis();
  void makeIGWASAnalysis();
  void makePredictCovarPhenotype();
  void makeMultiplePhenotypeGWAS();
  void makeMultiplePhenotypeResiduals();
  void makeFilterLabeledMatrix();
  void makeSNPStats();
};

#endif
