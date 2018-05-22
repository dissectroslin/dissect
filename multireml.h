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

#ifndef MULTIREML_H
#define MULTIREML_H

#include "matrix.h"
#include "genotype.h"
#include "covariate.h"
#include "reml.h"
#include "covariancematrix.h"
#include "kernel.h"

#include <vector>
#include <set>
#include <utility>
#include <fstream>

class MultiREML
{
public:
  REML * reml;
  
  bool writeResults;
  
  REMLResults subsampleREMLResults;
  
  MultiREML(bool argWriteResults = true);
  ~MultiREML();
  
  /**
   * Compute REML
   */
  void compute();
  
  /**
   * Compute regional REML
   * 
   * Compute regional REML. Each computation uses a regional GRM plus a global GRM with the regional part substracted.
   */
  void computeRegional();
  
  /**
  * Compute groups REML
  * 
  * Compute groups REML. One computation with SNPs divided into different GRMs and different variances all at same time.
  */
  void computeMultipleGroups();
  
  /**
   * Compare the results of two REMLs and store a p-value.
   */
  void compareREMLs(double baseLogLikelihood, double baseNVariances);
};

#endif
