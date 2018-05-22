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

#ifndef SINGLEREML_H
#define SINGLEREML_H

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

class SingleREML
{
public:
  REML * reml;
  
  bool writeResults;
  
  REMLResults subsampleREMLResults;
  
  SingleREML(bool argWriteResults = true);
  ~SingleREML();
  
  /**
   * Initialize REML parameters
   * 
   * \param grm The grm that will be used for preparing REML it will be deleted in this function..
   * \param genotypes Genotypes that will be used for computing SNPs BLUPs.
   */
  //void prepare(Kernel * grm, Genotype *genotypes = NULL);
  
  /**
   * Initialize REML parameters assuming multiple GRMs
   * 
   * \param grms Pointers to the grms that will be used for preparing REML. All these grms will be deleted and this vector cleared.
   * \param names Names that will be used for grms.
   * \param weights weights for the initial variances of each grm in grms. If empty, 1./grms.size() will be used for each grm.
   */
  //void prepare(std::vector<Kernel*> & grms, std::vector<std::string> names, std::vector<double> weights = std::vector<double>());
  
  /**
   * Compute REML initial values.
   * 
   * For computing REML initial values, this function can perform several REML analysis with subsets of full data or get initial values from a file.
   * Which option will be used depend on the options.
   */
  void precomputeInitialValues(int type, Kernel *srcGRM);
  
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
   * Compute just a RAW REML
   * 
   * Compute just a RAW REML without computing any BLUPs or BLUEs. ithou loading GRM from options. Without genotypes file. Just use srcGRM.
   * 
   * \param srcGRM The grm that will be used for computing REML.
   */
  void computeRaw(Kernel *srcGRM);

};

#endif
