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

#ifndef SIMULATEPHENOTYPE_H
#define SIMULATEPHENOTYPE_H

#include "matrix.h"

#include <string>
#include <vector>
#include <map>

class Genotype;

/**
 * Class for simulating the phenotypes of a population
 * 
 * Given genotypes and SNP effects, the phenotypes in a population will be simulated.
 * 
 */
class SimulatePhenotype
{
  Genotype * genotype;
  
  Genotype * adjustEffectsGenotype;
  
  Matrix * effectSizes;

public:
  std::vector<std::string> causalSNPIds;                ///< The names of the causal SNPs. The names are in the same order than in genotypes.
  std::map<std::string, double> causalSNPEffects;
  
  long seed;
  
  SimulatePhenotype(Genotype * srcGenotype, std::string fname, Genotype * argAdjustEffectsGenotype = NULL);
  ~SimulatePhenotype();
  
  void readCausalSNPs(std::string f);
  void simulatePhenotypes();
  
};

#endif
