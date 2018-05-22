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

#ifndef PREDICTPHENOTYPE_H
#define PREDICTPHENOTYPE_H

#include "matrix.h"
#include "genotype.h"

#include <string>
#include <vector>
#include <map>

class Genotype;

struct SNPEffect
{
  std::string name;
  std::string refAllele;
  double effect;
  double neffect;
  double stdev;
  double mean;
};

/**
 * Class for predicting the phenotypes from SNP effects
 * 
 * Given genotypes and SNP effects, the phenotypes in a population will be predicted.
 * 
 */
class PredictPhenotype
{
public:
  Genotype * genotype;
  
  std::map<std::string, SNPEffect> SNPEffects;
  
  std::vector<std::string> effectSNPIds;
  std::vector<std::string> individualIds;
  std::vector<Individual> individuals;
  
  std::vector<double> globalPhenotypes;
  std::vector<double> globalPhenotypeShift;

  std::vector<double> globalCovariateEffects;
  
  double backgroundSNPEffectShift;
  
  Matrix * matrixEffects;
  Matrix * matrixShift;
  
  bool effectSNPsEmpty;
  
  PredictPhenotype(Genotype * srcGenotype, std::string fname, bool errorWhenNoSNPEffects = true);
  ~PredictPhenotype();
  
  void readSNPsEffects(std::string f, bool errorWhenNoSNPEffects);
  SNPEffect loadREMLEffect(std::string f, std::string line);
  SNPEffect loadOldREMLEffect(std::string f, std::string line);
  SNPEffect loadGWASEffect(std::string f, std::string line);
  
  void predictPhenotypes();
  void addMoreEffects(PredictPhenotype *previousPredictions);
  
  void storePredictions();
};

#endif
