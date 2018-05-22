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

#ifndef ACCURACYBYSNP_H
#define ACCURACYBYSNP_H

#include "matrix.h"
#include "phenotype.h"
#include "predictphenotype.h"

class Genotype;

/*struct SNPEffect
{
  std::string name;
  std::string refAllele;
  double effect;
  double neffect;
  double stdev;
  double mean;
};*/

/**
 * Class for computing the prediction as a function of SNPs removed
 * 
 * 
 * 
 */
class AccuracyBySNP
{
  Genotype * genotype;
  PredictPhenotype * predictPhenotype;
  Phenotype *realPhenotypes;
  
public:
  AccuracyBySNP();
  ~AccuracyBySNP();
  
  void computeAccuracies();
  double computeAccuracy(PredictPhenotype * pp, Matrix * rp);
  void standardizeVector(Matrix * vector);
  double accuracyFilteringAt(std::vector<double> SNPaccuracies, double stdScaleThreshold);
};

#endif
