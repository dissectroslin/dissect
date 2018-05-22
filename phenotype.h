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

#ifndef PHENOTYPE_H
#define PHENOTYPE_H

#include <string>
#include <vector>
#include <map>
#include "matrix.h"

class Phenotype
{
public:
  Matrix *phenotypes;
  
  int nIndividuals;
  
  int nPhenotypesInFile;
  
  std::vector<std::string> individualIds;
  std::map<std::string, int> individualIdsIdx;
  
  std::vector<std::string> missings;
  
  

  
  Phenotype(DistributionType dist, std::string f, int column);
  ~Phenotype();
  
  void getNumberPhenotypesInFile(std::string f);
  
  void readPhenotype(DistributionType dist, std::string f, int column);
  
  void filterIndividuals(std::vector<std::string> keepIndividualIds);
  
  double computePhenotypeVariance();
  
  void replaceIndividualIds(std::vector<std::string> & newIndividualIds);
  
  void printPhenotype();
};

#endif
