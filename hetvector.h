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

#ifndef HETVECTOR_H
#define HETVECTOR_H

#include "communicator.h"
#include "global.h"
#include "message.h"
#include "genotype.h"
#include "covariate.h"

#include <string>
#include <map>

class HetVector
{
public:

  std::map<std::string, Genotype *> genotypes;
  std::map<std::string, Covariate *> covariates;
  
  std::set<std::string> elements;                       ///<All the keys in genotypes or covariates. There cannot be repeated keys.
  
  HetVector();
  ~HetVector();

  void addElement(std::string name, Genotype * genotype);
  void addElement(std::string name, Covariate * covariate);
  void deleteElement(std::string name);
  Matrix * getElementMatrix(std::string name, std::vector<std::string> individualIds);

};

#endif
