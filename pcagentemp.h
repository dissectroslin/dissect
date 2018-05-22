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

#ifndef PCAGENTEMP_H
#define PCAGENTEMP_H

#include "matrix.h"
#include "labeledmatrix.h"

#include <string>
#include <vector>
#include <map>

class Genotype;
class Kernel;

/**
 * Class for computing de PCA of a matrix;
 * 
 */
class PCAGenTemp
{
  std::vector<std::string> labels;
  
  std::vector<double> globalEigenValues;
  std::vector< std::vector<double> > globalEigenVectors;

public:
  PCAGenTemp(LabeledMatrix * lm);
};

#endif
