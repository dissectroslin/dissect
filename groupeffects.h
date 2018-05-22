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

#ifndef GROUPEFFECTS_H
#define GROUPEFFECTS_H

#include <math.h>
#include <string>
#include <map>

#include "labeledmatrix.h"

class GroupAttributes
{
public:
  std::string name;
  std::string chrom;
  double minpos;
  double maxpos;

  GroupAttributes() {};
  ~GroupAttributes() {};
  
  std::pair<bool, double> getDistance( GroupAttributes & group )
  {
    if( this->chrom == group.chrom )
    {
      if( (this->minpos > group.minpos && this->minpos < group.maxpos) ||
          (this->maxpos > group.minpos && this->maxpos < group.maxpos) ||
          (group.minpos > this->minpos && group.minpos < this->maxpos) ||
          (group.maxpos > this->minpos && group.maxpos < this->maxpos)
      )
      {
        return std::pair<bool, double>(true, 0);
      }
      else
      {
        double distance = std::min( fabs(this->maxpos - group.minpos),  fabs(this->minpos - group.maxpos));
        return std::pair<bool, double>(true, distance);
      }
    }
    else
    {
      return std::pair<bool, double>(false, -1);
    }
  }
};

class GroupEffects
{
public:
  LabeledMatrix * effects;
  
  GroupEffects(std::string fn);
  GroupEffects(std::vector<std::string> fns, RowColumn rowcolumn);
  ~GroupEffects();
  
  LabeledMatrix * computeCorrelations(RowColumn rowcolumn);
  LabeledMatrix * computeCovariances(RowColumn rowcolumn);
  
  void getPairsHighlyCorrelated(LabeledMatrix * correlations, double threshold, std::vector<int> & rowIdxs, std::vector<int> & colIdxs);
  std::map<std::string, GroupAttributes> getGroupPositions(std::string fn);
  void filterCorrelatedGroups(RowColumn rowcolumn, double threshold, std::string fn);
};

#endif