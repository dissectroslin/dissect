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

#ifndef RANGE_H
#define RANGE_H

#include <string>

class Range
{
private:
  bool downThresholdActive;
  bool upThresholdActive;
  
  int iDownThreshold;
  int iUpThreshold;
  
  double dDownThreshold;
  double dUpThreshold;

public:
  Range();
  
  Range(bool flag, int upThreshold);
  Range(int downThreshold, bool flag);
  Range(int downThreshold, int upThreshold);
  
  Range(bool flag, double upThreshold);
  Range(double downThreshold, bool flag);
  Range(double downThreshold, double upThreshold);
  
  bool checkRange(int x);
  bool checkRange(double x);
  
  std::string explainRange(int x);
  std::string explainRange(double x);
  
  ~Range();
};

#endif
