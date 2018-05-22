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

#ifndef LABELEDMATRIX_H
#define LABELEDMATRIX_H

#include <vector>
#include <string>
#include "matrix.h"

class LabeledMatrix
{
public:
  Matrix * matrix;
  std::vector<std::string> rowLabels;
  std::vector<std::string> colLabels;

  ////////////////////////////////////////////////////////////////

  LabeledMatrix();
  LabeledMatrix(int nRows, int nCols);
  LabeledMatrix(LabeledMatrix * lm);
  
  /**
   * Create a new labeled matrix from a file
   * 
   * \param lColumns if size == 0, load from a previously stored labeled matrix. Otherwise, load from a plain text file, using the first lColumns.size() as row labels and first row as column labels. If more then one column is used for labeling, the labels are concatenated using "@".
   */
  LabeledMatrix(std::string fn, std::vector<std::string> lColumns = std::vector<std::string>());
  ~LabeledMatrix();
  
  void clear();
  
  std::vector<std::string> getRowLabels(); //No longer required
  std::vector<std::string> getColLabels(); //No longer required
  Matrix * getMatrix();                    //No longer required
  
  void setRowLabels(std::vector<std::string> labels);
  void setColLabels(std::vector<std::string> labels);
  
  void insertRow(Matrix * row, int idx);
  void insertCol(Matrix * col, int idx);
  
  void filterRowsAndCols(std::vector<std::string> & rowLabelsKeep, std::vector<std::string> & colLabelsKeep);
  
  void appendVertically(LabeledMatrix * lm);
  void appendHorizontally(LabeledMatrix * lm);
  
  void save(std::string fn);
  void load(std::string fn);
  
  void loadRaw(std::string fn, std::vector<std::string> lColumns = std::vector<std::string>());
  
  void show(int padding = 5);
};

#endif