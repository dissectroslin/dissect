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

#include "pcagentemp.h"
#include "kernel.h"
#include "matrix.h"
#include "communicator.h"
#include "misc.h"
#include "options.h"
#include "labeledmatrix.h"

#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>

PCAGenTemp::PCAGenTemp(LabeledMatrix * lm)
{
  if(lm->getRowLabels() != lm->getColLabels())
  {
    misc.error("Error: An error has happened. The labels on rows and columns differ when performing a PCA.", 0);
  }
  
  this->labels = lm->getRowLabels();
  
  Matrix * eigenValues = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Matrix * eigenVectors = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Matrix * selectedEigenVectors = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  
  misc.setGetElapsedTime("PCA");
  misc.message << "Starting PCA analysis..." << std::endl;
  lm->getMatrix()->eigenDecomposition(eigenValues, eigenVectors);
  misc.message << "PCA analysis has ended after "  << misc.setGetElapsedTime("PCA") << ". Storing results." << std::endl;
  
  if(options.nEigenValues > eigenVectors->nGlobRows)
  {
    misc.message << "Warning: Due to matrix dimensions, only " << eigenVectors->nGlobRows << " eigenvectors can be stored." << std::endl;
    options.nEigenValues = eigenVectors->nGlobRows;
  }

  int rowsKeep[eigenVectors->nGlobRows];
  int colsKeep[options.nEigenValues];
  for(int i=0; i<eigenVectors->nGlobRows; i++)
  {
    rowsKeep[i] = i;
  }
  for(int i= 0; i<options.nEigenValues; i++)
  {
    int idx = (eigenVectors->nGlobCols - options.nEigenValues) + i;
    colsKeep[i] = idx;
  }

  eigenVectors->filterRowsAndColumns(selectedEigenVectors, rowsKeep, eigenVectors->nGlobRows, colsKeep, options.nEigenValues);
  
  this->globalEigenValues = eigenValues->diagonal();
  selectedEigenVectors->matrixToStandardVector(this->globalEigenVectors);
  
  delete eigenValues;
  delete eigenVectors;
  delete selectedEigenVectors;
  
  if(communicator->mpiRoot)
  {
    if(this->labels.size() != this->globalEigenVectors.size())
    {
      misc.error("Error: An internal error was happened when computing the pca.", 0);
    }
    
    Message message(options.outFile + ".pca.eigenvalues");
    for(int i = (int(globalEigenValues.size()) - 1); i>=0; i--)
    {
      message << this->globalEigenValues[i] << std::endl;
    }
    
    message.redirect(options.outFile + ".pca.eigenvectors");
    for(int i = 0; i<this->globalEigenVectors.size(); i++)
    {
      message << this->labels[i];
      for(int j = (options.nEigenValues - 1); j>=0; j--)
      {
        message << " " << this->globalEigenVectors[i][j];
      }
      message << std::endl;
    }
  }
}
