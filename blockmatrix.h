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

#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H

#include "communicator.h"
#include "global.h"
#include "matrix.h"

#include <vector>


class BlockMatrix
{
public:
  std::vector< std::vector<Matrix*> > m;        ///<Array with all blocks. It is assumed that all blocks in same row have same number of rows, and each block in same column hace same number of columns.

  int nGlobRows;
  int nGlobCols;
  
  int nBlockRows;       ///<Number of blocks in the rows
  int nBlockCols;       ///<Number of blocks in the columns
 
  /**
   * Creates an empty block matrix
   */
  BlockMatrix();
  
  /**
   * Creates a block matrix from matrices in mi
   * 
   * \param mi array of matrices that will create the block matrix. They must be un row-major order. Matrices will be copied, so after creating the block matrix, it must be cleared.
   */
  BlockMatrix(std::vector< std::vector<Matrix*> > & mi);
  
  /**
   * Creates a block matrix from a copy of another block matrix.
   * 
   * \param srcBlockMatrix Block Matrix to be copied.
   */
  BlockMatrix(BlockMatrix & srcBlockMatrix);
  
  ~BlockMatrix();
  
  /**
   * Delete all the infor in this block matrix transforming it to an empty block matrix.
   */
  void clear();
  
  /**
   * Adds a row of blocks after the last row of blocks i nthis block matrix
   */
  void addBlockRow(std::vector<Matrix*> & blockRow);
  
  /**
   * Multiply two block matrices
   * 
   * Multiply two block matrices and store the result in this matrix. Note that previous values in this matrix (if any) will be deleted.
   * 
   * \param m1 left block matrix to be multiplied.
   * \param m2 right block matrix to be multiplied.
   * \param factor A factor that be multiplied to matrix multiplication.
   */
  void multiply(BlockMatrix & m1, BlockMatrix & m2, double factor = 1.);
  
  /**
   * Multiply a block matrix with cyclic distributed matrix.
   * 
   * Multiply a block matrix with cyclic distributed matrix and store the result in this matrix. Note that previous values in this matrix (if any) will be deleted.
   * 
   * \param m1 left block matrix to be multiplied.
   * \param m2 right matrix to be multiplied.
   * \param factor A factor that be multiplied to matrix multiplication.
   */
  void multiply(BlockMatrix & m1, Matrix * m2, double factor = 1.);
  
  /**
   * Add another block matrix to this block matrix
   * 
   * Both block matrices must have same number of blocks with same dimensions.
   * 
   * \param m1 block matrix to be added to this.
   * \param thisFactor A factor multiplying this matrix before adition.
   * \param m1Factor A factor multiplying m1 matrix before adition.
   */  
  void add(BlockMatrix & m1, double thisFactor = 1., double m1Factor = 1.);
  
  /**
   * Invert this block matrix.
   * 
   * The inversion require that all matrices in the block must be symmetric, diagonal and positive definite.
   * 
   * \param logDeterminant If different than NULL, the log of the determinant of the block matrix will be returned.
   */
  bool invert(double * logDeterminant = NULL);
  
  /**
   * Get a cyclic distributed matrix from a block matrix
   * 
   * \return a pointer to the cyclic distributed matrix.
   */
  Matrix* block2distributed();
  
  void showGlobal();
};

#endif
