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

#include "matrix.h"
#include "options.h"
#include "communicator.h"
#include "global.h"
#include "misc.h"
#include "auxiliar.h"

#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm> 
#include <vector>
#include <limits>

#include <omp.h>
#include <string.h>


subMatrix::subMatrix(int ir, int ic, int nr, int nc)
{
  if(ir < 0 || ic < 0 || nr <= 0 || nc <= 0) {misc.error("Error: An internal error was happened", 0);}
  this->active=true; this->ir=ir; this->ic=ic; this->nr=nr; this->nc=nc;
};

subMatrix::subMatrix(Matrix * m)
{
  this->active=true; this->ir=0; this->ic=0; this->nr=m->nGlobRows; this->nc=m->nGlobCols;
}

subMatrix::subMatrix() {this->active=false;}

Matrix::Matrix(DistributionType dist, int ngr, int ngc, int nbr, int nbc)
{
  defaultMatrix(dist);
  this->initParameters(ngr, ngc, nbr, nbc, dist);
}

Matrix::Matrix(DistributionType dist)
{
  this->nGlobCols = 0;
  this->nGlobRows = 0;
  defaultMatrix(dist);
}

Matrix::Matrix()
{
  this->nGlobCols = 0;
  this->nGlobRows = 0;
  defaultMatrix(MATRIX_DEFAULT_DISTRIBUTION);
}

Matrix::Matrix(Matrix *srcMatrix)
{
  if(srcMatrix == NULL)
  {
    misc.error("Error: An internal error was happened. Trying to copy an inexistent matrix.", 0);
  }
  defaultMatrix(srcMatrix->distribution);
  duplicateMatrix(srcMatrix);
}

void Matrix::defaultMatrix(DistributionType dist)
{
  this->m = NULL;
  this->v = NULL;
  this->mSinglePrecision = NULL;
  this->distribution = dist;
  this->symmetric = false;
  this->uplo = 'B';
  this->vector = false;
}

Matrix::~Matrix()
{
  if(this->m != NULL)
  {
    delete [] this->m;
    this->m = NULL;
    misc.estimateMaxMemory(-(double(this->nGlobRows)*double(this->nGlobCols)*8.));
  }
  
  if(this->v != NULL)
  {
    delete [] this->v;
    this->v = NULL;
  }
  
  if(this->mSinglePrecision != NULL)
  {
    delete [] this->mSinglePrecision;
    this->mSinglePrecision = NULL;
    misc.estimateMaxMemory(-(double(this->nGlobRows)*double(this->nGlobCols)*4.));
  }
}


void Matrix::duplicateMatrix(Matrix *srcMatrix)
{
  this->initParameters(srcMatrix->nGlobRows, srcMatrix->nGlobCols, srcMatrix->nBlockRows, srcMatrix->nBlockCols, srcMatrix->distribution);
  this->distribution = srcMatrix->distribution;
  this->symmetric = srcMatrix->symmetric;
  this->uplo = srcMatrix->uplo;
  this->vector = srcMatrix->vector;
  
  #pragma omp parallel for
  for (int i = 0; i < (this->nRows*this->nCols); i++)
  {
    this->m[i] = srcMatrix->m[i];
  }
}

void Matrix::checkMatrixStructure(Matrix * m1)
{
  if(
    this->distribution != m1->distribution ||
    this->nGlobRows != m1->nGlobRows ||
    this->nGlobCols != m1->nGlobCols ||
    this->nBlockRows != m1->nBlockRows ||
    this->nBlockCols != m1->nBlockCols ||
    this->nRows != m1->nRows ||
    this->nCols != m1->nCols
  )
  {
    misc.error("Error: An internal error was happened. Tried to perform not allowed operations on matrices with different structure.", 0);
  }
}

void Matrix::scatterBlock(double *mBlockGlobal, int r, int c, int blockRowLength, int sourceProcessRow, int sourceProcessCol)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Scattering a block is not allowed on a non-distributed matrix.", 0);
  }
  
  if((r%this->nBlockRows) != 0 || (c%this->nBlockCols) != 0)
  {
    misc.error("Error: An internal error was happened", 0);
  }
  
  int sendProcessRow = 0, sendProcessColumn = 0, submatrixRowPosition = 0, submatrixColPosition = 0;

  //sendProcessRow and sendProcessColumn defines in which process this block will be stored
  sendProcessRow = r/this->nBlockRows;
  sendProcessRow = sendProcessRow % communicator->nProcRows;
  sendProcessColumn = c/this->nBlockCols;
  sendProcessColumn = sendProcessColumn % communicator->nProcCols;
  
  //submatrixRowPosition and submatrixColPosition defines where in the submatrix this block will be stored.
  submatrixRowPosition = r/this->nBlockRows;
  submatrixRowPosition = submatrixRowPosition/communicator->nProcRows;
  submatrixRowPosition = submatrixRowPosition*this->nBlockRows;
  submatrixColPosition = c/this->nBlockCols;
  submatrixColPosition = submatrixColPosition/communicator->nProcCols;
  submatrixColPosition = submatrixColPosition*this->nBlockCols;

  
  int nr = this->nBlockRows;
  if (this->nGlobRows-r < this->nBlockRows)
  {
    nr = this->nGlobRows-r;
  }
  int nc = this->nBlockCols;
  if (this->nGlobCols-c < this->nBlockCols)
  {
    nc = this->nGlobCols-c;
  }
  
  if (communicator->myRow == sourceProcessRow && communicator->myCol == sourceProcessCol) {
    // Send the submatrix to process (sendProcessRow, sendProcessColumn)
    Cdgesd2d(communicator->context, nr, nc, mBlockGlobal, blockRowLength, sendProcessRow, sendProcessColumn);
  }
  if (communicator->myRow == sendProcessRow && communicator->myCol == sendProcessColumn) {
    // Receive the same data. The leading dimension of the local matrix is this->nrows.
    Cdgerv2d(communicator->context, nr, nc, this->m+this->nRows*submatrixColPosition+submatrixRowPosition, this->nRows, sourceProcessRow, sourceProcessCol);
  }
  
}

void Matrix::gatherBlock(double *mBlockGlobal, int r, int c, int blockRowLength)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Gathering a block is not allowed on a non-distributed matrix.", 0);
  }
  
  if((r%this->nBlockRows) != 0 || (c%this->nBlockCols) != 0)
  {
    misc.error("Error: An internal error was happened", 0);
  }
  
  int sendProcessRow = 0, sendProcessColumn = 0, submatrixRowPosition = 0, submatrixColPosition = 0;
  
  //sendProcessRow and sendProcessColumn defines in which process this block will be stored
  sendProcessRow = r/this->nBlockRows;
  sendProcessRow = sendProcessRow % communicator->nProcRows;
  sendProcessColumn = c/this->nBlockCols;
  sendProcessColumn = sendProcessColumn % communicator->nProcCols;
  
  //submatrixRowPosition and submatrixColPosition defines where in the submatrix this block will be stored.
  submatrixRowPosition = r/this->nBlockRows;
  submatrixRowPosition = submatrixRowPosition/communicator->nProcRows;
  submatrixRowPosition = submatrixRowPosition*this->nBlockRows;
  submatrixColPosition = c/this->nBlockCols;
  submatrixColPosition = submatrixColPosition/communicator->nProcCols;
  submatrixColPosition = submatrixColPosition*this->nBlockCols;
  
  
  int nr = this->nBlockRows;
  if (this->nGlobRows-r < this->nBlockRows)
  {
    nr = this->nGlobRows-r;
  }
  int nc = this->nBlockCols;
  if (this->nGlobCols-c < this->nBlockCols)
  {
    nc = this->nGlobCols-c;
  }
  
  if (communicator->myRow == sendProcessRow && communicator->myCol == sendProcessColumn) {
    // Send a nr-by-nc submatrix to process (sendProcessRow, sendProcessColumn)
    Cdgesd2d(communicator->context, nr, nc, this->m+this->nRows*submatrixColPosition+submatrixRowPosition, this->nRows, 0, 0);
    submatrixColPosition = (submatrixColPosition+nc)%this->nCols;
  }
  
  if (communicator->mpiRoot) {
    // Receive the same data.The leading dimension of the local matrix is sendProcessRow!
    Cdgerv2d(communicator->context, nr, nc, mBlockGlobal, blockRowLength, sendProcessRow, sendProcessColumn);
  }
}

void Matrix::scatterMatrix(double *mGlobal, int sourceProcessRow, int sourceProcessCol)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Scattering a matrix is not allowed on a non-distributed matrix.", 0);
  }
  
  this->fillWithConstant(0.);
  
  for (int r = 0; r < this->nGlobRows; r += this->nBlockRows) {
    for (int c = 0; c < this->nGlobCols; c += this->nBlockCols) {
      scatterBlock(mGlobal+(this->nGlobRows*c)+r, r, c, this->nGlobRows, sourceProcessRow, sourceProcessCol);
    }
  }
}

void Matrix::scatterVector(double *vGlobal, RowColumn rowcolumn)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Scattering a vector is not allowed on a non-distributed matrix.", 0);
  }
  
  if(this->v != NULL)
  {
    delete [] this->v;
    this->v = NULL;
  }
  if(rowcolumn==row)
  {
    this->v = new double [this->nRows];
  }
  else
  {
    this->v = new double [this->nCols];
  }
  if(this->v == NULL)
  {
    misc.error("Error: There are not enough memory for performing the needed operation.", 0);
  }
  
  //sendProcessRow and sendProcessColumn defines in which process this block will be stored
  //submatrixRowPosition and submatrixColPosition defines where in the submatrix this block will be stored.
  int sendProcessRow = 0, sendProcessColumn = 0, submatrixRowPosition = 0, submatrixColPosition = 0;
  for (int r = 0; r < this->nGlobRows; r += this->nBlockRows, sendProcessRow=(sendProcessRow+1)%communicator->nProcRows) {
    sendProcessColumn = 0;
    // Number of rows to send.
    int nr = this->nBlockRows;
    if (this->nGlobRows-r < this->nBlockRows)
    {
      nr = this->nGlobRows-r;
    }
    
    for (int c = 0; c < this->nGlobCols; c += this->nBlockCols, sendProcessColumn=(sendProcessColumn+1)%communicator->nProcCols) {
      // Number of cols to send
      int nc = this->nBlockCols;
      if (this->nGlobCols-c < this->nBlockCols)
      {
	nc = this->nGlobCols-c;
      }

      if(rowcolumn==row)
      {
	if (communicator->mpiRoot) {
	  // Send the subvector to process (sendProcessRow, sendProcessColumn)
	  Cdgesd2d(communicator->context, nr, 1, vGlobal+r, this->nGlobRows, sendProcessRow, sendProcessColumn);
	}
	
	if (communicator->myRow == sendProcessRow && communicator->myCol == sendProcessColumn) {
	  // Receive the same data. The leading dimension of the local matrix is this->nRows
	  Cdgerv2d(communicator->context, nr, 1, this->v+submatrixRowPosition, this->nRows, 0, 0);
	  submatrixColPosition = (submatrixColPosition+nc)%this->nCols;
	}
      }
      else
      {
	if (communicator->mpiRoot) {
	  // Send the subvector to process (sendProcessRow, sendProcessColumn)
	  Cdgesd2d(communicator->context, nc, 1, vGlobal+c, this->nGlobCols, sendProcessRow, sendProcessColumn);
	}
	
	if (communicator->myRow == sendProcessRow && communicator->myCol == sendProcessColumn) {
	  // Receive the same data. The leading dimension of the local matrix is this->nCols
	  Cdgerv2d(communicator->context, nc, 1, this->v+submatrixColPosition, this->nCols, 0, 0);
	  submatrixColPosition = (submatrixColPosition+nc)%this->nCols;
	}
      }
      
    }
    
    if (communicator->myRow == sendProcessRow)
    {
      submatrixRowPosition = (submatrixRowPosition+nr)%this->nRows;
    }
  }

}

int * Matrix::scatterVectorRet(int *vGlobal, RowColumn rowcolumn)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Scattering a vector is not allowed on a non-distributed matrix.", 0);
  }
  
  int * vr;

  if(rowcolumn==row)
  {
    vr = new int [this->nRows];
  }
  else
  {
    vr = new int [this->nCols];
  }
  if(vr == NULL)
  {
    misc.error("Error: There are not enough memory for performing the needed operation.", 0);
  }
  
  //sendProcessRow and sendProcessColumn defines in which process this block will be stored
  //submatrixRowPosition and submatrixColPosition defines where in the submatrix this block will be stored.
  int sendProcessRow = 0, sendProcessColumn = 0, submatrixRowPosition = 0, submatrixColPosition = 0;
  for (int r = 0; r < this->nGlobRows; r += this->nBlockRows, sendProcessRow=(sendProcessRow+1)%communicator->nProcRows) {
    sendProcessColumn = 0;
    // Number of rows to send.
    int nr = this->nBlockRows;
    if (this->nGlobRows-r < this->nBlockRows)
    {
      nr = this->nGlobRows-r;
    }
    
    for (int c = 0; c < this->nGlobCols; c += this->nBlockCols, sendProcessColumn=(sendProcessColumn+1)%communicator->nProcCols) {
      // Number of cols to send
      int nc = this->nBlockCols;
      if (this->nGlobCols-c < this->nBlockCols)
      {
        nc = this->nGlobCols-c;
      }

      if(rowcolumn==row)
      {
        if (communicator->mpiRoot) {
          // Send the subvector to process (sendProcessRow, sendProcessColumn)
          Cigesd2d(communicator->context, nr, 1, vGlobal+r, this->nGlobRows, sendProcessRow, sendProcessColumn);
        }
        
        if (communicator->myRow == sendProcessRow && communicator->myCol == sendProcessColumn) {
          // Receive the same data. The leading dimension of the local matrix is this->nRows
          Cigerv2d(communicator->context, nr, 1, vr+submatrixRowPosition, this->nRows, 0, 0);
          submatrixColPosition = (submatrixColPosition+nc)%this->nCols;
        }
      }
      else
      {
        if (communicator->mpiRoot) {
          // Send the subvector to process (sendProcessRow, sendProcessColumn)
          Cigesd2d(communicator->context, nc, 1, vGlobal+c, this->nGlobCols, sendProcessRow, sendProcessColumn);
        }
        
        if (communicator->myRow == sendProcessRow && communicator->myCol == sendProcessColumn) {
          // Receive the same data. The leading dimension of the local matrix is this->nCols
          Cigerv2d(communicator->context, nc, 1, vr+submatrixColPosition, this->nCols, 0, 0);
          submatrixColPosition = (submatrixColPosition+nc)%this->nCols;
        }
      }
      
    }
    
    if (communicator->myRow == sendProcessRow)
    {
      submatrixRowPosition = (submatrixRowPosition+nr)%this->nRows;
    }
  }

  return vr;
}

int Matrix::local2global(RowColumn rc, int n)
{
  if( (n >= this->nCols && rc == column) || (n >= this->nRows && rc == row) || n < 0 )
  {
    misc.error("Error: An internal error was happened.", 0);
  }
  if(this->distribution == cyclicDistribution)  
  {
    int nProc;
    int myProc;
    int blockSize;
    if(rc == row)
    {
      nProc = communicator->nProcRows;
      myProc = communicator->myRow;
      blockSize = this->nBlockRows;
    }
    else if(rc == column)
    {
      nProc = communicator->nProcCols;
      myProc = communicator->myCol;
      blockSize = this->nBlockCols;
    }
    else
    {
      misc.error("Error: An internal error was happened.", 0);
    }
    int s = n/blockSize;
    return (nProc*s + myProc)*blockSize + (n%blockSize);
    
  }
  else if(this->distribution == global)
  {
    return n;
  }
  else
  {
    misc.error("Error: An internal error was happened. local to global transformation not allowed on this kind of distributed matrix.", 0);
    return -1;
  }
}

LocalPosition Matrix::global2local(int globalPosition, int sizeBlock, int nProcs)
{
  LocalPosition localPosition;
  
  int globalBlock = globalPosition / sizeBlock;
  localPosition.proc = globalBlock % nProcs;
  int block = globalBlock / nProcs;
  localPosition.position = ( globalPosition % sizeBlock ) + ( block * sizeBlock );
  
  return localPosition;
}

void Matrix::gatherMatrix(double *mGlobal)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Gathering a matrix is not allowed on a non-distributed matrix.", 0);
  }
  
  int sendProcessRow = 0, sendProcessColumn = 0, submatrixRowPosition = 0, submatrixColPosition = 0;
  
  sendProcessRow = 0;
  for (int r = 0; r < this->nGlobRows; r += this->nBlockRows) {
    for (int c = 0; c < this->nGlobCols; c += this->nBlockCols) {
      gatherBlock(mGlobal+this->nGlobRows*c+r, r, c, this->nGlobRows);
    }
  }
}

void Matrix::fillWithConstant(double value)
{
  #pragma omp parallel for
  for (int i = 0; i < (this->nRows*this->nCols); i++)
  {
    this->m[i] = value;
  }
}

void Matrix::fillDiagonal(double diagonal, double background)
{
  if(this->nGlobRows != this->nGlobCols)
  {
    misc.error("Error: Tried to fill the diagonal elements of a non quadratic matrix.", 0);
  }
  this->symmetric = true;
  this->uplo = 'B';
  
  if( this->distribution == cyclicDistribution )
  {
    int ia = 1;
    int ja = 1;
    pdlaset_(&this->uplo, &this->nGlobRows, &this->nGlobCols, &background, &diagonal, this->m, &ia, &ja, this->descriptor);
  }
  else if( this->distribution == diagonalDistribution )
  {
    if( background != 0. )
    {
      misc.error("Error: An internal error was happened. The background of a diagonal distributed matrix cannot be different than 0.", 0);
    }
    for(int i = 0; i < this->nRows; i++)
    {
      this->m[i] = diagonal;
    }
  }
  else
  {
    misc.error("Error: An internal error was happened. fill diagonal not allowed on this kind of distributed matrix.", 0);
  }
}

double Matrix::matrixElement(int r, int c, double defaultValue, bool *local, bool useSinglePrecision)
{
  if(this->distribution == cyclicDistribution)  
  {
    if(r >= this->nGlobRows || c >= this->nGlobCols)
    {
      misc.error("Error: An internal error was happened", 0);
    }
    
    int processRow = 0, processColumn = 0, submatrixRowPosition = 0, submatrixColPosition = 0;
    
    //sendProcessRow and sendProcessColumn defines in which process this block will be stored
    processRow = r/this->nBlockRows;
    processRow = processRow % communicator->nProcRows;
    processColumn = c/this->nBlockCols;
    processColumn = processColumn % communicator->nProcCols;
    
    if (communicator->myRow != processRow || communicator->myCol != processColumn)
    {
      *local = false;
      return defaultValue;
    }
    
    //submatrixRowPosition and submatrixColPosition defines where in the submatrix this block will be stored.
    submatrixRowPosition = r/this->nBlockRows;
    submatrixRowPosition = submatrixRowPosition/communicator->nProcRows;
    submatrixRowPosition = submatrixRowPosition*this->nBlockRows + r%this->nBlockRows;
    submatrixColPosition = c/this->nBlockCols;
    submatrixColPosition = submatrixColPosition/communicator->nProcCols;
    submatrixColPosition = submatrixColPosition*this->nBlockCols + c%this->nBlockCols;
    
    *local = true;
    if( useSinglePrecision == false )
    {
      return this->m[this->nRows*submatrixColPosition+submatrixRowPosition];
    }
    else
    {
      return double(this->mSinglePrecision[this->nRows*submatrixColPosition+submatrixRowPosition]);
    }
  }
  else if(this->distribution == global)
  {
    *local = true;
    if( useSinglePrecision == false )
    {
      return this->m[this->nRows*c+r];
    }
    else
    {
      return double(this->mSinglePrecision[this->nRows*c+r]);
    }
  }
  else
  {
    misc.error("Error: An internal error was happened. Matrix element cannot be retrieved on this kind of matrix distribution.", 0);
  }
}

void Matrix::gatherRowBlock(int r, double *blockRow, int * nr)
{
  if(this->distribution == cyclicDistribution)
  {
    for (int c = 0; c < this->nGlobCols; c += this->nBlockCols) {
      gatherBlock(blockRow+this->nGlobRows*c, r, c, this->nBlockRows);
    }
    *nr = this->nBlockRows;
    if (this->nGlobRows-r < this->nBlockRows)
    {
      *nr = this->nGlobRows-r;
    }
  }
  else
  {
    misc.error("Error: An internal error was happened. Function still not implemented.", 0);
  }
}

void Matrix::gatherColBlock(int c, double *blockCol, int * nc)
{
  if(this->distribution == cyclicDistribution)
  {
    for (int r = 0; r < this->nGlobRows; r += this->nBlockRows) {
      gatherBlock(blockCol+r, r, c, this->nGlobRows);
    }
    *nc = this->nBlockCols;
    if (this->nGlobCols-c < this->nBlockCols)
    {
      *nc = this->nGlobCols-c;
    }
  }
  else
  {
    misc.error("Error: An internal error was happened. Function still not implemented.", 0);
  }
}


void Matrix::generalResorting(Matrix *resultantMatrix, std::map<int, int> & rowsOriginDestination, std::map<int, int> & colsOriginDestination, bool unallocateThisMatrix)
{
  if( this->distribution == global || this->distribution == diagonalDistribution || this->distribution != cyclicDistribution )
  {
    misc.error("Error: An internal error was happened. Matrix resorting is not allowed on a non-distributed matrix.", 0);
  }
  
  if( rowsOriginDestination.size() == 0 || colsOriginDestination.size() == 0 )
  {
    misc.error("Error: Resorting an empty matrix or without specified indexs.", 0);
  }
  
  if( this->symmetric == true )
  {
    this->symmetrizeTriangularMatrix();
  }
  
  resultantMatrix->initParameters(rowsOriginDestination.size(), colsOriginDestination.size(), communicator->nDefaultBlockRows, communicator->nDefaultBlockCols, this->distribution, false);
  
  //Get the correspondence between the local rows/cols that must be sent and the local destination rows/cols
  
  std::vector< std::vector<int> > origProcessorRowLocalIndexs( communicator->mpiNumTasks );
  std::vector< std::vector<int> > destProcessorRowLocalIndexs( communicator->mpiNumTasks );
  int nSendingRows = 0;
  for(int lor = 0; lor < this->nRows; lor++) //lor -> local original row
  {
    int gor = local2global(row, lor); //global origin row
    if( rowsOriginDestination.count(gor) != 1 )
    {
      continue;
    }
    if( rowsOriginDestination[gor] >= resultantMatrix->nGlobRows )
    {
      misc.error("Error: An internal error was happeed when resorting a matrix. A row index is greater than resultant matrix dimensions.", 0);
    }
    LocalPosition localPosition = global2local(rowsOriginDestination[gor], resultantMatrix->nBlockRows, communicator->nProcRows);
    origProcessorRowLocalIndexs[localPosition.proc].push_back( lor );
    destProcessorRowLocalIndexs[localPosition.proc].push_back( localPosition.position );
    nSendingRows++;
  }

  std::vector< std::vector<int> > destProcessorColLocalIndexs( communicator->mpiNumTasks );
  std::vector< std::vector<int> > origProcessorColLocalIndexs( communicator->mpiNumTasks );
  int nSendingCols = 0;
  for(int loc = 0; loc < this->nCols; loc++) //loc -> local original column
  {
    int goc = local2global(column, loc); //global origin column
    if( colsOriginDestination.count(goc) != 1 )
    {
      continue;
    }
    if( colsOriginDestination[goc] >= resultantMatrix->nGlobCols )
    {
      misc.error("Error: An internal error was happeed when resorting a matrix. A column index is greater than resultant matrix dimensions.", 0);
    }
    LocalPosition localPosition = global2local(colsOriginDestination[goc], resultantMatrix->nBlockCols, communicator->nProcCols);
    origProcessorColLocalIndexs[localPosition.proc].push_back( loc );
    destProcessorColLocalIndexs[localPosition.proc].push_back( localPosition.position );
    nSendingCols++;
  }
  
  //Create vectors for sending. Vector structures are:
  
  double * sendValues = new double [ nSendingCols*nSendingRows ]; //Vector for sending array values. Values are grouped by destination processor.;
  misc.estimateMaxMemory((double(resultantMatrix->nGlobRows)*double(resultantMatrix->nGlobCols)*8.)); //Rough estimate of space used for sendValues
  std::vector< int > nSendValues( communicator->mpiNumTasks ); //Vector with the number of values for each processor.
  std::vector< int > sendIndexes; //Vector with the destination local indices. For each group of processors the structure is: number of columns, column indexes, number of rows, row indexes.
  std::vector< int > nSendIndexes( communicator->mpiNumTasks ); //Vector with the number of indexs for each processor.
  
  int sendValuesIdx = 0;
  for(int irp = 0; irp < communicator->nProcRows; irp++)
  {
    for(int icp = 0; icp < communicator->nProcCols; icp++)
    {
      int mpiRankTemp = icp + irp * communicator->nProcCols;
      nSendValues[ mpiRankTemp ] = origProcessorRowLocalIndexs[ irp ].size() * origProcessorColLocalIndexs[ icp ].size();
      if(nSendValues[ mpiRankTemp ] != 0)
      {
        nSendIndexes[ mpiRankTemp ] = origProcessorRowLocalIndexs[ irp ].size() + origProcessorColLocalIndexs[ icp ].size() + 2;
      }
      else
      {
        nSendIndexes[ mpiRankTemp ] = 2;
        sendIndexes.push_back(0);
        sendIndexes.push_back(0);
        continue;
      }
      
      sendIndexes.push_back( origProcessorColLocalIndexs[ icp ].size() );
      for(int ic = 0; ic < origProcessorColLocalIndexs[ icp ].size(); ic++)
      {
        int loc = origProcessorColLocalIndexs[ icp ][ ic ];
        sendIndexes.push_back( destProcessorColLocalIndexs[ icp ][ ic ] );
        for(int ir = 0; ir < origProcessorRowLocalIndexs[ irp ].size(); ir++)
        {
          int lor = origProcessorRowLocalIndexs[ irp ][ ir ];
          sendValues[ sendValuesIdx ] = this->m[loc*this->nRows + lor];
          sendValuesIdx++;
        }
      }
      sendIndexes.push_back( origProcessorRowLocalIndexs[ irp ].size() );
      for(int ir = 0; ir < origProcessorRowLocalIndexs[ irp ].size(); ir++)
      {
        sendIndexes.push_back( destProcessorRowLocalIndexs[ irp ][ ir ] );
      }
    }
  }
  
  //Send/receive data (indexes)
  
  std::vector< int > sendValuesDisplacements( communicator->mpiNumTasks );
  std::vector< int > sendIndexesDisplacements( communicator->mpiNumTasks );
  sendValuesDisplacements[0] = 0;
  sendIndexesDisplacements[0] = 0;
  for(int i = 1; i<communicator->mpiNumTasks; i++)
  {
    sendValuesDisplacements[ i ] = sendValuesDisplacements[ i - 1 ] + nSendValues[ i - 1 ];
    sendIndexesDisplacements[ i ] = sendIndexesDisplacements[ i - 1 ] + nSendIndexes[ i - 1 ];
  }
  
  int * nReceivedIndexes = new int [ communicator->mpiNumTasks ];
  MPI_Alltoall( &(nSendIndexes[0]), 1, MPI_INT, nReceivedIndexes, 1, MPI_INT, communicator->mpiCommunicator);
  
  int totalReceivedIndexes = 0;
  std::vector< int > receivedIndexesDisplacements( communicator->mpiNumTasks );
  totalReceivedIndexes += nReceivedIndexes[0];
  receivedIndexesDisplacements[0] = 0;
  for(int i = 1; i<communicator->mpiNumTasks; i++)
  {
    totalReceivedIndexes += nReceivedIndexes[ i ];
    receivedIndexesDisplacements[ i ] = receivedIndexesDisplacements[ i - 1 ] + nReceivedIndexes[ i - 1 ];
  }
  int * receivedIndexes = new int [ totalReceivedIndexes ];
  MPI_Alltoallv( &(sendIndexes[0]), &(nSendIndexes[0]), &(sendIndexesDisplacements[0]), MPI_INT, receivedIndexes, nReceivedIndexes, &(receivedIndexesDisplacements[0]), MPI_INT, communicator->mpiCommunicator);
  
  //Send/receive data (values)
  
  std::vector< std::vector<int> > receivedProcessorColLocalIndexs( communicator->mpiNumTasks );
  std::vector< std::vector<int> > receivedProcessorRowLocalIndexs( communicator->mpiNumTasks );
  std::vector< int > receivedValuesDisplacements( communicator->mpiNumTasks );
  std::vector< int > nReceivedValues( communicator->mpiNumTasks );
  int test = 0;
  for(int ip = 0; ip<communicator->mpiNumTasks; ip++)
  {
    int displacement = receivedIndexesDisplacements[ ip ];
    int shift = displacement;
    for(int ic = 0; ic < receivedIndexes[ shift ]; ic++)
    {
      receivedProcessorColLocalIndexs[ ip ].push_back(receivedIndexes[ shift + ic + 1 ]);
    }
    shift += receivedProcessorColLocalIndexs[ ip ].size() + 1;
    for(int ir = 0; ir < receivedIndexes[ shift ]; ir++)
    {
      receivedProcessorRowLocalIndexs[ ip ].push_back(receivedIndexes[ shift + ir + 1 ]);
    }
    
    nReceivedValues[ ip ] = receivedProcessorColLocalIndexs[ ip ].size()*receivedProcessorRowLocalIndexs[ ip ].size();
    if(ip == 0)
    {
      receivedValuesDisplacements[ ip ] = 0;
    }
    else
    {
      receivedValuesDisplacements[ ip ] = receivedValuesDisplacements[ ip - 1 ] + nReceivedValues[ ip - 1 ];
    }
    test += nReceivedValues[ ip ];
  }
  
  if( test != (resultantMatrix->nRows * resultantMatrix->nCols) )
  {
    misc.error("Error: An internal error was happened while resorting a matrix.", 0);
  }
  
  //Remove current matrix data if it is no longer needed.
  
  if( unallocateThisMatrix == true )
  {
    this->unallocateMemory();
  }
  
  //Send data
  
  double * receivedValues = new double [ resultantMatrix->nRows * resultantMatrix->nCols ];
  misc.estimateMaxMemory((double(resultantMatrix->nGlobRows)*double(resultantMatrix->nGlobCols)*8.)); //Rough estimate of space used for receivedValues
  
  MPI_Alltoallv( sendValues, &(nSendValues[0]), &(sendValuesDisplacements[0]), MPI_DOUBLE, receivedValues, &(nReceivedValues[0]), &(receivedValuesDisplacements[0]), MPI_DOUBLE, communicator->mpiCommunicator);
  
  //Allocate and free memory
  
  delete [] sendValues;
  misc.estimateMaxMemory(-(double(resultantMatrix->nGlobRows)*double(resultantMatrix->nGlobCols)*8.)); //Rough estimate of space used for sendValues
  resultantMatrix->allocateMemory();
  
  //Put received values in the local matrix
  
  for(int ip = 0; ip<communicator->mpiNumTasks; ip++)
  {
    int shift = receivedValuesDisplacements[ ip ];
    int receivedValuesIdx = 0;
    for(int ic = 0; ic < receivedProcessorColLocalIndexs[ ip ].size(); ic++)
    {
      int loc = receivedProcessorColLocalIndexs[ ip ][ ic ];
      for(int ir = 0; ir < receivedProcessorRowLocalIndexs[ ip ].size(); ir++)
      {
        int lor = receivedProcessorRowLocalIndexs[ ip ][ ir ];
        resultantMatrix->m[loc*resultantMatrix->nRows + lor] = receivedValues[ shift + receivedValuesIdx ];
        receivedValuesIdx++;
      }
    }
  }
  
  delete [] receivedIndexes;
  delete [] nReceivedIndexes;
  
  delete [] receivedValues;
  misc.estimateMaxMemory(-(double(resultantMatrix->nGlobRows)*double(resultantMatrix->nGlobCols)*8.)); //Rough estimate of space used for receivedValues
}

void Matrix::filterRowsAndColumns(Matrix *resultantMatrix, int *rows, int nElemRows, int *cols, int nElemCols, bool keep, bool unallocateThisMatrix)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Matrix filtering is not allowed on a non-distributed matrix.", 0);
  }
  
  std::map<int, int> rowsOriginDestination;
  std::map<int, int> colsOriginDestination;

  if(nElemRows > this->nGlobRows || nElemCols > this->nGlobCols)
  {
    misc.error("Error: An internal error was happened in function filterRowsAndColumns(). nElemRows > this->nGlobRows || nElemCols > this->nGlobCols.", 0);
  }

  int temp = -1;  
  for(int ir = 0; ir < nElemRows; ir++)
  {
    if(rows[ir] <= temp || rows[ir] >= this->nGlobRows || rows[ir] < 0)
    {
      misc.error("Error: An internal error was hapened when filtering a matrix. Invalid row indices.", 0);
    }
    temp = rows[ir];
  }
  temp = -1;
  for(int ic = 0; ic < nElemCols; ic++)
  {
    if(cols[ic] <= temp || cols[ic] >= this->nGlobCols || cols[ic] < 0)
    {
      misc.error("Error: An internal error was hapened when filtering a matrix. Invalid column indices.", 0);
    }
    temp = cols[ic];
  }


  if(keep)
  {
    for(int ir = 0; ir < nElemRows; ir++)
    {
      rowsOriginDestination[ rows[ir] ] = ir;
    }
    for(int ic = 0; ic < nElemCols; ic++)
    {
      colsOriginDestination[ cols[ic] ] = ic;
    }
  }
  else
  {
    int nRowsKeep = 0;
    int j = 0;
    for(int ir = 0; ir < this->nGlobRows; ir++)
    {
      if( j != nElemRows )
      {
        if( ir == rows[j] )
        {
          j++;
          continue;
        }
      }
      rowsOriginDestination[ ir ] = nRowsKeep;
      nRowsKeep++;
      if(nRowsKeep > this->nGlobRows - nElemRows)
      {
        misc.error("Error: An internal error was happened. keep = false. When filtering matrix rows/columns, there are rows in the resulting matrix not present in the original matrix.", 0);
      }
    }
    int nColsKeep = 0;
    j = 0;
    for(int ic = 0; ic < this->nGlobCols; ic++)
    {
      if( j != nElemCols )
      {
        if( ic == cols[j] )
        {
          j++;
          continue;
        }
      }
      colsOriginDestination[ ic ] = nColsKeep;
      nColsKeep++;
      if(nColsKeep > this->nGlobCols - nElemCols)
      {
        misc.error("Error: An internal error was happened. keep = false. When filtering matrix rows/columns, there are columns in the resulting matrix not present in the original matrix.", 0);
      }
    }
    
    if( (nRowsKeep != this->nGlobRows - nElemRows) || (nColsKeep != this->nGlobCols - nElemCols) )
    {
      misc.error("Error: An internal error was happened when resorting indices prior matrix filtration.", 0);
    }
  }
  
  generalResorting(resultantMatrix, rowsOriginDestination, colsOriginDestination, unallocateThisMatrix);
}

void Matrix::filterRowsAndColumns(Matrix *resultantMatrix, std::vector<int> & rows, std::vector<int> & cols, bool keep, bool unallocateThisMatrix)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Matrix filtering is not allowed on a non-distributed matrix.", 0);
  }
  
  int *rowIndxs;
  int *colIndxs;
  int nRows;
  int nCols;
  
  if(communicator->mpiRoot)
  {
    nRows = rows.size();
    nCols = cols.size();
    rowIndxs = new int [nRows];
    colIndxs = new int [nCols];
    
    for(int i = 0; i < nRows; i++)
    {
      rowIndxs[i] = rows[i];
    }
    for(int i = 0; i < nCols; i++)
    {
      colIndxs[i] = cols[i];
    }
  }
  
  communicator->broadcast(&nRows, 1);
  communicator->broadcast(&nCols, 1);
  if(!communicator->mpiRoot)
  {
    rowIndxs = new int [nRows];
    colIndxs = new int [nCols];
  }
  communicator->broadcast(rowIndxs, nRows);
  communicator->broadcast(colIndxs, nCols);
  
  this->filterRowsAndColumns(resultantMatrix, rowIndxs, nRows, colIndxs, nCols, keep, unallocateThisMatrix);
  
  delete [] rowIndxs;
  delete [] colIndxs;
}

void Matrix::filterRowsAndColumnsOld(Matrix *resultantMatrix, int *rows, int nElemRows, int *cols, int nElemCols, bool keep)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Matrix filtering is not allowed on a non-distributed matrix.", 0);
  }
  
  double * originalMatrixColsBlock;
  double * resultantMatrixColsBlock;
  int * rowsKeep;
  int * colsKeep;
  int nRowsKeep;
  int nColsKeep;
  
  if(nElemRows > this->nGlobRows || nElemCols > this->nGlobCols)
  {
    misc.error("Error: An internal error was happened in function filterRowsAndColumns(). nElemRows > this->nGlobRows || nElemCols > this->nGlobCols.", 0);
  }
  
  if(keep)
  {
    rowsKeep = rows;
    colsKeep = cols;
    nRowsKeep = nElemRows;
    nColsKeep = nElemCols;
  }
  else
  {
    nRowsKeep = 0;
    rowsKeep = new int [this->nGlobRows - nElemRows];
    int j = 0;
    for(int i = 0; i < this->nGlobRows; i++)
    {
      if( j != nElemRows )
      {
        if(i==rows[j])
        {
          j++;
          continue;
        }
      }
      rowsKeep[nRowsKeep] = i;
      nRowsKeep++;
      if(nRowsKeep > this->nGlobRows - nElemRows)
      {
	misc.error("Error: An internal error was happened. keep = false. When filtering matrix rows/columns, there are rows in the resulting matrix not present in the original matrix.", 0);
      }
    }
    nColsKeep = 0;
    colsKeep = new int [this->nGlobCols - nElemCols];
    j = 0;
    for(int i = 0; i < this->nGlobCols; i++)
    {
      if( j != nElemCols )
      {
        if(i==cols[j])
        {
          j++;
          continue;
        }
      }
      colsKeep[nColsKeep] = i;
      nColsKeep++;
      if(nColsKeep > this->nGlobCols - nElemCols)
      {
	misc.error("Error: An internal error was happened. keep = false. When filtering matrix rows/columns, there are columns in the resulting matrix not present in the original matrix.", 0);
      }
    }
    
    if( (nRowsKeep != this->nGlobRows - nElemRows) || (nColsKeep != this->nGlobCols - nElemCols) )
    {
      misc.error("Error: An internal error was happened when resorting indices prior matrix filtration.", 0);
    }
  }
  
  
  if(rowsKeep[nRowsKeep-1] >= this->nGlobRows || colsKeep[nColsKeep-1] >= this->nGlobCols || rowsKeep[0] < 0 || colsKeep[0] < 0)
  {
    misc.error("Error: An internal error was happened. When filtering matrix rows/columns, there are rows/cols in the resulting matrix not present in the original matrix.", 0);
  }
  
  resultantMatrix->initParameters(nRowsKeep, nColsKeep, this->nBlockRows, this->nBlockCols);
  
  if(communicator->mpiRoot)
  {
    originalMatrixColsBlock = new double [this->nGlobRows*this->nBlockCols];
    resultantMatrixColsBlock = new double [resultantMatrix->nGlobRows*resultantMatrix->nBlockCols];
  }
  
  int c;
  int r;
  int currentOriginalColBlock = -1;
  int resultBlockColPosition = 0;
  //This bucle collects recursively blocks of data of size this->nGlobRows*this->nBlockCols and then stores the block on the file.
  for(int ic=0; ic<nColsKeep; ic++)
  {
    c = colsKeep[ic];
    if(c/this->nBlockCols != currentOriginalColBlock) //I use this instead of c%this->nBlockCols == 0 because if an entire block is filtered, then this must not be gathered.
    {
      int originalBlockColPosition = c/this->nBlockCols;
      originalBlockColPosition *= this->nBlockCols;
      for(int r=0; r<this->nGlobRows; r+=this->nBlockRows)
      {
	this->gatherBlock((originalMatrixColsBlock+r), r, originalBlockColPosition, this->nGlobRows);
      }
      currentOriginalColBlock = c/this->nBlockCols;
    }
    
    if(communicator->mpiRoot)
    {
      for(int ir=0; ir<nRowsKeep; ir++)
      {
	r = rowsKeep[ir];
	resultantMatrixColsBlock[(ic%resultantMatrix->nBlockCols)*resultantMatrix->nGlobRows + ir] = originalMatrixColsBlock[(c%this->nBlockCols)*this->nGlobRows + r];
      }
    }
    
    if(ic%resultantMatrix->nBlockCols == resultantMatrix->nBlockCols - 1)
    {
      for(int resultBlockRowPosition=0; resultBlockRowPosition<resultantMatrix->nGlobRows; resultBlockRowPosition+=resultantMatrix->nBlockRows)
      {
	resultantMatrix->scatterBlock((resultantMatrixColsBlock+resultBlockRowPosition), resultBlockRowPosition, resultBlockColPosition, resultantMatrix->nGlobRows);
      }
      resultBlockColPosition += resultantMatrix->nBlockCols;
    }
  }
  
  if(nColsKeep%resultantMatrix->nBlockCols != 0)
  {
    for(int resultBlockRowPosition=0; resultBlockRowPosition<resultantMatrix->nGlobRows; resultBlockRowPosition+=resultantMatrix->nBlockRows)
    {
      resultantMatrix->scatterBlock((resultantMatrixColsBlock+resultBlockRowPosition), resultBlockRowPosition, resultBlockColPosition, resultantMatrix->nGlobRows);
    }
  }
  
  if(communicator->mpiRoot)
  {
    delete [] originalMatrixColsBlock;
    delete [] resultantMatrixColsBlock;
  }
  if(!keep)
  {
    delete [] rowsKeep;
    delete [] colsKeep;
  }
  
}

void Matrix::filterRowsAndColumnsOld(Matrix *resultantMatrix, std::vector<int> & rows, std::vector<int> & cols, bool keep)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Matrix filtering is not allowed on a non-distributed matrix.", 0);
  }
  
  int *rowIndxs;
  int *colIndxs;
  int nRows;
  int nCols;
  
  if(communicator->mpiRoot)
  {
    nRows = rows.size();
    nCols = cols.size();
    rowIndxs = new int [nRows];
    colIndxs = new int [nCols];
    
    for(int i = 0; i < nRows; i++)
    {
      rowIndxs[i] = rows[i];
    }
    for(int i = 0; i < nCols; i++)
    {
      colIndxs[i] = cols[i];
    }
  }
  
  communicator->broadcast(&nRows, 1);
  communicator->broadcast(&nCols, 1);
  if(!communicator->mpiRoot)
  {
    rowIndxs = new int [nRows];
    colIndxs = new int [nCols];
  }
  communicator->broadcast(rowIndxs, nRows);
  communicator->broadcast(colIndxs, nCols);
  
  this->filterRowsAndColumnsOld(resultantMatrix, rowIndxs, nRows, colIndxs, nCols, keep);
  
  delete [] rowIndxs;
  delete [] colIndxs;
}

Matrix * Matrix::redistributionToGroupedCommunicatorMatrices(Communicator * newCommunicator, std::vector<int> nGlobalRowsInGroup, std::vector<int> nGlobalColsInGroup, std::map<int, std::pair<int, int> > & rowsOriginDestination, std::map<int, std::pair<int, int> > & colsOriginDestination, bool unallocateThisMatrix)
{
  if( this->symmetric == true )
  {
    this->symmetrizeTriangularMatrix();
  }
  
  if( this->distribution == global || this->distribution == diagonalDistribution || this->distribution != cyclicDistribution )
  {
    misc.error("Error: An internal error was happened. Matrix resorting is not allowed on a non-distributed matrix.", 0);
  }
  if(newCommunicator->type != basicGroupedCommunicator)
  {
    misc.error("Error: An internal error was happened in redistributionToGroupedCommunicatorMatrices(). Wrong communicator type.", 0);
  }
  if( rowsOriginDestination.size() == 0 || colsOriginDestination.size() == 0 )
  {
    misc.error("Error: Redistributing an empty number of elements.", 0);
  }
  if( nGlobalRowsInGroup.size() != newCommunicator->nGroups || nGlobalColsInGroup.size() != newCommunicator->nGroups )
  {
    misc.error("Error: An internal error was happened in redistributionToGroupedCommunicatorMatrices(). nGlobalRowsInGroup or nGlobalColsInGroup size is different than the number of distributed communicator groups.", 0);
  }

  std::set< std::pair<int, int> > testRepetitions;
  for(std::map<int, std::pair<int, int> >::iterator it = rowsOriginDestination.begin(); it != rowsOriginDestination.end(); ++it)
  {
    if( testRepetitions.find(it->second) != testRepetitions.end() )
    {
      misc.error("Error: An internal error was happened in redistributionToGroupedCommunicatorMatrices(). Row destination coordinates repeated.", 0);
    }
    testRepetitions.insert(it->second);
    if(it->second.second != -1)
    {
      if( it->second.second >= newCommunicator->nGroups )
      {
        misc.error("Error: An internal error was happened in redistributionToGroupedCommunicatorMatrices(). Group index bigger than the number of groups.", 0);
      }
      if( it->second.first >= nGlobalRowsInGroup[it->second.second] )
      {
        misc.error("Error: An internal error was happened in redistributionToGroupedCommunicatorMatrices(). Matrix destination row bigger than the predefined dimensions.", 0);
      }
    }
    else
    {
      for(int idxGroup = 0; idxGroup < newCommunicator->nGroups; idxGroup++)
      {
        if( it->second.first >= nGlobalRowsInGroup[idxGroup] )
        {
          misc.error("Error: An internal error was happened in redistributionToGroupedCommunicatorMatrices(). Matrix destination row bigger than the predefined dimensions.", 0);
        }
      }
    }
  }
  testRepetitions.clear();
  for(std::map<int, std::pair<int, int> >::iterator it = colsOriginDestination.begin(); it != colsOriginDestination.end(); ++it)
  {
    if( testRepetitions.find(it->second) != testRepetitions.end() )
    {
      misc.error("Error: An internal error was happened in redistributionToGroupedCommunicatorMatrices(). Column destination coordinates repeated.", 0);
    }
    testRepetitions.insert(it->second);
    if(it->second.second != -1)
    {
      if( it->second.second >= newCommunicator->nGroups )
      {
        misc.error("Error: An internal error was happened in redistributionToGroupedCommunicatorMatrices(). Group index bigger than the number of groups.", 0);
      }
      if( it->second.first >= nGlobalColsInGroup[it->second.second] )
      {
        misc.error("Error: An internal error was happened in redistributionToGroupedCommunicatorMatrices(). Matrix destination column bigger than the predefined dimensions.", 0);
      }
    }
    else
    {
      for(int idxGroup = 0; idxGroup < newCommunicator->nGroups; idxGroup++)
      {
        if( it->second.first >= nGlobalColsInGroup[idxGroup] )
        {
          misc.error("Error: An internal error was happened in redistributionToGroupedCommunicatorMatrices(). Matrix destination column bigger than the predefined dimensions.", 0);
        }
      }
    }
  }
  testRepetitions.clear();
  
  Communicator * globalCommunicator = communicator;
  communicator = newCommunicator;
  Matrix *resultantMatrix = new Matrix();
  resultantMatrix->initParameters(nGlobalRowsInGroup[newCommunicator->group], nGlobalColsInGroup[newCommunicator->group], communicator->nDefaultBlockRows, communicator->nDefaultBlockCols, this->distribution, false);
  communicator = globalCommunicator;
  
  //Get the correspondence between the local rows/cols that must be sent and the local destination rows/cols
  
  std::map< std::pair<int, int>, std::vector<int> > origProcessorRowLocalIndexs; //The map key indicates the destination processor for the new communicator in terms of (dest. proc. row, group)
  std::map< std::pair<int, int>, std::vector<int> > destProcessorRowLocalIndexs; //The map key indicates the destination processor for the new communicator in terms of (dest. proc. row, group)
  
  for(int newGroup = 0; newGroup < newCommunicator->nGroups; newGroup++)
  {
    for(int proc = 0; proc < newCommunicator->groupNProcRows[newGroup]; proc++)
    {
      origProcessorRowLocalIndexs[std::pair<int, int>(proc, newGroup)].clear();
      destProcessorRowLocalIndexs[std::pair<int, int>(proc, newGroup)].clear();
    }
  }
  
  for(int lor = 0; lor < this->nRows; lor++) //lor -> local original row
  {
    int gor = local2global(row, lor); //global origin row
    if( rowsOriginDestination.count(gor) != 1 )
    {
      continue;
    }
    int gdr = rowsOriginDestination[gor].first;
    int dgroup = rowsOriginDestination[gor].second;
    if( dgroup == newCommunicator->group && gdr >= resultantMatrix->nGlobRows )
    {
      misc.error("Error: An internal error was happened when resorting a matrix. A row index is greater than resultant matrix dimensions.", 0);
    }
    
    if(dgroup != -1)
    {
      LocalPosition localPosition = global2local(gdr, resultantMatrix->nBlockRows, newCommunicator->groupNProcRows[dgroup]);
      origProcessorRowLocalIndexs[std::pair<int, int>(localPosition.proc, dgroup)].push_back( lor );
      destProcessorRowLocalIndexs[std::pair<int, int>(localPosition.proc, dgroup)].push_back( localPosition.position );
    }
    else
    {
      for(int idxGroup = 0; idxGroup < newCommunicator->nGroups; idxGroup++)
      {
        LocalPosition localPosition = global2local(gdr, resultantMatrix->nBlockRows, newCommunicator->groupNProcRows[idxGroup]);
        origProcessorRowLocalIndexs[std::pair<int, int>(localPosition.proc, idxGroup)].push_back( lor );
        destProcessorRowLocalIndexs[std::pair<int, int>(localPosition.proc, idxGroup)].push_back( localPosition.position );
      }
    }
  }

  std::map< std::pair<int, int>, std::vector<int> > destProcessorColLocalIndexs; //The map key indicates the destination processor for the new communicator in terms of (dest. proc. col., group)
  std::map< std::pair<int, int>, std::vector<int> > origProcessorColLocalIndexs; //The map key indicates the destination processor for the new communicator in terms of (dest. proc. col., group)
  
  for(int newGroup = 0; newGroup < newCommunicator->nGroups; newGroup++)
  {
    for(int proc = 0; proc < newCommunicator->groupNProcCols[newGroup]; proc++)
    {
      origProcessorColLocalIndexs[std::pair<int, int>(proc, newGroup)].clear();
      destProcessorColLocalIndexs[std::pair<int, int>(proc, newGroup)].clear();
    }
  }
  
  for(int loc = 0; loc < this->nCols; loc++) //loc -> local original column
  {
    int goc = local2global(column, loc); //global origin column
    if( colsOriginDestination.count(goc) != 1 )
    {
      continue;
    }
    int gdc = colsOriginDestination[goc].first;
    int dgroup = colsOriginDestination[goc].second;
    if( dgroup == newCommunicator->group && gdc >= resultantMatrix->nGlobCols )
    {
      misc.error("Error: An internal error was happened when resorting a matrix. A column index is greater than resultant matrix dimensions.", 0);
    }
    
    if(dgroup != -1)
    {
      LocalPosition localPosition = global2local(gdc, resultantMatrix->nBlockCols, newCommunicator->groupNProcCols[dgroup]);
      origProcessorColLocalIndexs[std::pair<int, int>(localPosition.proc, dgroup)].push_back( loc );
      destProcessorColLocalIndexs[std::pair<int, int>(localPosition.proc, dgroup)].push_back( localPosition.position );
    }
    else
    {
      for(int idxGroup = 0; idxGroup < newCommunicator->nGroups; idxGroup++)
      {
        LocalPosition localPosition = global2local(gdc, resultantMatrix->nBlockCols, newCommunicator->groupNProcCols[idxGroup]);
        origProcessorColLocalIndexs[std::pair<int, int>(localPosition.proc, idxGroup)].push_back( loc );
        destProcessorColLocalIndexs[std::pair<int, int>(localPosition.proc, idxGroup)].push_back( localPosition.position );
      }
    }
  }
  
  //Map global mpi processes to newCommunicator (proc, group) coords
  std::map<int, std::pair<int, int> > mpiToNewCommProcRowCoords;
  std::map<int, std::pair<int, int> > mpiToNewCommProcColCoords;
  int nTotalSendValues = 0;
  for(int g = 0; g<newCommunicator->nGroups; g++)
  {
    for(int pr = 0; pr<newCommunicator->groupNProcRows[g]; pr++)
    {
      for(int pc = 0; pc<newCommunicator->groupNProcCols[g]; pc++)
      {
        int globalMPIRank = (g*newCommunicator->basicGroupSize) + (pr*newCommunicator->groupNProcCols[g]) + pc;
        mpiToNewCommProcRowCoords[globalMPIRank] = std::pair<int, int>(pr, g);
        mpiToNewCommProcColCoords[globalMPIRank] = std::pair<int, int>(pc, g);
        
        nTotalSendValues += origProcessorColLocalIndexs[ mpiToNewCommProcColCoords[globalMPIRank] ].size() * origProcessorRowLocalIndexs[ mpiToNewCommProcRowCoords[globalMPIRank] ].size();
      }
    }
  }
  
  //Create vectors for sending. Vector structures are:
  
  double * sendValues = new double [ nTotalSendValues ]; //Vector for sending array values. Values are grouped by destination processor.;
  misc.estimateMaxMemory((double(resultantMatrix->nGlobRows)*double(resultantMatrix->nGlobCols)*8.)); //Rough estimate of space used for sendValues
  std::vector< int > nSendValues( communicator->mpiNumTasks ); //Vector with the number of values for each processor.
  std::vector< int > sendIndexes; //Vector with the destination local indices. For each group of processors the structure is: number of columns, column indexes, number of rows, row indexes.
  std::vector< int > nSendIndexes( communicator->mpiNumTasks ); //Vector with the number of indexs for each processor.
  
  int sendValuesIdx = 0;
  for(int irp = 0; irp < communicator->nProcRows; irp++)
  {
    for(int icp = 0; icp < communicator->nProcCols; icp++)
    {
      int mpiRankTemp = icp + irp * communicator->nProcCols;
      if( mpiToNewCommProcRowCoords.count(mpiRankTemp) == 0 || mpiToNewCommProcColCoords.count(mpiRankTemp) == 0 )
      {
        misc.error("Error: An internal error was happened when spliting a global matrix between grouped communicator matrices. Invalid MPI rank.", 0);
      }
      std::pair<int, int> irpNewCoords = mpiToNewCommProcRowCoords[mpiRankTemp];
      std::pair<int, int> icpNewCoords = mpiToNewCommProcColCoords[mpiRankTemp];
      
      nSendValues[ mpiRankTemp ] = origProcessorRowLocalIndexs[ irpNewCoords ].size() * origProcessorColLocalIndexs[ icpNewCoords ].size();
      if(nSendValues[ mpiRankTemp ] != 0)
      {
        nSendIndexes[ mpiRankTemp ] = origProcessorRowLocalIndexs[ irpNewCoords ].size() + origProcessorColLocalIndexs[ icpNewCoords ].size() + 2;
      }
      else
      {
        nSendIndexes[ mpiRankTemp ] = 2;
        sendIndexes.push_back(0);
        sendIndexes.push_back(0);
        continue;
      }
      
      sendIndexes.push_back( origProcessorColLocalIndexs[ icpNewCoords ].size() );
      for(int ic = 0; ic < origProcessorColLocalIndexs[ icpNewCoords ].size(); ic++)
      {
        int loc = origProcessorColLocalIndexs[ icpNewCoords ][ ic ];
        sendIndexes.push_back( destProcessorColLocalIndexs[ icpNewCoords ][ ic ] );
        for(int ir = 0; ir < origProcessorRowLocalIndexs[ irpNewCoords ].size(); ir++)
        {
          int lor = origProcessorRowLocalIndexs[ irpNewCoords ][ ir ];
          sendValues[ sendValuesIdx ] = this->m[loc*this->nRows + lor];
          sendValuesIdx++;
        }
      }
      sendIndexes.push_back( origProcessorRowLocalIndexs[ irpNewCoords ].size() );
      for(int ir = 0; ir < origProcessorRowLocalIndexs[ irpNewCoords ].size(); ir++)
      {
        sendIndexes.push_back( destProcessorRowLocalIndexs[ irpNewCoords ][ ir ] );
      }
    }
  }
  
  if( sendValuesIdx != nTotalSendValues )
  {
    misc.error("Error: An internal error was happened when redistributing matrices in idfferent communicators. Expected sizes of the arrays differ.", 0);
  }
  
  //Send/receive data (indexes)
  
  std::vector< int > sendValuesDisplacements( communicator->mpiNumTasks );
  std::vector< int > sendIndexesDisplacements( communicator->mpiNumTasks );
  sendValuesDisplacements[0] = 0;
  sendIndexesDisplacements[0] = 0;
  for(int i = 1; i<communicator->mpiNumTasks; i++)
  {
    sendValuesDisplacements[ i ] = sendValuesDisplacements[ i - 1 ] + nSendValues[ i - 1 ];
    sendIndexesDisplacements[ i ] = sendIndexesDisplacements[ i - 1 ] + nSendIndexes[ i - 1 ];
  }
  
  int * nReceivedIndexes = new int [ communicator->mpiNumTasks ];
  MPI_Alltoall( &(nSendIndexes[0]), 1, MPI_INT, nReceivedIndexes, 1, MPI_INT, communicator->mpiCommunicator);
  
  int totalReceivedIndexes = 0;
  std::vector< int > receivedIndexesDisplacements( communicator->mpiNumTasks );
  totalReceivedIndexes += nReceivedIndexes[0];
  receivedIndexesDisplacements[0] = 0;
  for(int i = 1; i<communicator->mpiNumTasks; i++)
  {
    totalReceivedIndexes += nReceivedIndexes[ i ];
    receivedIndexesDisplacements[ i ] = receivedIndexesDisplacements[ i - 1 ] + nReceivedIndexes[ i - 1 ];
  }
  int * receivedIndexes = new int [ totalReceivedIndexes ];
  MPI_Alltoallv( &(sendIndexes[0]), &(nSendIndexes[0]), &(sendIndexesDisplacements[0]), MPI_INT, receivedIndexes, nReceivedIndexes, &(receivedIndexesDisplacements[0]), MPI_INT, communicator->mpiCommunicator);
  
  //Send/receive data (values)
  
  std::vector< std::vector<int> > receivedProcessorColLocalIndexs( communicator->mpiNumTasks );
  std::vector< std::vector<int> > receivedProcessorRowLocalIndexs( communicator->mpiNumTasks );
  std::vector< int > receivedValuesDisplacements( communicator->mpiNumTasks );
  std::vector< int > nReceivedValues( communicator->mpiNumTasks );
  int test = 0;
  for(int ip = 0; ip<communicator->mpiNumTasks; ip++)
  {
    int displacement = receivedIndexesDisplacements[ ip ];
    int shift = displacement;
    for(int ic = 0; ic < receivedIndexes[ shift ]; ic++)
    {
      receivedProcessorColLocalIndexs[ ip ].push_back(receivedIndexes[ shift + ic + 1 ]);
    }
    shift += receivedProcessorColLocalIndexs[ ip ].size() + 1;
    for(int ir = 0; ir < receivedIndexes[ shift ]; ir++)
    {
      receivedProcessorRowLocalIndexs[ ip ].push_back(receivedIndexes[ shift + ir + 1 ]);
    }
    
    nReceivedValues[ ip ] = receivedProcessorColLocalIndexs[ ip ].size()*receivedProcessorRowLocalIndexs[ ip ].size();
    if(ip == 0)
    {
      receivedValuesDisplacements[ ip ] = 0;
    }
    else
    {
      receivedValuesDisplacements[ ip ] = receivedValuesDisplacements[ ip - 1 ] + nReceivedValues[ ip - 1 ];
    }
    test += nReceivedValues[ ip ];
  }
  
  if( test != (resultantMatrix->nRows * resultantMatrix->nCols) )
  {
    misc.error("Error: An internal error was happened while resorting a matrix.", 0);
  }
  
  //Remove current matrix data if it is no longer needed.
  
  if( unallocateThisMatrix == true )
  {
    this->unallocateMemory();
  }
  
  //Send data
  
  double * receivedValues = new double [ resultantMatrix->nRows * resultantMatrix->nCols ];
  misc.estimateMaxMemory((double(resultantMatrix->nGlobRows)*double(resultantMatrix->nGlobCols)*8.)); //Rough estimate of space used for receivedValues
  
  MPI_Alltoallv( sendValues, &(nSendValues[0]), &(sendValuesDisplacements[0]), MPI_DOUBLE, receivedValues, &(nReceivedValues[0]), &(receivedValuesDisplacements[0]), MPI_DOUBLE, communicator->mpiCommunicator);
  
  //Allocate and free memory
  
  delete [] sendValues;
  misc.estimateMaxMemory(-(double(resultantMatrix->nGlobRows)*double(resultantMatrix->nGlobCols)*8.)); //Rough estimate of space used for sendValues
  resultantMatrix->allocateMemory();
  resultantMatrix->fillWithConstant(0.);
  
  //Put received values in the local matrix
  
  for(int ip = 0; ip<communicator->mpiNumTasks; ip++)
  {
    int shift = receivedValuesDisplacements[ ip ];
    int receivedValuesIdx = 0;
    for(int ic = 0; ic < receivedProcessorColLocalIndexs[ ip ].size(); ic++)
    {
      int loc = receivedProcessorColLocalIndexs[ ip ][ ic ];
      for(int ir = 0; ir < receivedProcessorRowLocalIndexs[ ip ].size(); ir++)
      {
        int lor = receivedProcessorRowLocalIndexs[ ip ][ ir ];
        resultantMatrix->m[loc*resultantMatrix->nRows + lor] = receivedValues[ shift + receivedValuesIdx ];
        receivedValuesIdx++;
      }
    }
  }
  
  delete [] receivedIndexes;
  delete [] nReceivedIndexes;
  
  delete [] receivedValues;
  misc.estimateMaxMemory(-(double(resultantMatrix->nGlobRows)*double(resultantMatrix->nGlobCols)*8.)); //Rough estimate of space used for receivedValues
  
  return resultantMatrix;
}

Matrix* Matrix::copyToGroupedCommunicator(Communicator * groupedCommunicator)
{
  if(groupedCommunicator->type != basicGroupedCommunicator)
  {
    misc.error("Error: An internal error was happened in redistributionToGroupedCommunicatorMatrices(). Wrong communicator type.", 0);
  }
  
  Matrix * distributedCopy;
  
  if( this->distribution == cyclicDistribution )
  {
    std::vector<int> nGlobalRowsInGroup;
    std::vector<int> nGlobalColsInGroup;
    for(int ig = 0; ig < groupedCommunicator->nGroups; ig++)
    {
      nGlobalRowsInGroup.push_back(this->nGlobRows);
      nGlobalColsInGroup.push_back(this->nGlobCols);
    }
    
    std::map<int, std::pair<int, int> > rowsOriginDestination;
    for(int r = 0; r < this->nGlobRows; r++)
    {
      rowsOriginDestination[r] = std::pair<int, int>(r,-1);
    }
    std::map<int, std::pair<int, int> > colsOriginDestination;
    for(int c = 0; c < this->nGlobCols; c++)
    {
      colsOriginDestination[c] = std::pair<int, int>(c,-1);
    }
    
    distributedCopy = this->redistributionToGroupedCommunicatorMatrices(groupedCommunicator, nGlobalRowsInGroup, nGlobalColsInGroup, rowsOriginDestination, colsOriginDestination, false);
  }
  else if( this->distribution == diagonalDistribution )
  {
    int dimensions;
    if( communicator->mpiRoot == true )
    {
      dimensions = this->nGlobRows;
    }
    communicator->broadcast(&dimensions);
    
    double * values;
    if( communicator->mpiRoot != true )
    {
      values = new double [dimensions];
    }
    else
    {
      values = this->m;
    }
    communicator->broadcast(values, dimensions);
    
    Communicator * globalCommunicator = communicator;  
    communicator = groupedCommunicator;
    distributedCopy = new Matrix(diagonalDistribution, dimensions, dimensions);
    distributedCopy->setDiagonal(values, dimensions);
    communicator = globalCommunicator;
    
    if( communicator->mpiRoot != true )
    {
      delete [] values;
    }
    values = NULL;
  }
  else
  {
    misc.error("Error: An internal error was happened. Matrix copy a non-distributed or diagonal distributed matrix over grouped communicators.", 0);
  }
  
  return distributedCopy;
}

void Matrix::allocateMemory()
{
  unallocateMemory();
  
  if( this->distribution != diagonalDistribution )
  {
    this->m = new(std::nothrow) double[this->nRows*this->nCols];
    
    misc.estimateMaxMemory((double(this->nGlobRows)*double(this->nGlobCols)*8.));
  }
  else if( this->distribution == diagonalDistribution )
  {
    if(communicator->mpiRoot == true)
    {
      this->m = new(std::nothrow) double[this->nRows*this->nCols];
    }
  }
  else
  {
    misc.error("Error: An internal error was happened when allocating memory. Invalid Matrix distribution.", 0);
  }
  
  if( (this->m == NULL) && (this->nRows*this->nCols != 0) )
  {
    if( this->distribution != diagonalDistribution )
    {
      misc.estimateMaxMemory(-(double(this->nGlobRows)*double(this->nGlobCols)*8.));
    }
    
    misc.error("Error: Sorry, I am unable to allocate enough memory. I was using (rough estimate) " + getString(misc.currentMemory/(1024.*1024.*double(communicator->mpiNumTasks))) + " MB and I required: " + getString(double(this->nRows*this->nCols*sizeof(double))/(1024.*1024.)) + " additional MB for a matrix of size " + i2s(this->nGlobRows) + "x" + i2s(this->nGlobCols) + ".", 0);
  }
}

void Matrix::unallocateMemory()
{
  if(this->m != NULL)
  {
    delete [] this->m;
    this->m = NULL;
    if( this->distribution != diagonalDistribution )
    {
      misc.estimateMaxMemory(-(double(this->nGlobRows)*double(this->nGlobCols)*8.));
    }
  }
  if(this->v != NULL)
  {
    delete [] this->v;
    this->v = NULL;
  }
  if(this->mSinglePrecision != NULL)
  {
    delete [] this->mSinglePrecision;
    this->mSinglePrecision = NULL;
    if( this->distribution != diagonalDistribution )
    {
      misc.estimateMaxMemory(-(double(this->nGlobRows)*double(this->nGlobCols)*4.));
    }
  }
}

void Matrix::initParameters(int ngr, int ngc, int nbr, int nbc, DistributionType newDist, bool allocateMemoryNow)
{
  if(ngr == 0 || ngc == 0)
  {
    misc.error("Error: An internal error was happened. A matrix with zero dimensions can not be created.", 0);
  }
  if(ngr < 0 || ngc < 0)
  {
    misc.error("Error: An internal error was happened. A matrix with negative dimensions can not be created.", 0);
  }
  
  unallocateMemory();

  this->distribution = newDist;
  //reinitialize the parameters to their default values.
  defaultMatrix(this->distribution);
  
  if(this->distribution == global)
  {
    this->nGlobRows = ngr;
    this->nGlobCols = ngc;
    this->nBlockRows = ngr;
    this->nBlockCols = ngc;
    this->nRows = ngr;
    this->nCols = ngc;

    long int rTest = this->nRows;
    long int cTest = this->nCols;
    long int maxIntTest = std::numeric_limits<int>::max();
    if(rTest*cTest > maxIntTest)
    {
      misc.error("Error: Sorry, to analyse this amount of data using DISSECT, you have to increase the number of MPI processes.", 0);
    }
    
    if(allocateMemoryNow == true)
    {
      allocateMemory();
    }
  }
  else if(this->distribution == cyclicDistribution)
  {
    this->nGlobRows = ngr;
    this->nGlobCols = ngc;
    this->nBlockRows = nbr;
    this->nBlockCols = nbc;
    
    int iZero = 0;
    
    this->nRows = numroc_(&this->nGlobRows, &this->nBlockRows, &communicator->myRow, &iZero, &communicator->nProcRows);
    this->nCols = numroc_(&this->nGlobCols, &this->nBlockCols, &communicator->myCol, &iZero, &communicator->nProcCols);
    
    int firstProcRow = 0;  /// We are assuming all distributed arrays start in the process (0, 0)
    int firstProcCol = 0;  
    int errorInfo;

    
    int leadingDim = (this->nRows>1?this->nRows:1);
    descinit_(this->descriptor, &this->nGlobRows, &this->nGlobCols, &this->nBlockRows, &this->nBlockCols,
	      &firstProcRow, &firstProcCol, &communicator->context, &leadingDim, &errorInfo);
    
    if (errorInfo != 0) {
      std::cerr << communicator->myRow << " " << communicator->myCol << "******* loc:" << this->nRows << "x" << this->nCols << " g:" << this->nGlobRows << "x" << this->nGlobCols << " proc" << communicator->myRow << "x" << communicator->myCol << std::endl; std::cout.flush();
      for(int jj=0; jj< DESCRIPTOR_SIZE; jj++)
      {
	std::cerr << "- " << jj<< " " << descriptor[jj] << std::endl;
      }
      std::cerr.flush();
      std::stringstream temp;
      temp << "Process " << communicator->mpiRank << " descinit_ failed." << std::endl;
      temp << "error_info = " << errorInfo << std::endl;
      misc.error(temp.str(), 0);
    }
    
    long int rTest = this->nRows;
    long int cTest = this->nCols;
    long int maxIntTest = std::numeric_limits<int>::max();
    if(rTest*cTest > maxIntTest)
    {
      misc.error("Error: Sorry, the size of a local matrix exceeds the integer limits of this program. Maybe by increasing the number of MPI processes the problem could be solved.", 0);
    }
    
    Cblacs_barrier(communicator->context, "All");
   
    if(allocateMemoryNow == true)
    {
      allocateMemory();
    }
  }
  else if( this->distribution == diagonalDistribution )
  {
    if( ngr != ngc )
    {
      misc.error("Error: An internal error was happened. A diagonal distributed matrix must be a square matrix.", 0);
    }
    this->nGlobRows = ngr;
    this->nGlobCols = ngc;
    this->nBlockRows = 0;
    this->nBlockCols = 0;
    if( communicator->mpiRoot == true )
    {
      this->nRows = this->nGlobRows;
      this->nCols = 1;
      
      if(allocateMemoryNow == true)
      {
        allocateMemory();
      }
    }
    else
    {
      this->nRows = 0;
      this->nCols = 0;
    }
    this->symmetric = true;
  }
  else
  {
    misc.error("Error: An internal error was happened. Invalid Matrix distribution.", 0);
  }
  
  if( (allocateMemoryNow == true) && (this->m == NULL) && (this->nRows*this->nCols != 0) )
  {
    misc.error("Error: Sorry, unable to allocate enough memory.", 0);
  }
  
  if( (this->nGlobCols == 1 || this->nGlobRows == 1) && (this->distribution != diagonalDistribution) )
  {
    this->vector = true;
  }
  else
  {
    this->vector = false;
  }
}

void Matrix::writeMatrixFile(std::ofstream & file)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. This function for storing a matrix is not allowed on a non-distributed matrix.", 0);
  }
  
  double * matrixColsBlock;
  
  if(communicator->mpiRoot)
  {
    char *header = new char [4];
    header[0] = 0x1;
    header[1] = ((this->symmetric==true)?'S':'N');
    header[2] = this->uplo;
    header[3] = ((this->vector==true)?'V':'N');
    file.write((char*)header, 4);
    delete [] header;
    
    matrixColsBlock = new double [this->nGlobRows*this->nBlockCols];
  }
  
  //This bucle collects recursively blocks of data of size this->nGlobRows*this->nBlockCols and then stores the block on the file.
  for(int c=0; c<this->nGlobCols; c += this->nBlockCols)
  {
    for(int r=0; r<this->nGlobRows; r+=this->nBlockRows)
    {
      this->gatherBlock((matrixColsBlock+r), r, c, this->nGlobRows);
    }
    
    if(communicator->mpiRoot)
    {
      int nBlockColsTemp = this->nBlockCols;
      if((this->nGlobCols - c) < this->nBlockCols)
      {
	nBlockColsTemp = (this->nGlobCols - c);
      }
      file.write((char*)matrixColsBlock,sizeof(double)*nBlockColsTemp*this->nGlobRows);
    }
    
  }
  
  if(communicator->mpiRoot)
  {
    delete [] matrixColsBlock;
  }
}

void Matrix::readMatrixFile(std::ifstream & file, int ngr, int ngc, int nbr, int nbc)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. This function for reading a matrix is not allowed on a non-distributed matrix.", 0);
  }
  
  double * matrixColsBlock;
  
  this->initParameters(ngr, ngc, nbr, nbc);
  
  if(communicator->mpiRoot)
  {
    char *header = new char [4];
    file.read((char*)header, 4);
    if( header[0] != 0x1 || (header[1] != 'S' && header[1] != 'N') || (header[2] != 'L' && header[2] != 'U' && header[2] != 'B') || (header[3] != 'V' && header[3] != 'N') )
    {
      misc.error("Error: Trying to read a matrix file with a malformed header.", 0);
    }
    this->symmetric = ((header[1]=='S')?true:false);
    this->uplo = header[2];
    this->vector = ((header[3]=='V')?true:false);
    
    delete [] header;
    
    matrixColsBlock = new double [this->nGlobRows*this->nBlockCols];
  }
  communicator->broadcast(&this->symmetric);
  communicator->broadcast(&this->uplo);
  communicator->broadcast(&this->vector);
  
  //This bucle reads recursively blocks of data of size this->nGlobRows*this->nBlockCols from the file and then distribute it between nodes.
  for(int c=0; c<this->nGlobCols; c += this->nBlockCols)
  {
    if(communicator->mpiRoot)
    {
      int nBlockColsTemp = this->nBlockCols;
      if((this->nGlobCols - c) < this->nBlockCols)
      {
	nBlockColsTemp = (this->nGlobCols - c);
      }
      file.read((char*)matrixColsBlock,sizeof(double)*nBlockColsTemp*this->nGlobRows);
    }
    for(int r=0; r<this->nGlobRows; r+=this->nBlockRows)
    {
      this->scatterBlock((matrixColsBlock+r), r, c, this->nGlobRows);
    }
  }
  
  
  
  if(communicator->mpiRoot)
  {
    delete [] matrixColsBlock;
  }
}

void Matrix::writeMatrixFilev2(std::ofstream & file)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. This function for storing a matrix is not allowed on a non-distributed matrix.", 0);
  }
  
  double * matrixColsBlock;
  double * storeMatrixColsBlock;
  char storedUplo;
  
  if(communicator->mpiRoot)
  {
    if(this->symmetric == true && this->uplo == 'B')
    {
      storedUplo = 'L';
    }
    else
    {
      storedUplo = this->uplo;
    }
    
    char *header = new char [4];
    header[0] = 0x1;
    header[1] = ((this->symmetric==true)?'S':'N');
    header[2] = storedUplo;
    header[3] = ((this->vector==true)?'V':'N');
    file.write((char*)header, 4);
    delete [] header;
    
    matrixColsBlock = new double [this->nGlobRows*this->nBlockCols];
    storeMatrixColsBlock = new double [this->nGlobRows*this->nBlockCols];
  }
  
  //This bucle collects recursively blocks of data of size this->nGlobRows*this->nBlockCols and then stores the block on the file.
  for(int c=0; c<this->nGlobCols; c += this->nBlockCols)
  {
    for(int r=0; r<this->nGlobRows; r+=this->nBlockRows)
    {
      this->gatherBlock((matrixColsBlock+r), r, c, this->nGlobRows);
    }
    
    if(communicator->mpiRoot)
    {
      int nBlockColsTemp = this->nBlockCols;
      if((this->nGlobCols - c) < this->nBlockCols)
      {
        nBlockColsTemp = (this->nGlobCols - c);
      }
      
      if(this->symmetric == true) //If matrix is symmetric, store only the upper or lower part.
      {
        if(this->nGlobRows != this->nGlobCols)
        {
          misc.error("Error: An internal error was happened when writing a matrix. A symmetric nonsquare matrix is not a valid matrix.", 0);
        }
        int size;
        int idx;
        size = nBlockColsTemp*(nBlockColsTemp+1);
        size /= 2;
        if(storedUplo == 'L')
        {
          size += nBlockColsTemp*(this->nGlobCols - c - nBlockColsTemp);
          idx = 0;
          for(int lc = 0; lc < nBlockColsTemp; lc++) //Copy to storeMatrixColsBlock the lower part of matrixColsBlock.
          {
            for(int lr = c + lc; lr < this->nGlobRows; lr++)
            {
              storeMatrixColsBlock[idx] = matrixColsBlock[lc*this->nGlobRows + lr];
              idx++;
            }
          }
        }
        else if(storedUplo == 'U')
        {
          size += nBlockColsTemp*c;
          idx = 0;
          for(int lc = 0; lc < nBlockColsTemp; lc++) //Copy to storeMatrixColsBlock the upper part of matrixColsBlock.
          {
            for(int lr = 0; lr < c + lc + 1; lr++)
            {
              storeMatrixColsBlock[idx] = matrixColsBlock[lc*this->nGlobRows + lr];
              idx++;
            }
          }
        }
        else
        {
          misc.error("Error: An internal error was happened when writing a matrix. Invalid uplo for a symmetric matrix.", 0);
        }
        if(idx != size)
        {
          misc.error("Error: An internal error was happened when writing a matrix. Discordant sizes.", 0);
        }
        
        file.write((char*)storeMatrixColsBlock,sizeof(double)*size);
      }
      else
      {
        file.write((char*)matrixColsBlock,sizeof(double)*nBlockColsTemp*this->nGlobRows);
      }
    }
  }
  
  if(communicator->mpiRoot)
  {
    delete [] matrixColsBlock;
    delete [] storeMatrixColsBlock;
  }
}

void Matrix::readMatrixFilev2(std::ifstream & file, int ngr, int ngc, int nbr, int nbc)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. This function for reading a matrix is not allowed on a non-distributed matrix.", 0);
  }
  
  double * matrixColsBlock;
  double * storeMatrixColsBlock;
  
  this->initParameters(ngr, ngc, nbr, nbc);
  
  if(communicator->mpiRoot)
  {
    char *header = new char [4];
    file.read((char*)header, 4);
    if( header[0] != 0x1 || (header[1] != 'S' && header[1] != 'N') || (header[2] != 'L' && header[2] != 'U' && header[2] != 'B') || (header[3] != 'V' && header[3] != 'N') )
    {
      misc.error("Error: Trying to read a matrix file with a malformed header.", 0);
    }
    this->symmetric = ((header[1]=='S')?true:false);
    this->uplo = header[2];
    this->vector = ((header[3]=='V')?true:false);
    
    delete [] header;
    
    matrixColsBlock = new double [this->nGlobRows*this->nBlockCols];
    storeMatrixColsBlock = new double [this->nGlobRows*this->nBlockCols];
  }
  communicator->broadcast(&this->symmetric);
  communicator->broadcast(&this->uplo);
  communicator->broadcast(&this->vector);
  
  //This bucle reads recursively blocks of data of size this->nGlobRows*this->nBlockCols from the file and then distribute it between nodes.
  for(int c=0; c<this->nGlobCols; c += this->nBlockCols)
  {
    if(communicator->mpiRoot)
    {
      int nBlockColsTemp = this->nBlockCols;
      if((this->nGlobCols - c) < this->nBlockCols)
      {
        nBlockColsTemp = (this->nGlobCols - c);
      }
      
      if(this->symmetric == true) //If matrix is symmetric, store only the upper or lower part
      {
        if(this->nGlobRows != this->nGlobCols)
        {
          misc.error("Error: An internal error was happened when reading a matrix. A symmetric nonsquare matrix is not a valid matrix.", 0);
        }
        int size;
        size = nBlockColsTemp*(nBlockColsTemp+1);
        size /= 2;
        if(this->uplo == 'L')
        {
          size += nBlockColsTemp*(this->nGlobCols - c - nBlockColsTemp);
        }
        else if(this->uplo == 'U')
        {
          size += nBlockColsTemp*c;
        }
        else
        {
          misc.error("Error: An internal error was happened when reading a matrix. Invalid uplo for a symmetric matrix.", 0);
        }
        
        file.read((char*)storeMatrixColsBlock,sizeof(double)*size);
        
        int idx;
        if(this->uplo == 'L')
        {
          idx = 0;
          for(int lc = 0; lc < nBlockColsTemp; lc++) //Copy to the the lower part of matrixColsBlock the data in storeMatrixColsBlock.
          {
            for(int lr = c + lc; lr < this->nGlobRows; lr++)
            {
              matrixColsBlock[lc*this->nGlobRows + lr] = storeMatrixColsBlock[idx];
              idx++;
            }
          }
        }
        else if(this->uplo == 'U')
        {
          idx = 0;
          for(int lc = 0; lc < nBlockColsTemp; lc++) //Copy to the the upper part of matrixColsBlock the data in storeMatrixColsBlock.
          {
            for(int lr = 0; lr < c + lc + 1; lr++)
            {
              matrixColsBlock[lc*this->nGlobRows + lr] = storeMatrixColsBlock[idx];
              idx++;
            }
          }
        }
        if(idx != size)
        {
          misc.error("Error: An internal error was happened when reading a matrix. Discordant sizes.", 0);
        }
      }
      else
      {
        file.read((char*)matrixColsBlock,sizeof(double)*nBlockColsTemp*this->nGlobRows);
      }
    }
    
    for(int r=0; r<this->nGlobRows; r+=this->nBlockRows)
    {
      this->scatterBlock((matrixColsBlock+r), r, c, this->nGlobRows);
    }
  }
  
  if(communicator->mpiRoot)
  {
    delete [] matrixColsBlock;
    delete [] storeMatrixColsBlock;
  }
}

void Matrix::writeMatrixMPI(std::string fname, char * header, int offset)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. This function for storing a matrix is not allowed on a non-distributed matrix.", 0);
  }
  
  int dims[] = {this->nGlobRows, this->nGlobCols};
  int dargs[] = {this->nBlockRows, this->nBlockCols};
  int distribs[] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
  int dim[] = {communicator->nProcRows, communicator->nProcCols};
  char nat[] = "native";
  int rc;
  MPI_Datatype dcarray; 
  MPI_File cFile;
  MPI_Status status;

  MPI_Type_create_darray(communicator->mpiNumTasks, communicator->mpiRank, 2, dims, distribs, dargs, dim, MPI_ORDER_FORTRAN, MPI_DOUBLE, &dcarray); 
  MPI_Type_commit(&dcarray);
  
  std::vector<char> fn(fname.begin(), fname.end());
  fn.push_back('\0');
  rc = MPI_File_open(communicator->mpiCommunicator, &fn[0], MPI_MODE_EXCL | MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &cFile);
  if(rc)
  {
    MPI_File_delete (&fn[0], MPI_INFO_NULL);
    communicator->barrier();
    rc = MPI_File_open(communicator->mpiCommunicator, &fn[0], MPI_MODE_EXCL | MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &cFile);
  }
  if(rc){
    std::stringstream ss;
    ss << "Error: Failed to open file: " << rc;
    misc.error(ss.str(), 0);
  }
  else
  {
    int err0 = MPI_SUCCESS;
    if(communicator->mpiRoot && offset != 0)
    {
      err0 = MPI_File_write(cFile, header, offset, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    communicator->barrier();
    int err1 = MPI_File_set_view(cFile, MPI_Offset(offset), MPI_DOUBLE, dcarray, nat, MPI_INFO_NULL);
    int err2 = MPI_File_write_all(cFile, this->m, this->nRows*this->nCols, MPI_DOUBLE, &status);    
    if(err0 != MPI_SUCCESS || err1 != MPI_SUCCESS || err2 != MPI_SUCCESS)
    {
      misc.error("Error: An error was happened when writing a matrix. Is there enough space in the hard disk?", 0);
    }
  }
  MPI_File_close(&cFile);
  MPI_Type_free(&dcarray);
}

void Matrix::readMatrixMPI(std::string fname, int offset)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. This function for reading a matrix is not allowed on a non-distributed matrix.", 0);
  }
  
  int dims[] = {this->nGlobRows, this->nGlobCols};
  int dargs[] = {this->nBlockRows, this->nBlockCols};
  int distribs[] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
  int dim[] = {communicator->nProcRows, communicator->nProcCols};
  char nat[] = "native";
  int rc;
  MPI_Datatype dcarray; 
  MPI_File cFile;
  MPI_Status status;
  
  MPI_Type_create_darray(communicator->mpiNumTasks, communicator->mpiRank, 2, dims, distribs, dargs, dim, MPI_ORDER_FORTRAN, MPI_DOUBLE, &dcarray); 
  MPI_Type_commit(&dcarray);
  
  std::vector<char> fn(fname.begin(), fname.end());
  fn.push_back('\0');
  rc = MPI_File_open(communicator->mpiCommunicator, &fn[0], MPI_MODE_RDONLY, MPI_INFO_NULL, &cFile);
  if(rc){
    std::stringstream ss;
    ss << "Error: Failed to open file: " << rc;
    misc.error(ss.str(), 0);
  }
  else
  {
    int err1 = MPI_File_set_view(cFile, MPI_Offset(offset), MPI_DOUBLE, dcarray, nat, MPI_INFO_NULL);
    int err2 = MPI_File_read_all(cFile, this->m, this->nRows*this->nCols, MPI_DOUBLE, &status);    
    if(err1 != MPI_SUCCESS || err2 != MPI_SUCCESS)
    {
      misc.error("Error: An error was happened when reading a matrix. Could the file be corrupt?", 0);
    }
  }
  MPI_File_close(&cFile);
  MPI_Type_free(&dcarray);
}

void Matrix::packMatrices(Matrix * m1, Matrix * m2)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Matrix packing is not allowed on a non-distributed matrix.", 0);
  }
  
  double * m1ColsBlock;
  double * m2ColsBlock;
  double * packedMatrixColsBlock;
  
  if(m1->symmetric == false || m2->symmetric == false)
  {
    misc.error("Error: An internal error was happened. Only symmetric matrices can be packed.", 0);
  }
  
  if(m1->nGlobRows != m2->nGlobRows || m1->nGlobCols != m2->nGlobCols || m1->nGlobRows != m1->nGlobCols)
  {
    misc.error("Error: An internal error was happened. Matrices must be symmetric and have same dimensions.", 0);
  }
  
  if(m1->nBlockRows != m2->nBlockRows || m1->nBlockCols != m2->nBlockCols)
  {
    misc.error("Error: An internal error was happened. Matrices must be distributed using same block structure.", 0);
  }
  
  if(m1->uplo != 'L' && m1->uplo != 'B')
  {
    m1->symmetrizeTriangularMatrix();
  }
  if(m2->uplo != 'U' && m2->uplo != 'B')
  {
    m2->symmetrizeTriangularMatrix();
  }
  
  this->initParameters(m1->nGlobRows + 1, m1->nGlobCols, m1->nBlockRows, m1->nBlockCols);
  
  if(communicator->mpiRoot)
  {
    m1ColsBlock = new double [m1->nGlobRows*m1->nBlockCols];
    m2ColsBlock = new double [m2->nGlobRows*m2->nBlockCols];
    packedMatrixColsBlock = new double [this->nGlobRows*this->nBlockCols];
  }
  
  
  //This bucle merges m1 and m2 o this.
  for(int c=0; c<this->nGlobCols; c += this->nBlockCols)
  {
    for(int r=0; r<m1->nGlobRows; r+=m1->nBlockRows)
    {
      m1->gatherBlock((m1ColsBlock+r), r, c, m1->nGlobRows);
      m2->gatherBlock((m2ColsBlock+r), r, c, m2->nGlobRows);
    }
    
    if(communicator->mpiRoot)
    {
      int nBlockColsTemp = this->nBlockCols;
      if((this->nGlobCols - c) < this->nBlockCols)
      {
        nBlockColsTemp = (this->nGlobCols - c);
      }
      
      for(int lc = 0; lc < nBlockColsTemp; lc++) 
      {
        for(int lr = 0; lr <= lc+c; lr++)
        {
          packedMatrixColsBlock[lc*this->nGlobRows + lr] = m2ColsBlock[lc*m2->nGlobRows + lr];
        }
        for(int lr = lc+c; lr < m1->nGlobRows; lr++)
        {
          packedMatrixColsBlock[lc*this->nGlobRows + (lr + 1)] = m1ColsBlock[lc*m1->nGlobRows + lr];
        }
      }
    }
    
    for(int r=0; r<this->nGlobRows; r+=this->nBlockRows)
    {
      this->scatterBlock((packedMatrixColsBlock+r), r, c, this->nGlobRows);
    }
  }
  
  if(communicator->mpiRoot)
  {
    delete [] m1ColsBlock;
    delete [] m2ColsBlock;
    delete [] packedMatrixColsBlock;
  }
}

void Matrix::unpackMatrices(Matrix * m1, Matrix * m2, bool unallocateThisMatrix)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Matrix unpacking is not allowed on a non-distributed matrix.", 0);
  }
  
  if(this->nGlobRows - 1 != this->nGlobCols)
  {
    misc.error("Error: An internal error was happened. When unpacking two symmetric matrices: the relation of rows/columns must be: nGlobCols == nGlobRows - 1.", 0);
  }
  
  int filterIndex = 0;
  this->filterRowsAndColumns(m1, &filterIndex, 1, &filterIndex, 0, false);
  m1->symmetric = true;
  m1->uplo = 'L';
  if(m1->nGlobRows != this->nGlobCols || m1->nGlobCols != this->nGlobCols)
  {
    misc.error("Error: An internal error was happened when unpacking two matrices.", 0);
  }
  
  filterIndex = this->nGlobRows - 1;
  Matrix * m2U = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->filterRowsAndColumns(m2U, &filterIndex, 1, &filterIndex, 0, false, unallocateThisMatrix);
  m2U->symmetric = true;
  m2U->uplo = 'U';
  m2->transpose(m2U);
  delete m2U;
  
  if(m2->nGlobRows != this->nGlobCols || m2->nGlobCols != this->nGlobCols)
  {
    misc.error("Error: An internal error was happened when unpacking two matrices.", 0);
  }
}

void Matrix::unpackMatricesOld(Matrix * m1, Matrix * m2)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Matrix unpacking is not allowed on a non-distributed matrix.", 0);
  }
  
  double * m1ColsBlock;
  double * m2ColsBlock;
  double * packedMatrixColsBlock;
  
  if(this->nGlobRows - 1 != this->nGlobCols)
  {
    misc.error("Error: An internal error was happened. When unpacking two symmetric matrices: the relation of rows/columns must be: nGlobCols == nGlobRows - 1.", 0);
  }
  
  m1->initParameters(this->nGlobRows - 1, this->nGlobCols, this->nBlockRows, this->nBlockCols);
  m1->symmetric = true;
  m1->uplo = 'L';
  
  Matrix * m2U = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  m2U->initParameters(this->nGlobRows - 1, this->nGlobCols, this->nBlockRows, this->nBlockCols);
  m2U->symmetric = true;
  m2U->uplo = 'U';
  
  if(communicator->mpiRoot)
  {
    m1ColsBlock = new double [m1->nGlobRows*m1->nBlockCols];
    m2ColsBlock = new double [m2U->nGlobRows*m2U->nBlockCols];
    packedMatrixColsBlock = new double [this->nGlobRows*this->nBlockCols];
  }
  
  
  //This bucle merges splits this into m1 and m2.
  for(int c=0; c<this->nGlobCols; c += this->nBlockCols)
  {
    for(int r=0; r<this->nGlobRows; r+=this->nBlockRows)
    {
      this->gatherBlock((packedMatrixColsBlock+r), r, c, this->nGlobRows);
    }
    
    if(communicator->mpiRoot)
    {
      int nBlockColsTemp = this->nBlockCols;
      if((this->nGlobCols - c) < this->nBlockCols)
      {
        nBlockColsTemp = (this->nGlobCols - c);
      }
      
      for(int lc = 0; lc < nBlockColsTemp; lc++) 
      {
        for(int lr = 0; lr <= lc+c; lr++)
        {
          m2ColsBlock[lc*m2U->nGlobRows + lr] = packedMatrixColsBlock[lc*this->nGlobRows + lr];
        }
        for(int lr = lc+c; lr < m1->nGlobRows; lr++)
        {
          m1ColsBlock[lc*m1->nGlobRows + lr] = packedMatrixColsBlock[lc*this->nGlobRows + (lr + 1)];
        }
      }
    }
    
    for(int r=0; r<m1->nGlobRows; r+=m1->nBlockRows)
    {
      m1->scatterBlock((m1ColsBlock+r), r, c, m1->nGlobRows);
      m2U->scatterBlock((m2ColsBlock+r), r, c, m2U->nGlobRows);
    }
  }
  
  if(communicator->mpiRoot)
  {
    delete [] m1ColsBlock;
    delete [] m2ColsBlock;
    delete [] packedMatrixColsBlock;
  }
  
  m2->transpose(m2U);
  delete m2U;
}

void Matrix::multiply(Matrix * m1, char t1, Matrix * m2, char t2, double scale, subMatrix smr)
{
  if( m1->distribution == diagonalDistribution || m2->distribution == diagonalDistribution)
  {
    multiplyDiagonalMatrixMatrix(m1, t1, m2, t2, scale, smr);
  }
  else
  {
    if(m1->vector == true && m2->vector == true)
    {
      multiplyMatrixMatrix(m1, t1, m2, t2, scale, smr); //Needs to be implemented
    }
    else if(m1->vector == false && m2->vector == true)
    {
      multiplyMatrixVector(m1, t1, m2, t2, scale, smr); //Needs to be implemented
    }
    else if(m1->vector == true && m2->vector == false)
    {
      multiplyMatrixVector(m1, t1, m2, t2, scale, smr); //Needs to be implemented
    }
    else if(m1->vector == false && m2->vector == false)
    {
      multiplyMatrixMatrix(m1, t1, m2, t2, scale, smr);
    }
  }
}

void Matrix::multiplyMatrixVector(Matrix * m1, char t1, Matrix * m2, char t2, double scale, subMatrix smr)
{
  int m = m1->nGlobRows;
  int n = m2->nGlobCols;
  double alpha = scale;
  int ia = 1;
  int ja = 1;
  int ix = 1;
  int jx = 1;
  double beta = 0.;
  int iy = 1;
  int jy = 1;
  
  //To be implemented
  //pdgemv_(trans, &m, &n, &alpha, m1->m, &ia, &ja, m1->description, v1->m, &ix, &jx, x->description, incx, &beta, this->m, &iy, &jy, this->description, incy)
  multiplyMatrixMatrix(m1, t1, m2, t2, scale, smr);
}

void Matrix::multiplyMatrixMatrix(Matrix * m1, char t1, Matrix * m2, char t2, double scale, subMatrix smr)
{
  if( m1->distribution == global || m1->distribution == diagonalDistribution || m2->distribution == global || m2->distribution == diagonalDistribution || this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Wrong multiplication function for a on a non-distributed matrix.", 0);
  }
  
  int m1NGlobRows, m1NGlobCols, m1NBlockRows, m1NBlockCols;
  int m2NGlobRows, m2NGlobCols, m2NBlockRows, m2NBlockCols;
  
  //Switch the number of columns/rows of the matrices depending whether the transpose will be used.
  if( t1=='T' || t1=='t' )
  {
    m1NGlobRows = m1->nGlobCols, 
    m1NGlobCols = m1->nGlobRows;
    m1NBlockRows = m1->nBlockCols;
    m1NBlockCols = m1->nBlockRows;
    
    t1 = 'T';
  }
  else
  {
    if ( t1 != 'N' && t1 != 'n')
    {
      misc.error("Error: An internal error was happened", 0);
    }
    m1NGlobRows = m1->nGlobRows, 
    m1NGlobCols = m1->nGlobCols;
    m1NBlockRows = m1->nBlockRows;
    m1NBlockCols = m1->nBlockCols;
    
    t1 = 'N';
  }
  
  if( t2=='T' || t2=='t' )
  {
    m2NGlobRows = m2->nGlobCols, 
    m2NGlobCols = m2->nGlobRows;
    m2NBlockRows = m2->nBlockCols;
    m2NBlockCols = m2->nBlockRows;
    
    t2 = 'T';
  }
  else
  {
    if ( t2 != 'N' && t2 != 'n')
    {
      misc.error("Error: An internal error was happened", 0);
    }
    m2NGlobRows = m2->nGlobRows, 
    m2NGlobCols = m2->nGlobCols;
    m2NBlockRows = m2->nBlockRows;
    m2NBlockCols = m2->nBlockCols;
    
    t2 = 'N';
  }
  
  //The result will be stored in a submatrix? Then, it is supposed that the result matrix is already allocated.
  if(!smr.active)
  {
    //Initiate parameters for the result matrix
    this->initParameters(m1NGlobRows, m2NGlobCols, m1NBlockRows, m2NBlockCols);
  }
  else
  {
    if(this->m == NULL)
    {
      misc.error("Error: An internal error was happened", 0);
    }
    if(m1NGlobRows != smr.nr || m2NGlobCols != smr.nc)
    {
      misc.error("Error: An internal error was happened when multiplying two matrices. There is a discrepancy between the subMatrix dimensions and the resultant matrix.", 0);
    }
  }
  
  if(m1NGlobCols != m2NGlobRows)
  {
    misc.error("Error: An internal error was happened. The dimensions of the matrices do not match.", 0);
  }
  
  //char transa = 'N';
  //char transb = 'T';
  int m = m1NGlobRows;
  int n = m2NGlobCols;
  int k = m1NGlobCols;
  double alpha = scale;
  //a, 
  int ia = 1;
  int ja = 1;
  //desca
  //b
  int ib = 1;
  int jb = 1;
  //descb
  double beta = 0.;
  //c
  int ic = 1;
  int jc = 1;
  //descc
  
  //If submatrix is specified, change where the resulting matrix submatrix will be copied.
  if(smr.active)
  {
    this->symmetric = false;
    if(smr.ir+smr.nr>this->nGlobRows || smr.ic+smr.nc>this->nGlobCols)
    {
      misc.error("Error: An internal error was happened. Specified submatrix does not fit in the resultant matrix.", 0);
    }
    //Change ic and jc
    ic = smr.ir + 1;
    jc = smr.ic + 1;
  }
  
  if(m1->symmetric == true || m2->symmetric == true)
  {
    if(m1->uplo != 'B' && m2->uplo != 'B')
    {
      m2->symmetrizeTriangularMatrix();
    }
    if(t1 == 'T' || t2 == 'T')
    {
      misc.error("Error: Multiplication between a symmetric matrix with transposition. Still not implemented.", 0);
    }

    char side;
    char tempuplo;
    if(m1->uplo != 'B')
    {
      if(m1->symmetric == false) { misc.error("Error: An internal error was happened", 0); }
      side = 'L';
      tempuplo = m1->uplo;
    }
    else if(m2->uplo != 'B')
    {
      if(m2->symmetric == false) { misc.error("Error: An internal error was happened", 0); }
      side = 'R';
      tempuplo = m2->uplo;
    }
    else
    {
      if(m1->symmetric == true)
      {
	side = 'L';
	tempuplo = 'L';
      }
      else
      {
	side = 'R';
	tempuplo = 'L';
      }
    }
    //m1->showGlobal("m1");
    //std::cout << side << " " << m1->uplo << " " <<m<< " " <<n<< " " <<alpha<< " " <<m1->m<< " " <<ia<< " " <<ja<< " " <<m1->descriptor<< " " <<m2->m<< " " <<ib<< " " <<jb<< " " <<m2->descriptor<< " " <<beta<< " " <<this->m<< " " <<ic<< " " <<jc<< " " << this->descriptor << std::endl; std::cout.flush();
    //std::cout << this->nGlobRows << " " << this->nGlobCols << " " <<m<< " " <<n<< " " <<alpha<< " " <<m1->m<< " " <<ia<< " " <<ja<< " " <<m1->descriptor<< " " <<m2->m<< " " <<ib<< " " <<jb<< " " <<m2->descriptor<< " " <<beta<< " " <<this->m<< " " <<ic<< " " <<jc<< " " << this->descriptor << std::endl; std::cout.flush();
    //m2->showGlobal("m2");
    this->symmetric = false;
    this->uplo = 'B';
    if(side == 'L')
    {
      pdsymm_(&side, &tempuplo, &m, &n, &alpha, m1->m, &ia, &ja, m1->descriptor, m2->m, &ib, &jb, m2->descriptor, &beta, this->m, &ic, &jc, this->descriptor);
    }
    else
    {
      pdsymm_(&side, &tempuplo, &m, &n, &alpha, m2->m, &ia, &ja, m2->descriptor, m1->m, &ib, &jb, m1->descriptor, &beta, this->m, &ic, &jc, this->descriptor);
    }
  }
  else if(smr.active == false && m1==m2 && ((t1 =='T' && t2 == 'N') || (t1 =='N' && t2 == 'T')) )
  {
    this->symmetric = true;
    this->uplo = 'L';
    
    pdsyrk_(&this->uplo, &t1, &m, &k, &alpha, m1->m, &ia, &ja, m1->descriptor, &beta, this->m, &ic, &jc, this->descriptor);
  }
  else
  {
    this->symmetric = false;
    this->uplo = 'B';
    pdgemm_(&t1, &t2, &m, &n, &k, &alpha, m1->m, &ia, &ja, m1->descriptor, m2->m, &ib, &jb, m2->descriptor, &beta, this->m, &ic, &jc, this->descriptor);
  }
  
}

void Matrix::multiplyDiagonalMatrixMatrix(Matrix * m1, char t1, Matrix * m2, char t2, double scale, subMatrix smr)
{
  if( m1->distribution != diagonalDistribution && m2->distribution != diagonalDistribution )
  {
    misc.error("Error: An internal error was happened. Wrong multiplication function for a on a non-distributed matrix.", 0);
  }
  
  if( smr.active == true )
  {
    Matrix * m1Temp;
    Matrix * m2Temp;
    if(m1->distribution == diagonalDistribution)
    {
      m1Temp = new Matrix(cyclicDistribution, m1->nGlobRows, m1->nGlobCols);
      m1Temp->fillWithConstant(0.);
      m1Temp->setDiagonal(m1->m, m1->nGlobRows);
      m1Temp->symmetric = true;
    }
    else
    {
      m1Temp = m1;
    }
    if(m2->distribution == diagonalDistribution)
    {
      m2Temp = new Matrix(cyclicDistribution, m2->nGlobRows, m2->nGlobCols);
      m2Temp->fillWithConstant(0.);
      m2Temp->setDiagonal(m2->m, m2->nGlobRows);
      m2Temp->symmetric = true;
    }
    else
    {
      m2Temp = m2;
    }
    multiplyMatrixMatrix(m1Temp, t1, m2Temp, t2, scale, smr);
    if(m1Temp != m1)
    {
      delete m1Temp;
    }
    if(m2Temp != m2)
    {
      delete m2Temp;
    }
  }
  else
  {
    if( m1->distribution == diagonalDistribution && m2->distribution == diagonalDistribution ) //diagonal by diagonal matrix product
    {
      this->duplicateMatrix(m1);
      this->elementWiseMultiplication(m2, scale);
    }
    else //diagonal by non-diagonal matrix product
    {
      //Get diagonal and non-diagonal matrices
      Matrix *nonDiagonalMatrix;
      Matrix *diagonalMatrix;
      char transposeNonDiagonal;
      if( m1->distribution == diagonalDistribution )
      {
        diagonalMatrix = m1;
        nonDiagonalMatrix = m2;
        transposeNonDiagonal = t2;
      }
      else
      {
        diagonalMatrix = m2;
        nonDiagonalMatrix = m1;
        transposeNonDiagonal = t1;
      }
      
      //Copy non-diagonal matrix to this and transpose (whether necessary)
      if( transposeNonDiagonal=='T' || transposeNonDiagonal=='t' )
      {
        this->transpose(nonDiagonalMatrix);
      }
      else
      {
        if ( transposeNonDiagonal != 'N' && transposeNonDiagonal != 'n')
        {
          misc.error("Error: An internal error was happened.", 0);
        }
        this->duplicateMatrix(nonDiagonalMatrix);
      }
      
      if(this->symmetric == true)
      {
        this->symmetrizeTriangularMatrix();
      }
      this->symmetric = false;
      
      //Perform the product
      if( m1->distribution == diagonalDistribution )
      {
        if( this->nGlobRows != diagonalMatrix->nGlobCols )
        {
          misc.error("Error: An internal error was happened. The dimensions of the matrices do not match.", 0);
        }
        
        this->scatterVector(diagonalMatrix->m, row);
        
        #pragma omp parallel for
        for(int c = 0; c<this->nCols; c++)
        {
          for(int r = 0; r<this->nRows; r++)
          {
            this->m[c*this->nRows + r] *= scale*this->v[r];
          }
        }
      }
      else
      {
        if( this->nGlobCols != diagonalMatrix->nGlobRows )
        {
          misc.error("Error: An internal error was happened. The dimensions of the matrices do not match.", 0);
        }
        
        this->scatterVector(diagonalMatrix->m, column);
        
        #pragma omp parallel for
        for(int c = 0; c<this->nCols; c++)
        {
          for(int r = 0; r<this->nRows; r++)
          {
            this->m[c*this->nRows + r] *= scale*this->v[c];
          }
        }
      }
      
    }
  }
}

void Matrix::add(Matrix * m1, double beta, double alpha, subMatrix smt, subMatrix sm1)
{
  if( this->distribution == diagonalDistribution || m1->distribution == diagonalDistribution )
  {
    addDiagonalMatrix(m1, beta, alpha, smt, sm1);
  }
  else
  {
    addMatrix(m1, beta, alpha, smt, sm1);
  }
}

void Matrix::addDiagonalMatrix(Matrix * m1, double beta, double alpha, subMatrix smt, subMatrix sm1)
{
  if( this->distribution == diagonalDistribution && m1->distribution == diagonalDistribution )
  {
    if(smt.active == false || sm1.active == false)
    {
      if(this->nGlobRows != m1->nGlobRows || this->nGlobCols != m1->nGlobCols )
      {
        misc.error("Error: An internal error was happened", 0);
      }
      if(communicator->mpiRoot == true)
      {
        for(int i = 0; i < this->nGlobRows; i++)
        {
          this->m[ i ] = beta*this->m[ i ] + alpha*m1->m[ i ];
        }
      }
    }
    else
    {
      if(sm1.ir != sm1.ic || sm1.nr != sm1.nc || smt.ir != smt.ic || smt.nr != smt.nc)
      {
        misc.error("Error: An internal error was happened when adding two diagonal distributed matrices. Non-symmetric submatrices cannot be defined.", 0);
      }
      if(smt.ir+smt.nr>this->nGlobRows || smt.ic+smt.nc>this->nGlobCols)
      {
        misc.error("Error: An internal error was happened", 0);
      }
      if(sm1.ir+sm1.nr>m1->nGlobRows || sm1.ic+sm1.nc>m1->nGlobCols)
      {
        misc.error("Error: An internal error was happened", 0);
      }
      if(smt.nr != sm1.nr || smt.nc != sm1.nc)
      {
        misc.error("Error: An internal error was happened", 0);
      }
      if(communicator->mpiRoot == true)
      {
        for(int i = 0; i < sm1.nr; i++)
        {
          this->m[ smt.ir + i ] = beta*this->m[ smt.ir + i ] + alpha*m1->m[ sm1.ir + i ];
        }
      }
    }
  }
  else if( this->distribution == diagonalDistribution && m1->distribution != diagonalDistribution )
  {
    double * diagonalTemp;
    int dimension = this->nGlobRows;
    if(communicator->mpiRoot)
    {
      diagonalTemp = new double[dimension];
      for(int i = 0; i<dimension; i++)
      {
        diagonalTemp[i] = this->m[i];
      }
    }
    initParameters(dimension, dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    fillWithConstant(0.);
    setDiagonal(diagonalTemp, dimension);
    this->symmetric = true;
    
    addMatrix(m1, beta, alpha, smt, sm1);
    
    if(communicator->mpiRoot)
    {
      delete [] diagonalTemp;
    }
  }
  else if( this->distribution != diagonalDistribution && m1->distribution == diagonalDistribution )
  {
    Matrix * m1Temp = new Matrix(cyclicDistribution, m1->nGlobRows, m1->nGlobCols);
    m1Temp->fillWithConstant(0.);
    m1Temp->setDiagonal(m1->m, m1->nGlobRows);
    m1Temp->symmetric = true;
    
    addMatrix(m1Temp, beta, alpha, smt, sm1);
    
    delete m1Temp;
  }
  else
  {
    misc.error("Error: An internal error was happened.", 0);
  }
}

void Matrix::addMatrix(Matrix * m1, double beta, double alpha, subMatrix smt, subMatrix sm1)
{
  int ia, ja, ic, jc, m, n;
  char t1 = 'N';
  
  if(smt.active == false || sm1.active == false)
  {
    if(this->nGlobRows != m1->nGlobRows || this->nGlobCols != m1->nGlobCols )
    {
      misc.error("Error: An internal error was happened", 0);
    }
    ia = 1;
    ja = 1;
    ic = 1;
    jc = 1;
    m = this->nGlobRows;
    n = this->nGlobCols;
    
    if(this->symmetric == true && this->symmetric == true)
    {
      this->symmetric = true;
    }
    else
    {
      this->symmetric = false;
    }
    if(this->uplo != 'B' || m1->uplo != 'B')
    {
      if(this->uplo != 'B' && m1->uplo != 'B' && this->uplo != m1->uplo)
      {
	misc.error("Error: An internal error was happened", 0);
      }
      if(this->uplo == 'B' && m1->uplo != 'B')
      {
	if(this->symmetric)
	{
	  this->uplo = m1->uplo; //The resultant matrix has the uplo property according to the two original matrices.
	}
	else
	{
	  m1->symmetrizeTriangularMatrix();
	}
      }
    }
  }
  else
  {
    if(smt.ir+smt.nr>this->nGlobRows || smt.ic+smt.nc>this->nGlobCols)
    {
      misc.error("Error: An internal error was happened", 0);
    }
    if(sm1.ir+sm1.nr>m1->nGlobRows || sm1.ic+sm1.nc>m1->nGlobCols)
    {
      misc.error("Error: An internal error was happened", 0);
    }
    if(smt.nr != sm1.nr || smt.nc != sm1.nc)
    {
      misc.error("Error: An internal error was happened", 0);
    }
    if(this->uplo != 'B')
    {
      this->symmetrizeTriangularMatrix();
    }
    if(m1->uplo != 'B')
    {
      m1->symmetrizeTriangularMatrix();
    }
    if(this->uplo != 'B' || m1->uplo != 'B')
    {
      misc.error("Error: An internal error was happened", 0);
    }
    this->symmetric = false;
    this->uplo = 'B';
    
    
    ia = sm1.ir + 1;
    ja = sm1.ic + 1;
    ic = smt.ir + 1;
    jc = smt.ic + 1;
    m = sm1.nr;
    n = sm1.nc;
  }
  
  pdgeadd_(&t1, &m, &n, &alpha, m1->m, &ia, &ja, m1->descriptor, &beta, this->m, &ic, &jc, this->descriptor);
}

bool Matrix::symmetricInvert(double * logDeterminant, bool useSinglePrecision)
{
  int m1NGlobRows, m1NGlobCols, m1NBlockRows, m1NBlockCols;
  double logPartialDeterminant = 0.;
  
  if(this->symmetric == false || (this->uplo != 'U' && this->uplo != 'L' && this->uplo != 'B') )
  {
    misc.error("Error: An internal error was happened. symmetricInverse() expects a symmetric matrix", 0);
  }
  
  if( this->distribution == diagonalDistribution )
  {
    bool inverted = true;
    double logDeterminantTemp = 0;
    if(communicator->mpiRoot)
    {
      bool isNegativeDefinite = false;
      for(int i = 0; i<this->nGlobRows; i++)
      {
        logDeterminantTemp += log( fabs(this->m[i]) );
        if( this->m[i] == 0. )
        {
          std::stringstream ss;
          ss << i;
          misc.message << "Warning: Matrix inversion cannot be computed. The (" + ss.str()  +  ", " + ss.str()  + ") element is zero, and the inverse cannot be computed." << std::endl;
          inverted = false;
        }
        if( this->m[i] < 0. )
        {
          isNegativeDefinite = true;
        }
        this->m[i] = 1./this->m[i];
      }
      if( isNegativeDefinite == true )
      {
        misc.message << "Warning: The inverted matrix is not positive definite. The determinant is computed using diagonal absolute values." << std::endl;
      }
    }
    
    if( misc.gt(inverted == false) )
    {
      return false;
    }
    
    if(logDeterminant != NULL)
    {
      *logDeterminant = logDeterminantTemp;
      communicator->broadcast(logDeterminant, 1);
    }
  }
  else
  {
    //Initiate the parameters for the result matrix
    //this->initParameters(this->nGlobRows, this->nGlobCols, this->nBlockRows, this->nBlockCols);
    if(this->nGlobCols != this->nGlobRows)
    {
      misc.error("Error: An internal error has been happened while computing the inverse matrix. A square matrix is expected.", 0);
    }
    if(this->uplo == 'B')
    {
      this->uplo = 'L';
    }
    
    //char uplo = '?';
    int n = this->nGlobRows;
    //a, 
    int ia = 1;
    int ja = 1;
    //desca
    int info;
    
    //Compute Cholesky decomposition
    if(useSinglePrecision == false)
    {
      pdpotrf_(&this->uplo, &n, this->m, &ia, &ja, this->descriptor, &info);
    }
    else
    {
      copyDoubleToSingle();
      pspotrf_(&this->uplo, &n, this->mSinglePrecision, &ia, &ja, this->descriptor, &info);
    }
    if(info<0)
    {
      misc.error("Error: An internal error was happened (Illegal argument to pdpotrf function).", 0);
    }
    if(info>0){
      misc.message << "Warning: The Cholesky factorization cannot be completed because the matrix is not positive-definite." << std::endl;
      return false;
    }
    
    //Compute the determinant logarithm. I am not sure if this is the faster approach. This can be improved if block size is symmetric.
    bool local;
    if(logDeterminant != NULL)
    {
      logPartialDeterminant = 0.;
      *logDeterminant = 0.;
      for(int r=0; r<this->nGlobRows; r++)
      {
        double element = matrixElement(r, r, 1., &local, useSinglePrecision);
        logPartialDeterminant += log(element*element);
      }
      for (int pr = 0; pr < communicator->nProcRows; pr++) {
        for (int pc = 0; pc < communicator->nProcCols; pc++) {
          if (communicator->myRow == pr && communicator->myCol == pc) {
            //Send logPartialDeterminant to Root
            Cdgesd2d(communicator->context, 1, 1, &logPartialDeterminant, 1, 0, 0);
          }
          
          if (communicator->mpiRoot) {
            //Receive logPartialDeterminant
            double temp;
            Cdgerv2d(communicator->context, 1, 1, &temp, 1, pr, pc);
            *logDeterminant += temp;
          }
        }
      }
      communicator->barrier();
      communicator->broadcast(logDeterminant, 1);
    }
    
    
    //Compute the inverse matrix
    if(useSinglePrecision == false)
    {
      pdpotri_(&this->uplo, &n, this->m, &ia, &ja, this->descriptor, &info);
    }
    else
    {
      pspotri_(&this->uplo, &n, this->mSinglePrecision, &ia, &ja, this->descriptor, &info);
      copySingleToDouble();
    }
    if(info<0)
    {
      misc.error("Error: An internal error was happened (Illegal argument to pdpotri function).", 0);
    }
    if(info>0)
    {
      std::stringstream ss;
      ss << info;
      misc.message << "Warning: Matrix inversion cannot be computed. The (" + ss.str()  +  ", " + ss.str()  + ") element is zero, and the inverse cannot be computed." << std::endl;
      return false;
    }
  }
  
  return true;
}

bool Matrix::invert(double * logDeterminant)
{
  int m1NGlobRows, m1NGlobCols, m1NBlockRows, m1NBlockCols;
  double logPartialDeterminant = 0.;
  
  if( this->distribution == diagonalDistribution )
  {
    bool inverted = true;
    double logDeterminantTemp = 0;
    if(communicator->mpiRoot)
    {
      bool isNegativeDefinite = false;
      for(int i = 0; i<this->nGlobRows; i++)
      {
        logDeterminantTemp += log( fabs(this->m[i]) );
        if( this->m[i] == 0. )
        {
          std::stringstream ss;
          ss << i;
          misc.message << "Warning: Matrix inversion cannot be computed. The (" + ss.str()  +  ", " + ss.str()  + ") element is zero, and the inverse cannot be computed." << std::endl;
          inverted = false;
        }
        if( this->m[i] < 0. )
        {
          isNegativeDefinite = true;
        }
        this->m[i] = 1./this->m[i];
      }
      if( isNegativeDefinite == true )
      {
        misc.message << "Warning: The inverted matrix is not positive definite. The determinant is computed using diagonal absolute values." << std::endl;
      }
    }
    
    if( misc.gt(inverted == false) )
    {
      return false;
    }
    
    if(logDeterminant != NULL)
    {
      *logDeterminant = logDeterminantTemp;
      communicator->broadcast(logDeterminant, 1);
    }
  }
  else
  {
    if(this->uplo != 'B')
    {
      this->symmetrizeTriangularMatrix();
    }
    
    //Initiate the parameters for the result matrix
    //this->initParameters(this->nGlobRows, this->nGlobCols, this->nBlockRows, this->nBlockCols);
    if(this->nGlobCols != this->nGlobRows)
    {
      misc.error("Error: An internal error has been happened while computing the inverse matrix. A square matrix is expected.", 0);
    }
    
    //char uplo = '?';
    int m = this->nGlobRows;
    int n = this->nGlobCols;
    //a, 
    int ia = 1;
    int ja = 1;
    //desca
    int info;
    
    int ipiv[this->nRows + this->nBlockRows];
    
    //Compute LU factorization
    pdgetrf_(&m, &n, this->m, &ia, &ja, this->descriptor, ipiv, &info);
    if(info<0)
    {
      misc.error("Error: An internal error was happened (Illegal argument to pdgetrf function).", 0);
    }
    if(info>0){
      misc.message << "Warning: There are factors which are singular in the LU factorization." << std::endl;
      return false;
    }
    
    //Compute the determinant logarithm. I am not sure if this is the faster approach. This can be improved if block size is symmetric.
    bool local;
    if(logDeterminant != NULL)
    {
      logPartialDeterminant = 0.;
      *logDeterminant = 0.;
      for(int r=0; r<this->nGlobRows; r++)
      {
        double element = matrixElement(r, r, 1., &local); //std::cout << element << std::endl;
        logPartialDeterminant += log(fabs(element)); //Take the fabs() is correct???? This is the GCTA solution.... Maybe when the matrices are nearly positive definite??
      }
      for (int pr = 0; pr < communicator->nProcRows; pr++) {
        for (int pc = 0; pc < communicator->nProcCols; pc++) {
          if (communicator->myRow == pr && communicator->myCol == pc) {
            //Send logPartialDeterminant to Root
            Cdgesd2d(communicator->context, 1, 1, &logPartialDeterminant, 1, 0, 0);
          }
          
          if (communicator->mpiRoot) {
            //Receive logPartialDeterminant
            double temp;
            Cdgerv2d(communicator->context, 1, 1, &temp, 1, pr, pc);
            *logDeterminant += temp;
          }
        }
      }
      communicator->barrier();
      communicator->broadcast(logDeterminant, 1);
    }
    
    int lwork; //= this->nRows*this->nBlockCols;
    int liwork;
    /*
    //This uses the formula at the manual, but apparently do not work?
    if(communicator->nProcCols == communicator->nProcRows)
    {
      liwork = this->nCols + this->nBlockCols;
    }
    else
    {
      int lcm = leastCommonMultiple(communicator->nProcRows, communicator->nProcCols);
      double temp = ceil( double(ceil(double(this->nRows)/double(this->nBlockRows))) / double(lcm/communicator->nProcRows) );
      liwork = this->nCols + std::max( int(temp), this->nBlockCols );
    }
    std::cout << communicator->myRow << " " << communicator->myCol << " " << lwork << " " << liwork << std::endl;*/
    
    lwork = -1;
    liwork = -1;
    double qlwork;
    int qliwork;
    pdgetri_(&n, this->m, &ia, &ja, this->descriptor, ipiv, &qlwork, &lwork, &qliwork, &liwork, &info); //Query the lwork and liwork values.
    lwork = int(qlwork);
    liwork = int(qliwork);
    
    double *work = new double [lwork];
    int *iwork = new int [liwork];
    
    //Compute the inverse matrix
    pdgetri_(&n, this->m, &ia, &ja, this->descriptor, ipiv, work, &lwork, iwork, &liwork, &info);
    
    delete [] work;
    delete [] iwork;
    
    if(info<0)
    {
      misc.error("Error: An internal error was happened (Illegal argument to pdgetri function).", 0);
    }
    if(info>0)
    {
      std::stringstream ss;
      ss << info;
      misc.message << "Warning: Matrix inversion cannot be computed. The (" + ss.str()  +  ", " + ss.str()  + ") element is zero, and the inverse cannot be computed." << std::endl;
      return false;
    }
  }
  
  return true;
}

void Matrix::eigenDecomposition(Matrix * eigenValues, Matrix * eigenVectors)
{
  if(this->symmetric == false || (this->uplo != 'U' && this->uplo != 'L' && this->uplo != 'B') )
  {
    misc.error("Error: An internal error was happened. eigenDecomposition() expects a symmetric matrix", 0);
  }

  if(this->nGlobCols != this->nGlobRows)
  {
    misc.error("Error: An internal error has been happened while computing the eigen decomposition. A square matrix is expected.", 0);
  }
  
  if(this->distribution != cyclicDistribution)
  {
    misc.error("Error: Ops, eigendecomposition on a diagonal matrix is not implemented, yet. (Yes, it should be straightforward.)", 0);
  }
  
  if(this->uplo == 'B')
  {
    this->uplo = 'L';
  }
  
  double *globalEigenValues = new double [this->nGlobRows];
  eigenValues->initParameters(this->nGlobRows, this->nGlobCols, this->nBlockRows, this->nBlockCols, diagonalDistribution);
  eigenVectors->initParameters(this->nGlobRows, this->nGlobCols, this->nBlockRows, this->nBlockCols);
  
  char jobz = 'V';
  int n = this->nGlobRows;
  int ia = 1;
  int ja = 1;
  
  int iz = 1;
  int jz = 1;
  
  int info;
  
  int lwork = -1;
  double qlwork;
  
  //Query work size;
  pdsyev_(&jobz, &this->uplo, &n, this->m, &ia, &ja, this->descriptor, globalEigenValues, eigenVectors->m, &iz, &jz, eigenVectors->descriptor, &qlwork, &lwork, &info);
  
  lwork = int(qlwork);
  double * work = new double [lwork];
  
  //compute decomposition
  pdsyev_(&jobz, &this->uplo, &n, this->m, &ia, &ja, this->descriptor, globalEigenValues, eigenVectors->m, &iz, &jz, eigenVectors->descriptor, work, &lwork, &info);
  if(info<0)
  {
    misc.error("Error: An internal error was happened (Illegal argument to pdsyev function).", 0);
  }
  if(info>0)
  {
    misc.error("Error: The eigendecomposition cannot be performed. The computation did not converge or the system is not homogeneous.", 0);
  }
  
  eigenValues->setDiagonal(globalEigenValues, this->nGlobRows);
  
  //eigenValues->symmetric = false;
  //eigenValues->uplo = 'B';
  eigenVectors->symmetric = false;
  eigenVectors->uplo = 'B';
  
  delete [] globalEigenValues;
  delete [] work;
}

void Matrix::bendMatrix()
{
  if(this->symmetric == false || (this->uplo != 'U' && this->uplo != 'L' && this->uplo != 'B') )
  {
    misc.error("Error: An internal error was happened. bendMatrix() expects a symmetric matrix", 0);
  }
  
  if(this->nGlobCols != this->nGlobRows)
  {
    misc.error("Error: An internal error has been happened while bending a matrix. A square matrix is expected.", 0);
  }
  
  if(this->distribution != cyclicDistribution)
  {
    misc.error("Error: Matrix bending is not implemented on non-cyclic distributed matrices.)", 0);
  }
  
  Matrix * eigenValues = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Matrix * eigenVectors = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  
  this->eigenDecomposition(eigenValues, eigenVectors);
  
  std::vector<double> globalEigenValues = eigenValues->diagonal();
  if(communicator->mpiRoot)
  {
    double evalMean = 0.;
    double S = 0.;
    double P = 0.;
    bool flag = false;
    for(int i = 0; i < eigenValues->nGlobRows; i++)
    {
      evalMean += globalEigenValues[i];
      if(globalEigenValues[i] > 0)
      {
        continue;
      }
      if(globalEigenValues[i] == 0.)
      {
        misc.error("Error: Bending a singular matrix. The results could be wrong. Stopping.", 0);
      }
      S += globalEigenValues[i];
      P = ( (P<globalEigenValues[i])?P:globalEigenValues[i] );
      flag = true;
    }
    evalMean /= double(eigenValues->nGlobRows);
    if(flag == true)
    {
      double W = S*S*100.+1.;
      double correctedEvalMean = 0.;
      for(int i=0; i < eigenValues->nGlobRows; i++)
      {
        if(globalEigenValues[i] >= 0)
        {
          correctedEvalMean += globalEigenValues[i];
          continue;
        }
        globalEigenValues[i] = -P*( S - globalEigenValues[i] )*( S - globalEigenValues[i] )/W;
        correctedEvalMean += globalEigenValues[i];
      }
      correctedEvalMean /= double(eigenValues->nGlobRows);
      for(int i=0; i < eigenValues->nGlobRows; i++)
      {
        globalEigenValues[i] *= evalMean/correctedEvalMean;
        if(globalEigenValues[i] < 0.)
        {
          misc.error("Error: An error was happened when bending a matrix. Not all eigenvalues have been converted to positive values.", 0);
        }
      }
    }
  }
  eigenVectors->scatterVector(&(globalEigenValues[0]), column);
  
  #pragma omp parallel for
  for(int c = 0; c<eigenVectors->nCols; c++)
  {
    eigenVectors->v[c] = sqrt(eigenVectors->v[c]);
  }
  
  #pragma omp parallel for
  for(int c = 0; c<eigenVectors->nCols; c++)
  {
    for(int r = 0; r<eigenVectors->nRows; r++)
    {
      eigenVectors->m[c*eigenVectors->nRows + r] *= eigenVectors->v[c];
    }
  }
  
  this->multiply(eigenVectors, 'N', eigenVectors, 'T');
  
  delete eigenValues;
  delete eigenVectors;
}

void Matrix::QRDecomposition()
{
  if(this->distribution != cyclicDistribution)
  {
    misc.error("Error: Ops, QR on a diagonally distributed matrix is not implemented, yet.", 0);
  }
  
  if(this->symmetric == true)
  {
    symmetrizeTriangularMatrix();
  }

  int m = this->nGlobRows;  
  int n = this->nGlobCols;
  int ia = 1;
  int ja = 1;
  
  int info;

  int tauSize = (this->nGlobRows < this->nGlobCols)?this->nGlobRows:this->nGlobCols;
  double * tau = new double [tauSize];
  
  int lwork = -1;
  double qlwork;
  
  //Query work size;
  pdgeqrf_(&m, &n, this->m, &ia, &ja, this->descriptor, tau, &qlwork, &lwork, &info);
  
  lwork = int(qlwork);
  double * work = new double [lwork];
  
  //compute decomposition
  pdgeqrf_(&m, &n, this->m, &ia, &ja, this->descriptor, tau, work, &lwork, &info);
  if(info<0)
  {
    misc.error("Error: An internal error was happened (Illegal argument to pdgeqrf function).", 0);
  }
  if(info>0)
  {
    misc.error("Error: The QR decomposition cannot be performed.", 0);
  }
  
  this->symmetric = false;
  this->uplo = 'B';
  
  delete [] work;
  delete [] tau;
}

std::vector<int> Matrix::getDependentColumns(double threshold)
{
  if(this->distribution != cyclicDistribution)
  {
    misc.error("Error: Ops, Getting independent columns on a diagonally distributed matrix is not implemented, yet.", 0);
  }
  if(this->nGlobRows < this->nGlobCols)
  {
    misc.error("Error: Ops, Getting independent columns on a matrix with more columns than rows is not implemented, yet. (It should be straightforward at some point.)", 0);
  }
  
  QRDecomposition();
  
  std::vector<double> RDiagonal = diagonal(true);
  std::vector<int> dependentColumns;
  
  if(communicator->mpiRoot == true)
  {
    if(RDiagonal.size() != this->nGlobCols)
    {
      misc.error("Error: An internal error was happened when searching for independent columns. Diagonal and matrix dimensions disagree.", 0);
    }
    for(int i = 0; i<RDiagonal.size(); i++)
    {
      if( fabs(RDiagonal[i]) < threshold )
      {
        dependentColumns.push_back(i);
      }
    }
  }
  
  return dependentColumns;
}

double Matrix::trace()
{
  double totalTrace = 0.;
  double partialTrace = 0.;
  
  if( this->distribution == diagonalDistribution)
  {
    if( communicator->mpiRoot )
    {
      for(int i = 0; i<this->nGlobRows; i++)
      {
        totalTrace += this->m[i];
      }
    }
    communicator->broadcast(&totalTrace, 1);
  }
  else
  {
    //Compute the trace. I am not sure if this is the faster approach. This can be improved if block size is symmetric.
    bool local;
    for(int r=0; r<this->nGlobRows; r++)
    {
      double element = matrixElement(r, r, 0., &local);
      partialTrace += element;
    }
    
    for (int pr = 0; pr < communicator->nProcRows; pr++) {
      for (int pc = 0; pc < communicator->nProcCols; pc++) {
        if (communicator->myRow == pr && communicator->myCol == pc) {
          //Send partialTrace to Root
          Cdgesd2d(communicator->context, 1, 1, &partialTrace, 1, 0, 0);
        }
        
        if (communicator->mpiRoot) {
          //Receive partialTrace
          double temp;
          Cdgerv2d(communicator->context, 1, 1, &temp, 1, pr, pc);
          totalTrace += temp;
        }
      }
    }
    communicator->barrier();
    communicator->broadcast(&totalTrace, 1);
  }

  return totalTrace;
}

double Matrix::elementsAverage()
{
  if( this->symmetric == true )
  {
    this->symmetrizeTriangularMatrix();
  }
  
  double totalAverage = 0.;
  double partialAverage = 0.;
  
  if( this->distribution == diagonalDistribution)
  {
    if( communicator->mpiRoot )
    {
      for(int i = 0; i<this->nGlobRows; i++)
      {
        totalAverage += this->m[i];
      }
    }
    totalAverage = totalAverage/(double(this->nGlobRows)*double(this->nGlobCols));
    communicator->broadcast(&totalAverage, 1);
  }
  else
  {
    for(int c = 0; c<this->nCols; c++)
    {
      for(int r = 0; r<this->nRows; r++)
      {
        partialAverage += this->m[c*this->nRows + r]/double(this->nGlobRows);
      }
    }
    
    partialAverage /= double(this->nGlobCols);
    double * vectorPartialAverages = communicator->gather(&partialAverage, 1);
    
    if (communicator->mpiRoot)
    {
      for(int ip = 0; ip < communicator->mpiNumTasks; ip++)
      {
        totalAverage += vectorPartialAverages[ip];
      }
      delete [] vectorPartialAverages;
    }
    communicator->broadcast(&totalAverage, 1);
  }

  return totalAverage;
}

void Matrix::applyExponentialOperator(double alpha)
{
  if(this->distribution == diagonalDistribution)
  {
    std::vector<double> diagonalElements = diagonal();
    if(communicator->mpiRoot == true)
    {
      for( int i = 0; i < diagonalElements.size(); i++ )
      {
        diagonalElements[i] = exp(alpha*diagonalElements[i]);
      }
    }
    initParameters(this->nGlobRows, this->nGlobCols, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols, cyclicDistribution);
    fillWithConstant(1.);
    setDiagonal(&(diagonalElements[0]), this->nGlobRows);
  }
  else
  {
    #pragma omp parallel for
    for(int c = 0; c<this->nCols; c++)
    {
      for(int r = 0; r<this->nRows; r++)
      {
        this->m[c*this->nRows + r] = exp(alpha*this->m[c*this->nRows + r]);
      }
    }
  }
}

void Matrix::applyInverseLogistic()
{
  if(this->distribution == diagonalDistribution)
  {
    std::vector<double> diagonalElements = diagonal();
    if(communicator->mpiRoot == true)
    {
      for( int i = 0; i < diagonalElements.size(); i++ )
      {
        diagonalElements[i] = 1./(1. + exp( -diagonalElements[i] ));
      }
    }
    initParameters(this->nGlobRows, this->nGlobCols, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols, cyclicDistribution);
    fillWithConstant(0.5);
    setDiagonal(&(diagonalElements[0]), this->nGlobRows);
  }
  else
  {
    #pragma omp parallel for
    for(int c = 0; c<this->nCols; c++)
    {
      for(int r = 0; r<this->nRows; r++)
      {
        this->m[c*this->nRows + r] = 1./(1. + exp( -this->m[c*this->nRows + r] ));
      }
    }
  }
}

std::vector<double> Matrix::diagonal(bool allowNonSquareMatrices)
{
  if(this->nGlobRows != this->nGlobCols && allowNonSquareMatrices == false)
  {
    misc.error("Error: An internal error was happened. Diagonal entries can only be retrieved from square matrices.", 0);
  }
  if(this->nBlockRows != this->nBlockCols)
  {
    misc.error("Error: An internal error was happened. Diagonal entries can only be retrieved from matrices with symmetric block sizes.", 0);
  }

  std::vector<double> result;
  if( this->distribution == diagonalDistribution)
  {
    if( communicator->mpiRoot )
    {
      for(int i = 0; i<this->nGlobRows; i++)
      {
        result.push_back( this->m[i] );
      }
    }
  }
  else
  {
    double * mGlobalBlock;
    if(communicator->mpiRoot)
    {
      mGlobalBlock = new double [this->nBlockRows*this->nBlockCols];
    }
    
    int minimumDimension = 0;
    if(this->nGlobCols < this->nGlobRows)
    {
      minimumDimension = this->nGlobCols;
    }
    else
    {
      minimumDimension = this->nGlobRows;
    }
    
    for (int r = 0; r < minimumDimension; r += this->nBlockRows) {
      gatherBlock(mGlobalBlock, r, r, this->nBlockRows);
      
      int nr = this->nBlockRows;
      if (minimumDimension-r < this->nBlockRows)
      {
        nr = minimumDimension-r;
      }
      
      if(communicator->mpiRoot)
      {
        for(int id = 0; id < nr; id++)
        {
          result.push_back(mGlobalBlock[id*this->nBlockRows + id]);
        }
      }
    }
    
    if(communicator->mpiRoot)
    {
      delete [] mGlobalBlock;
    }
  }
  
  return result;
}

void Matrix::setDiagonal(double * diagonal, int nElements)
{
  if(this->nGlobRows != this->nGlobCols)
  {
    misc.error("Error: An internal error was happened. Diagonal entries can only be set on square matrices.", 0);
  }
  if(communicator->mpiRoot == true && (this->nGlobRows != nElements || this->nGlobCols != nElements) )
  {
    misc.error("Error: An internal error was happened. Number of diagonal entries differ from matrix dimensions.", 0);
  }
  
  if(this->distribution == diagonalDistribution)
  {
    if(communicator->mpiRoot == true)
    {
      for (int r = 0; r < this->nRows; r += 1)
      {
        this->m[r] = diagonal[r];
      }
    }
  }
  else if(this->distribution == cyclicDistribution)
  {
    communicator->broadcast(&nElements);
    if(communicator->mpiRoot == false)
    {
      diagonal = new double[nElements];
    }
    communicator->broadcast(diagonal, nElements);
    
    for (int r = 0; r < this->nGlobRows; r += 1)
    {
      LocalPosition lrp = global2local(r, this->nBlockRows, communicator->nProcRows);
      LocalPosition lcp = global2local(r, this->nBlockCols, communicator->nProcCols);
      if( communicator->myRow == lrp.proc && communicator->myCol == lcp.proc)
      {
        if(lrp.position >= this->nRows || lcp.position >= this->nCols)
        {
          misc.error("Error: An internal error was happened. Indices out of rank when setting the diagonal on a distributed matrix.", 0);
        }
        this->m[lcp.position*this->nRows + lrp.position] = diagonal[r];
      }
    }
    
    if(communicator->mpiRoot == false)
    {
      delete [] diagonal;
    }
  }
  else
  {
    misc.error("Error: An internal error was happened. Matrix diagonal cannot be set on a matrix with this distribution.", 0);
  }
}

double Matrix::traceOfMatrixProduct(Matrix * m1)
{
  if( this->nGlobRows != this->nGlobCols || m1->nGlobRows != m1->nGlobCols )
  {
    misc.error("Error: An internal error was happened: Trying to compute the trace of a product of a two non-square matrices.", 0);
  }
  if( this->nGlobRows != m1->nGlobRows || this->nGlobCols != m1->nGlobCols )
  {
    misc.error("Error: An internal error was happened: Trying to compute the trace of a product of a two matrices with different dimensions.", 0);
  }

  double trace = 0.;
  if( this->distribution == diagonalDistribution || m1->distribution == diagonalDistribution )
  {
    std::vector<double> diagonal1 = this->diagonal();
    std::vector<double> diagonal2 = m1->diagonal();
    if( diagonal1.size() != diagonal2.size() )
    {
      misc.error("Error: An internal error was happened when computing the trace of a matrix product.", 0);
    }
    if( communicator->mpiRoot )
    {
      for(int i = 0; i<diagonal1.size(); i++)
      {
        trace += diagonal1[i]*diagonal2[i];
      }
    }
    
    communicator->broadcast(&trace, 1);
  }
  else
  {
    if( this->symmetric == true )
    {
      this->symmetrizeTriangularMatrix();
    }
    Matrix * m1t;
    if( m1->symmetric == true )
    {
      m1->symmetrizeTriangularMatrix();
      m1t = m1;
    }
    else
    {
      m1t = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
      m1t->transpose(m1);
    }
    
    checkMatrixStructure(m1t);

    if( (this->nBlockRows != m1t->nBlockRows) || (this->nBlockCols != m1t->nBlockCols) || (this->nRows != m1t->nRows) || (this->nCols != m1t->nCols))
    {
      misc.error("Error: An internal error was happened. The trace of a matrix product can not be performed on matrices with different structural properties.", 0);
    }
    
    double sum = 0.;
    for(int c = 0; c<this->nCols; c++)
    {
      for(int r = 0; r<this->nRows; r++)
      {
        sum += this->m[c*this->nRows + r] * m1t->m[c*this->nRows + r];
      }
    }
    
    double * result = communicator->gather(&sum, 1);
    if(communicator->mpiRoot)
    {
      for(int i = 0; i<communicator->mpiNumTasks; i++)
      {
        trace += result[i];
      }
    }
    delete [] result;
    
    communicator->broadcast(&trace);
    
    if( m1t != m1 )
    {
      delete m1t;
    }
  }
  
  return trace;
}

double Matrix::diagonalOfABAt(Matrix * A, Matrix * B)
{
  if(B->nGlobCols != B->nGlobRows)
  {
    misc.error("Error: An internal error was happened. On function diagonalOfABAt, B has to be a square matrix.", 0);
  }
  
  Matrix * temp = new Matrix();
  temp->multiply(A, 'N', B, 'N');
  
  if(A->symmetric == true)
  {
    A->symmetrizeTriangularMatrix();
  }
  if(temp->symmetric == true)
  {
    temp->symmetrizeTriangularMatrix();
  }
  
  temp->elementWiseMultiplication(A);
  
  Matrix * tempColumnWithOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, A->nGlobCols, 1);
  tempColumnWithOnes->fillWithConstant(1.);
  
  Matrix * temp2 = new Matrix();
  temp2->multiply(temp, 'N', tempColumnWithOnes, 'N');
  
  double * temp2Global;
  if( communicator->mpiRoot == true )
  {
    temp2Global = new double [A->nGlobRows];
  }
  temp2->gatherMatrix(temp2Global);
  
  initParameters(A->nGlobRows, A->nGlobRows, A->nBlockRows, A->nBlockCols, diagonalDistribution);
  setDiagonal(temp2Global, A->nGlobRows);
  
  delete temp;
  delete temp2;
  delete tempColumnWithOnes;
  if( communicator->mpiRoot == true )
  {
    delete [] temp2Global;
  }
}

void Matrix::transpose(Matrix * m1)
{
  if(m1->distribution == diagonalDistribution)
  {
    if(m1->nGlobRows != m1->nGlobCols)
    {
      misc.error("Error: An internal error was happened when transposing a diagonal Matrix.", 0);
    }
    this->duplicateMatrix(m1);
    return;
  }
  
  this->initParameters(m1->nGlobCols, m1->nGlobRows, m1->nBlockRows, m1->nBlockCols);
  
  int m = m1->nGlobCols;
  int n = m1->nGlobRows;
  double alpha = 1.;
  int ia = 1;
  int ja = 1;
  double beta = 0.;
  int ic = 1;
  int jc = 1;
  
  pdtran_(&m, &n, &alpha, m1->m, &ia, &ja, m1->descriptor, &beta, this->m, &ic, &jc, this->descriptor);
  
  this->symmetric = m1->symmetric;
  this->vector = m1->vector;
  if(m1->uplo == 'B')
  {
    this->uplo = 'B';
  }
  else if(m1->uplo == 'L' && m1->symmetric == true)
  {
    this->uplo = 'U';
  }
  else if(m1->uplo == 'U' && m1->symmetric == true)
  {
    this->uplo = 'L';
  }
  else
  {
    misc.error("Error: An internal error was happened when transposing a matrix.", 0);
  }
}

void Matrix::symmetrizeTriangularMatrix()
{
  if(this->nGlobRows != this->nGlobCols || this->symmetric == false)
  {
    misc.error("Error: An internal error was happened. A not square or not triangular/symmetric matrix cannot be symmetrized.", 0);
  }
  if(this->distribution == diagonalDistribution)
  {
    return;
  }
  if(this->uplo == 'B')
  {
    return;
  }
  
  Matrix * temp = new Matrix(this->distribution);
  temp->transpose(this);
  
  int m = temp->nGlobRows;
  int n = temp->nGlobCols;
  int ia = 1;
  int ja = 1;
  int ib = 1;
  int jb = 1;
  
  if(temp->nRows != this->nRows || temp->nCols != this->nCols || temp->nBlockRows != this->nBlockRows || temp->nBlockCols != this->nBlockCols || temp->nGlobRows != this->nGlobRows || temp->nGlobCols != this->nGlobCols)
  {
    misc.error("Error: An internal error was happened. Invalid matrix specifications when symmetrizing matrix.", 0);
  }
  
  pdlacpy_(&(temp->uplo), &m, &n, temp->m, &ia, &ja, temp->descriptor, this->m, &ib, &jb, this->descriptor);
  this->uplo = 'B';
  
  delete temp;
}

void Matrix::elementWiseMultiplication(Matrix * m1, double scale)
{
  if( (this->nGlobRows != m1->nGlobRows) || (this->nGlobCols != m1->nGlobCols) || (this->nRows != m1->nRows) || (this->nCols != m1->nCols) || (this->distribution != m1->distribution))
  {
    misc.error("Error: An internal error was happened. Element wise multiplication can not be performed on matrices with different structural properties.", 0);
  }
  if( (this->nBlockRows != m1->nBlockRows) || (this->nBlockCols != m1->nBlockCols))
  {
    misc.error("Error: An internal error was happened. Element wise multiplication can not be performed on matrices with different structural properties.", 0);
  }
  if(this->symmetric != m1->symmetric || this->uplo != m1->uplo)
  {
    if( (this->symmetric == true && this->uplo != 'B') || (m1->symmetric == true && m1->uplo != 'B') )
    {
      misc.error("Error: An internal error was happened. Element wise multiplication can not be performed on matrices with different symmetry properties.", 0);
    }
  }
  #pragma omp parallel for
  for(int c = 0; c<this->nCols; c++)
  {
    for(int r = 0; r<this->nRows; r++)
    {
      this->m[c*this->nRows + r] *= scale*m1->m[c*this->nRows + r];
    }
  }
}

void Matrix::elementWiseDivision(Matrix * m1, double scale)
{
  if( (this->nGlobRows != m1->nGlobRows) || (this->nGlobCols != m1->nGlobCols) || (this->nRows != m1->nRows) || (this->nCols != m1->nCols) || (this->distribution != m1->distribution))
  {
    misc.error("Error: An internal error was happened. Element wise division can not be performed on matrices with different structural properties.", 0);
  }
  if( (this->nBlockRows != m1->nBlockRows) || (this->nBlockCols != m1->nBlockCols))
  {
    misc.error("Error: An internal error was happened. Element wise division can not be performed on matrices with different structural properties.", 0);
  }
  if(this->symmetric != m1->symmetric || this->uplo != m1->uplo)
  {
    if( (this->symmetric == true && this->uplo != 'B') || (m1->symmetric == true && m1->uplo != 'B') )
    {
      misc.error("Error: An internal error was happened. Element wise division can not be performed on matrices with different symmetry properties.", 0);
    }
  }
  #pragma omp parallel for
  for(int c = 0; c<this->nCols; c++)
  {
    for(int r = 0; r<this->nRows; r++)
    {
      this->m[c*this->nRows + r] *= scale / m1->m[c*this->nRows + r];
    }
  }
}

void Matrix::scaleBy(double scale)
{
  #pragma omp parallel for
  for(int c = 0; c<this->nCols; c++)
  {
    for(int r = 0; r<this->nRows; r++)
    {
      this->m[c*this->nRows + r] *= scale;
    }
  }
}

void Matrix::makeIntersectionMatrix(std::vector<int> categories, double onIntersection, double onNoIntersection)
{
  int dimension;
  
  dimension = categories.size();
  communicator->broadcast(&dimension, 1);
  
  initParameters(dimension, dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  fillWithConstant(onNoIntersection);
  
  int * rv = scatterVectorRet(&(categories[0]), row);
  int * cv = scatterVectorRet(&(categories[0]), column);

  for(int c = 0; c < nCols; c++)
  {
    for(int r = 0; r < nRows; r++)
    {
      this->m[c*this->nRows + r] = ((rv[r]==cv[c])?onIntersection:onNoIntersection);
    }
  }
  
  this->symmetric = true;
  this->uplo = 'B';
  
  delete [] rv;
  delete [] cv;
}

void Matrix::makeDifferenceMatrix(std::vector<double> v1, std::vector<double> v2)
{
  if(communicator->mpiRoot == true && (v1.size() == 0 || v2.size() == 0))
  {
    misc.error("Error: An internal error was happened. The makeDifferenceMatrix() method requires two non-empty vectors.", 0);
  }
  int newGlobalRows = v1.size();
  int newGlobalCols = v2.size();
  communicator->broadcast(&newGlobalRows);
  communicator->broadcast(&newGlobalCols);
  
  initParameters(newGlobalRows, newGlobalCols, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  scatterVector(&(v1[0]), row);
  double *dv1 = this->v;
  this->v = NULL;
  scatterVector(&(v2[0]), column);
  double *dv2 = this->v;
  this->v = NULL;
  
  #pragma omp parallel for
  for(int c = 0; c<this->nCols; c++)
  {
    for(int r = 0; r<this->nRows; r++)
    {
      this->m[c*this->nRows + r] = dv1[r] - dv2[c];
    }
  }
  
  delete [] dv1;
  delete [] dv2;
}

void Matrix::copyDoubleToSingle()
{
  if(this->mSinglePrecision == NULL)
  {
    this->mSinglePrecision = new(std::nothrow) float[this->nRows*this->nCols];
    if(this->mSinglePrecision == NULL)
    {
      misc.error("Error: Sorry, unable to allocate enough memory.", 0);
    }
    
    misc.estimateMaxMemory((double(this->nGlobRows)*double(this->nGlobCols)*4.));
  }
  #pragma omp parallel for
  for(int c = 0; c<this->nCols; c++)
  {
    for(int r = 0; r<this->nRows; r++)
    {
      this->mSinglePrecision[c*this->nRows + r] = float(this->m[c*this->nRows + r]);
    }
  }
}

void Matrix::copySingleToDouble(bool unallocateSinglePrecisionMatrix)
{
  if(this->mSinglePrecision == NULL)
  {
    misc.error("Error: An internal error was happened. There is not a single precission version of the matrix.", 0);
  }
  #pragma omp parallel for
  for(int c = 0; c<this->nCols; c++)
  {
    for(int r = 0; r<this->nRows; r++)
    {
      this->m[c*this->nRows + r] = double(this->mSinglePrecision[c*this->nRows + r]);
    }
  }
  
  if(unallocateSinglePrecisionMatrix == true)
  {
    delete [] this->mSinglePrecision;
    this->mSinglePrecision = NULL;
    misc.estimateMaxMemory(-(double(this->nGlobRows)*double(this->nGlobCols)*4.));
  }
}

void Matrix::getGlobalIndexElementsGreatherThan(double threshold, std::vector<int> & rows, std::vector<int> & cols)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: Get global index elements greather than a threshold is not allowed on a non-distributed matrix such as a diagonal matrix. If you need this function, contact us.", 0);
  }
  if( this->symmetric == true && ( this->uplo != 'B' && this->uplo != 'L' && this->uplo != 'U' ) )
  {
    misc.error("Error: An internal error was happened. Symetric matrix with invalid properties in getGlobalIndexElementsGreatherThan() method.", 0);
  }
  
  rows.clear();
  cols.clear();
  for(int c = 0; c < this->nCols; c++)
  {
    for(int r = 0; r < this->nRows; r++)
    {
      if(this->m[c*this->nRows + r] > threshold)
      {
        int gr = local2global(row, r);
        int gc = local2global(column, c);
        if( this->symmetric == true && this->uplo != 'B' )
        {
          if( ( this->uplo == 'L' && gc <= gr ) || ( this->uplo == 'U' && gr <= gc ) )
          {
            rows.push_back(gr);
            cols.push_back(gc);
            
            if(gc != gr)
            {
              rows.push_back(gc);
              cols.push_back(gr);
            }
          }
        }
        else
        {
          rows.push_back(gr);
          cols.push_back(gc);
        }
      }
    }
  }
}

void Matrix::getGlobalIndexOutsideRange(double lowerThreshold, double upperThreshold, std::vector<int> & idxGlobalRows, std::vector<int> & idxGlobalCols)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: Get global index elements greather than a threshold is not allowed on a non-distributed matrix such as a diagonal matrix. If you need this function, contact us.", 0);
  }
  if( this->symmetric == true && ( this->uplo != 'B' && this->uplo != 'L' && this->uplo != 'U' ) )
  {
    misc.error("Error: An internal error was happened. Symetric matrix with invalid properties in getGlobalIndexElementsGreatherThan() method.", 0);
  }
  
  std::vector<int> rows;
  std::vector<int> cols;

  for(int c = 0; c < this->nCols; c++)
  {
    for(int r = 0; r < this->nRows; r++)
    {
      if(this->m[c*this->nRows + r] < lowerThreshold || this->m[c*this->nRows + r] > upperThreshold)
      {
        int gr = local2global(row, r);
        int gc = local2global(column, c);
        if( this->symmetric == true && this->uplo != 'B' )
        {
          if( ( this->uplo == 'L' && gc <= gr ) || ( this->uplo == 'U' && gr <= gc ) )
          {
            rows.push_back(gr);
            cols.push_back(gc);
            
            if(gc != gr)
            {
              rows.push_back(gc);
              cols.push_back(gr);
            }
          }
        }
        else
        {
          rows.push_back(gr);
          cols.push_back(gc);
        }
      }
    }
  }

  //Combine the local data to global arrays in the root process.
  int * globalRows;
  int * globalCols;
  int nGlobalRows;
  int nGlobalCols;
  globalRows = communicator->asymmetricGather(&(rows[0]), rows.size(), &nGlobalRows);
  globalCols = communicator->asymmetricGather(&(cols[0]), cols.size(), &nGlobalCols);
  if(nGlobalRows != nGlobalCols)
  {
    misc.error("Error: An internal error was happened when searching for matrix elements outside a range.", 0);
  }

  idxGlobalRows.clear();
  idxGlobalCols.clear();
  if(communicator->mpiRoot == true)
  {
    idxGlobalRows.assign(globalRows, globalRows + nGlobalRows);
    idxGlobalCols.assign(globalCols, globalCols + nGlobalCols);
  }
  
  delete [] globalRows;
  delete [] globalCols;
}

void Matrix::getGlobalIndexInsideRange(double lowerThreshold, double upperThreshold, std::vector<int> & idxGlobalRows, std::vector<int> & idxGlobalCols)
{
  if( this->distribution == global || this->distribution == diagonalDistribution )
  {
    misc.error("Error: Get global index elements greather than a threshold is not allowed on a non-distributed matrix such as a diagonal matrix. If you need this function, contact us.", 0);
  }
  if( this->symmetric == true && ( this->uplo != 'B' && this->uplo != 'L' && this->uplo != 'U' ) )
  {
    misc.error("Error: An internal error was happened. Symetric matrix with invalid properties in getGlobalIndexElementsGreatherThan() method.", 0);
  }
  
  std::vector<int> rows;
  std::vector<int> cols;

  for(int c = 0; c < this->nCols; c++)
  {
    for(int r = 0; r < this->nRows; r++)
    {
      if(this->m[c*this->nRows + r] > lowerThreshold && this->m[c*this->nRows + r] < upperThreshold)
      {
        int gr = local2global(row, r);
        int gc = local2global(column, c);
        if( this->symmetric == true && this->uplo != 'B' )
        {
          if( ( this->uplo == 'L' && gc <= gr ) || ( this->uplo == 'U' && gr <= gc ) )
          {
            rows.push_back(gr);
            cols.push_back(gc);
            
            if(gc != gr)
            {
              rows.push_back(gc);
              cols.push_back(gr);
            }
          }
        }
        else
        {
          rows.push_back(gr);
          cols.push_back(gc);
        }
      }
    }
  }

  //Combine the local data to global arrays in the root process.
  int * globalRows;
  int * globalCols;
  int nGlobalRows;
  int nGlobalCols;
  globalRows = communicator->asymmetricGather(&(rows[0]), rows.size(), &nGlobalRows);
  globalCols = communicator->asymmetricGather(&(cols[0]), cols.size(), &nGlobalCols);
  if(nGlobalRows != nGlobalCols)
  {
    misc.error("Error: An internal error was happened when searching for matrix elements outside a range.", 0);
  }

  idxGlobalRows.clear();
  idxGlobalCols.clear();
  if(communicator->mpiRoot == true)
  {
    idxGlobalRows.assign(globalRows, globalRows + nGlobalRows);
    idxGlobalCols.assign(globalCols, globalCols + nGlobalCols);
  }
  
  delete [] globalRows;
  delete [] globalCols;
}

void Matrix::joinMatrices(Matrix * m1, subMatrix sm1, Matrix * m2, subMatrix sm2, double backgroundValue)
{
  if(m1->distribution == diagonalDistribution && m2->distribution == diagonalDistribution && backgroundValue == 0. && sm1.ir == sm1.ic && sm2.ir == sm2.ic)
  {
    joinDiagonalMatrices(m1, sm1, m2, sm2, backgroundValue);
  }
  else
  {
    Matrix * m1Temp;
    Matrix * m2Temp;
    if(m1->distribution == diagonalDistribution)
    {
      m1Temp = new Matrix(cyclicDistribution, m1->nGlobRows, m1->nGlobCols);
      m1Temp->fillWithConstant(backgroundValue);
      m1Temp->setDiagonal(m1->m, m1->nGlobRows);
    }
    else
    {
      m1Temp = m1;
    }
    if(m2->distribution == diagonalDistribution)
    {
      m2Temp = new Matrix(cyclicDistribution, m2->nGlobRows, m2->nGlobCols);
      m2Temp->fillWithConstant(backgroundValue);
      m2Temp->setDiagonal(m2->m, m2->nGlobRows);
    }
    else
    {
      m2Temp = m2;
    }
    
    joinGeneralMatrices(m1Temp, sm1, m2Temp, sm2, backgroundValue);
    
    if(m1Temp != m1)
    {
      delete m1Temp;
    }
    if(m2Temp != m2)
    {
      delete m2Temp;
    }
  }
}

void Matrix::joinDiagonalMatrices(Matrix * m1, subMatrix sm1, Matrix * m2, subMatrix sm2, double backgroundValue)
{
  if(m1->distribution != diagonalDistribution || m2->distribution != diagonalDistribution || backgroundValue != 0. || sm1.ir != sm1.ic || sm2.ir != sm2.ic)
  {
    misc.error("Error: An internal error was happened. Joining diagonal matrices can not be performed using these matrices.", 0);
  }
  if(sm1.nr != m1->nGlobRows || sm1.nc != m1->nGlobCols)
  {
    misc.error("Error: An internal error was happened. subMatrix 1 is not properly defined.", 0);
  }
  if(sm2.nr != m2->nGlobRows || sm2.nc != m2->nGlobCols)
  {
    misc.error("Error: An internal error was happened. subMatrix 2 is not properly defined.", 0);
  }

  int nNewGlobRows = (sm1.ir+sm1.nr)>=(sm2.ir+sm2.nr)?(sm1.ir+sm1.nr):(sm2.ir+sm2.nr);
  int nNewGlobCols = (sm1.ic+sm1.nc)>=(sm2.ic+sm2.nc)?(sm1.ic+sm1.nc):(sm2.ic+sm2.nc);
  
  initParameters(nNewGlobRows, nNewGlobCols, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols, diagonalDistribution);
  fillWithConstant(backgroundValue);
  
  if( communicator->mpiRoot == true )
  {
    for(int i = 0; i<sm1.nr; i++)
    {
      this->m[sm1.ir + i] = m1->m[i];
    }
    for(int i = 0; i<sm2.nr; i++)
    {
      this->m[sm2.ir + i] = m2->m[i];
    }
  }
}

void Matrix::joinGeneralMatrices(Matrix * m1, subMatrix sm1, Matrix * m2, subMatrix sm2, double backgroundValue)
{
  int nNewGlobRows;
  int nNewGlobCols;
  
  if(m1->distribution != cyclicDistribution)
  {
    misc.error("Error: An internal error was happened. General matrix joining connot be performed on diagonal matrices.", 0);
  }
  if(sm1.active == false || sm2.active == false)
  {
    misc.error("Error: An internal error was happened. Submatrices are no defined.", 0);
  }
  if(sm1.nr != m1->nGlobRows || sm1.nc != m1->nGlobCols)
  {
    misc.error("Error: An internal error was happened. subMatrix 1 is not properly defined.", 0);
  }
  if(sm2.nr != m2->nGlobRows || sm2.nc != m2->nGlobCols)
  {
    misc.error("Error: An internal error was happened. subMatrix 2 is not properly defined.", 0);
  }
  
  nNewGlobRows = (sm1.ir+sm1.nr)>=(sm2.ir+sm2.nr)?(sm1.ir+sm1.nr):(sm2.ir+sm2.nr);
  nNewGlobCols = (sm1.ic+sm1.nc)>=(sm2.ic+sm2.nc)?(sm1.ic+sm1.nc):(sm2.ic+sm2.nc);
  
  initParameters(nNewGlobRows, nNewGlobCols, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  fillWithConstant(backgroundValue);
  
  add(m1, 0., 1., sm1, subMatrix(m1));
  add(m2, 0., 1., sm2, subMatrix(m2));
}

void Matrix::joinMatricesVertically(Matrix * m1, Matrix * m2)
{
  if(m1->nGlobCols != m2->nGlobCols)
  {
    misc.error("Error: An internal error was happened. Matrices with different number of columns can not be joined vertically.", 0);
  }
  joinMatrices(m1, subMatrix(0, 0, m1->nGlobRows, m1->nGlobCols), m2, subMatrix(m1->nGlobRows, 0, m2->nGlobRows, m2->nGlobCols));
}

void Matrix::joinMatricesHorizontally(Matrix * m1, Matrix * m2)
{
  if(m1->nGlobRows != m2->nGlobRows)
  {
    misc.error("Error: An internal error was happened. Matrices with different number of rows can not be joined horizontally.", 0);
  }
  joinMatrices(m1, subMatrix(0, 0, m1->nGlobRows, m1->nGlobCols), m2, subMatrix(0, m1->nGlobCols, m2->nGlobRows, m2->nGlobCols));
}

void Matrix::matrixToStandardVector(std::vector< std::vector<double> > & globalVector)
{
  if( this->distribution == diagonalDistribution )
  {
    globalVector.clear();
    if( communicator->mpiRoot == true )
    {
      globalVector.resize(this->nGlobRows, std::vector<double>(this->nGlobCols, 0));
      for(int i = 0; i<this->nGlobRows; i++)
      {
        globalVector[i][i] = this->m[i];
      }
    }
  }
  else
  {
    double * globalMatrix;
    globalVector.clear();
    
    if(communicator->mpiRoot)
    {
      globalMatrix = new double [this->nGlobRows*this->nGlobCols];
    }
    
    this->gatherMatrix(globalMatrix);
    
    if(communicator->mpiRoot)
    {
      globalVector.resize(this->nGlobRows, std::vector<double>(this->nGlobCols, 0));
      for(int c = 0; c<this->nGlobCols; c++)
      {
        for(int r = 0; r<this->nGlobRows; r++)
        {
          globalVector[r][c] = globalMatrix[c*this->nGlobRows + r];
        }
      }
    }
    
    if(communicator->mpiRoot)
    {
      delete [] globalMatrix;
    }
  }
}

void Matrix::matrixToStandardVector(std::vector<double> & globalVector)
{
  if( this->distribution == diagonalDistribution )
  {
    globalVector.clear();
    if( communicator->mpiRoot == true )
    {
      globalVector.resize(this->nGlobRows*this->nGlobCols, 0);
      for(int i = 0; i<this->nGlobRows; i++)
      {
        globalVector[i*this->nGlobRows + i] = this->m[i];
      }
    }
  }
  else
  {
    double * globalMatrix;
    globalVector.clear();
    
    if(communicator->mpiRoot)
    {
      globalMatrix = new double [this->nGlobRows*this->nGlobCols];
    }
    
    this->gatherMatrix(globalMatrix);
    
    if(communicator->mpiRoot)
    {
      globalVector.resize(this->nGlobRows*this->nGlobCols, 0);
      for(int c = 0; c<this->nGlobCols; c++)
      {
        for(int r = 0; r<this->nGlobRows; r++)
        {
          globalVector[c*this->nGlobRows + r] = globalMatrix[c*this->nGlobRows + r];
        }
      }
    }
    
    if(communicator->mpiRoot)
    {
      delete [] globalMatrix;
    }
  }
}

void Matrix::standardizeMatrix(RowColumn rowcolumn, int ddof)
{
  Matrix * colOnes;
  char matrixTransposeStatus;
  if(rowcolumn == column)
  {
    colOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->nGlobRows, 1);
    matrixTransposeStatus = 'T';
  }
  else
  {
    colOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->nGlobCols, 1);
    matrixTransposeStatus = 'N';
  }
  colOnes->fillWithConstant(1.);
  
  Matrix * mMeans = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  mMeans->multiply(this, matrixTransposeStatus, colOnes, 'N', 1./double(colOnes->nGlobRows));
  double * means = new double [mMeans->nGlobRows];
  mMeans->gatherMatrix(means);
  scatterVector(means, rowcolumn);
  delete mMeans;
  delete [] means;
  
  if(rowcolumn == column)
  {
    #pragma omp parallel for
    for(int c = 0; c<this->nCols; c++)
    {
      for(int r = 0; r<this->nRows; r++)
      {
        this->m[c*this->nRows + r] -= this->v[c];
      }
    }
  }
  else
  {
    #pragma omp parallel for
    for(int c = 0; c<this->nCols; c++)
    {
      for(int r = 0; r<this->nRows; r++)
      {
        this->m[c*this->nRows + r] -= this->v[r];
      }
    }
  }
  delete [] this->v;
  this->v = NULL;

  Matrix * matrixSquared = new Matrix(this);
  matrixSquared->elementWiseMultiplication(this);
  Matrix * mVars = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  mVars->multiply(matrixSquared, matrixTransposeStatus, colOnes, 'N', 1./(double(colOnes->nGlobRows) - double(ddof)));
  double * stds = new double [mVars->nGlobRows];
  mVars->gatherMatrix(stds);
  for(int i =0; i < mVars->nGlobRows; i++)
  {
    stds[i] = sqrt(stds[i]);
  }
  scatterVector(stds, rowcolumn);
  delete matrixSquared;
  delete mVars;
  delete [] stds;
  
  if(rowcolumn == column)
  {
    #pragma omp parallel for
    for(int c = 0; c<this->nCols; c++)
    {
      for(int r = 0; r<this->nRows; r++)
      {
        this->m[c*this->nRows + r] /= this->v[c];
      }
    }
  }
  else
  {
    #pragma omp parallel for
    for(int c = 0; c<this->nCols; c++)
    {
      for(int r = 0; r<this->nRows; r++)
      {
        this->m[c*this->nRows + r] /= this->v[r];
      }
    }
  }
  delete [] this->v;  
  this->v = NULL;
  
  delete colOnes;
}

void Matrix::centerMatrixRowsColumns(RowColumn rowcolumn)
{
  Matrix * colOnes;
  char matrixTransposeStatus;
  if(rowcolumn == column)
  {
    colOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->nGlobRows, 1);
    matrixTransposeStatus = 'T';
  }
  else
  {
    colOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->nGlobCols, 1);
    matrixTransposeStatus = 'N';
  }
  colOnes->fillWithConstant(1.);
  
  Matrix * mMeans = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  mMeans->multiply(this, matrixTransposeStatus, colOnes, 'N', 1./double(colOnes->nGlobRows));
  double * means = new double [mMeans->nGlobRows];
  mMeans->gatherMatrix(means);
  scatterVector(means, rowcolumn);
  delete mMeans;
  delete [] means;
  
  if(rowcolumn == column)
  {
    #pragma omp parallel for
    for(int c = 0; c<this->nCols; c++)
    {
      for(int r = 0; r<this->nRows; r++)
      {
        this->m[c*this->nRows + r] -= this->v[c];
      }
    }
  }
  else
  {
    #pragma omp parallel for
    for(int c = 0; c<this->nCols; c++)
    {
      for(int r = 0; r<this->nRows; r++)
      {
        this->m[c*this->nRows + r] -= this->v[r];
      }
    }
  }
  delete [] this->v;
  this->v = NULL;

  delete colOnes;
}


void Matrix::copySubMatrix(Matrix * m1, subMatrix smt, subMatrix sm1)
{
  int ia, ja, ib, jb, m, n;
  char t1 = 'N';
  
  misc.error("Error: An internal error was happened. This function cannot be used. See header file.", 0);
  
  if(smt.active == false)
  {
    if(this->nGlobRows != m1->nGlobRows || this->nGlobCols != m1->nGlobCols )
    {
      misc.error("Error: An internal error was happened", 0);
    }
  }
  if(sm1.active == false)
  {
    sm1 = subMatrix(0, 0, m1->nGlobRows, m1->nGlobCols);
  }
  
  if(smt.ir+smt.nr>this->nGlobRows || smt.ic+smt.nc>this->nGlobCols)
  {
    misc.error("Error: An internal error was happened", 0);
  }
  if(sm1.ir+sm1.nr>m1->nGlobRows || sm1.ic+sm1.nc>m1->nGlobCols)
  {
    misc.error("Error: An internal error was happened", 0);
  }
  if(smt.nr != sm1.nr || smt.nc != sm1.nc)
  {
    misc.error("Error: An internal error was happened", 0);
  }
  if(this->uplo != 'B')
  {
    this->symmetrizeTriangularMatrix();
  }
  if(m1->uplo != 'B')
  {
    m1->symmetrizeTriangularMatrix();
  }
  if(this->uplo != 'B' || m1->uplo != 'B')
  {
    misc.error("Error: An internal error was happened", 0);
  }
  this->symmetric = false;
  this->uplo = 'B';
  
  ia = sm1.ir + 1;
  ja = sm1.ic + 1;
  ib = smt.ir + 1;
  jb = smt.ic + 1;
  m = sm1.nr;
  n = sm1.nc;

  pdlacpy_(&(this->uplo), &m, &n, m1->m, &ia, &ja, m1->descriptor, this->m, &ib, &jb, this->descriptor);
  communicator->barrier();
}

void Matrix::showPartial(RowColumn rowcolumn, bool showv, bool showContext)
{
  for (int id = 0; id < communicator->nProc; ++id) {
    if (id == communicator->myId) {
      if(showContext == false)
      {
	std::cout << "Matrix " << communicator->myRow << " " << communicator->myCol  << " (" << this->nRows << ", " << this->nCols << ")" << std::endl;
        std::cout << "Flags: symm: " << this->symmetric << " pos: " << this->uplo  << " vec: " << this->vector << std::endl;
	
        std::cout << "     ";
        for (int c = 0; c < this->nCols; ++c) {
          std::cout << std::setw(4) << local2global(column, c) + 1;
        }
        std::cout << std::endl;
        
	if(this->v != NULL && rowcolumn == column && showv == true)
	{
	  for (int c = 0; c < this->nCols; c++)
	  {
	    std::cout << std::setw(4) << *(this->v + c);
	  }
	  std::cout << std::endl;
	  
	  for (int c = 0; c < this->nCols; c++)
	  {
	    std::cout << "----";
	  }
	  std::cout << std::endl;
	}
	
	for (int r = 0; r < this->nRows; r++)
	{
          std::cout << std::setw(4) << local2global(row, r) + 1 << ": ";
	  if(this->v != NULL && rowcolumn == row  && showv == true)
	  {
	    std::cout << std::setw(4) << *(this->v + r) << " | ";
	  }
	  for (int c = 0; c < this->nCols; ++c) {
	    std::cout << std::setw(4) << *(this->m + this->nRows*c + r);
	  }
	  std::cout << std::endl;
	}
	std::cout << std::endl;
      }
      else
      {
	bool local;
	double value;
	std::cout << "Blocks in process (" << communicator->myRow << ", " << communicator->myCol << "): " << std::endl;
	for (int r = 0; r < this->nGlobRows; r++)
	{
	  for (int c = 0; c < this->nGlobCols; ++c) {
	    value = matrixElement(r, c, 0., &local);
	    if(local)
	    {
	      std::cout << std::setw(3) << value << " ";
	    }
	    else
	    {
	      std::cout << std::setw(3) << "n" << " ";
	    }
	  }
	  std::cout << std::endl;
	}
	std::cout << std::endl;
      }
    } //end id == communicator->myId
    std::cout.flush();
    Cblacs_barrier(communicator->context, "All");
  } //end for loop
  
}

void Matrix::showGlobal(std::string name, bool symmetrize, int precision, double zeroThreshold)
{
  if(this->distribution == cyclicDistribution)
  {
    double *g;
  
    if (communicator->mpiRoot) {
      g = new double [this->nGlobRows*this->nGlobCols];
    }
    
    gatherMatrix(g);
    
    if (communicator->mpiRoot)
    {
      std::cout << "Global Matrix (" << this->nGlobRows << ", " << this->nGlobCols << "): " << name << std::endl;
      std::cout << "Cyclic distribution. Blocks (" << this->nBlockRows << ", " << this->nBlockCols << ")" << std::endl;
      std::cout << "Symmetric? " << (this->symmetric == true?"yes":"no") << std::endl;
      if(this->symmetric)
      {
        std::cout << "Symmetric matrix: " << this->uplo << std::endl;
        if(symmetrize)
        {
          for (int r = 0; r < this->nGlobRows; r++)
          {
            for (int c = 0; c < this->nGlobCols; c++)
            {
              if(this->uplo == 'L')
              {
                g[c*this->nGlobRows + r] = g[r*this->nGlobRows + c];
              }
              else
              {
                g[r*this->nGlobRows + c] = g[c*this->nGlobRows + r];
              }
            }
          }
        }
      }
    
      std::cout << "[" << std::endl;
      for (int r = 0; r < this->nGlobRows; r++)
      {
        std::cout << "[" << std::setprecision(precision) << std::setw(precision + 4) << ( (fabs(g[r])<zeroThreshold)?0:g[r] );
        for (int c = 1; c < this->nGlobCols; c++)
        {
          std::cout  << ", " << std::setprecision(precision) << std::setw(precision + 4) << ((fabs(g[c*this->nGlobRows + r])<zeroThreshold)?0:g[c*this->nGlobRows + r]);
        }
        std::cout << (r!=(this->nGlobRows-1)?"],":"]") << std::endl;
      }
      std::cout << "]" << std::endl;
      std::cout.flush();
    }
    
    if (communicator->mpiRoot) {
      delete [] g;
    }
  }
  else if(this->distribution == diagonalDistribution)
  {
    if(communicator->mpiRoot)
    {
      std::cout << "Global Matrix (" << this->nGlobRows << ", " << this->nGlobCols << "): " << name << std::endl;
      std::cout << "Diagonal distribution." << std::endl << std::endl;
      
      for (int i = 0; i<this->nGlobRows; i++)
      {
        std::cout << std::setprecision(precision) << this->m[i] << " ";
      }
      std::cout << std::endl << std::endl;
    }
  }
}

void Matrix::debugRead(std::string fname, int ttnr, int ttnc, int sourceProcessRow, int sourceProcessCol, int ttnblockr, int ttnblockc)
{
  misc.checkFileExists(fname);
  if(this->distribution != diagonalDistribution)
  {
    initParameters(ttnr, ttnc, ttnblockr, ttnblockc);
    fillWithConstant(0.);
    double *mGlob = NULL;
    if (communicator->myRow == sourceProcessRow && communicator->myCol == sourceProcessCol) {
      std::cout << "Scatering matrix from " << communicator->myRow << " " << communicator->myCol << std::endl;
      mGlob = new double[this->nGlobRows*this->nGlobCols];
      std::ifstream file(fname.c_str());
      for (int r = 0; r < this->nGlobRows; ++r) {
        for (int c = 0; c < this->nGlobCols; ++c) {
          file >> *(mGlob + this->nGlobRows*c + r);
        }
      }
    }
    //std::cout << "Pointer: " << mGlob << " in " << communicator->myRow << " " << communicator->myCol << std::endl;
    
    scatterMatrix(mGlob, sourceProcessRow, sourceProcessCol);
    if (communicator->myRow == sourceProcessRow && communicator->myCol == sourceProcessCol) {
      delete [] mGlob;
    }
  }
  else
  {
    initParameters(ttnr, ttnc, ttnblockr, ttnblockc, diagonalDistribution);
    if(communicator->mpiRoot == true)
    {
      std::ifstream file(fname.c_str());
      for (int r = 0; r < this->nGlobRows; ++r) {
        file >> this->m[r];
      }
      file.close();
    }
  }
  communicator->barrier();
  
}

void Matrix::debugWrite(std::string fname)
{
  double *g;
  
  if (communicator->mpiRoot) {
    g = new double [this->nGlobRows*this->nGlobCols];
  }
  
  if(this->symmetric)
  {
    symmetrizeTriangularMatrix();
  }
  
  gatherMatrix(g);
  
  if (communicator->mpiRoot)
  {
    Message message(fname);
    for (int r = 0; r < this->nGlobRows; ++r) {
      for (int c = 0; c < this->nGlobCols; ++c) {
        message << *(g + this->nGlobRows*c + r) << " ";
      }
      message << std::endl;
    }
  }
  
  if (communicator->mpiRoot) {
    delete [] g;
  }
}

void Matrix::fillWithRandom(double min, double max, long * seed)
{
  //srand(time(NULL));

  for(int c = 0; c<this->nCols; c++)
  {
    for(int r = 0; r<this->nRows; r++)
    {
      if(seed == NULL)
      {
        this->m[c*this->nRows + r] = min + ( (max-min)*double(rand() % 10000)/10000. );
      }
      else
      {
        this->m[c*this->nRows + r] = min + ( (max-min)*double(unif_rand_dbl(seed)) );
      }
    }
  }
}

