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

#include "blockmatrix.h"
#include "matrix.h"
#include "misc.h"
#include "communicator.h"
#include "auxiliar.h"

#include <vector>


BlockMatrix::BlockMatrix()
{
  this->m.clear();
  
  this->nGlobRows = 0;
  this->nGlobCols = 0;
  
  this->nBlockRows = 0;
  this->nBlockCols = 0;
}

BlockMatrix::BlockMatrix(std::vector< std::vector<Matrix*> > & mi)
{
  this->m.clear();
  
  this->nGlobRows = 0;
  this->nGlobCols = 0;
  
  this->nBlockRows = mi.size();
  if( this->nBlockRows < 1)
  {
    misc.error("Error: An internal error was happened. A block matrix cannot be created from an empty array.", 0);
  }
  
  this->m = std::vector< std::vector<Matrix*> >(this->nBlockRows, std::vector<Matrix*>());
  
  this->nBlockCols = 0;
  for(int r = 0; r < this->nBlockRows; r++)
  {
    if( mi[r].size() == 0 )
    {
      misc.error("Error: An internal error was happened. A block matrix cannot be created from an empty array.", 0);
    }
    if( (this->nBlockCols != 0) && (mi[r].size() != this->nBlockCols) )
    {
      misc.error("Error: An internal error was happened. The block matrix is not properly formated.", 0);
    }
    this->nBlockCols = mi[r].size();
    
    for(int c = 0; c < this->nBlockCols; c++)
    {
      if( r != 0 )
      {
        if( this->m[r-1][c]->nGlobCols != mi[r][c]->nGlobCols )
        {
          misc.error("Error: An internal error was happened. Matrices in a block matrix have discordant number of columns.", 0);
        }
      }
      else
      {
        this->nGlobCols += mi[0][c]->nGlobCols;
      }
      if( c != 0 )
      {
        if( mi[r][c - 1]->nGlobRows != mi[r][c]->nGlobRows )
        {
          misc.error("Error: An internal error was happened. Matrices in a block matrix have discordant number of rows.", 0);
        }
      }
        
      Matrix *mTemp = new Matrix(mi[r][c]);
      this->m[r].push_back(mTemp);
    }
    
    this->nGlobRows += mi[r][0]->nGlobRows;
  }
}

BlockMatrix::BlockMatrix(BlockMatrix & srcBlockMatrix)
{
  this->m.clear();
  
  this->nGlobRows = srcBlockMatrix.nGlobRows;
  this->nGlobCols = srcBlockMatrix.nGlobCols;
  
  this->nBlockRows = srcBlockMatrix.nBlockRows;
  this->nBlockCols = srcBlockMatrix.nBlockCols;
  
  this->m = std::vector< std::vector<Matrix*> >(this->nBlockRows, std::vector<Matrix*>());
  
  for(int r = 0; r < this->nBlockRows; r++)
  {
    for(int c = 0; c < this->nBlockCols; c++)
    {
      Matrix *mTemp = new Matrix(srcBlockMatrix.m[r][c]);
      this->m[r].push_back(mTemp);
    }
  }
}

BlockMatrix::~BlockMatrix()
{
  clear();
}

void BlockMatrix::clear()
{
  for(int r = 0; r < this->nBlockRows; r++)
  {
    for(int c = 0; c < this->nBlockCols; c++)
    {
      if( this->m[r][c] != NULL )
      {
        delete this->m[r][c];
      }
    }
  }
  
  this->m.clear();
  
  this->nGlobRows = 0;
  this->nGlobCols = 0;
  
  this->nBlockRows = 0;
  this->nBlockCols = 0;
}

void BlockMatrix::addBlockRow(std::vector<Matrix*> & blockRow)
{
  if( this->nBlockRows != 0 && this->nBlockCols != blockRow.size() )
  {
    misc.error("Error: An internal error was happened. A row of blocks cannot be added in a block matrix. Discordant number of block columns.", 0);
  }
  if( blockRow.size() < 1)
  {
    misc.error("Error: An internal error was happened. A row of 0 blocks cannot be added in a block matrix.", 0);
  }
  
  if( this->nBlockRows != 0 )
  {
    for(int c = 0; c < this->nBlockCols; c++)
    {
      if( this->m[0][c]->nGlobCols != blockRow[c]->nGlobCols )
      {
        misc.error("Error: An internal error was happened. Block row cannot be added. Matrices in row have discordant number of columns with current block matrix.", 0);
      }
      if( c != 0 )
      {
        if( blockRow[c - 1]->nGlobRows != blockRow[c]->nGlobRows )
        {
          misc.error("Error: An internal error was happened. Block row cannot be added. Matrices in row have discordant number of rows.", 0);
        }
      }
    }
  }
  else
  {
    if(this->nBlockCols != 0)
    {
      misc.error("Error: An internal error was happened when adding a row of blocks in a Block Matrix.", 0);
    }
  }
  
  int nParcialRows = 0;
  for(int c = 0; c < blockRow.size(); c++)
  {
    if( c != 0 && nParcialRows != blockRow[c]->nGlobRows )
    {
      misc.error("Error: An internal error was happened when adding a row of blocks in a Block Matrix. Inconsistent number of rows in each block", 0);
    }
    if(this->nBlockRows == 0)
    {
      this->nGlobCols += blockRow[c]->nGlobCols;
    }
    nParcialRows = blockRow[c]->nGlobRows;
  }
  this->nGlobRows += nParcialRows;
  
  this->m.push_back(blockRow);
  if( this->nBlockRows == 0)
  {
    this->nBlockCols = blockRow.size();
  }
  this->nBlockRows++;
}

void BlockMatrix::multiply(BlockMatrix & m1, BlockMatrix & m2, double factor)
{
  if( m1.nBlockCols != m2.nBlockRows )
  {
    misc.error("Error: An internal error was happened. Block matrices with a discrepant number of block rows and block columns cannot be multiplied.", 0);
  }
  if( m1.nGlobCols != m2.nGlobRows )
  {
    misc.error("Error: An internal error was happened. Block matrices with a discrepant number of rows and columns cannot be multiplied.", 0);
  }
  if( m1.nBlockRows == 0 || m1.nBlockCols == 0 || m2.nBlockRows == 0 || m2.nBlockCols == 0)
  {
    misc.error("Error: An internal error was happened. Empty block matrices cannot be multiplied.", 0);
  }
  if( m1.nGlobRows == 0 || m1.nGlobCols == 0 || m2.nGlobRows == 0 || m2.nGlobCols == 0)
  {
    misc.error("Error: An internal error was happened. Empty block matrices cannot be multiplied.", 0);
  }
  
  clear();
  
  for(int m1r = 0; m1r < m1.nBlockRows; m1r++)
  {
    std::vector<Matrix*> blockRow;
    for(int m2c = 0; m2c < m2.nBlockCols; m2c++)
    {
      Matrix * resultBlock = NULL;
      for(int m1c = 0; m1c < m1.nBlockCols; m1c++)
      {
        Matrix * temp = new Matrix();
        temp->multiply(m1.m[m1r][m1c], 'N', m2.m[m1c][m2c], 'N', factor);
        if(resultBlock == NULL)
        {
          resultBlock = temp;
        }
        else
        {
          resultBlock->add(temp);
          delete temp;
        }
      }
      blockRow.push_back(resultBlock);
    }
    addBlockRow(blockRow);
  }
}

void BlockMatrix::multiply(BlockMatrix & m1, Matrix * m2, double factor)
{
  if( m1.nGlobCols != m2->nGlobRows )
  {
    misc.error("Error: An internal error was happened. Block matrices with a discrepant number of rows and columns cannot be multiplied.", 0);
  }
  if( m1.nBlockRows == 0 || m1.nBlockCols == 0 )
  {
    misc.error("Error: An internal error was happened. Empty block matrices cannot be multiplied.", 0);
  }
  if( m1.nGlobRows == 0 || m1.nGlobCols == 0 || m2->nGlobRows == 0 || m2->nGlobCols == 0)
  {
    misc.error("Error: An internal error was happened. Empty block matrices cannot be multiplied.", 0);
  }
  
  clear();
  
  std::vector<Matrix*> m2Blocks;
  int rowShift = 0;
  for(int m1c = 0; m1c < m1.nBlockCols; m1c++)
  {
    Matrix * m2Block = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, m1.m[0][m1c]->nGlobCols, m2->nGlobCols);
    m2Block->fillWithConstant(0.);
    subMatrix smDest(m2Block);
    subMatrix smSrc(rowShift, 0, m1.m[0][m1c]->nGlobCols, m2->nGlobCols);
    m2Block->add(m2, 1., 1., smDest, smSrc);
    m2Blocks.push_back(m2Block);
    rowShift += m1.m[0][m1c]->nGlobCols;
  }
  
  for(int m1r = 0; m1r < m1.nBlockRows; m1r++)
  {
    std::vector<Matrix*> blockRow;
    Matrix * resultBlock = NULL;
    for(int m1c = 0; m1c < m1.nBlockCols; m1c++)
    {
      Matrix * temp = new Matrix();
      temp->multiply(m1.m[m1r][m1c], 'N', m2Blocks[m1c], 'N', factor);
      if(resultBlock == NULL)
      {
        resultBlock = temp;
      }
      else
      {
        resultBlock->add(temp);
        delete temp;
      }
    }
    blockRow.push_back(resultBlock);
    addBlockRow(blockRow);
  }
  
  for(int m1c = 0; m1c < m1.nBlockCols; m1c++)
  {
    delete m2Blocks[m1c];
  }
  m2Blocks.clear();
}

void BlockMatrix::add(BlockMatrix & m1, double thisFactor, double m1Factor)
{
  if( m1.nBlockRows != this->nBlockRows || m1.nBlockCols != this->nBlockCols )
  {
    misc.error("Error: An internal error was happened. Block matrices with a discrepant number of block rows and block columns. They cannot be added.", 0);
  }
  if( m1.nGlobRows != this->nGlobRows || m1.nGlobCols != this->nGlobCols )
  {
    misc.error("Error: An internal error was happened. Block matrices with a discrepant number of rows and columns. They cannot be added.", 0);
  }
  if( m1.nBlockRows == 0 || m1.nBlockCols == 0 || this->nBlockRows == 0 || this->nBlockCols == 0)
  {
    misc.error("Error: An internal error was happened. Empty block matrices cannot be added.", 0);
  }
  
  for(int r = 0; r < this->nBlockRows; r++)
  {
    for(int c = 0; c < this->nBlockCols; c++)
    {
      this->m[r][c]->add(m1.m[r][c], thisFactor, m1Factor);
    }
  }
}

bool BlockMatrix::invert(double * logDeterminant)
{
  if( this->nBlockRows != this->nBlockCols )
  {
    misc.error("Error: An internal error was happened. A block matrix with different number of blocks in the rows and in the columns cannot be inverted.", 0);
  }
  
  if( this->nBlockRows == 1 && this->nBlockCols == 1)
  {
    if( this->m[0][0]->symmetric == false || this->m[0][0]->distribution != diagonalDistribution )
    {
      misc.error("Error: An internal error was happened. A block matrix can only be inverted if all blocks are diagonal.", 0);
    }
    this->m[0][0]->symmetricInvert(logDeterminant);
  }
  else
  {
    //Divide this matrix in four blocks
    // [ A B ]
    // [ C D ]
    // Where D contains only the last block row and column, amd A, B, C are defined accordingly. Then use the block inversion when matrix is formed by four blocks.

    std::vector< std::vector<Matrix*> > ABlocks;
    std::vector< std::vector<Matrix*> > BBlocks;
    std::vector< std::vector<Matrix*> > CBlocks;
    std::vector< std::vector<Matrix*> > DBlocks;
    
    for(int r = 0; r < this->nBlockRows; r++)
    {
      std::vector<Matrix*> blockRowLeft;
      std::vector<Matrix*> blockRowRight;
      for(int c = 0; c < this->nBlockCols; c++)
      {
        if( c < (this->nBlockCols - 1) )
        {
          blockRowLeft.push_back(this->m[r][c]);
        }
        else
        {
          blockRowRight.push_back(this->m[r][c]);
        }
      }
      if( r < (this->nBlockRows - 1) )
      {
        ABlocks.push_back(blockRowLeft);
        BBlocks.push_back(blockRowRight);
      }
      else
      {
        CBlocks.push_back(blockRowLeft);
        DBlocks.push_back(blockRowRight);
      }
    }
    
    BlockMatrix A(ABlocks);
    BlockMatrix B(BBlocks);
    BlockMatrix C(CBlocks);
    BlockMatrix D(DBlocks);
    
    bool inverted;
    
    double DLogDeterminant;
    BlockMatrix D_1(D);
    inverted = D_1.invert(&DLogDeterminant);
    if( inverted == false )
    {
      return false;
    }
    BlockMatrix D_1C;
    D_1C.multiply(D_1, C);
    BlockMatrix A_BD_1C_1;
    A_BD_1C_1.multiply(B, D_1C);
    A_BD_1C_1.add(A, -1., 1.);
    double A_BD_1CLogDeterminant;
    inverted = A_BD_1C_1.invert(&A_BD_1CLogDeterminant);
    if( inverted == false )
    {
      return false;
    }
    BlockMatrix minD_1CA_BD_1C_1;
    minD_1CA_BD_1C_1.multiply(D_1C, A_BD_1C_1, -1.);
    D_1C.clear();
    
    BlockMatrix A_1(A);
    inverted = A_1.invert();    
    if( inverted == false )
    {
      return false;
    }
    BlockMatrix A_1B;
    A_1B.multiply(A_1, B);
    BlockMatrix D_CA_1B_1;
    D_CA_1B_1.multiply(C, A_1B);
    D_CA_1B_1.add(D, -1., 1.);
    inverted = D_CA_1B_1.invert();
    if( inverted == false )
    {
      return false;
    }
    BlockMatrix minA_1BD_CA_1B_1;
    minA_1BD_CA_1B_1.multiply(A_1B, D_CA_1B_1, -1.);
    A_1B.clear();
    
    clear();
    
    for(int r = 0; r < A.nBlockRows; r++)
    {
      std::vector<Matrix*> blockRow;
      for(int c = 0; c < A.nBlockCols; c++)
      {
        blockRow.push_back(A_BD_1C_1.m[r][c]);
        A_BD_1C_1.m[r][c] = NULL;
      }
      blockRow.push_back(minA_1BD_CA_1B_1.m[r][0]);
      minA_1BD_CA_1B_1.m[r][0] = NULL;
      addBlockRow(blockRow);
    }
    std::vector<Matrix*> blockRow;
    for(int c = 0; c < A.nBlockCols; c++)
    {
      blockRow.push_back(minD_1CA_BD_1C_1.m[0][c]);
      minD_1CA_BD_1C_1.m[0][c] = NULL;
    }
    blockRow.push_back(D_CA_1B_1.m[0][0]);
    D_CA_1B_1.m[0][0] = NULL;
    addBlockRow(blockRow);
    
    A_BD_1C_1.clear();
    minA_1BD_CA_1B_1.clear();
    D_CA_1B_1.clear();
    minD_1CA_BD_1C_1.clear();
    
    if(logDeterminant != NULL)
    {
      *logDeterminant = DLogDeterminant + A_BD_1CLogDeterminant;
      communicator->broadcast(logDeterminant, 1);
    }
  }
  
  return true;
}

Matrix* BlockMatrix::block2distributed()
{
  if(this->nBlockRows < 1 || this->nBlockCols < 1)
  {
    misc.error("Error: An internal error was happened. A cyclic matrix cannot be generated from an empty block matrix.", 0);
  }
  
  Matrix* result = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->nGlobRows, this->nGlobCols);
  result->fillWithConstant(0.);
  
  int rowShift = 0;
  for(int r = 0; r < this->nBlockRows; r++)
  {
    int colShift = 0;
    for(int c = 0; c < this->nBlockCols; c++)
    {
      result->add( this->m[r][c], 1., 1., subMatrix(rowShift, colShift ,this->m[r][c]->nGlobRows, this->m[r][c]->nGlobCols), subMatrix(this->m[r][c]) );
      colShift += this->m[r][c]->nGlobCols;
    }
    rowShift += this->m[r][0]->nGlobRows;
  }
  
  return result;
}

void BlockMatrix::showGlobal()
{
  misc.message << "Blocks: (" << this->nBlockRows << ", " << this->nBlockCols << ")" << std::endl;
  misc.message << "Global dimensions: (" << this->nGlobRows << ", " << this->nGlobCols << ")" << std::endl;
  
  Matrix * temp = block2distributed();
  temp->showGlobal();
  delete temp;
}