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

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include "labeledmatrix.h"
#include "misc.h"
#include "communicator.h"
#include "options.h"
#include "matrix.h"
#include "auxiliar.h"

LabeledMatrix::LabeledMatrix()
{
  this->matrix = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  
  this->rowLabels.clear();
  this->colLabels.clear();
}

LabeledMatrix::LabeledMatrix(int nRows, int nCols)
{
  if(nRows <= 1 || nCols <= 1)
  {
    misc.error("Error: An internal error was happened. The effects matrix cannot have dimensions smaller than 1.", 0);
  }
  
  this->matrix = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nRows, nCols);
  this->matrix->fillWithConstant(0.);
  
  this->rowLabels.clear();
  this->colLabels.clear();
}

LabeledMatrix::LabeledMatrix(LabeledMatrix * lm)
{
  this->rowLabels = lm->rowLabels;
  this->colLabels = lm->colLabels;
  
  this->matrix = new Matrix(lm->matrix);
}

LabeledMatrix::LabeledMatrix(std::string fn, std::vector<std::string> lColumns)
{
  this->matrix = NULL;
  
  this->rowLabels.clear();
  this->colLabels.clear();
  
  if( misc.gt(lColumns.size() == 0) )
  {
    load(fn);
  }
  else
  {
    loadRaw(fn, lColumns);
  }
}

LabeledMatrix::~LabeledMatrix()
{
  if( this->matrix != NULL )
  {
    delete this->matrix;
    this->matrix = NULL;
  }
}

void LabeledMatrix::clear()
{
  if(this->matrix != NULL)
  {
    delete this->matrix;
    this->matrix = NULL;
  }
  this->rowLabels.clear();
  this->colLabels.clear();
}

void LabeledMatrix::loadRaw(std::string fn, std::vector<std::string> lColumns)
{
  clear();
  
  misc.setGetElapsedTime("LabeledMatrixLoad");
  misc.message << "Reading data from file [ " << fn << " ] ..." << std::endl;
  
  double *valuesGlob = NULL;
  int nRows = 0;
  int nCols = 0;
  if(communicator->mpiRoot == true)
  {
    std::vector<std::string> header = getHeader(fn);
    std::set<std::string> stemp(header.begin(), header.end());
    if(header.size() != stemp.size())
    {
      misc.error("Error: An error has happened. The header of the file [ " + fn + " ] has repeated elements.", 0);
    }
    if(lColumns.size() >= header.size() )
    {
      misc.error("Error: An error has happened. The file [ " + fn + " ] has to have more than " + i2s(lColumns.size()) + " columns.", 0);
    }
    for(int i = 0; i<lColumns.size(); i++)
    {
      if(header[i] != lColumns[i])
      {
        misc.error("Error: An error has happened. The header element in column " + i2s(i + 1) + " ( '" + header[i] + "' ) of the file [ " + fn + " ] has to be '" + lColumns[i] + "'.", 0);
      }
    }
    
    for(int i = lColumns.size(); i<header.size(); i++)
    {
      this->colLabels.push_back(header[i]);
    }
    
    misc.checkFileExists(fn);
    std::ifstream file;
    file.open(fn.c_str());
    
    std::string line;
    getline(file,line); //Remove header we already have read.
    std::vector< std::vector<double> > values;
    while(getline(file,line))
    {
      std::istringstream sstemp(line);
      
      std::string rowLabel = "";
      for(int i = 0; i < lColumns.size(); i++)
      {
        std::string element;
        if((sstemp >> element).fail())
        {
          misc.error(" Error: Line " + i2s(values.size() + 2) + " in file [ " + fn + " ] does not have " + i2s(this->colLabels.size()) + " data columns.", 0);
        }
        rowLabel += ((i==0)?"":"@") + element;
      }
      this->rowLabels.push_back(rowLabel);
      
      std::vector<double> row;
      for(int i = 0; i < this->colLabels.size(); i++)
      {
        double value;
        if( (sstemp >> value).fail() )
        {
          misc.error("Error: The row " + i2s(values.size() + 2) + " in file [ " + fn + " ] contains an element which is an invalid number, or does not contain enough elements in the row.", 0);
        }
        row.push_back(value);
      }
      if(row.size() != this->colLabels.size())
      {
        misc.error(" Error: An internal error has happened when reading the file file [ " + fn + " ]. Not enough row elements.", 0);
      }
      std::string stemp;
      if( (sstemp >> stemp).fail() == false )
      {
        misc.error("Error: The row " + i2s(values.size() + 2) + " in file [ " + fn + " ] contains more elements than expected.", 0);
      }
      values.push_back(row);
    }
    file.close();
    
    if(values.size() == 0)
    {
      misc.error("Error: File [ " + fn + " ] is empty.", 0);
    }
    if(this->rowLabels.size() != values.size())
    {
      misc.error("Error: An internal error has happened when reading the file [ " + fn + " ]. Unexpected number of rows.", 0);
    }
    
    
    nRows = this->rowLabels.size();
    nCols = this->colLabels.size();
    valuesGlob = new double[nRows*nCols];
    for (int r = 0; r < nRows; ++r)
    {
      for (int c = 0; c < nCols; ++c)
      {
        valuesGlob[nRows*c + r] = values[r][c];
      }
    }
  }
  communicator->broadcast(&nRows);
  communicator->broadcast(&nCols);
  
  this->matrix = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nRows, nCols);
  this->matrix->scatterMatrix(valuesGlob);
  
  if(communicator->mpiRoot == true)
  {
    delete [] valuesGlob;
  }
  
  if(communicator->mpiRoot == true)
  {
    std::vector<std::string> temp = this->rowLabels;
    std::set<std::string> stemp(temp.begin(), temp.end());
    if(temp.size() != stemp.size())
    {
      misc.error("Error: An error has happened. The row labels of the file [ " + fn + " ] have repeated elements.", 0);
    }
  }
  
  misc.message << "Data loaded after " << misc.setGetElapsedTime("LabeledMatrixLoad", true) << std::endl;
}

std::vector<std::string> LabeledMatrix::getRowLabels()
{
  return this->rowLabels;
}

std::vector<std::string> LabeledMatrix::getColLabels()
{
  return this->colLabels;
}

Matrix * LabeledMatrix::getMatrix()
{
  return this->matrix;
}

void LabeledMatrix::setRowLabels(std::vector<std::string> labels)
{
  if( communicator->mpiRoot == true )
  {
    if( this->matrix->nGlobRows != labels.size() )
    {
      misc.error("Error: An internal error was happened. The labels for the effects matrix does not have the proper dimensions.", 0);
    }
    
    this->rowLabels = labels;
  }
}

void LabeledMatrix::setColLabels(std::vector<std::string> labels)
{
  if( communicator->mpiRoot == true )
  {
    if( this->matrix->nGlobCols != labels.size() )
    {
      misc.error("Error: An internal error was happened. The labels for the effects matrix does not have the proper dimensions.", 0);
    }
    
    this->colLabels = labels;
  }
}

void LabeledMatrix::insertRow(Matrix * row, int idx)
{
  if(row->nGlobCols != this->matrix->nGlobCols || row->nGlobRows != 1)
  {
    misc.error("Error: An internal error was happened. The row for inserting in effects matrix does not have the proper dimensions.", 0);
  }
  if(this->matrix->nGlobRows <= idx)
  {
    misc.error("Error: An internal error was happened. The row index for inserting in effects matrix is out of range.", 0);
  }
  
  this->matrix->add(row, 0., 1., subMatrix(idx, 0, 1, this->matrix->nGlobCols), subMatrix(0, 0, 1, row->nGlobCols));
}

void LabeledMatrix::insertCol(Matrix * col, int idx)
{
  if(col->nGlobRows != this->matrix->nGlobRows || col->nGlobCols != 1)
  {
    misc.error("Error: An internal error was happened. The col for inserting in effects matrix does not have the proper dimensions.", 0);
  }
  if(this->matrix->nGlobCols <= idx)
  {
    misc.error("Error: An internal error was happened. The col index for inserting in effects matrix is out of range.", 0);
  }
  
  this->matrix->add(col, 0., 1., subMatrix(0, idx, this->matrix->nGlobRows, 1), subMatrix(0, 0, col->nGlobRows, 1));
}

void LabeledMatrix::filterRowsAndCols(std::vector<std::string> & rowLabelsKeep, std::vector<std::string> & colLabelsKeep)
{
  if( communicator->mpiRoot == true && (this->matrix->nGlobCols != this->colLabels.size() || this->matrix->nGlobRows != this->rowLabels.size()) )
  {
    misc.error("Error: An internal error was happened when filtering effects matrix. The labels for the effects matrix does not have the proper dimensions.", 0);
  }
  
  
  std::vector<std::string> newRowLabels;
  int * keepRowIdxs;
  int nRowIdxs = 0;
  
  std::vector<std::string> newColLabels;
  int * keepColIdxs;
  int nColIdxs = 0;
  
  if( communicator->mpiRoot == true )
  {
    std::map<std::string, int> rowLabelIdxs;
    std::map<std::string, int> colLabelIdxs;
    
    for(int i = 0; i < this->rowLabels.size(); i++)
    {
      rowLabelIdxs[ this->rowLabels[i] ] = i;
    }
    for(int i = 0; i < this->colLabels.size(); i++)
    {
      colLabelIdxs[ this->colLabels[i] ] = i;
    }
    
    keepRowIdxs = new int [rowLabelsKeep.size()];
    for(int i = 0; i<rowLabelsKeep.size(); i++)
    {
      std::string label = rowLabelsKeep[i];
      if(rowLabelIdxs.count(label) == 0)
      {
        misc.error("Error: An internal error was happened when filtering effects matrix. Inexistent row label.", 0);
      }
      newRowLabels.push_back(label);
      keepRowIdxs[ i ] = rowLabelIdxs[label];
      nRowIdxs++;
    }
    this->rowLabels = newRowLabels;
    
    keepColIdxs = new int [colLabelsKeep.size()];
    for(int i = 0; i<colLabelsKeep.size(); i++)
    {
      std::string label = colLabelsKeep[i];
      if(colLabelIdxs.count(label) == 0)
      {
        misc.error("Error: An internal error was happened when filtering effects matrix. Inexistent col label.", 0);
      }
      newColLabels.push_back(label);
      keepColIdxs[ i ] = colLabelIdxs[label];
      nColIdxs++;
    }
    this->colLabels = newColLabels;
  }
  
  communicator->broadcast(&nRowIdxs);
  communicator->broadcast(&nColIdxs);
  
  if( communicator->mpiRoot == false )
  {
    keepRowIdxs = new int [nRowIdxs];
    keepColIdxs = new int [nColIdxs];
  }
  
  communicator->broadcast(keepRowIdxs, nRowIdxs);
  communicator->broadcast(keepColIdxs, nColIdxs);
  
  std::map<int, int> rowsOriginDestination;
  for(int i = 0; i < nRowIdxs; i++)
  {
    rowsOriginDestination[ keepRowIdxs[i] ] = i;
  }
  
  std::map<int, int> colsOriginDestination;
  for(int i = 0; i < nColIdxs; i++)
  {
    colsOriginDestination[ keepColIdxs[i] ] = i;
  }
  
  Matrix * newMatrix = new Matrix(cyclicDistribution);
  this->matrix->generalResorting(newMatrix, rowsOriginDestination, colsOriginDestination, true);
  newMatrix->symmetric = false;
  newMatrix->uplo = 'B';
  delete this->matrix;
  this->matrix = newMatrix;
  
  delete [] keepRowIdxs;
  delete [] keepColIdxs;
}

void LabeledMatrix::appendVertically(LabeledMatrix * lm)
{
  Matrix * temp = new Matrix();
  
  if(this->colLabels != lm->getColLabels())
  {
    misc.error("Error: When joining two labeled matrices, the column labels of the matrices do not match.", 0);
  }
  
  temp->joinMatricesVertically(this->matrix, lm->getMatrix());
  delete this->matrix;
  this->matrix = temp;
  
  std::vector<std::string> labeltemp = lm->getRowLabels();
  this->rowLabels.insert(this->rowLabels.end(), labeltemp.begin(), labeltemp.end());
  
  if(communicator->mpiRoot == true && this->rowLabels.size() != this->matrix->nGlobRows)
  {
    misc.error("Error: An internal error has happened when joining two labeled matrices, the row labels size differs from the number of matrix rows.", 0);
  }
}

void LabeledMatrix::appendHorizontally(LabeledMatrix * lm)
{
  Matrix * temp = new Matrix();
  
  if(this->rowLabels != lm->getRowLabels())
  {
    misc.error("Error: When joining two labeled matrices, the row labels of the matrices do not match.", 0);
  }
  
  temp->joinMatricesHorizontally(this->matrix, lm->getMatrix());
  delete this->matrix;
  this->matrix = temp;
  
  std::vector<std::string> labeltemp = lm->getColLabels();
  this->colLabels.insert(this->colLabels.end(), labeltemp.begin(), labeltemp.end());
  
  if(communicator->mpiRoot == true && this->colLabels.size() != this->matrix->nGlobCols)
  {
    misc.error("Error: An internal error has happened when joining two labeled matrices, the column labels size differs from the number of matrix columns.", 0);
  }
}

void LabeledMatrix::save(std::string fn)
{
  if( communicator->mpiRoot == true && (this->matrix->nGlobCols != this->colLabels.size() || this->matrix->nGlobRows != this->rowLabels.size()) )
  {
    misc.error("Error: An internal error was happened when saving effects matrix. The labels for the effects matrix does not have the proper dimensions.", 0);
  }
  
  misc.setGetElapsedTime("EffectsWrite");
  misc.message << "Writing effects to [ " << fn << ".(rowids/colids/dat) ]..." << std::endl;
  
  if(communicator->mpiRoot)
  {
    Message messageRows(fn + ".rowids");
    for(int i = 0; i < rowLabels.size(); i++)
    {
      messageRows << rowLabels[i] << std::endl;
    }
    
    Message messageCols(fn + ".colids");
    for(int i = 0; i < colLabels.size(); i++)
    {
      messageCols << colLabels[i] << std::endl;
    }
  }
  
  //Define header
  unsigned char *header = new unsigned char [14];
  header[0] = 'E';
  header[1] = 'F';
  header[2] = 'F';
  header[3] = 'E';
  header[4] = 'C';
  header[5] = 'T';
  header[6] = 'S';
  header[7] = 0x5A; //Magic number 1
  header[8] = 0x99; //Magic number 2
  header[9] = 0x1;  //Version
  header[10] = 0x1; //1 if stored doubles
  header[11] = sizeof(double); //double size
  header[12] = 0x0; //unused
  header[13] = 0x0; //unused
  
  if(options.useMPIForWriting == true)
  {
    this->matrix->writeMatrixMPI(fn + ".dat", (char*)header, 14);
  }
  else
  {
    misc.error("Error: Effects currently cannot be stored using this options.", 0);
  }
  
  delete [] header;
  
  misc.message << "Effects file stored after " << misc.setGetElapsedTime("EffectsWrite", true) << std::endl;
}


void LabeledMatrix::load(std::string fn)
{
  misc.setGetElapsedTime("EffectsLoad");
  misc.message << "Loading effects file from [ " << fn << ".(rowids/colids/dat) ]..." << std::endl;
  
  this->rowLabels.clear();
  this->colLabels.clear();
  
  if(this->matrix != NULL)
  {
    delete this->matrix;
    this->matrix = NULL;
  }

  std::string fdataname = fn + ".dat";
  if(communicator->mpiRoot)
  {
    getListFromFile(fn + ".rowids", this->rowLabels);
    getListFromFile(fn + ".colids", this->colLabels);
    
    //Read data
    std::ifstream file;
    misc.checkFileExists(fdataname);
    file.open(fdataname.c_str(), std::ifstream::in | std::ifstream::binary);
    
    //Read header
    unsigned char *header = new unsigned char [14];
    file.read((char*)header, 14);
    if (
      header[0] != 'E' ||
      header[1] != 'F' ||
      header[2] != 'F' ||
      header[3] != 'E' ||
      header[4] != 'C' ||
      header[5] != 'T' ||
      header[6] != 'S' ||
      header[7] != 0x5A ||
      header[8] != 0x99 ||
      header[9] != 0x1 ||
      header[10] != 0x1
    )
    {
      misc.error("Error: The file [ " + fdataname + " ] do not have a proper header or is not a valid file.", 0);
    }
    if( header[11] != sizeof(double) )
    {
      misc.error("Error: Data types in file [ " + fdataname + " ] differ with this system datatypes. Are you loading a grm computed in a different system?", 0);
    }
    delete [] header;
    file.close();
  }
  
  int nRows = this->rowLabels.size();
  int nCols = this->colLabels.size();
  
  communicator->broadcast(&nRows);
  communicator->broadcast(&nCols);
  
  communicator->barrier();
  
  this->matrix = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nRows, nCols, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->matrix->readMatrixMPI(fdataname, 14);
  
  misc.message << "Effects file loaded after " << misc.setGetElapsedTime("EffectsLoad", true) << std::endl;
}

void LabeledMatrix::show(int padding)
{
  if( communicator->mpiRoot == true && (this->matrix->nGlobCols != this->colLabels.size() || this->matrix->nGlobRows != this->rowLabels.size()) )
  {
    misc.error("Error: An internal error was happened when showing effects matrix. The labels for the effects matrix does not have the proper dimensions.", 0);
  }
  
  double *g;
  
  if (communicator->mpiRoot) {
    g = new double [this->matrix->nGlobRows*this->matrix->nGlobCols];
  }

  this->matrix->gatherMatrix(g);
  if (communicator->mpiRoot) {
    misc.message << "Effects:" << std::endl;
    misc.message << this->rowLabels.size() << " " << this->colLabels.size() << std::endl;
    misc.message << "     ";
    for (int i = 0; i < this->colLabels.size(); ++i)
    {
      misc.message << std::setw(padding) << this->colLabels[i] << " ";
    }
    misc.message << std::endl;
    for (int r = 0; r < this->matrix->nGlobRows; ++r)
    {
      misc.message << this->rowLabels[r] << ": " << std::setw(10);
      for (int c = 0; c < this->matrix->nGlobCols; ++c)
      {
        misc.message << std::setw(padding) << *(g + this->matrix->nGlobRows*c + r) << " ";
      }
      misc.message << "\n";
    }
    misc.message << std::endl;
  }
      
  if (communicator->mpiRoot) {
    delete [] g;
  }
}