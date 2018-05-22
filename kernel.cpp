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

#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

#if defined(BOOSTLIB) && defined(ZLIB)
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#endif

#include "matrix.h"
#include "kernel.h"
#include "genotype.h"
#include "global.h"
#include "options.h"
#include "communicator.h"
#include "misc.h"
#include "auxiliar.h"

Kernel::Kernel()
{
  this->type = kernelGRM;
  
  this->name = "default";
  
  this->kernel = NULL;
  this->N = NULL;
  this->eigenValues = NULL;
  this->eigenVectors = NULL;
  
  this->diagonalized = false;
  this->normalized = false;
  this->asymmetric = false;
  
  this->genotypes = NULL;
  this->covariates = NULL;
}

Kernel::Kernel(Genotype * genotype, bool normalizeGRM)
{
  this->type = kernelGRM;
  
  this->name = "default";
  
  this->asymmetric = false;
  this->diagonalized = false;
  
  if(communicator->mpiRoot)
  {
    this->individuals.clear();
    for(int i = 0; i<genotype->nIndividuals; i++)
    {
      this->individuals.push_back(genotype->individuals[i]);
      std::string key = genotype->individuals[i].familyID + "@" + genotype->individuals[i].individualID;
      this->individualIds.push_back(key);
      if(this->individualIdsIdx.count(key) != 0)
      {
	misc.error("Error: The individual with family Id: " + genotype->individuals[i].familyID + " and individual Id: " + genotype->individuals[i].individualID + " appears more than one time in the genotypes.", 0);
      }
      this->individualIdsIdx[key] = i;
    }
    this->nIndividuals = this->individuals.size();
    this->randomVarNames = genotype->SNPIds;
  }
  communicator->broadcast(&this->nIndividuals, 1);
  
  misc.setGetElapsedTime("ComputeGRM");
  misc.message << "Starting GRM computation..." << std::endl;
  
  this->kernel = new Matrix(cyclicDistribution);
  genotype->normalizeGenotypes();
  this->kernel->multiply(genotype->genotypes, 'T', genotype->genotypes, 'N');
  
  if( options.flatGRMNormalization == false )
  {
    misc.message << "GRM will be normalized using the individuals shared non-missing SNPs." << std::endl;
    this->N = new Matrix(cyclicDistribution);
    this->N->multiply(genotype->missings, 'T', genotype->missings, 'N');
  }
  else
  {
    misc.message << "GRM will be normalized for all individuals by the total number of SNPs: " << double(genotype->nSNPs) << "." << std::endl;
    this->N = new Matrix(cyclicDistribution, this->kernel->nGlobRows, this->kernel->nGlobCols);
    this->N->fillWithConstant( double(genotype->nSNPs) );
    this->N->uplo = this->kernel->uplo;
    this->N->symmetric = true;
  }


  this->normalized = false;
  if( normalizeGRM  == true)
  {
    normalize();
  }
  
  this->eigenValues = NULL;
  this->eigenVectors = NULL;
  
  this->genotypes = NULL;
  this->covariates = NULL;
  
  misc.message << "GRM computation finished after " << misc.setGetElapsedTime("ComputeGRM", true) << "." << std::endl;
}

Kernel::Kernel(std::string f, KernelType readType, int useColumn, std::set<std::string> keepIndividuals)
{
  this->type = readType;
  
  this->name = "default";
  
  this->kernel = NULL;
  this->N = NULL;
  this->eigenValues = NULL;
  this->eigenVectors = NULL;
  
  this->asymmetric = false;
  this->diagonalized = false;
  this->normalized = true;
  
  this->genotypes = NULL;
  this->covariates = NULL;
  
  if(this->type == kernelGRM)
  {
    readKernel(f);
  }
  else if(this->type == kernelFromDiscreteCovariates)
  {
    createKernelFromDiscreteCovariates(f, useColumn, keepIndividuals);
  }
  else if(this->type == kernelFromMultiDiscreteCovariates)
  {
    createKernelFromMultipleDiscreteCovariates(f, useColumn, keepIndividuals);
  }
  else if(this->type == kernelFromContinuousCovariates)
  {
    createKernelFromContinuousCovariates(f);
  }
  else if(this->type == kernelSquaredExponential)
  {
    createKernelSquaredExponential(f, keepIndividuals);
  }
  else if(this->type == kernelGCTAGRM)
  {
    this->type = kernelGRM;
    readGCTAGRM(f);
  }
  else
  {
    misc.error("Error: An internal error was happened. Invalid Kernel type when creating a new kernel.", 0);
  }
}

Kernel::Kernel(Kernel * k1, Kernel * k2)
{
  if(k1->diagonalized == true)
  {
    misc.error("Error: An error has happened. Diagonal covariances ( " + k1->name + " ) cannot be used for computing interaction covariances.", 0);
  }
  if(k2->diagonalized == true)
  {
    misc.error("Error: An error has happened. Diagonal covariances ( " + k2->name + " ) cannot be used for computing interaction covariances.", 0);
  }

  this->kernel = new Matrix(k1->kernel);
  if(k1->N != NULL)
  {
    this->N = new Matrix(k1->N);
  }
  else
  {
    this->N = NULL;
  }
  this->eigenValues = NULL;
  this->eigenVectors = NULL;
  
  copyParameters(k1);
  
  this->genotypes = NULL;
  this->covariates = NULL;
  
  std::vector<std::string> tempIndividuals = orderVectorAsTemplate(this->individualIds, k2->individualIds);
  if( communicator->mpiRoot == true && tempIndividuals.size() == 0 )
  {
    misc.error("Error: An error has happened. There are not overlapping individuals when computing the interaction between " + this->name + " and " + k2->name + ".", 0);
  }
  
  filterIndividuals(tempIndividuals);
  
  Kernel * k2temp = NULL;
  if( misc.gt(k2->individualIds != this->individualIds) )
  {
    k2temp = new Kernel(k2);
    k2temp->filterIndividuals(tempIndividuals);
  }
  else
  {
    k2temp = k2;
  }

  if( this->individualIds != k2temp->individualIds )
  {
    misc.error("Error: An internal error has happened. The individuals when computing the interaction between " + this->name + " and " + k2temp->name + " covariances do not match.", 0);
  }
  
  this->kernel->symmetrizeTriangularMatrix();
  k2temp->kernel->symmetrizeTriangularMatrix();
  
  this->kernel->elementWiseMultiplication(k2temp->kernel);
  
  this->name = this->name + "x" + k2->name;
  this->type = kernelInteraction;
  
  if( k2temp != k2 )
  {
    delete k2temp;
  }
  if( this->N != NULL )
  {
    delete this->N;
    this->N = NULL;
  }
  
  this->randomVarNames.clear();
}

Kernel::Kernel(Kernel * kernel)
{
  if(kernel->diagonalized == false)
  {
    this->kernel = new Matrix(kernel->kernel);
    if(kernel->N != NULL)
    {
      this->N = new Matrix(kernel->N);
    }
    else
    {
      this->N = NULL;
    }
    this->eigenValues = NULL;
    this->eigenVectors = NULL;
  }
  else
  {
    this->kernel = NULL;
    this->N = NULL;
    this->eigenValues = new Matrix(kernel->eigenValues);
    this->eigenVectors = new Matrix(kernel->eigenVectors);
  }
  
  copyParameters(kernel);
  
  this->genotypes = NULL;
  this->covariates = NULL;
}

Kernel::Kernel(Kernel * kernel, KernelType createType, Genotype* genotype)
{
  if( kernel->type == kernelGRM && createType == kernelEpistaticGRM && kernel->diagonalized == false )
  {
    kernel->normalize();
    this->kernel = new Matrix(kernel->kernel);
    this->N = NULL;
    this->eigenValues = NULL;
    this->eigenVectors = NULL;
    
    copyParameters(kernel);
  
    this->genotypes = NULL;
    this->covariates = NULL;
    
    this->kernel->elementWiseMultiplication(kernel->kernel);
    this->type = kernelEpistaticGRM;
    
    /*genotype->normalizeGenotypes();
    Matrix * temp = new Matrix(genotype->genotypes);
    temp->elementWiseMultiplication(genotype->genotypes);
    Matrix * temp2 = new Matrix();
    temp2->multiply(temp, 'T', temp, 'N');
    temp2->symmetrizeTriangularMatrix();
    delete temp;
    temp = new Matrix(kernel->N);
    temp->elementWiseMultiplication(kernel->N);
    temp->symmetrizeTriangularMatrix();
    temp2->elementWiseDivision(temp);
    this->kernel->add(temp2, 0.5, -0.5); //After this substraction, the kernel have to be divided by 2. This can be arranged with the factors.
    delete temp;
    delete temp2;*/
  }
  else
  {
    misc.error("Error: An internal error was happened. A new kernel cannot be created using this types in src or dest kernels..", 0);
  }
}

Kernel::~Kernel()
{
  if(this->kernel != NULL)
  {
    delete this->kernel;
    this->kernel = NULL;
  }
  if(this->N != NULL)
  {
    delete this->N;
    this->N = NULL;
  }
  if(this->eigenValues != NULL)
  {
    delete this->eigenValues;
    this->eigenValues = NULL;
  }
  if(this->eigenVectors != NULL)
  {
    delete this->eigenVectors;
    this->eigenVectors = NULL;
  }
  if(this->genotypes != NULL)
  {
    delete this->genotypes;
    this->genotypes = NULL;
  }
  if(this->covariates != NULL)
  {
    delete this->covariates;
    this->covariates = NULL;
  }
}

void Kernel::copyParameters(Kernel *kernel)
{
  this->type = kernel->type;
  
  this->name = kernel->name;
  
  this->individuals = kernel->individuals;
  this->individualIds = kernel->individualIds;
  this->individualIdsIdx = kernel->individualIdsIdx;
  this->nIndividuals = kernel->nIndividuals;

  this->randomVarNames = kernel->randomVarNames;
  
  this->normalized = kernel->normalized;
  
  this->diagonalized = kernel->diagonalized;
  
  this->asymmetric = kernel->asymmetric;
  
  this->individualsRows = kernel->individualsRows;
  this->individualIdsRows = kernel->individualIdsRows;
  this->individualIdsIdxRows = kernel->individualIdsIdxRows;
  this->nIndividualsRows = kernel->nIndividualsRows;
  
  this->individualsCols = kernel->individualsCols;
  this->individualIdsCols = kernel->individualIdsCols;
  this->individualIdsIdxCols = kernel->individualIdsIdxCols;
  this->nIndividualsCols = kernel->nIndividualsCols;
}

void Kernel::normalize()
{
  if(this->normalized == true)
  {
    return;
  }
  
  if(this->diagonalized == true)
  {
    misc.error("Error: Diagonalized Kernel cannot be normalized.", 0);
  }
  
  if(this->type == kernelGRM)
  {
    if(this->N == NULL)
    {
      misc.error("Error: An internal error was happened. GRM cannot be normalized because normalization matrix is not present.", 0);
    }
    if(this->N->uplo != this->kernel->uplo)
    {
      misc.error("Error: An internal error was happened when normalizing a GRM. GRM and normalization matrix should have same uplo properties.", 0);
    }
    
    this->kernel->elementWiseDivision(this->N);
    
    this->normalized = true;
  }
  else
  {
    misc.error("Error: Still not supported normalization of non GRM kernels.", 0);
    this->normalized = true;
  }
}

void Kernel::denormalize()
{
  if(this->diagonalized == true)
  {
    misc.error("Error: Diagonalized Kernel cannot be unnormalized.", 0);
  }
  
  if(this->normalized == false)
  {
    return;
  }
  
  if(this->type == kernelGRM)
  {
    if(this->N == NULL)
    {
      misc.error("Error: An internal error was happened. GRM cannot be unnormalized because normalization matrix is not present.", 0);
    }
    if(this->N->uplo != this->kernel->uplo)
    {
      misc.error("Error: An internal error was happened when normalizing a GRM. GRM and normalization matrix should have same uplo properties.", 0);
    }
    
    this->kernel->elementWiseMultiplication(this->N);
    
    this->normalized = false;
  }
  else
  {
    misc.error("Error: Still not supported denormalization of non GRM kernels.", 0);
    this->normalized = true;
  }
}

Matrix * Kernel::getNormalizedKernel()
{
  normalize();
  if(this->diagonalized == false)
  {
    return this->kernel;
  }
  else
  {
    return this->eigenValues;
  }
}

void Kernel::createKernelFromDiscreteCovariates(std::string fn, int useColumn, std::set<std::string> keepIndividuals)
{
  std::vector<int> categories;
  
  this->randomVarNames.clear();
  this->individuals.clear();
  this->individualIds.clear();
  this->individualIdsIdx.clear();
  
  if( communicator->mpiRoot == true )
  {
    std::vector< std::vector<std::string> > tableWithIndividualIds;
    getTableFromFile(fn, tableWithIndividualIds, 2 + useColumn + 1);
    
    if( tableWithIndividualIds[0][0] != "FID" || tableWithIndividualIds[0][1] != "IID" )
    {
      misc.error("Error: The file [ " + fn + " ] is not properly formated. The two first header labels have to be FID and IID.", 0);
    }
    this->name = tableWithIndividualIds[0][ 2 + useColumn ];
    misc.setGetElapsedTime("ComputeRandomKernel");
    misc.message << "Starting computation of covariance from categories defined in the column " << useColumn + 1 << " of the file [ " << fn << " ]: " << this->name << "..." << std::endl;
    
    int indIdx = 0;
    std::vector<std::string> strCategories;
    std::map<std::string, int> categoriesIndex;
    int nCategories = 0;
    
    for(int idx = 1; idx<tableWithIndividualIds.size(); idx++)
    {
      Individual individual;
      
      individual.familyID = tableWithIndividualIds[idx][0];
      individual.individualID = tableWithIndividualIds[idx][1];
      individual.paternalID = "";
      individual.maternalID = "";
      individual.sex = "";
      individual.phenotype = -9;
      
      std::string key = individual.familyID + "@" + individual.individualID;
      std::string category = tableWithIndividualIds[idx][ 2 + useColumn ];
      
      if( keepIndividuals.size() != 0 && keepIndividuals.find(key) == keepIndividuals.end() )
      {
        continue;
      }
      if( category == "NA" )
      {
        continue;
      }
      
      this->individuals.push_back(individual);
      this->individualIds.push_back(key);
      if(this->individualIdsIdx.count(key) != 0)
      {
        misc.error("Error: The individual with family Id: " + individual.familyID + " and individual Id: " + individual.individualID + " appears more than one time in the file [ " + fn + " ].", 0);
      }
      this->individualIdsIdx[key] = indIdx;
      indIdx++;
      
      strCategories.push_back(category);
      if(categoriesIndex.count(category) == 0)
      {
        this->randomVarNames.push_back(category);
        categoriesIndex[ category ] = nCategories;
        nCategories++;
      }
    }
    this->nIndividuals = this->individuals.size();
    
    for(int i = 0; i < strCategories.size(); i++)
    {
      int idx = categoriesIndex[ strCategories[i] ];
      categories.push_back(idx);
    }
  }
  communicator->broadcast(&this->nIndividuals, 1);
  communicator->broadcast(this->name);
  
  if( this->nIndividuals == 0 )
  {
    misc.error("Error: Only 0 individuals kept after filtering missings from random effects defined in file [ " + fn + " ]. Please, check that in this file there are enough overlapping individuals with other data such as GRMs, covariates, and/or phenotypes.", 0);
  }
  
  misc.message << "Creating covariance matrix using categories for " << this->nIndividuals << " individuals..." << std::endl;
  
  if(this->kernel != NULL)
  {
    delete this->kernel;
  }
  if(this->N != NULL)
  {
    delete this->N;
    this->N = NULL;
  }
  if(this->eigenValues != NULL)
  {
    delete this->eigenValues;
    this->eigenValues = NULL;
  }
  if(this->eigenVectors != NULL)
  {
    delete this->eigenVectors;
    this->eigenVectors = NULL;
  }
  
  this->kernel = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->kernel->makeIntersectionMatrix(categories, 1., 0.);
  
  this->normalized = true;
  this->asymmetric = false;
  this->diagonalized = false;
  
  misc.message << "Covariance computation finished after " << misc.setGetElapsedTime("ComputeRandomKernel", true) << "." << std::endl;
}

void Kernel::createKernelFromMultipleDiscreteCovariates(std::string fn, int useColumn, std::set<std::string> keepIndividuals)
{
  std::vector<int> categories;
  
  this->randomVarNames.clear();
  this->individuals.clear();
  this->individualIds.clear();
  this->individualIdsIdx.clear();
  
  std::vector<int> nCatsByIndividual;
  std::vector<int> allCatsFlattened;
  if( communicator->mpiRoot == true )
  {
    std::vector< std::vector<std::string> > tableWithIndividualIds;
    getTableFromFile(fn, tableWithIndividualIds, 2 + useColumn + 1);
    
    if( tableWithIndividualIds[0][0] != "FID" || tableWithIndividualIds[0][1] != "IID" )
    {
      misc.error("Error: The file [ " + fn + " ] is not properly formated. The two first header labels have to be FID and IID.", 0);
    }
    this->name = tableWithIndividualIds[0][ 2 + useColumn ];
    misc.setGetElapsedTime("ComputeRandomKernel");
    misc.message << "Starting computation of covariance from multiple categories defined in the column " << useColumn + 1 << " of the file [ " << fn << " ]: " << this->name << "..." << std::endl;
    
    int indIdx = 0;
    std::map<std::string, int> categoriesIndex;
    int nCategories = 0;
    for(int idx = 1; idx<tableWithIndividualIds.size(); idx++)
    {
      Individual individual;
      
      individual.familyID = tableWithIndividualIds[idx][0];
      individual.individualID = tableWithIndividualIds[idx][1];
      individual.paternalID = "";
      individual.maternalID = "";
      individual.sex = "";
      individual.phenotype = -9;
      
      std::string key = individual.familyID + "@" + individual.individualID;
      std::string categories = tableWithIndividualIds[idx][ 2 + useColumn ];
      
      if( keepIndividuals.size() != 0 && keepIndividuals.find(key) == keepIndividuals.end() )
      {
        continue;
      }
      if( categories == "NA" )
      {
        continue;
      }
      
      this->individuals.push_back(individual);
      this->individualIds.push_back(key);
      if(this->individualIdsIdx.count(key) != 0)
      {
        misc.error("Error: The individual with family Id: " + individual.familyID + " and individual Id: " + individual.individualID + " appears more than one time in the file [ " + fn + " ].", 0);
      }
      this->individualIdsIdx[key] = indIdx;
      indIdx++;
      
      
      std::vector<std::string> catsTemp = splitString(categories, ",");
      std::vector<std::string> catsTempIdxs;
      for(int catidx = 0; catidx < catsTemp.size(); catidx++)
      {
        std::string category = catsTemp[ catidx ];
        if(categoriesIndex.count(category) == 0)
        {
          this->randomVarNames.push_back(category);
          categoriesIndex[ category ] = nCategories;
          nCategories++;
        }
        allCatsFlattened.push_back( categoriesIndex[ category ] );
      }
      nCatsByIndividual.push_back(catsTemp.size());
      //indCategories.push_back(catsTemp);
    }
    this->nIndividuals = this->individuals.size();
  }
  communicator->broadcast(&this->nIndividuals, 1);
  communicator->broadcast(this->name);
  
  communicator->broadcast(nCatsByIndividual);
  communicator->broadcast(allCatsFlattened);

  std::vector< std::set<int> > catsByIndividual;
  std::vector< double > normByIndividual;
  std::vector< int > indIdxs;
  int shift = 0;
  for(int indidx = 0; indidx < nCatsByIndividual.size(); indidx++)
  {
    std::set<int> temp;
    for(int idx = shift; idx < (shift + nCatsByIndividual[indidx]); idx++)
    {
      temp.insert(allCatsFlattened[idx]);
    }
    shift += nCatsByIndividual[indidx];
    catsByIndividual.push_back(temp);
    normByIndividual.push_back(sqrt(double(temp.size())));
    indIdxs.push_back(indidx);
  }
  
  if( this->nIndividuals == 0 )
  {
    misc.error("Error: Only 0 individuals kept after filtering missings from random effects defined in file [ " + fn + " ]. Please, check that in this file there are enough overlapping individuals with other data such as GRMs, covariates, and/or phenotypes.", 0);
  }
  
  misc.message << "Creating covariance matrix using categories for " << this->nIndividuals << " individuals..." << std::endl;
  
  if(this->kernel != NULL)
  {
    delete this->kernel;
  }
  if(this->N != NULL)
  {
    delete this->N;
    this->N = NULL;
  }
  if(this->eigenValues != NULL)
  {
    delete this->eigenValues;
    this->eigenValues = NULL;
  }
  if(this->eigenVectors != NULL)
  {
    delete this->eigenVectors;
    this->eigenVectors = NULL;
  }
  
  this->kernel = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->nIndividuals, this->nIndividuals);
  this->kernel->fillWithConstant(0.);
  
  //This part should be integrated in matrix.cpp in the future. It assumes a particular matrix distribution.
  int * iridx = this->kernel->scatterVectorRet(&(indIdxs[0]), row);
  int * icidx = this->kernel->scatterVectorRet(&(indIdxs[0]), column);
  
  for(int c = 0; c < this->kernel->nCols; c++)
  {
    for(int r = 0; r < this->kernel->nRows; r++)
    {
      int i1 = iridx[r];
      int i2 = icidx[c];
      std::set<int> intersectElements;
      std::set_intersection(catsByIndividual[i1].begin(),catsByIndividual[i1].end(),catsByIndividual[i2].begin(),catsByIndividual[i2].end(), std::inserter(intersectElements,intersectElements.begin()));
      
      this->kernel->m[c*this->kernel->nRows + r] = double(intersectElements.size())/(normByIndividual[i1]*normByIndividual[i2]);
    }
  }
  this->kernel->symmetric = true;
  this->kernel->uplo = 'B';
  
  delete [] iridx;
  delete [] icidx;
  
  this->normalized = true;
  this->asymmetric = false;
  this->diagonalized = false;
  
  misc.message << "Covariance computation finished after " << misc.setGetElapsedTime("ComputeRandomKernel", true) << "." << std::endl;
}

void Kernel::createKernelFromContinuousCovariates(std::string fn)
{
}

void Kernel::createKernelSquaredExponential(std::string fn, std::set<std::string> keepIndividuals)
{
  misc.setGetElapsedTime("ComputeRandomKernel");
  misc.message << "Starting computation of squared exponential covariance from coordinates defined in the file [ " << fn << " ]..." << std::endl;
  
  this->randomVarNames.clear();
  this->individuals.clear();
  this->individualIds.clear();
  this->individualIdsIdx.clear();
  this->nIndividuals = 0;
  
  this->randomVarNames.clear();
  
  int nIndividualsInFile = 0;
  
  std::vector< std::vector<double> > vectorsCoordinates;
  if( communicator->mpiRoot == true )
  {
    this->randomVarNames.push_back("distance");
    
    int fNColumns = getNumberOfFileColumns(fn);
    if(fNColumns < 3)
    {
      misc.error("Error: Error when computing a squared exponential kernel. The file [ " + fn + " ] does not have enough columns of data.", 0);
    }
    std::vector< std::vector<std::string> > rawData;
    getTableFromFile(fn, rawData, fNColumns);
    if( rawData[0][0] != "FID" || rawData[0][1] != "IID" )
    {
      misc.error("Error: The file [ " + fn + " ] is not properly formated. The two first header labels have to be FID and IID.", 0);
    }
    this->name = rawData[0][ 2 ];
    
    vectorsCoordinates = std::vector< std::vector<double> >(fNColumns - 2, std::vector<double>());
    
    nIndividualsInFile = rawData.size();
    
    for( int ir = 1; ir < rawData.size(); ir++)
    {
      Individual individual;
      
      individual.familyID = rawData[ir][0];
      individual.individualID = rawData[ir][1];
      individual.paternalID = "";
      individual.maternalID = "";
      individual.sex = "";
      individual.phenotype = -9;
      
      std::string key = individual.familyID + "@" + individual.individualID;
      
      bool skipIndividual = false;
      for(int ic = 2; ic<fNColumns; ic++)
      {
        if(rawData[ir][ic] == "NA")
        {
          skipIndividual = true;
        }
      }
      if( keepIndividuals.size() != 0 && keepIndividuals.find(key) == keepIndividuals.end() )
      {
        skipIndividual = true;
      }
      if( skipIndividual == true )
      {
        continue;
      }
      
      this->individuals.push_back(individual);
      this->individualIds.push_back(key);
      if(this->individualIdsIdx.count(key) != 0)
      {
        misc.error("Error: The individual with family Id: " + individual.familyID + " and individual Id: " + individual.individualID + " appears more than one time in the file [ " + fn + " ].", 0);
      }
      this->individualIdsIdx[key] = this->nIndividuals;
      this->nIndividuals++;
      
      for(int ic = 2; ic<fNColumns; ic++)
      {
        double value = string2Number<double>(rawData[ir][ic], "Error: The individual " + rawData[ir][0] + " " + rawData[ir][1] + " has not a valid value in one of their covariables: " + rawData[ir][ic]);
        vectorsCoordinates[ic - 2].push_back(value);
      }
    }
  }
  communicator->broadcast(&(this->nIndividuals));
  int nDimensions = vectorsCoordinates.size();
  communicator->broadcast(&nDimensions);
  
  if(this->nIndividuals == 0)
  {
    misc.error("Error: The file [ " + fn + " ] is empty, or all individuals have at least a missing value.", 0);
  }
  
  misc.message << "Read data for " << nIndividualsInFile << " individuals. " << this->nIndividuals << " will be kept." << std::endl;
  
  if(this->kernel != NULL)
  {
    delete this->kernel;
  }
  if(this->N != NULL)
  {
    delete this->N;
    this->N = NULL;
  }
  if(this->eigenValues != NULL)
  {
    delete this->eigenValues;
    this->eigenValues = NULL;
  }
  if(this->eigenVectors != NULL)
  {
    delete this->eigenVectors;
    this->eigenVectors = NULL;
  }
  this->kernel = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->nIndividuals, this->nIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->kernel->fillWithConstant(0);
  for(int id = 0; id<nDimensions; id++)
  {
    std::vector<double> tempVector;
    if(communicator->mpiRoot == true)
    {
      tempVector = vectorsCoordinates[id];
    }
    Matrix * temp = new Matrix();
    temp->makeDifferenceMatrix(tempVector, tempVector);
    temp->elementWiseMultiplication(temp);
    this->kernel->add(temp);
    delete temp;
  }
  
  this->kernel->symmetric = true;
  this->kernel->uplo = 'B';
  
  this->normalized = true;
  this->asymmetric = false;
  this->diagonalized = false;
  
  misc.message << "Covariance computation finished after " << misc.setGetElapsedTime("ComputeRandomKernel", true) << "." << std::endl;
}

void Kernel::writeKernel(std::string f)
{
  if(this->type == kernelGRM)
  { 
    writeGRM(f);
  }
  else
  {
    misc.error("Error: Still not supported storing of non GRM kernels.", 0);
  }
}

void Kernel::writeGRM(std::string f)
{
  std::ofstream file;
  
  if(this->asymmetric == true)
  {
    misc.error("Error: An internal error was happened. This operation can not be performed on an asymmetric GRM", 0);
  }
  if(this->normalized == false)
  {
    misc.error("Error: An internal error was happened. GRM can only be stored in normalized form.", 0);
  }
  if(this->N == NULL && this->diagonalized == false)
  {
    misc.error("Error: An internal error was happened. GRM can only be stored when N information is stored.", 0);
  }
  
  misc.setGetElapsedTime("GRMWrite");
  if(this->diagonalized == false)
  {
    misc.message << "Writing GRM on [ " << f << ".grm.(ids/snps/dat) ]..." << std::endl;
  }
  else
  {
    misc.message << "Writing GRM on [ " << f << ".grm.(ids/snps/dat/diag) ]..." << std::endl;
  }
  
  if(communicator->mpiRoot)
  {
    std::string fname;
    
    //Write Individual info
    fname = f + ".grm.ids";
    file.open(fname.c_str(), std::ofstream::out);
    for(int i = 0; i<this->nIndividuals; i++)
    {
      file << this->individuals[i].familyID << " " << this->individuals[i].individualID << std::endl;
    }
    file.close();
    
    //Write SNP info
    fname = f + ".grm.snps";
    file.open(fname.c_str(), std::ofstream::out);
    for(int i = 0; i<this->randomVarNames.size(); i++)
    {
      file << this->randomVarNames[i] << std::endl;
    }
    file.close();
  }

  //Define header
  unsigned char *header = new unsigned char [14];
  header[0] = 'G';
  header[1] = 'R';
  header[2] = 'M';
  header[3] = '\0';
  header[4] = 0x5A; //Magic number 1
  header[5] = 0x99; //Magic number 2
  header[6] = 0x2; //Version
  header[7] = 0x1; //1 if stored doubles
  header[8] = sizeof(double); //double size
  if(this->diagonalized == false)
  {
    header[9] = 0x1; //1 normalized 2 non normalized
  }
  else
  {
    header[9] = 0x3; //3 diagonalized
  }
  
  header[10] = 0x0; //non used
  header[11] = 0x0; //non used
  header[12] = 0x0; //non used
  header[13] = 0x0; //non used
  
  if(this->diagonalized == false)
  {
    Matrix * grmPacked = new Matrix(cyclicDistribution);
    grmPacked->packMatrices(this->kernel, this->N);
    if(options.useMPIForWriting == true)
    {
      grmPacked->writeMatrixMPI(f + ".grm.dat", (char*)header, 14);
    }
    else
    {
      if(communicator->mpiRoot)
      {
        std::string fname = f + ".grm.dat";
        file.open(fname.c_str(), std::ofstream::out | std::ofstream::binary);
        file.write((char*)header, 10);
      }
      grmPacked->writeMatrixFile(file);
      if(communicator->mpiRoot)
      {
        file.close();
      }
    }
    delete grmPacked;
  }
  else
  {
    this->eigenVectors->writeMatrixMPI(f + ".grm.dat", (char*)header, 14);
    //Write diagonal
    if(communicator->mpiRoot)
    {
      std::string fname = f + ".grm.diag";
      file.open(fname.c_str(), std::ofstream::out | std::ofstream::binary);
      file.write((char*)this->eigenValues->m,sizeof(double)*this->eigenValues->nGlobRows);
      file.close();
    }
  }
  
  delete [] header;
  
  misc.message << "GRM stored after " << misc.setGetElapsedTime("GRMWrite", true) << std::endl;
  
}

void Kernel::readKernel(std::string f)
{
  if(this->type == kernelGRM)
  { 
    readGRM(f);
  }
  else
  {
    misc.error("Error: Still not supported reading of non GRM kernels.", 0);
  }
}

void Kernel::readGRM(std::string f)
{
  std::string fname;
  std::ifstream file;
  
  if(this->asymmetric == true)
  {
    misc.error("Error: An internal error was happened. This operation can not be performed on an asymmetric GRM", 0);
  }
  this->normalized = true;
  this->asymmetric = false;
  
  misc.setGetElapsedTime("GRMRead");
  misc.message << "Reading GRM from [ " << f << ".grm.(ids/snps/dat) ]..." << std::endl;
  
  if(communicator->mpiRoot)
  {
    std::string line;
    
    //read Individual info
    fname = f + ".grm.ids";
    misc.message << "Reading Individuals Ids data from GRM file [ " << fname << " ] ..." << std::endl;
    misc.checkFileExists(fname);
    file.open(fname.c_str(), std::ifstream::in);
    this->individuals.clear();
    this->individualIds.clear();
    this->individualIdsIdx.clear();
    int idx = 0;
    while(getline(file,line))
    {
      std::istringstream sstemp(line); //->Comprova que te 2 elements.
      
      Individual individual;
      
      sstemp >> individual.familyID;
      sstemp >> individual.individualID;
      individual.paternalID = "";
      individual.maternalID = "";
      individual.sex = "";
      individual.phenotype = -9;
      
      this->individuals.push_back(individual);
      std::string key = individual.familyID + "@" + individual.individualID;
      this->individualIds.push_back(key);
      if(this->individualIdsIdx.count(key) != 0)
      {
	misc.error("Error: The individual with family Id: " + individual.familyID + " and individual Id: " + individual.individualID + " appears more than one time in the genotypes.", 0);
      }
      this->individualIdsIdx[key] = idx;
      idx++;
    }
    this->nIndividuals = this->individuals.size();
    file.close();
    misc.message << this->nIndividuals << " individuals found." << std::endl;
    
    //Read SNP info
    fname = f + ".grm.snps";
    misc.message << "Reading SNP Ids data from GRM file [ " << fname << " ] ..." << std::endl;
    misc.checkFileExists(fname);
    file.open(fname.c_str(), std::ifstream::in);
    this->randomVarNames.clear();
    while(getline(file,line))
    {
      std::istringstream sstemp(line);
      std::string SNPId;
      sstemp >> SNPId;
      this->randomVarNames.push_back(SNPId);
    }
    file.close();
    
    //Read data
    fname = f + ".grm.dat";
    misc.checkFileExists(fname);
    file.open(fname.c_str(), std::ifstream::in | std::ifstream::binary);
    
    //Read header
    unsigned char *header = new unsigned char [10];
    file.read((char*)header, 10);
    if (
      header[0] != 'G' ||
      header[1] != 'R' ||
      header[2] != 'M' ||
      header[3] != '\0' ||
      header[4] != 0x5A ||
      header[5] != 0x99 ||
      header[6] != 0x2 ||
      header[7] != 0x1
    )
    {
      misc.error("Error: The file [ " + fname + " ] do not have a proper header or is not a valid file.", 0);
    }
    if( header[8] != sizeof(double) )
    {
      misc.error("Error: Data types in file [ " + fname + " ] differ with this system datatypes. Are you loading a grm computed in a different system?", 0);
    }
    if( header[9] != 0x1 && header[9] != 0x3 )
    {
      misc.error("Error: Error in file [ " + fname + " ]. Non normalized matrices still not supported.", 0);
    }
    if(header[9] == 0x3)
    {
      this->diagonalized = true;
      misc.message << "Reading diagonalized GRM data from file [ " << f << ".grm.dat/diag ] ..." << std::endl;
    }
    else
    {
      this->diagonalized = false;
      misc.message << "Reading GRM data from file [ " << f << ".grm.dat ] ..." << std::endl;
    }
    delete [] header;
    file.close();
  }
  
  communicator->broadcast(&(this->diagonalized));
  
  communicator->broadcast(&this->nIndividuals, 1);
  communicator->barrier();
  
  if(this->kernel != NULL)
  {
    delete this->kernel;
    this->kernel = NULL;
  }
  if(this->N != NULL)
  {
    delete this->N;
    this->N = NULL;
  }
  if(this->eigenValues != NULL)
  {
    delete this->eigenValues;
    this->eigenValues = NULL;
  }
  if(this->eigenVectors != NULL)
  {
    delete this->eigenVectors;
    this->eigenVectors = NULL;
  }
  
  if(this->diagonalized == false)
  {
    this->kernel = new Matrix(cyclicDistribution);
    this->N = new Matrix(cyclicDistribution);

    Matrix * grmPacked = new Matrix(cyclicDistribution, this->nIndividuals + 1, this->nIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    grmPacked->readMatrixMPI(f + ".grm.dat", 14);
    grmPacked->unpackMatrices(this->kernel, this->N, true);
    delete grmPacked;
  }
  else
  {
    this->eigenVectors = new Matrix(cyclicDistribution, this->nIndividuals, this->nIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    this->eigenVectors->readMatrixMPI(f + ".grm.dat", 14);
    //Read diagonal elements
    this->eigenValues = new Matrix(diagonalDistribution, this->nIndividuals, this->nIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    if(communicator->mpiRoot == true)
    {
      fname = f + ".grm.diag";
      misc.checkFileExists(fname);
      file.open(fname.c_str(),  std::ifstream::in | std::ifstream::binary);
      file.seekg (0, file.end);
      int fileSize = file.tellg();
      file.seekg (0, file.beg);
      if(fileSize != sizeof(double)*this->eigenValues->nGlobRows)
      {
        misc.error("Error: The file [ " + f + ".grm.diag ] is not properly formated.", 0);
      }
      file.read((char*)this->eigenValues->m,sizeof(double)*this->eigenValues->nGlobRows);
      file.close();
    }
  }
  
  misc.message << "GRM loaded after " << misc.setGetElapsedTime("GRMRead", true) << std::endl;
}

void Kernel::readGCTAGRM(std::string f)
{
#if defined(BOOSTLIB) && defined(ZLIB)
  if(this->asymmetric == true)
  {
    misc.error("Error: An internal error was happened. This operation can not be performed on an asymmetric GRM", 0);
  }
  this->normalized = true;
  this->asymmetric = false;
  
  misc.setGetElapsedTime("GRMRead");
  misc.message << "Reading GRM from [ " << f << ".grm.(id/gz) ]..." << std::endl;
  
  if(communicator->mpiRoot)
  {
    std::string line;
    
    //read Individual info
    misc.message << "Reading Individuals Ids data from GRM file [ " << f + ".grm.id" << " ] ..." << std::endl;
    
    std::vector< std::vector<std::string> > table;
    getTableFromFile(f + ".grm.id", table, 2);
    
    this->individuals.clear();
    this->individualIds.clear();
    this->individualIdsIdx.clear();
    for(int iline = 0; iline < table.size(); iline++)
    {
      Individual individual;
      
      individual.familyID = table[iline][0];
      individual.individualID = table[iline][1];
      individual.paternalID = "";
      individual.maternalID = "";
      individual.sex = "";
      individual.phenotype = -9;
      
      this->individuals.push_back(individual);
      std::string key = individual.familyID + "@" + individual.individualID;
      this->individualIds.push_back(key);
      if(this->individualIdsIdx.count(key) != 0)
      {
        misc.error("Error: The individual with family Id: " + individual.familyID + " and individual Id: " + individual.individualID + " appears more than one time in the genotypes.", 0);
      }
      this->individualIdsIdx[key] = iline;
    }
    this->nIndividuals = this->individuals.size();
    misc.message << this->nIndividuals << " individuals found." << std::endl;
    
    //Read SNP info
    misc.message << "No SNPs data for this type of file. Skipping reading SNPs." << std::endl;
    this->randomVarNames.clear();
    //this->randomVarNames.push_back("nosnpdata");
    
    
    this->diagonalized = false;
  }
  
  misc.message << "Reading GRM data from file [ " << f << ".grm.gz ] ..." << std::endl;
  
  
  communicator->broadcast(&(this->diagonalized));
  
  communicator->broadcast(&this->nIndividuals, 1);
  communicator->barrier();
  
  if(this->kernel != NULL)
  {
    delete this->kernel;
    this->kernel = NULL;
  }
  if(this->N != NULL)
  {
    delete this->N;
    this->N = NULL;
  }
  if(this->eigenValues != NULL)
  {
    delete this->eigenValues;
    this->eigenValues = NULL;
  }
  if(this->eigenVectors != NULL)
  {
    delete this->eigenVectors;
    this->eigenVectors = NULL;
  }
  
  double * gGRM = NULL;
  double * gN = NULL;
  if(communicator->mpiRoot)
  {
    gGRM = new double [this->nIndividuals*this->nIndividuals];
    gN = new double [this->nIndividuals*this->nIndividuals];
    
    std::string fname = f + ".grm.gz";
    std::ifstream file(fname.c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);
    
    std::istream is(&in);
    
    std::string line;
    int nLoadedValues = 0;
    while(std::getline(is, line))
    {
      std::istringstream sstemp(line);
      int idx1;
      int idx2;
      double normval;
      double grmval;
      if( (sstemp >> idx1).fail() )
      {
        misc.error("Error: An error has happened when getting a matrix index value from the file [ " + f + ".grm.gz ]. The file may not be properly formated.", 0);
      }
      if( (sstemp >> idx2).fail() )
      {
        misc.error("Error: An error has happened when getting a matrix index value from the file [ " + f + ".grm.gz ]. The file may not be properly formated.", 0);
      }
      if( (sstemp >> normval).fail() )
      {
        misc.error("Error: An error has happened when getting a number of shared SNPs value from the file [ " + f + ".grm.gz ]. The file may not be properly formated.", 0);
      }
      if( (sstemp >> grmval).fail() )
      {
        misc.error("Error: An error has happened when getting a relationship value from the file [ " + f + ".grm.gz ]. The file may not be properly formated.", 0);
      }
      
      idx1--;
      idx2--;
      
      if(idx1 < 0 || idx1 >= this->nIndividuals || idx2 < 0 || idx2 >= this->nIndividuals)
      {
        misc.error("Error: An error has happened when getting a matrix index value from the file [ " + f + ".grm.gz ]. Invalid index(s). The file may not be properly formated.", 0);
      }
      
      gGRM[ idx2*this->nIndividuals + idx1 ] = grmval;
      gGRM[ idx1*this->nIndividuals + idx2 ] = grmval;
      
      gN[ idx2*this->nIndividuals + idx1 ] = normval;
      gN[ idx1*this->nIndividuals + idx2 ] = normval;
      
      nLoadedValues++;
    }
    
    if( 2*nLoadedValues < (this->nIndividuals*(this->nIndividuals + 1)) )
    {
      misc.error("Error: An error has happened when loading the file [ " + f + ".grm.gz ]. The number of entries is not enough to define a complete matrix.", 0);
    }
    if( 2*nLoadedValues != (this->nIndividuals*(this->nIndividuals + 1)) && nLoadedValues != (this->nIndividuals*this->nIndividuals) )
    {
      misc.message << "WARNING: Unexpected number of entries on file [ " + f + ".grm.gz ]: " + i2s(nLoadedValues) + ". Please, check that the file is properly formated." << std::endl;
    }
  }

  this->kernel = new Matrix(cyclicDistribution, this->nIndividuals, this->nIndividuals);
  this->N = new Matrix(cyclicDistribution, this->nIndividuals, this->nIndividuals);

  this->kernel->scatterMatrix(gGRM);
  this->N->scatterMatrix(gN);
  
  this->kernel->symmetric = true;
  this->kernel->uplo = 'B';
  
  this->N->symmetric = true;
  this->N->uplo = 'B';
  
  misc.message << "GRM loaded after " << misc.setGetElapsedTime("GRMRead", true) << std::endl;
  
  if(communicator->mpiRoot)
  {
    delete [] gGRM;
    delete [] gN;
  }
#else
  misc.error("Error: The current version of DISSECT is compiled without zlib and boost support. These libraries are required for loading GCTA gz GRM matrices. Please, recompile it using these libraries or ask for support.", 0);
#endif

}

void Kernel::filterIndividuals(std::vector<std::string> & keepIndividualIds, bool filterN)
{
  int *keepIndxs;

  if(this->asymmetric == true)
  {
    misc.error("Error: An internal error was happened. This operation cannot be performed on an asymmetric GRM.", 0);
  }
  if(filterN == true && this->N == NULL && this->type == kernelGRM)
  {
    misc.error("Error: An internal error was happened when filtering a Kernel. N matrix is not stored and cannot be filtered.", 0);
  }
  
  if( this->type != kernelGRM )
  {
    filterN = false;
  }
  
  //If the individuals not changed, return.
  int noChangesFlag = 0;
  if(communicator->mpiRoot)
  {
    noChangesFlag = ((this->individualIds == keepIndividualIds)?1:0);
  }
  communicator->broadcast(&noChangesFlag, 1);
  if(noChangesFlag == 1)
  {
    return;
  }
  
  if(this->diagonalized == true)
  {
    misc.error("Error: An error was happened. Diagonalized GRM cannot be filtered. Please, check that there are the same non-missing individuals in your grm, genotype, phenotypes, and covariates files.", 0);
  }
  
  this->kernel->symmetrizeTriangularMatrix();
  if(this->N != NULL)
  {
    this->N->symmetrizeTriangularMatrix();
  }
  
  if(communicator->mpiRoot)
  {
    std::map<std::string, int> newIndividualIdsIdx;
    std::vector<Individual> newIndividuals;
    
    keepIndxs = new int [keepIndividualIds.size()];
    int previousOldIdx = -1;
    
    for(int r=0; r<keepIndividualIds.size(); r++) //Is it better access rows contiguously and search in map each step? Or inverse?
    {
      if(this->individualIdsIdx.count(keepIndividualIds[r]) == 0)
      {
	misc.error("Error: An internal error was happened. The individual '" + keepIndividualIds[r]  + "' is not in the GRM matrix.", 0);
      }
      if(newIndividualIdsIdx.count(keepIndividualIds[r]) != 0)
      {
        misc.error("Error: An internal error was happened. The individual '" + keepIndividualIds[r]  + "' id is repeated.", 0);
      }
      int oldIdx = this->individualIdsIdx[keepIndividualIds[r]];
      
      keepIndxs[r] = oldIdx;
      newIndividualIdsIdx[keepIndividualIds[r]] = r;
      newIndividuals.push_back(this->individuals[oldIdx]);
      previousOldIdx = oldIdx;
    }
    this->nIndividuals = keepIndividualIds.size();
    this->individualIds = keepIndividualIds;
    this->individualIdsIdx = newIndividualIdsIdx;
    this->individuals = newIndividuals;
  }
  
  communicator->broadcast(&this->nIndividuals, 1);
  if(!communicator->mpiRoot)
  {
    keepIndxs = new int [this->nIndividuals];
  }
  communicator->broadcast(keepIndxs, this->nIndividuals);
  
  std::map<int, int> rowsOriginDestination;
  std::map<int, int> colsOriginDestination;
  for(int i = 0; i<this->nIndividuals; i++)
  {
    rowsOriginDestination[ keepIndxs[i] ] = i;
    colsOriginDestination[ keepIndxs[i] ] = i;
  }
  
  Matrix * newKernelMatrix = new Matrix(cyclicDistribution);
  this->kernel->generalResorting(newKernelMatrix, rowsOriginDestination, colsOriginDestination, true);
  newKernelMatrix->symmetric = this->kernel->symmetric;
  newKernelMatrix->uplo = this->kernel->uplo;
  delete this->kernel;
  this->kernel = newKernelMatrix;
  
  if(filterN)
  {
    Matrix * newNMatrix = new Matrix(cyclicDistribution);
    this->N->generalResorting(newNMatrix, rowsOriginDestination, colsOriginDestination, true);
    newNMatrix->symmetric = this->N->symmetric;
    newNMatrix->uplo = this->N->uplo;
    delete this->N;
    this->N = newNMatrix;
  }
  else
  {
    if(this->N != NULL)
    {
      delete this->N;
      this->N = NULL;
    }
  }
  
  delete [] keepIndxs;
}

void Kernel::filterIndividualsAsymmetric(std::vector<std::string> & keepIndividualIdsRows, std::vector<std::string> & keepIndividualIdsCols, bool filterN)
{
  int *keepIndxsRows;
  int *keepIndxsCols;
  
  if(this->asymmetric == true)
  {
    misc.error("Error: An internal error was happened. This operation can not be performed on an asymmetric GRM", 0);
  }
  if(filterN == true && this->N == NULL)
  {
    misc.error("Error: An internal error was happened when filtering a Kernel. N matrix is not stored and cannot be filtered.", 0);
  }
  
  if( this->type != kernelGRM )
  {
    filterN = false;
  }
  
  if( misc.gt( (keepIndividualIdsRows == this->individualIds) && (keepIndividualIdsCols == this->individualIds) ) )
  {
    this->nIndividualsRows = this->nIndividuals;
    this->individualIdsRows = this->individualIds;
    this->individualIdsIdxRows = this->individualIdsIdx;
    this->individualsRows = this->individuals;

    this->nIndividualsCols = this->nIndividuals;
    this->individualIdsCols = this->individualIds;
    this->individualIdsIdxCols = this->individualIdsIdx;
    this->individualsCols = this->individuals;

    this->nIndividuals = 0;
    this->individualIds.clear();
    this->individualIdsIdx.clear();
    this->individuals.clear();
    
    communicator->broadcast(&this->nIndividuals, 1);
    communicator->broadcast(&this->nIndividualsRows, 1);
    communicator->broadcast(&this->nIndividualsCols, 1);
    
    this->asymmetric = true;
    
    return;
  }
  
  if(this->diagonalized == true)
  {
    misc.error("Error: An internal error was happened. Diagonalized grm cannot be filtered.", 0);
  }
  
  this->asymmetric = true;
  this->kernel->symmetrizeTriangularMatrix();
  this->kernel->symmetric = false;
  if(this->N != NULL)
  {
    this->N->symmetrizeTriangularMatrix();
    this->N->symmetric = false;
  }
  
  if(communicator->mpiRoot)
  {
    for(int i=0; i<2; i++) //0 = rows, 1 = cols
    {
      std::vector<std::string> keepIndividualIds;
      if(i == 0)
      {
        keepIndividualIds = keepIndividualIdsRows;
      }
      else
      {
        keepIndividualIds = keepIndividualIdsCols;
      }
      
      std::map<std::string, int> newIndividualIdsIdx;
      std::vector<Individual> newIndividuals;
      
      int * keepIndxs = new int [keepIndividualIds.size()];
      int previousOldIdx = -1;
      
      for(int r=0; r<keepIndividualIds.size(); r++) //Is it better access rows contiguously and search in map each step? Or inverse?
      {
        if(this->individualIdsIdx.count(keepIndividualIds[r]) == 0)
        {
          misc.error("Error: An internal error was happened. The individual '" + keepIndividualIds[r]  + "' is not in the GRM matrix.", 0);
        }
        if(newIndividualIdsIdx.count(keepIndividualIds[r]) != 0)
        {
          misc.error("Error: An internal error was happened. The individual '" + keepIndividualIds[r]  + "' id is repeated.", 0);
        }
        int oldIdx = this->individualIdsIdx[keepIndividualIds[r]];
        
        keepIndxs[r] = oldIdx;
        newIndividualIdsIdx[keepIndividualIds[r]] = r;
        newIndividuals.push_back(this->individuals[oldIdx]);
        previousOldIdx = oldIdx;
      }
      
      if(i == 0)
      {
        this->nIndividualsRows = keepIndividualIds.size();
        this->individualIdsRows = keepIndividualIds;
        this->individualIdsIdxRows = newIndividualIdsIdx;
        this->individualsRows = newIndividuals;
        keepIndxsRows = keepIndxs;
      }
      else
      {
        this->nIndividualsCols = keepIndividualIds.size();
        this->individualIdsCols = keepIndividualIds;
        this->individualIdsIdxCols = newIndividualIdsIdx;
        this->individualsCols = newIndividuals;
        keepIndxsCols = keepIndxs;
      }
    } // End for(int i=0; i<2; i++)
    this->nIndividuals = 0;
    this->individualIds.clear();
    this->individualIdsIdx.clear();
    this->individuals.clear();
  } // End communicator->mpiRoot
  
  communicator->broadcast(&this->nIndividuals, 1);
  communicator->broadcast(&this->nIndividualsRows, 1);
  communicator->broadcast(&this->nIndividualsCols, 1);
  if(!communicator->mpiRoot)
  {
    keepIndxsRows = new int [this->nIndividualsRows];
    keepIndxsCols = new int [this->nIndividualsCols];
  }
  communicator->broadcast(keepIndxsRows, this->nIndividualsRows);
  communicator->broadcast(keepIndxsCols, this->nIndividualsCols);
  
  std::map<int, int> rowsOriginDestination;
  std::map<int, int> colsOriginDestination;
  for(int i = 0; i<this->nIndividualsRows; i++)
  {
    rowsOriginDestination[ keepIndxsRows[i] ] = i;
  }
  for(int i = 0; i<this->nIndividualsCols; i++)
  {
    colsOriginDestination[ keepIndxsCols[i] ] = i;
  }
  
  Matrix * newKernelMatrix = new Matrix(cyclicDistribution);
  this->kernel->generalResorting(newKernelMatrix, rowsOriginDestination, colsOriginDestination, true);
  delete this->kernel;
  this->kernel = newKernelMatrix;
  
  if(filterN)
  {
    Matrix * newNMatrix = new Matrix(cyclicDistribution);
    this->N->generalResorting(newNMatrix, rowsOriginDestination, colsOriginDestination, true);
    delete this->N;
    this->N = newNMatrix;
  }
  else
  {
    if(this->N != NULL)
    {
      delete this->N;
      this->N = NULL;
    }
  }
  
  delete [] keepIndxsRows;
  delete [] keepIndxsCols;
}

void Kernel::replaceIndividualIds(std::vector<Individual> & newIndividuals)
{
  if( this->asymmetric == true )
  {
    misc.error("Error: An internal error was happened. Individual Ids cannot be replaced on an asymmetric kernel.", 0);
  }
  
  this->individuals.clear();
  this->individualIds.clear();
  this->individualIdsIdx.clear();
  
  if(communicator->mpiRoot)
  {
    if(this->nIndividuals != newIndividuals.size())
    {
      misc.error("Error: An internal error was happened. Individual Ids cannot be replaced. Size of new individual Ids differ from Kernel dimensions.", 0);
    }
      
    this->individuals = newIndividuals;
      
    for(int i = 0; i<this->nIndividuals; i++)
    {
      std::string key = this->individuals[i].familyID + "@" + this->individuals[i].individualID;
      this->individualIds.push_back(key);
      if(this->individualIdsIdx.count(key) != 0)
      {
        misc.error("Error: The individual with family Id: " + this->individuals[i].familyID + " and individual Id: " + this->individuals[i].individualID + " appears more than one time in the Kernel replacement keys.", 0);
      }
      this->individualIdsIdx[key] = i;
    }
  }
}

void Kernel::addKernels(double scalingFactor1, Kernel *kernel1, double scalingFactor2, Kernel *kernel2)
{
  if(this->type == kernelGRM)
  { 
    addGRMs(scalingFactor1, kernel1, scalingFactor2, kernel2);
  }
  else
  {
    misc.error("Error: Adding non GRM kernels still not supported.", 0);
  }
}

void Kernel::addGRMs(double scalingFactor1, Kernel *grm1, double scalingFactor2, Kernel *grm2)
{
  if(grm1 == this || grm2 == this)
  {
    misc.error("Error: An internal error was happened when adding two GRMs. GRMs to add must be different instances than the resulting GRM.", 0);
  }
  if(grm2 == NULL)
  {
    grm2 = this;
  }
  if( (scalingFactor1 != 1. && scalingFactor1 != -1.) || (scalingFactor2 != 1. && scalingFactor2 != -1.) )
  {
    misc.error("Error: An internal error was happened. Scaling factor when adding or substracting GRMs can only be 1 or -1.", 0);
  }
  if(grm1->diagonalized == true || grm2->diagonalized == true || this->diagonalized == true)
  {
    misc.error("Error: An internal error was happened. Diagonalized grm's cannot be added.", 0);
  }
  if(grm1->N == NULL || grm2->N == NULL)
  {
    misc.error("Error: An internal error was happened. GRMs can not be added without the N matrix.", 0);
  }
  if(grm1->asymmetric != grm2->asymmetric)
  {
    misc.error("Error: An internal error was happened. Symmetric GRM can not be added to a asymetric GRM.", 0);
  }
  if(communicator->mpiRoot)
  {
    if(this->asymmetric == false)
    {
      if( (grm1->individualIds != grm2->individualIds) ||
          (grm1->individualIdsIdx != grm2->individualIdsIdx) ||
          (grm1->nIndividuals != grm2->nIndividuals)
        )
      {
        misc.error("Error: An internal error was happened. GRMs with different individuals can not be added.", 0);
      }
    }
    else
    {
      if(
        (grm1->individualIdsRows != grm2->individualIdsRows) ||
        (grm1->individualIdsIdxRows != grm2->individualIdsIdxRows) ||
        (grm1->nIndividualsRows != grm2->nIndividualsRows)
      )
      {
        misc.error("Error: An internal error was happened. GRMs with different individuals can not be added.", 0);
      }
      
      if(
        (grm1->individualIdsCols != grm2->individualIdsCols) ||
        (grm1->individualIdsIdxCols != grm2->individualIdsIdxCols) ||
        (grm1->nIndividualsCols != grm2->nIndividualsCols)
      )
      {
        misc.error("Error: An internal error was happened. GRMs with different individuals can not be added.", 0);
      }
    }
  }
  
  if(grm2 != this)
  {
    copyParameters(grm1);
    
    if(this->kernel != NULL)
    {
      delete this->kernel;
      this->kernel = NULL;
    }
    if(this->N != NULL)
    {
      delete this->N;
      this->N = NULL;
    }
    
    this->kernel = new Matrix(grm1->kernel);
    this->N = new Matrix(grm1->N);
    
    bool thisNormalized = this->normalized;
    bool grm2Normalized = grm2->normalized;
    this->denormalize();
    grm2->denormalize();
    
    this->kernel->add(grm2->kernel, scalingFactor1, scalingFactor2);
    this->N->add(grm2->N, scalingFactor1, scalingFactor2);
    
    if(thisNormalized)
    {
      this->normalize();
    }
    if(grm2Normalized)
    {
      grm2->normalize();
    }
  }
  else
  {
    bool grm1Normalized = grm1->normalized;
    bool grm2Normalized = grm2->normalized;
    grm1->denormalize();
    grm2->denormalize();
    
    grm2->kernel->add(grm1->kernel, scalingFactor2, scalingFactor1);
    grm2->N->add(grm1->N, scalingFactor2, scalingFactor1);
    
    if(grm1Normalized)
    {
      grm1->normalize();
    }
    if(grm2Normalized)
    {
      grm2->normalize();
    }
  }
  
  std::vector<std::string> resultantSNPs;
  if(scalingFactor1 == scalingFactor2)
  {
    std::vector<std::string> test = intersectionStringVectors(2, &grm1->randomVarNames, &grm2->randomVarNames);
    if( test.size() != 0 )
    {
      misc.error("Error: Invalid operation. Two GRMs cannot be added together if the SNPs used for computing them intersect.", 0);
    }
    resultantSNPs = grm2->randomVarNames;
    resultantSNPs.insert(resultantSNPs.end(), grm1->randomVarNames.begin(), grm1->randomVarNames.end());
  }
  else
  {
    if(scalingFactor1 == -1. && scalingFactor2 == 1.)
    {
      std::vector<std::string> test = orderVectorAsTemplate(grm1->randomVarNames, grm2->randomVarNames);
      if(test != grm1->randomVarNames)
      {
        misc.error("Error: Invalid operation. Two GRMs cannot be substracted if all the SNPs of one are not contained into the other.", 0);
      }
      resultantSNPs = differenceBetweenTwoVectors(grm2->randomVarNames, grm1->randomVarNames);
    }
    else if(scalingFactor1 == 1. && scalingFactor2 == -1.)
    {
      std::vector<std::string> test = orderVectorAsTemplate(grm2->randomVarNames, grm1->randomVarNames);
      if(test != grm2->randomVarNames)
      {
        misc.error("Error: Invalid operation. Two GRMs cannot be substracted if all the SNPs of one are not contained into the other.", 0);
      }
      resultantSNPs = differenceBetweenTwoVectors(grm1->randomVarNames, grm2->randomVarNames);
    }
    else
    {
      misc.error("Error: An internal error was happened. Unexpected scaling factors.", 0);
    }
  }
  this->randomVarNames = resultantSNPs;
  
  this->name = grm1->name + " " + grm2->name;
}

std::vector<std::string> Kernel::minimumIndividualsRemovePairs(std::vector<int> globalIdxs1, std::vector<int> globalIdxs2)
{
  int * individualFrequency;
  std::set<int> idxsToDelete;
  std::vector<std::string> individualIdsToKeep;
  
  if(communicator->mpiRoot)
  {
    //Search the number of times each individuals (i.e. each index) appears.
    individualFrequency = new int[this->nIndividuals];
    for(int i = 0; i < this->nIndividuals; i++)
    {
      individualFrequency[i] = 0;
    }
    for(int i = 0; i < globalIdxs1.size(); i++)
    {
      if(globalIdxs1[i] != globalIdxs2[i])
      {
        individualFrequency[ globalIdxs1[i] ]++;
        individualFrequency[ globalIdxs2[i] ]++;
      }
    }
    //Which indices should be deleted?
    idxsToDelete.clear();
    for(int i = 0; i < globalIdxs1.size(); i++)
    {
      if(globalIdxs1[i] == globalIdxs2[i])
      {
        continue;
      }
      //From the combination of two indices (referring each one to one individual) that define an element in the GRM greather
      //than the cutoff, delete only the index that appeared more times.
      if(individualFrequency[ globalIdxs1[i] ] < individualFrequency[ globalIdxs2[i] ])
      {
        //If the other individual is already to eleiminate, this don't need to be eliminated.
        //Otherwise, when there are individuals with the same number of hits and share a hit, then both will be eliminated.
        if( idxsToDelete.find(globalIdxs1[i]) == idxsToDelete.end() )
        {
          idxsToDelete.insert(globalIdxs2[i]);
        }
      }
      else
      {
        //If the other individual is already to eleiminate, this don't need to be eliminated.
        //Otherwise, when there are individuals with the same number of hits and share a hit, then both will be eliminated.
        if( idxsToDelete.find(globalIdxs2[i]) == idxsToDelete.end() )
        {
          idxsToDelete.insert(globalIdxs1[i]);
        }
      }
    }
    //Which individuals must be kept?
    individualIdsToKeep.clear();
    for(int i = 0; i < this->nIndividuals; i++)
    {
      if(idxsToDelete.find(i) == idxsToDelete.end())
      {
        individualIdsToKeep.push_back(this->individualIds[i]);
      }
    }
   
    delete [] individualFrequency;
  }
  
  return individualIdsToKeep;
}

std::vector<std::string> Kernel::searchNoHighRelatedIndividuals(double cutoff)
{
  if(this->diagonalized == true)
  {
    misc.error("Error: An internal error was happened. Diagonalized grm cannot be pruned.", 0);
  }
  
  if(this->asymmetric == true)
  {
    misc.error("Error: An internal error was happened. This operation can not be performed on an asymmetric GRM", 0);
  }

  //Get the global indices (r, c) for all elements in the GRM > cutoff
  std::vector<int> localIdxs1;
  std::vector<int> localIdxs2;
  std::map<int, int> idxCount;
  this->kernel->getGlobalIndexElementsGreatherThan(cutoff, localIdxs2, localIdxs1);
  
  int * globalIdxs1;
  int * globalIdxs2;
  int nGlobalIdxs1;
  int nGlobalIdxs2;
  globalIdxs1 = communicator->asymmetricGather(&(localIdxs1[0]), localIdxs1.size(), &nGlobalIdxs1);
  globalIdxs2 = communicator->asymmetricGather(&(localIdxs2[0]), localIdxs2.size(), &nGlobalIdxs2);
  if(nGlobalIdxs1 != nGlobalIdxs2)
  {
    misc.error("Error: An internal error was happened when searching for individuals to prune in the GRM.", 0);
  }
  
  std::vector<std::string> individualIdsToKeep;
  if(communicator->mpiRoot == true)
  {
    std::vector<int> temp1(globalIdxs1, globalIdxs1 + nGlobalIdxs1);
    std::vector<int> temp2(globalIdxs2, globalIdxs2 + nGlobalIdxs2);
    individualIdsToKeep = minimumIndividualsRemovePairs(temp1, temp2);
  }

  if(communicator->mpiRoot)
  {
    delete [] globalIdxs1;
    delete [] globalIdxs2;
  }
  
  return individualIdsToKeep;
}

void Kernel::pruneKernel(double cutoff)
{
  if(this->type == kernelGRM)
  {
    misc.setGetElapsedTime("PruneGRM");
    misc.message << "Pruning GRM..." << std::endl;

    int previousNIndividuals = this->individualIds.size();
    std::vector<std::string> individualIdsToKeep = this->searchNoHighRelatedIndividuals(cutoff);
    this->filterIndividuals( individualIdsToKeep );
    
    misc.message << "GRM pruned after " << misc.setGetElapsedTime("ComputeGRM", true) << ". " <<  previousNIndividuals - this->individualIds.size() << " individuals filtered using " << cutoff << " as a threshold." << std::endl;
  }
  else
  {
    misc.error("Error: Pruning non GRM kernels still not supported.", 0);
  }
}

bool Kernel::sanitizeKernel()
{
  if(this->type == kernelGRM)
  {
    if(options.pruneGRM)
    {
      if( this->diagonalized == false )
      {
        pruneKernel(options.grmCutoff);
      }
      else
      {
        misc.error("Error: A diagonalized GRM cannot be prunned.", 0);
      }
    }
    if(this->diagonalized == false)
    {
      if(this->N == NULL)
      {
        misc.error("Error: An internal error was happened. Kernel cannot be checked without the normalization matrix.", 0);
      }
      std::vector<int> idxInds1;
      std::vector<int> idxInds2;
      this->N->getGlobalIndexInsideRange(-0.1, options.minimumNumberOverlapingSNPs - 0.1, idxInds1, idxInds2);
      std::vector<std::string> individualIdsToKeep = minimumIndividualsRemovePairs(idxInds1, idxInds2);
      double proportion = double(individualIdsToKeep.size())/double(this->nIndividuals);
      if( misc.gt( proportion < options.maximumProportionOfElimitaedIndividuals ) )
      {
        misc.message << "More than 10% of the individuals have to be removed for not having an overlaping of at least " << int(options.minimumNumberOverlapingSNPs) << " SNPs on all pairs. Skipping this GRM." << std::endl;
        return false;
      }
      if( misc.gt(individualIdsToKeep.size() == 1) ) //If any individual pass the test, at least one individual will remain since we are avoiding the diagonal. This solves the problem.
      {
        misc.message << "Sorry, not enough individuals have an overlaping of at least " << int(options.minimumNumberOverlapingSNPs) << " SNPs. Skipping this GRM." << std::endl;
        return false;
      }
      if( misc.gt(individualIdsToKeep.size() != this->nIndividuals) )
      {
        misc.message << this->nIndividuals - individualIdsToKeep.size() << " individuals will be removed since they do not have at least " << options.minimumNumberOverlapingSNPs << " overlapping non-missing SNPs." << std::endl;
        filterIndividuals( individualIdsToKeep );
      }
    }
  }
  return true;
}

void Kernel::keepWithRelatednessOutside(double lowerThreshold, double upperThreshold)
{
  if(this->diagonalized == true)
  {
    misc.error("Error: The GRM cannot be filtered when it is diagonal.", 0);
  }
  
  misc.message << "Filtering all individuals that have not at least one relationship outside the range (" << lowerThreshold << ", " << upperThreshold << ")." << std::endl;
  
  std::vector<int> idxInds1;
  std::vector<int> idxInds2;
  this->kernel->getGlobalIndexOutsideRange(lowerThreshold, upperThreshold, idxInds1, idxInds2);
  
  std::set<int> idxsToKeep;
  for(int i = 0; i<idxInds1.size(); i++)
  {
    if(idxInds1[i] != idxInds2[i])
    {
      idxsToKeep.insert(idxInds1[i]);
      idxsToKeep.insert(idxInds2[i]);
    }
  }
  std::vector<std::string> individualIdsToKeep;
  for(std::set<int>::iterator it = idxsToKeep.begin(); it != idxsToKeep.end(); ++it)
  {
    individualIdsToKeep.push_back(this->individualIds[*it]);
  }
  
  misc.message << "Filtering individuals " << this->nIndividuals - int(individualIdsToKeep.size()) << " from " << this->nIndividuals << "." << std::endl;
  
  filterIndividuals( individualIdsToKeep );
}

void Kernel::randomSubSample(double fraction, int minimum)
{
  if(this->diagonalized == true)
  {
    misc.error("Error: An internal error was happened. Diagonalized grm cannot be randomly subsampled.", 0);
  }
  
  misc.setGetElapsedTime("GRMRandomSubsampling");
  
  std::vector<std::string> randomIndividualIds;
  if(communicator->mpiRoot)
  {
    int nElements = (int)floor(fraction*double(this->nIndividuals));
    if(nElements < minimum)
    {
      if(this->nIndividuals <= nElements)
      {
        randomIndividualIds = this->individualIds;
      }
      else
      {
        randomIndividualIds = getRandomSample(this->individualIds, minimum);
      }
    }
    else
    {
      randomIndividualIds = getRandomSample(this->individualIds, nElements);
    }
  }
  filterIndividuals(randomIndividualIds);
  
  misc.message << "Created a GRM random subsample with " << this->nIndividuals << " after " << misc.setGetElapsedTime("GRMRandomSubsampling", true) << std::endl;
}

void Kernel::diagonalizeKernel()
{
  if(this->asymmetric == true)
  {
    misc.error("Error: An internal error was happened. Kernel diagonalization cannot be performed on an asymmetric GRM", 0);
  }
  if(this->diagonalized == true)
  {
    return;
  }
  
  misc.message << "Diagonalizing Kernel..." << std::endl;
  misc.setGetElapsedTime("KernelDiagonalization");
  
  normalize();
  
  this->eigenValues = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->eigenVectors = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  
  this->kernel->eigenDecomposition(this->eigenValues, this->eigenVectors);

  this->diagonalized = true;
  
  if(this->kernel != NULL)
  {
    delete this->kernel;
    this->kernel = NULL;
  }
  if(this->N != NULL)
  {
    delete this->N;
    this->N = NULL;
  }
  
  misc.message << "Kernel diagonalized after " << misc.setGetElapsedTime("KernelDiagonalization", true) << std::endl;
}

void Kernel::recoverKernelFromEigenDecomposition()
{
  if(this->diagonalized == false)
  {
    return;
  }
  
  misc.message << "Restoring Kernel to non-diagonalized form..." << std::endl;
  misc.setGetElapsedTime("KernelUnDiagonalization");
  
  Matrix *temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  temp->multiply(this->eigenVectors, 'N', this->eigenValues, 'N');
  this->kernel = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->kernel->multiply(temp, 'N',  this->eigenVectors, 'T');
  this->kernel->symmetric = true;

  this->diagonalized = false;
  
  if(this->eigenVectors != NULL)
  {
    delete this->eigenVectors;
    this->eigenVectors = NULL;
  }
  if(this->eigenValues != NULL)
  {
    delete this->eigenValues;
    this->eigenValues = NULL;
  }
  delete temp;
  
  misc.message << "Kernel restored after " << misc.setGetElapsedTime("KernelUnDiagonalization", true) << std::endl;
}

void Kernel::printGRM(bool showSNPs, int padding)
{
  double *g;
  
  if(this->diagonalized == false)
  {
    if (communicator->mpiRoot) {
      g = new double [this->kernel->nGlobRows*this->kernel->nGlobCols];
    }
    
    if(!this->asymmetric)
    {
      this->kernel->symmetrizeTriangularMatrix();
      if(this->N != NULL)
      {
        this->N->symmetrizeTriangularMatrix();
      }
      
      this->kernel->gatherMatrix(g);
      if (communicator->mpiRoot) {
        misc.message << "GRM ( " << (this->normalized?"normalized":"denormalized") << " ):" << std::endl;
        misc.message << this->nIndividuals << " " << this->individuals.size() << " " << this->individualIds.size() << " " << this->individualIdsIdx.size() << std::endl;
        for (int r = 0; r < this->kernel->nGlobRows; ++r)
        {
          misc.message << this->individuals[r].familyID << " " << this->individuals[r].individualID << " " << this->individualIds[r] << " " << this->individualIdsIdx[this->individualIds[r]] << ": ";
          for (int c = 0; c < this->kernel->nGlobCols; ++c)
          {
            misc.message << std::setw(padding) << *(g + this->kernel->nGlobRows*c + r) << " ";
          }
          misc.message << "\n";
        }
        misc.message << std::endl;
      }
      
      if (communicator->mpiRoot) {
        misc.message << "Cutoff pattern ( " << options.grmCutoff << " ):" << std::endl;
        for (int r = 0; r < this->kernel->nGlobRows; ++r)
        {
          misc.message << this->individualIds[r] << ": ";
          for (int c = 0; c < this->kernel->nGlobCols; ++c)
          {
            if( (*(g + this->kernel->nGlobRows*c + r) > options.grmCutoff) && (c!=r) )
            {
              misc.message << "_ ";
            }
            else
            {
              misc.message << ". ";
            }
          }
          misc.message << "\n";
        }
        misc.message << std::endl;
      }
      
      if(this->N != NULL)
      {
        this->N->gatherMatrix(g);
        if (communicator->mpiRoot) {
          misc.message << "N GRM:\n";
          for (int r = 0; r < this->kernel->nGlobRows; ++r)
          {
            misc.message << this->individuals[r].familyID << " " << this->individuals[r].individualID << " " << this->individualIds[r] << " " << this->individualIdsIdx[this->individualIds[r]] << ": ";
            for (int c = 0; c < this->kernel->nGlobCols; ++c)
            {
              misc.message << *(g + this->kernel->nGlobRows*c + r) << " ";
            }
            misc.message << "\n";
          }
          misc.message << std::endl;
        }
      }
    }
    else
    {
      this->kernel->gatherMatrix(g);
      if (communicator->mpiRoot) {
        misc.message << "GRM ( " << (this->normalized?"normalized":"denormalized") << " ):" << std::endl;
        misc.message << this->nIndividualsRows << " " << this->individualsRows.size() << " " << this->individualIdsRows.size() << " " << this->individualIdsIdxRows.size() << std::endl;
        misc.message << this->nIndividualsCols << " " << this->individualsCols.size() << " " << this->individualIdsCols.size() << " " << this->individualIdsIdxCols.size() << std::endl;
        misc.message << this->nIndividuals << " " << this->individuals.size() << " " << this->individualIds.size() << " " << this->individualIdsIdx.size() << std::endl;
        for (int c = 0; c < this->kernel->nGlobCols; ++c)
        {
          misc.message << this->individualsCols[c].familyID << ":" << this->individualsCols[c].individualID << ":" << this->individualIdsCols[c] << ":" << this->individualIdsIdxCols[this->individualIdsCols[c]] << " ";
        }
        misc.message << std::endl;
        for (int r = 0; r < this->kernel->nGlobRows; ++r)
        {
          misc.message << this->individualsRows[r].familyID << " " << this->individualsRows[r].individualID << " " << this->individualIdsRows[r] << " " << this->individualIdsIdxRows[this->individualIdsRows[r]] << ": ";
          for (int c = 0; c < this->kernel->nGlobCols; ++c)
          {
            misc.message << std::setw(11) << *(g + this->kernel->nGlobRows*c + r) << " ";
          }
          misc.message << "\n";
        }
        misc.message << std::endl;
      }
      
      if(this->N != NULL)
      {
        this->N->gatherMatrix(g);
        if (communicator->mpiRoot) {
          misc.message << "N GRM:\n";
          for (int c = 0; c < this->N->nGlobCols; ++c)
          {
            misc.message << this->individualsCols[c].familyID << ":" << this->individualsCols[c].individualID << ":" << this->individualIdsCols[c] << ":" << this->individualIdsIdxCols[this->individualIdsCols[c]] << " ";
          }
          misc.message << std::endl;
          for (int r = 0; r < this->N->nGlobRows; ++r)
          {
            misc.message << this->individualsRows[r].familyID << " " << this->individualsRows[r].individualID << " " << this->individualIdsRows[r] << " " << this->individualIdsIdxRows[this->individualIdsRows[r]] << ": ";
            for (int c = 0; c < this->N->nGlobCols; ++c)
            {
              misc.message << *(g + this->N->nGlobRows*c + r) << " ";
            }
            misc.message << "\n";
          }
          misc.message << std::endl;
        }
      }
    }
    
    if (communicator->mpiRoot && showSNPs == true)
    {
      for(int i = 0; i<this->randomVarNames.size(); i++)
      {
        std::cout << this->randomVarNames[i] << std::endl;
      }
    }
    
    if (communicator->mpiRoot) {
      delete [] g;
    }
  }
  else
  {
    this->eigenValues->showGlobal("EigenValues");
    this->eigenVectors->showGlobal("EigenVectors");
  }
  
  misc.message << this->kernel << " " << this->N << " " << this->eigenValues << " " << this->eigenVectors << std::endl;
}
