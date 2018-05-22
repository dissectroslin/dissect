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
#include "blockmatrix.h"
#include "communicator.h"
#include "misc.h"
#include "global.h"
#include "reml.h"
#include "kernel.h"
#include "options.h"
#include "covariancematrix.h"

#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <iomanip>
#include <algorithm>

CovarianceMatrix::CovarianceMatrix(int dimension, bool useDiagonalDistribution, int tBlockDimension)
{
  this->dimension = dimension;
  
  this->nCovarianceSubMatrices = 0;
  
  this->blockDimension = tBlockDimension;
  if( tBlockDimension >= 1)
  {
    this->byBlocks = true;
  }
  
  this->allCovarianceSubMatricesDiagonal = useDiagonalDistribution;
  
  DistributionType baseMatrixDistribution;
  if( useDiagonalDistribution == true && tBlockDimension == 1 )
  {
    baseMatrixDistribution = diagonalDistribution;
  }
  else
  {
    baseMatrixDistribution = MATRIX_DEFAULT_DISTRIBUTION;
  }
  
  this->m = new Matrix(baseMatrixDistribution, this->dimension, this->dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->mStatus = undefined;
  this->mDerivate = new Matrix(baseMatrixDistribution, this->dimension, this->dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  this->mInBlocks = BlockMatrix();
}

CovarianceMatrix::~CovarianceMatrix()
{
  if(this->m != NULL)
  {
    delete this->m;
    this->m = NULL;
  }
  
  for(int i = 0; i<this->nCovarianceSubMatrices; i++)
  {
    if(this->covarianceSubMatrices[i] != NULL)
    {
      delete this->covarianceSubMatrices[i];
      this->covarianceSubMatrices[i] = NULL;
    }
  }
  
  if(this->mDerivate != NULL)
  {
    delete this->mDerivate;
    this->mDerivate = NULL;
  }
}

void CovarianceMatrix::insertCovarianceMatrix(std::string name, Matrix * m)
{
  if(name == "")
  {
    misc.error("Error: covariance matrix can not be created with name: '" + name + "'.", 0);
  }
  for(int i = 0; i<this->nCovarianceSubMatrices; i++)
  {
    if(this->covarianceSubMatrixNames[i] == name)
    {
      misc.error("Error: There are two covariance matrices with the same name: '" + name + "'", 0);
    }
  }
  if( (this->allCovarianceSubMatricesDiagonal == true && m->distribution != diagonalDistribution) || 
      (this->allCovarianceSubMatricesDiagonal == false && m->distribution == diagonalDistribution)  )
  {
    misc.error("Error: An internal error was happened. Diagonal matrices cannot be mixed with non diagonal matrices in the covariance matrix.", 0);
  }
  
  if(m->symmetric == true)
  {
    m->symmetrizeTriangularMatrix();
  }
  this->covarianceSubMatrixNames.push_back(name);
  this->covarianceSubMatrices.push_back(m);
  this->nCovarianceSubMatrices = this->covarianceSubMatrices.size();
}

void CovarianceMatrix::insertCovarianceMatrix(std::string name, Kernel * kernel)
{
  insertCovarianceMatrix(name, kernel->getNormalizedKernel());
  
  kernel->kernel = NULL;
  kernel->eigenValues = NULL;
  delete kernel;
}

void CovarianceMatrix::insertVarianceGroup(std::string name, double groupExpectedMagnitude)
{
  if(name == "")
  {
    misc.error("Error: this group can not be created with name: '" + name + "'.", 0);
  }
  if(this->varianceGroupExpectedMagnitudes.count(name) != 0)
  {
    misc.error("Error: There are two groups with the same name: '" + name + "'", 0);
  }
  this->varianceGroupExpectedMagnitudes[name] = groupExpectedMagnitude;
}

void CovarianceMatrix::insertVariance(std::string name, std::string groupName, ParameterAttributes::VarianceType type, ParameterAttributes::VarianceTypeEffect typeEffect, double initialValue, std::set<std::string> constrainedDependingOnProductOf)
{
  if(name == "")
  {
    misc.error("Error: Invalid variance name.", 0);
  }
  if(type != ParameterAttributes::variance && type != ParameterAttributes::covariance && type != ParameterAttributes::correlation && type != ParameterAttributes::parameter)
  {
    misc.error("Error: Invalid variance type.", 0);
  }
  for(int i = 0; i<this->variances.size(); i++)
  {
    if(this->variances[i].name == name)
    {
      misc.error("Error: There are two variances with the same name: '" + name + "'", 0);
    }
  }
  if(this->varianceGroupExpectedMagnitudes.count(groupName) != 1)
  {
    misc.error("Error: When inserting a new variance, there are not any group with the name: '" + groupName + "'", 0);
  }
  
  Variance newVariance;
  newVariance.name = name;
  newVariance.group = groupName;
  newVariance.type = type;
  newVariance.typeEffect = typeEffect;
  newVariance.fixed = false;
  newVariance.variance = initialValue;
  newVariance.initialVariance = initialValue;
  newVariance.inElements = std::set<int>();
  newVariance.constrainedDependingOnProductOfName = constrainedDependingOnProductOf;

  this->variances.push_back(newVariance);
  
}

void CovarianceMatrix::insertElement(std::string subCovarianceMatrixId, std::string elementName, KernelType type, std::string insertedCovarianceMatrixName, double constantFactor, std::pair<int, int> blockPosition, subMatrix outcomeSubMatrix, subMatrix sourceSubMatrix)
{
  bool flag;
  Element newElement;
  
  if(getElement(elementName) >= 0)
  {
    misc.error("Error: Invalid element name. The element names must be unique.", 0);
  }
  
  newElement.subCovarianceId = subCovarianceMatrixId;
  newElement.name = elementName;
  newElement.type = type;
  
  flag = false;
  for(int i = 0; i<this->nCovarianceSubMatrices; i++)
  {
    if(this->covarianceSubMatrixNames[i] == insertedCovarianceMatrixName)
    {
      newElement.m = this->covarianceSubMatrices[i];
      flag = true;
    }
  }
  if(!flag)
  {
    misc.error("Error: Invalid covariance matrix name when constructing an Element: '" + insertedCovarianceMatrixName + "'", 0);
  }
  
  newElement.iVariances.clear();
  newElement.nameVariances.clear();
  newElement.transformationTypes.clear();
  newElement.parametersPosition.clear();
  
  newElement.constantFactor = constantFactor;

  if( this->byBlocks == true )
  {
    if( blockPosition.first < 0 || blockPosition.second < 0 || blockPosition.first >= this->blockDimension || blockPosition.second >= this->blockDimension )
    {
      misc.error("Error: An internal error was happened. Inserting a matrix with no block position properly specified in a blocked covariance matrix.", 0);
    }
  }
  
  newElement.blockPosition = blockPosition;
  newElement.outcomeSubMatrix = outcomeSubMatrix;
  newElement.sourceSubMatrix = sourceSubMatrix;
  
  this->elements.push_back(newElement);
}

void CovarianceMatrix::appendVarianceToElement(std::string elementName, std::string varianceName, ParameterAttributes::TransformationType transformationType, ParameterAttributes::ParameterPosition parameterPosition)
{
  int idxElement = getElement(elementName);
  if(idxElement < 0)
  {
    misc.error("Error: Variance cannot be appended to a non existing element.", 0);
  }
  for(int j = 0; j<this->variances.size(); j++)
  {
    if(this->variances[j].name == varianceName)
    {
      this->elements[idxElement].iVariances.push_back(j);
      this->elements[idxElement].nameVariances.push_back(varianceName);
      this->elements[idxElement].transformationTypes.push_back( transformationType );
      this->elements[idxElement].parametersPosition.push_back( parameterPosition );
      this->variances[j].inElements.insert(idxElement);
      return;
    }
  }
  misc.error("Error: Invalid variance name when adding variances to an Element: " + varianceName, 0);
}

bool CovarianceMatrix::checkMoreThanOneVariance()
{
  bool moreThanOneVariance = false;
  for(int i = 0; i < elements.size(); i++)
  {
    if(elements[i].iVariances.size() > 1)
    {
      moreThanOneVariance = true;
    }
  }
  
  communicator->broadcast(&moreThanOneVariance);
  
  return moreThanOneVariance;
}

int CovarianceMatrix::getElement(std::string elementName)
{
  for(int i = 0; i<this->elements.size(); i++)
  {
    if(this->elements[i].name == elementName)
    {
      return i;
    }
  }
  return -1;
}

int CovarianceMatrix::getVariance(std::string varianceName)
{
  for(int i = 0; i<this->variances.size(); i++)
  {
    if(this->variances[i].name == varianceName)
    {
      return i;
    }
  }
  misc.error("Error: An internal error was happened. Malformed Covariance Matrix. The variance name do no exist.", 0);
  return -1;
}

bool CovarianceMatrix::checkVarianceExists(std::string varianceName)
{
  for(int i = 0; i<this->variances.size(); i++)
  {
    if(this->variances[i].name == varianceName)
    {
      return true;
    }
  }
  return false;
}

void CovarianceMatrix::reinitializeVariances()
{
  for(int i = 0; i<this->variances.size(); i++)
  {
    this->variances[i].variance = this->variances[i].initialVariance;
  }
  broadcastVariances();
}

void CovarianceMatrix::setInitialVariancesWithCurrentValues()
{
  broadcastVariances();
  for(int i = 0; i<this->variances.size(); i++)
  {
    this->variances[i].initialVariance = this->variances[i].variance;
  }
}

void CovarianceMatrix::backupFullModel()
{
  if( this->variancesFullModel.size() != 0 ||  this->varianceGroupExpectedMagnitudesFullModel.size() != 0 ||  this->elementsFullModel.size() != 0 )
  {
    misc.error("Error: An internal error was happened. Full model cannot be stored. There is already a model stored.", 0);
  }
  
  this->variancesFullModel = this->variances;
  this->varianceGroupExpectedMagnitudesFullModel = this->varianceGroupExpectedMagnitudes;
  this->elementsFullModel = this->elements;
}

void CovarianceMatrix::restoreBackupFullModel()
{
  if( this->variancesFullModel.size() == 0 ||  this->varianceGroupExpectedMagnitudesFullModel.size() == 0 ||  this->elementsFullModel.size() == 0 )
  {
    misc.error("Error: An internal error was happened. Full model cannot be recovered. There stored model is empty.", 0);
  }
  
  this->variances = this->variancesFullModel;
  this->varianceGroupExpectedMagnitudes = this->varianceGroupExpectedMagnitudesFullModel;
  this->elements = this->elementsFullModel;
  
  this->variancesFullModel.clear();
  this->varianceGroupExpectedMagnitudesFullModel.clear();
  this->elementsFullModel.clear();
}

void CovarianceMatrix::clearElementVariances(std::string elementName)
{
  int idxElement = getElement(elementName);
  if(idxElement < 0)
  {
    misc.error("Error: An internal error was happened. Inexistent element when trying to clear element.", 0);
  }
  for(int i = 0; i<this->elements[idxElement].iVariances.size(); i++)
  {
    int idxVariance = this->elements[idxElement].iVariances[i];
    int removed = this->variances[idxVariance].inElements.erase(idxElement);
    if(removed != 1)
    {
      misc.error("Error: An internal error was happened. Malformed Covariance Matrix.", 0);
    }
  }
  this->elements[idxElement].iVariances.clear();
  this->elements[idxElement].nameVariances.clear();
  this->elements[idxElement].transformationTypes.clear();
  this->elements[idxElement].parametersPosition.clear();
}

void CovarianceMatrix::changeElementConstantFactor(std::string elementName, double constantFactor)
{
  int idxElement = getElement(elementName);
  if(idxElement < 0)
  {
    misc.error("Error: An internal error was happened. Inexistent element when trying to change the constant factor.", 0);
  }
  this->elements[idxElement].constantFactor = constantFactor;
}

void CovarianceMatrix::deleteVariance(std::string name)
{
  int idxVariance = -1;
  for(int i = 0; i<this->variances.size(); i++)
  {
    if(this->variances[i].name == name)
    {
      idxVariance = i;
    }
  }
  if(idxVariance < 0)
  {
    misc.error("Error: An internal error was happened. Invalid variance name when deleting the variance: " + name, 0);
  }
  if(this->variances[idxVariance].inElements.size() != 0)
  {
    misc.error("Error: An internal error was happened. Trying to remove a variance before removing it from elements.", 0);
  }
  
  std::string group = this->variances[idxVariance].group;
  double varianceValue = this->variances[idxVariance].variance;
  ParameterAttributes::VarianceType varianceType = this->variances[idxVariance].type;
  this->variances.erase(this->variances.begin() + idxVariance);
  
  bool flag = false;
  int nVariancesInGroupToDistribute = 0;
  double totalVarianceGroup = 0.;
  for(int i = 0; i<this->variances.size(); i++)
  {
    if(this->variances[i].group == group)
    {
      flag = true;
      if( varianceType == ParameterAttributes::variance && this->variances[i].type == ParameterAttributes::variance )
      {
        nVariancesInGroupToDistribute++;
        totalVarianceGroup += this->variances[i].variance;
      }
    }
  }
  if(flag == false)
  {
    int erased = varianceGroupExpectedMagnitudes.erase(group);
    if(erased != 1)
    {
      misc.error("Error: An internal error was happened. Malformed Covariance Matrix.", 0);
    }
  }
  else
  {
    //If the deleted variance is of type variance and there are more variances in the group of type variance,
    //then distribute the effect of the deleted variance between them
    if(nVariancesInGroupToDistribute != 0)
    {
      for(int i = 0; i<this->variances.size(); i++)
      {
        if(this->variances[i].group == group && this->variances[i].type == ParameterAttributes::variance )
        {
          double factor = this->variances[i].variance/totalVarianceGroup;
          this->variances[i].variance += factor*varianceValue/double(nVariancesInGroupToDistribute);
        }
      }
    }
  }
  
  for(int i = 0; i < this->elements.size(); i++)
  {
    for(int j = 0; j < this->elements[i].nameVariances.size(); j++)
    {
      if( (this->elements[i].nameVariances[j] == name) || (this->elements[i].nameVariances.size() != this->elements[i].iVariances.size()) )
      {
        misc.error("Error: An internal error was happened. Malformed Covariance Matrix.", 0);
      }
      this->elements[i].iVariances[j] = getVariance(this->elements[i].nameVariances[j]);
    }
  }
  
  updateVariancesDependenceIndexs();
}

void CovarianceMatrix::deleteElementsWithSubCovarianceId(std::string subCovarianceMatrixId)
{
  std::set<std::string> varianceToDelete;
  std::vector<int> idxElementsToDelete;
  
  //Search elements to delete and recopile variances in these elements.
  for(int i = 0; i < this->elements.size(); i++)
  {
    if(this->elements[i].subCovarianceId == subCovarianceMatrixId)
    {
      for(int j = 0; j < this->elements[i].nameVariances.size(); j++)
      {
        varianceToDelete.insert(elements[i].nameVariances[j]);
      }
      idxElementsToDelete.push_back(i);
    }
  }
  
  if(idxElementsToDelete.size() < 1)
  {
    misc.error("Error: An internal error was happened. There are not any element to delete with subCovarianceId  equal to " + subCovarianceMatrixId + ".", 0);
  }
  
  //Remove elements.
  for(int i = idxElementsToDelete.size() - 1; i >= 0; i--)
  {
    this->elements.erase( this->elements.begin() + idxElementsToDelete[i] );
  }
  
  //Check if the variances in removed elements are in other elements. In this case these cariances have to be kept
  for(int i = 0; i < this->elements.size(); i++)
  {
    for(int j = 0; j < this->elements[i].nameVariances.size(); j++)
    {
      if(varianceToDelete.find( elements[i].nameVariances[j] ) != varianceToDelete.end())
      {
        varianceToDelete.erase( elements[i].nameVariances[j] );
      }
    }
  }
  
  //Remove variances
  for (std::set<std::string>::iterator it = varianceToDelete.begin(); it != varianceToDelete.end(); ++it) //To avoid errors when deleting variances due a dependences no longer satisfied.
  {
    int temp = getVariance(*it);
    this->variances[temp].constrainedDependingOnProductOfName.clear();
    this->variances[temp].constrainedDependingOnProductOfi.clear();
  }
  for (std::set<std::string>::iterator it = varianceToDelete.begin(); it != varianceToDelete.end(); ++it)
  {
    this->variances[getVariance(*it)].inElements.clear();
    deleteVariance(*it);
  }
  
  //Recode inElements for remaining variances
  for(int i = 0; i<this->variances.size(); i++)
  {
    this->variances[i].inElements.clear();
  }
  for(int i = 0; i < this->elements.size(); i++)
  {
    for(int j = 0; j < this->elements[i].iVariances.size(); j++)
    {
      int vidx = this->elements[i].iVariances[j];
      this->variances[vidx].inElements.insert(i);
    }
  }
}

void CovarianceMatrix::updateVariancesDependenceIndexs()
{
  for(int i = 0; i<this->variances.size(); i++)
  {
    this->variances[i].constrainedDependingOnProductOfi.clear();
    for(std::set<std::string>::iterator it = this->variances[i].constrainedDependingOnProductOfName.begin(); it != this->variances[i].constrainedDependingOnProductOfName.end(); ++it)
    {
      int temp = getVariance(*it);
      this->variances[i].constrainedDependingOnProductOfi.insert(temp);
    }
  }
}

void CovarianceMatrix::computeCovariance()
{
  this->m->fillWithConstant(0.);
  this->m->symmetric = true;
  this->m->uplo = 'B';
  
  for(int i=0; i<this->elements.size(); i++)
  {
    Element currentElement = this->elements[i];
    
    std::pair<double, Matrix *> mf = computeElementMatrixAndFactor(currentElement);
    double currentVariance = mf.first;
    Matrix * melement = mf.second;
    
    //Add to the covariance a subcovariance matrix multiplied by the previous product of variances.
    if(currentElement.outcomeSubMatrix.active == false || currentElement.sourceSubMatrix.active == false)
    {
      this->m->add(melement, (i==0)?0.:1., currentVariance);
    }
    else
    {
      this->m->add(melement, (i==0)?0.:1., currentVariance, currentElement.outcomeSubMatrix, currentElement.sourceSubMatrix);
    }
    
    if( currentElement.m != melement && melement != NULL )
    {
      delete melement;
    }
  }
  this->m->symmetric = true;
  this->m->uplo = 'U';
  this->mStatus = covarianceMatrix;
}

void CovarianceMatrix::computeBlockCovariance()
{
  if( this->byBlocks == false )
  {
    misc.error("Error: An internal error was happened. Covariance matrix cannot be computed by blocks. Blocks undefined.", 0);
  }
  
  int nrBlocks = 0;
  int ncBlocks = 0;
  
  std::vector< std::vector<Matrix*> > blocks = std::vector< std::vector<Matrix*> >(this->blockDimension, std::vector<Matrix*>(this->blockDimension));
    
  for(int i=0; i<this->elements.size(); i++)
  {
    Element currentElement = this->elements[i];
    
    //Compute the product of the variance elements
    std::pair<double, Matrix *> mf = computeElementMatrixAndFactor(currentElement);
    double currentVariance = mf.first;
    Matrix * melement = mf.second;
    
    
    int br = currentElement.blockPosition.first;
    int bc = currentElement.blockPosition.second;
    
    if( br > bc )
    {
      misc.error("Error: An internal error was happened wen computing the covariance matrix in block form. No elements expected below the diagonal.", 0);
    }
    
    //Add to the covariance a subcovariance matrix multiplied by the previous product of variances.
    if(blocks[br][bc] == NULL)
    {
      Matrix * temp = new Matrix(melement->distribution, melement->nGlobRows, melement->nGlobCols);
      temp->fillWithConstant(0.);
      blocks[br][bc] = temp;
    }
    blocks[br][bc]->add(melement, 1., currentVariance);
    
    if( currentElement.m != melement && melement != NULL )
    {
      delete melement;
    }
  }
  
  this->mInBlocks.clear();
  
  for(int br=0; br<this->blockDimension; br++)
  {
    std::vector<Matrix*> rowBlocks;
    for(int bc=0; bc<this->blockDimension; bc++)
    {
      if( bc >= br && blocks[br][bc] == NULL )
      {
        misc.error("Error: An internal error was happened wen computing the covariance matrix in block form. Empty block.", 0);
      }
      if( bc < br)
      {
        if(blocks[br][bc] != NULL)
        {
          misc.error("Error: An internal error was happened wen computing the covariance matrix in block form. Empty block.", 0);
        }
        blocks[br][bc] = new Matrix(blocks[bc][br]);
      }
      rowBlocks.push_back(blocks[br][bc]);
    }
    this->mInBlocks.addBlockRow(rowBlocks);
  }
  blocks.clear();
  
  this->mStatus = covarianceMatrixByBlocks;
}

Matrix * CovarianceMatrix::computeDerivateCovariance(int idxDerivateVariance)
{
  bool flag = false;
  
  if(variances[idxDerivateVariance].inElements.size() == 0)
  {
    misc.error("Error: An internal error was happened on method computeDerivateCovariance().", 0);
  }
  
  this->mDerivate->fillWithConstant(0.);
  this->mDerivate->symmetric = true;
  this->mDerivate->uplo = 'B';
  
  //Iterate over the elements that has the variance with index idxDerivateVariance in their components
  int i = 0;
  for (std::set<int>::iterator it=this->variances[idxDerivateVariance].inElements.begin(); it!=this->variances[idxDerivateVariance].inElements.end(); ++it)
  {
    flag = true;
    Element currentElement = this->elements[ *it ];
    
    //Compute the derivate of the current element variances product with respect the variance idxDerivateVariance.
    std::pair<double, Matrix *> mf = computeElementDerivateMatrixAndFactor(currentElement, idxDerivateVariance);
    double currentVariance = mf.first;
    Matrix * melement = mf.second;
    
    if(melement != NULL)
    {
      if(currentElement.outcomeSubMatrix.active == false || currentElement.sourceSubMatrix.active == false)
      {
        this->mDerivate->add(melement, (i==0)?0.:1., currentVariance);
      }
      else
      {
        this->mDerivate->add(melement, (i==0)?0.:1., currentVariance, currentElement.outcomeSubMatrix, currentElement.sourceSubMatrix);
      }
      i++;
    }
    
    if( currentElement.m != melement && melement != NULL )
    {
      delete melement;
    }
  }
  
  if(flag == false)
  {
    misc.error("Error: An internal error was happened on method computeDerivateCovariance().", 0);
  }
  
  this->mDerivate->symmetric = true;
  this->mDerivate->uplo = 'U';
  
  return this->mDerivate;
}

Matrix * CovarianceMatrix::computeDerivateCovariance(int idxDerivateVariance1, int idxDerivateVariance2)
{
  if(variances[idxDerivateVariance1].inElements.size() == 0 || variances[idxDerivateVariance2].inElements.size() == 0)
  {
    misc.error("Error: An internal error was happened on method computeDerivateCovariance().", 0);
  }
  
  Matrix * result = NULL;
  
  std::set<int> s1 = this->variances[idxDerivateVariance1].inElements;
  std::set<int> s2 = this->variances[idxDerivateVariance2].inElements;
  std::set<int> intersectElements;
  std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(), std::inserter(intersectElements,intersectElements.begin()));
  
  int i = 0;
  for (std::set<int>::iterator it=intersectElements.begin(); it!=intersectElements.end(); ++it)
  {
    Element currentElement = this->elements[ *it ];
    
    //Compute the derivate of the current element variances product with respect the variances idxDerivateVariance1 and idxDerivateVariance2.
    std::pair<std::pair<bool, double>, Matrix *> mf = computeElementDerivateMatrixAndFactor(currentElement, idxDerivateVariance1, idxDerivateVariance2);
    std::pair<bool, double> currentVarianceResult = mf.first;
    double currentVariance = currentVarianceResult.second;
    Matrix * melement = mf.second;
    
    
    if(currentVarianceResult.first == true)
    {
      if(result == NULL)
      {
        result = new Matrix(this->mDerivate->distribution, this->dimension, this->dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
        result->fillWithConstant(0.);
        result->symmetric = true;
        result->uplo = 'B';
      }
      
      if(currentElement.outcomeSubMatrix.active == false || currentElement.sourceSubMatrix.active == false)
      {
        result->add(melement, (i==0)?0.:1., currentVariance);
      }
      else
      {
        result->add(melement, (i==0)?0.:1., currentVariance, currentElement.outcomeSubMatrix, currentElement.sourceSubMatrix);
      }
      i++;
    }
    if( currentElement.m != melement && melement != NULL)
    {
      delete melement;
    }
  }
  
  if( result != NULL )
  {
    result->symmetric = true;
    result->uplo = 'U';
  }
  
  return result;
}

std::pair<double, Matrix *> CovarianceMatrix::computeElementMatrixAndFactor(Element element)
{
  Matrix * m = NULL;
  double factor = 1.;
  
  factor = computeElementMultiplyingVariance(element);
  
  if( element.type == kernelGRM || 
      element.type == kernelEpistaticGRM ||
      element.type == kernelFromDiscreteCovariates ||
      element.type == kernelFromMultiDiscreteCovariates ||
      element.type == kernelCovarianceMatrix ||
      element.type == kernelEnvirontmental ||
      element.type == kernelInteraction
    ) 
  {
    m = element.m;
  }
  else if( element.type == kernelSquaredExponential )
  {
    int iParam = -1;
    int nParams = 0;
    for(int i = 0; i<element.iVariances.size(); i++)
    {
      int itemp = element.iVariances[i];
      if( this->variances[itemp].type == ParameterAttributes::parameter )
      {
        nParams++;
        iParam = itemp;
      }
    }
    if(nParams != 1)
    {
      misc.error("Error: An internal error was happened. Unexpected element number of parameters parameters when computing a covariance element.", 0);
    }
    double param = this->variances[iParam].variance;
    m = new Matrix(element.m);
    m->applyExponentialOperator(-param);
  }
  else
  {
    misc.error("Error: An internal error was happened. Invalid kernel type when computing a covariance element.", 0);
  }
  
  
  std::pair<double, Matrix *> result = std::pair<double, Matrix *>(factor, m);
  return result;
}

std::pair<double, Matrix *> CovarianceMatrix::computeElementDerivateMatrixAndFactor(Element element, int idxDerivateVariance)
{
  Matrix * m = NULL;
  double factor = 1.;
  
  if(this->variances[idxDerivateVariance].fixed == true)
  {
    return std::pair<double, Matrix *>(0, NULL);
  }
  
  if( element.type == kernelGRM || 
      element.type == kernelEpistaticGRM ||
      element.type == kernelFromDiscreteCovariates ||
      element.type == kernelFromMultiDiscreteCovariates ||
      element.type == kernelCovarianceMatrix ||
      element.type == kernelEnvirontmental ||
      element.type == kernelInteraction
    )
  {
    factor = computeElementMultiplyingVarianceDerivate(element, idxDerivateVariance);
    m = element.m;
  }
  else if( element.type == kernelSquaredExponential )
  {
    int iParam = -1;
    int nParams = 0;
    for(int i = 0; i<element.iVariances.size(); i++)
    {
      int itemp = element.iVariances[i];
      if( this->variances[itemp].type == ParameterAttributes::parameter )
      {
        nParams++;
        iParam = itemp;
      }
    }
    if(nParams != 1)
    {
      misc.error("Error: An internal error was happened. Unexpected element number of parameters parameters when computing a covariance element.", 0);
    }
    double param = this->variances[iParam].variance;
    m = new Matrix(element.m);
    m->applyExponentialOperator(-param);
    if(idxDerivateVariance == iParam)
    {
      m->elementWiseMultiplication(element.m, -1.);
      factor = computeElementMultiplyingVariance(element);
    }
    else
    {
      factor = computeElementMultiplyingVarianceDerivate(element, idxDerivateVariance);
    }
  }
  else
  {
    misc.error("Error: An internal error was happened. Invalid kernel type when computing a covariance element derivate.", 0);
  }
  
  std::pair<double, Matrix *> result = std::pair<double, Matrix *>(factor, m);
  return result;
}


std::pair<std::pair<bool, double>, Matrix *> CovarianceMatrix::computeElementDerivateMatrixAndFactor(Element element, int idxDerivateVariance1, int idxDerivateVariance2)
{
  Matrix * m = NULL;
  //double factor = 1.;
  std::pair<bool, double> factor = std::pair<bool, double>(false, 0.);
  
  if(this->variances[idxDerivateVariance1].fixed == true || this->variances[idxDerivateVariance2].fixed == true)
  {
    return std::pair<std::pair<bool, double>, Matrix *>(std::pair<bool, double>(false, 0.), NULL);
  }
  
  if( element.type == kernelGRM || 
      element.type == kernelEpistaticGRM ||
      element.type == kernelFromDiscreteCovariates ||
      element.type == kernelFromMultiDiscreteCovariates ||
      element.type == kernelCovarianceMatrix ||
      element.type == kernelEnvirontmental ||
      element.type == kernelInteraction
    )
  {
    factor = computeElementMultiplyingVarianceDerivate(element, idxDerivateVariance1, idxDerivateVariance2);
    m = element.m;
  }
  else if( element.type == kernelSquaredExponential )
  {
    int iParam = -1;
    int nParams = 0;
    for(int i = 0; i<element.iVariances.size(); i++)
    {
      int itemp = element.iVariances[i];
      if( this->variances[itemp].type == ParameterAttributes::parameter )
      {
        nParams++;
        iParam = itemp;
      }
    }
    if(nParams != 1)
    {
      misc.error("Error: An internal error was happened. Unexpected element number of parameters parameters when computing a covariance element.", 0);
    }
    double param = this->variances[iParam].variance;
    m = new Matrix(element.m);
    m->applyExponentialOperator(-param);
    
    if(idxDerivateVariance1 == iParam && idxDerivateVariance2 == iParam)
    {
      m->elementWiseMultiplication(element.m);
      m->elementWiseMultiplication(element.m);

      double temp = computeElementMultiplyingVariance(element);
      factor = std::pair<bool, double>(true, temp);
    }
    else if(idxDerivateVariance1 == iParam && idxDerivateVariance2 != iParam)
    {
      m->elementWiseMultiplication(element.m, -1.);
      double temp = computeElementMultiplyingVarianceDerivate(element, idxDerivateVariance2);
      factor = std::pair<bool, double>(true, temp);
    }
    else if(idxDerivateVariance1 != iParam && idxDerivateVariance2 == iParam)
    {
      m->elementWiseMultiplication(element.m, -1.);
      double temp = computeElementMultiplyingVarianceDerivate(element, idxDerivateVariance1);
      factor = std::pair<bool, double>(true, temp);
    }
    else
    {
      factor = computeElementMultiplyingVarianceDerivate(element, idxDerivateVariance1, idxDerivateVariance2);
    }
  }
  else
  {
    misc.error("Error: An internal error was happened. Invalid kernel type when computing a covariance element derivates.", 0);
  }
  
  std::pair<std::pair<bool, double>, Matrix *> result = std::pair<std::pair<bool, double>, Matrix *>(factor, m);

  return result;
}

double CovarianceMatrix::computeElementMultiplyingVariance(Element & element)
{
  if(element.iVariances.size() == 0)
  {
    misc.error("Error: An internal error was happened. Covariance matrix malformed. There is an element without variances.", 0);
  }

  double currentVariance = 1.;
  for(int j = 0; j<element.iVariances.size(); j++)
  {
    if(element.parametersPosition[j] != ParameterAttributes::multiplyingMatrix)
    {
      continue;
    }
    
    if(element.transformationTypes[j] == ParameterAttributes::nochange)
    {
      currentVariance *= this->variances[element.iVariances[j]].variance;
    }
    else if(element.transformationTypes[j] == ParameterAttributes::squareRoot)
    {
      currentVariance *= sqrt(this->variances[element.iVariances[j]].variance); //The sqrt if this variance is, for instance, a standard deviation.
    }
    else
    {
      misc.error("Error: An internal error was happened", 0);
    }
  }
  return currentVariance*element.constantFactor;
}

double CovarianceMatrix::computeElementMultiplyingVarianceDerivate(Element & element, int idxDerivateVariance)
{
  bool flag = false;
  double currentVariance = 1.;
  for(int j = 0; j<element.iVariances.size(); j++) //Creates the derivate of the product of variances with respect the variance idxDerivateVariance
  {
    int idx  = element.iVariances[j];
    
    if(element.parametersPosition[j] != ParameterAttributes::multiplyingMatrix)
    {
      if( idx == idxDerivateVariance )
      {
        misc.error("Error: An internal error was happened. The derivate factor of a non multiplying parameter cannot be computed inside computeElementVarianceDerivate() method.", 0);
      }
      continue;
    }
    
    if(idx != idxDerivateVariance)
    {
      if(element.transformationTypes[j] == ParameterAttributes::nochange)
      {
        currentVariance *= this->variances[idx].variance;
      }
      else if(element.transformationTypes[j] == ParameterAttributes::squareRoot)
      {
        currentVariance *= sqrt(this->variances[idx].variance);
      }
      else
      {
        misc.error("Error: An internal error was happened.", 0);
      }
    }
    else
    {
      flag = true; //Just check everything is ok.
      if(element.transformationTypes[j] == ParameterAttributes::nochange)
      {
        currentVariance *= 1.;
      }
      else if(element.transformationTypes[j] == ParameterAttributes::squareRoot)
      {
        currentVariance *= 0.5/sqrt(this->variances[idx].variance); //If the product depends not on the variance but on the standard deviation, the derivate is 1/(2*std) instead of 1.
      }
      else
      {
        misc.error("Error: An internal error was happened.", 0);
      }
    }
  }
  if(!flag)
  {
    misc.error("Error: An internal error was happened when computing a variance derivate.", 0);
  }
  return currentVariance*element.constantFactor;
}

std::pair<bool, double> CovarianceMatrix::computeElementMultiplyingVarianceDerivate(Element & element, int idxDerivateVariance1, int idxDerivateVariance2)
{
  std::pair<bool, double> result = std::pair<bool, double>(false, 0.);
  
  int flag1 = 0;
  int flag2 = 0;
  for(int j = 0; j<element.iVariances.size(); j++) //Creates the derivate of the product of variances with respect the variance idxDerivateVariance
  {
    int idx  = element.iVariances[j];
    
    if(element.parametersPosition[j] != ParameterAttributes::multiplyingMatrix)
    {
      if( idx == idxDerivateVariance1 || idx == idxDerivateVariance2 )
      {
        misc.error("Error: An internal error was happened. The derivate factor of a non multiplying parameter cannot be computed inside computeElementVarianceDerivate() method.", 0);
      }
      continue;
    }
    
    if(idx == idxDerivateVariance1)
    {
      flag1++;
    }
    if(idx == idxDerivateVariance2)
    {
      flag2++;
    }
  }
  if( flag1 > 1 || flag2 > 1)
  {
    misc.error("Error: An internal error was happened when computing the second derivative of the covariance.", 0);
  }
  if( flag1 == 0 || flag2 == 0)
  {
    return result;
  }
  
  double currentVariance = 1.;
  for(int j = 0; j<element.iVariances.size(); j++) //Creates the derivate of the product of variances with respect the variance idxDerivateVariance1, and idxDerivateVariance2
  {
    int idx  = element.iVariances[j];
    if(idx != idxDerivateVariance1 && idx != idxDerivateVariance2)
    {
      if(element.transformationTypes[j] == ParameterAttributes::nochange)
      {
        currentVariance *= this->variances[idx].variance;
      }
      else if(element.transformationTypes[j] == ParameterAttributes::squareRoot)
      {
        currentVariance *= sqrt(this->variances[idx].variance);
      }
      else
      {
        misc.error("Error: An internal error was happened.", 0);
      }
    }
    else
    {
      if( idxDerivateVariance1 != idxDerivateVariance2 )
      {
        if(element.transformationTypes[j] == ParameterAttributes::nochange)
        {
          currentVariance *= 1.;
        }
        else if(element.transformationTypes[j] == ParameterAttributes::squareRoot)
        {
          currentVariance *= 0.5/sqrt(this->variances[idx].variance); //If the product depends not on the variance but on the standard deviation, the derivate is 1/(2*std) instead of 1.
        }
        else
        {
          misc.error("Error: An internal error was happened.", 0);
        }
      }
      else
      {
        if(element.transformationTypes[j] == ParameterAttributes::nochange)
        {
          currentVariance = 0.;
          result = std::pair<bool, double>(false, 0.);
          return result;
        }
        else if(element.transformationTypes[j] == ParameterAttributes::squareRoot)
        {
          double temp = this->variances[idx].variance;
          currentVariance *= -0.25/sqrt(temp*temp*temp); //If the product depends not on the variance but on the standard deviation, the second derivate is -0.25/(std^3) instead of 1.
        }
        else
        {
          misc.error("Error: An internal error was happened.", 0);
        }
      }
    }
  }
  
  result = std::pair<bool, double>(true, currentVariance*element.constantFactor);
  
  return result;
}

bool CovarianceMatrix::invertCovariance(double * logDeterminant, bool useSinglePrecision)
{
  if( this->allCovarianceSubMatricesDiagonal == false || this->byBlocks == false || this->blockDimension <= 1)
  {
    if(this->mStatus != covarianceMatrix)
    {
      computeCovariance();
    }
    bool inverted = this->m->symmetricInvert(logDeterminant, useSinglePrecision);
    if(inverted == false)
    {
      misc.message << "Warning: The Cholesky factorization cannot be completed because the covariance matrix is not positive-definite." << std::endl;
      computeCovariance();
      inverted = this->m->invert(logDeterminant);
      if(inverted == false)
      {
        misc.message << "Error: V matrix can not be inverted. REML iterations cannot continue." << std::endl;
        return false;
      }
    }
    this->mStatus = inverseCovarianceMatrix;
    return true;
  }
  else //When all matrices are diagonal and block dimension is gretaher than one, use blockMatrix for inversion.
  {
    computeBlockCovariance();
    bool inverted = this->mInBlocks.invert(logDeterminant);
    if(inverted == false)
    {
      misc.message << "Error: V matrix can not be inverted. REML iterations cannot continue." << std::endl;
      return false;
    }
    
    this->mStatus = inverseCovarianceMatrixByBlocks;
    
    return true;
  }
}

int CovarianceMatrix::constrainVariancesM1(std::set<std::string> & constrainedVariances)
{
  constrainedVariances.clear();
  
  int nVariancesConstrained = 0;
  std::map<std::string, int> nVariancesConstrainedGroup;
  std::map<std::string, double> correctedMagnitude;
  std::map<std::string, double> nonCorrectedMagnitude;
  for(std::map<std::string, double>::iterator it = varianceGroupExpectedMagnitudes.begin(); it != varianceGroupExpectedMagnitudes.end(); ++it)
  {
    nVariancesConstrainedGroup[it->first] = 0;
    correctedMagnitude[it->first] = 0.;
    nonCorrectedMagnitude[it->first] = 0.;
  }
  std::vector<bool> constrained(this->variances.size(), false);
  
  for(int i = 0; i < this->variances.size(); i++)
  {
    std::string idGroup = this->variances[i].group;
    double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
    
    if(this->variances[i].variance < 0 && this->variances[i].type == ParameterAttributes::variance)
    {
      correctedMagnitude[idGroup] += baseMagnitude*options.varianceConstrainProportion - this->variances[i].variance;
      this->variances[i].variance = baseMagnitude*options.varianceConstrainProportion;
      constrained[i] = true;
      nVariancesConstrainedGroup[idGroup]++;
      nVariancesConstrained++;
      constrainedVariances.insert(this->variances[i].name);
    }
    else if(this->variances[i].variance > 0 && this->variances[i].type == ParameterAttributes::variance)
    {
      nonCorrectedMagnitude[idGroup] += this->variances[i].variance;
    }
  }
  
  //Use GCTA methods?
  if(options.remlGCTAMode)
  {
    for(std::map<std::string, double>::iterator it = varianceGroupExpectedMagnitudes.begin(); it != varianceGroupExpectedMagnitudes.end(); ++it)
    {
      if(nVariancesConstrainedGroup[it->first] != 0)
      {
        correctedMagnitude[it->first] /= double(nVariancesConstrainedGroup[it->first]);
      }
    }
    
    for(int i = 0; i < this->variances.size(); i++)
    {
      std::string idGroup = this->variances[i].group;
      double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
      if(constrained[i] == false && this->variances[i].type == ParameterAttributes::variance)
      {
        this->variances[i].variance -= correctedMagnitude[idGroup];
        /*this->variances[i].variance -= correctedMagnitude[idGroup]*this->variances[i].variance/nonCorrectedMagnitude[idGroup];
        if(this->variances[i].variance < 0.)
        {
          this->variances[i].variance = baseMagnitude*options.varianceConstrainProportion;
        }*/
      }
    }
  }
  
  //Constrain covariances.
  for(int i = 0; i < this->variances.size(); i++)
  {
    
    if( this->variances[i].type == ParameterAttributes::covariance && this->variances[i].constrainedDependingOnProductOfi.size()>0 )
    {
      double constrainedTo = options.maximumCorrelationCovarianceConstrain;
      for(std::set<int>::iterator it = this->variances[i].constrainedDependingOnProductOfi.begin(); it != this->variances[i].constrainedDependingOnProductOfi.end(); ++it)
      {
        constrainedTo *= this->variances[*it].variance;
      }
      constrainedTo = sqrt(fabs(constrainedTo));
      if( fabs(this->variances[i].variance) > constrainedTo )
      {
        this->variances[i].variance = copysign(constrainedTo, this->variances[i].variance);
	nVariancesConstrained++;
        constrainedVariances.insert(this->variances[i].name);
      }
    }
  }
  
  //Constrain correlations.
  for(int i = 0; i < this->variances.size(); i++)
  {
    if( this->variances[i].type == ParameterAttributes::correlation )
    {
      if( fabs(this->variances[i].variance) > options.maximumCorrelationCovarianceConstrain )
      {
        this->variances[i].variance = copysign(options.maximumCorrelationCovarianceConstrain, this->variances[i].variance);
	nVariancesConstrained++;
        constrainedVariances.insert(this->variances[i].name);
      }
    }
  }
  
  communicator->broadcast(&nVariancesConstrained, 1); 
  broadcastVariances();
  
  /*if( options.varianceConstrainProportion > options.minimumVarianceConstrainProportion )
  {
    options.varianceConstrainProportion = options.varianceConstrainProportion/10.;
    if( options.varianceConstrainProportion < options.minimumVarianceConstrainProportion )
    {
      options.varianceConstrainProportion = options.minimumVarianceConstrainProportion;
    }
  }*/
  
  return nVariancesConstrained;
}

double CovarianceMatrix::constrainVariancesM2(std::vector<double> & oldVariances, std::vector<double> & delta)
{
  int nVariancesToConstrain = 0;
  std::vector<double> scalingFactors;
  std::vector<bool> infiniteScalingFactor;
  
  if(delta.size() != this->variances.size())
  {
    misc.error("Error: An internal error was happened. The dimensions of variances and delta differ.", 0);
  }
  
  for(int i = 0; i < this->variances.size(); i++)
  {
    std::string idGroup = this->variances[i].group;
    double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
    infiniteScalingFactor.push_back(false);
    if((oldVariances[i] + delta[i]) < 0 && this->variances[i].type == ParameterAttributes::variance)
    {
      double tempScalingFactor = -delta[i]/(oldVariances[i] - baseMagnitude*options.varianceConstrainProportion);
      if( isinf(tempScalingFactor) || isnan(tempScalingFactor) )
      {
        infiniteScalingFactor[i] = true;
      }
      else
      {
        scalingFactors.push_back( tempScalingFactor );
      }
      nVariancesToConstrain++;
    }
  }
  
  if(nVariancesToConstrain == 0)
  {
    return 1.;
  }
  
  if(scalingFactors.size() == 0)
  {
    misc.error("Error: When constraining variances, an error was happened. Too many negative variances and a scaling factor can not be computed.", 0);
  }
  
  //Search the minimum scaling factor that don't produce negative variances.
  double minimumScalingFactor = -1.;
  bool flagCheck = false; //Just for checking
  for(int i = 0; i < scalingFactors.size(); i++)
  {
    int flag = true;
    for(int j = 0; j < this->variances.size(); j++)
    {
      std::string idGroup = this->variances[j].group;
      double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
      if((oldVariances[j] + delta[j]) < 0 && this->variances[j].type == ParameterAttributes::variance)
      {
        double temp = oldVariances[j] + (delta[j]/scalingFactors[i]);
        if(temp < 0. && infiniteScalingFactor[j] == false)
        {
          flag = false;
        }
      }
    }
    if(flag && ( (scalingFactors[i] <= minimumScalingFactor) || minimumScalingFactor < 0. ))
    {
      minimumScalingFactor = scalingFactors[i];
      flagCheck = true;
    }
  }
  if(!flagCheck)
  {
    misc.error("Error: An internal error was happened when constraining variances (Method 2).", 0);
  }
  
  for(int i = 0; i < this->variances.size(); i++)
  {
    std::string idGroup = this->variances[i].group;
    double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
    this->variances[i].variance = oldVariances[i] + (delta[i]/minimumScalingFactor);
    if(this->variances[i].variance < 0 && this->variances[i].type == ParameterAttributes::variance)
    {
      if(infiniteScalingFactor[i] == true)
      {
        this->variances[i].variance = baseMagnitude*options.varianceConstrainProportion;
      }
      else
      {
        misc.error("Error: An internal error was happened when constraining variances (Method 2). Variances still negative.", 0);
      }
    }
  }
  
  broadcastVariances();
  
  return minimumScalingFactor;
}

double CovarianceMatrix::constrainVariancesM3(std::vector<double> & oldVariances, std::vector<double> & delta)
{
  int nVariancesToConstrain = 0;
  
  if(delta.size() != this->variances.size())
  {
    misc.error("Error: An internal error was happened. The dimensions of variances and delta differ.", 0);
  }
  
  for(int i = 0; i < this->variances.size(); i++)
  {
    if((oldVariances[i] + delta[i]) < 0. && this->variances[i].type == ParameterAttributes::variance)
    {
      nVariancesToConstrain++;
    }
  }
  
  if( misc.gt(nVariancesToConstrain == 0) )
  {
    return 1.;
  }
  
  //Search the minimum scaling factor that don't produce negative variances.
  double scalingFactor = 1.;
  if(communicator->mpiRoot == true)
  {
    bool variancesToConstrain = true;
    while(variancesToConstrain == true)
    {
      scalingFactor *= options.stepWeightingConstant;
      variancesToConstrain = false;
      for(int i = 0; i < this->variances.size(); i++)
      {
        this->variances[i].variance = oldVariances[i] + (delta[i]*scalingFactor);
        if( this->variances[i].variance < 0. && this->variances[i].type == ParameterAttributes::variance )
        {
          variancesToConstrain = true;
        }
      }
      if( scalingFactor == 0. )
      {
        misc.error("Error: An internal error was happened when constraining variances (Method 3).", 0);
      }
    }
  }
  
  broadcastVariances();
  communicator->broadcast(&scalingFactor);
  
  return scalingFactor;
}

int CovarianceMatrix::fixVariancesToZero()
{
  int nVariancesFixed = 0;
  for(int i = 0; i < this->variances.size(); i++)
  {
    std::string idGroup = this->variances[i].group;
    double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
    if((this->variances[i].variance*(1e20)) < baseMagnitude && this->variances[i].type == ParameterAttributes::variance)
    {
      this->variances[i].variance = 0.;
      nVariancesFixed++;
    }
  }
  
  broadcastVariances();
  communicator->broadcast(&nVariancesFixed);
  
  return nVariancesFixed;
}

std::string CovarianceMatrix::getStringVarianceNames(std::string separator)
{
  std::stringstream result;
  for(int i = 0;  i< this->variances.size(); i++)
  {
    if(options.logFieldWidth > this->variances[i].name.size())
    {
      result << std::setw(options.logFieldWidth) << this->variances[i].name;
    }
    else
    {
      result << (i==0?"":separator) << this->variances[i].name;
    }
  }
  return result.str();
}

std::string CovarianceMatrix::getStringVarianceValues(std::string separator)
{
  std::stringstream result;
  for(int i = 0;  i< this->variances.size(); i++)
  {
    result << std::setprecision(options.logOutputPrecision) << std::setw(options.logFieldWidth) << this->variances[i].variance << ((this->variances[i].fixed == true)?"(f)":"");
  }
  return result.str();
}

std::string CovarianceMatrix::getMatrixName(Matrix *m)
{
  for(int i = 0; i<this->covarianceSubMatrices.size(); i++)
  {
    if(this->covarianceSubMatrices[i] == m)
    {
      return this->covarianceSubMatrixNames[i];
    }
  }
  misc.error("Error: An internal error was happened. In covariance matrix: Matrix not found.", 0);
}

void CovarianceMatrix::broadcastVariances()
{
  double * vv = new double [this->variances.size()];
  
  for(int i = 0; i < this->variances.size(); i++)
  {
    vv[i] = this->variances[i].variance;
  }
  communicator->broadcast(vv, this->variances.size());
  for(int i = 0; i < this->variances.size(); i++)
  {
    this->variances[i].variance = vv[i];
  }
  
  delete [] vv;
}

void CovarianceMatrix::broadcastInitialVariances()
{
  double * vv = new double [this->variances.size()];
  
  for(int i = 0; i < this->variances.size(); i++)
  {
    vv[i] = this->variances[i].initialVariance;
  }
  communicator->broadcast(vv, this->variances.size());
  for(int i = 0; i < this->variances.size(); i++)
  {
    this->variances[i].initialVariance = vv[i];
  }
  
  delete [] vv;
}

std::vector<double> CovarianceMatrix::getVectorVariances()
{
  std::vector<double> vv;
  for(int i = 0; i < this->variances.size(); i++)
  {
    vv.push_back(this->variances[i].variance);
  }
  return vv;
}

void CovarianceMatrix::storeVectorVariances(std::vector<double> & vv)
{
  if(vv.size() != this->variances.size())
  {
    misc.error("Error: An internal error was happened. The size of vector variances is not equal to the argument vector.", 0);
  }
  for(int i = 0; i < this->variances.size(); i++)
  {
    this->variances[i].variance = vv[i];
  }
}

/*
void CovarianceMatrix::getElementMatrixInContext(int i, Matrix * m1)
{
  m1->fillWithConstant(0.);
  
  Element currentElement = this->elements[i];
  
  //Add to the covariance a subcovariance matrix multiplied by the previous product of variances.
  if(currentElement.outcomeSubMatrix.active == false || currentElement.sourceSubMatrix.active == false)
  {
    m1->add(currentElement.m, 0., 1.);
  }
  else
  {
    m1->add(currentElement.m, 0., 1., currentElement.outcomeSubMatrix, currentElement.sourceSubMatrix);
  }
  
  m1->symmetric = true;
  m1->uplo = 'U';
}
*/

Matrix * CovarianceMatrix::getSubCovariance(std::string subCovarianceId)
{
  Matrix * subCovarianceMatrix = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->dimension, this->dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  subCovarianceMatrix->fillWithConstant(0.);
  subCovarianceMatrix->symmetric = true;
  subCovarianceMatrix->uplo = 'B';
  
  bool nameExist = false;
  for(int i=0; i<this->elements.size(); i++)
  {
    if(elements[i].subCovarianceId != subCovarianceId)
    {
      continue;
    }
    nameExist = true;
    
    Element currentElement = this->elements[i];
    
    //Compute the product of the variance elements
    std::pair<double, Matrix *> mf = computeElementMatrixAndFactor(currentElement);
    double currentVariance = mf.first;
    Matrix * melement = mf.second;
    
    //Add to the covariance a subcovariance matrix multiplied by the previous product of variances.
    if(currentElement.outcomeSubMatrix.active == false || currentElement.sourceSubMatrix.active == false)
    {
      subCovarianceMatrix->add(melement, (i==0)?0.:1., currentVariance);
    }
    else
    {
      subCovarianceMatrix->add(melement, (i==0)?0.:1., currentVariance, currentElement.outcomeSubMatrix, currentElement.sourceSubMatrix);
    }
    
    if( currentElement.m != melement && melement != NULL )
    {
      delete melement;
    }
  }
  
  if(nameExist == false)
  {
    misc.error("Error: An internal error was happened. Inexistent subCovarianceId.", 0);
  }
  
  subCovarianceMatrix->symmetric = true;
  subCovarianceMatrix->uplo = 'U';
  
  return subCovarianceMatrix;
}

Matrix * CovarianceMatrix::multiply(std::string subCovarianceId, Matrix * m, double scale)
{
  Matrix * mResult;
  if( this->allCovarianceSubMatricesDiagonal == false )
  {
    mResult = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    Matrix * subCovariance = getSubCovariance(subCovarianceId);
    mResult->multiply(subCovariance, 'N', m, 'N', scale);
    delete subCovariance;
  }
  else
  {
    mResult = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    Matrix * subCovariance = getSubCovariance(subCovarianceId);
    mResult->multiply(subCovariance, 'N', m, 'N', scale);
    delete subCovariance;
  }
  return mResult;
}

void CovarianceMatrix::setVarianceInitialValues(std::vector<Variance> newVariances)
{
  int noInitialize;
  if(communicator->mpiRoot && newVariances.size() == 0)
  {
    noInitialize = 1;
  }
  communicator->broadcast(&noInitialize);
  if(noInitialize == 1)
  {
    return;
  }
  
  if(communicator->mpiRoot)
  {
    if(newVariances.size() != this->variances.size())
    {
      misc.error("Error: The defined initial variances do not agree with REML variances.", 0);
    }
    for(int i = 0; i < newVariances.size(); i++)
    {
      int check = 0;
      for(int j = 0; j < this->variances.size(); j++)
      {
        if(newVariances[i].name == this->variances[j].name)
        {
          this->variances[j].initialVariance = newVariances[i].variance;
          check++;
        }
      }
      if(check != 1)
      {
        misc.error("Error: Error when setting initial variances. Variance names do not agree.", 0);
      }
    }
  }
  broadcastInitialVariances();
}

void CovarianceMatrix::setVarianceInitialValuesFromFile(std::string fileName)
{
  std::vector<Variance> initialVariances;
  
  if(fileName != "")
  {
    misc.message << "Reading initial variance values from file [ " << options.initialVariancesFile << " ]" << std::endl;
    if( communicator->mpiRoot == true )
    {
      std::ifstream file;
      std::string line;
      
      misc.checkFileExists(options.initialVariancesFile);
      file.open(options.initialVariancesFile.c_str());
      
      while(getline(file,line))
      {
        if(!file)
        {
          break;
        }
        if( line.find_first_not_of("\t\r\n ") == std::string::npos )
        {
          continue;
        }
        
        std::istringstream sstemp(line);
        
        Variance variance;
        
        sstemp >> variance.name;
        if( (sstemp >> variance.variance).fail() )
        {
          misc.error("Error: The variance " + variance.name + " has not a valid value in file [ " + options.initialVariancesFile + " ].", 0);
        }
        initialVariances.push_back(variance);
      }
      file.close();
    }
    
    setVarianceInitialValues(initialVariances);
  }
}

void CovarianceMatrix::showCovarianceMatrix(bool showMatrices, bool computedMatrix)
{
  if(showMatrices)
  {
    Matrix * temp = new Matrix(cyclicDistribution, this->dimension, this->dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    for(int i=0; i<this->elements.size(); i++)
    {
      Element currentElement = this->elements[i];
      temp->fillWithConstant(0.);
      temp->add(currentElement.m, 0., 1, currentElement.outcomeSubMatrix, currentElement.sourceSubMatrix);
      temp->showGlobal(getMatrixName(currentElement.m));
    }
    delete temp;
  }
 
  misc.message << "*********************************************" << std::endl;
  misc.message << getStringVarianceNames() << std::endl;
  misc.message << getStringVarianceValues() << std::endl << std::endl;
  
  for(int i=0; i<this->elements.size(); i++)
  {
    Element currentElement = this->elements[i];
    
    std::string mod1 = "";
    std::string mod2 = "";
    
    
    misc.message << currentElement.constantFactor << "*M(" << getMatrixName(currentElement.m) << ")M ";
    for(int j = 0; j<currentElement.iVariances.size(); j++)
    {
      if(currentElement.transformationTypes[j] == ParameterAttributes::squareRoot)
      {
        mod1 = "sqrt(";
        mod2 = ")";
      }
      else if(currentElement.transformationTypes[j] == ParameterAttributes::nochange)
      {
        mod1 = "";
        mod2 = "";
      }
      else if(currentElement.transformationTypes[j] == ParameterAttributes::nochange)
      {
        mod1 = "";
        mod2 = "";
      }
      else
      {
        mod1 = "";
        mod2 = "";
      }
      misc.message << "*" << mod1 << this->variances[currentElement.iVariances[j]].name << mod2;
    }
    misc.message << ((i!=(elements.size()-1))?" + ":"");
  }
  misc.message << std::endl;
  misc.message << "*********************************************" << std::endl;
  
  if(computedMatrix == true)
  {
    computeCovariance();
    this->m->showGlobal("Full Covariance", false);
  }
}
