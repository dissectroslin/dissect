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

#include "covariate.h"
#include "matrix.h"
#include "misc.h"
#include "global.h"
#include "options.h"
#include "auxiliar.h"
#include "communicator.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>

Covariate::Covariate(std::string fcov, std::string fqcov, std::vector<std::string> emptyIndividualIds, bool constructMatrix, std::set<std::string> prefilterIndividuals)
{
  this->covariates = NULL;
  
  this->nIndividuals = 0;
  this->nCovariates = 0;
  this->nQuantitativeCovariates = 0;
  this->nDiscreteCovariates = 0;
  
  readRawCovariate(fqcov, this->rawQCovars, prefilterIndividuals);
  readRawCovariate(fcov, this->rawCovars, prefilterIndividuals);
  this->rawCleared = false;
  getDiscreteCovariateCategories(this->rawCovars);
  
  if(constructMatrix == true)
  {
    parseRawCovariates(emptyIndividualIds);
  }
}

Covariate::~Covariate()
{
  if(this->covariates!=NULL)
  {
    delete this->covariates;
    this->covariates = NULL;
  }
}

void Covariate::parseRawCovariates(std::vector<std::string> emptyIndividualIds, int nMeans, int idxThisMean)
{
  if(this->rawCleared == true)
  {
    misc.error("Error: An internal error was happened. The covariates can not be parsed when raw vectors are empty.", 0);
  }

  std::vector< std::vector<double> > qCovars;
  std::vector< std::vector<double> > Covars;
  
  reestructureQuantitativeCovariate(this->rawQCovars, qCovars);
  reestructureDiscreteCovariateUsingDifferences(this->rawCovars, Covars);
  
  //Compute the total number of covariates and individuals
  //this->nIndividuals = Covars.size()!=0?Covars.size():qCovars.size();
  this->nCovariates = 0;
  this->nCovariates +=nMeans;
  this->nCovariates += this->nDiscreteCovariates;
  this->nCovariates += this->nQuantitativeCovariates;
  
  this->meanNames.clear();
  for(int i = 0; i < nMeans; i++)
  {
    std::stringstream ss;
    ss << "mean-" << i;
    this->meanNames.push_back(ss.str());
  }
  
  //If no covariates are specified, then, the number of individuals is 0.
  //However the mean matrix must be created. To this end, the individuals in emptyIndividuals are used.
  if(this->nIndividuals == 0)
  {
    if(communicator->mpiRoot == true)
    {
      this->individualIds = emptyIndividualIds;
      for(int i = 0; i < this->individualIds.size(); i++)
      {
        individualIdsIdx[ this->individualIds[i] ] = i;
      }
      this->nIndividuals = this->individualIds.size();
    }
    communicator->broadcast(&this->nIndividuals, 1);
  }
  
  //Create the covariance matrix
  constructCovariateMatrix(Covars, qCovars, nMeans, idxThisMean);
  
  clearRaw();
}

void Covariate::clearRaw()
{
  this->rawCovars.clear();
  this->rawQCovars.clear();
  this->rawCleared = true;
}

void Covariate::readRawCovariate(std::string f, std::vector< std::vector<std::string> > & covars, std::set<std::string> prefilterIndividuals)
{
  std::ifstream file;
  std::string line;
  
  if(f == "")
  {
    covars.clear();
    return;
  }
  
  if(communicator->mpiRoot)
  {
    misc.message << "Reading covariates data from file [ " << f << " ] ..." << std::endl;
    
    misc.checkFileExists(f);
    file.open(f.c_str());
    
    covars.clear();
    int nCovars = -1;
    int idx = 0;
    while(getline(file,line))
    {
      if(!file)
      {
	break;
      }
      
      std::istringstream sstemp(line);
      
      std::string familyId, individualId;
      sstemp >> familyId >> individualId;
      std::string key = familyId + "@" + individualId;
      
      if( prefilterIndividuals.size() != 0 && prefilterIndividuals.find(key) == prefilterIndividuals.end() )
      {
        continue;
      }
      
      if(this->nIndividuals != 0) //If a file is readed before, nIndividuals should be != 0 and the individualIds between two files must be the same.
      {
	if(idx>=this->individualIds.size())
	{
	  misc.error("Error: There are not the same individuals in covariate and quantitative covariate files.", 0);
	}
	if(this->individualIds[idx] != familyId + "@" + individualId)
	{
	  misc.error("Error: There are not the same individuals in covariate and quantitative covariate files.", 0);
	}
      }
      else
      {
	
	this->individualIds.push_back(key);
	if(this->individualIdsIdx.count(key) != 0)
	{
	  misc.error("Error: The individual with family Id: " + familyId + " and individual Id: " + individualId + " appears more than one time in file [ " + f + " ].", 0);
	}
	this->individualIdsIdx[key] = idx;
      }

      std::vector<std::string> covarsLine;
      std::string sCovar;
      
      //Read all de covariates for the current individual
      while(sstemp >> sCovar)
      {
	covarsLine.push_back(sCovar);
        if(sCovar == "NA")
        {
          this->individualIdsWithMissingData.insert(key);
        }
      }
      
      //Check that each individual has the same number of covariates
      if(nCovars == -1)
      {
	nCovars = covarsLine.size();
      }
      if(nCovars != covarsLine.size())
      {
	misc.error("Error: The covariables file is not properly formated. The individual " + familyId + " " + individualId + " has less/more covariables than previous individuals.", 0);
      }

      covars.push_back(covarsLine);
      idx++;
    }
    file.close();
    
    if(covars.size()==0)
    {
      misc.error("Error: The covariates file [ " + f + " ] is empty.", 0);
    }
    misc.message << covars[0].size() << " covariates have been read for " << covars.size() << " individuals." << std::endl;
    
    if(this->nIndividuals != 0 && this->nIndividuals != this->individualIds.size())
    {
      misc.error("Error: There are not the same individuals in covariate and quantitative covariate files.", 0);
    }
    this->nIndividuals = this->individualIds.size();
  }
  communicator->broadcast(&this->nIndividuals, 1);
}

void Covariate::getDiscreteCovariateCategories(std::vector< std::vector<std::string> > & rawCovars)
{
  if(this->rawCleared == true)
  {
    misc.error("Error: An internal error was happened. The covariate categories can not be parsed when raw vectors are empty.", 0);
  }
  
  //Get all the posible values of the discrete covars.
  if(communicator->mpiRoot)
  {
    if(rawCovars.size() == 0)
    {
      return;
    }
    
    int nTotalValuesCovar = 0;
    this->covarCategories.clear();
    for(int j=0; j<rawCovars[0].size(); j++)
    {
      std::map<std::string, int> temp;
      int nValuesCurrentCovar = 0;
      for(int i=0; i<rawCovars.size(); i++)
      {
        if(temp.count(rawCovars[i][j]) != 0)
        {
          continue;
        }
        temp[ rawCovars[i][j] ] = nValuesCurrentCovar;
        nValuesCurrentCovar++;
        nTotalValuesCovar++;
      }
      if(double(temp.size()) > double(this->nIndividuals)*options.maximumProportionOfCategories)
      {
        misc.error("Error: There are too many categories in the " + i2s( j + 1 ) + " column of covariates file.", 0);
      }
      this->covarCategories.push_back(temp);
    }
    this->nDiscreteCovariates = nTotalValuesCovar;
  }
}

void Covariate::syncronizeDiscreteCovariateCategoriesWith(Covariate * covariate)
{
  if(this->rawCleared == true)
  {
    misc.error("Error: An internal error was happened. Covariates can not be syncronized when raw vectors are empty.", 0);
  }
  
  if(communicator->mpiRoot)
  {
    std::vector< std::map<std::string, int> > combinedCovarCategories;
    if( this->covarCategories.size() != covariate->covarCategories.size() )
    {
      misc.error("Error: Error when combining two covariate files. It seems that their columns are not the same or not represent the same thing.", 0);
    }
    
    int nTotalValuesCovar = 0;
    for(int i = 0; i < this->covarCategories.size(); i++)
    {
      //get the categories in both Covariate objects
      std::set<std::string> categories;
      for(std::map<std::string, int>::iterator it = this->covarCategories[i].begin(); it != this->covarCategories[i].end(); ++it)
      {
        categories.insert(it->first);
      }
      for(std::map<std::string, int>::iterator it = covariate->covarCategories[i].begin(); it != covariate->covarCategories[i].end(); ++it)
      {
        categories.insert(it->first);
      }
      //index the categories
      std::map<std::string, int> indexedCategories;
      int idx = 0;
      for(std::set<std::string>::iterator it = categories.begin(); it != categories.end(); ++it)
      {
        indexedCategories[*it] = idx;
        idx++;
        nTotalValuesCovar++;
      }
      //Append to the new vector of covarCategories;
      combinedCovarCategories.push_back(indexedCategories);
    }
    
    this->covarCategories = combinedCovarCategories;
    covariate->covarCategories = combinedCovarCategories;
    this->nDiscreteCovariates = nTotalValuesCovar;
    covariate->nDiscreteCovariates = nTotalValuesCovar;
  }
}

void Covariate::reestructureQuantitativeCovariate(std::vector< std::vector<std::string> > & rawCovars, std::vector< std::vector<double> > & covars)
{
  covars.clear();
  if(communicator->mpiRoot)
  {
    this->quantitativeCovarNames.clear();
    if(rawCovars.size() != 0)
    {
      this->nQuantitativeCovariates = rawCovars[0].size();
    }
    else
    {
      return;
    }
    
    //Parse data in rawCovars
    std::vector<double> means(rawCovars[0].size(), 0.);
    std::vector<double> nNonMissing(rawCovars[0].size(), 0.);
    std::vector<int> allMissing(rawCovars[0].size(), 1);
    for(int i=0; i<rawCovars.size(); i++)
    {
      std::vector<double> covarsLine;
      for(int j=0; j<rawCovars[i].size(); j++)
      {
        double covar;
        if(rawCovars[i][j] != "NA")
        {
          std::istringstream sstemp2(rawCovars[i][j]);
          if( (sstemp2 >> covar).fail() )
          {
            misc.error("Error: The individual " + this->individualIds[i] + " has not a valid value in one of their covariables: " + rawCovars[i][j], 0);
          }
          covarsLine.push_back(covar);
          means[j] += covar;
          nNonMissing[j] += 1;
          allMissing[j] = 0;
        }
        else
        {
          covarsLine.push_back(0.);
        }
      }
      covars.push_back(covarsLine);
    }
    
    //Check there are no columns with all missings.
    for(int j=0; j<rawCovars[0].size(); j++)
    {
      if(allMissing[j] == 1)
      {
        misc.error("Error: A column in the quantitative covariates file have all values missing.", 0);
      }
    }
    
    //Replace NA values with the mean
    for(int i=0; i<rawCovars.size(); i++)
    {
      for(int j=0; j<rawCovars[i].size(); j++)
      {
        if(rawCovars[i][j] == "NA")
        {
          covars[i][j] = means[j]/nNonMissing[j];
        }
      }
    }
    
    //Store names in order
    for(int j=0; j<rawCovars[0].size(); j++)
    {
      std::stringstream ss;
      ss << "q-col-" << j;
      this->quantitativeCovarNames.push_back(ss.str());
    }
  }
}


void Covariate::reestructureDiscreteCovariate(std::vector< std::vector<std::string> > & rawCovars, std::vector< std::vector<double> > & covars)
{
  covars.clear();
  if(communicator->mpiRoot)
  {
    if(rawCovars.size() == 0)
    {
      return;
    }
   
    //Generate the matrix of dummy variables (each value for each covariate has a column. Each individual has a 1 on the columns that match the their value and category, 0 otherwise)
    for(int i=0; i<rawCovars.size(); i++)
    {
      std::vector<double> temp(this->nDiscreteCovariates, 0);
      int shift = 0;
      for(int j=0; j<rawCovars[0].size(); j++)
      {
	if(this->covarCategories[j].count(rawCovars[i][j]) == 0)
	{
	  misc.error("Error: An internal error was happened.", 0);
	}
	temp[ shift + this->covarCategories[j][ rawCovars[i][j] ] ] = 1.;
	shift += this->covarCategories[j].size();
      }
      covars.push_back(temp);
    }
    
    //Store names in order
    int namesShift = 0;
    this->discreteCovarNames = std::vector<std::string>(this->nDiscreteCovariates, "");
    for(int j=0; j<rawCovars[0].size(); j++)
    {
      for(std::map<std::string, int>::iterator it = this->covarCategories[j].begin(); it != this->covarCategories[j].end(); ++it)
      {
        std::stringstream ss;
        ss << "col-" << j << "(" << it->first << ")";
        this->discreteCovarNames[ namesShift + it->second ] = ss.str();
      }
      namesShift += this->covarCategories[j].size();
    }
  }
}

void Covariate::reestructureDiscreteCovariateUsingDifferences(std::vector< std::vector<std::string> > & rawCovars, std::vector< std::vector<double> > & covars)
{
  covars.clear();
  if(communicator->mpiRoot)
  {
    if(rawCovars.size() == 0)
    {
      return;
    }
    
    this->nDiscreteCovariates = this->nDiscreteCovariates - rawCovars[0].size(); //The first category of each covariate is removed.
    
    //Generate the matrix of dummy variables (each value for each covariate has a column. Each individual has a 1 on the columns that match the their value and category, 0 otherwise)
    for(int i=0; i<rawCovars.size(); i++)
    {
      std::vector<double> temp(this->nDiscreteCovariates, 0);
      int shift = 0;
      for(int j=0; j<rawCovars[0].size(); j++)
      {
        if(this->covarCategories[j][ rawCovars[i][j] ] == 0) //The first category will not be included.
        {
          shift += this->covarCategories[j].size() - 1;
          continue;
        }
        if(this->covarCategories[j].count(rawCovars[i][j]) == 0)
        {
          misc.error("Error: An internal error was happened.", 0);
        }
        temp[ shift + this->covarCategories[j][ rawCovars[i][j] ] - 1 ] = 1.;
        shift += this->covarCategories[j].size() - 1;
      }
      covars.push_back(temp);
    }
    
    //Store names in order
    int namesShift = 0;
    this->discreteCovarNames = std::vector<std::string>(this->nDiscreteCovariates, "");
    for(int j=0; j<rawCovars[0].size(); j++)
    {
      std::string baseCategoryName = "";
      for(std::map<std::string, int>::iterator it = this->covarCategories[j].begin(); it != this->covarCategories[j].end(); ++it)
      {
        if(it->second == 0)
        {
          baseCategoryName = it->first;
        }
      }
      for(std::map<std::string, int>::iterator it = this->covarCategories[j].begin(); it != this->covarCategories[j].end(); ++it)
      {
        if(it->second == 0)
        {
          continue;
        }
        std::stringstream ss;
        ss << "col-" << j << "(" << it->first << "/" << baseCategoryName << ")";
        this->discreteCovarNames[ namesShift + it->second - 1 ] = ss.str();
      }
      namesShift += this->covarCategories[j].size() - 1;
    }
  }
}

void Covariate::constructCovariateMatrix(std::vector< std::vector<double> > & Covars, std::vector< std::vector<double> > & qCovars, int nMeans, int idxThisMean)
{
  double * covariatesTempMatrix;
  
  if (communicator->mpiRoot) {
    if(Covars.size()!=0 && qCovars.size()!=0 && Covars.size() != qCovars.size())
    {
      misc.error("Error: The number of individuals in the covariates file differ from the number of individuals in the quantitative covariate file.", 0);
    }
    
    //Copy Covars and qCovars to covariatesTempMatrix
    covariatesTempMatrix = new double [this->nCovariates*this->nIndividuals];
    
    int columnShift = 0;
    for(int i = 0; i < nMeans; i++)
    {
      double value = ((i == idxThisMean)?1.:0.);
      for(int r = 0; r < this->nIndividuals; r++)
      {
        for(int c = 0; c < 1; c++)
        {
          covariatesTempMatrix[(c+columnShift)*this->nIndividuals + r] = value;
        }
      }
      columnShift++;
    }
    if(Covars.size() != 0)
    {
      for(int r = 0; r < Covars.size(); r++)
      {
	for(int c = 0; c < Covars[r].size(); c++)
	{
	  covariatesTempMatrix[(c+columnShift)*this->nIndividuals + r] = Covars[r][c];
	}
      }
      columnShift += Covars[0].size();
    }
    if(qCovars.size() != 0)
    {
      for(int r = 0; r < qCovars.size(); r++)
      {
	for(int c = 0; c < qCovars[r].size(); c++)
	{
	  covariatesTempMatrix[(c+columnShift)*this->nIndividuals + r] = qCovars[r][c];
	}
      }
      columnShift += qCovars[0].size();
    }
  }
  
  communicator->broadcast(&this->nIndividuals, 1);
  communicator->broadcast(&this->nCovariates, 1);
  
  //Distribute covariatesTempMatrix between nodes.
  if(this->covariates!=NULL)
  {
    delete this->covariates;
    this->covariates = NULL;
  }
  
  this->covariates = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->nIndividuals, this->nCovariates, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->covariates->scatterMatrix(covariatesTempMatrix);
  
  if (communicator->mpiRoot) {
    delete [] covariatesTempMatrix;
  }
}

void Covariate::filterIndividuals(std::vector<std::string> keepIndividualIds)
{
  double * originalGlobMatrix;
  double * filteredGlobMatrix;
  if(communicator->mpiRoot)
  {
    originalGlobMatrix = new double [this->covariates->nGlobRows*this->covariates->nGlobCols];
    filteredGlobMatrix = new double [keepIndividualIds.size()*this->nCovariates];
  }
  
  this->covariates->gatherMatrix(originalGlobMatrix);
  
  if(communicator->mpiRoot)
  {
    std::map<std::string, int> newIndividualIdsIdx;
    std::set<std::string> newIndividualIdsWithMissingData;
    
    for(int r=0; r<keepIndividualIds.size(); r++) //Is it better access rows contiguously and searxh in map each step? Or inverse?
    {
      if(this->individualIdsIdx.count(keepIndividualIds[r]) == 0)
      {
	misc.error("Error: An internal error was happened. The individual '" + keepIndividualIds[r]  + "' is not in the covariance matrix.", 0);
      }
      if(newIndividualIdsIdx.count(keepIndividualIds[r]) != 0)
      {
	misc.error("Error: An internal error was happened. The individual '" + keepIndividualIds[r]  + "' is repeated.", 0);
      }
      if(this->individualIdsWithMissingData.find(keepIndividualIds[r]) != this->individualIdsWithMissingData.end())
      {
        newIndividualIdsWithMissingData.insert(keepIndividualIds[r]);
      }
      
      int oldRowIdx = this->individualIdsIdx[keepIndividualIds[r]];
      for(int c=0; c<this->nCovariates; c++)
      {
	filteredGlobMatrix[c*keepIndividualIds.size() + r] = originalGlobMatrix[c*this->nIndividuals + oldRowIdx];
      }
      newIndividualIdsIdx[keepIndividualIds[r]] = r;
    }
    this->nIndividuals = keepIndividualIds.size();
    this->individualIds = keepIndividualIds;
    this->individualIdsIdx = newIndividualIdsIdx;
    
    this->individualIdsWithMissingData = newIndividualIdsWithMissingData;
  }
  communicator->broadcast(&this->nIndividuals, 1);
  communicator->barrier();
  
  this->covariates->initParameters(this->nIndividuals, this->nCovariates, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->covariates->scatterMatrix(filteredGlobMatrix);
  
  if(communicator->mpiRoot)
  {
    delete [] originalGlobMatrix;
    delete [] filteredGlobMatrix;
  }
}

std::map<std::string, double> Covariate::loadEffectPrediction(std::string discretefname, std::string quantitativefname)
{
  if(this->rawCleared == true)
  {
    misc.error("Error: An internal error was happened. The covariate effects can not be estimated when raw vectors are empty.", 0);
  }
  
  
  std::map< int, std::map<std::string, double> > covarEffects = loadDiscreteEffects(discretefname);
  std::vector<double> qCovarEffects = loadQuantitativeEffects(quantitativefname);
  
  if( this->rawCovars.size() != 0 )
  {
    if( covarEffects.size() != this->rawCovars[0].size() )
    {
      misc.error("Error: An error was happened. The number of columns in the discrete covariate files disagree with the number of columns in the effect files.", 0);
    }
  }
  if( this->rawQCovars.size() != 0 )
  {
    if( qCovarEffects.size() != this->rawQCovars[0].size() )
    {
      misc.error("Error: An error was happened. The number of columns in the quantitative covariate files disagree with the number of columns in the effect files.", 0);
    }
  }
  
  std::map<std::string, double> result;
  
  if( misc.gt( covarEffects.size() == 0 && qCovarEffects.size() == 0 ))
  {
    return result;
  }
  
  if(communicator->mpiRoot == true)
  {
    std::vector<double> values;
    //Add effects from discrete covariates
    if( this->rawCovars.size() != 0 )
    {
      for( int i = 0; i<this->nIndividuals; i++ )
      {
        double value = 0.;
        for ( int col = 0; col<this->rawCovars[i].size(); col++ )
        {
          if( covarEffects[col].count( this->rawCovars[i][col] ) != 0 )
          {
            value += covarEffects[col][ this->rawCovars[i][col] ];
          }
          else
          {
            if( options.forceUseNonEstimatedCovarKeys == false )
            {
              misc.error("Error: In the file of discrete covariates there is the key " + this->rawCovars[i][col] + " which is not in the effects file. You can force the effect of this key to be 0. by using the flag --force-use-unestimated-values.", 0);
            }
          }
        }
        values.push_back(value);
      }
    }
    else
    {
      values.assign(this->nIndividuals, 0.);
    }
    
    //Add effects from quantitative covariates
    std::vector< std::vector<double> > qCovars;
    reestructureQuantitativeCovariate(this->rawQCovars, qCovars);
    if( qCovars.size() != 0 )
    {
      for( int i = 0; i<this->nIndividuals; i++ )
      {
        for ( int col = 0; col<qCovars[i].size(); col++ )
        {
          values[i] += qCovars[i][col]*qCovarEffects[col];
        }
      }
    }
    
    //Sort results by individual Id
    for( int i = 0; i<this->individualIds.size(); i++ )
    {
      result[ this->individualIds[i] ] = values[i];
    }
  }
  
  this->nQuantitativeCovariates = 0;
  this->quantitativeCovarNames.clear();
  
  return result;
}

std::map< int, std::map<std::string, double> > Covariate::loadDiscreteEffects(std::string fname)
{
  std::ifstream file;
  std::string line;
  std::map< int, std::map<std::string, double> > effects;
  std::map< int, std::string > baseEffect;
  
  if(fname == "")
  {
    return effects;
  }
  
  if(communicator->mpiRoot)
  {
    misc.message << "Reading Quantitative covariate effects from file [ " << fname << " ]..." << std::endl;
    
    std::vector< std::vector<std::string> > table;
    getTableFromFile(fname, table, 3);
    
    if(table.size() <= 1)
    {
      return effects;
    }
    if(table[0][0] != "NAME" ||  table[0][1] != "BETA" || table[0][2] != "STD")
    {
      return loadDiscreteEffectsOld(fname);
    }
    
    int nEffects = 0;
    for(int ei = 1; ei < table.size(); ei++)
    {
      //std::istringstream sstemp(line); //Check that number of words of line <= that column + 2.
      
      std::string covarId = table[ei][0];
      std::vector<std::string> substrings = splitString(covarId, "(");
      if(substrings.size() != 2)
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + covarId + " is not a valid key.", 0);
      }
      std::vector<std::string> temp = splitString(substrings[0], "-");
      if(temp.size() < 2)
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + covarId + " is not a valid key.", 0);
      }
      std::istringstream sscolumn( temp[ temp.size() - 1 ] );
      int column;
      if( (sscolumn >> column).fail() )
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + covarId + " is not a valid key.", 0);
      }
      std::vector<std::string> subkeys = splitString(substrings[1], ")");
      if( subkeys.size() != 2 )
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + covarId + " is not a valid key.", 0);
      }
      subkeys = splitString(subkeys[0], "/");
      if( subkeys.size() != 2 )
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + covarId + " is not a valid key.", 0);
      }
      std::string baseKey = subkeys[1];
      std::string currentKey = subkeys[0];
      
      if(baseEffect.count(column) == 0)
      {
        baseEffect[column] = baseKey;
        effects[column] = std::map<std::string, double>();
      }
      else
      {
        if(baseEffect[column] != baseKey)
        {
          misc.error("Error: An error was happened when reading [ " + fname + " ]. Column " + i2s(column) + " has defined more than one base key.", 0);
        }
      }
      
      double effect = string2Number<double>(table[ei][1], "Error: The discrete covariate file is not properly formatted. The key: " + covarId + " has not a valid value.");
      
      if(effects[column].count(currentKey) != 0)
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. Column " + i2s(column) + " has a duplicated key.", 0);
      }
      effects[column][currentKey] = effect;
      nEffects++;
    }
    
    for(int col = 0; col<baseEffect.size(); col++)
    {
      if( baseEffect.count(col) == 0 || effects.count(col) == 0 )
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. Column " + i2s(col) + " is not present.", 0);
      }
      if( effects[col].count(baseEffect[col]) != 0 )
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. There is a key in column " + i2s(col) + " which is equal to the base effect.", 0);
      }
      effects[col][baseEffect[col]] = 0.;
    }
    
    misc.message << nEffects << " effects have been read from file [ " + fname + " ]." << std::endl;
  }
  
  return effects;
}

std::vector<double> Covariate::loadQuantitativeEffects(std::string fname)
{
  std::ifstream file;
  std::string line;
  std::vector<double> effects;
  
  if(fname == "")
  {
    return effects;
  }
  
  if(communicator->mpiRoot)
  {
    misc.message << "Reading Quantitative covariate effects from file [ " << fname << " ]..." << std::endl;
    
    std::vector< std::vector<std::string> > table;
    getTableFromFile(fname, table, 3);
    
    if(table.size() <= 1)
    {
      return effects;
    }
    if(table[0][0] != "NAME" ||  table[0][1] != "BETA" || table[0][2] != "STD")
    {
      return loadQuantitativeEffectsOld(fname);
    }
    
    int idx = 0;
    for(int ei = 1; ei < table.size(); ei++)
    {
      //std::istringstream sstemp(line); //Check that number of words of line <= that column + 2.
      
      std::string qCovarId = table[ei][0];
      std::vector<std::string> temp = splitString(qCovarId, "-");
      if( temp.size() < 3)
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + qCovarId + " is not a valid key.", 0);
      }
      if(temp[temp.size()-3] != "q" || temp[temp.size()-2] != "col" || temp[temp.size()-1] != i2s(idx))
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + qCovarId + " is not a valid key.", 0);
      }
      
      double effect = string2Number<double>(table[ei][1], "Error: The quantitative covariate file is not properly formatted. The key: " + qCovarId + " has not a valid value.");
      
      effects.push_back(effect);
      idx++;
    }
    file.close();
    
    misc.message << effects.size() << " effects have been read from file [ " + fname + " ]." << std::endl;
  }
  
  return effects;
}

std::map< int, std::map<std::string, double> > Covariate::loadDiscreteEffectsOld(std::string fname)
{
  std::ifstream file;
  std::string line;
  std::map< int, std::map<std::string, double> > effects;
  std::map< int, std::string > baseEffect;
  
  if(fname == "")
  {
    return effects;
  }
  
  if(communicator->mpiRoot)
  {
    misc.message << "Reading Quantitative covariate effects from file [ " << fname << " ]..." << std::endl;
    
    misc.checkFileExists(fname);
    file.open(fname.c_str());
    
    int nEffects = 0;
    while(getline(file,line))
    {
      if(!file)
      {
        break;
      }
      
      std::istringstream sstemp(line); //Check that number of words of line <= that column + 2.
      
      std::string covarId;
      sstemp >> covarId;
      std::vector<std::string> substrings = splitString(covarId, "(");
      if(substrings.size() != 2)
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + covarId + " is not a valid key.", 0);
      }
      std::vector<std::string> temp = splitString(substrings[0], "-");
      if(temp.size() < 2)
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + covarId + " is not a valid key.", 0);
      }
      std::istringstream sscolumn( temp[ temp.size() - 1 ] );
      int column;
      if( (sscolumn >> column).fail() )
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + covarId + " is not a valid key.", 0);
      }
      std::vector<std::string> subkeys = splitString(substrings[1], ")");
      if( subkeys.size() != 2 )
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + covarId + " is not a valid key.", 0);
      }
      subkeys = splitString(subkeys[0], "/");
      if( subkeys.size() != 2 )
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + covarId + " is not a valid key.", 0);
      }
      std::string baseKey = subkeys[1];
      std::string currentKey = subkeys[0];
      
      if(baseEffect.count(column) == 0)
      {
        baseEffect[column] = baseKey;
        effects[column] = std::map<std::string, double>();
      }
      else
      {
        if(baseEffect[column] != baseKey)
        {
          misc.error("Error: An error was happened when reading [ " + fname + " ]. Column " + i2s(column) + " has defined more than one base key.", 0);
        }
      }
      
      double effect;
      std::string seffect;
      if((sstemp >> seffect).fail())
      {
        misc.error("Error: The discrete covariate file is not properly formatted for key: " + covarId + ".", 0);
      }
      
      sstemp.clear();
      sstemp.str(seffect);
      if( (sstemp >> effect).fail() )
      {
        misc.error("Error: The discrete covariate file is not properly formatted. The key: " + covarId + " has not a valid value.", 0);
      }
      
      if(effects[column].count(currentKey) != 0)
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. Column " + i2s(column) + " has a duplicated key.", 0);
      }
      effects[column][currentKey] = effect;
      nEffects++;
    }
    file.close();
    
    for(int col = 0; col<baseEffect.size(); col++)
    {
      if( baseEffect.count(col) == 0 || effects.count(col) == 0 )
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. Column " + i2s(col) + " is not present.", 0);
      }
      if( effects[col].count(baseEffect[col]) != 0 )
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. There is a key in column " + i2s(col) + " which is equal to the base effect.", 0);
      }
      effects[col][baseEffect[col]] = 0.;
    }
    
    misc.message << nEffects << " effects have been read from file [ " + fname + " ]." << std::endl;
  }
  
  return effects;
}

std::vector<double> Covariate::loadQuantitativeEffectsOld(std::string fname)
{
  std::ifstream file;
  std::string line;
  std::vector<double> effects;
  
  if(fname == "")
  {
    return effects;
  }
  
  if(communicator->mpiRoot)
  {
    misc.message << "Reading Quantitative covariate effects from file [ " << fname << " ]..." << std::endl;
    
    misc.checkFileExists(fname);
    file.open(fname.c_str());
    
    int idx = 0;
    while(getline(file,line))
    {
      if(!file)
      {
        break;
      }
      
      std::istringstream sstemp(line); //Check that number of words of line <= that column + 2.
      
      std::string qCovarId;
      sstemp >> qCovarId;
      std::vector<std::string> temp = splitString(qCovarId, "-");
      if( temp.size() < 3)
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + qCovarId + " is not a valid key.", 0);
      }
      if(temp[temp.size()-3] != "q" || temp[temp.size()-2] != "col" || temp[temp.size()-1] != i2s(idx))
      {
        misc.error("Error: An error was happened when reading [ " + fname + " ]. " + qCovarId + " is not a valid key.", 0);
      }
      
      double effect;
      std::string seffect;
      if((sstemp >> seffect).fail())
      {
        misc.error("Error: The quantitative covariate file is not properly formatted for key: " + qCovarId + ".", 0);
      }
      
      sstemp.clear();
      sstemp.str(seffect);
      if( (sstemp >> effect).fail() )
      {
        misc.error("Error: The quantitative covariate file is not properly formatted. The key: " + qCovarId + " has not a valid value.", 0);
      }
      effects.push_back(effect);
      idx++;
    }
    file.close();
    
    misc.message << effects.size() << " effects have been read from file [ " + fname + " ]." << std::endl;
  }
  
  return effects;
}


void Covariate::printCovariate(int nSetw)
{
  double *cov;
  
  if (communicator->mpiRoot) {
    cov = new double [this->nIndividuals*this->nCovariates];
  }
  
  this->covariates->gatherMatrix(cov);
  if (communicator->mpiRoot) {
    misc.message << "Covariates Matrix (" << this->covariates->nGlobRows << ", " << this->covariates->nGlobCols << "):\n";
    misc.message << this->individualIds.size() << " " << this->individualIdsIdx.size() << std::endl ;
    for (int r = 0; r < this->nIndividuals; r++)
    {
      misc.message << this->individualIds[r] << " " << this->individualIdsIdx[this->individualIds[r]] << ": " ;
      for (int c = 0; c < this->nCovariates; c++)
      {
	misc.message << cov[c*this->nIndividuals + r] << " ";
      }
      misc.message << std::endl;
    }
    misc.message << std::endl;
    
    for (int r = 0; r < this->nIndividuals; r++)
    {
      misc.message << this->individualIds[r] << " " << this->individualIdsIdx[this->individualIds[r]] << ": mu " ;
      int shift = 0;
      for (int c = 0; c < this->meanNames.size(); c++)
      {
        misc.message << std::setw(nSetw) << cov[(c + shift)*this->nIndividuals + r] << " ";
      }
      shift += this->meanNames.size();
      for (int c = 0; c < this->nDiscreteCovariates; c++)
      {
        if(cov[(c + shift)*this->nIndividuals + r] == 1)
        {
          misc.message << std::setw(nSetw) << discreteCovarNames[c];
        }
        else
        {
          misc.message << std::setw(nSetw) << "-";
        }
      }
      shift += this->nDiscreteCovariates;
      for (int c = 0; c < this->nQuantitativeCovariates; c++)
      {
        misc.message << std::setw(nSetw) << cov[(c + shift)*this->nIndividuals + r] << " ";
      }
      misc.message << std::endl;
    }
    misc.message << std::endl;
    
    misc.message.flush();
    
    misc.message << "With missings:" << std::endl;
    for(std::set<std::string>::iterator it = this->individualIdsWithMissingData.begin(); it != this->individualIdsWithMissingData.end(); ++it)
    {
      misc.message << *it << std::endl;
    }
  }
  
  if (communicator->mpiRoot) {
    delete [] cov;
  }
}
