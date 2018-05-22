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

#include "phenotype.h"
#include "matrix.h"
#include "misc.h"
#include "global.h"
#include "options.h"
#include "communicator.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

Phenotype::Phenotype(DistributionType dist, std::string f, int column)
{
  this->phenotypes = NULL;
  this->nIndividuals = 0;
  this->nPhenotypesInFile = 0;

  getNumberPhenotypesInFile(f);
  
  if(this->nPhenotypesInFile < column)
  {
    misc.error("Error: Invalid phenotype column. Phenotypes file [ " + f + " ] does not have enough columns.", 0);
  }
  
  if(column < 1)
  {
    misc.error("Error: Invalid column when reading phenotypes file [ " + f + " ].", 0);
  }
  
  readPhenotype(dist, f, column);
}

Phenotype::~Phenotype()
{
  if(this->phenotypes!=NULL)
  {
    delete this->phenotypes;
    this->phenotypes = NULL;
  }
}

void Phenotype::getNumberPhenotypesInFile(std::string f)
{
  if(communicator->mpiRoot)
  {
    std::ifstream file;
    std::string line;
  
    misc.checkFileExists(f);
    file.open(f.c_str());
    
    getline(file,line);
    if(!file)
    {
      misc.error("Error: Phenotypes file [ " + f + " ] is empty.", 0);
    }
    
    std::istringstream sstemp(line);
    std::string dummy;
    this->nPhenotypesInFile = 0;
    while(sstemp >> dummy)
    { 
      this->nPhenotypesInFile++;
    }
    this->nPhenotypesInFile -= 2;
    
    if(this->nPhenotypesInFile < 1)
    {
      misc.error("Error: Phenotypes file [ " + f + " ] is not properly formated.", 0);
    }
    
    file.close();
  }
  communicator->broadcast(&this->nPhenotypesInFile);
}

void Phenotype::readPhenotype(DistributionType dist, std::string f, int column)
{
  std::ifstream file;
  std::string line;
  std::vector<double> phenosArrayTemp;
  
  if(communicator->mpiRoot)
  {
    int nRemovedPhenotypesFromOptions = 0;
    misc.message << "Reading Phenotype data from column " << column << " of file [ " << f << " ]..." << std::endl;
    if(options.keysIndividualsToKeep.size() != 0)
    {
      misc.message << "All individuals not present in file [ " << options.fileIndividualsToKeep << " ] will be removed." << std::endl;
    }
    
    misc.checkFileExists(f);
    file.open(f.c_str());
    
    this->individualIds.clear();
    this->missings.clear();
    this->individualIdsIdx.clear();
    phenosArrayTemp.clear();
    int idx = 0;
    while(getline(file,line))
    {
      if(!file)
      {
	break;
      }
      
      std::istringstream sstemp(line); //Check that number of words of line <= that column + 2.
      
      std::string familyId, individualId;
      sstemp >> familyId >> individualId;
      std::string key = familyId + "@" + individualId;
      if( familyId == "FID" && individualId == "IID" && idx == 0 )
      {
        misc.message << "First row starting with FID and IID. Assuming the first row is a header." << std::endl;
        continue;
      }
      if( options.keysIndividualsToKeep.size() != 0 && options.keysIndividualsToKeep.find(key) == options.keysIndividualsToKeep.end() )
      {
        nRemovedPhenotypesFromOptions++;
        continue;
      }
      
      std::string dummy;
      for(int i=0; i<(column-1); i++)
      {
	sstemp >> dummy;
      }
      
      double phenotype;
      std::string sphenotype;
      if((sstemp >> sphenotype).fail())
      {
        misc.error("Error: The phenotype file is not properly formated for the individual: " + familyId + " " + individualId, 0);
      }
      
      if(sphenotype != "NA")
      {
	sstemp.clear();
	sstemp.str(sphenotype);
	if( (sstemp >> phenotype).fail() )
	{
	  misc.error("Error: The phenotype of the individual " + familyId + " " + individualId + " has not a valid value: " + sphenotype, 0);
	}
	phenosArrayTemp.push_back(phenotype);
	
	this->individualIds.push_back(key);
	if(this->individualIdsIdx.count(key) != 0)
	{
	  misc.error("Error: The individual with family Id: " + familyId + " and individual Id: " + individualId + " appears more than one time in file [ " + f + " ].", 0);
	}
	this->individualIdsIdx[key] = idx;
	idx++;
      }
      else
      {
	this->missings.push_back(key);
      }
    }
    file.close();
    
    this->nIndividuals = phenosArrayTemp.size();

    misc.message << phenosArrayTemp.size() << " phenotypes have been read from column " << column << " of file [ " + f + " ]. " << this->missings.size() << " missings found." << std::endl;
    if( nRemovedPhenotypesFromOptions != 0 )
    {
      misc.message << nRemovedPhenotypesFromOptions << " individuals specified in file [ " + options.fileIndividualsToKeep + " ] have been removed." << std::endl;
    }
  }
  
  communicator->broadcast(&this->nIndividuals);
  
  if(this->phenotypes!=NULL)
  {
    delete this->phenotypes;
    this->phenotypes = NULL;
  }
  this->phenotypes = new Matrix(dist, this->nIndividuals, 1, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->phenotypes->scatterMatrix(&phenosArrayTemp[0]);
}

void Phenotype::filterIndividuals(std::vector<std::string> keepIndividualIds)
{
  double * originalGlobMatrix;
  double * filteredGlobMatrix;
  if(communicator->mpiRoot)
  {
    originalGlobMatrix = new double [this->phenotypes->nGlobRows*this->phenotypes->nGlobCols];
    filteredGlobMatrix = new double [keepIndividualIds.size()*1];
  }
  
  this->phenotypes->gatherMatrix(originalGlobMatrix);
  
  if(communicator->mpiRoot)
  {
    std::map<std::string, int> newIndividualIdsIdx;
    
    for(int r=0; r<keepIndividualIds.size(); r++) //Is it better access rows contiguously and search in map each step? Or inverse?
    {
      if(this->individualIdsIdx.count(keepIndividualIds[r]) == 0)
      {
	misc.error("Error: An internal error was happened. The individual '" + keepIndividualIds[r]  + "' is not in the phenotypes matrix.", 0);
      }
      if(newIndividualIdsIdx.count(keepIndividualIds[r]) != 0)
      {
	misc.error("Error: An internal error was happened. The individual '" + keepIndividualIds[r]  + "' is repeated.", 0);
      }
      int oldRowIdx = this->individualIdsIdx[keepIndividualIds[r]];
      for(int c=0; c<1; c++)
      {
	filteredGlobMatrix[c*keepIndividualIds.size() + r] = originalGlobMatrix[c*this->nIndividuals + oldRowIdx];
      }
      newIndividualIdsIdx[keepIndividualIds[r]] = r;
    }
    this->nIndividuals = keepIndividualIds.size();
    this->individualIds = keepIndividualIds;
    this->individualIdsIdx = newIndividualIdsIdx;
  }
  communicator->broadcast(&this->nIndividuals, 1);
  communicator->barrier();
  
  this->phenotypes->initParameters(this->nIndividuals, 1, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->phenotypes->scatterMatrix(filteredGlobMatrix);
  
  if(communicator->mpiRoot)
  {
    delete [] originalGlobMatrix;
    delete [] filteredGlobMatrix;
  }
}

double Phenotype::computePhenotypeVariance()
{
  double *globalPhenotype;
  double variance;
  
  if(this->nIndividuals <= 1)
  {
    misc.error("Error: Unable to compute the phenotypic variance of less than two individuals.", 0);
  }
  
  if (communicator->mpiRoot) {
    globalPhenotype = new double [this->nIndividuals];
  }
  
  this->phenotypes->gatherMatrix(globalPhenotype);
  
  if (communicator->mpiRoot)
  {
    double avg = 0.;
    for(int r = 0; r<this->nIndividuals; r++)
    {
      avg += globalPhenotype[r];
    }
    avg /= double(this->nIndividuals);
    variance = 0.;
    for(int r = 0; r<this->nIndividuals; r++)
    {
      double temp = globalPhenotype[r] - avg;
      variance += temp * temp;
    }
    variance /= (double(this->nIndividuals) - 1.);
  }
  
  communicator->broadcast(&variance, 1);
  
  if (communicator->mpiRoot) {
    delete [] globalPhenotype;
  }
  
  return variance;
}

void Phenotype::replaceIndividualIds(std::vector<std::string> & newIndividualIds)
{
  this->individualIds.clear();
  this->individualIdsIdx.clear();
  
  if(communicator->mpiRoot)
  {
    if(this->nIndividuals != newIndividualIds.size())
    {
      misc.error("Error: An internal error was happened. Individual Ids cannot be replaced. Size of new individual Ids differ from number of current phenotypes.", 0);
    }
      
    this->individualIds = newIndividualIds;
      
    for(int i = 0; i<this->nIndividuals; i++)
    {
      std::string key = newIndividualIds[i];
      if(this->individualIdsIdx.count(key) != 0)
      {
        misc.error("Error: The individual " + key + " appears more than one time in the Phenotypes replacement keys.", 0);
      }
      this->individualIdsIdx[key] = i;
    }
  }
}

void Phenotype::printPhenotype()
{
  double *p;
  
  if (communicator->mpiRoot) {
    p = new double [this->nIndividuals];
  }
  
  this->phenotypes->gatherMatrix(p);
  if (communicator->mpiRoot) {
    misc.message << "Phenotype Matrix:\n";
    misc.message << this->nIndividuals << " " << this->individualIds.size() << " " << this->individualIdsIdx.size() << std::endl;
    //int nextMissing = 0;
    for (int r = 0; r < this->individualIds.size(); r++)
    {
      misc.message << this->individualIds[r] << " " << this->individualIdsIdx[this->individualIds[r]] << ": " ;
      /*if(nextMissing < this->missings.size())
      {
	if(this->missings[nextMissing] == r)
	{
	  misc.message << "NA" << std::endl;
	  nextMissing++;
	  continue;
	}
      }*/
      misc.message << p[r] << std::endl;
    }
    misc.message << std::endl;
    for (int r = 0; r < this->missings.size(); r++)
    {
      misc.message << this->missings[r] << " NA" << std::endl;
    }
  }
  
  if (communicator->mpiRoot) {
    delete [] p;
  }
}
