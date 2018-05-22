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

#include "misc.h"
#include "global.h"
#include "options.h"
#include "communicator.h"
#include "message.h"
#include "hetvector.h"

#include <string>

HetVector::HetVector()
{
}

HetVector::~HetVector()
{
  for(std::map<std::string, Genotype *>::iterator it = this->genotypes.begin(); it != this->genotypes.end(); ++it)
  {
    delete it->second;
  }
  for(std::map<std::string, Covariate *>::iterator it = this->covariates.begin(); it != this->covariates.end(); ++it)
  {
    delete it->second;
  }
  genotypes.clear();
  covariates.clear();
}


void HetVector::addElement(std::string name, Genotype * genotype)
{
  if( this->elements.find(name) != this->elements.end() )
  {
    misc.error("Error: An error was happened. Repeated name on HetVector.", 0);
  }
  
  this->elements.insert(name);
  
  this->genotypes[ name ] = genotype;
}

void HetVector::addElement(std::string name, Covariate * covariate)
{
  if( this->elements.find(name) != this->elements.end() )
  {
    misc.error("Error: An error was happened. Repeated name on HetVector.", 0);
  }
  
  this->elements.insert(name);
  
  this->covariates[ name ] = covariate;
}

void HetVector::deleteElement(std::string name)
{
  if( this->genotypes.count(name) != 0 )
  {
    delete this->genotypes[name];
    this->genotypes.erase(name);
  }
  else if( this->covariates.count(name) != 0 )
  {
    delete this->covariates[name];
    this->covariates.erase(name);
  }
  else
  {
    misc.error("Error: An error was happened. Element name does not exist in HetVector. It cannot be deleted.", 0);
  }
}

Matrix * HetVector::getElementMatrix(std::string name, std::vector<std::string> individualIds)
{
  Matrix * mResult;
  if( this->genotypes.count(name) != 0 )
  {
    Genotype * tempGenotype;
    this->genotypes[name]->normalizeGenotypes();
    //****************************************************
    misc.exit("Error: This function has to be finished. Take into account that the order of individualIds can be changed by the order present in the genotypes, if it is not the same.")
    //****************************************************
    this->genotypes[name]->filterSNPsAndIndividuals(this->genotypes[name]->SNPIds, individualIds, false, tempGenotype);
    mResult = tempGenotype->genotypes;
    tempGenotype->genotypes = NULL;
    delete tempGenotype;
  }
  else if( this->covariates.count(name) != 0 )
  {
    misc.error("Error: An error was happened. This function must still be implemented.", 0);
  }
  else
  {
    misc.error("Error: An error was happened. Inexistent element in HetVector.", 0);
  }
  return mResult;
}

