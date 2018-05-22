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

#include "predictphenotype.h"
#include "genotype.h"
#include "matrix.h"
#include "communicator.h"
#include "misc.h"
#include "options.h"
#include "auxiliar.h"
#include "covariate.h"

#include <vector>
#include <map>
#include <cmath>
#include <sstream>
#include <algorithm>

PredictPhenotype::PredictPhenotype(Genotype * srcGenotype, std::string fname, bool errorWhenNoSNPEffects)
{
  this->matrixEffects = NULL;
  this->matrixShift = NULL;
  
  this->genotype = srcGenotype;
  this->individualIds = this->genotype->individualIds;
  this->individuals = this->genotype->individuals;

  
  this->effectSNPsEmpty = true;
  
  readSNPsEffects(fname, errorWhenNoSNPEffects);
}

PredictPhenotype::~PredictPhenotype()
{
  if(this->genotype != NULL)
  {
    delete this->genotype;
  }
  if(this->matrixEffects != NULL)
  {
    delete this->matrixEffects;
  }
  if(this->matrixShift != NULL)
  {
    delete this->matrixShift;
  }
}

void PredictPhenotype::readSNPsEffects(std::string f, bool errorWhenNoSNPEffects)
{
  std::ifstream file;
  std::string line;
  
  this->SNPEffects.clear();
  this->effectSNPIds.clear();
  
  std::vector<std::string> tempSNPIds;
  if(communicator->mpiRoot)
  {
    misc.message << "Reading SNP effects from file [ " << f << " ] ..." << std::endl;
    
    misc.checkFileExists(f);
    file.open(f.c_str());
    
    getline(file,line);
    if( misc.checkFileHeader(line, "SNP ALLELE BLUP STDEV MEAN NBLUP") == true ) //Load predictions from REML analysis
    {
      misc.message << "REML effects detected." << std::endl;
      while(getline(file,line))
      {
        if(!file)
        {
          break;
        }
        SNPEffect effect = loadREMLEffect(f, line);
        this->SNPEffects[effect.name] = effect;
        tempSNPIds.push_back(effect.name);
      }
    }
    else if( misc.checkFileHeader(line, "GROUP SNP ALLELE MEAN STDEV BETA NBETA SE PV GROUPPV") == true ) //Load predictions from GWAS analysis
    {
      misc.message << "GWAS effects detected." << std::endl;
      while(getline(file,line))
      {
        if(!file)
        {
          break;
        }
        SNPEffect effect = loadGWASEffect(f, line);
        this->SNPEffects[effect.name] = effect;
        tempSNPIds.push_back(effect.name);
      }
    }
    else if( misc.checkFileHeader(line, "SNP ALLELE BLUP STDEV MEAN") == true ) //Load predictions from REML analysis (old file version)
    {
      misc.message << "REML effects detected." << std::endl;
      while(getline(file,line))
      {
        if(!file)
        {
          break;
        }
        SNPEffect effect = loadOldREMLEffect(f, line);
        this->SNPEffects[effect.name] = effect;
        tempSNPIds.push_back(effect.name);
      }
    }
    else
    {
      misc.error("Error: The file [ " + f + " ] has a bad header.", 0);
    }
    
    file.close();
    
    if(this->SNPEffects.size() == 0)
    {
      misc.error("Error: The file [ " + f + " ] is empty.", 0);
    }
    
    misc.message << this->SNPEffects.size() << " SNP effects have been read properly." << std::endl;
  }
  
  //Order effectSNPIds as in genotype.
  this->effectSNPsEmpty = false;
  if(communicator->mpiRoot)
  {
    this->effectSNPIds = orderVectorAsTemplate(this->genotype->SNPIds, tempSNPIds);
    if(this->effectSNPIds.size() == 0)
    {
      this->effectSNPsEmpty = true;
    }
  }
  communicator->broadcast(&this->effectSNPsEmpty);
  if( this->effectSNPsEmpty == true )
  {
    if ( errorWhenNoSNPEffects == true)
    {
      misc.error("Error: There is no overlap between SNPs in effects file and SNPs in genotypes file.", 0);
    }
    else
    {
      misc.message << "Warning: There is no overlap between SNPs in effects file and SNPs in genotypes file." << std::endl;
      return;
    }
  }
  misc.message << this->effectSNPIds.size() << " SNPs will be used for predicting phenotypes." << std::endl;
  
  this->genotype->filterSNPsAndIndividuals(this->effectSNPIds, this->genotype->individualIds);

  //Correct for reference allele.
  misc.message << "Correcting for reference allele." << std::endl;  
  if(communicator->mpiRoot)
  {
    if(this->effectSNPIds.size() != genotype->SNPs.size())
    {
      misc.error("Error: An internal error was happened when phenotype computation.", 0);
    }
    
    for(int i = 0; i < this->effectSNPIds.size(); i++)
    {
      std::string name = this->effectSNPIds[i];
      if(this->SNPEffects[name].name != genotype->SNPs[i].name)
      {
        misc.error("Error: An internal error was happened when correcting for reference allele in predicting phenotype computation.", 0);
      }
      
      if(this->SNPEffects[name].refAllele == genotype->SNPs[i].allele2)
      {
        //Nothing to do
      }
      else if(this->SNPEffects[name].refAllele == genotype->SNPs[i].allele1)
      {
        this->SNPEffects[name].mean = 2. - this->SNPEffects[name].mean;
        this->SNPEffects[name].effect = -this->SNPEffects[name].effect;
        this->SNPEffects[name].neffect = -this->SNPEffects[name].neffect;
      }
      else
      {
        misc.error("Error: The SNP " + name + " have different allele deffinitions on genotypes file and SNP effects file.", 0);
      }
    }
  }
  
  //Populate slope and shift matrices
  double * globalMatrixEffects;
  double * globalMatrixShift;
  if(communicator->mpiRoot)
  {
    globalMatrixEffects = new double [this->effectSNPIds.size()];
    globalMatrixShift = new double [this->effectSNPIds.size()];
    for(int i = 0; i<this->effectSNPIds.size(); i++)
    {
      std::string name = this->effectSNPIds[i];
      globalMatrixEffects[i] = this->SNPEffects[ name ].neffect;
      globalMatrixShift[i] = -this->SNPEffects[ name ].mean*this->SNPEffects[ name ].neffect;
    }
  }
  this->matrixEffects = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, genotype->nSNPs, 1, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->matrixShift = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, genotype->nSNPs, 1, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->matrixEffects->scatterMatrix(&globalMatrixEffects[0]);
  this->matrixShift->scatterMatrix(&globalMatrixShift[0]);
  
  if(communicator->mpiRoot)
  {
    delete [] globalMatrixEffects;
    delete [] globalMatrixShift;
  }
}

SNPEffect PredictPhenotype::loadREMLEffect(std::string f, std::string line)
{
  std::istringstream sstemp(line);
      
  SNPEffect effect;

  sstemp >> effect.name;
  sstemp >> effect.refAllele;
  
  if( (sstemp >> effect.effect).fail() )
  {
    misc.error("Error: The SNP " + effect.name + " has not a valid number in column 3 of file [ " + f + " ].", 0);
  }
  if( (sstemp >> effect.stdev).fail() )
  {
    misc.error("Error: The SNP " + effect.name + " has not a valid number in column 4 of file [ " + f + " ].", 0);
  }
  if( (sstemp >> effect.mean).fail() )
  {
    misc.error("Error: The SNP " + effect.name + " has not a valid number in column 5 of file [ " + f + " ].", 0);
  }
  if( (sstemp >> effect.neffect).fail() )
  {
    misc.error("Error: The SNP " + effect.name + " has not a valid number in column 6 of file [ " + f + " ].", 0);
  }
  
  if(this->SNPEffects.count(effect.name) != 0)
  {
    misc.error("Error: The SNP " + effect.name + " appears repeated in file [ " + f + " ].", 0);
  }
  
  return effect;
}

SNPEffect PredictPhenotype::loadOldREMLEffect(std::string f, std::string line)
{
  std::istringstream sstemp(line);
      
  SNPEffect effect;

  sstemp >> effect.name;
  sstemp >> effect.refAllele;
  
  if( (sstemp >> effect.effect).fail() )
  {
    misc.error("Error: The SNP " + effect.name + " has not a valid number in column 3 of file [ " + f + " ].", 0);
  }
  if( (sstemp >> effect.stdev).fail() )
  {
    misc.error("Error: The SNP " + effect.name + " has not a valid number in column 4 of file [ " + f + " ].", 0);
  }
  if( (sstemp >> effect.mean).fail() )
  {
    misc.error("Error: The SNP " + effect.name + " has not a valid number in column 5 of file [ " + f + " ].", 0);
  }
  effect.neffect = effect.effect/effect.stdev;
  
  if(this->SNPEffects.count(effect.name) != 0)
  {
    misc.error("Error: The SNP " + effect.name + " appears repeated in file [ " + f + " ].", 0);
  }
  
  return effect;
}

SNPEffect PredictPhenotype::loadGWASEffect(std::string f, std::string line)
{
  std::istringstream sstemp(line);
  std::string dummy;
      
  SNPEffect effect;

  sstemp >> dummy; //Get the group
  
  sstemp >> effect.name;
  sstemp >> effect.refAllele;
  
  if( (sstemp >> effect.mean).fail() )
  {
    misc.error("Error: The SNP " + effect.name + " has not a valid number in column 4 of file [ " + f + " ].", 0);
  }
  if( (sstemp >> effect.stdev).fail() )
  {
    misc.error("Error: The SNP " + effect.name + " has not a valid number in column 5 of file [ " + f + " ].", 0);
  }
  if( (sstemp >> effect.effect).fail() )
  {
    misc.error("Error: The SNP " + effect.name + " has not a valid number in column 6 of file [ " + f + " ].", 0);
  }
  if( (sstemp >> effect.neffect).fail() )
  {
    misc.error("Error: The SNP " + effect.name + " has not a valid number in column 7 of file [ " + f + " ].", 0);
  }
  
  if(this->SNPEffects.count(effect.name) != 0)
  {
    misc.error("Error: The SNP " + effect.name + " appears repeated in file [ " + f + " ].", 0);
  }
  
  return effect;
}

void PredictPhenotype::predictPhenotypes()
{
  if( this->genotype == NULL )
  {
    misc.error("Error: An internal error was happened when trying to predict the phenotypes.", 0);
  }
  if ( this->effectSNPsEmpty == true )
  {
    return;
  }
  
  Matrix * phenotypes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Matrix * phenotypeShift = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  
  Matrix * correction = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  
  phenotypes->multiply(this->genotype->genotypes, 'T', this->matrixEffects, 'N');
  phenotypeShift->multiply(this->genotype->missings, 'T', this->matrixShift, 'N');
  phenotypes->add(phenotypeShift, 1., 1.);
  correction->multiply(this->genotype->missings, 'T', this->matrixEffects, 'N', -1.); //Correction because our genotypes are coded 1, 2, and 3 instead of 0, 1, 2.
  phenotypes->add(correction, 1., 1.);
  
  phenotypes->matrixToStandardVector(this->globalPhenotypes);
  phenotypeShift->matrixToStandardVector(this->globalPhenotypeShift);
  
  delete this->genotype;
  this->genotype = NULL;
  
  delete phenotypes;
  delete phenotypeShift;

  delete correction;
}

void PredictPhenotype::storePredictions()
{
  if ( this->effectSNPsEmpty == true )
  {
    misc.error("Error: There is no overlap between SNPs in effects file(s) and SNPs in genotypes file(s)", 0);
  }
  
  Covariate covariate(options.covarsFile, options.qCovarsFile, std::vector<std::string>(), false);
  std::map<std::string, double> covariateEffects = covariate.loadEffectPrediction(options.fileCovarEffects, options.fileQCovarEffects);
  
  if(communicator->mpiRoot)
  {
    Message message(options.outFile + ".predicted.phenos");
    for(int i = 0; i < this->individuals.size(); i++)
    {
      message << this->individuals[i].familyID << " " << this->individuals[i].individualID;
      message << " " << this->globalPhenotypes[i]; //Genotype computed prediction
      if( covariateEffects.count(individualIds[i]) != 0 )
      {
        message << " " << covariateEffects[individualIds[i]]; //Covariate computed prediction
        message << " " << this->globalPhenotypes[i] + covariateEffects[individualIds[i]]; //Total prediction
      }
      else
      {
        message << " NA"; //Covariate prediction is not computed
        message << " NA"; //Total prediction is not computed
      }
      message << " " << this->globalPhenotypeShift[i];
      message << " " << this->globalPhenotypes[i] - this->globalPhenotypeShift[i];
      message << std::endl;
    }
  }
}

void PredictPhenotype::addMoreEffects(PredictPhenotype *previousPredictions)
{
  //Same individuals in two predictions?
  if(previousPredictions->individualIds != this->individualIds || previousPredictions->individuals != this->individuals)
  {
    misc.error("Error: When computing predictions using multiple files. All files must contain the same exact individuals.", 0);
  }
  
  //Is previousPredictions empty?
  if(previousPredictions->effectSNPsEmpty == true)
  {
    misc.message << "No effect SNPs, skipping effect addition." << std::endl;
    return;
  }
  
  //Is this prediction empty? If so, copy previousPrediction and return.
  if(this->effectSNPsEmpty == true)
  {
    if(communicator->mpiRoot)
    {
      for(int i = 0; i < previousPredictions->globalPhenotypes.size(); i++)
      {
        this->globalPhenotypes.push_back(previousPredictions->globalPhenotypes[i]);
        this->globalPhenotypeShift.push_back(previousPredictions->globalPhenotypeShift[i]);
      }
    }
    this->effectSNPIds = previousPredictions->effectSNPIds;
    this->effectSNPsEmpty = false;
    misc.message << "Effects added." << std::endl;
    return;
  }
  
  //Sizes on previous predictions agree with this?
  if(
    previousPredictions->globalPhenotypes.size() != this->globalPhenotypes.size() ||
    previousPredictions->globalPhenotypeShift.size() != this->globalPhenotypeShift.size() ||
    previousPredictions->globalPhenotypes.size() != previousPredictions->globalPhenotypeShift.size()
  )
  {
    misc.error("Error: An internal error was happened when joining two predictions.", 0);
  }
  
  //Are there repeated markers between two predictions?
  std::vector<std::string> SNPIntersection = intersectionStringVectors(2, &this->effectSNPIds, &previousPredictions->effectSNPIds);
  if(SNPIntersection.size() != 0)
  {
    misc.error("Error: When computing predictions using multiple files. There are at least one SNPs repeated in more than one file. This could account for the same effect more than one time.", 0);
  }
  
  //Add two predictions
  if(communicator->mpiRoot)
  {
    for(int i = 0; i < this->globalPhenotypes.size(); i++)
    {
      this->globalPhenotypes[i] += previousPredictions->globalPhenotypes[i];
      this->globalPhenotypeShift[i] += previousPredictions->globalPhenotypeShift[i];
    }
  }
  this->effectSNPIds.insert( this->effectSNPIds.end(), previousPredictions->effectSNPIds.begin(), previousPredictions->effectSNPIds.end() );
  misc.message << "Effects added." << std::endl;
}

