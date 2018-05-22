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

#include "simulatephenotype.h"
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

SimulatePhenotype::SimulatePhenotype(Genotype * srcGenotype, std::string fname, Genotype * argAdjustEffectsGenotype)
{
  this->genotype = srcGenotype;
  this->effectSizes = NULL;
  this->seed = options.randomSeed;
  
  this->adjustEffectsGenotype = argAdjustEffectsGenotype;
  
  readCausalSNPs(fname);
  
  this->genotype->normalizeGenotypes();
  this->genotype->filterSNPsAndIndividuals(this->causalSNPIds, this->genotype->individualIds, true);
  
}

SimulatePhenotype::~SimulatePhenotype()
{
  if(this->genotype != NULL)
  {
    delete this->genotype;
  }
  if(this->adjustEffectsGenotype != NULL)
  {
    delete this->adjustEffectsGenotype;
  }
  if(this->effectSizes != NULL)
  {
    delete this->effectSizes;
  }
}

void SimulatePhenotype::readCausalSNPs(std::string f)
{
  std::ifstream file;
  std::string line;
  std::vector<std::string> unorderedCausalSNPIds;
  
  if(communicator->mpiRoot)
  {
    int nRandom = 0;
    
    misc.message << "Reading causal loci effect sizes from file [ " << f << " ] ..." << std::endl;
    
    misc.checkFileExists(f);
    file.open(f.c_str());
    Message message(options.outFile + ".simulated.effects");
    
    unorderedCausalSNPIds.clear();
    this->causalSNPEffects.clear();
    int idx = 0;
    while(getline(file,line))
    {
      if(!file)
      {
        break;
      }
      
      std::istringstream sstemp(line); //Check that number of words of line <= that column + 2.
      
      std::string snp;
      
      sstemp >> snp;
      
      if( this->genotype->SNPIdsIdx.count(snp) == 0 )
      {
        misc.message << "Warning: SNP " << snp << " is in effects file but not in genotypes file." << std::endl;
        continue;
      }
      if( this->adjustEffectsGenotype != NULL )
      {
        if( this->adjustEffectsGenotype->SNPIdsIdx.count(snp) == 0 )
        {
          misc.message << "Warning: SNP " << snp << " is in effects file but not in genotypes file for adjusting SNP effects." << std::endl;
          continue;
        }
      }
      
      std::string sEffect;
      double effect;
      
      bool flag = (sstemp >> sEffect).fail();
      if( flag == true )
      {
        effect = box_muller(0., 1., &this->seed);
        nRandom++;
      }
      else
      {
        sstemp.clear();
        sstemp.str(sEffect);
        if( (sstemp >> effect).fail() )
        {          
          misc.error("Error: The effect size of the SNP " + snp + " has not a valid value.", 0);
        }
        
        //Correct the effect using other genotype file frequencies if specified
        if( this->adjustEffectsGenotype != NULL )
        {
          int SrcSNPIdx = this->genotype->SNPIdsIdx[snp];
          int AdjustSNPIdx = this->adjustEffectsGenotype->SNPIdsIdx[snp];
          effect = effect * (this->genotype->SNPs[ SrcSNPIdx ].standardDev / this->adjustEffectsGenotype->SNPs[ AdjustSNPIdx ].standardDev );
        }
      }
      
      if(this->causalSNPEffects.count(snp) != 0)
      {
        misc.error("Error: The SNP " + snp + " appears repeated in file [ " + f + " ].", 0);
      }
      unorderedCausalSNPIds.push_back(snp);
      this->causalSNPEffects[snp] = effect;
      message << snp << " " << effect << std::endl;
      idx++;
    }
    file.close();
    
    if(unorderedCausalSNPIds.size() == 0)
    {
      misc.error("Error: The file [ " + f + " ] is empty or there is not overlapping with genotype file(s).", 0);
    }
    
    misc.message << "Effect sizes have been read properly:" << std::endl;
    misc.message << "  " << unorderedCausalSNPIds.size() - nRandom << " fixed effect sizes." << std::endl;
    misc.message << "  " << nRandom << " random generated effect sizes." << std::endl;
    
  }
  
  if(communicator->mpiRoot) //Order causalSNPIds that appear in genotypes.
  {
    this->causalSNPIds = orderVectorAsTemplate(this->genotype->SNPIds, unorderedCausalSNPIds);
    if(this->causalSNPIds.size() != unorderedCausalSNPIds.size())
    {
      misc.error("Error: An internal error was happened. There is a SNP in the causal SNP file which is not in the BIM file.", 0);
    }
  }
  
  //Populate effectSizes matrix
  double * tempEffectSizes;
  if(communicator->mpiRoot)
  {
    tempEffectSizes = new double [this->causalSNPIds.size()];
    for(int i = 0; i<this->causalSNPIds.size(); i++)
    {
      tempEffectSizes[i] = this->causalSNPEffects[ this->causalSNPIds[i] ];
    }
  }
  
  int nCausalSNPs = this->causalSNPIds.size();
  communicator->broadcast(&nCausalSNPs);
  this->effectSizes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nCausalSNPs, 1, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->effectSizes->scatterMatrix(&tempEffectSizes[0]);
  
  if(communicator->mpiRoot)
  {
    delete [] tempEffectSizes;
  }
}

void SimulatePhenotype::simulatePhenotypes()
{
  int nCases = 0;
  int nControls = 0;
  
  Matrix * geneticEffects = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  geneticEffects->multiply(this->genotype->genotypes, 'T', this->effectSizes, 'N');
  
  double geneticVariance = computeVariance(geneticEffects);
  double environmentVariance = geneticVariance*((1.0 - options.simulationh2)/options.simulationh2);
  double environmentStd = sqrt(environmentVariance);

  double * globalGeneticEffects;
  double * globalPhenotypes;
  double * globalEnvironmental;
  if(communicator->mpiRoot)
  {
    globalGeneticEffects = new double [this->genotype->nIndividuals];
    globalPhenotypes = new double [this->genotype->nIndividuals];
    globalEnvironmental = new double [this->genotype->nIndividuals];
  }
  
  geneticEffects->gatherMatrix(globalGeneticEffects);
  
  if(communicator->mpiRoot)
  {
    for(int i = 0; i < this->genotype->nIndividuals; i++)
    {
      globalPhenotypes[i] = globalGeneticEffects[i] + box_muller(0., environmentStd, &this->seed);
      globalEnvironmental[i] = globalPhenotypes[i] - globalGeneticEffects[i];
    }
    if(options.simulateBinaryTrait)
    {
      std::vector<double> sortedEffects(globalPhenotypes, globalPhenotypes + this->genotype->nIndividuals);
      std::sort(sortedEffects.begin(), sortedEffects.end());
      int nEstimatedControls = int( (this->genotype->nIndividuals*(1.0 - options.prevalence)) );
      if(nEstimatedControls == 0)
      {
        misc.error("Error: The prevalence is too high or the population too small and there are not controls. Please, increase the population size or decrease the prevalence.", 0);
      }
      double effectThreshold = 0.5*( sortedEffects[ nEstimatedControls ] + sortedEffects[ nEstimatedControls - 1 ] );
      
      for(int i = 0; i < this->genotype->nIndividuals; i++)
      {
        if(globalPhenotypes[i] > effectThreshold)
        {
          globalPhenotypes[i] = 2;
          nCases++;
        }
        else
        {
          globalPhenotypes[i] = 1;
          nControls++;
        }
      }
    }
    
    Message message(options.outFile + ".simulated.phenos");
    if( options.fileCovarEffects == "" && options.fileQCovarEffects == "")
    {
      for(int i = 0; i < this->genotype->nIndividuals; i++)
      {
        message << genotype->individuals[i].familyID << " " << genotype->individuals[i].individualID << " " << globalPhenotypes[i] << std::endl;
      }
    }
    else
    {
      Covariate covariate(options.covarsFile, options.qCovarsFile, std::vector<std::string>(), false);
      std::map<std::string, double> covariateEffects = covariate.loadEffectPrediction(options.fileCovarEffects, options.fileQCovarEffects);
      for(int i = 0; i < this->genotype->nIndividuals; i++)
      {
        if( covariateEffects.count(genotype->individualIds[i]) != 0 )
        {
          message << genotype->individuals[i].familyID << " " << genotype->individuals[i].individualID << " " << globalPhenotypes[i] << " " << covariateEffects[genotype->individualIds[i]] << " " << globalPhenotypes[i] + covariateEffects[genotype->individualIds[i]] << std::endl;
        }
      }
    }
    
    message.redirect(options.outFile + ".simulated.blups");
    for(int i = 0; i < this->genotype->nIndividuals; i++)
    {
      message << genotype->individuals[i].familyID << " " << genotype->individuals[i].individualID << " " << globalGeneticEffects[i] << " " << globalEnvironmental[i] << std::endl;
    }
  }
  
  delete geneticEffects;
  if(communicator->mpiRoot)
  {
    delete [] globalGeneticEffects;
    delete [] globalPhenotypes;
    delete [] globalEnvironmental;
  }
}
