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

#include "multireml.h"
#include "genotype.h"
#include "kernel.h"
#include "covariancematrix.h"
#include "reml.h"
#include "covariate.h"
#include "phenotype.h"
#include "matrix.h"
#include "options.h"
#include "auxiliar.h"
#include "message.h"
#include "misc.h"
#include "global.h"
#include "results.h"
#include "singlereml.h"

#include <cmath>
#include <set>

MultiREML::MultiREML(bool argWriteResults)
{
  this->writeResults = argWriteResults;
  
  this->reml = NULL;
}

MultiREML::~MultiREML()
{
  if(this->reml != NULL)
  {
    delete this->reml;
    this->reml = NULL;
  }
}

void MultiREML::compute()
{
  if(this->reml != NULL)
  {
    delete this->reml;
  }
  
  Genotype *genotypes = NULL;
  Kernel *grm = loadGRMUsingOptions(true, &genotypes);
  if(options.computeSNPsBLUP)
  {
    if( genotypes == NULL && options.SNPBLUPGenotypeFiles.size() == 0 )
    {
      genotypes = loadGenotypeUsingOptions();
    }
  }
  else
  {
    if(genotypes != NULL)
    {
      delete genotypes;
      genotypes = NULL;
    }
  }
  
  std::vector<std::string> reducedModels;
  std::vector<std::string> individualBLUPNames;
  std::vector<Kernel*> kernels;
  grm->name = options.baseVarianceNames["GRM"];
  kernels.push_back(grm);
  reducedModels.push_back(grm->name);
  individualBLUPNames.push_back(grm->name);
  introduceResortedGRMsByCouples(kernels, reducedModels, individualBLUPNames);
  addKernelsUsingOptions(kernels, reducedModels, individualBLUPNames);
  
  std::vector<double> weights;
  weights.assign( kernels.size() , 1./double(kernels.size()) );
  
  std::vector<std::pair<std::string, std::string> > covariateFiles;
  if(options.covarsFiles.size() != options.qCovarsFiles.size())
  {
    misc.error("Error: An internal error was happened. Inconsistent number of covars and qcovars.", 0);
  }
  for(int i = 0; i<options.covarsFiles.size(); i++)
  {
    covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFiles[i], options.qCovarsFiles[i]));
  }
  
  this->reml = new REML(this->writeResults);
  
  if(genotypes != NULL && options.computeSNPsBLUP)
  {
    this->reml->mSNPIds[grm->name] = grm->randomVarNames;
    this->reml->SNPsBLUPGenotypes[grm->name] = genotypes;
  }
  else if(options.computeSNPsBLUP)
  {
    this->reml->mSNPIds[grm->name] = grm->randomVarNames;
    this->reml->SNPsBLUPGenotypeFiles[grm->name] = options.SNPBLUPGenotypeFiles;
  }
  if(options.computeIndividualsBLUP == true)
  {
    this->reml->individualBLUPNames = individualBLUPNames;
  }
  
  this->reml->elementsToTest = reducedModels;
  
  bool prepared = this->reml->prepare(bivariateREMLType, kernels, weights, options.phenotypeColumns, options.initialh2Traits, covariateFiles);
  
  if( prepared == true )
  {
    this->reml->computeREMLWithReducedModels();
  }
  else
  {
    misc.message << "Sorry, a problem was happened while preparing data for performing REML. The MLM cannot be fitted. Please, check the logs." << std::endl;
  }
  
  delete this->reml;
  this->reml = NULL;
}

void MultiREML::computeRegional()
{
  if(this->reml != NULL)
  {
    delete this->reml;
  }
  
  //Load genotypes and GRM
  Genotype * genotype = loadGenotypeUsingOptions();
  Kernel * grm;
  if(options.grmFile == "")
  {
    grm = new Kernel(genotype);
  }
  else
  {
    grm = new Kernel(options.grmFile);
    genotype->normalizeGenotypes();
  }
  
  if(grm->individualIds != genotype->individualIds)
  {
    misc.error("Error: The individuals in GRM and genotypes are not the same. Maybe you are not using the proper GRM or genotype file?", 0);
  }

  if( grm->diagonalized == true )
  {
    if(options.forceUseDiagonalizedKernels == true)
    {
      grm->recoverKernelFromEigenDecomposition();
    }
    else
    {
      misc.error("Error: Sorry, this analysis cannot be performed with diagonal GRMs. You can use the option --force-use-diag-grms to force converting diagonal GRMs to their non-diagonalized form before starting the analysis.", 0);
    }
  }
  
  //Create groups
  genotype->groupSNPs(options.regionBy);
  
  //Start iterations over groups
  std::string baseOutFile = options.outFile;
  int nGroups = genotype->groupedSNPs.size();
  communicator->broadcast(&nGroups);
  std::map<std::string, std::set<std::string> >::iterator it = genotype->groupedSNPs.begin();
  for(int i = 0; i < nGroups; i++)
  {
    std::string group = "";
    if(communicator->mpiRoot)
    {
      group = it->first;
      it++;
    }
    communicator->broadcast(group);
    
    misc.message << "\nAnalysing region " << group << "..." << std::endl;
    
    options.outFile = baseOutFile + "." + group;
    
    this->reml = new REML();
    
    Genotype *regionalGenotype = new Genotype();
    genotype->genotypeOfSNPsGroup(group, regionalGenotype);
    double proportionRegionalTotal = double(regionalGenotype->nSNPs)/double(genotype->nSNPs);
    Kernel *regionalGRM = new Kernel(regionalGenotype);
    Kernel *globalGRM = new Kernel();
    globalGRM->addKernels(1., grm, -1., regionalGRM);
    delete regionalGenotype;
    globalGRM->name = options.baseVarianceNames["Global-GRM"];
    regionalGRM->name = options.baseVarianceNames["Regional-GRM"];
    
    std::vector<std::string> snpIntersection = orderVectorAsTemplate(regionalGRM->randomVarNames, grm->randomVarNames);
    if( snpIntersection != regionalGRM->randomVarNames )
    {
      misc.error("Error: SNPs in regional GRM are not present in global GRM. This could be due that some SNPs in genotype file are not present in the the GRM (i.e. they are not used for computing the global GRM).", 0);
    }
    
    std::vector<Kernel*> kernels;
    kernels.push_back(globalGRM);
    kernels.push_back(regionalGRM);
    
    std::vector<double> weights;
    weights.push_back(1. - proportionRegionalTotal);
    weights.push_back(proportionRegionalTotal);
    
    std::vector<std::pair<std::string, std::string> > covariateFiles;
    if(options.covarsFiles.size() != options.qCovarsFiles.size())
    {
      misc.error("Error: An internal error was happened. Inconsistent number of covars and qcovars.", 0);
    }
    for(int i = 0; i<options.covarsFiles.size(); i++)
    {
      covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFiles[i], options.qCovarsFiles[i]));
    }
    
    if(options.computeIndividualsBLUP == true)
    {
      this->reml->individualBLUPNames.push_back(globalGRM->name);
      this->reml->individualBLUPNames.push_back(regionalGRM->name);
    }
    
    bool prepared = this->reml->prepare(multipleBivariateREMLType, kernels, weights, options.phenotypeColumns, options.initialh2Traits, covariateFiles);
    
    if( prepared == true )
    {
      this->reml->computeREML();
    }
    else
    {
      misc.message << "WARNING: Sorry, a problem was happened while preparing data for performing REML on region " << group << ". The MLM cannot be fitted on this region. Please, check the logs." << std::endl;
    }
    
    delete this->reml;
    this->reml = NULL;
  }
  
  delete grm;
  delete genotype;
}

void MultiREML::computeMultipleGroups()
{
  SingleREML singleREML;
  singleREML.computeMultipleGroups();
}


void MultiREML::compareREMLs(double baseLogLikelihood, double baseNVariances)
{
  if(communicator->mpiRoot)
  {
    double LogRatio = 2.0*(baseLogLikelihood - this->reml->logLikelihood);
    if(LogRatio < 0.)
    {
      LogRatio = 0.;
    }
    int df = baseNVariances - this->reml->V->variances.size();
    
    Message message(options.outFile + ".reml.fixed");
    
    message << "LRT\t" << std::setprecision(3) << LogRatio << std::endl;
    message << "df\t" << std::setprecision(1) << df << std::endl;
    message << "Pval\t" << std::setprecision(4) << 0.5*chi1_CDF(df, LogRatio) << std::endl;
  }
}