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

#include "singlereml.h"
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
#include "glmm.h"

#include <sstream>
#include <vector>

SingleREML::SingleREML(bool argWriteResults)
{
  this->writeResults = argWriteResults;
  
  this->reml = NULL;
}

SingleREML::~SingleREML()
{
  if(this->reml != NULL)
  {
    delete this->reml;
    this->reml = NULL;
  }
}

void SingleREML::compute()
{
  if(this->reml != NULL)
  {
    delete this->reml;
  }
  
  //Prepare data
  std::map<std::string, Genotype*> genotypesList;
  std::vector<Kernel *> grmsBase;
  std::map<std::string, std::vector<std::string> > genotypesFilesList;
  std::map<std::string, std::vector<std::string> > SNPIds;
  loadGRMUsingOptions(grmsBase, true, genotypesList, genotypesFilesList, SNPIds);
  
  //precomputeInitialValues(0, grmBase);
  
  if( options.computeEpistasisVariance == true )
  {
    Kernel *epistaticGRMBase = NULL;
    std::vector<Kernel *> tempList;
    for(int gi = 0; gi<grmsBase.size(); gi++)
    {
      epistaticGRMBase = new Kernel(grmsBase[gi], kernelEpistaticGRM, NULL);
      if(grmsBase.size() == 1)
      {
        epistaticGRMBase->name = options.baseVarianceNames["epi"];
      }
      else
      {
        epistaticGRMBase->name = options.baseVarianceNames["epi"] + "-" + grmsBase[gi]->name;
      }
      tempList.push_back(epistaticGRMBase);
    }
    grmsBase.insert(grmsBase.end(), tempList.begin(), tempList.end());
  }
  
  //Add other covariances
  std::vector<std::string> reducedModels;
  std::vector<std::string> individualBLUPNames;
  std::vector<Kernel*> kernelsBase;
  for(int gi = 0; gi<grmsBase.size(); gi++)
  {
    kernelsBase.push_back(grmsBase[gi]);
    reducedModels.push_back(grmsBase[gi]->name);
    individualBLUPNames.push_back(grmsBase[gi]->name);
  }
  grmsBase.clear();
  introduceResortedGRMsByCouples(kernelsBase, reducedModels, individualBLUPNames);
  addKernelsUsingOptions(kernelsBase, reducedModels, individualBLUPNames);
  
  //Compute REML
  std::vector<int> phenotypesForAnalyze = getPhenotyesForAnalysis();
  std::vector<Kernel*> kernels;
  std::string baseOutFile = options.outFile;
  for(int fni = 0; fni < options.phenotypesFiles.size(); fni++ )
  {
    options.phenotypesFile = options.phenotypesFiles[fni];
    std::vector<std::string> header = getHeader(options.phenotypesFile);
    if(header[0] == "FID" && header[1] == "IID")
    {
      misc.message << "Using phenotype names from file [ " + options.phenotypesFile + " ] header." << std::endl;
    }
    else
    {
        header.clear();
    }
    
    for(int i = 0; i<phenotypesForAnalyze.size(); i++)
    {
      for(int ki = 0; ki<kernels.size(); ki++)
      {
        delete kernels[ki];
      }
      kernels.clear();
      if(phenotypesForAnalyze.size() != 1 || options.phenotypesFiles.size() != 1)
      {
        for(int ik = 0; ik < kernelsBase.size(); ik++)
        {
          Kernel * temp = new Kernel(kernelsBase[ik]);
          kernels.push_back(temp);
        }
        options.outFile = baseOutFile;
        if(options.phenotypesFiles.size() != 1)
        {
          options.outFile += "." + getFileName(options.phenotypesFiles[fni]);
        }
        if(phenotypesForAnalyze.size() != 1)
        {
          if(header.size() == 0)
          {
            options.outFile += "." + i2s(phenotypesForAnalyze[i]);
          }
          else
          {
            options.outFile += "." + header[phenotypesForAnalyze[i] + 1];
          }
        }
        options.phenotypeColumn = phenotypesForAnalyze[i];
      }
      else
      {
        for(int ik = 0; ik < kernelsBase.size(); ik++)
        {
          kernels.push_back(kernelsBase[ik]);
        }
        kernelsBase.clear();
      }
    
      std::vector<double> weights;
      weights.assign( kernels.size() , 1./double(kernels.size()) );
      
      std::vector<int> phenotypeColumns;
      phenotypeColumns.push_back(options.phenotypeColumn);

      std::vector<double> heritabilities;
      heritabilities.push_back(options.initialh2Trait);
      
      std::vector<std::pair<std::string, std::string> > covariateFiles;
      covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFile, options.qCovarsFile));
      
      this->reml = new REML(this->writeResults);

      if( options.computeSNPsBLUP )
      {
        this->reml->SNPsBLUPGenotypes = genotypesList;
        this->reml->SNPsBLUPGenotypeFiles = genotypesFilesList;
        this->reml->mSNPIds = SNPIds;
      }
      if(options.computeIndividualsBLUP == true)
      {
        this->reml->individualBLUPNames = individualBLUPNames;
      }
      this->reml->elementsToTest = reducedModels;
      
      bool prepared = this->reml->prepare(singleREMLType, kernels, weights, phenotypeColumns, heritabilities, covariateFiles);
      if( prepared == true )
      {
        if( options.analysis != GLMMAnalysis )
        {
          this->reml->computeREMLWithReducedModels();
        }
        else
        {
          GLMM * glmm = new GLMM(this->reml);
          glmm->fit();
        }
      }
      else
      {
        misc.message << "Sorry, a problem was happened while preparing data for performing REML. The MLM cannot be fitted. Please, check the logs." << std::endl;
      }

      this->reml->SNPsBLUPGenotypes.clear();    
      delete this->reml;
      this->reml = NULL;

    } //End of pheno file column iterations
  } //End of pheno file iterations
  
  options.outFile = baseOutFile;
  
  options.phenotypesFile = options.phenotypesFiles[0];
  options.phenotypeColumn = phenotypesForAnalyze[0];
  
  for(std::map<std::string, Genotype *>::iterator it = genotypesList.begin(); it != genotypesList.end(); ++it)
  {
    delete it->second;
  }
  for(int ki = 0; ki<kernelsBase.size(); ki++)
  {
    delete kernelsBase[ki];
  }
}

void SingleREML::computeRegional()
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
    std::vector<std::string> individualBLUPNames;
    
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
    
    std::vector<std::string> snpIntersection = orderVectorAsTemplate(regionalGRM->randomVarNames, grm->randomVarNames);
    if( snpIntersection != regionalGRM->randomVarNames )
    {
      misc.error("Error: SNPs in regional GRM are not present in global GRM. This could be due that some SNPs in genotype file are not present in the the GRM (i.e. they are not used for computing the global GRM).", 0);
    }
    
    std::vector<std::string> reducedModels;
    std::vector<Kernel*> kernels;
    globalGRM->name = options.baseVarianceNames["Global-GRM"];
    kernels.push_back(globalGRM);
    if(options.skipTestGlobalGRMRegionalAnalysis == false)
    {
      reducedModels.push_back(globalGRM->name);
    }
    individualBLUPNames.push_back(globalGRM->name);
    regionalGRM->name = options.baseVarianceNames["Regional-GRM"];
    kernels.push_back(regionalGRM);
    reducedModels.push_back(regionalGRM->name);
    individualBLUPNames.push_back(regionalGRM->name);
    addKernelsUsingOptions(kernels, reducedModels, individualBLUPNames);
    
    std::vector<double> weights;
    weights.push_back( (1. - proportionRegionalTotal)*2./double(kernels.size()) );
    weights.push_back( proportionRegionalTotal*2./double(kernels.size()) );
    for(int widx = 2; widx < kernels.size(); widx++ )
    {
      weights.push_back(1./kernels.size());
    }
    
    std::vector<int> phenotypeColumns;
    phenotypeColumns.push_back(options.phenotypeColumn);

    std::vector<double> heritabilities;
    heritabilities.push_back(options.initialh2Trait);
    
    std::vector<std::pair<std::string, std::string> > covariateFiles;
    covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFile, options.qCovarsFile));
    
    if(options.computeIndividualsBLUP == true)
    {
      this->reml->individualBLUPNames = individualBLUPNames;
    }
    
    this->reml->elementsToTest = reducedModels;
    
    bool prepared = this->reml->prepare(regionalSingleREMLType, kernels, weights, phenotypeColumns, heritabilities, covariateFiles);
    if( prepared == true )
    {
      this->reml->computeREMLWithReducedModels();
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

void SingleREML::computeMultipleGroups()
{
  if(this->reml != NULL)
  {
    delete this->reml;
  }
  
  //Load genotypes and GRM
  Genotype * genotype = loadGenotypeUsingOptions();
  genotype->normalizeGenotypes();
  
  //Create groups
  genotype->groupSNPs(options.regionBy);
  
  //Start iterations over groups
  misc.message << "Computing all regional GRMs..." << std::endl;
  misc.setGetElapsedTime("ConstructRegionalGRMs");
  
  this->reml = new REML();
  
  std::vector<Kernel*> kernels;
  std::vector<std::string> reducedModels;
  std::vector<std::string> individualBLUPNames;
  std::vector<double> weights;
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
    misc.message << "Computing group " << group << "..." << std::endl;
    
    Genotype *regionalGenotype = new Genotype();
    genotype->genotypeOfSNPsGroup(group, regionalGenotype);
    double proportionRegionalTotal = double(regionalGenotype->nSNPs)/double(genotype->nSNPs);
    Kernel *regionalGRM = new Kernel(regionalGenotype);
    regionalGRM->name = group;
    
    if(options.computeSNPsBLUP)
    {
      this->reml->mSNPIds[regionalGRM->name] = regionalGRM->randomVarNames;
      this->reml->SNPsBLUPGenotypes[regionalGRM->name] = regionalGenotype;
    }
    else
    {
      delete regionalGenotype;
    }
    
    kernels.push_back(regionalGRM);
    reducedModels.push_back(regionalGRM->name);
    individualBLUPNames.push_back(regionalGRM->name);
    weights.push_back( proportionRegionalTotal );
  }
  misc.message << "Regional GRMs computed after " << misc.setGetElapsedTime("ConstructRegionalGRMs", true) << std::endl;
  
  delete genotype;
  genotype = NULL;
  
  addKernelsUsingOptions(kernels, reducedModels, individualBLUPNames);
  for(int widx = nGroups; widx < kernels.size(); widx++ )
  {
    weights.push_back(1./double(kernels.size()));
  }
  for(int widx = 0; widx < nGroups; widx++ )
  {
    weights[widx] *= double(nGroups)/double(kernels.size());
  }
  
  
  std::vector<int> phenotypeColumns = options.phenotypeColumns;

  std::vector<double> heritabilities = options.initialh2Traits;
  
  std::vector<std::pair<std::string, std::string> > covariateFiles;
  if(options.covarsFiles.size() != options.qCovarsFiles.size())
  {
    misc.error("Error: An internal error was happened. Inconsistent number of covars and qcovars.", 0);
  }
  for(int i = 0; i<options.covarsFiles.size(); i++)
  {
    covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFiles[i], options.qCovarsFiles[i]));
  }
  
  this->reml->elementsToTest = reducedModels;
  if(options.computeIndividualsBLUP == true)
  {
    this->reml->individualBLUPNames = individualBLUPNames;
  }
  
  bool prepared = false;
  if( phenotypeColumns.size() == 1 )
  {
    prepared = this->reml->prepare(multipleSingleREMLType, kernels, weights, phenotypeColumns, heritabilities, covariateFiles);
  }
  else
  {
    prepared = this->reml->prepare(multipleBivariateREMLType, kernels, weights, phenotypeColumns, heritabilities, covariateFiles);
  }
  
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

void SingleREML::computeRaw(Kernel *srcGRM)
{
  std::vector<Kernel*> grms;
  srcGRM->name = options.baseVarianceNames["GRM"];
  grms.push_back(srcGRM);
  
  std::vector<double> weights;
  weights.push_back(1.);
  
  std::vector<int> phenotypeColumns;
  phenotypeColumns.push_back(options.phenotypeColumn);

  std::vector<double> heritabilities;
  heritabilities.push_back(options.initialh2Trait);
  
  std::vector<std::pair<std::string, std::string> > covariateFiles;
  covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFile, options.qCovarsFile));
  
  bool prepared = this->reml->prepare(singleREMLType, grms, weights, phenotypeColumns, heritabilities, covariateFiles);
  if( prepared == true )
  {
    this->reml->computeREML();
  }
  else
  {
    misc.message << "Sorry, a problem was happened while preparing data for performing REML. The MLM cannot be fitted. Please, check the logs." << std::endl;
  }
}

void SingleREML::precomputeInitialValues(int type, Kernel *srcGRM)
{
  REMLResults results;
  results.variances.clear();
  results.logLikelihood = 0;
  
  if(options.initialVariancesFile != "" )
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
        
        std::istringstream sstemp(line);
        
        Variance variance;
        
        sstemp >> variance.name;
        if( (sstemp >> variance.variance).fail() )
        {
          misc.error("Error: The variance " + variance.name + " has not a valid value in file [ " + options.initialVariancesFile + " ].", 0);
        }
        results.variances.push_back(variance);
      }
      file.close();
    }
  }
  else if(options.computeREMLInSubsample == true)
  {
    misc.error("Error: An error was happened. Sorry for the inconvenience, subsampling function has to be updated.", 0);
    
    //If number of individuals is not enough, dissable subsampling.
    int avoidPrecomputing = 0;
    if(3*options.minimumSubsample > srcGRM->nIndividuals)
    {
      avoidPrecomputing = 1;
    }
    communicator->broadcast(&avoidPrecomputing);
    if(avoidPrecomputing == 1)
    {
      misc.message << "Random subsampling is disabled due the low number of individuals." << std::endl;
      this->subsampleREMLResults = results;
      return;
    }
    
    misc.message << "Estimating starting variance values from REML on random subsampling." << std::endl;
    misc.message.tab = "  ";
    misc.setGetElapsedTime("REMLSubsampling");
    
    //Start random subsampling
    std::vector<REMLResults> partialREMLResults;
    for(int i = 0; i < options.nSubSampleIterations; i++)
    {
      SingleREML singleREML(false);
      if(type == 0)
      {
        Kernel *grm = new Kernel(srcGRM);
        grm->randomSubSample(options.initialSubsampleFraction, options.minimumSubsample);
        singleREML.computeRaw(grm);
      }
      
      if(singleREML.reml != NULL)
      {
        if(singleREML.reml->success == true) //If REML finished successful, store results.
        {
          REMLResults temp;
          temp.variances = singleREML.reml->V->variances;
          temp.logLikelihood = singleREML.reml->logLikelihood;
          partialREMLResults.push_back(temp);
        }
      }
    }
    
    if(partialREMLResults.size() > 0) //Average results from all REML analysis
    {
      results.variances = partialREMLResults[0].variances;
      results.logLikelihood = partialREMLResults[0].logLikelihood;
      for(int i = 1; i < partialREMLResults.size(); i++) //Average variances
      {
        std::vector<Variance> variancesToAdd = partialREMLResults[i].variances;
        if(variancesToAdd.size() != results.variances.size())
        {
          misc.error("Error: An internal error was happened when estimating initial variances using random subsampling. Different number of variances.", 0);
        }
        for(int j = 0; j < results.variances.size(); j++)
        {
          if(results.variances[j].name != variancesToAdd[j].name)
          {
            misc.error("Error: An internal error was happened when estimating initial variances using random subsampling. Different variance names.", 0);
          }
          results.variances[j].variance += variancesToAdd[j].variance;
        }
      }
      
      for(int j = 0; j < results.variances.size(); j++)
      {
        results.variances[j].variance /= double(partialREMLResults.size());
      }
    } //End of results averaging
    
    misc.message.tab = "";
    misc.message << "Initial variance values have been estimated after " << misc.setGetElapsedTime("REMLSubsampling", true) << std::endl;
  }
  
  this->subsampleREMLResults = results;
}

