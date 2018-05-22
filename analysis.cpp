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

#include "analysis.h"
#include "singlereml.h"
#include "multireml.h"
#include "accuracybysnp.h"
#include "groupeffects.h"
#include "labeledmatrix.h"
#include "pcagentemp.h"
#include "mpresiduals.h"

#include <cstdlib>

Analysis::Analysis()
{
  unsigned int seed = options.randomSeed;
  std::srand ( seed );
}

Analysis::~Analysis()
{
}

void Analysis::makeGRM()
{
  if(options.regionalAnalysis == false)
  {
    Kernel *grm = loadGRMUsingOptions();
    bool sanitized = grm->sanitizeKernel();
    if( sanitized == true )
    {
      if( options.diagonalizeGRM == true )
      {
        if( options.writeAlsoNormalGRM == true )
        {
          grm->writeKernel(options.outFile + ".nondiagonal");
        }
        grm->diagonalizeKernel();
      }
      grm->writeKernel(options.outFile);
    }
    else
    {
      misc.message << "Sorry, a problem was happened while computing the GRM. Please, check the logs." << std::endl;
    }
    delete grm;
  }
  else
  {
    Genotype * genotype = loadGenotypeUsingOptions();
    genotype->groupSNPs(options.regionBy);
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

      options.outFile = baseOutFile + "." + group;
      Genotype *regionalGenotype = new Genotype();
      genotype->genotypeOfSNPsGroup(group, regionalGenotype);
      Kernel *grm = new Kernel(regionalGenotype);
      bool sanitized = grm->sanitizeKernel();
      if( sanitized == true )
      {
        if( options.diagonalizeGRM == true )
        {
          if( options.writeAlsoNormalGRM == true )
          {
            grm->writeKernel(options.outFile + ".nondiagonal");
          }
          grm->diagonalizeKernel();
        }
        grm->writeKernel(options.outFile);
      }
      else
      {
        misc.message << "WARNING: Sorry, a problem was happened while computing the GRM on region " << group << ". Please, check the logs." << std::endl;
      }
      delete grm;
      delete regionalGenotype;
    }
    options.outFile = baseOutFile;
  }
}

void Analysis::makeGRMAndStoreMostRelated()
{
  Kernel *grm = loadGRMUsingOptions();
  grm->writeKernel(options.outFile);
  
  Kernel *grmMostRelated = new Kernel(grm);
  grmMostRelated->keepWithRelatednessOutside(options.mostRelatedLowerThreshold, options.mostRelatedUpperThreshold);
  grmMostRelated->writeKernel(options.outFile + ".mostRelated");
  delete grmMostRelated;
  
  for(int i = 0; i<options.pruneThresholdsCheck.size(); i++)
  {
    Kernel *grmBase = new Kernel(grm);
    options.pruneGRM = true;
    options.grmCutoff = options.pruneThresholdsCheck[i];
    grmBase->sanitizeKernel();
    misc.message << int(grm->nIndividuals) - int(grmBase->nIndividuals) << " individuals have been filtered from " << grm->nIndividuals << " when cutoff is " << options.grmCutoff << ". (" << double(int(grm->nIndividuals) - int(grmBase->nIndividuals))/double(grm->nIndividuals) << ")" << std::endl;
    delete grmBase;
  }
  
  delete grm;
  
}

void Analysis::makeREML()
{
  if(options.regionalAnalysis == false)
  {
    SingleREML singleREML;
    singleREML.compute();
  }
  else
  {
    if( options.allRegionsTogether == false )
    {
      SingleREML singleREML;
      singleREML.computeRegional();
    }
    else
    {
      SingleREML singleREML;
      singleREML.computeMultipleGroups();
    }
  }
}

void Analysis::makeMultivarREML()
{
  if(options.regionalAnalysis == false) //Regional analysis?
  {
    MultiREML multiREML;
    multiREML.compute();
  }
  else
  {
    if( options.allRegionsTogether == false )
    {
      MultiREML multiREML;
      multiREML.computeRegional();
    }
    else
    {
      MultiREML multiREML;
      multiREML.computeMultipleGroups();
    }
  }
}

void Analysis::makeSimulatePhenotype()
{
  Genotype * genotype = loadGenotypeUsingOptions();
  Genotype * adjustEffectsGenotype = NULL;
  if( options.adjustEffectsGenotypeListFile != "" )
  {
    adjustEffectsGenotype = new Genotype();
    adjustEffectsGenotype->loadList(options.adjustEffectsGenotypeListFile);
  }
  SimulatePhenotype simulatePhenotype(genotype, options.effectsSizeFile, adjustEffectsGenotype);
  simulatePhenotype.simulatePhenotypes();
}

void Analysis::makePredictPhenotype()
{
  Genotype *genotype;
  if(options.genotypeFile != "") //Read single genotype file
  {
    genotype = new Genotype(options.genotypeFile);
    PredictPhenotype predictPhenotype(genotype, options.snpEffectsFile);
    predictPhenotype.predictPhenotypes();
    predictPhenotype.storePredictions();
  }
  else if (options.genotypeListFile != "") //Read multiple genotype files
  {
    misc.message << "Predicting phenotypes from genotypes specified in [ " << options.genotypeListFile << " ]..." << std::endl;
    misc.message.tab = "  ";
    misc.setGetElapsedTime("AddingEffects");
    std::vector<std::string> fileList;
    getListFromFile(options.genotypeListFile, fileList);

    genotype = new Genotype(fileList[0]);
    PredictPhenotype predictPhenotype(genotype, options.snpEffectsFile, false);
    predictPhenotype.predictPhenotypes();

    for(int i = 1; i < fileList.size(); i++)
    {
      genotype = new Genotype(fileList[i]);
      PredictPhenotype phenotypeToAdd(genotype, options.snpEffectsFile, false);
      phenotypeToAdd.predictPhenotypes();
      predictPhenotype.addMoreEffects(&phenotypeToAdd);
    }
    misc.message.tab = "";
    misc.message << "Genotype effects added successfully after " << misc.setGetElapsedTime("AddingEffects", true) << ". Storing results..." << std::endl;
    predictPhenotype.storePredictions();
  }
  else
  {
    misc.error("Error: No genotype file(s) specified.", 0);
  }
}

void Analysis::makePCA()
{
  Kernel *grm = loadGRMUsingOptions();
  bool sanitized = grm->sanitizeKernel();
  if( sanitized != true )
  {
    misc.message << "WARNING: Sorry, a problem was happened while computing GRM for performing a PCA. Please, check the logs." << std::endl;
  }
  PCA pca(grm);
  delete grm;
}

void Analysis::makeGWAS()
{
  GWAS gwas;
  gwas.computeGWAS();
}

void Analysis::makeRecursiveGWAS()
{
  GWAS gwas;
  gwas.computeGWAS();
}

void Analysis::makeComputeAccuracyBySNP()
{
  AccuracyBySNP accuracyBySNP;
}

void Analysis::makeEffectsAnalysis()
{
  LabeledMatrix * correlations;
  PCAGenTemp * pca;
  
  std::string baseOutFile = options.outFile;
  
  if(options.crossedEffectCorrelations == false)
  {
    //currently testing function
    //GroupEffects groupEffects("/home/ocanela/big/ukbiobank/gwas.byGenes.corrected2/results.byGene.withvariance/gwas.byGenes.corrected.50-0.0.Standing_height.outThrNone.rep.chr8.effects");
    //GroupEffects groupEffects("/home/ocanela/ris-lx01/ownsoftware/dissect/src/tmp/effects");
    //correlations->save("/home/ocanela/big/tmp/test.correlations");  
    
    GroupEffects groupEffects(options.groupEffectsFileList, row);
    
    std::vector<std::string> rowLabelsToKeep = groupEffects.effects->getRowLabels();
    std::vector<std::string> colLabelsToKeep = groupEffects.effects->getColLabels();
    if(options.fileGroupsToKeep != "")
    {
      misc.message << "Removing groups not specified in file [ " << options.fileGroupsToKeep << " ]..." << std::endl;
      if(communicator->mpiRoot == true)
      {
        std::vector<std::string> temp;
        getListFromFile(options.fileGroupsToKeep, temp);
        std::set<std::string> test(temp.begin(), temp.end());
        if( test.size() != temp.size() )
        {
          misc.error("Error: There are repeated groups in file [ " + options.fileGroupsToKeep + " ].", 0);
        }
        
        std::vector<std::string> colLabelsTemp = groupEffects.effects->getColLabels();
        colLabelsToKeep = intersectionStringVectors(2, &temp, &(colLabelsTemp));
        if( colLabelsToKeep.size() != temp.size() )
        {
          misc.message << "Not all groups specified in [ " << options.fileGroupsToKeep << " ] are in the current effects file." << std::endl;
        }
      }
    }
    if(options.fileIndividualsToKeep != "")
    {
      misc.message << "Removing individuals not specified in file [ " << options.fileIndividualsToKeep << " ]..." << std::endl;
      if(communicator->mpiRoot == true)
      {
        std::vector< std::vector<std::string> > temp0;
        getTableFromFile(options.fileIndividualsToKeep, temp0, 2);
        std::vector<std::string> temp;
        for(int i = 0; i < temp0.size(); i++)
        {
          temp.push_back(temp0[i][0] + "@" + temp0[i][1]);
        }
        std::set<std::string> test(temp.begin(), temp.end());
        if( test.size() != temp.size() )
        {
          misc.error("Error: There are repeated individuals in file [ " + options.fileIndividualsToKeep + " ].", 0);
        }
        
        std::vector<std::string> rowLabelsTemp = groupEffects.effects->getRowLabels();
        rowLabelsToKeep = intersectionStringVectors(2, &temp, &(rowLabelsTemp));
        if( rowLabelsToKeep.size() != temp.size() )
        {
          misc.message << "Not all individuals specified in [ " << options.fileIndividualsToKeep << " ] are in the current effects file." << std::endl;
        }
      }
    }
    if(options.fileGroupsToKeep != "" || options.fileIndividualsToKeep != "")
    {
      groupEffects.effects->filterRowsAndCols(rowLabelsToKeep, colLabelsToKeep);
      misc.message << groupEffects.effects->getMatrix()->nGlobCols << " groups after filtering." << std::endl;
      misc.message << groupEffects.effects->getMatrix()->nGlobRows << " individuals after filtering." << std::endl;
    }
    
    correlations = groupEffects.computeCorrelations(column);
    //correlations = groupEffects.computeCovariances(column);
    correlations->save(options.outFile + ".gene.correlations.unfiltered");
    delete correlations;
    
    //correlations = groupEffects.computeCorrelations(row);
    correlations = groupEffects.computeCovariances(row);
    //correlations = new LabeledMatrix();
    //correlations->getMatrix()->multiply(groupEffects.effects->getMatrix(), 'N', groupEffects.effects->getMatrix(), 'T');
    //correlations->save(options.outFile + ".indiv.covariances.unfiltered");
    options.outFile = baseOutFile + ".indiv.covariances.unfiltered";
    pca = new PCAGenTemp(correlations);
    delete correlations;
    options.outFile = baseOutFile;
    delete pca;
    
    
    //options.groupDistanceForDiscarding = 500000;
    groupEffects.filterCorrelatedGroups(column, 0.1, options.groupsPositions);
    correlations = groupEffects.computeCorrelations(column);
    correlations->save(options.outFile + ".gene.correlations." + i2s(options.groupDistanceForDiscarding));
    delete correlations;
    correlations = groupEffects.computeCovariances(column);
    options.outFile = baseOutFile + ".gene.covariances." + i2s(options.groupDistanceForDiscarding);
    pca = new PCAGenTemp(correlations);
    options.outFile = baseOutFile;
    delete pca;
    delete correlations;
    
    //correlations = groupEffects.computeCorrelations(row);
    correlations = groupEffects.computeCovariances(row);
    //correlations = new LabeledMatrix();
    //correlations->getMatrix()->multiply(groupEffects.effects->getMatrix(), 'N', groupEffects.effects->getMatrix(), 'T');
    //correlations->save(options.outFile + ".indiv.covariances." + i2s(options.groupDistanceForDiscarding));
    options.outFile = baseOutFile + ".indiv.covariances." + i2s(options.groupDistanceForDiscarding);
    pca = new PCAGenTemp(correlations);
    delete correlations;
    delete pca;
    options.outFile = baseOutFile;
    /*
    options.groupDistanceForDiscarding = 1000000;
    groupEffects.filterCorrelatedGroups(column, 0.1, "/home/ocanela/big/ukbiobank/gwas.byGenes.corrected2/notebooks/genes.info.quoted");
    correlations = groupEffects.computeCorrelations(column);
    correlations->save(options.outFile + ".correlations.1000000");
    delete correlations;
    
    options.groupDistanceForDiscarding = 2000000;
    groupEffects.filterCorrelatedGroups(column, 0.1, "/home/ocanela/big/ukbiobank/gwas.byGenes.corrected2/notebooks/genes.info.quoted");
    correlations = groupEffects.computeCorrelations(column);
    correlations->save(options.outFile + ".correlations.2000000");
    delete correlations;  
    */
  }
  else
  {
    std::vector<std::string> list1;
    std::vector<std::string> list2;
    for(int i = 0; i < options.groupEffectsPairFileList.size(); i+=2)
    {
      list1.push_back(options.groupEffectsPairFileList[ i ]);
      list2.push_back(options.groupEffectsPairFileList[ i + 1 ]);
    }
    GroupEffects groupEffects1(list1, row);
    GroupEffects groupEffects2(list2, row);
    
    Matrix * standardizedEffects1 = new Matrix( groupEffects1.effects->getMatrix() );
    standardizedEffects1->standardizeMatrix(column);
    Matrix * standardizedEffects2 = new Matrix( groupEffects2.effects->getMatrix() );
    standardizedEffects2->standardizeMatrix(column);
    
    LabeledMatrix * correlations = new LabeledMatrix();
    correlations->getMatrix()->multiply(standardizedEffects1, 'T', standardizedEffects2, 'N', 1./double(standardizedEffects1->nGlobRows));
    
    correlations->setRowLabels( groupEffects1.effects->getColLabels() );
    correlations->setColLabels( groupEffects2.effects->getColLabels() );
    correlations->save(options.outFile + ".gene.crossed.correlations");
    
    delete standardizedEffects1;
    delete standardizedEffects2;
    delete correlations;
  }
}

void Analysis::makeGLMMAnalysis()
{
  if(options.regionalAnalysis == false)
  {
    SingleREML singleREML;
    singleREML.compute();
  }
  else
  {
    misc.error("Error: An internal error has happened. Multi GLMM analysis not implemented, yet.", 0);
  }
}

void Analysis::makeIGWASAnalysis()
{
  IGWAS igwas;
  igwas.computeGWAS();
}

void Analysis::makePredictCovarPhenotype()
{
  Covariate covariate(options.covarsFile, options.qCovarsFile, std::vector<std::string>(), false);
  std::map<std::string, double> covariateEffects = covariate.loadEffectPrediction(options.fileCovarEffects, options.fileQCovarEffects);
  
  if(communicator->mpiRoot)
  {
    Message message(options.outFile + ".covars.predicted.phenos");
    message << "FID IID CPHENO" << std::endl;
    for(std::map<std::string, double>::iterator it = covariateEffects.begin(); it != covariateEffects.end(); ++it)
    {
      std::vector<std::string> individualIdPairs = splitString(it->first, "@");
      if(individualIdPairs.size() != 2)
      {
        misc.error("Error: An internal error has happened. Unexpected id: " + it->first + ".", 0);
      }
      message << individualIdPairs[0] << " " << individualIdPairs[1] << " " << it->second << std::endl;
    }
  }

}

void Analysis::makeMultiplePhenotypeGWAS()
{
  GWAS gwas;
  if(options.originalRedistributionMethod == true)
  {
    gwas.computeMultiplePhenotypeGWAS();
  }
  else
  {
    gwas.computeMultiplePhenotypeGWASGroupedCommunicator();
  }
}

void Analysis::makeMultiplePhenotypeResiduals()
{
  MPResiduals mpresiduals;
  LabeledMatrix * residuals = mpresiduals.computeResiduals(options.grmFile);
  residuals->save(options.outFile + ".residuals");
  delete residuals;
}

void Analysis::makeFilterLabeledMatrix()
{
  if(options.inputLabeledMatrix == options.outFile)
  {
    misc.error("Error: The inpute and output files are the same.", 0);
  }
  
  std::vector<std::string> keepRows;
  std::vector<std::string> keepColumns;
  
  getListFromFile(options.rowLabelsFile, keepRows);
  getListFromFile(options.columnLabelsFile, keepColumns);
  
  
  LabeledMatrix * labeledMatrix = NULL;
  labeledMatrix = new LabeledMatrix(options.inputLabeledMatrix);
  labeledMatrix->filterRowsAndCols(keepRows, keepColumns);
  labeledMatrix->save(options.outFile);
  
  delete labeledMatrix;
}

void Analysis::makeSNPStats()
{
  std::set<std::string> keepSNPs;
  if( options.fileSNPsToKeep != "" )
  {
    std::vector<std::string> keepSNPsVector;
    getListFromFile(options.fileSNPsToKeep, keepSNPsVector);
    keepSNPs = std::set<std::string>(keepSNPsVector.begin(), keepSNPsVector.end());
    misc.message << "Only SNPs specified in file [ " << options.fileSNPsToKeep << " ] will be kept." << std::endl;
  }
  
  std::set<std::string> keepIndividualIds = getIndividualIdsSpecifiedByOptionKeep();
  
  Genotype * genotype = NULL;
  if(options.genotypeFile != "")
  {
    genotype = new Genotype(options.genotypeFile, keepSNPs, keepIndividualIds);
  }
  else if(options.bgenGenotypeFile != "")
  {
    genotype = new Genotype(options.bgenGenotypeFile, keepSNPs, keepIndividualIds, GenotypeAttributes::probabilities);
  }
  else
  {
    misc.error("Error: A genotype file have been not specified.", 0);
  }
  
  if( communicator->mpiRoot == true )
  {
    Message message(options.outFile + ".freq");
    message << "CHROM SNP A1 A2 FREQ1 FREQ2 MAF" << std::endl;
    for(int i = 0; i < genotype->SNPs.size(); i++)
    {
      double maf = (genotype->SNPs[i].p1<genotype->SNPs[i].p2)?genotype->SNPs[i].p1:genotype->SNPs[i].p2;
      
      message << genotype->SNPs[i].chr;
      message << " " << genotype->SNPs[i].name;
      message << " " << genotype->SNPs[i].allele1;
      message << " " << genotype->SNPs[i].allele2;
      message << " " << genotype->SNPs[i].p1;
      message << " " << genotype->SNPs[i].p2;
      message << " " << maf;
      message << std::endl;
    }
  }
  
  delete genotype;
}
