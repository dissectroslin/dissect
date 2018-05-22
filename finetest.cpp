#include "finetest.h"
#include "genotype.h"
#include "phenotype.h"

#include <vector>
#include <string>
#include <set>


FineTest::FineTest(std::string genotypesFile, std::string phenotypeFile, std::vector< std::vector<std::string> > SNPPairs)
{
  Genotype genotypes(genotypesFile);
  Phenotype phenotypes(phenotypeFile);
  
  std::set<std::string> setSNPsList;
  for(int i = 0; i<SNPPairs.size(); i++)
  {
    setSNPsList.insert(SNPPairs[i][0]);
    setSNPsList.insert(SNPPairs[i][1]);
  }
  std::vector<std::string> SNPsList(setSNPsList.begin(), setSNPsList.end());
  
  std::vector<std::string> tempOrderedIndividuals = orderVectorAsTemplate(genotypes.individualIds, phenotypes.individualIds);
  genotypes.filterSNPsAndIndividuals(SNPsList, tempOrderedIndividuals);
  phenotypes.filterIndividuals(tempOrderedIndividuals);
  
  
  std::vector<std::vector<double> > globalGenotypes;
  genotypes.genotypes->matrixToStandardVector(globalGenotypes);
  
  std::vector<double> globalPhenotypes;
  phenotypes.phenotypes->matrixToStandardVector(globalPhenotypes);
  
  if(communicator->mpiRoot == true)
  {
    for(int i = 0; i<SNPPairs.size(); i++)
    {
      snp1 = SNPPairs[i][0];
      snp2 = SNPPairs[i][1];
      
      idx1 = genotypes.SNPIdsIdx[snp1];
      idx2 = genotypes.SNPIdsIdx[snp2];
      
      //Order snp1, snp2 as a funtion of the one that has the larger pvalue.
    }
  }
}

FineTest::~FineTest()
{
}

std::vector<std::vector<double>> FineTest::createContingencyTables(std::vector<std::vector<double> > globalGenotypes, std::vector<double> globalPhenotypes, Genotype & genotypes, Phenotype & phenotypes, SNP1idx, SNP2idx)
{
  if(communicator->mpiRoot == true)
  {
    std::vector< std::vector<double> > s1s2Full(3, std::vector<double>(3, 0.)); //The full contingency tables including the heterozygot and 2 homozygots
    std::vector< std::vector<double> > s1affFull(3, std::vector<double>(2, 0.));
    std::vector< std::vector<double> > s2affFull(3, std::vector<double>(2, 0.));
    
    for(int indGenoIdx = 0; indGenoIdx < genotypes.nIndividuals; indGenoIdx++)
    {
      indPhenoIdx = phenotypes.individualIdsIdx[ genotypes.individualIds[indGenoIdx] ];
      
      double temp = globalPhenotypes[ indPhenoIdx ];
      if(affectionStatus != 1 or affectionStatus != 2)
      {
        temp = 0;
      }
      int affectionStatus = int(temp) - 1;
      int s1a = int(globalGenotypes[SNP1Idx][indPhenoIdx])-1;
      int s2a = int(globalGenotypes[SNP2Idx][indPhenoIdx])-1;
      
      //s1s2 table
      if(s1a != -1 && s2a != -1)
      {
        s1s2Full[s1a][s2b]++;
      }
      if(s1a != -1 && affectionStatus != -1)
      {
        s1affFull[s1a][affectionStatus]++;
      }
      if(s2a != -1 && affectionStatus != -1)
      {
        s2affFull[s2a][affectionStatus]++;
      }
    }
    
    std::vector< std::vector<double> > s1s2(2, std::vector<double>(2, 0.));
    s1s2[0][0] = 4*s1s2Full[0][0] + 2*s1s2Full[1][0] + 2*s1s2Full[0][1] +  s1s2Full[1][1];
    s1s2[0][1] = 4*s1s2Full[0][2] + 2*s1s2Full[1][2] + 2*s1s2Full[0][1] +  s1s2Full[1][1];
    s1s2[1][0] = 4*s1s2Full[2][0] + 2*s1s2Full[1][0] + 2*s1s2Full[2][1] +  s1s2Full[1][1];
    s1s2[1][1] = 4*s1s2Full[2][2] + 2*s1s2Full[1][2] + 2*s1s2Full[2][1] +  s1s2Full[1][1];
    
    std::vector< std::vector<double> > s1aff(2, std::vector<double>(2, 0.));
    s1aff[0][0] = 2*s1affFull[0][0] + s1affFull[1][0];
    s1aff[0][1] = 2*s1affFull[0][1] + s1affFull[1][1];
    s1aff[1][0] = 2*s1affFull[2][0] + s1affFull[1][0];
    s1aff[1][1] = 2*s1affFull[2][1] + s1affFull[1][1];
    
    std::vector< std::vector<double> > s2aff(2, std::vector<double>(2, 0.));
    s2aff[0][0] = 2*s2affFull[0][0] + s2affFull[1][0];
    s2aff[0][1] = 2*s2affFull[0][1] + s2affFull[1][1];
    s2aff[1][0] = 2*s2affFull[2][0] + s2affFull[1][0];
    s2aff[1][1] = 2*s2affFull[2][1] + s2affFull[1][1];
    
  }
}