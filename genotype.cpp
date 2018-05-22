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

#include "genotype.h"
#include "matrix.h"
#include "options.h"
#include "misc.h"
#include "global.h"
#include "communicator.h"
#include "auxiliar.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <cmath>

SNP::SNP()
{
  this->chr = "NA";
  this->name = "NA";
  this->position = -1;
  
  this->allele1 = "X";
  this->allele2 = "X";
  
  this->p1 = 0.;
  this->p2 = 0.;
  this->standardDev = 0.;
  this->nNonMissing = 0;
  for(int i=0; i<4; i++)
  {
    this->frequencies[i] = 0;
  }
}


SNP::SNP(std::string chr, std::string name, std::string position, std::string allele1, std::string allele2)
{
  this->chr = chr;
  this->name = name;
  
  std::stringstream ss;
  ss << position;
  if( (ss >> this->position).fail() )
  {
    misc.error("The SNP " + name + " has an invalid coordinate: " + position, 0);
  }
  
  this->allele1 = allele1;
  this->allele2 = allele2;
  
  this->p1 = 0.;
  this->p2 = 0.;
  this->standardDev = 0.;
  this->nNonMissing = 0;
  for(int i=0; i<4; i++)
  {
    this->frequencies[i] = 0;
  }
}

Genotype::Genotype(DistributionType dist, int ngr, int ngc, int nbr, int nbc)
{
  this->type = GenotypeAttributes::empty;
  
  this->genotypes = new Matrix(dist, ngr, ngc, nbr, nbc);
  this->missings = new Matrix(dist, ngr, ngc, nbr, nbc);
  
  this->normalized = false;
  
  this->groupedBy = ungrouped;
  this->groupedSNPs.clear();
}

Genotype::Genotype()
{
  this->type = GenotypeAttributes::empty;
  
  this->genotypes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->missings = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  
  this->normalized = false;
  
  this->groupedBy = ungrouped;
  this->groupedSNPs.clear();
}

Genotype::Genotype(std::string fname, bool loadGenotypes)
{
  this->type = GenotypeAttributes::empty;
  
  this->genotypes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->missings = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  
  this->normalized = false;
  
  this->groupedBy = ungrouped;
  this->groupedSNPs.clear();
  
  load(fname, loadGenotypes);
}

Genotype::Genotype( std::string fbase, std::set<std::string> & keepSNPs, std::set<std::string> & keepIndividualIds )
{
  this->type = GenotypeAttributes::empty;
  
  Genotype * baseGenotype = new Genotype(fbase, false);
  load(fbase, baseGenotype, keepSNPs, keepIndividualIds);
  delete baseGenotype;
}

Genotype::Genotype( std::string fbase, Genotype * baseGenotype, std::set<std::string> & keepSNPs, std::set<std::string> & keepIndividualIds )
{
  this->type = GenotypeAttributes::empty;
  
  load(fbase, baseGenotype, keepSNPs, keepIndividualIds);
}

Genotype::Genotype( std::string fbase, std::set<std::string> & keepSNPs, std::set<std::string> & keepIndividualIds,  GenotypeAttributes::Type formatType)
{
#ifdef BGEN
  this->originalFileBEDSNPIdxs.clear();
  this->originalFileBedIndividualsIdxs.clear();
  
  this->genotypes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->missings = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  
  this->normalized = false;
  
  this->groupedBy = ungrouped;
  this->groupedSNPs.clear();
  
  if( formatType == GenotypeAttributes::probabilities )
  {
    this->type = GenotypeAttributes::probabilities;
    readBGENFile(fbase, keepSNPs, keepIndividualIds);
  }
  else
  {
    misc.error("Error: An internal error has happened. Invalid format type when loading genotypes.", 0);
  }
#endif
}

void Genotype::load( std::string fbase, Genotype * baseGenotype, std::set<std::string> & keepSNPs, std::set<std::string> & keepIndividualIds )
{
  if(baseGenotype->originalFileBEDSNPIdxs.find(fbase) == originalFileBEDSNPIdxs.end() || baseGenotype->originalFileBedIndividualsIdxs.find(fbase) == originalFileBedIndividualsIdxs.end())
  {
    misc.error("Error: An internal error has happened. Unexpected file key when loading filtered genotypes.", 0);
  }
  
  this->type = GenotypeAttributes::calls;
  
  this->originalFileBEDSNPIdxs.clear();
  this->originalFileBedIndividualsIdxs.clear();
  
  this->genotypes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->missings = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  
  this->normalized = false;
  
  this->groupedBy = ungrouped;
  this->groupedSNPs.clear();
  
  //Get SNPs
  std::vector<int> SNPsBEDIdxsToKeep;
  if(communicator->mpiRoot)
  {
    this->SNPs.clear();
    this->SNPIds.clear();
    this->SNPIdsIdx.clear();
    int idx = 0;
    for(int i = 0; i < baseGenotype->SNPs.size(); i++)
    {
      if( keepSNPs.size() != 0 && keepSNPs.find(baseGenotype->SNPs[i].name) == keepSNPs.end() )
      {
        continue;
      }
      
      SNP snp = SNP(baseGenotype->SNPs[i].chr, baseGenotype->SNPs[i].name, i2s(baseGenotype->SNPs[i].position), baseGenotype->SNPs[i].allele1, baseGenotype->SNPs[i].allele2);
      
      SNPsBEDIdxsToKeep.push_back( getMapValue(baseGenotype->originalFileBEDSNPIdxs[fbase], baseGenotype->SNPs[i].name) );
      
      this->SNPs.push_back(snp);
      this->SNPIds.push_back(snp.name);
      if(this->SNPIdsIdx.count(snp.name) != 0)
      {
        misc.error("Error: An internal error has happened. The SNP with name: " + snp.name + " appears more than once in the source genotype.", 0);
      }
      this->SNPIdsIdx[snp.name] = idx;
      idx++;
    }
    this->nSNPs = SNPs.size();
  }
  communicator->broadcast(&this->nSNPs);
  communicator->broadcast(SNPsBEDIdxsToKeep);
  
  
  //Get Individuals
  std::vector<int> individualBEDIdxsMaskToKeep(baseGenotype->originalFileBedIndividualsIdxs[fbase].size(), 0);
  if(communicator->mpiRoot)
  {
    this->individuals.clear();
    this->individualIds.clear();
    this->individualIdsIdx.clear();
    int idx = 0;
    for(int i = 0; i < baseGenotype->individuals.size(); i++)
    {
      Individual individual = baseGenotype->individuals[i];
      std::string key = individual.familyID + "@" + individual.individualID;
      
      if( keepIndividualIds.size() != 0 && keepIndividualIds.find(key) == keepIndividualIds.end() )
      {
        continue;
      }
      individualBEDIdxsMaskToKeep[ getMapValue(baseGenotype->originalFileBedIndividualsIdxs[fbase], key) ] = 1;
      
      this->individuals.push_back(individual);
      this->individualIds.push_back(key);
      if(this->individualIdsIdx.count(key) != 0)
      {
        misc.error("Error: The individual with family Id: " + individual.familyID + " and individual Id: " + individual.individualID + " appears more than one time in the genotypes FAM file.", 0);
      }
      this->individualIdsIdx[key] = idx;
      idx++;
    }
    this->nIndividuals = this->individuals.size();
  }
  communicator->broadcast(&this->nIndividuals);
  communicator->broadcast(individualBEDIdxsMaskToKeep);
  
  if(SNPsBEDIdxsToKeep.size() == 0 || this->nIndividuals == 0) // individualBEDIdxsMaskToKeep.size() == this->nIndividuals || 
  {
    misc.error("Error: An internal error has happened. The number of SNPs and/or individuals to load is empty.", 0);
  }
  
  readBEDFile(fbase, SNPsBEDIdxsToKeep, individualBEDIdxsMaskToKeep);
}

Genotype::~Genotype()
{
  if(this->genotypes != NULL)
  {
    delete this->genotypes;
    this->genotypes = NULL;
  }
  if(this->missings != NULL)
  {
    delete this->missings;
    this->missings = NULL;
  }
}

void Genotype::load(std::string f, bool loadGenotypes)
{
  this->type = GenotypeAttributes::calls;
  
  std::set<std::string> keepSNPs;
  if(options.fileSNPsToKeep != "")
  {
    misc.message << "Only SNPs specified in file [ " << options.fileSNPsToKeep << " ] will be kept." << std::endl;
    std::vector<std::string> temp;
    getListFromFile(options.fileSNPsToKeep, temp);
    keepSNPs = std::set<std::string>(temp.begin(), temp.end());
  }
  std::set<std::string> keepIndividualIds;
  if(options.fileIndividualsToKeep != "")
  {
    misc.message << "Only individuals specified in file [ " << options.fileIndividualsToKeep << " ] will be kept." << std::endl;
    std::vector< std::vector<std::string> > temp;
    getTableFromFile(options.fileIndividualsToKeep, temp, 2);
    for(int ir = 0; ir < temp.size(); ir++)
    {
      keepIndividualIds.insert(temp[ir][0] + "@" + temp[ir][1]);
    }
  }
  
  std::vector<int> SNPsBEDIdxsToKeep;
  std::vector<int> individualBEDIdxsMaskToKeep;
  readBIMFile(f, keepSNPs, SNPsBEDIdxsToKeep);
  readFAMFile(f, keepIndividualIds, individualBEDIdxsMaskToKeep);
  
  if(this->nSNPs == 0 || this->nIndividuals == 0)
  {
    misc.error("Error: There are 0 SNPs and/or 0 individuals kept for the analysis. Aborting...", 0);
  }
  
  if( loadGenotypes == true )
  {
    readBEDFile(f, SNPsBEDIdxsToKeep, individualBEDIdxsMaskToKeep);
  }
  else
  {
    if(this->genotypes != NULL)
    {
      delete this->genotypes;
      this->genotypes = NULL;
    }
    if(this->missings != NULL)
    {
      delete this->missings;
      this->missings = NULL;
    }
  }
}

void Genotype::filterGenotypeUsingOptions()
{
  if( this->spaceChangedKernelName != "" )
  {
    misc.error("Error: An internal error has happened. Operation not permitted on a transformed genotype.", 0);
  }
  
  if(options.fileSNPsToKeep != "")
  {
    misc.message << "Removing SNPs not specified in file [ " << options.fileSNPsToKeep << " ]..." << std::endl;
    std::vector<std::string> SNPsToKeep;
    if(communicator->mpiRoot == true)
    {
      std::vector<std::string> temp;
      getListFromFile(options.fileSNPsToKeep, temp);
      std::set<std::string> test(temp.begin(), temp.end());
      if( test.size() != temp.size() )
      {
        misc.error("Error: There are repeated SNPs in file [ " + options.fileSNPsToKeep + " ].", 0);
      }
      
      SNPsToKeep = intersectionStringVectors(2, &temp, &(this->SNPIds));
      if( SNPsToKeep.size() != temp.size() )
      {
        misc.message << "Not all SNPs specified in [ " << options.fileSNPsToKeep << " ] are in the current genotype file." << std::endl;
      }
    }
    std::vector<std::string> currentIndividuals = this->individualIds;
    filterSNPsAndIndividuals(SNPsToKeep, currentIndividuals);
    misc.message << SNPIds.size() << " SNPs after filtering." << std::endl;
  }
}

void Genotype::loadList(std::string f)
{
  misc.message << "Reading genotypes from files specified in [ " << f << " ]..." << std::endl;
  misc.message.tab = "  ";

  std::vector<std::string> fileNames;
  getListFromFile(f, fileNames);
  
  for(int idx = 0; idx<fileNames.size(); idx++)
  {
    std::string partialFileName = fileNames[idx];
    
    if( idx == 0 )
    {
      this->load(partialFileName);
    }
    else
    {
      misc.setGetElapsedTime("AddGenotype");
      misc.message << "Adding genotype..." << std::endl;
      Genotype partialGenotype = Genotype(partialFileName);
      appendGenotype(&partialGenotype);
      misc.message << "Adding genotype needed " << misc.setGetElapsedTime("AddGenotype", true) << std::endl;
    }
  }

  misc.message.tab = "";
  misc.message << "Reading genotype files ended." << std::endl;
}

void Genotype::readBIMFile(std::string f, std::set<std::string> & keepSNPs, std::vector<int> & SNPsBEDIdxsToKeep)
{
  std::ifstream file;
  std::string line;
  
  if(this->originalFileBEDSNPIdxs.count(f) != 0)
  {
    misc.error("Error: The file [ " + f + ".bim/bed/fam ] appears defined twice.", 0);
  }
  this->originalFileBEDSNPIdxs[f].clear();
  
  SNPsBEDIdxsToKeep.clear();
  
  if(communicator->mpiRoot)
  {
    misc.message << "Reading SNPs data from file [ " << f << ".bim ] ..." << std::endl;
    
    misc.checkFileExists(f + ".bim");
    file.open((f + ".bim").c_str());

    this->SNPs.clear();
    this->SNPIds.clear();
    this->SNPIdsIdx.clear();
    int idx = 0;
    int ibed = 0;
    while(getline(file,line))
    {
      if(!file)
      {
	break;
      }
      std::istringstream sstemp(line); //->Comprova que te 6 elements.
      
      std::string temp[6];
      sstemp >> temp[0] >> temp[1] >> temp[2] >> temp[3] >> temp[4] >> temp[5];
      
      SNP snp = SNP(temp[0], temp[1], temp[3], temp[4], temp[5]);
      
      this->originalFileBEDSNPIdxs[f][snp.name] = ibed;
      if( keepSNPs.size() != 0 && keepSNPs.find(snp.name) == keepSNPs.end() )
      {
        ibed++;
        continue;
      }
      SNPsBEDIdxsToKeep.push_back(ibed);
      ibed++;
      
      this->SNPs.push_back(snp);
      this->SNPIds.push_back(snp.name);
      if(this->SNPIdsIdx.count(snp.name) != 0)
      {
        misc.error("Error: The SNP with name: " + snp.name + " appears more than one time in the genotypes BIM file.", 0);
      }
      this->SNPIdsIdx[snp.name] = idx;
      idx++;
    }
    file.close();

    this->nSNPs = SNPs.size();
  
    if(keepSNPs.size() == 0)
    {
      misc.message << this->nSNPs << " SNPs found." << std::endl;
      SNPsBEDIdxsToKeep.clear();
    }
    else
    {
      misc.message << ibed << " SNPs found. " << this->nSNPs << " kept for the analysis." << std::endl;
    }
  }

  communicator->broadcast(&this->nSNPs);
  communicator->broadcast(SNPsBEDIdxsToKeep);
}

void Genotype::readFAMFile(std::string f, std::set<std::string> & keepIndividualIds, std::vector<int> & individualBEDIdxsMaskToKeep)
{
  std::ifstream file;
  std::string line;
  
  if(this->originalFileBedIndividualsIdxs.count(f) != 0)
  {
    misc.error("Error: The file [ " + f + ".bim/bed/fam ] appears defined twice.", 0);
  }
  this->originalFileBedIndividualsIdxs[f].clear();
  
  individualBEDIdxsMaskToKeep.clear();
  
  if(communicator->mpiRoot)
  {
    misc.message << "Reading Individuals data from file [ " << f << ".fam ] ..." << std::endl;
    misc.checkFileExists(f + ".fam");
    file.open((f + ".fam").c_str());
    
    this->individuals.clear();
    this->individualIds.clear();
    this->individualIdsIdx.clear();
    int idx = 0;
    int ibed = 0;
    while(getline(file,line))
    {
      std::istringstream sstemp(line); //->Comprova que te 6 elements.
      
      Individual individual;
      
      sstemp >> individual.familyID;
      sstemp >> individual.individualID;
      sstemp >> individual.paternalID;
      sstemp >> individual.maternalID;
      sstemp >> individual.sex;
      
      std::string key = individual.familyID + "@" + individual.individualID;
      
      if( (sstemp >> individual.phenotype).fail() )
      {
	misc.error("Error. The phenotype of the individual " + individual.familyID + " " + individual.individualID + " is not an integer.", 0);
      }
      
      this->originalFileBedIndividualsIdxs[f][key] = ibed;
      if( keepIndividualIds.size() != 0 && keepIndividualIds.find(key) == keepIndividualIds.end() )
      {
        individualBEDIdxsMaskToKeep.push_back(0);
        ibed++;
        continue;
      }
      individualBEDIdxsMaskToKeep.push_back(1);
      ibed++;
      
      this->individuals.push_back(individual);
      this->individualIds.push_back(key);
      if(this->individualIdsIdx.count(key) != 0)
      {
        misc.error("Error: The individual with family Id: " + individual.familyID + " and individual Id: " + individual.individualID + " appears more than one time in the genotypes FAM file.", 0);
      }
      this->individualIdsIdx[key] = idx;
      idx++;
    }
    
    file.close();
    
    this->nIndividuals = this->individuals.size();
    
    if(keepIndividualIds.size() == 0)
    {
      misc.message << this->nIndividuals << " individuals found." << std::endl;
    }
    else
    {
      misc.message << individualBEDIdxsMaskToKeep.size() << " individuals found. " << this->nIndividuals << " kept for the analysis." << std::endl;
    }
  }
  
  communicator->broadcast(&this->nIndividuals);
  communicator->broadcast(individualBEDIdxsMaskToKeep);
}

void Genotype::readBEDFile(std::string f, std::vector<int> & SNPsBEDIdxsToKeep, std::vector<int> & individualBEDIdxsMaskToKeep)
{
  std::ifstream file;
  unsigned char * header; //Temporal array where header information will be stored.
  unsigned char * memBlock; //Temporal array where a file row will be stored
  double * genotypeMatrixBlock; //Temporal matrix where a block of genotypes will be stored
  double * missingMatrixBlock; //Temporal matrix where a block of missing information will be stored
  int nrf; //Number of bytes that store 4 genotypes
  int rest; //Number of genotypes in the last byte.
  int inc; //1 if rest!=0, 0 otherwise

  misc.setGetElapsedTime("GenotypeLoad");
  misc.message << "Reading genotype data from file [ " << f << ".bed ] ..." << std::endl;
  
  this->genotypes->initParameters(this->nSNPs, this->nIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockRows);
  this->missings->initParameters(this->nSNPs, this->nIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockRows);
  
  if(individualBEDIdxsMaskToKeep.size() == 0)
  {
    misc.error("Error: An internal error has happened. Empty individuals mask when loading gneotype data.", 0);
  }
  
  //Load the BED header. Check if the file is valid BED file and if it is in SNP-Major mode
  if(communicator->myCol == 0)
  {
    misc.checkFileExists(f + ".bed");
    file.open((f + ".bed").c_str(),std::ios::binary);
    
    //Checking the BED header.
    header=new unsigned char [3];
    file.read((char*)( &header[0]), 3);
    
    //Is a BED file?
    if(header[0] != 0x6C || header[1]!=0x1B)
    {
      misc.error("Error: The data file is not a valid binary BED > v1.0", 0);
    }
    
    //The file is in SNP-Major mode?
    if(header[2]==0x1)
    {
      int temp = individualBEDIdxsMaskToKeep.size();
      nrf = temp/4;
      rest = temp%4;
      if(((nrf*4)+rest)!=temp)
      {
	misc.error("Error: An internal error was happened.", 0);
      }
      misc.message << "Reading the binary file in SNP-major mode." << std::endl;
    }
    else if(header[2]==0x0)
    {
      misc.error("Error: The BED file is in Individual-major mode. This mode is not currently supported. Please, read PLINK manual.", 0);
    }
    else
    {
      misc.error("Error: The data file is not a valid binary BED.\nIt is not specified whether the file is in SNP-major o Individual-major mode.", 0);
    }
    delete [] header;
    
    inc=0;
    if(rest!=0)
    {
      inc = 1;
    }
    memBlock = new unsigned char [nrf + inc + 1];
    genotypeMatrixBlock = new double [this->nIndividuals*this->genotypes->nBlockRows];
    missingMatrixBlock = new double [this->nIndividuals*this->genotypes->nBlockRows];
    
    if(memBlock==NULL || genotypeMatrixBlock==NULL || missingMatrixBlock==NULL)
    {
      misc.error("Error: There are not enough memory for reading the binary file.", 0);
    }
  }
  
  //Resort the individual and SNP indexs to be loaded.
  if(SNPsBEDIdxsToKeep.size() != 0)
  {
    if(this->nSNPs != SNPsBEDIdxsToKeep.size())
    {
      misc.error("Error: An internal error has happened when loading a BED file. Unexpected number of BED SNP indexs to load.", 0);
    }
  }
  
  //Read genotype data
  std::vector<int> loadedSNPsByThisProcess;
  loadedSNPsByThisProcess.clear();
  if(communicator->mpiRoot == false)
  {
    this->SNPs = std::vector<SNP>(this->nSNPs, SNP());
  }
  for(int i=0; i<this->nSNPs; i += this->genotypes->nBlockRows)
  {
    //Which process will read and scatter the data?
    int sendProcessRow;
    int sendProcessCol;
    sendProcessRow = i/this->genotypes->nBlockRows;
    sendProcessRow = sendProcessRow % communicator->nProcRows;
    sendProcessCol = 0;
    
    //Read a block of genotypes of size this->genotypes->nBlockRows*this->nIndividuals
    if(communicator->myCol == 0)
    {
      int nFileRows = this->genotypes->nBlockRows;
      if((this->nSNPs-i) < nFileRows)
      {
        nFileRows = this->nSNPs-i;
      }
      
      bool seekFailed = false;
      if(communicator->myRow == sendProcessRow)
      {
        if(SNPsBEDIdxsToKeep.size() == 0) //Evaluate whether read sequentally all SNPs, or only those preselected indexs.
        {
          long offset = long(i)*(long(nrf) + long(inc)) + long(3);
          file.seekg(offset);
          seekFailed = seekFailed || file.fail();
          for(int j=0; j<nFileRows; j++)
          {
            file.read((char*)(&memBlock[0]), nrf + inc);
            parseSNP(nrf, rest, memBlock, &(genotypeMatrixBlock[j]), &(missingMatrixBlock[j]), &(SNPs[i+j]), individualBEDIdxsMaskToKeep); //Change from row major to col major.
            loadedSNPsByThisProcess.push_back(i+j);
          }
        }
        else
        {
          long offset = long(SNPsBEDIdxsToKeep[i])*(long(nrf) + long(inc)) + long(3);
          file.seekg(offset);
          seekFailed = seekFailed || file.fail();
          for(int j=0; j<nFileRows; j++)
          {
            if( j != 0 )
            {
              if( SNPsBEDIdxsToKeep[i + j] != SNPsBEDIdxsToKeep[i + j - 1] + 1 )
              {
                offset = long(SNPsBEDIdxsToKeep[i + j])*(long(nrf) + long(inc)) + long(3);
                file.seekg(offset);
                seekFailed = seekFailed || file.fail();
              }
            }
            file.read((char*)(&memBlock[0]), nrf + inc);
            parseSNP(nrf, rest, memBlock, &(genotypeMatrixBlock[j]), &(missingMatrixBlock[j]), &(SNPs[i+j]), individualBEDIdxsMaskToKeep); //Change from row major to col major.
            loadedSNPsByThisProcess.push_back(i+j);
          }
        }
        if(file.good() == false || file.fail() == true || seekFailed == true)
        {
          misc.error("Error: An internal error has happened when reading the genotypes file.", 0);
        }
      }
    }
    
    //Distribute the read block of genotypes between the genotypes matrix
    for(int k=0; k<this->nIndividuals; k+=this->genotypes->nBlockCols)
    {
      this->genotypes->scatterBlock((genotypeMatrixBlock+k*this->genotypes->nBlockRows), i, k, this->genotypes->nBlockRows, sendProcessRow, sendProcessCol);
      this->missings->scatterBlock((missingMatrixBlock+k*this->genotypes->nBlockRows), i, k, this->genotypes->nBlockRows, sendProcessRow, sendProcessCol);
    }
    
  }
  
  //Release temporal allocated arrays
  if(communicator->myCol == 0)
  {
    delete [] memBlock;
    delete [] genotypeMatrixBlock;
    delete [] missingMatrixBlock;
    
    file.close();
  }
  
  gatherSNPData(loadedSNPsByThisProcess);

  Cblacs_barrier(communicator->context, "All");
  
  misc.message << "Genotype loaded after " << misc.setGetElapsedTime("GenotypeLoad", true) << std::endl;
}

void Genotype::parseSNP(int nrf, int rest, unsigned char * memBlock, double *genotypes, double *missings, SNP* snp, const std::vector<int> & individualBEDIdxsMaskToKeep)
{
  unsigned int keepIdxs = 0;
  for(int j=0; j<nrf; j++)
  {
    parseSNPbyte(memBlock[j], 4, 0, genotypes, missings, snp, keepIdxs, individualBEDIdxsMaskToKeep, j*4);
  }
  parseSNPbyte(memBlock[nrf], rest, 0, genotypes, missings, snp, keepIdxs, individualBEDIdxsMaskToKeep, nrf*4);
  
  snp->nNonMissing = snp->frequencies[1] + snp->frequencies[2] + snp->frequencies[3];
  snp->p1 = ( 2.*double(snp->frequencies[1]) + double(snp->frequencies[2]) )/(2.*double(snp->nNonMissing));
  snp->p2 = ( 2.*double(snp->frequencies[3]) + double(snp->frequencies[2]) )/(2.*double(snp->nNonMissing));
  snp->standardDev = sqrt(2.*snp->p1*(1.-snp->p1));
}

void Genotype::parseSNPbyte(unsigned char g, int n, int d, double * genotypes, double * missings, SNP* snp, unsigned int & keepIdxs, const std::vector<int> & individualBEDIdxsMaskToKeep, int maskShift)
{
  unsigned char temp1, temp2;
  
  temp1 = g;
  temp1 = temp1 >> d;
  
  for(unsigned int ibed=0; ibed<n; ibed++)
  {
    if(individualBEDIdxsMaskToKeep[maskShift + ibed] != 0)
    {
      temp2 = temp1 & 0x3;
      if( temp2 == 0x0 )  //homozygot
      {
        genotypes[keepIdxs*this->genotypes->nBlockRows] = 1;
        missings[keepIdxs*this->genotypes->nBlockRows] = 1;
        snp->frequencies[ 1 ]++;
      }
      else if( temp2 == 0x3 ) //homozygot
      {
        genotypes[keepIdxs*this->genotypes->nBlockRows] = 3;
        missings[keepIdxs*this->genotypes->nBlockRows] = 1;
        snp->frequencies[ 3 ]++;
      }
      else if( temp2 == 0x2 ) //heterozygot
      {
        genotypes[keepIdxs*this->genotypes->nBlockRows] = 2;
        missings[keepIdxs*this->genotypes->nBlockRows] = 1;
        snp->frequencies[ 2 ]++;
      }
      else if( temp2 == 0x1 ) //missing
      {
        genotypes[keepIdxs*this->genotypes->nBlockRows] = 0;
        missings[keepIdxs*this->genotypes->nBlockRows] = 0;
        snp->frequencies[ 0 ]++;
      }
      else
      {
        misc.error("Error: An error happened when reading the binary file.", 0);
      }
      keepIdxs++;
    }

    temp1 = temp1 >> 2;
    
  }
}

void Genotype::gatherSNPData(std::vector<int> & loadedSNPsByThisProcess)
{
  std::vector<int> procNNonMissing;
  std::vector<int> procFrequencies0;
  std::vector<int> procFrequencies1;
  std::vector<int> procFrequencies2;
  std::vector<int> procFrequencies3;
  std::vector<double> procP1;
  std::vector<double> procP2;
  std::vector<double> procStandardDev;
  
  for(int i = 0; i < loadedSNPsByThisProcess.size(); i++)
  {
    int idx = loadedSNPsByThisProcess[i];
    procNNonMissing.push_back(this->SNPs[idx].nNonMissing);
    procFrequencies0.push_back(this->SNPs[idx].frequencies[0]);
    procFrequencies1.push_back(this->SNPs[idx].frequencies[1]);
    procFrequencies2.push_back(this->SNPs[idx].frequencies[2]);
    procFrequencies3.push_back(this->SNPs[idx].frequencies[3]);
    procP1.push_back(this->SNPs[idx].p1);
    procP2.push_back(this->SNPs[idx].p2);
    procStandardDev.push_back(this->SNPs[idx].standardDev);
  }
  
  int * globalLoadedSNPs;
  int * globalNNonMissing;
  int * globalFrequencies0;
  int * globalFrequencies1;
  int * globalFrequencies2;
  int * globalFrequencies3;
  double * globalP1;
  double * globalP2;
  double * globalStandardDev;
  
  int partialSize;
  int totalSize = 0;
  globalLoadedSNPs = communicator->asymmetricGather(&(loadedSNPsByThisProcess[0]), loadedSNPsByThisProcess.size(), &partialSize);
  totalSize += partialSize;
  
  globalNNonMissing = communicator->asymmetricGather(&(procNNonMissing[0]), procNNonMissing.size(), &partialSize);
  totalSize += partialSize;
  
  globalFrequencies0 = communicator->asymmetricGather(&(procFrequencies0[0]), procFrequencies0.size(), &partialSize);
  totalSize += partialSize;
  
  globalFrequencies1 = communicator->asymmetricGather(&(procFrequencies1[0]), procFrequencies1.size(), &partialSize);
  totalSize += partialSize;
  
  globalFrequencies2 = communicator->asymmetricGather(&(procFrequencies2[0]), procFrequencies2.size(), &partialSize);
  totalSize += partialSize;
  
  globalFrequencies3 = communicator->asymmetricGather(&(procFrequencies3[0]), procFrequencies3.size(), &partialSize);
  totalSize += partialSize;
  
  globalP1 = communicator->asymmetricGather(&(procP1[0]), procP1.size(), &partialSize);
  totalSize += partialSize;
  
  globalP2 = communicator->asymmetricGather(&(procP2[0]), procP2.size(), &partialSize);
  totalSize += partialSize;
  
  globalStandardDev = communicator->asymmetricGather(&(procStandardDev[0]), procStandardDev.size(), &partialSize);
  totalSize += partialSize;
  
  if(communicator->mpiRoot)
  {
    if(totalSize != 9*this->nSNPs)
    {
      misc.error("Error: An internal error was happened while gathering SNP information between threads.", 0);
    }
    
    for(int i = 0; i < this->nSNPs; i++)
    {
      int idx = globalLoadedSNPs[i];
      this->SNPs[idx].nNonMissing = globalNNonMissing[i];
      this->SNPs[idx].frequencies[0] = globalFrequencies0[i];
      this->SNPs[idx].frequencies[1] = globalFrequencies1[i];
      this->SNPs[idx].frequencies[2] = globalFrequencies2[i];
      this->SNPs[idx].frequencies[3] = globalFrequencies3[i];
      this->SNPs[idx].p1 = globalP1[i];
      this->SNPs[idx].p2 = globalP2[i];
      this->SNPs[idx].standardDev = globalStandardDev[i];
    }
    
    delete [] globalLoadedSNPs;
    delete [] globalNNonMissing;
    delete [] globalFrequencies0;
    delete [] globalFrequencies1;
    delete [] globalFrequencies2;
    delete [] globalFrequencies3;
    delete [] globalP1;
    delete [] globalP2;
    delete [] globalStandardDev;
  }
  else
  {
    this->SNPs.clear();
  }
}

void Genotype::normalizeGenotypes()
{
  double *vSnpsP, *vSnpsStdDev;
  
  if(this->genotypes == NULL)
  {
    misc.error("Error: An internal error has happened. Empty genotype data cannot be normalized.", 0);
  }
  
  if(this->normalized)
  {
    return;
  }
  this->normalized = true;
  
  if( this->spaceChangedKernelName != "" )
  {
    misc.error("Error: An internal error has happened. Operation not permitted on a transformed genotype.", 0);
  }
  
  if( this->type != GenotypeAttributes::calls )
  {
    misc.error("Error: An internal error has happened. Unexpected non-normalized imputed genotypes.", 0);
  }
  
  //Put loci standard deviation and mean in a vector
  if(communicator->mpiRoot)
  {
    std::vector<std::string> WrongSNPs;
    
    vSnpsP = new double [this->nSNPs];
    vSnpsStdDev = new double [this->nSNPs];
    for(int i=0; i<this->nSNPs; i++)
    {
      vSnpsP[i] = 2.*this->SNPs[i].p2;
      vSnpsStdDev[i] = this->SNPs[i].standardDev;
      if(vSnpsStdDev[i] == 0)
      {
        WrongSNPs.push_back(this->SNPs[i].name);
      }
    }
    if(WrongSNPs.size() != 0)
    {
      misc.message << "Error: There are SNPs with same genotype for all samples or insignificant variability. Please, filter them before performing this analysis." << std::endl;
      Message message(options.outFile + ".badsnps");
      for(int i = 0; i < WrongSNPs.size(); i++)
      {
        message << WrongSNPs[i] << std::endl;
      }
      misc.message << "The SNPs are stored in the file [ " << options.outFile + ".badsnps" << " ]." << std::endl;
      misc.error("Error: There are SNPs with same genotype for all samples or insignificant variability. Please, filter them before performing this analysis.", 0);
    }
  }
  
  this->genotypes->scatterVector(vSnpsP, row);
  #pragma omp parallel for
  for(int r = 0; r<this->genotypes->nRows; r++)
  {
    for(int c = 0; c<this->genotypes->nCols; c++)
    {
      int temp = int(this->genotypes->m[c*this->genotypes->nRows + r]!=0.);
      this->genotypes->m[c*this->genotypes->nRows + r] -= double(temp)*(1. + this->genotypes->v[r]); //Substracts one to pass from our coding to the number of reference alleles. allele2 used as the reference allele.
    }
  }
  
  this->genotypes->scatterVector(vSnpsStdDev, row);
  #pragma omp parallel for
  for(int r = 0; r<this->genotypes->nRows; r++)
  {
    for(int c = 0; c<this->genotypes->nCols; c++)
    {
      this->genotypes->m[c*this->genotypes->nRows + r] /= this->genotypes->v[r];
    }
  }
  
  delete [] this->genotypes->v;
  this->genotypes->v = NULL;
  if(communicator->mpiRoot)
  {
    delete [] vSnpsP;
    delete [] vSnpsStdDev;
  }
}

void Genotype::filterSNPsAndIndividuals(std::vector<std::string> & unorderedKeepSNPIds, std::vector<std::string> & unorderedKeepIndividualIds, bool filterMissings, Genotype *newGenotype)
{
  int *keepIndividualIndxs;
  int *keepSNPIndxs;
  
  std::vector<std::string> keepSNPIds;
  std::vector<std::string> keepIndividualIds;
  
  if(this->genotypes == NULL)
  {
    misc.error("Error: An internal error has happened. Empty genotype data cannot be filtered.", 0);
  }
  
  if(this->normalized == false)
  {
    //misc.error("Error: An internal error was happened. Filtering of unnormalized genotypes still not implemented.", 0);
  }
  if( filterMissings == true && this->missings == NULL )
  {
    misc.error("Error: An internal error was happened. Missing information is not present. It can not be filtered.", 0);
  }
  
  //We filter in place?
  Genotype * filteredGenotype;
  if(newGenotype == NULL)
  {
    filteredGenotype = this;
    
    filteredGenotype->groupedBy = ungrouped;
    filteredGenotype->groupedSNPs.clear();
  }
  else
  {
    filteredGenotype = newGenotype;
    filteredGenotype->normalized = this->normalized;
    
    filteredGenotype->groupedBy = ungrouped;
    filteredGenotype->groupedSNPs.clear();
  }
  
  if(communicator->mpiRoot)
  {
    if(unorderedKeepSNPIds.size() == 0)
    {
      misc.error("Error: Error when filtering genotypes. The resultant number of SNPs is empty.", 0);
    }
    if(unorderedKeepIndividualIds.size() == 0)
    {
      misc.error("Error: Error when filtering genotypes. The resultant number of individuals is empty.", 0);
    }
  }
  
  if(communicator->mpiRoot) //Order unorderedKeepSNPIds and unorderedKeepIndividualIds in the same order that appear here.
  {
    keepSNPIds = orderVectorAsTemplate(this->SNPIds, unorderedKeepSNPIds);
    keepIndividualIds = orderVectorAsTemplate(this->individualIds, unorderedKeepIndividualIds);
    if(keepSNPIds.size() != unorderedKeepSNPIds.size())
    {
      misc.error("Error: Searching SNPs in genotypes that are not present in the BIM file.", 0);
    }
    if(keepIndividualIds.size() != unorderedKeepIndividualIds.size())
    {
      misc.error("Error: Searching individuals in genotypes that are not present in the FAM file.", 0);
    }
  }
  
  //If the individuals and SNPs not changed, return.
  int noChangesFlag = 0;
  if(communicator->mpiRoot)
  {
    noChangesFlag = ( ( (this->individualIds == keepIndividualIds) && (this->SNPIds == keepSNPIds) )?1:0 );
  }
  communicator->broadcast(&noChangesFlag, 1);
  if(noChangesFlag == 1 && newGenotype == NULL)
  {
    return;
  }
  
  if( this->spaceChangedKernelName != "" && noChangesFlag != 1 )
  {
    misc.error("Error: An internal error has happened. Operation not permitted on a transformed genotype.", 0);
  }
  if( newGenotype != NULL )
  {
    newGenotype->spaceChangedKernelName = this->spaceChangedKernelName;
  }
  
  if(communicator->mpiRoot)
  {
    //Individuals
    std::map<std::string, int> newIndividualIdsIdx;
    std::vector<Individual> newIndividuals;
    
    keepIndividualIndxs = new int [keepIndividualIds.size()];
   
    int previousOldIndividualIdx = -1;
    
    for(int r=0; r<keepIndividualIds.size(); r++) //Is it better access rows contiguously and search in map each step? Or inverse?
    {
      if(this->individualIdsIdx.count(keepIndividualIds[r]) == 0)
      {
        misc.error("Error: An error was happened. The individual '" + keepIndividualIds[r]  + "' is not in the genotypes file.", 0);
      }
      int oldIdx = this->individualIdsIdx[keepIndividualIds[r]];
      if(oldIdx<=previousOldIndividualIdx)
      {
        misc.error("Error: An internal error was happened. The individuals passed to filter the genotypes are not in the proper order.", 0);
      }
      keepIndividualIndxs[r] = oldIdx;
      newIndividualIdsIdx[keepIndividualIds[r]] = r;
      newIndividuals.push_back(this->individuals[oldIdx]);
      previousOldIndividualIdx = oldIdx;
    }
    filteredGenotype->nIndividuals = keepIndividualIds.size();
    filteredGenotype->individualIds = keepIndividualIds;
    filteredGenotype->individualIdsIdx = newIndividualIdsIdx;
    filteredGenotype->individuals = newIndividuals;
    
    //SNPs
    std::map<std::string, int> newSNPIdsIdx;
    std::vector<SNP> newSNPs;
    
    keepSNPIndxs = new int [keepSNPIds.size()];
    
    int previousOldSNPIdx = -1;
    
    for(int r=0; r<keepSNPIds.size(); r++) //Is it better access rows contiguously and search in map each step? Or inverse?
    {
      if(this->SNPIdsIdx.count(keepSNPIds[r]) == 0)
      {
        misc.error("Error: An error was happened. The SNP '" + keepSNPIds[r]  + "' is not in the genotypes file.", 0);
      }
      int oldIdx = this->SNPIdsIdx[keepSNPIds[r]];
      if(oldIdx<=previousOldSNPIdx)
      {
        misc.error("Error: An internal error was happened. The SNPs passed to filter the genotypes matrix are not in the proper order.", 0);
      }
      keepSNPIndxs[r] = oldIdx;
      newSNPIdsIdx[keepSNPIds[r]] = r;
      newSNPs.push_back(this->SNPs[oldIdx]);
      previousOldSNPIdx = oldIdx;
    }
    filteredGenotype->nSNPs = keepSNPIds.size();
    filteredGenotype->SNPIds = keepSNPIds;
    filteredGenotype->SNPIdsIdx = newSNPIdsIdx;
    filteredGenotype->SNPs = newSNPs;
  }
  
  communicator->broadcast(&filteredGenotype->nIndividuals, 1);
  communicator->broadcast(&filteredGenotype->nSNPs, 1);
  if(!communicator->mpiRoot)
  {
    keepIndividualIndxs = new int [filteredGenotype->nIndividuals];
    keepSNPIndxs = new int [filteredGenotype->nSNPs];
  }
  communicator->broadcast(keepIndividualIndxs, filteredGenotype->nIndividuals);
  communicator->broadcast(keepSNPIndxs, filteredGenotype->nSNPs);
  
  Matrix * newGenotypes = new Matrix(cyclicDistribution);
  this->genotypes->filterRowsAndColumns(newGenotypes, keepSNPIndxs, filteredGenotype->nSNPs, keepIndividualIndxs, filteredGenotype->nIndividuals);
  delete filteredGenotype->genotypes;
  filteredGenotype->genotypes = newGenotypes;
  
  if(filterMissings)
  {
    Matrix * newMissings = new Matrix(cyclicDistribution);
    this->missings->filterRowsAndColumns(newMissings, keepSNPIndxs, filteredGenotype->nSNPs, keepIndividualIndxs, filteredGenotype->nIndividuals);
    delete filteredGenotype->missings;
    filteredGenotype->missings = newMissings;
  }
  else
  {
    delete filteredGenotype->missings;
    filteredGenotype->missings = NULL;
  }
  
  delete [] keepIndividualIndxs;
  delete [] keepSNPIndxs;
}

void Genotype::appendGenotype(Genotype * newGenotype)
{
  if(this->genotypes == NULL || newGenotype->genotypes == NULL || this->missings == NULL || newGenotype->missings == NULL)
  {
    misc.error("Error: An internal error has happened. Empty genotype data cannot be joined.", 0);
  }
  
  if( this->spaceChangedKernelName != "" )
  {
    misc.error("Error: An internal error has happened. Operation not permitted on a transformed genotype.", 0);
  }
  
  if( this->normalized == true || newGenotype->normalized == true )
  {
    misc.error("Error: An internal error was happened. Normalized genotypes can not be joined.", 0);
  }
  this->groupedSNPs.clear();
  this->groupedBy = ungrouped;
  
  std::vector<std::string> SNPsIntersection = orderVectorAsTemplate(this->SNPIds, newGenotype->SNPIds);
  std::vector<std::string> individualsIntersection = orderVectorAsTemplate(this->individualIds, newGenotype->individualIds);
  if(SNPsIntersection != this->SNPIds && individualsIntersection != this->individualIds)
  {
    misc.error("Error: When joining genotypes from multiple files, one of these conditions must be met: All files have same SNPs but different individuals or all files have different SNPs but same individuals. There are at least two files that have different SNPs and different individuals.", 0);
  }
  
  int method;
  if(communicator->mpiRoot)
  {
    method = -1;
    if( this->SNPIds == newGenotype->SNPIds && individualsIntersection.size() == 0 )
    {
      method = 1;
    }
    else if( SNPsIntersection.size() == 0 && this->individualIds == newGenotype->individualIds )
    {
      method = 2;
    }
  }
  communicator->broadcast(&method);
  
  if( method == 1 )
  {
    Matrix gTemp(this->genotypes);
    Matrix mTemp(this->missings);
    this->genotypes->joinMatricesHorizontally(&gTemp, newGenotype->genotypes);
    this->missings->joinMatricesHorizontally(&mTemp, newGenotype->missings);
    
    if(communicator->mpiRoot)
    {
      //Join individual data
      for(int i = 0; i < newGenotype->individuals.size(); i++)
      {
        this->individuals.push_back(newGenotype->individuals[i]);
        std::string newIndividualId = newGenotype->individualIds[i];
        this->individualIds.push_back(newIndividualId);
        this->individualIdsIdx[newIndividualId] = newGenotype->individualIdsIdx[newIndividualId] + this->nIndividuals;
      }
      this->nIndividuals = this->individuals.size();
      if(nIndividuals != this->individualIds.size() || nIndividuals != this->individualIdsIdx.size())
      {
        misc.error("Error: An internal error was happened when joining two genotypes.", 0);
      }
      
      //Check agreement in SNP data
      for(int i = 0; i < newGenotype->SNPs.size(); i++)
      {
        if( newGenotype->SNPs[i].allele1 != this->SNPs[i].allele1 || newGenotype->SNPs[i].allele2 != this->SNPs[i].allele2 )
        {
          misc.error("Error: An error was happened when joining two genotypes. Alleles of SNP " + newGenotype->SNPs[i].name + " do not agree between two files.", 0);
        }
        if( newGenotype->SNPs[i].chr != this->SNPs[i].chr )
        {
          misc.error("Error: An error was happened when joining two genotypes. Chromosome of SNP " + newGenotype->SNPs[i].name + " do not agree between two files.", 0);
        }
        if( newGenotype->SNPs[i].position != this->SNPs[i].position )
        {
          misc.error("Error: An error was happened when joining two genotypes. Position of SNP " + newGenotype->SNPs[i].name + " do not agree between two files.", 0);
        }
      }
      
      //Update SNPs data
      for(int i = 0; i < newGenotype->SNPs.size(); i++)
      {
        this->SNPs[i].nNonMissing += newGenotype->SNPs[i].nNonMissing;
        for(int j=0; j<4; j++)
        {
          this->SNPs[i].frequencies[j] += newGenotype->SNPs[i].frequencies[j];
        }
        this->SNPs[i].p1 = ( 2.*double(this->SNPs[i].frequencies[1]) + double(this->SNPs[i].frequencies[2]) )/(2.*double(this->SNPs[i].nNonMissing));
        this->SNPs[i].p2 = ( 2.*double(this->SNPs[i].frequencies[3]) + double(this->SNPs[i].frequencies[2]) )/(2.*double(this->SNPs[i].nNonMissing));
        this->SNPs[i].standardDev = sqrt(2*this->SNPs[i].p1*(1.-this->SNPs[i].p1));
      }
    }
  }
  else if( method == 2 )
  {
    Matrix gTemp(this->genotypes);
    Matrix mTemp(this->missings);
    this->genotypes->joinMatricesVertically(&gTemp, newGenotype->genotypes);
    this->missings->joinMatricesVertically(&mTemp, newGenotype->missings);
    
    if(communicator->mpiRoot)
    {
      //Join SNP data
      for(int i = 0; i < newGenotype->SNPs.size(); i++)
      {
        this->SNPs.push_back(newGenotype->SNPs[i]);
        for(int ttii = 0; ttii < 4; ttii++) //Just check the copy of the constant array.
        {
          if(this->SNPs.back().frequencies[ttii] != newGenotype->SNPs[i].frequencies[ttii]) { misc.error("Error: An internal error was happened.", 0); }
        }
        std::string newSNPId = newGenotype->SNPIds[i];
        this->SNPIds.push_back(newSNPId);
        this->SNPIdsIdx[newSNPId] = newGenotype->SNPIdsIdx[newSNPId] + this->nSNPs;
      }
      this->nSNPs = this->SNPs.size();
      if(nSNPs != this->SNPIds.size() || nSNPs != this->SNPIdsIdx.size())
      {
        misc.error("Error: An internal error was happened when joining two genotypes.", 0);
      }
      
      //Check agreement in individual data
      for(int i = 0; i < newGenotype->individuals.size(); i++)
      {
        if(this->individuals[i].phenotype != newGenotype->individuals[i].phenotype)
        {
          misc.error("Error: The phenotypes when joining two genotypes, do not agree.", 0);
        }
      }
    }
  }
  else
  {
    misc.error("Error: When joining genotypes from multiple files, one of these conditions must be met: All files have same SNPs but different individuals or all files have different SNPs but same individuals. The SNPs and Individuals must be in the same order. There are at least two files that share some individuals and some SNPs or the order between files is not the same.", 0);
  }
  
  communicator->broadcast(&this->nSNPs);
  communicator->broadcast(&this->nIndividuals);
}

void Genotype::groupSNPs(GroupBy groupBy, std::vector<std::string> SNPIdsSubset)
{
  if(groupBy == byPosition)
  {
    misc.message << "Grouping SNPs in regions by position..." << std::endl;
    groupSNPsByPosition();
  }
  else if(groupBy == byGene || groupBy == byGroup)
  {
    misc.message << "Grouping SNPs in regions by group..." << std::endl;
    groupSNPsByGene();
  }
  else if(groupBy == byOrderedFixedSize)
  {
    misc.message << "Grouping SNPs in regions by ordered groups of size " << options.fixedGroupSize << "..." << std::endl;
    groupSNPsByOrderedFixedSize(options.fixedGroupSize, SNPIdsSubset);
  }
  else if(groupBy == byAll)
  {
    this->groupedSNPs.clear();
    this->groupedBy = byAll;
    std::set<std::string> temp(this->SNPIds.begin(), this->SNPIds.end());
    groupedSNPs["all"] = temp;
  }
  else
  {
    this->groupedSNPs.clear();
    this->groupedBy = ungrouped;
    return;
  }
  misc.message << "Created " << this->groupedSNPs.size() << " SNP group(s)" << std::endl;

  //Remove groups with not enough SNPs inside.
  std::vector<std::string> groupsToErase;
  for(std::map<std::string, std::set<std::string> >::iterator it = this->groupedSNPs.begin(); it != this->groupedSNPs.end(); ++it)
  {
    if(it->second.size() < options.minSNPsInRegions)
    {
      groupsToErase.push_back(it->first);
    }
  }
  if(groupsToErase.size() != 0)
  {
    Message message(options.outFile + ".filteredregions");
    for(int i = 0; i < groupsToErase.size(); i++)
    {
      message << groupsToErase[i] << " " << this->groupedSNPs[groupsToErase[i]].size();
      this->groupedSNPs.erase( groupsToErase[i] );
    }
    misc.message << groupsToErase.size() << " regions with less than " << options.minSNPsInRegions << " SNPs have been filtered. Regions saved on file [ " << options.outFile << ".filteredregions ]."<< std::endl;
  }
}

void Genotype::groupSNPsByPosition()
{
  this->groupedSNPs.clear();
  this->groupedBy = byPosition;
  
  if(communicator->mpiRoot)
  {
    int initialShift = options.regionSize; //Distance between the begining of the chromosome and the end of the first region.
    int regionEndDistance = options.regionSize - options.regionOverlapping; //Distance between the beginings of each region.
    int regionStartDistance = regionEndDistance; //Distance between the ends of each region.
    int regionSize = options.regionSize; //Region size.
    
    
    std::map<std::string, std::map<int, std::vector <std::string> > > chrGroups; //Group SNPs by position in each chromosome.
    for(int i = 0; i<this->nSNPs; i++)
    {
      std::string chr = SNPs[i].chr;
      std::string name = SNPs[i].name;
      int position = SNPs[i].position;
      
      int lastRegionStartBefore = position/regionStartDistance;  //The last region that begins before SNP position.
      int firstRegionEndAfter;  //The first region that ends before the snp position.
      if(position < initialShift)
      {
        firstRegionEndAfter = 0;
      }
      else
      {
        firstRegionEndAfter = ((position - initialShift)/regionEndDistance) + 1;
      }
      for(int j = firstRegionEndAfter; j <= lastRegionStartBefore; j++)
      {
        chrGroups[chr][j].push_back(name);
      }
    }
    
    //Join groups of all chromosomes together
    std::map<std::string, std::map<int, std::vector <std::string> > >::iterator it;
    for(it = chrGroups.begin(); it != chrGroups.end(); ++it)
    {
      std::string chr = it->first;
      std::map<int, std::vector <std::string> > groupsInChromosome = it->second;
      std::map<int, std::vector <std::string> >::iterator it2;
      for(it2 = groupsInChromosome.begin(); it2 != groupsInChromosome.end(); ++it2)
      {
        int groupId = it2->first;
        std::stringstream groupName;
        groupName << chr << "_" << groupId*regionStartDistance << "-" << groupId*regionStartDistance + regionSize;
        this->groupedSNPs[groupName.str()] = std::set<std::string>(it2->second.begin(), it2->second.end());
      }
    }
  }
}

void Genotype::groupSNPsByGene()
{
  this->groupedSNPs.clear();
  this->groupedBy = byGene;
  
  std::string f = options.regionsFile;
  std::ifstream file;
  std::string line;
  
  if(communicator->mpiRoot)
  {
    misc.message << "Reading SNP grouping data from file [ " << f << " ] ..." << std::endl;
    
    misc.checkFileExists(f);
    file.open(f.c_str());
    
    while(getline(file,line))
    {
      if(!file)
      {
        break;
      }
      
      std::istringstream sstemp(line);
      
      std::string snp;
      sstemp >> snp;
      if(this->SNPIdsIdx.count(snp) == 0)
      {
        //misc.error("Error: There are a SNP in the regions file not present in the genotypes file.", 0);
        continue;
      }

      std::string group;
      while( sstemp >> group )
      {
        if(groupedSNPs[group].find(snp) != groupedSNPs[group].end())
        {
          misc.message << "Warning: The group " << group << " is defined more than once for SNP " << snp << "." << std::endl;
        }
        this->groupedSNPs[group].insert(snp);
      }
    }
    file.close();
  }
}

void Genotype::groupSNPsByOrderedFixedSize(int groupsSize, std::vector<std::string> SNPIdsSubsetArg)
{
  this->groupedSNPs.clear();
  this->groupedBy = byOrderedFixedSize;
  
  if(groupsSize < 1)
  {
    misc.error("Error: An internal error was happened when grouping SNPs using fixed group sizes. Invalid group size.", 0);
  }
  
  if(communicator->mpiRoot)
  {
    std::vector<std::string> SNPIdsSubset;
    
    if(SNPIdsSubsetArg.size() == 0)
    {
      SNPIdsSubset = this->SNPIds;
    }
    else
    {
      SNPIdsSubset = SNPIdsSubsetArg;
    }
    
    std::map<std::string, std::map<int, std::vector <std::string> > > chrOrdered;
    
    //Organize SNPs in chromosomes and order in a map within each chromosome.
    //Accepts SNPs within same chromosome and same coordinates. In this case, order of SNPs with same coordinates is indetermined.
    for(int i = 0; i < SNPIdsSubset.size(); i++)
    {
      if(this->SNPIdsIdx.count( SNPIdsSubset[i] ) == 0)
      {
        continue;
      }
      int idx = this->SNPIdsIdx[ SNPIdsSubset[i] ];
      
      std::string chr = SNPs[idx].chr;
      std::string name = SNPs[idx].name;
      int position = SNPs[idx].position;
      
      chrOrdered[chr][position].push_back(name);
    }
    
    //Make the groups
    std::map<std::string, std::map<int, std::vector <std::string> > >::iterator it;
    std::vector<std::string> currentGroup;
    int groupId = 0;
    for(it = chrOrdered.begin(); it != chrOrdered.end(); ++it)
    {
      std::string chr = it->first;
      std::map<int, std::vector <std::string> > orderedSNPsInChromosome = it->second;
      std::map<int, std::vector <std::string> >::iterator it2;
      int previousPosition = -1;
      for(it2 = orderedSNPsInChromosome.begin(); it2 != orderedSNPsInChromosome.end(); ++it2)
      {
        if(previousPosition > it2->first && it2 != orderedSNPsInChromosome.begin())
        {
          misc.error("Error: An internal error was happened when grouping SNPs using fixed group sizes. Error in ordering.", 0);
        }
        previousPosition = it2->first;
        std::vector<std::string> positionSNPs = it2->second;
        for(int i = 0; i < positionSNPs.size(); i++)
        {
          currentGroup.push_back(positionSNPs[i]);
          if(currentGroup.size() == groupsSize)
          {
            std::stringstream groupName;
            groupName << chr << "_" << groupId;
            this->groupedSNPs[groupName.str()] = std::set<std::string>(currentGroup.begin(), currentGroup.end());
            currentGroup.clear();
            groupId++;
          }
        }
      }
      if(currentGroup.size() != 0)
      {
        std::stringstream groupName;
        groupName << chr << "_" << groupId;
        this->groupedSNPs[groupName.str()] = std::set<std::string>(currentGroup.begin(), currentGroup.end());
      }
      currentGroup.clear();
      groupId = 0;
    }
  }
}

void Genotype::groupSNPsByFileOrderedWindows(int windowSize)
{
  this->groupedSNPs.clear();
  this->groupedBy = byFileOrderedWindows;
  
  if(windowSize < 1)
  {
    misc.error("Error: An internal error was happened when grouping SNPs using fixed file ordered windows. Invalid group size.", 0);
  }
  
  if(communicator->mpiRoot)
  {
    int idx = 0;
    int ig = 0;
    while(idx < this->SNPIds.size())
    {
      int istart = idx;
      int iend = idx + windowSize;
      if( iend >= this->SNPIds.size() )
      {
        iend = this->SNPIds.size();
      }
      
      std::set<std::string> gsnps;
      for(int i = istart; i < iend; i++)
      {
        gsnps.insert(this->SNPIds[i]);
        idx++;
      }
      this->groupedSNPs[ "fw" + i2s(windowSize) + "_" + i2s(ig) ] = gsnps;
      ig++;
    }
  }
}

std::vector<std::string> Genotype::getGroup(std::string group, bool writeGroup)
{
  if( this->groupedBy == ungrouped )
  {
    misc.error("Error: An internal error was happened. The SNPs within a group can not be retrieved. Groups are not defined.", 0);
  }
  if( communicator->mpiRoot == true && this->groupedSNPs.count(group) == 0 )
  {
    misc.error("Error: An internal error was happened. The group " + group + " is not a valid group.", 0);
  }
  
  std::vector<std::string> groupList;
  if(communicator->mpiRoot)
  {
    groupList.assign(this->groupedSNPs[group].begin(), this->groupedSNPs[group].end());
    
    if(writeGroup == true)
    {
      Message message(options.outFile + ".region." + group);
      for(int i = 0; i < groupList.size(); i++)
      {
        message << groupList[i] << "\n";
      }
    }
  }
  
  return groupList;
}

void Genotype::genotypeOfSNPsGroup(std::string group, Genotype *newGenotype, bool writeGroup)
{
  if(this->groupedBy == ungrouped)
  {
    misc.error("Error: An internal error was happened. The genotype of a group can not be retrieved. Groups are not defined.", 0);
  }
  if( communicator->mpiRoot == true && this->groupedSNPs.count(group) == 0 )
  {
    misc.error("Error: An internal error was happened. The group " + group + " is not a valid group.", 0);
  }
  if(newGenotype == this)
  {
    misc.error("Error: An internal error was happened. When getting a genotpe group, the resultat genotype must be different than the source genotype.", 0);
  }
  if( this->spaceChangedKernelName != "" )
  {
    misc.error("Error: An internal error has happened. Operation not permitted on a transformed genotype.", 0);
  }
  
  std::vector<std::string> SNPGroup = getGroup(group, writeGroup);
  if(communicator->mpiRoot == true && SNPGroup.size() == 0)
  {
    misc.error("Error: The group " + group + " is not a valid group. There are not SNPs inside this group.", 0);
  }
  filterSNPsAndIndividuals(SNPGroup, this->individualIds, true, newGenotype);
}

Matrix* Genotype::genotypesRedistributionToGroupedCommunicatorMatrices(Communicator * groupedCommunicator, std::map< int, std::vector<int> > & SNPidxs)
{
  if(this->genotypes == NULL)
  {
    misc.error("Error: An internal error has happened. Empty genotype data cannot be redistributed.", 0);
  }
  
  if( this->spaceChangedKernelName != "" )
  {
    misc.error("Error: An internal error has happened. Operation not permitted on a transformed genotype.", 0);
  }
  
  int isShift = 0;
  SNPidxs.clear();
  std::map<int, std::pair<int, int> > rowsOriginDestination;
  std::vector<int> nGlobalRowsInGroup;
  std::vector<int> nGlobalColsInGroup;
  
  int nSNPsBaseInGroup = this->nSNPs/groupedCommunicator->nGroups;
  int remainingSNPs = this->nSNPs - (nSNPsBaseInGroup*groupedCommunicator->nGroups);
  if(remainingSNPs < 0 || remainingSNPs >= groupedCommunicator->nGroups)
  {
    misc.error("Error: An internal error was happened. Unexpected remainder when redistributing genotypes.", 0);
  }
  for(int ig = 0; ig < groupedCommunicator->nGroups; ig++)
  {
    int nSNPsInGroup = nSNPsBaseInGroup;
    if( ig < remainingSNPs )
    {
      nSNPsInGroup++;
    }
    for(int is = 0; is < nSNPsInGroup; is++)
    {
      int globsi = is + isShift; //Global SNP index.
      rowsOriginDestination[globsi] = std::pair<int, int>(is,ig);
      SNPidxs[ig].push_back(globsi);
    }
    isShift += nSNPsInGroup;
    
    nGlobalRowsInGroup.push_back(nSNPsInGroup);
    nGlobalColsInGroup.push_back(this->nIndividuals);
  }
  if(isShift != this->nSNPs)
  {
    misc.error("Error: An internal error has happened. Not all SNPs have been redistributed.", 0);
  }
  
  std::map<int, std::pair<int, int> > colsOriginDestination;
  for(int ii = 0; ii < this->nIndividuals; ii++)
  {
    colsOriginDestination[ii] = std::pair<int, int>(ii,-1);
  }
  
  Matrix * groupedGenotypes = this->genotypes->redistributionToGroupedCommunicatorMatrices(groupedCommunicator, nGlobalRowsInGroup, nGlobalColsInGroup, rowsOriginDestination, colsOriginDestination, false);
  
  return groupedGenotypes;
}

void Genotype::changeSpace(std::string name)
{
  if(this->spaceChangedKernelName != "")
  {
    misc.error("Error: An internal error has happened. The space was already changed.", 0);
  }
  if(name == "")
  {
    misc.error("Error: An internal error has happened. Invalid kernel name.", 0);
  }
  
  this->spaceChangedKernelName = name;
}

void Genotype::printGenotype(bool showGroups, bool unnormalizeProbabilities)
{
  double *g;
  double *m;
  
  if (communicator->mpiRoot) {
    g = new double [this->nSNPs*this->nIndividuals];
    m = new double [this->nSNPs*this->nIndividuals];
  }
  
  this->genotypes->gatherMatrix(g);
  this->missings->gatherMatrix(m);
  if (communicator->mpiRoot) {
    if(showGroups)
    {
      misc.message << "********" << std::endl;
      misc.message << "Genotype grouped by: " << this->groupedBy << std::endl;
      misc.message << this->groupedSNPs.size() << " groups." << std::endl;
      misc.message << "********" << std::endl;
      
      std::map<std::string, std::set<std::string> >::iterator it;
      for(it = this->groupedSNPs.begin(); it != this->groupedSNPs.end(); ++it)
      {
        misc.message << it->first << std::endl;
        std::set<std::string>::iterator it2;
        for(it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
          misc.message << *it2 << std::endl;
        }
      }
    }
    
    misc.message << "Genotype Matrix:\n";
    misc.message << this->nSNPs << " " << this->SNPs.size() << " " << this->SNPIds.size() << " " << this->SNPIdsIdx.size() << std::endl;
    misc.message << this->nIndividuals << " " << this->individuals.size() << " " << this->individualIds.size() << " " << this->individualIdsIdx.size() << std::endl;
    
    for (int r = 0; r < this->genotypes->nGlobRows; r++)
    {
      misc.message << this->SNPIds[r] << " " << this->SNPIdsIdx[this->SNPIds[r]] << " ";
    }
    misc.message << std::endl;
    for (int c = 0; c < this->genotypes->nGlobCols; c++)
    {
      misc.message << this->individuals[c].familyID << " " << this->individuals[c].individualID << " " << this->individuals[c].paternalID << " " << this->individuals[c].maternalID << " " << this->individuals[c].sex << " " << this->individuals[c].phenotype << " " << this->individualIds[c] << " " << this->individualIdsIdx[this->individualIds[c]] << " ";
      for (int r = 0; r < this->genotypes->nGlobRows; r++)
      {
        if(this->type == GenotypeAttributes::calls)
        {
          int temp = int(*(g + this->genotypes->nGlobRows*c + r));
          if(temp == 1)
          {
            misc.message << SNPs[r].allele1 << " " << SNPs[r].allele1 << " ";
          }
          else if(temp == 2)
          {
            misc.message << SNPs[r].allele1 << " " << SNPs[r].allele2 << " ";
          }
          else if(temp == 3)
          {
            misc.message << SNPs[r].allele2 << " " << SNPs[r].allele2 << " ";
          }
          else
          {
            misc.message << "0" << " " << "0" << " ";
          }
        }
        else if(this->type == GenotypeAttributes::probabilities)
        {
          if(int(*(m + this->genotypes->nGlobRows*c + r)) == 1)
          {
            if(unnormalizeProbabilities == false)
            {
              misc.message << *(g + this->genotypes->nGlobRows*c + r) << " ";
            }
            else
            {
              misc.message << ((*(g + this->genotypes->nGlobRows*c + r))*this->SNPs[r].standardDev) + 2.*this->SNPs[r].p2 << " ";
            }
          }
          else
          {
            misc.message << "-1.0000000 ";
          }
        }
        else
        {
          misc.error("Error: Invalid type for genotype printing.", 0);
        }
      }
      misc.message << "\n";
    }
    misc.message << std::endl;
  }
  
  if (communicator->mpiRoot) {
    misc.message << "Missings Matrix:\n";
    for (int c = 0; c < this->genotypes->nGlobCols; c++)
    {
      misc.message << this->individuals[c].familyID << " " << this->individuals[c].individualID << " " << this->individuals[c].paternalID << " " << this->individuals[c].maternalID << " " << this->individuals[c].sex << " " << this->individuals[c].phenotype << " ";
      for (int r = 0; r < this->genotypes->nGlobRows; r++)
      {
	int temp = int(*(m + this->genotypes->nGlobRows*c + r));
	if(temp == 1)
	{
	  misc.message << " " << "-" << " " << " ";
	}
	else
	{
	  misc.message << "0" << "-" << "0" << " ";
	}
      }
      misc.message << std::endl;
    }
    misc.message << std::endl;
  }
  
  if (communicator->mpiRoot) {
    delete [] g;
    delete [] m;
  }
}

void Genotype::printGenotypeT(bool unnormalizeProbabilities)
{
  double *g;
  double *m;
  
  if (communicator->mpiRoot) {
    g = new double [this->nSNPs*this->nIndividuals];
    m = new double [this->nSNPs*this->nIndividuals];
  }
  
  this->genotypes->gatherMatrix(g);
  this->missings->gatherMatrix(m);
  if (communicator->mpiRoot) {
    
    misc.message << "Genotype Matrix:\n";
    misc.message << this->nSNPs << " " << this->SNPs.size() << " " << this->SNPIds.size() << " " << this->SNPIdsIdx.size() << std::endl;
    misc.message << this->nIndividuals << " " << this->individuals.size() << " " << this->individualIds.size() << " " << this->individualIdsIdx.size() << std::endl;
    
    for (int c = 0; c < this->genotypes->nGlobCols; c++)
    {
      misc.message << this->individualIds[c] << " ";
    }
    misc.message << std::endl;
    for (int r = 0; r < this->genotypes->nGlobRows; r++)
    {
      misc.message << this->SNPIds[r] << " ";
      for (int c = 0; c < this->genotypes->nGlobCols; c++)
      {
        if(this->type == GenotypeAttributes::calls)
        {
          int temp = int(*(g + this->genotypes->nGlobRows*c + r));
          if(temp == 1)
          {
            misc.message << SNPs[r].allele1 << " " << SNPs[r].allele1 << " ";
          }
          else if(temp == 2)
          {
            misc.message << SNPs[r].allele1 << " " << SNPs[r].allele2 << " ";
          }
          else if(temp == 3)
          {
            misc.message << SNPs[r].allele2 << " " << SNPs[r].allele2 << " ";
          }
          else
          {
            misc.message << "0" << " " << "0" << " ";
          }
        }
        else if(this->type == GenotypeAttributes::probabilities)
        {
          if(int(*(m + this->genotypes->nGlobRows*c + r)) == 1)
          {
            if(unnormalizeProbabilities == false)
            {
              misc.message << *(g + this->genotypes->nGlobRows*c + r) << " ";
            }
            else
            {
              misc.message << ((*(g + this->genotypes->nGlobRows*c + r))*this->SNPs[r].standardDev) + 2.*this->SNPs[r].p2 << " ";
            }
          }
          else
          {
            misc.message << "-1.0000000(" << *(g + this->genotypes->nGlobRows*c + r) << ") ";
          }
        }
        else
        {
          misc.error("Error: Invalid type for genotype printing.", 0);
        }
      }
      misc.message << "\n\n";
    }
    misc.message << std::endl;
  }
  
  if (communicator->mpiRoot) {
    misc.message << "Missings Matrix:\n";
    for (int r = 0; r < this->genotypes->nGlobRows; r++)
    {
      misc.message << this->SNPIds[r] << " ";
      for (int c = 0; c < this->genotypes->nGlobCols; c++)
      {
        int temp = int(*(m + this->genotypes->nGlobRows*c + r));
        if(temp == 1)
        {
          misc.message << " " << "-" << " " << " ";
        }
        else
        {
          misc.message << "0" << "-" << "0" << " ";
        }
      }
      misc.message << std::endl;
    }
    misc.message << std::endl;
  }
  
  if (communicator->mpiRoot) {
    delete [] g;
    delete [] m;
  }
}


/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
// Old functions
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/


void Genotype::readBEDFileOld(std::string f, std::vector<int> & SNPsIdxsToKeep, std::vector<int> & individualIdxsToKeep)
{/*
  std::ifstream file;
  unsigned char * header; //Temporal array where header information will be stored.
  unsigned char * memBlock; //Temporal array where a file row will be stored
  double * genotypeMatrixBlock; //Temporal matrix where a block of genotypes will be stored
  double * missingMatrixBlock; //Temporal matrix where a block of missing information will be stored
  int nrf; //Number of bytes that store 4 genotypes
  int rest; //Number of genotypes in the last byte.
  int inc; //1 if rest!=0, 0 otherwise

  misc.setGetElapsedTime("GenotypeLoad");
  misc.message << "Reading genotype data from file [ "<< f <<" ] ..." << std::endl;
  
  this->genotypes->initParameters(this->nSNPs, this->nIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockRows);
  this->missings->initParameters(this->nSNPs, this->nIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockRows);
  
  //Load the BED header. Check if the file is valid BED file and if it is in SNP-Major mode
  if(communicator->myCol == 0)
  {
    misc.checkFileExists(f);
    file.open(f.c_str(),std::ios::binary);
    
    //Checking the BED header.
    header=new unsigned char [3];
    file.read((char*)( &header[0]), 3);
    
    //Is a BED file?
    if(header[0] != 0x6C || header[1]!=0x1B)
    {
      misc.error("Error: The data file is not a valid binary BED > v1.0", 0);
    }
    
    //The file is in SNP-Major mode?
    if(header[2]==0x1)
    {
      nrf = this->nIndividuals/4;
      rest = this->nIndividuals%4;
      if(((nrf*4)+rest)!=this->nIndividuals)
      {
        misc.error("Error: An internal error was happened.", 0);
      }
      misc.message << "Reading the binary file in SNP-major mode." << std::endl;
    }
    else if(header[2]==0x0)
    {
      misc.error("Error: The BED file is in Individual-major mode. This mode is not currently supported. Please, read PLINK manual.", 0);
    }
    else
    {
      misc.error("Error: The data file is not a valid binary BED.\nIt is not specified whether the file is in SNP-major o Individual-major mode.", 0);
    }
    delete [] header;
    
    inc=0;
    if(rest!=0)
    {
      inc = 1;
    }
    memBlock = new unsigned char [nrf + inc];
    genotypeMatrixBlock = new double [this->nIndividuals*this->genotypes->nBlockRows];
    missingMatrixBlock = new double [this->nIndividuals*this->genotypes->nBlockRows];
    
    if(memBlock==NULL || genotypeMatrixBlock==NULL || missingMatrixBlock==NULL)
    {
      misc.error("Error: There are not enough memory for reading the binary file.", 0);
    }
  }
  
  //Resort the individual and SNP indexs to be loaded.
  std::vector<int> keepSNPsMask;
  if(SNPsIdxsToKeep.size() != 0)
  {
    keepSNPsMask = std::vector<int>();
  }
  
  //Read genotype data
  std::vector<int> loadedSNPsByThisProcess;
  loadedSNPsByThisProcess.clear();
  if(communicator->mpiRoot == false)
  {
    this->SNPs = std::vector<SNP>(this->nSNPs, SNP());
  }
  for(int i=0; i<this->nSNPs; i += this->genotypes->nBlockRows)
  {
    //Which process will read and scatter the data?
    int sendProcessRow;
    int sendProcessCol;
    sendProcessRow = i/this->genotypes->nBlockRows;
    sendProcessRow = sendProcessRow % communicator->nProcRows;
    sendProcessCol = 0;
    
    //Read a block of genotypes of size this->genotypes->nBlockRows*this->nIndividuals
    if(communicator->myCol == 0)
    {
      int nFileRows = this->genotypes->nBlockRows;
      if((this->nSNPs-i) < nFileRows)
      {
        nFileRows = this->nSNPs-i;
      }
      if(communicator->myRow == sendProcessRow)
      {
        for(int j=0; j<nFileRows; j++)
        {
          //std::cout << i << " " << j << " " << i+j << std::endl; std::cout.flush();
          file.read((char*)(&memBlock[0]), nrf + inc);
          //parseSNP(nrf, rest, memBlock, &(genotypeMatrixBlock[j*this->nIndividuals]), &(genotypeMissingMatrixBlock[j*this->nIndividuals]), &(SNPs[i+j])); //Change from row major to col major.
          parseSNP(nrf, rest, memBlock, &(genotypeMatrixBlock[j]), &(missingMatrixBlock[j]), &(SNPs[i+j])); //Change from row major to col major.
          loadedSNPsByThisProcess.push_back(i+j);
        }
      }
      else if(communicator->myRow != sendProcessRow)
      {
        file.seekg(nFileRows*(nrf + inc), std::ios_base::cur);
      }
    }
    
    //Distribute the read block of genotypes between the genotypes matrix
    for(int k=0; k<this->nIndividuals; k+=this->genotypes->nBlockCols)
    {
      this->genotypes->scatterBlock((genotypeMatrixBlock+k*this->genotypes->nBlockRows), i, k, this->genotypes->nBlockRows, sendProcessRow, sendProcessCol);
      this->missings->scatterBlock((missingMatrixBlock+k*this->genotypes->nBlockRows), i, k, this->genotypes->nBlockRows, sendProcessRow, sendProcessCol);
    }
    
  }
  
  //Release temporal allocated arrays
  if(communicator->myCol == 0)
  {
    delete [] memBlock;
    delete [] genotypeMatrixBlock;
    delete [] missingMatrixBlock;
    
    file.close();
  }
  
  gatherSNPData(loadedSNPsByThisProcess);

  Cblacs_barrier(communicator->context, "All");
  
  misc.message << "Genotype loaded after " << misc.setGetElapsedTime("GenotypeLoad", true) << std::endl;
*/}