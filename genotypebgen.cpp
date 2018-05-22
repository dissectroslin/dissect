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
// The structs ProbabilitiesSetter and ReadBGENData are developed from an example
// created by Gavin Band with Copyright:
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifdef BGEN

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <memory>
#include <vector>
#include <utility>
#include "genfile/bgen/bgen.hpp"
#include "genotype.h"
#include "communicator.h"
#include "misc.h"
#include "auxiliar.h"
#include "options.h"

struct ProbabilitiesSetter
{
private:
  int layout;
  
  double* rowGenotypes;
  int* rowMissings;
  
  double layout1Missing;
  
  int nTotalSamples;
  
  int indFileIdx;
  int snpIdx;
  
  std::vector<int> * indFileIdxsMaskToKeep;
  bool avoidThis;
  
  int indIdx;
  int nEntries;
  
  double *mean;
  double *std;
  
  bool *withMissings;

public:
  ProbabilitiesSetter( int bgenlayout, double* buffer, int* mbuffer, double * rmean, double * rstd, bool * rwithMissings, int nKeepSamples, std::vector<int> * individualIdxsMaskToKeep )
  {
    this->layout = bgenlayout;
    
    this->rowGenotypes = buffer;
    this->rowMissings = mbuffer;
    
    this->nTotalSamples = nKeepSamples;
    
    this->indFileIdx = -1;
    this->snpIdx = 0;
    
    this->indFileIdxsMaskToKeep = individualIdxsMaskToKeep;
    this->avoidThis = false;
    
    this->indIdx = -1;
    this->nEntries = 0;
    
    this->mean = rmean;
    this->std = rstd;
    *this->mean = 0;
    *this->std = 0;
    
    this->withMissings = rwithMissings;
    *this->withMissings = false;
  }
  
  // Called once allowing us to set storage.
  void initialise( std::size_t sampleNumber, std::size_t allelesNumber ) {
    //m_result->clear();
    //m_result->resize( number_of_samples );
    if( sampleNumber < this->nTotalSamples)
    {
      misc.error("Error: An internal error has happened. Unexpected small number of samples.", 0);
    }
    if(allelesNumber != 2)
    {
      misc.error("Error: The bgen file contains non-biallelic genotypes. DISSECT currently can only work with biallelic genotypes.", 0);
    }
  }
  
  // If present with this signature, called once after initialise()
  // to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
  // This enables us to set up storage for the data ahead of time.
  void set_min_max_ploidy( uint32_t minPloidy, uint32_t maxPloidy, uint32_t minEntries, uint32_t maxEntries ) {
    /*for( std::size_t i = 0; i < m_result->size(); ++i ) {
      m_result->at( i ).reserve( max_entries );
    }*/
    if( minPloidy != 2 || maxPloidy != 2 )
    {
      misc.error("Error: The bgen file contains non-diploid genotypes. DISSECT currently can only work with diploid genotypes.", 0);
    }
    if( minEntries != maxEntries )
    {
      misc.error("Error: Unexpected different number of entries for different individuals.", 0);
    }
    this->nEntries = minEntries;
  }
  
  bool set_sample( std::size_t i )
  {
    if(this->indFileIdx != -1 && this->snpIdx != 3 && this->avoidThis == false)
    {
      misc.error("Error: An internal error has happened. Unexpected number of individual genotypes when reading a bgen file.", 0);
    }
    if((this->indFileIdx + 1) != i) //From the documentation, I think i has to be consecutive, but check just in case.
    {
      misc.error("Error: An internal error has happened. Unexpected sequence of indices when reading bgen file.", 0);
    }
    this->indFileIdx = i;
    
    if((*this->indFileIdxsMaskToKeep)[this->indFileIdx] == 1)
    {
      this->avoidThis = false;
      this->indIdx++;
      this->rowMissings[this->indIdx] = 1;
    }
    else
    {
      this->avoidThis = true;
    }
    
    this->snpIdx = 0;
    
    this->layout1Missing = 0;
    
    return true;
  }
  
  // Called once per sample to set the number of probabilities that are present.
  void set_number_of_entries(std::size_t ploidy, std::size_t numberOfEntries, genfile::OrderType orderType, genfile::ValueType valueType)
  {
    /*if( (this->nEntries != 3 && orderType == genfile::ePerUnorderedGenotype) || 
         (this->nEntries != 4 && orderType == genfile::ePerPhasedHaplotypePerAllele)
         || valueType == genfile::eProbability ) 
    ) 
    {
      misc.error("Error: Unexpected number of probabilities when reading bgen file.", 0);
    }*/
    if( this->nEntries != 3 || numberOfEntries != 3 || orderType != genfile::ePerUnorderedGenotype || valueType != genfile::eProbability ) 
    {
      misc.error("Error: Currently can only load unphased genotypes. If you are interested on loading phased genotypes, please, contact us.", 0);
    }
  }
  
  // Called once for each genotype (or haplotype) probability per sample.
  void set_value( uint32_t i, double value )
  {
    if(avoidThis == false)
    {
      if(this->snpIdx != i) //From the documentation, I am not sure this has to be consecutive, check just in case.
      {
        misc.error("Error: An internal error has happened. Unexpected sequence of indices when reading genotype individual of bgen file.", 0);
      }
      this->snpIdx++;
      this->layout1Missing += value;
      
      if(i == 0)
      {
      }
      else if(i == 1)
      {
        rowGenotypes[this->indIdx] = value;
      }
      else if(i == 2)
      {
        rowGenotypes[this->indIdx] += 2.*value;
        *this->mean += rowGenotypes[this->indIdx];
        
        if(this->layout1Missing == 0.) //Layout 1 codes missings with all probabilities == 0. Check this.
        {
          if(this->layout == 1)
          {
            this->rowMissings[this->indIdx] = 0;
            *this->withMissings = true;
          }
          else
          {
            misc.error("Error: The sum of all probabilities is 0. This is unexpected for a BGEN file coded using Layout 2.", 0);
          }
        }
      }
      else if(i > 2)
      {
        misc.error("Error: An internal error has happened. Unexpected index when reading individual genotype.", 0);
      }
    }
  }

  void set_value( uint32_t i, genfile::MissingValue value )
  {
    if(avoidThis == false)
    {
      if(this->snpIdx != i) //From the documentation, I am not sure this has to be consecutive, check just in case.
      {
        misc.error("Error: An internal error has happened. Unexpected sequence of indices when reading genotype individual of bgen file.", 0);
      }
      this->snpIdx++;
      
      this->rowMissings[this->indIdx] = 0;
      *this->withMissings = true;
    }
  }
  
  // If present with this signature, called once after all data has been set.
  void finalise()
  {
    if(this->nTotalSamples != (this->indIdx + 1))
    {
      misc.error("Error: An internal error has happened. Read more individuals than expected.", 0);
    }
    
    int nNonMissings = 0;
    if(*this->withMissings == true)
    {
      for(int i = 0; i < this->nTotalSamples; i++)
      {
        nNonMissings += this->rowMissings[i];
      }
    }
    else
    {
      nNonMissings = this->nTotalSamples;
    }
    
    *this->mean = *this->mean/double(nNonMissings);
    
    *this->std = 0;
    for(int i = 0; i < this->nTotalSamples; i++)
    {
      if(this->rowMissings[i] == 1)
      {
        this->rowGenotypes[i] = this->rowGenotypes[i] - *this->mean;
        *this->std += this->rowGenotypes[i]*this->rowGenotypes[i];
      }
    }
    *this->std = sqrt(*this->std/double(this->nTotalSamples));
    
    for(int i = 0; i < this->nTotalSamples; i++)
    {
      if(this->rowMissings[i] == 1)
      {
        this->rowGenotypes[i] /= *this->std;
      }
      else
      {
        this->rowGenotypes[i] = 0.;
      }
    }
  }
};


struct ReadBGENData
{
private:
  bool isActive;
  
  int layout;
  
  std::string fileName;
  std::istream * stream;
  
  genfile::bgen::Context context;
  
  uint32_t offset;
  
  enum State { fileNotOpen, fileOpen, streamReadyForVariant, streamReadyForProbs };
  State state;
  
  std::vector< std::pair<std::string, std::string> > individualRawIds;
  std::vector< int > individualBGENIdxsMaskToKeep;
  
  std::vector< genfile::byte_t > tempBuffer1;
  std::vector< genfile::byte_t > tempBuffer2;
  
  bool beQuiet;
  
public:
  ReadBGENData( std::string const& fn, std::set<std::string> & keepIndividualIds, bool active, bool quiet = false )
  {
    this->isActive = active;
    this->beQuiet = quiet;
    if(this->isActive)
    {
      this->fileName = fn;
      this->state = fileNotOpen;
      
      // Open the stream
      this->stream = new std::ifstream( (fn + ".bgen").c_str(), std::ifstream::binary );
      if( !*this->stream ) {
        misc.error( "Error: Unable to read the file [ " + fn + ".bgen ].", 0 );
      }
      this->state = fileOpen;
      
      // Read the offset, header, and sample IDs if present.
      std::vector< std::pair<std::string, std::string> > bgenIndividualIds;
      genfile::bgen::read_offset( *this->stream, &this->offset );
      genfile::bgen::read_header_block( *this->stream, &this->context );
      
      if((this->context.flags & genfile::bgen::e_Layout) == genfile::bgen::e_Layout0)
      {
        this->layout = 0;
        misc.error("Error: The BGEN file is in Layout 0. Currently DISSECT only supports Layout 1 and Layout 2.", 0);
      }
      else if((this->context.flags & genfile::bgen::e_Layout) == genfile::bgen::e_Layout1)
      {
        this->layout = 1;
        if(options.allowBGENLayout1 == false)
        {
          misc.error("Error: The BGEN file is using Layout 1. Currently DISSECT only supports Layout 1 for testing purposes. You can use the option --bgen-l1 to enable reading BGEN layout 1 files. But we recommend to update the file to Layout 2.", 0);
        }
      }
      else if((this->context.flags & genfile::bgen::e_Layout) == genfile::bgen::e_Layout2)
      {
        this->layout = 2;
      }
      else
      {
        misc.error("Error: Unrecognized BGEN layout type. Supported types are Layout 1 and Layout 2.", 0);
      }
      
      if( this->context.flags & genfile::bgen::e_SampleIdentifiers )
      {
        genfile::bgen::read_sample_identifier_block(
          *this->stream, this->context,
          [&bgenIndividualIds]( std::string id ) { bgenIndividualIds.push_back( std::pair<std::string, std::string>(id, id) ); }
        );
      }
      else
      {
        // Load ids from sample file (if exists)
        misc.message << "There are no individual identifiers in the [ " + fn + ".bgen ] file. Trying to read them from [ " + fn + ".sample file ]." << std::endl;
        bgenIndividualIds = readSampleFile(fn + ".sample");
      }
      if( this->context.number_of_samples != bgenIndividualIds.size() )
      {
        misc.error("Error: An error has happened whe reading bgen file. Unexpected number of individual Ids.", 0);
      }
      
      for (int i = 0; i < bgenIndividualIds.size(); i++)
      {
        std::string key = bgenIndividualIds[i].first + "@" + bgenIndividualIds[i].second;
        if( keepIndividualIds.size() != 0 && keepIndividualIds.find(key) == keepIndividualIds.end() )
        {
          this->individualBGENIdxsMaskToKeep.push_back(0);
        }
        else
        {
          this->individualBGENIdxsMaskToKeep.push_back(1);
          this->individualRawIds.push_back(bgenIndividualIds[i]);
        }
      }
      
      // Jump to the first variant data block.
      this->stream->seekg( this->offset + 4 );
      
      // We keep track of state (though it's not really needed for this implementation.)
      this->state = streamReadyForVariant;
    }
    else
    {
      this->stream = NULL;
    }
  }
  
  ~ReadBGENData()
  {
    if(this->stream != NULL)
    {
      delete this->stream;
    }
  }
  
  void showFileInfo( )
  {
    if(this->isActive == true && this->beQuiet == false)
    {
      misc.message << "Reading bgen file [ " << this->fileName <<" ] ("
      << ( this->context.flags & genfile::bgen::e_Layout2 ? "v1.2 layout" : "v1.1 layout" )
      << ", " << this->layout << ", "
      << ( this->context.flags & genfile::bgen::e_CompressedSNPBlocks ? "compressed" : "uncompressed" ) << ")"
      << " with " << this->context.number_of_samples << " individuals and " << this->context.number_of_variants << " variants." << std::endl;
    }
  }
  
  std::vector< std::pair<std::string, std::string> > readSampleFile(std::string fn)
  {
    std::vector< std::pair<std::string, std::string> > indIds;
    if( this->isActive == true )
    {
      misc.checkFileExists(fn);
      indIds.clear();
      std::vector< std::vector<std::string> > tempIdsTable;
      getTableFromFile(fn, tempIdsTable, 2);
      if(tempIdsTable.size() <= 2)
      {
        misc.error("Error: The file [ " + fn + " ] is empty.", 0);
      }
      if(tempIdsTable[0][0] != "ID_1" || tempIdsTable[0][1] != "ID_2" || tempIdsTable[1][0] != "0" || tempIdsTable[1][1] != "0")
      {
        misc.error("Error: Unexpected header in the file [" + fn + "].", 0);
      }
      for(int rid = 2; rid<tempIdsTable.size(); rid++)
      {
        indIds.push_back(std::pair<std::string, std::string>(tempIdsTable[rid][0], tempIdsTable[rid][1]));
      }
    }
    else
    {
      indIds.clear();
    }
    return indIds;
  }
  
  std::vector< std::pair<std::string, std::string> > getIndividualRawIds()
  {
    if(this->isActive)
    {
      return this->individualRawIds;
    }
    else
    {
      return std::vector< std::pair<std::string, std::string> >();
    }
  }
  int getNSNPs()
  {
    if(this->isActive)
    {
      return this->context.number_of_variants;
    }
    else
    {
      return 0;
    }
  }
  
  bool readVariantInfo(std::string* chromosome, uint32_t* position, std::string* rsid, std::vector< std::string >* alleles)
  {
    if(this->isActive)
    {
      if( this->state != streamReadyForVariant )
      {
        misc.error("Error: An internal error has happened.", 0);
      }
      std::string dummy;
      
      if(
        genfile::bgen::read_snp_identifying_data(
          *this->stream, this->context,
          &dummy, rsid, chromosome, position,
          [&alleles]( std::size_t n ) { alleles->resize( n ); },
                                                [&alleles]( std::size_t i, std::string const& allele ) { alleles->at(i) = allele; }
        )
      )
      {
        this->state = streamReadyForProbs;
        
        return true;
      }
      else
      {
        return false; //Data has not been properly read, this indicates end of file.
      }
    }
  }
  
  void readProbabilities( double* probs, int *miss, double *mean, double *std, bool * withMissings )
  {
    if(this->isActive)
    {
      if( this->state != streamReadyForProbs )
      {
        misc.error("Error: An internal error has happened", 0);
      }
      
      ProbabilitiesSetter setter( this->layout, probs, miss, mean, std, withMissings, this->individualRawIds.size(), &this->individualBGENIdxsMaskToKeep);
      genfile::bgen::read_and_parse_genotype_data_block< ProbabilitiesSetter >(
        *this->stream,
        this->context,
        setter,
        &this->tempBuffer1,
        &this->tempBuffer2
      );
      this->state = streamReadyForVariant;
    }
  }
  
  void ignoreProbabilities()
  {
    if(this->isActive)
    {
      if( this->state != streamReadyForProbs )
      {
        misc.error("Error: An internal error has happened", 0);
      }
      
      genfile::bgen::ignore_genotype_data_block( *this->stream, this->context );
      this->state = streamReadyForVariant;
    }
  }
};

std::vector<std::string> Genotype::getBGENFileSNPs(std::string f)
{
  std::vector<std::string> fileSNPs;
  
  misc.setGetElapsedTime("GenotypeGetSNPs");
  misc.message << "Reading SNPs info from file [ " << f << " ] ..." << std::endl;
  
  std::set<std::string> dummySet;
  ReadBGENData bgenParser( f, dummySet, true, true);
  int tempNSNPs = bgenParser.getNSNPs();
  
  for(int i=0; i<tempNSNPs; i++)
  {
    std::string chromosome; uint32_t position; std::string rsid; std::vector< std::string > alleles(2, "");
    bgenParser.readVariantInfo( &chromosome, &position, &rsid, &alleles );
    
    fileSNPs.push_back(rsid);
    
    bgenParser.ignoreProbabilities();
  }
  
  misc.message << "SNPs info loaded after " << misc.setGetElapsedTime("GenotypeGetSNPs", true) << std::endl;
  
  return fileSNPs;
}

void Genotype::readBGENFile(std::string f, std::set<std::string> & keepSNPs, std::set<std::string> & keepIndividualIds)
{
  misc.setGetElapsedTime("GenotypeLoad");
  misc.message << "Reading genotype data from file [ " << f << " ] ..." << std::endl;
  
  this->normalized = true;
  
  //Load the BGEN header, individual ids, and number of SNPs.
  ReadBGENData bgenParser( f, keepIndividualIds, communicator->myCol == 0 ); //bgenParser can only be accessed in the processes satisfying communicator->myCol == 0.
  if(communicator->myCol == 0)
  {
    if( communicator->mpiRoot == true )
    {
      bgenParser.showFileInfo();
      
      std::vector< std::pair<std::string, std::string> > rawIds = bgenParser.getIndividualRawIds();
      for(int i = 0; i < rawIds.size(); i++)
      {
        Individual individual;
        
        individual.familyID = rawIds[i].first;
        individual.individualID = rawIds[i].second;
        individual.paternalID = "NA";
        individual.maternalID = "NA";
        individual.sex = "NA";
        individual.phenotype = -1;
        
        std::string key = individual.familyID + "@" + individual.individualID;
        
        this->originalFileBedIndividualsIdxs[f][key] = i;
        
        this->individuals.push_back(individual);
        this->individualIds.push_back(key);
        if(this->individualIdsIdx.count(key) != 0)
        {
          misc.error("Error: The individual with id: " + individual.individualID + " appears more than one time in the genotypes BGEN file.", 0);
        }
        this->individualIdsIdx[key] = i;
      }
      this->nIndividuals = this->individuals.size();
      if(keepSNPs.size() == 0)
      {
        this->nSNPs = bgenParser.getNSNPs();
      }
      else
      {
        this->nSNPs = keepSNPs.size();
      }
    }
  }
  communicator->broadcast(&this->nIndividuals);
  communicator->broadcast(&this->nSNPs);
  
  this->genotypes->initParameters(this->nSNPs, this->nIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockRows);
  this->missings->initParameters(this->nSNPs, this->nIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockRows);
  
  double * probabilities; //Temporal array where a file row will be stored
  int * missingsArray; //Temporal array where a file row missings will be stored
  double * genotypeMatrixBlock; //Temporal matrix where a block of genotypes will be stored
  double * missingMatrixBlock; //Temporal matrix where a block of missing information will be stored
  if(communicator->myCol == 0)
  {
    probabilities = new double [this->nIndividuals];
    missingsArray = new int [this->nIndividuals];
    
    genotypeMatrixBlock = new double [this->nIndividuals*this->genotypes->nBlockRows];
    missingMatrixBlock = new double [this->nIndividuals*this->genotypes->nBlockRows];
  }
  
  //Read genotype data
  std::vector<int> loadedSNPsByThisProcess;
  loadedSNPsByThisProcess.clear();
  
  this->SNPs = std::vector<SNP>(this->nSNPs, SNP());
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
      
      int j = 0;
      while(j < nFileRows)
      {
        std::string chromosome; uint32_t position; std::string rsid; std::vector< std::string > alleles(2, "");
        bgenParser.readVariantInfo( &chromosome, &position, &rsid, &alleles );
        if(keepSNPs.size() != 0 && keepSNPs.find(rsid) == keepSNPs.end())
        {
          bgenParser.ignoreProbabilities();
          continue;
        }
        if(communicator->mpiRoot == true)
        {
          SNP snp = SNP(chromosome, rsid, i2s(position), alleles.at(0), alleles.at(1));
          this->originalFileBEDSNPIdxs[f][snp.name] = i+j;
          this->SNPs[i+j] = snp;
          this->SNPIds.push_back(snp.name);
          if(this->SNPIdsIdx.count(snp.name) != 0)
          {
            misc.error("Error: The SNP with name: " + snp.name + " appears more than one time in the genotypes BGEN file.", 0);
          }
          this->SNPIdsIdx[snp.name] = i+j;
        }
        
        if(communicator->myRow == sendProcessRow)
        {
          double mean;
          double std;
          bool withMissings;
          
          bgenParser.readProbabilities( probabilities, missingsArray, &mean, &std, &withMissings );
          parseBGENSNP(probabilities, missingsArray, &(genotypeMatrixBlock[j]), &(missingMatrixBlock[j]), &(SNPs[i+j]), mean, std); //Change from row major to col major.
          loadedSNPsByThisProcess.push_back(i+j);
        }
        else
        {
          bgenParser.ignoreProbabilities();
        }
        
        j++;
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
    delete [] probabilities;
    delete [] missingsArray;
    delete [] genotypeMatrixBlock;
    delete [] missingMatrixBlock;
    
  }
  
  gatherSNPData(loadedSNPsByThisProcess);

  Cblacs_barrier(communicator->context, "All");
  
  removeSNPsWithZeroVariance();
  
  misc.message << "Genotype loaded after " << misc.setGetElapsedTime("GenotypeLoad", true) << std::endl;
}

void Genotype::parseBGENSNP(double* probabilities, int* missingsArray, double *genotypes, double *missings, SNP* snp, double mean, double std)
{
  unsigned int keepIdxs = 0;
  for(int j=0; j<this->nIndividuals; j++)
  {
    genotypes[j*this->genotypes->nBlockRows] = probabilities[j];
    missings[j*this->genotypes->nBlockRows] = missingsArray[j];
  }
  
  snp->nNonMissing = -1;
  snp->p1 = 1. - (mean/2.);
  snp->p2 = mean/2.;
  snp->standardDev = std;
}


void Genotype::removeSNPsWithZeroVariance()
{
  std::vector<std::string> keepSNPs;
  std::vector<std::string> removeSNPs;
  if( communicator->mpiRoot == true )
  {
    for( int isnp = 0; isnp < this->SNPs.size(); isnp++ )
    {
      if(this->SNPs[isnp].standardDev != 0.)
      {
        keepSNPs.push_back(this->SNPIds[isnp]);
      }
      else
      {
        removeSNPs.push_back(this->SNPIds[isnp]);
      }
    }
  }
  
  if( misc.gt(removeSNPs.size() != 0) )
  {
    if( communicator->mpiRoot == true )
    {
      Message message(options.outFile + ".snps.novariance");
      for(int i = 0; i < removeSNPs.size(); i++)
      {
        message << removeSNPs[i] << std::endl;
      }
    }
    
    if(options.keepSNPsWithZeroVariance == false)
    {
      misc.message << "Warning: There are " << removeSNPs.size() << " SNPs which have zero variance and will be removed. They are stored in file [ " << (options.outFile + ".snps.novariance") << " ]." << std::endl;
      filterSNPsAndIndividuals(keepSNPs, this->individualIds);
    }
    else
    {
      misc.message << "Warning: There are " << removeSNPs.size() << " SNPs which have zero variance BUT WILL BE KEPT. They are stored in file [ " << (options.outFile + ".snps.novariance") << " ]." << std::endl;
    }
  }
}
#endif
