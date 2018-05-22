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

#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include "matrix.h"


namespace GenotypeAttributes
{
  enum Type
  {
    empty,
    calls,
    probabilities
  };
};

enum GroupBy
{
  ungrouped,            ///<SNPs not grouped
  byPosition,           ///<SNPs grouped by regions defined by their position.
  byGene,               ///<SNPs grouped by their gene defined in options.regionsFile
  byGroup,              ///<SNPs grouped by their group defined in options.regionsFile
  byOrderedFixedSize,   ///<SNPs grouped by ordered groups of fixed size. (i.e. order SNPs and then divide in groups of fixed size). Two SNPs in different chromosomes cannot be in the same group.
  byAll,                ///<Only one group with all SNPs
  byFileOrderedWindows  ///<Grouped by fixed size windows with SNPs and groups in the order of the original genotypes file.
};

/**
 * Store SNP information
 */
class SNP
{
public:
  std::string chr;	///< SNP Chromosome
  std::string name;	///< SNP name
  int position;		///< position in the chromosome
  std::string allele1;	///< Allele 1 description
  std::string allele2;	///< Allele 2 description
  
  int nNonMissing;
  int frequencies[4]; 	///< Number of individuals as a function of their genotype [0] num. missings [1] num. homozygote 1 [2] num. heterozygote [3] num. homozygote 2
  double p1; 		///< allele1 frequency
  double p2; 		///< allele2 frequency
  double standardDev;	///< sqrt(2*p1*(1-p1))
  
  SNP();
  SNP(std::string chr, std::string name, std::string position, std::string allele1, std::string allele2);
};

inline bool operator==(const SNP& lhs, const SNP& rhs){
  return (
    (lhs.chr == rhs.chr) && (lhs.name == rhs.name)  && (lhs.position == rhs.position) && (lhs.allele1 == rhs.allele1) && (lhs.allele2 == rhs.allele2) &&
    (lhs.nNonMissing == rhs.nNonMissing) && (lhs.frequencies[0] == rhs.frequencies[0]) && (lhs.frequencies[1] == rhs.frequencies[1]) && (lhs.frequencies[2] == rhs.frequencies[2]) &&
    (lhs.frequencies[3] == rhs.frequencies[3]) && (lhs.p1 == rhs.p1) && (lhs.p2 == rhs.p2) && (lhs.standardDev == rhs.standardDev) 
  );
}
inline bool operator!=(const SNP& lhs, const SNP& rhs){ return !operator==(lhs,rhs); }

/**
 * Store individual information
 */
class Individual
{
public:
  std::string familyID;
  std::string individualID;
  std::string paternalID;
  std::string maternalID;
  std::string sex;
  int phenotype;		///< The individual phenotype
};
inline bool operator==(const Individual& lhs, const Individual& rhs){
  return (
    (lhs.familyID == rhs.familyID) && (lhs.individualID == rhs.individualID)  && (lhs.paternalID == rhs.paternalID) && (lhs.maternalID == rhs.maternalID) && (lhs.sex == rhs.sex) &&
    (lhs.phenotype == rhs.phenotype)
  );
}
inline bool operator!=(const Individual& lhs, const Individual& rhs){ return !operator==(lhs,rhs); }

/**
 * Class for loading genotype data
 * This class contains the genotype and missing data. It also contains summary statistics for each SNP
 */
class Genotype
{
public:
  GenotypeAttributes::Type type;                ///<It is from calls or probabilities.
  
  Matrix *genotypes;                            ///<The matrix of genotypes. Values could be 0: missing, 1: homozygote 1, 2: heterozygote, 3: homozygote 2, before normalizing. SNPs in rows Individuals in columns
  Matrix *missings;                             ///<The matrix of missing genotypes. Values could be 0: missing 1: non-missing. SNPs in rows Individuals in columns.
  
  bool normalized;                              ///<The genotypes are normalized? true = yes, false = no.
  
  int nSNPs;                                    ///<number of SNPs
  std::vector<SNP> SNPs;                        ///<list of SNP classes. Element i has the information of the SNP in row i of genotypes Matrix
  std::vector<std::string> SNPIds;              ///<List of the SNP names: snp.name. Only in the mpiRoot process.
  std::map<std::string, int> SNPIdsIdx;         ///<Map between the SNP name and their position in the genotypes matrix. Only in the mpiRoot process.
  
  int nIndividuals;                             ///<number of individuals
  std::vector<Individual> individuals;          ///<List of the Individuals Ids. Only in the mpiRoot process.
  std::vector<std::string> individualIds;       ///<List of the individual keys: familyId + @ + individualId. Only in the mpiRoot process.
  std::map<std::string, int> individualIdsIdx;  ///<Map between the individual keys: familyId + @ + individualId and their position in the genotypes matrix. Only in the mpiRoot process.
  
  std::map<std::string, std::set<std::string> > groupedSNPs;    ///<The SNPs grouped. only in mpiRoot;
  GroupBy groupedBy;                                            ///<The grouping strategy used for grouping the SNPs;
  
  
  std::map< std::string, std::map<std::string, int> > originalFileBEDSNPIdxs;                ///<The dict of the original position of SNPs in the genotypes file before any kind of filtering.
  std::map< std::string, std::map<std::string, int> > originalFileBedIndividualsIdxs;        ///<The dict of the original position of Individuals in the genotypes file before any kind of filtering.
  
  std::string spaceChangedKernelName;           ///<The genotype matrix has been rotated using the kernel eigenvetors with name stored here.
  
  /**
   * Creates the Genotype and allocates memory for the genotypes and missings matrices
   */
  Genotype(DistributionType dist, int ngr, int ngc, int nbr, int nbc);
  
  /**
   * Creates the Genotype without allocating memory for the matrices
   */
  Genotype();
  
  /**
   * Creates the Genotype and loads data from file fname
   */
  Genotype(std::string fname, bool loadGenotypes = true);
  
  /**
   * Creates the Genotype and loads data from file fbase keeping only keepSNPs and keepIndividualIds
   * 
   * The pourpose of this function is avoid loading full genotype files which may require huge ammounts of memory.
   */
  Genotype( std::string fbase, std::set<std::string> & keepSNPs, std::set<std::string> & keepIndividualIds );
  
  /**
   * Creates the Genotype and loads data from file fbase, using SNPs and Inds from baseGenotypes and keeping only keepSNPs and keepIndividualIds
   * 
   * The pourpose of this function is avoid loading full genotype files which may require huge ammounts of memory.
   */
  Genotype( std::string fbase, Genotype * baseGenotype, std::set<std::string> & keepSNPs, std::set<std::string> & keepIndividualIds );
  
  
  /**
   * Creates the Genotype and loads data from file fbase keeping only keepSNPs and keepIndividualIds
   * 
   * The pourpose of this function is avoid loading full genotype files which may require huge ammounts of memory.
   * This constructor loads from another non-plink formats of genotype data (currently keepSNPs is ignored.).
   */
  Genotype( std::string fbase, std::set<std::string> & keepSNPs, std::set<std::string> & keepIndividualIds,  GenotypeAttributes::Type formatType);
  
  /**
   * Destructor
   */
  ~Genotype();
  
  
  /**
   * Load genotype data from file
   * \param f string with the base name of the files (without the .bed, .bim or .fam extensions)
   */
  void load(std::string f, bool loadGenotypes = true);
  
  /**
   * Load genotype data from file
   * \param f string with the base name of the files (without the .bed, .bim or .fam extensions)
   * \param baseGenotype baseGenotype with the Individual and SNP info preloaded.
   * \param keepSNPs Ids of the SNPs to load
   * \param keepIndividualIds Ids of the individuals to load
   */
  void load(std::string fbase, Genotype * baseGenotype, std::set<std::string> & keepSNPs, std::set<std::string> & keepIndividualIds);
  
  /**
   * Filter genotype as a function of option parameters.
   * 
   * Filter genotype as a function of option parameters. Currently this option only filters SNPs if a file with a list of SNPs to keep is specified.
   */
  void filterGenotypeUsingOptions();
  
  /**
   * Load genotype data from files specified in a list
   * \param f string with the base name of the file where the list is stored
   */
  void loadList(std::string f);
  
  /**
   * Read the SNPs info (BIM) file
   * \param f string with the name of the BIM file
   * \param keepSNPs SNPs to keep from the file.
   * \param SNPsBEDIdxsToKeep SNPs bed idxs to load.
   */
  void readBIMFile(std::string f, std::set<std::string> & keepSNPs, std::vector<int> & SNPsBEDIdxsToKeep);
  
  /**
   * Read the individuals info (FAM) file
   * \param f string with the name of the FAM file
   * \param keepIndividualIds Individual ids to keep from the file.
   * \param individualBEDIdxsMaskToKeep Mask of individual ids to keep from bed.
   */
  void readFAMFile(std::string f, std::set<std::string> & keepIndividualIds, std::vector<int> & individualBEDIdxsMaskToKeep);
  
  /**
   * Read the genotype data and store in matrix genotypes and missings
   * \param f string with the base name of the BED file
   * \param SNPsIdxsToKeep BED Index of the SNPs that have to be loaded from the file. If empty, all individuals will be loaded.
   * \param individualBEDIdxsMaskToKeep Mask of the BED individuals indexs that have to be loaded from the file (0, individual is going to be filtered, !=0 it is going to be kept). This has to be especified. If empty, an error is raised.
   */
  void readBEDFile(std::string f, std::vector<int> & SNPsBEDIdxsToKeep, std::vector<int> & individualBEDIdxsMaskToKeep);
  void readBEDFileOld(std::string f, std::vector<int> & SNPsBEDIdxsToKeep, std::vector<int> & individualBEDIdxsMaskToKeep);
  
  /**
   * Read a block of genotype data from the genotype file
   * 
   * This is a helping function for readBEDFile(). Reads a block of genotypes of size nIndividuals*BlockRows.
   * Where BlockRows is the row block size of the genotypes Matrix. The genotypes and missings are stored in
   * temporal arrays pointed by genotypes and missings function parameters.
   * 
   * \param nrf the number of bytes to read in a row. Each byte stores 4 genotypes. nrf = floor(nIndividuals/4)
   * \param rest the number of genotypes at the end of the row that not fill an entire byte. rest = nIndividuals%4
   * \param memBlock a block of memory for storing a row of the genotype file
   * \param genotypes a pointer to a temporal genotype matrix of size nIndividuals*BlockRows. BlockRows is the row block size of the genotypes Matrix
   * \param missings a pointer to a temporal missings matrix of size nIndividuals*BlockRows. BlockRows is the row block size of the missings Matrix
   * \param snp pointer to the current SNP class for updating allele frequencies
   * 
   * \sa readBEDFile()
   */
  void parseSNP(int nrf, int rest, unsigned char * memBlock, double *genotypes, double *missings, SNP* snp, const std::vector<int> & individualBEDIdxsMaskToKeep);
  
  /**
   * Parse a byte with the information of four genotypes
   * 
   * This is a helping function for parseSNP()
   * 
   * \param g a byte which can store up to four genotypes that will be parsed by the function.
   * \param n The number of genotypes stored in the byte.
   * \param d a parameter indicating if a logical shift must be performed on g
   * \param genotypes a pointer to a temporal genotype matrix of size nIndividuals*BlockRows where genotypes of g will be stored. BlockRows is the row block size of the genotypes Matrix
   * \param missings a pointer to a temporal missings matrix of size nIndividuals*BlockRows where genotypes of g will be stored. BlockRows is the row block size of the missings Matrix
   * \param snp pointer to the current SNP class for updating allele frequencies
   * 
   * \sa readBEDFile(), parseSNPbyte()
   */
  void parseSNPbyte(unsigned char g, int n, int d, double * genotypes, double * missings, SNP* snp, unsigned int & keepIdxs, const std::vector<int> & individualBEDIdxsMaskToKeep, int maskShift);
  
  
  /**
   * Gather SNP data in the root process
   * 
   * This is a helping function for readBEDFile. A portion of genotype data are readed in different threads. This function gathers this information on the root node.
   * 
   * \param loadedSNPsByThisProcess A vector with the indices of the SNPs loaded by the current thread.
   */
  void gatherSNPData(std::vector<int> & loadedSNPsByThisProcess);
  
  /**
   * Normalize genotype data using the GCTA method
   * 
   * Normalize genotype data using the GCTA method. The flag this->normalized will be set to true;
   */
  void normalizeGenotypes();
  
  /**
   * Filter SNPs and individuals from genotypes matrix
   * 
   * Take into account that if individuals are filtered, the current version of the function does not actualize SNP info regarding frequencies, mean values, etc.
   * 
   * \param unorderedKeepSNPIds List of SNPs that will be kept. Their value could be modified by this function (not anymore).
   * \param unorderedKeepIndividualIds List of individuals that will be kept. Their value could be modified by this function (not anymore).
   * \param filterMissings if true, this->missings matrix will be filtered too. Otherwise, it will be deleted.
   * \param newGenotype if NULL, genotypes will be filtered in place. Otherwise, it will be stored in the Genotype() pointed by newGenotype.
   */
  void filterSNPsAndIndividuals(std::vector<std::string> & unorderedKeepSNPIds, std::vector<std::string> & unorderedKeepIndividualIds, bool filterMissings = true, Genotype *newGenotype = NULL);
  
  /**
   * Append a new genotype to the current genotype.
   * 
   * The new genotype must met one of these two conditions: new genotype have same SNPs but different individuals or new genotype have different SNPs but same individuals.
   * 
   * \param newGenotype A pointer to the genotype to append.
   */
  void appendGenotype(Genotype * newGenotype);
  
  /**
   * Group SNPs using different criteria.
   * 
   * The SNPs are grouped and the groups stored in this->groupedSNPs. The grouping method is stored in this->groupedBy.
   * 
   * \param groupBy The criteria used for grouping the SNPs.
   * \param SNPIdsSubsetArg For the moment, This option only applies when grouping byOrderedFixedSize. If specified, groups would be made by using only SNPs in this list of SNPIds. 
   */
  void groupSNPs(GroupBy groupBy, std::vector<std::string> SNPIdsSubset = std::vector<std::string>());
  
  /**
   * Auxiliar method used by groupSNPs method.
   */
  void groupSNPsByPosition();
  
  /**
   * Auxiliar method used by groupSNPs method.
   */
  void groupSNPsByGene();
  
  /**
   * Auxiliar method used by groupSNPs method.
   * 
   * It groups SNPs in an ordered manner assuming
   * 
   * \param groupsSize The size of the groups.
   * \param SNPIdsSubsetArg If specified, groups would be made by using only SNPs in this list of SNPIds.
   */
  void groupSNPsByOrderedFixedSize(int groupsSize, std::vector<std::string> SNPIdsSubsetArg = std::vector<std::string>());
  
  /**
   * Grouped by fixed size windows with SNPs and groups in the order of the original genotypes file.
   * 
   * \param windowSize The size of the window.
   */
  void groupSNPsByFileOrderedWindows(int windowSize);
  
  /**
   * Get a vector of SNPs in a particular group
   * 
   * \param group Name of group of SNPs. Only needs to be defined in root process.
   * \param writeGroup Store the group elements in a file?
   * \return Unordered list of SNPs in group. Only in root process.
   */
  std::vector<std::string> getGroup(std::string group, bool writeGroup = true);
  
  /**
   * Get the genotype of SNPs in a particular group
   * 
   * \param group Name of group of SNPs. Only needs to be defined in root process.
   * \param newGenotype Pointer to Genotype where tha genotype data for the selected group will be stored.
   * \param writeGroup Store the group elements in a file?
   */
  void genotypeOfSNPsGroup(std::string group, Genotype *newGenotype, bool writeGroup = true);
  
  /**
   * This function splits and redistributes the genotypes matrix to grouped genotype matrices defined with a new communicator.
   * 
   * \param newCommunicator The new grouped communicator which will be used for redistributing the matrix. Its type has to be basicGroupedCommunicator.
   * \param SNPidxs Return the index of SNPs in each communicator group.
   * \return Pointer to the redistributed genotype matrices.
   */
  Matrix* genotypesRedistributionToGroupedCommunicatorMatrices(Communicator * newCommunicator, std::map< int, std::vector<int> > & SNPidxs);
  
  
  /**
   * beta functions for loading from BGEN file
   */
  std::vector<std::string> getBGENFileSNPs(std::string f);
  void parseBGENSNP(double* probabilities, int* missingsArray, double *genotypes, double *missings, SNP* snp, double mean, double std);
  void readBGENFile(std::string f, std::set<std::string> & keepSNPs, std::set<std::string> & keepIndividualIds);
  
  /**
   * Remove SNPs with zero variance.
   */
  void removeSNPsWithZeroVariance();
  
  /**
   * Change the space of the genotypes matrix using the eigenvectors from a kernel.
   * 
   * ATTENTION: Currtenly the function just sets the name of the kernel used to chenge the space. It has to be done manually before calling the function. This has to be updated.
   */
  void changeSpace(std::string name);
  
  /**
   * Function for debugging pourposes thhat prints the current genotypes and missing matrices
   */
  void printGenotype(bool showGroups = false, bool unnormalizeProbabilities = false);
  /**
   * Function for debugging pourposes thhat prints the current genotypes and missing matrices
   */
  void printGenotypeT(bool unnormalizeProbabilities = false);
};

#endif
