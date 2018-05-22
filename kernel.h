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

#ifndef KERNEL_H
#define KERNEL_H

#include "matrix.h"
#include "genotype.h"
#include "covariate.h"

#include <string>
#include <vector>
#include <map>

class Options;

enum KernelType
{
  kernelGRM,
  kernelEpistaticGRM,
  kernelFromDiscreteCovariates,
  kernelFromMultiDiscreteCovariates,
  kernelFromContinuousCovariates,
  kernelSquaredExponential,
  kernelCovarianceMatrix,
  kernelEnvirontmental,
  kernelInteraction,
  kernelGCTAGRM
};

class Kernel
{
public:
  KernelType type;                              ///<The type of the Kernel. i.e. GRM
  
  std::string name;                             ///<The name of the kernel
  
  Matrix *kernel;			        ///<The Kernel matrix. e.g. a grm. Set to NULL when diagonalized == true.
  Matrix *N;					///<An auxiliar normalization matrix which could be used for computing the Kernel. For instance, in grm cases, a matrix that stores the number of genotypes shared between two individuals.  Set to NULL when diagonalized == true.
  
  Matrix *eigenValues;                          ///<The eigenvalues of the diagonalized kernel. Set to NULL when diagonalized == false.
  Matrix *eigenVectors;                         ///<The eigenvectors after kenel diagonalization. Set to NULL when diagonalized == false.
  
  std::vector<Individual> individuals;  	///<List of the Individuals Ids. Only in the mpiRoot process. Only valid if asymmetric == false.
  std::vector<std::string> individualIds;	///<List of the individual keys: familyId + @ + individualId. Only in the mpiRoot process. Only valid if asymmetric == false.
  std::map<std::string, int> individualIdsIdx;	///<Map between the individual keys: familyId + @ + individualId and their position in the grm. Only in the mpiRoot process. Only valid if asymmetric == false.
  int nIndividuals;				///<Total number of individuals. Only valid if asymmetric == false.
  
  std::vector<std::string> randomVarNames;      ///<Names of the variables used for computing the Kernel. For instance, for a GRM it will be the SNP ids used for computing the grm.
  
  bool normalized;                                      ///<true if GRM is normalized (i.e. each grm element is divided by the corresponding N element).
  
  bool asymmetric;                                      ///<true if individuals in rows and columns are different. false otherwise.
  
  bool diagonalized;                                    ///<If grm diagonalized, grm is stored as their diagonal (in eigenValues variable) an their eigenvectors (in eigenVectors).  grm = NULL and N = NULL.
  
  std::vector<Individual> individualsRows;              ///<List of the Individuals Ids in rows. Only in the mpiRoot process. Only valid if asymmetric == true.
  std::vector<std::string> individualIdsRows;           ///<List of the individual keys in rows: familyId + @ + individualId. Only in the mpiRoot process. Only valid if asymmetric == true.
  std::map<std::string, int> individualIdsIdxRows;      ///<Map between the individual keys in rows: familyId + @ + individualId and their position in the grm rows. Only in the mpiRoot process. Only valid if asymmetric == true.
  int nIndividualsRows;                                 ///<Total number of individuals in rows. Only valid if asymmetric == true.
  
  std::vector<Individual> individualsCols;              ///<List of the Individuals Ids in columns. Only in the mpiRoot process. Only valid if asymmetric == true.
  std::vector<std::string> individualIdsCols;           ///<List of the individual keys in columns: familyId + @ + individualId. Only in the mpiRoot process. Only valid if asymmetric == true.
  std::map<std::string, int> individualIdsIdxCols;      ///<Map between the individual keys in columns: familyId + @ + individualId and their position in the grm columns. Only in the mpiRoot process. Only valid if asymmetric == true.
  int nIndividualsCols;                                 ///<Total number of individuals in columns. Only valid if asymmetric == true.
  
  Genotype * genotypes;
  Covariate * covariates;
  
  
  /**
   * Constructor that creates an empty kernel of type grm
   */
  Kernel();
  
  /**
   * Constructor that generates a grm from a Genotype class
   * 
   * \param genotype pointer to the genotype class
   */
  Kernel(Genotype * genotype, bool normalizeGRM = true);
  
  /**
   * Constructor that creates a Kernel of grm type and loads a grm from a file
   * 
   * \param f the base name of the file that contains the grm (without the extensions .grm.id or .grm.dat)
   * \param readType the type of the Kernel to load
   * \param useColumn If the type of the kernel is kernelFromDiscreteCovariates, then thie var indicates which column of the file will be used for creating the intersection matrix.
   * \param keepIndividuals If keepIndividuals.size() != 0, the kernel will be created using only the intersection between individuals in file and this set. Only required in root. (Currently, only used for kernelSquaredExponential, pending implementation for others).
   */
  Kernel(std::string f, KernelType readType = kernelGRM, int useColumn = 0, std::set<std::string> keepIndividuals = std::set<std::string>());
  
  /**
   * Constructor that creates a new Kernel from an existent Kernel
   * 
   * \param kernel The kernel that will be copied
   */
  Kernel(Kernel * kernel);
  
  /**
   * Constructor that generates a new kernel from an existent kernel
   * 
   * Constructor that generates a new kernel from an existent kernel. For instance, an epistatic kernel can be created from standard GRM.
   * 
   * \param kernel The kernel used as base kernel
   * \param createType Tht type for the new kernel
   */
  Kernel(Kernel * kernel, KernelType createType, Genotype* genotype);
  
  /**
   * Constructor that generates an interaction kernel between two previous kernels.
   * 
   * This functions create a new kernel which is the elementwise product of the components of both kernels. The individuals used is the overlapping individuals between k1 and k2.
   *
   * \param k1 The kernel used as base kernel
   * \param k2 The second kernel used for computing the elementwise product.
   */
  Kernel(Kernel * k1, Kernel * k2);
  
  /**
   * Destructor
   */
  ~Kernel();
  
  /**
   * Copy parameters from existent kernel
   * 
   * \param kernel Kernel from which parameters will be copied.
   */
  void copyParameters(Kernel *kernel);
  
  /**
   * Normalize GRM
   * 
   * Normalizes the GRM. Namely, element-wise division: grm/N. normalized flag is set to true.
   */
  void normalize();
  
  /**
   * Denormalize GRM
   * 
   * Denormalizes the GRM. Namely, element-wise multiplication: grm*N. normalized flag is set to false.
   */
  void denormalize();
  
  /**
   * Returns the normalized Kernel
   * 
   * Returns the normalized this->kernel if this->diagonalized == false or this->eigenValues if this->diagonalized == true
   */
  Matrix * getNormalizedKernel();
  
  /**
   * Create kernel from discrete covariates.
   * 
   * This function must only be called from the constructor.
   * 
   * \param fn The file name where the covariates are stored.
   * \param useColumn Use this column from the discrete covariates file for constructing the kernel.
   * \param keepIndividuals If non-empty set, then only keep the individuals present in the set.
   */
  void createKernelFromDiscreteCovariates(std::string fn, int useColumn, std::set<std::string> keepIndividuals);
  
  /**
   * Create kernel from multiple discrete covariates.
   * 
   * This function must only be called from the constructor.
   * 
   * \param fn The file name where the covariates are stored.
   * \param useColumn Use this column from the discrete covariates file for constructing the kernel.
   * \param keepIndividuals If non-empty set, then only keep the individuals present in the set.
   */
  void createKernelFromMultipleDiscreteCovariates(std::string fn, int useColumn, std::set<std::string> keepIndividuals);

  /**
   * Create kernel from continuous covariates.
   * 
   * This function must only be called from the constructor.
   * 
   * \param fn The file name where the covariates are stored.
   */
  void createKernelFromContinuousCovariates(std::string fn);
  
  /**
   * Create a squared exponential kernel.
   * 
   * This function must only be called from the constructor.
   * 
   * \param fn The file name where the values for creating the kernel are stored.
   * \param keepIndividuals If keepIndividuals.size() != 0, the kernel will be created using only the intersection between individuals in file and this set. Only required in root.
   */
  void createKernelSquaredExponential(std::string fn, std::set<std::string> keepIndividuals);
  
  /**
   * Compute the GRM from genotype data
   * 
   * \param genotype pointer to genotype data used for computation
   */
  void computeGRM(Genotype * genotype);
  
  /**
   * Save Kernel to a file
   * 
   * \param fn the base name of the file (without the extensions, e.g. .grm.id or .grm.dat)
   */
  void writeKernel(std::string fn);
  
  /**
   * Save GRM to a file
   * 
   * \param fn the base name of the file (without the extensions .grm.id or .grm.dat)
   */
  void writeGRM(std::string fn);
  
  /**
   * Load Kernel from a file
   * 
   * \param fn the base name of the file (without the extensions, e.g. .grm.id or .grm.dat)
   */
  void readKernel(std::string f);
  
  /**
   * Load GRM from a file
   * 
   * \param fn the base name of the file (without the extensions .grm.id or .grm.dat)
   */
  void readGRM(std::string f);
  
  /**
   * Load GRM from a file in GCTA gz format.
   * 
   * \param fn the base name of the file (without the extensions .grm.id or .grm.gz)
   */
  void readGCTAGRM(std::string f);
  
  /**
   * Filter individuals from GRM
   * 
   * \param keepIndividualIds List of individuals that will be kept. The order of this list will determine the order of the individuals in the new Kernel.
   * \param filterN if true, this->N matrix will be filtered too (if apply for kernel type). Otherwise, it will be deleted.
   */
  void filterIndividuals(std::vector<std::string> & keepIndividualIds, bool filterN = true);
  
  /**
   * Filter individuals from GRM. Different individuals in rows and columns.
   * 
   * Filter individuals uwing different individuals in rows and columns. this->asymmetric flag fill be set to true.
   * 
   * \param keepIndividualIdsRows List of individuals that will be kept in rows. The order of this list will determine the order of the individuals in the new Kernel.
   * \param keepIndividualIdsCols List of individuals that will be kept in columns. The order of this list will determine the order of the individuals in the new Kernel.
   * \param filterN if true, this->N matrix will be filtered too (if apply for kernel type). Otherwise, it will be deleted.
   */
  void filterIndividualsAsymmetric(std::vector<std::string> & keepIndividualIdsRows, std::vector<std::string> & keepIndividualIdsCols, bool filterN = true);
  
  /**
   * Replace the individual data for new individual data and indexes.
   * 
   * Only for symmetric kernels.
   * 
   * \param newIndividuals The new individual data. New individualIds indexes will be created from this. The vector elements has to be in the proper order for the current Kernel matrix.
   */
  void replaceIndividualIds(std::vector<Individual> & newIndividuals);
  
  /**
   * Add two Kernels
   * 
   *  
   * Performs the operation of adding two kernels: Performs the operation this = scalingFactor1*grm1 + scalingFactor2*grm2. How this is done, depends on the Kernel type.
   * For GRMs, see addGRMs() method description. kernel1 and kernel2 must be different than this. If kernel2 == NULL, then kernel1 will be added to this:
   * this = scalingFactor2*this + scalingFactor1*kernel1. randomVarNames behaviour strongly depends on Kernel type. Resultant Kernel name is the joined names of grm1
   * and grm2.
   * 
   * \param scalingFactor1 Scaling factor for kernel1.
   * \param grm1 the first GRM for the sum.
   * \param scalingFactor2 Scaling factor for kernel2.
   * \param grm2 the second GRM for the sum.
   */
  void addKernels(double scalingFactor1, Kernel *kernel1, double scalingFactor2 = 1., Kernel *kernel2 = NULL);
  
  /**
   * Add two GRMs
   * 
   * Performs the operation this = scalingFactor1*grm1 + scalingFactor2*grm2. Before the sum, the GRMs are denormalized. grm1 and grm2 must be different than this.
   * If grm2 == NULL, then grm1 will be added to this: this = scalingFactor2*this + scalingFactor1*grm1. If scalingFactors are equal, SNPIds from two GRMs are joined, if
   * scaling factors are different, SNPIds are substracted according to scalingFactors signs.
   * 
   * \param scalingFactor1 Scaling factor for grm1. Now it can only be 1. or -1. indicating the GRM is substracted or added.
   * \param grm1 the first GRM for the sum.
   * \param scalingFactor2 Scaling factor for grm2
   * \param grm2 the second GRM for the sum. Now it can only be 1. or -1. indicating the GRM is substracted or added.
   */
  void addGRMs(double scalingFactor1, Kernel *grm1, double scalingFactor2 = 1., Kernel *grm2 = NULL);
  
  /**
   * Returns a list of individuals to keep after removing the minimum number of individuals needed for eliminating all pairs of the input.
   * 
   * e.g. If globalIdxs1 = {2, 3, 5} and globalIdxs1 = {3, 2, 6}, i.e., pairs (3, 2), (2,3), (5,6). By eliminating only indexes 2 and 5,
   * we are eliminating all pairs.
   * Only runs on root node.
   * 
   * \param globalIdxs1 list of first individuals of each pair
   * \param globalIdxs2 list of second individuals of each pair
   */
  std::vector<std::string> minimumIndividualsRemovePairs(std::vector<int> globalIdxs1, std::vector<int> globalIdxs2);
  
  /**
   * Search the ids of the maximum number of individuals can be kept to eliminate all elements in the GRM greather than a particular cutoff.
   * 
   * \param cutoff the threshold over which the GRM elements must be filtered.
   */
  std::vector<std::string> searchNoHighRelatedIndividuals(double cutoff);
  
  /**
   * Prune a kernel removing the highest related individuals
   * 
   * Remove the minimum number of individuals such there are not relations with a value above cutoff
   * 
   * \param cutoff threshold above which relations will be pruned.
   */
  void pruneKernel(double cutoff);
  
  /**
   * Apply different filterings and controls to the Kernel.
   *
   * Apply different filterings and controls to the Kernel based on the options passed in the command line.
   * 
   * \return true if everythink is ok. false if there is a problem with the kernel or the kernel is diagonal.
   */
  bool sanitizeKernel();
  
  /**
   * Keep only all those individuals which have at least a relationship smaller than lowerThreshold or upper than upperThreshold.
   * 
   * \param lowerThreshold lower boundary
   * \param upperThreshold upper boundary
   */
  void keepWithRelatednessOutside(double lowerThreshold, double upperThreshold);
  
  /**
   * Randomly filter this grm
   * 
   * \param fraction fraction of individuals that will be kept.
   * \param limit minimum number of individuals that must be kept.
   */
  void randomSubSample(double fraction, int minimum);
  
  /**
   * Diagonalize this kernel and store their eigenvalues and their eigenvectors.
   */
  void diagonalizeKernel();
  
  /**
   * Recover original kernel from their eigenvalues and their eigenvectors.
   */
  void recoverKernelFromEigenDecomposition();
  
  /**
   * Function for debugging pourposes thhat prints the current GRM
   */
  void printGRM(bool showSNPs = false, int padding = 11);
};

#endif
