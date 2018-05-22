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

#ifndef REML_H
#define REML_H

#include "matrix.h"
#include "genotype.h"
#include "covariate.h"
#include "covariancematrix.h"

#include <vector>
#include <set>
#include <utility>
#include <fstream>

class Matrix;
class GRM;
class Genotype;
class REML;
class CovarianceMatrix;
class Message;

namespace REMLParameters
{
  enum SecondDerivativesMatrixTypes
  {
    AI,
    MLF,
    REMLF,
    undefined
  };
};

enum REMLType
{
  rawREMLType,
  singleREMLType,
  multipleSingleREMLType,
  regionalSingleREMLType,
  bivariateREMLType,
  multipleBivariateREMLType,
  regionalBivariateREMLType
};

struct REMLResults
{
  std::vector<Variance> variances;
  double logLikelihood;
};

struct CovariateNames
{
  std::vector< std::pair<std::string, int> > meanNames;                           ///<Names of the means and their column position on X matrix. Only on root process.
  std::vector< std::pair<std::string, int> > quantitativeCovarNames;              ///<Names of the quantitative covariates and their column position on X matrix. Only on root process.
  std::vector< std::pair<std::string, int> > discreteCovarNames;                  ///<Names of the discrete covariates and their column position on X matrix. Only on root process.
};

class REML
{
public:
  REMLType type;
  
  bool useMLinsteadOfREML;
  
  Matrix * y;					///<Phenotypes matrix.
  std::vector<Matrix *> yList;                  ///<Vector of phenotypes matrices (for averaging likelihoods with equal incidence and covariances over different uncorrelated samples). ATTENTION: Currently it is not used on EM steps (it has to be implemented).
  Matrix * X;					///<Covariates matrix.
  CovarianceMatrix * V;				///<The total covariance matrix.
  
  double logDetXtViX;
  double logDetV;
  double logLikelihood;
  double logLikelihoodDifference;
  double logLikelihoodRelativeDifference;
  
  bool singlePrecisionInversion;                ///<The invert of the covariance matrix will be performed using single precision?
  
  int currentStep;

  Matrix * ViX;  
  Matrix * XtViX_i;
  Matrix * P;
  Matrix * PDiagonal;
  Matrix * Py;
  Matrix * subVPy;
  Matrix * AI;                                  ///<Stores the AI, ML F or REML F inverses, depending on the analysis. It is defined by the variable secondDerivativesMatrixType
  Matrix * yPsubVPy_trPsubV;                    ///<This is only used when this->useMLinsteadOfREML == false
  Matrix * yPsubVPy_trVisubV;                   ///<This is only used when this->useMLinsteadOfREML == true
  Matrix * VisubVTraces;                        ///<This is only used when this->useMLinsteadOfREML == true
  double ytPy;                                  
  
  REMLParameters::SecondDerivativesMatrixTypes secondDerivativesMatrixType;   ///<Indicates what is stored in AI matrix.
  
  std::vector<double> oldVariances;             ///<The old variances before a REML step
  std::vector<double> delta;                    ///<The variance changes after a REML step: variances = oldVariances + delta
  
  int dimension;
  
  int nPhenotypes;
  
  bool usingDiagonalKernels;                                             ///<Computation is based on diagonal kernels?
  
  std::vector<CovariateNames> covariateNames;                            ///<Vector of covariate names for each phenotype (e.g. for bivariate, vector with two components.)

  std::vector< std::vector<Individual> > vIndividuals;                   ///<individual Ids for the different phenotypes analyzed. They are in the same order than in the covariance matrix.
  std::vector<std::string> individualBLUPNames;                          ///<List with the covarianceIds of the individual BLUPs to compute
  std::map<std::string, std::vector<std::string> > mSNPIds;              ///<Map of list of SNPs used for computing GRMs. Only defined for singleREMLType and regionalSingleREMLType analysis.  The key must be the same that the GRM name used for preparing the covariance matrix.
  std::map<std::string, Genotype *> SNPsBLUPGenotypes;                   ///<Map of pointer to genotypes which could be used for computing SNP BLUPs. The key must be the same that the GRM name used for preparing the covariance matrix.
  std::map<std::string, std::vector<std::string> > SNPsBLUPGenotypeFiles;///<Map of pointer to genotype files which could be used for computing SNP BLUPs. The key must be the same that the GRM name used for preparing the covariance matrix.
  std::map<std::string, Matrix *> GRMEigenVectors;                       ///<Pointers to GRM EigenVectors in case REML is computed with diagonalized GRM. NULL otherwise.
  
  bool success;                                 ///<computeREML converged?
  bool logLikelihoodConverged;                  ///<log Likelihood converged?
  bool variancesConverged;                      ///<variances converged?
  bool gradientConverged;                       ///<gradient converged?
  
  bool writeResults;                            ///<write results?
  bool writeBLUE;                               ///<Can override writeResults. If writeResults is false, but this is true, BLUE results will be stored.
  bool writeSummary;                            ///<Can override writeResults. If writeResults is false, but this is true, Summary results will be stored.
  
  std::set<std::string> constrainedVariances;   ///<Stores the constrained variances in a REML step.
  
  std::set<std::string> warnings;               ///<Stores the warnings produced during REML.
  
  std::string stepModifications;
  
  std::vector<std::string> elementsToTest;                               ///<List of Element.subCovarianceId's to test when removed.
  std::vector< std::pair<std::string, double> > covarsToFix;             ///<List of covars to fix and the value to be fixed.
  
  std::string sREMLML;
  
  REML(bool argWriteResults = true);
  ~REML();
  
  /**
   * Delete intermediate matrices used in REML computation
   */
  void deleteIntermediateMatrices();
  
  /**
   * Initialize REML parameters assuming multiple GRMs
   * 
   * All parameters must be defined on all processes.
   * 
   * \param type A type defining the REML
   * \param kernels Pointers to the kernels that will be used for preparing REML. All these grms will be deleted and this vector cleared.
   * \param weights The initial variance weights for each kernel. I empty the weights will be balanced equally on all kernels.
   * \param phenotypeColumns A list of the columns of the phenotypes file ised for fitting the REML.
   * \param heritabilities A list of initial heritability assumed for each phenotype column. If it is empty, an heritability of 0.5 is assumed.
   * \param covariateFiles The covariate file that should be used for each phenotype.
   * \return true if REML computation can be started, false otherwise.
   */
  bool prepare(REMLType type, std::vector<Kernel*> & kernels, std::vector<double> weights, std::vector<int> phenotypeColumns, std::vector<double> heritabilities, std::vector<std::pair<std::string, std::string> > covariateFiles);
  
  /**
   * Initialize REML parameters from raw matrices for performing raw REML
   * 
   * \param yparam matrix with phenotypes
   * \param Xparam matrix with covariates
   * \param kernels vector of matrices representing the kernels to use.
   * \param heritability initial heritability to be used
   * \param weights The initial variance weights for each kernel. I empty the weights will be balanced equally on all kernels.
   * \return true if REML computation can be started, false otherwise.
   */
  bool prepare(Matrix* yparam, Matrix* Xparam, std::vector<Matrix*> & kernels, std::vector<KernelType> kernelTypes, double heritability, std::vector<double> weights);
  
  /**
   * Create a covar matrix from another matrix.
   * 
   * Suppose that srcCovarMatrix = M and transpose(srcCovarMatrix) = Mt. This function creates a new matrix R as:
   * 
   *       [ 0  M ]
   *       [ Mt 0 ]
   * 
   * Assumes that M is filled on upper and lower parts of the diagonal. (i.e. if M is symmetric and uplo != 'B', the resultant matrix could be wrong).
   * 
   * \param srcCovarMatrix The source matrix (M matrix of description).
   * \param resultCovarMatrix The resulting matrix (R matrix of description).
   */
  void createCovarMatrix(Matrix * srcCovarMatrix, Matrix * resultCovarMatrix);
  
  /**
   * Compute an initial estimation of the variance, fitting a OLS
   * 
   * if Xt*X is invertible returns var( X * ((Xt*X)-1) * Xt*y ). var( y ) otherwise.
   * 
   * \param yMatrix Matrix of phenotypes (y in the description).
   * \param XMatrix Design matrix (X in the tescription).
   * \return The computed variance
   */ 
  double computeInitialVariance(Matrix* yMatrix, Matrix* XMatrix);
  
  /**
   * Transform the matrix m using the kernel eigenvectors defined by name.
   * 
   * Transform the matrix m using the kernel eigenvectors defined by name. If the eigenvectors are A, we can define
   * 
   *        [ At 0  0  ] 
   *   AB = [ 0  At 0  ] 
   *        [ 0  0  At ]
   * 
   * where At is the transposed of A if t == 'T' (At==A if t=='N') and AB is a block matrix with number of blocks defined by numberBlocks (in the example, numberBlocks = 3).
   * The m is transformed doint the product:
   * 
   *  AB*m
   * 
   * is performed. numberBlocks usually depends on the number of phenotypes being analyzed. 
   * 
   * \param name The name of the diagonalized kernel, used for get their eigenvectors.
   * \param m a pointer to a Matrix pointer to be transformed.
   * \param numberBlocks the number of blocks used for transform m. It usually depends on the number of phenotypes analyzed. m number rows have to be equal to numberBlocks*(eigenvalues->nGlobCols)
   */
  void transformUsingEigenVectors(std::string name, char t, Matrix ** m, int numberBlocks);
  
  /**
   * Fill this meanNames, discreteCovarNames, and quantitativeCovarName.
   * 
   * Fills this meanNames, discreteCovarNames, and quantitativeCovarName respectively with names on cm. The index is an increasing index starting with shift value.
   * A prefix string can be added on the names.
   * 
   * \param cm covariates matrix from where getting the names.
   * \param prefix add this prefix on the names
   * \param shift Number of columns that contain means in the covariate matrices.
   */
  void addCovariatesNames(Covariate * cm, std::string prefix = "", int shift = 0);
  
  /**
   * If this->nPhenotypes > 1 returns sufix. Returns "" otherwise.
   */
  std::string addSuffix(std::string suffix);
  
  void setyList(std::vector<Matrix *> & ys);
  void deleteyList();
  
  /**
   * Performs REML computation and test nested models
   * 
   * Before starting REML computation, this->V is reinitialized. This function can be called consecutively for computing REML again.
   */
  void computeREMLWithReducedModels();
  
  /**
   * Performs REML computation
   * 
   * Before starting REML computation, this->V is reinitialized. This function can be called consecutively for computing REML again.
   * 
   * \return this->logLikelihood
   */
  double computeREML();
  
  void aiREMLStep();
  void emREMLStep();
  void emPartialREMLStep(Matrix * newVariances, Matrix * delta, double scale);
  
  void computePMatrix();
  //void computePyMatrix();
  void computePxMatrix(Matrix * Px, Matrix * x);
  void computeytPy();
  void computeSubVPyMatrix();
  void computeAIMatrix();
  void computeAIMatrixCrossedDerivatesCorrection();
  void computeMLFMatrix();
  void computeREMLFMatrix();
  void computeyPsubVPy_trPsubVVector(double scale = 1.);
  void computeyPsubVPy_trVisubVVector(double scale = 1.);
  
  double computeLogLikelihood();
  
  double varianceDifferenceNorm();
  bool allVariancesRelativeDifferencesLowerThan(double threshold);
  bool gradientLowerThanThreshold(double threshold);
  void showNonConvergenceCauses();
  void fixVariancesAndParameters();
  void unfixVariancesAndParameters();
  int getNVariancesFixed();
  
  void computeSummary();
  void computeBLUE();
  void computeIndividualsBLUP();
  std::vector<double> computeBLUPErrors(std::string name);
  void computeSNPsBLUP();
  void computeSNPsBLUP(std::string name, Genotype* genotypes);
  
  void verboseLogOutput1();
  
  void experimentalEpistasisPrediction();
  
  double logistic(double x);
  double logisticInv(double x);
};

#endif
