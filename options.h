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

#ifndef OPTIONS_H
#define OPTIONS_H

#include "range.h"
#include "genotype.h"

#include <string>
#include <sstream>
#include <map>

#define BASEGRM "GRM"

enum AnalysisToPerform
{
  withoutAnalysis,                         ///<No analysis specified
  justCheckAnalysis,                       ///<Just check if the options can be parsed properly
  showHelpAnalysis,                        ///<Show help?
  makeGRMAnalysis,                         ///<Compute the GRM?
  makeGRMKeepMostRelated,                  ///<Compute the GRM and keep only the most related.
  REMLAnalysis,                            ///<Perform REML analysis?
  GWASAnalysis,                            ///<Perform GWAS analysis?
  recursiveGWASAnalysis,                   ///<Perform a recurive GWAS analysis.
  iGWASAnalysis,                           ///<Perform an inverse GWAS analysis?
  bivarREMLAnalysis,                       ///<Perform bivariate REML analysis?
  multivarREMLAnalysis,                    ///<Perform multi REML analysis?
  simulatePhenotypeAnalysis,               ///<SimulatePhenotypes?
  predictPhenotypeAnalysis,                ///<PredictPhenotypes?
  PCAAnalysis,                             ///<Perform a PCA?
  accuracyBySNPAnalysis,                   ///<Compute the accuracy as a function of SNPs included in the predictor. Experimental analysis.
  effectsAnalysis,                         ///<Perform analysis from gene effects matrices.
  GLMMAnalysis,                            ///<Perform GLMM analysis?
  predictCovarPhenotypeAnalysis,           ///<Predict phenotypes but only from covars.
  multiplePhenotypeGWASAnalysis,           ///<Perform Multiple phenotype GWAS analysis.
  multiplePhenotypeResiduals,              ///<Compute residuals from multiple phenotypes for performing Multiple phenotype GWAS analysis (multiplePhenotypeGWASAnalysis).
  filterLabeledMatrix,                     ///<Filter LabeledMatrices.
  computeSNPStatsAnalysis                  ///<Filter LabeledMatrices.
};

struct REMLCovarianceRestricted
{
  std::string baseId;                      ///<The base Id for the variance.
  int p1;                                  ///<The phenotype 1.
  int p2;                                  ///<The phenotype 2.
  double correlation;                      ///<The fixed correlation.
};

class Options
{
public:
  ////////////////////////////
  // Actions
  AnalysisToPerform analysis;
  
  ////////////////////////////
  // Base input files
  std::string genotypeFile;             ///<Genotype file
  std::string genotypeListFile;         ///<File with a list of genotype files
  std::string bgenGenotypeFile;         ///<BGEN Genotype file.
  std::string grmFile;                  ///<GRM file
  std::string grmListFile;              ///<File with a list of GRMs
  std::string phenotypesFile;           ///<Phenotypes file
  std::vector<std::string> phenotypesFiles;     ///<List of phenotype files
  std::string covarsFile;               ///<Covariates file
  std::vector<std::string> covarsFiles; ///<Covariates files list (for bivariate multi REML analysis)
  std::string qCovarsFile;              ///<Quantitative covariates file
  std::vector<std::string> qCovarsFiles;///<Quantitative covariates files list (for bivariate multi REML analysis)
  std::string GCTAGRMsFile;             ///<File defining GRMs in GCTA format.
  
  ////////////////////////////
  // Output files
  std::string outFile;                  ///<Output file
  bool compressedOutput;                ///<Compress output?
  
  ////////////////////////////
  // GRM
  int GRMJoinMethod;                    ///<How GRM will be constructed from multiple genotype files. 0: computing grm from each file and then add all together. 1: join all genotypes, then compute grm.
  int useMPIForWriting;                 ///<GRM writing will be serialized if this is false. It could be useful for some file systems.
  bool diagonalizeGRM;                  ///<Diagonalize computed GRMs?
  bool writeAlsoNormalGRM;              ///<If diagonalizeGRMs == true, store also undiagonalized GRM
  double minimumNumberOverlapingSNPs;   ///<The minimum number of overlapping SNPs for considering correct the correlation between two individuals.
  double maximumProportionOfElimitaedIndividuals; ///<The maximum proportion of individuals allowed to be filtered due to not enouth SNP overlapping before disregarding a GRM.
  bool flatGRMNormalization;            ///<If true, the GRM is normalized assuming there are no missings.
  
  double mostRelatedLowerThreshold;     ///<All individuals with a relation below this threshold will be kept (makeGRMKeepMostRelated analysis).
  double mostRelatedUpperThreshold;     ///<All individuals with a relation above this threshold will be kept (makeGRMKeepMostRelated analysis).
  std::vector<double> pruneThresholdsCheck; ///<Check the number of individuals which remains after applying this prune threshold.
  
  bool addAllGRMsInOne;                 ///<When true, if several GRMs are loaded, all are added in a single GRM.
  
  ////////////////////////////
  // REML
  bool computeBLUE;                     ///<Compute the best linear unbiased estimator?
  bool computeIndividualsBLUP;          ///<Compute the best linear unbiased predictor for the total genetic effect?
  bool computeIndividualsBLUPErrors;    ///<Compute the best linear unbiased predictor errors for the total genetic effect?
  bool computeSNPsBLUP;                 ///<Compute the best linear unbiased predictor for each individual SNP?
  bool fixCorrelation;                  ///<Repeat REML with a fixed correlation?
  double fixedCorrelation;              ///<Value of the fixed correlation.
  int phenotypeColumn;                  ///<Column of the phenotype used in RMEL in file phenotypesFile.
  std::vector<int> phenotypeColumns;    ///<Columns of the phenotypes in file phenotypesFile to analyze in bivar/multi REML.
  bool analyzeAllPhenos;                ///<Iteratively analyze all phenotypes.
  int REMLMethod;                       ///<REML method 0->AI, 1->EM.
  bool environmentalCovariance;         ///<Environmental covariance must be computed on bivariate REML?
  bool computeREMLInSubsample;          ///<Compute REML in a subsample before REML with all samples.
  int nSubSampleIterations;             ///<Number of REML subsample computations.
  double initialSubsampleFraction;      ///<Proportion of individuals that will be used as a first REML approximation.
  int minimumSubsample;                 ///<Minimum number of individuals on random subsample.
  double initialh2Trait;                ///<The initial estimated heretability used for estimating initial variance values for trait in REML
  std::vector<double> initialh2Traits;  ///<The initial estimated heretabilities used for estimating initial variance values for traits when computing bivar/multi REML;
  double varianceConvergenceThreshold;  ///<The threshold for relative variance differences criteria for convergence.
  double gradientConvergenceThreshold;  ///<The threshold for gradient criteria for convergence.
  double changeAIStepThreshold;         ///<Threshold for switching from AI to short AI or EM
  double allowSwitchFromAItoEM;         ///<When loglikelihood relative difference > changeAIStepThreshold, a switch from AI to EM is performed when this is true. (needs testing)
  double firstStepEM;                   ///<Use the EM for the first step?
  std::string initialVariancesFile;     ///<Define a file with initial variances in.
  bool correctLinkageDisequilibrium;    ///<Correct for linkage disequillibrium in single REML?
  bool joinCovariatesVertically;        ///<If true, on bivar REML analysis, joinCovariates vertically.
  bool forceUseDiagonalizedKernels;     ///<Force the use of diagonalized GRMs for analysis that do not allow them. GRMs will be undiagonalized before use.
  bool computeEpistasisVariance;        ///<Compute the variance due to epistatic effects?
  double stepWeightingConstant;         ///<Constant used for weighting AI REML step when previous change in likelihood is large.
  std::string discreteRandomEffectsFile;///<The file containing discrete randome effects
  std::vector<int> discreteRandomEffectsFileColumns;          ///<A vector with the columns used from discreteRandomEffectsFile file.
  std::string multiDiscreteRandomEffectsFile;///<The file containing multiple discrete randome effects
  std::vector<int> multiDiscreteRandomEffectsFileColumns;          ///<A vector with the columns used from multiDiscreteRandomEffectsFile file.
  std::vector<std::string> continuousRandomEffectsFile;       ///<The files containing discrete randome effects
  bool skipComputeReducedModels;        ///<When true, reduced models will not be tested.
  bool skipTestGlobalGRMRegionalAnalysis;                     ///<In a regional analysis, test the Global GRM?
  std::vector<REMLCovarianceRestricted> restrictedCovariances;///<List and definitions of restricted covariances to use.
  std::string epistatitcPredictionGRMFile;  ///
  bool allowConvergenceWithConstrainedVars; ///<Allow REML convergence even when some variances are constrained.
  double maximumCorrelationCovarianceConstrain;               ///<The maximum correlation which a covariance can reach.
  bool useCorrelations;                 ///<On the multi REML, adjust the correlations instead of the covariances.
  bool useLogLogisticScale;             ///<On the REML optimization, use the log scale for variances and logistic for correlations.
  bool writeBLUEInReducedModels;        ///<Write BLUE results on reduced REML models?
  std::string indirectEffectsCouplesFile;   ///<The fie indicating the couples for computing indirect effects.
  bool computeDominanceGRM;             ///<Create a dominance GRM?
  std::vector<std::string> SNPBLUPGenotypeFiles;              ///<A file with a list of genotypes used for SNP BLUP computation.
  bool removeNonOvarlapingBLUPSNPs;     ///<Remove the SNPs present in the genotype file but not in the GRM when computing the BLUPs?
  bool useMLinsteadOfREML;              ///<Use ML instead of REML.
  std::vector<std::string> squaredExponentialKernelFiles;     ///<The files containing the coordinates used for computing distances for the squared exponential kernels.
  bool forceUseREMLAIWhenNoLinearCovariance;    ///<When non linear covariance matrices are used, avoid swith to use F matrix for iteration steps.
  double expKernelParameterInitialFactor;       ///<The factor used to multiply the initial estimated value for the parameter of an exponential covariance.
  int remlStepsToUnfixExpKernelParameter;       ///<Number of REML steps before the kernel parameter will be unfixed.
  std::string environmentalWeightsFile;         ///<File used to define environmental weights.
  int environmentalWeightsCol;                  ///<Column on the files from where environmental weights will be read.
  bool scaleEnvironmentalWeightTrace;           ///<If true, divide the weights so that the trace of the environmental matrix is equal to the dimension of the matrix.
  bool includeGxEInteractions;                  ///<Include GxE interactions?
  bool switchToEMLargeParameterChange;          ///<Switch to EM step when the parameter of a squared exponential kernel changes over a threshold.
  bool switchToEMParameterChangeThreshold;      ///<Switch to EM step when the parameter of a squared exponential kernel changes over a this threshold and switchToEMLargeParameterChange == true.
  std::vector<std::string> reducedModelsOnlyCovariances;       ///<A list of covariances used when testing reduced models. Reduced models were tested were only this covariances with uncorrelated environment ("E") are fited.
  
  ////////////////////////////
  // GWAS
  bool gwasWithMixedModel;              ///<Perform GWAS using a mixed model correcting for population structure.
  double correlatedSNPsThreshold;       ///<Threshold above which two SNPs are considered to be correlated enough for filtering one of two.
  std::string genotypeAndGRMsListFile;  ///<File with a list of genotype files and GRM files used for GWAS correcting for population structure.
  std::string genotypeBGENAndGRMsListFile; ///<File with a list of genotype (BGEN) files and GRM files used for GWAS correcting for population structure.
  std::string genotypeAndResidualsListFile;  ///<File with a list of genotype files and residuals files used for multiple phenotype GWAS.
  std::string genotypeBGENAndResidualsListFile; ///<File with a list of genotype (BGEN) files and GRM files used for multiple phenotype GWAS.
  bool parallelGWAS;                    ///<If true, use a grouped communicator, each group with a subgroup of SNPs.
  bool gwasUseAlwaysDiagonalCovariances;///<In a corrected GWAS analysis, diagonalize the covariance matrix if it is not previously diagonalized.
  bool computeGroupVariance;            ///<Compute the variance of the effects of each fitted group.
  bool saveGroupEffects;                ///<Compute and save the effects of the groups for each individual.
  
  ////////////////////////////
  // IGWAS
  std::string igwasTestCovariates;      ///<Covariates to test in a igwas as fixed effects
  std::string igwasTestQCovariates;     ///<Quantitative covariates to test in a igwas as fixed effects
  
  ////////////////////////////
  // MPGWAS
  bool originalRedistributionMethod;   ///<Specifies the redistribution method used in mpgwas.
  
  ////////////////////////////
  // Recursive GWAS
  double relationFitSNPsIndividuals;      ///<When performing a grouped recursive GWAS, the ratio of the number of SNPs to fit related to the number of samples.
  std::vector<double> significanceThresholdsFilterSNPs;         ///<Significance thresholds for filtering SNPs in a recursive GWAS.
  int recursiveGWASMaxIterations;         ///<The maximum number of iterations for recursive GWAS. If negative, there are no maximum.
  
  ////////////////////////////misc.message
  // Limits
  int maxREMLIterations;                ///<Max number of REML iterations
  double varianceConstrainProportion;   ///<If a variance is negative in a REML step. It will be constrained to a base variance multiplied by this variable.
  double minimumVarianceConstrainProportion;  ///<The REML algorithm can change the varianceConstrainProportion option. This sets a lower limit to this change.
  long randomSeed;                      ///<Seed used for random number generation.
  int logOutputPrecision;
  int resultsOutputPrecision;
  int logFieldWidth;
  int resultsFieldWidth;
  double grmCutoff;                     ///<Threshold over which individuals will be pruned in the grm
  bool pruneGRM;                        ///<Prune GRMs?
  double maximumProportionOfCategories; ///<Var that defines the maximum proportion fo categories for fixed effects compared to sample size.
  
  ////////////////////////////
  // Phenotype simulation vars
  std::string effectsSizeFile;          ///<file of SNP effects
  std::string adjustEffectsGenotypeListFile; ///<genotype file for adjusting SNP effects using their frequencies.
  double simulationh2;                  ///<h square for the simulation.
  double prevalence;                    ///<Prevalence of the binary trait
  bool simulateBinaryTrait;             ///<The trait is binary?
  bool simulateQuantitativeTrait;       ///<The trait is quantitative?
  
  ////////////////////////////
  // Phenotype prediction vars
  std::string snpEffectsFile;           ///<file with SNP effects.
  std::string fileCovarEffects;         ///<file with discrete covariate effects.
  std::string fileQCovarEffects;        ///<file with quantitative covariate effects.
  bool forceUseNonEstimatedCovarKeys;   ///<When estimating covariate effects, allow using non-estimated effect keys as 0.
  
  ////////////////////////////
  // PCA
  int nEigenValues;                     ///<The number of eigenvalues will be stored.
  
  ////////////////////////////
  // Regional analysis
  bool regionalAnalysis;                ///<Perform regional analysis?
  GroupBy regionBy;                     ///<How regions will be defined?
  std::string regionsFile;              ///<SNPs are grouped by groups defined in this file. This could be used in REML regional analysis or GWAS.
  int regionSize;                       ///<In a regional analysis, the regions (defined by position) are of this size.
  int regionOverlapping;                ///<In a regional analysis, the regions (defined by position) have this overlapping.
  int minSNPsInRegions;                 ///<The minimum SNPs ina region for computing the region.
  int fixedGroupSize;                   ///<The size of the group when SNPs grouped using ordered fixed sized groups.
  bool allRegionsTogether;              ///<In a regional analysis, if true, all regions are fitted together.
  
  ////////////////////////////
  // Group effects analysis
  std::string fileGroupsToKeep;         ///<File indicating groups to keep. The others will be filtered.
  std::string groupsPositions;          ///<File with positions of the groups.
  double groupDistanceForDiscarding;    ///<Maximum distance for discarding one group of a correlated group pair.
  std::vector<std::string> groupEffectsFileList;  ///<List of files of group effects.
  std::vector<std::string> groupEffectsPairFileList;  ///<List of paired files of group effects.
  bool crossedEffectCorrelations;
  
  ////////////////////////////
  // Filter Labeled Matrix
  std::string inputLabeledMatrix;
  std::string rowLabelsFile;
  std::string columnLabelsFile;
  
  ////////////////////////////
  // Divers
  bool removeMissings;                  ///<If true, individuals with a missing value on any fixed effect will be removed.
  bool remlGCTAMode;
  int defaultBlockSize;                 ///<Set the default matrix distribution block sizes.
  bool allowSinglePrecisionInversion;   ///<Allow single precision on matrix inversions.
  std::string fileSNPsToKeep;           ///<When loading a genotype file, keep only SNPs in this file.
  std::string fileIndividualsToKeep;    ///<Keep only those individuals (currently only used for the testing effects analysis, and filtering phenotypes).
  std::set<std::string> keysIndividualsToKeep;     ///<If it is not empty, keep only those individuals in this set (currently only used for the testing effects analysis, and filtering phenotypes).
  bool allowFixingVariancesToZero;      ///<Allow variances to be fixed to zero when performing a REML.
  std::map< std::string, std::string> baseVarianceNames;  ///<Used output base variance names.
  int numberOfProcessesDefaultMPIGroup; ///<Number of processes on the default MPI group communicator (basicGroupedCommunicator).
  bool verboseLog;                      ///<Increase log outputs.
  bool allowBGENLayout1;                ///<Allow reading files with BGEN Layout 1.
  bool keepSNPsWithZeroVariance;        ///<Keep SNPs with zero variance when only computing allele frequencies in imputed variants.
  
  ////////////////////////////
  // Debug options
  bool mpiDebug;
  bool debug;
  
  ////////////////////////////
  // Parsed options
  std::map<std::string, std::string> parsedOptions;     ///<Parsed options
  std::string warnings;                                 ///<warnings
  
  Options();
  Options(int argc, char **argv);
  ~Options();
  
  void showOptions();
  
  void defaultOptions();
  void parseOptions(int argc, char **argv);
  void checkIncompatibilities();
  void fillMissingOptions();
  
  void setAnalysis(AnalysisToPerform analysisToSet);
  std::string getString(int argc, char **argv, int i);
  std::string getFileName(int argc, char **argv, int i);
  int getInt(int argc, char **argv, int i, Range range = Range());
  double getDouble(int argc, char **argv, int i, Range range = Range());
  std::vector<std::string> getStringList(int argc, char **argv, int i);
  std::vector<std::string> getFileNameList(int argc, char **argv, int i, bool checkExistence = true);
  std::vector<int> getIntList(int argc, char **argv, int i, Range range);
  std::vector<double> getDoubleList(int argc, char **argv, int i, Range range);
  
  std::vector<REMLCovarianceRestricted> getRestrictedCovariancesList( std::vector<std::string> listOfArgs);
  
  void appendParsedOptions(std::string option, std::string parameters);
  void showParsedOptions();
  
  /**
   * Check if passed GRMs are Diagonalized
   * 
   * \param fn file name of the grm
   * \return true if grm is diagonal, false otherwise.
   */
  bool isGRMDiagonal(std::string fn);
};

#endif
