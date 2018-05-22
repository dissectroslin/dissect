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

#include "options.h"
#include "misc.h"
#include "auxiliar.h"

#include <iostream>
#include <sstream>
#include <map>
#include <ctime>
#include <cstring>

Options::Options()
{
  defaultOptions();
}

Options::Options(int argc, char **argv)
{
  defaultOptions();
  
  parseOptions(argc, argv);
  
  checkIncompatibilities();
}

Options::~Options()
{
}

void Options::defaultOptions()
{
  //////////////////////////////////////////////////
  // Actions
  this->analysis = withoutAnalysis;

  //////////////////////////////////////////////////
  // Base input files
  this->genotypeFile = "";
  this->genotypeListFile = "";
  this->bgenGenotypeFile = "";
  this->grmFile = "";
  this->grmListFile = "";
  this->phenotypesFile = "";
  this->qCovarsFile = "";
  this->qCovarsFiles.clear();
  this->covarsFile = "";
  this->covarsFiles.clear();
  
  //////////////////////////////////////////////////
  // Output files
  this->outFile = "dissect";
  this->compressedOutput = false;
  
  //////////////////////////////////////////////////
  // GRM
  this->GRMJoinMethod = 0;
  this->useMPIForWriting = true;
  this->diagonalizeGRM = false;
  this->writeAlsoNormalGRM = false;
  this->minimumNumberOverlapingSNPs = 10.;
  this->maximumProportionOfElimitaedIndividuals = 0.9;
  this->flatGRMNormalization = false;
  
  this->mostRelatedLowerThreshold = -1.;
  this->mostRelatedUpperThreshold = 0.025;
  this->pruneThresholdsCheck.clear();
  
  this->addAllGRMsInOne = false;
  
  this->GCTAGRMsFile = "";
  
  //////////////////////////////////////////////////
  // REML
  this->computeBLUE = false;
  this->computeIndividualsBLUP = false;
  this->computeIndividualsBLUPErrors = false;
  this->computeSNPsBLUP = false;
  this->fixCorrelation = false;
  this->fixedCorrelation = 0.;
  this->phenotypeColumn = 1;
  this->phenotypeColumns.clear();
  this->REMLMethod = 0;
  this->environmentalCovariance = true;
  this->computeREMLInSubsample = false;
  this->nSubSampleIterations = 1;
  this->initialSubsampleFraction = 0.3;
  this->minimumSubsample = 1000;
  this->initialh2Trait = 0.5;
  this->initialh2Traits.clear();
  this->varianceConvergenceThreshold = 1e-5;
  this->gradientConvergenceThreshold = 1e-2;
  this->changeAIStepThreshold = 0.001;
  this->allowSwitchFromAItoEM = false;
  this->firstStepEM = true;
  this->initialVariancesFile = "";
  this->correctLinkageDisequilibrium = false;
  this->joinCovariatesVertically = false; //The oposite situation is untested
  this->forceUseDiagonalizedKernels = false;
  this->computeEpistasisVariance = false;
  this->stepWeightingConstant = 0.3;
  this->discreteRandomEffectsFile = "";
  this->discreteRandomEffectsFileColumns.clear();
  this->multiDiscreteRandomEffectsFile = "";
  this->multiDiscreteRandomEffectsFileColumns.clear();
  this->continuousRandomEffectsFile.clear();
  this->skipComputeReducedModels = false;
  this->skipTestGlobalGRMRegionalAnalysis = false;
  this->restrictedCovariances.clear();
  this->epistatitcPredictionGRMFile = "";
  this->allowConvergenceWithConstrainedVars = true;
  this->maximumCorrelationCovarianceConstrain = 1.;
  this->useCorrelations = false;
  this->useLogLogisticScale = false;
  this->writeBLUEInReducedModels = false;
  this->indirectEffectsCouplesFile = "";
  this->computeDominanceGRM = false;
  this->SNPBLUPGenotypeFiles.clear();
  this->removeNonOvarlapingBLUPSNPs = true;
  this->useMLinsteadOfREML = false;
  this->squaredExponentialKernelFiles.clear();
  this->forceUseREMLAIWhenNoLinearCovariance = true;
  this->expKernelParameterInitialFactor = 1.;
  this->remlStepsToUnfixExpKernelParameter = 8;
  this->environmentalWeightsFile = "";
  this->environmentalWeightsCol = 1;
  this->scaleEnvironmentalWeightTrace = true;
  this->includeGxEInteractions = false;
  this->switchToEMLargeParameterChange = false;
  this->switchToEMParameterChangeThreshold = 0.5;
  this->reducedModelsOnlyCovariances.clear();
  
  //////////////////////////////////////////////////
  // Recursive GWAS
  this->relationFitSNPsIndividuals = 0.01;
  this->significanceThresholdsFilterSNPs = std::vector<double>(1, 0.05);
  this->recursiveGWASMaxIterations = -1;

  //////////////////////////////////////////////////
  // GWAS
  this->correlatedSNPsThreshold = -1.;
  this->genotypeAndGRMsListFile = "";
  this->genotypeBGENAndGRMsListFile = "";
  this->parallelGWAS = true;
  this->gwasUseAlwaysDiagonalCovariances = true;
  this->computeGroupVariance = false;
  this->saveGroupEffects = false;
  
  //////////////////////////////////////////////////
  // IGWAS
  this->igwasTestCovariates = "";
  this->igwasTestQCovariates = "";
  
  ////////////////////////////
  // MPGWAS
  this->originalRedistributionMethod = true;
  
  //////////////////////////////////////////////////
  // Limits
  this->maxREMLIterations = 40;
  this->varianceConstrainProportion = 1e-6;
  this->minimumVarianceConstrainProportion = 1e-10;
  this->randomSeed = time(NULL);
  this->logOutputPrecision = 4;
  this->resultsOutputPrecision = 6;
  this->logFieldWidth = this->logOutputPrecision + 10;
  this->logFieldWidth = (this->logFieldWidth < 15)?15:this->logFieldWidth;
  this->resultsFieldWidth = this->resultsOutputPrecision + 10;
  this->resultsFieldWidth = (this->resultsFieldWidth < 15)?15:this->resultsFieldWidth;
  this->grmCutoff = 0.025;
  this->pruneGRM = false;
  this->maximumProportionOfCategories = 0.1;
  
  //////////////////////////////////////////////////
  // PCA
  this->nEigenValues = 20;
  
  //////////////////////////////////////////////////
  // Phenotype simulation vars
  this->effectsSizeFile = "";
  this->adjustEffectsGenotypeListFile = "";
  this->simulationh2 = 0.1;
  this->prevalence = 0.1;
  this->simulateBinaryTrait = false;
  this->simulateQuantitativeTrait = false;
  
  //////////////////////////////////////////////////
  // Phenotype prediction vars
  this->snpEffectsFile = "";
  this->fileCovarEffects = "";
  this->fileQCovarEffects = "";
  this->forceUseNonEstimatedCovarKeys = false;
  
  //////////////////////////////////////////////////
  // Regional analysis
  this->regionalAnalysis = false;
  this->regionsFile = "";
  this->regionBy = ungrouped;
  this->regionSize = -1;
  this->regionOverlapping = 0;
  this->minSNPsInRegions = 1;
  this->fixedGroupSize = -1;
  this->allRegionsTogether = false;
  
  ////////////////////////////
  // Filter Labeled Matrix
  this->inputLabeledMatrix = "";
  this->rowLabelsFile = "";
  this->columnLabelsFile = "";
  
  //////////////////////////////////////////////////
  // Group effects analysis
  this->fileGroupsToKeep = "";
  this->groupsPositions = "";
  this->groupDistanceForDiscarding = 500000;    ///<Maximum distance for discarding one group of a correlated group pair.
  this->groupEffectsFileList = std::vector<std::string>();
  this->groupEffectsPairFileList.clear();  ///<List of paired files of group effects.
  this->crossedEffectCorrelations = false;
  
  //////////////////////////////////////////////////
  // Divers
  this->analyzeAllPhenos = false;
  this->removeMissings = false;
  this->remlGCTAMode = false;
  this->defaultBlockSize = 256;
  this->allowSinglePrecisionInversion = true;
  this->fileSNPsToKeep = "";
  this->fileIndividualsToKeep = "";
  this->keysIndividualsToKeep.clear();
  this->allowFixingVariancesToZero = false;
  this->numberOfProcessesDefaultMPIGroup = 1;
  this->verboseLog = false;
  this->allowBGENLayout1 = false;
  this->keepSNPsWithZeroVariance = false;
  
  this->baseVarianceNames["GRM"] = "GRM";
  this->baseVarianceNames["epi"] = "epi";
  this->baseVarianceNames["Global-GRM"] = "Global-GRM";
  this->baseVarianceNames["Regional-GRM"] = "Regional-GRM";
  this->baseVarianceNames["cov-"] = "cov-";
  this->baseVarianceNames["SEK-"] = "SEK-";
  
  //////////////////////////////////////////////////
  // Debugging vars
  this->mpiDebug = false;
  this->debug = false;
}

void Options::parseOptions(int argc, char **argv)
{
  std::vector<std::string> tempRestrictedCovariances;
  
  for (int i = 1; i < argc; i++) 
  {
    //////////////////////////////////////////////////
    // Analysis
    
    //Just check the options?
    if( !strcmp(argv[i],"--check") )
    {
      setAnalysis(justCheckAnalysis);
    }
    //Show help?
    else if( !strcmp(argv[i],"--help") || !strcmp(argv[i],"-h") )
    {
      setAnalysis(showHelpAnalysis);
    }
    //Only compute grm
    else if( !strcmp(argv[i],"--make-grm") )
    {
      setAnalysis(makeGRMAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //
    else if( !strcmp(argv[i],"--make-grm-mr") )
    {
      setAnalysis(makeGRMKeepMostRelated);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform univariate REML
    else if( !strcmp(argv[i],"--reml") )
    {
      setAnalysis(REMLAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform bivariate REML
    else if( !strcmp(argv[i],"--bivar-reml") )
    {
      setAnalysis(bivarREMLAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform bivariate REML
    else if( !strcmp(argv[i],"--multi-reml") )
    {
      setAnalysis(multivarREMLAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Simulate Phenotypes
    else if( !strcmp(argv[i],"--simulate") )
    {
      setAnalysis(simulatePhenotypeAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Predict phenotypes
    else if( !strcmp(argv[i],"--predict") )
    {
      setAnalysis(predictPhenotypeAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform PCA
    else if( !strcmp(argv[i],"--pca") )
    {
      setAnalysis(PCAAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform GWAS
    else if( !strcmp(argv[i],"--gwas") )
    {
      setAnalysis(GWASAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform recursive grouped GWAS
    else if( !strcmp(argv[i],"--rgwas") )
    {
      setAnalysis(recursiveGWASAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform inverse GWAS
    else if( !strcmp(argv[i],"--igwas") )
    {
      setAnalysis(iGWASAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Compute the accuracy as a function of SNPs included in the predictor. Experimental analysis.
    else if( !strcmp(argv[i],"--accuracy-by-snp") )
    {
      setAnalysis(accuracyBySNPAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform analysis from gene effects matrices
    else if( !strcmp(argv[i],"--effects") )
    {
      setAnalysis(effectsAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform GLMM analysis
    else if( !strcmp(argv[i],"--glmm") )
    {
      setAnalysis(GLMMAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform Phenotype prediction from covars
    else if( !strcmp(argv[i],"--cov-predict") )
    {
      setAnalysis(predictCovarPhenotypeAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform MultiPhenotype GWAS
    else if( !strcmp(argv[i],"--mpgwas") )
    {
      setAnalysis(multiplePhenotypeGWASAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Perform MultiPhenotype residuals
    else if( !strcmp(argv[i],"--mpresiduals") )
    {
      setAnalysis(multiplePhenotypeResiduals);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //
    else if( !strcmp(argv[i],"--filter-matrix") )
    {
      setAnalysis(filterLabeledMatrix);
      appendParsedOptions(std::string(argv[i]), "");
    }
    //
    else if( !strcmp(argv[i],"--snp-stats") )
    {
      setAnalysis(computeSNPStatsAnalysis);
      appendParsedOptions(std::string(argv[i]), "");
    }
    
    //////////////////////////////////////////////////
    // Base input files
    
    //Genotype file
    else if ( !strcmp(argv[i],"--bfile") )
    {
      this->genotypeFile = getString(argc, argv, i + 1);
      misc.checkFileExists(this->genotypeFile + ".bed");
      misc.checkFileExists(this->genotypeFile + ".bim");
      misc.checkFileExists(this->genotypeFile + ".fam");
    }
    //File with a list of genotype files
    else if( !strcmp(argv[i],"--bfile-list") )
    {
      this->genotypeListFile = getFileName(argc, argv, i + 1);
    }
    //bgen genotype file
    else if( !strcmp(argv[i],"--bgen") )
    {
      this->bgenGenotypeFile = getString(argc, argv, i + 1);
      misc.checkFileExists(this->bgenGenotypeFile + ".bgen");
    }
    //GRM file
    else if( !strcmp(argv[i],"--grm") )
    {
      this->grmFile = getString(argc, argv, i + 1);
      misc.checkFileExists(this->grmFile + ".grm.ids");
      misc.checkFileExists(this->grmFile + ".grm.dat");
      misc.checkFileExists(this->grmFile + ".grm.snps");
    }
    //File with a list of GRM files
    else if( !strcmp(argv[i],"--grm-list") )
    {
      this->grmListFile = getFileName(argc, argv, i + 1);
    }
    //Phenotypes file
    else if( !strcmp(argv[i],"--pheno") )
    {
      this->phenotypesFile = getFileName(argc, argv, i + 1);
    }
    //Phenotypes files
    else if( !strcmp(argv[i],"--phenos") )
    {
      this->phenotypesFiles = getFileNameList(argc, argv, i + 1);
    }
    //Discrete covariance file 1
    else if( !strcmp(argv[i],"--covar") )
    {
      this->covarsFile = getFileName(argc, argv, i + 1);
    }
    //Discrete covariance file 2
    else if( !strcmp(argv[i],"--covars") )
    {
      this->covarsFiles = getFileNameList(argc, argv, i + 1);
    }
    //Quantitative covariance file 1
    else if( !strcmp(argv[i],"--qcovar") )
    {
      this->qCovarsFile = getFileName(argc, argv, i + 1);
    }
    //Quantitative covariance file 2
    else if( !strcmp(argv[i],"--qcovars") )
    {
      this->qCovarsFiles = getFileNameList(argc, argv, i + 1);
    }
    //File defining GRMs in GCTA format.
    else if( !strcmp(argv[i],"--gcta-grms-gz") )
    {
#if defined(BOOSTLIB) && defined(ZLIB)
      this->GCTAGRMsFile = getFileName(argc, argv, i + 1);
#else
      misc.error("Error: The current version of DISSECT is compiled without zlib and boost support. These libraries are required for using the option --gcta-grms-gz. Please, recompile it using these libraries or ask for support.", 0);
#endif
    }
    
    //////////////////////////////////////////////////
    // Output files
    
    //Output file
    else if( !strcmp(argv[i],"--out") )
    {
      this->outFile = getString(argc, argv, i + 1);
    }
    else if( !strcmp(argv[i],"--zout") )
    {
#if defined(BOOSTLIB) && defined(ZLIB)      
      this->outFile = getString(argc, argv, i + 1);
      this->compressedOutput = true;
#else
      misc.error("Error: The current version of DISSECT is compiled without zlib and boost support. These libraries are required for using the --zout option. Please, recompile it using these libraries or ask for support.", 0);
#endif
    }
    
    //////////////////////////////////////////////////
    // GRM
    else if( !strcmp(argv[i],"--grm-join-method") )
    {
      this->GRMJoinMethod = getInt(argc, argv, i + 1, Range(0, 1));
    }
    else if( !strcmp(argv[i],"--grm-no-mpi-write") )
    {
      this->useMPIForWriting = false;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--diagonalize") )
    {
      this->diagonalizeGRM = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--store-both") )
    {
      this->writeAlsoNormalGRM = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--min-overlap-snps") )
    {
      this->minimumNumberOverlapingSNPs = double( getInt(argc, argv, i + 1, Range(0, false)) );
    }
    
    else if( !strcmp(argv[i],"--mostr-lower-thr") )
    {
      this->mostRelatedLowerThreshold = getDouble(argc, argv, i + 1, Range());
    }
    else if( !strcmp(argv[i],"--mostr-upper-thr") )
    {
      this->mostRelatedUpperThreshold = getDouble(argc, argv, i + 1, Range());
    }
    else if( !strcmp(argv[i],"--cutoff-thrs") )
    {
      this->pruneThresholdsCheck = getDoubleList(argc, argv, i + 1, Range());
    }
    else if( !strcmp(argv[i],"--grm-flat-norm") )
    {
      this->flatGRMNormalization = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--add-grms") )
    {
      this->addAllGRMsInOne = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Specify the maximum proportion of individuals of the GRM that can be filtered due to missing genotypes without discarding the GRM.
    else if( !strcmp(argv[i],"--min-prop-grm-inds-kept") )
    {
      this->maximumProportionOfElimitaedIndividuals = getDouble(argc, argv, i + 1, Range(0, 1));
    }
    
    //////////////////////////////////////////////////
    // REML
    else if( !strcmp(argv[i],"--blue") )
    {
      this->computeBLUE = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--indiv-blup") )
    {
      this->computeIndividualsBLUP = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--indiv-blup-error") )
    {
      this->computeIndividualsBLUPErrors = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--snp-blup") )
    {
      this->computeSNPsBLUP = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--fix-correlation") )
    {
      tempRestrictedCovariances = getStringList(argc, argv, i + 1);
      //this->fixedCorrelation = getDouble(argc, argv, i + 1, Range(-1., 1.));
      //this->fixCorrelation = true;
    }
    //Get the columns should be used for phenotypes
    else if( !strcmp(argv[i],"--pheno-col") )
    {
      this->phenotypeColumn = getInt(argc, argv, i + 1, Range(1, false));
    }
    else if( !strcmp(argv[i],"--pheno-cols") )
    {
      this->phenotypeColumns = getIntList(argc, argv, i + 1, Range(1, false));
    }
    //Get the REML method
    else if( !strcmp(argv[i],"--reml-method-em") )
    {
      //this->REMLMethod = 1;
      //appendParsedOptions(std::string(argv[i]), "");
      misc.error("Error: Sorry, --reml-method-em is temporarily not available.", 0);
    }
    else if( !strcmp(argv[i],"--reml-method-ai") )
    {
      this->REMLMethod = 0;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--no-environment-cov") )
    {
      this->environmentalCovariance = false;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--reml-subsample") )
    {
      this->computeREMLInSubsample = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Get the number of subsample iterations to perform
    else if( !strcmp(argv[i],"--subsample-replicates") )
    {
      this->nSubSampleIterations = getInt(argc, argv, i + 1, Range(1, 100));
    }
    else if( !strcmp(argv[i],"--initial-h2") )
    {
      this->initialh2Trait = getDouble(argc, argv, i + 1, Range(0., 1.));
    }
    else if( !strcmp(argv[i],"--initial-h2s") )
    {
      this->initialh2Traits = getDoubleList(argc, argv, i + 1, Range(0., 1.));
    }
    //The threshold for relative variance differences criteria for convergence.
    else if( !strcmp(argv[i],"--variance-threshold") )
    {
      this->varianceConvergenceThreshold = getDouble(argc, argv, i + 1, Range(0., 0.1));
    }
    //The threshold for gradient criteria for convergence.
    else if( !strcmp(argv[i],"--gradient-threshold") )
    {
      this->gradientConvergenceThreshold = getDouble(argc, argv, i + 1, Range(0., false));
    }
    //AI reml step is modified when loglikelihood relative difference > this->changeAIStepThreshold.
    else if( !strcmp(argv[i],"--ai-switch-threshold") )
    {
      this->changeAIStepThreshold = getDouble(argc, argv, i + 1, Range(0., 1.));
    }
    //When loglikelihood relative difference > this->changeAIStepThreshold, a switch from AI to EM is performed when this is true.
    else if( !strcmp(argv[i],"--ai-em-switch") )
    {
      this->allowSwitchFromAItoEM = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Use EM in the first step?
    else if( !strcmp(argv[i],"--no-first-em") )
    {
      this->firstStepEM = false;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Give a file with initial variances.
    else if( !strcmp(argv[i],"--initial-variances") )
    {
      this->initialVariancesFile = getFileName(argc, argv, i + 1);
    }
    //Try on single REML for correcting for LD. Experimental.
    else if( !strcmp(argv[i],"--correct-ld") )
    {
      this->correctLinkageDisequilibrium = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Force the use of diagonalized GRMs for analysis that do not allow them. GRMs will be undiagonalized before use.
    else if( !strcmp(argv[i],"--force-use-diag-kernels") )
    {
      this->forceUseDiagonalizedKernels = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Compute epistasis?
    else if( !strcmp(argv[i],"--epistasis-var") )
    {
      this->computeEpistasisVariance = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Define the scaling constant for the delta of a q step.
    else if( !strcmp(argv[i],"--reml-qstep-scale") )
    {
      this->stepWeightingConstant = getDouble(argc, argv, i + 1, Range(0., 1.));
    }
    //Define file with discrete values for computing random effects
    else if( !strcmp(argv[i],"--random-effects") )
    {
      this->discreteRandomEffectsFile = getFileName(argc, argv, i + 1);
    }
    //Define which columns of file this->discreteRandomEffectsFile will be used
    else if( !strcmp(argv[i],"--random-effects-cols") )
    {
      this->discreteRandomEffectsFileColumns = getIntList(argc, argv, i + 1, Range(1, false));
    }
    //Define file with multiple discrete values for computing random effects
    else if( !strcmp(argv[i],"--multirandom-effects") )
    {
      this->multiDiscreteRandomEffectsFile = getFileName(argc, argv, i + 1);
    }
    //Define which columns of file this->multiDiscreteRandomEffectsFile will be used
    else if( !strcmp(argv[i],"--multirandom-effects-cols") )
    {
      this->multiDiscreteRandomEffectsFileColumns = getIntList(argc, argv, i + 1, Range(1, false));
    }
    else if( !strcmp(argv[i],"--skip-test-reduced-models") )
    {
      this->skipComputeReducedModels = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--grm-epi") )
    {
      this->epistatitcPredictionGRMFile = getString(argc, argv, i + 1);
      misc.checkFileExists(this->epistatitcPredictionGRMFile + ".grm.ids");
      misc.checkFileExists(this->epistatitcPredictionGRMFile + ".grm.dat");
      misc.checkFileExists(this->epistatitcPredictionGRMFile + ".grm.snps");
    }
    //The maximum correlation to which covariances will be constrained
    else if( !strcmp(argv[i],"--max-correlation") )
    {
      this->maximumCorrelationCovarianceConstrain = getDouble(argc, argv, i + 1, Range(0., false));
    }
    //Use correlations instead of covariances when optimizing REML
    else if( !strcmp(argv[i],"--use-correlations") )
    {
      this->useCorrelations = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Optimize REML on the log and logistic scale
    else if( !strcmp(argv[i],"--use-log-logistic") )
    {
      this->useLogLogisticScale = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--write-blue-reduced") )
    {
      this->writeBLUEInReducedModels = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Indirect effects file
    else if( !strcmp(argv[i],"--indirect-effects-couples") )
    {
      this->indirectEffectsCouplesFile = getFileName(argc, argv, i + 1);
    }
    //A file with a list of genotypes used for SNP BLUP computation.
    else if( !strcmp(argv[i],"--blup-bfile-list") )
    {
      std::string tempFileName = getFileName(argc, argv, i + 1);
      getListFromFile(tempFileName, this->SNPBLUPGenotypeFiles);
    }
    //Remove the SNPs present in the genotype file but not in the GRM when computing the BLUPs?
    else if( !strcmp(argv[i],"--blup-no-filter-snps") )
    {
      this->removeNonOvarlapingBLUPSNPs = false;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Use ML instead of REML?
    else if( !strcmp(argv[i],"--use-ml") )
    {
      this->useMLinsteadOfREML = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Squared exponential kernel coordinate files
    else if( !strcmp(argv[i],"--sqrt-exp-coord-files") )
    {
      this->squaredExponentialKernelFiles = getFileNameList(argc, argv, i + 1);
      this->switchToEMLargeParameterChange = true;
    }
    //The factor used to multiply the initial estimated value for the parameter of an exponential covariance.
    else if( !strcmp(argv[i], "--param-init-fac") )
    {
      this->expKernelParameterInitialFactor = getDouble(argc, argv, i + 1, Range());
    }
    //Number of REML steps before the kernel parameter will be unfixed.
    else if( !strcmp(argv[i],"--steps-to-unfix") )
    {
      this->remlStepsToUnfixExpKernelParameter = getInt(argc, argv, i + 1, Range(0, false));
    }
    //File used to define environmental weights.
    else if( !strcmp(argv[i],"--weights") )
    {
      this->environmentalWeightsFile = getFileName(argc, argv, i + 1);
    }
    //Column on the files from where environmental weights will be read.
    else if( !strcmp(argv[i],"--weights-col") )
    {
      this->environmentalWeightsCol = getInt(argc, argv, i + 1, Range(1, false));
    } 
    //If true, divide the weights so that the trace of the environmental matrix is equal to the dimension of the matrix.
    else if( !strcmp(argv[i],"--no-scale-weights") )
    {
      this->scaleEnvironmentalWeightTrace = false;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Include GxE interactions?
    else if( !strcmp(argv[i],"--gxe") )
    {
      this->includeGxEInteractions = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Include covariances ids used to fit reduced models with only this variances individually (together with environmental variance).
    else if( !strcmp(argv[i],"--reduced-with-only") )
    {
      this->reducedModelsOnlyCovariances = getStringList(argc, argv, i + 1);
    }
    
    //////////////////////////////////////////////////
    // Recursive GWAS
    else if( !strcmp(argv[i],"--rgwas-ratio") )
    {
      this->relationFitSNPsIndividuals = getDouble(argc, argv, i + 1, Range(0.0000001, 0.1));
    }
    else if( !strcmp(argv[i],"--rgwas-thresholds") )
    {
      this->significanceThresholdsFilterSNPs = getDoubleList(argc, argv, i + 1, Range(0., 1.));
    }
    else if( !strcmp(argv[i],"--rgwas-maxit") )
    {
      this->recursiveGWASMaxIterations = getInt(argc, argv, i + 1, Range(1, false));
    }
    
    //////////////////////////////////////////////////
    // GWAS
    else if( !strcmp(argv[i],"--snp-corr-threshold") )
    {
      this->correlatedSNPsThreshold = getDouble(argc, argv, i + 1, Range(0., 1.));
    }
    else if( !strcmp(argv[i],"--bfile-grm-list") )
    {
      this->genotypeAndGRMsListFile = getFileName(argc, argv, i + 1);
    }
    else if( !strcmp(argv[i],"--bgen-grm-list") )
    {
#if defined(BGEN) && defined(ZLIB)
      this->genotypeBGENAndGRMsListFile = getFileName(argc, argv, i + 1);
#else
      misc.error("Error: The current version of DISSECT is compiled without zlib and bgen support. These libraries are required for using the option --bgen-grm-list. Please, recompile it using these libraries or ask for support.", 0);
#endif
    }
    else if( !strcmp(argv[i],"--bfile-residuals-list") )
    {
      this->genotypeAndResidualsListFile = getFileName(argc, argv, i + 1);
    }
    else if( !strcmp(argv[i],"--bgen-residuals-list") )
    {
#if defined(BGEN) && defined(ZLIB)
      this->genotypeBGENAndResidualsListFile = getFileName(argc, argv, i + 1);
#else
      misc.error("Error: The current version of DISSECT is compiled without zlib and bgen support. These libraries are required for using the option --bgen-residuals-list. Please, recompile it using these libraries or ask for support.", 0);
#endif
    }
    else if( !strcmp(argv[i],"--parallel-gwas") )
    {
      this->parallelGWAS = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--nonparallel-gwas") )
    {
      this->parallelGWAS = false;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Compute the variance of the effects of each fitted group.
    else if( !strcmp(argv[i],"--group-var") )
    {
      this->computeGroupVariance = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Compute and save the effects of the groups for each individual.
    else if( !strcmp(argv[i],"--group-effects") )
    {
      this->saveGroupEffects = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    
    //////////////////////////////////////////////////
    // IGWAS
    else if( !strcmp(argv[i],"--igwas-covar") )
    {
      this->igwasTestCovariates = getFileName(argc, argv, i + 1);
    }
    else if( !strcmp(argv[i],"--igwas-qcovar") )
    {
      this->igwasTestQCovariates = getFileName(argc, argv, i + 1);
    }
    
    //////////////////////////////////////////////////
    // MPGWAS
    else if( !strcmp(argv[i],"--redist-meth2") )
    {
      this->originalRedistributionMethod = false;
      appendParsedOptions(std::string(argv[i]), "");
    }
    
    //////////////////////////////////////////////////
    // Limits
    
    //How many REML iterations can be done without convergence?
    else if( !strcmp(argv[i],"--reml-maxit") )
    {
      this->maxREMLIterations = getInt(argc, argv, i + 1);
    }
    else if( !strcmp(argv[i],"--variance-constrain") )
    {
      this->varianceConstrainProportion = getDouble(argc, argv, i + 1);
    }
    else if( !strcmp(argv[i],"--random-seed") )
    {
      this->randomSeed = long(getInt(argc, argv, i + 1, Range(0, false)));
    }
    //Which elements of GRM should be deleted?
    else if( !strcmp(argv[i],"--grm-cutoff") )
    {
      this->grmCutoff = getDouble(argc, argv, i + 1, Range(-1., 2.));
      this->pruneGRM = true;
    }
    
    //////////////////////////////////////////////////
    // PCA
    
    //The number of eigenvectors that will be stored
    else if( !strcmp(argv[i],"--num-eval") )
    {
      this->nEigenValues = getInt(argc, argv, i + 1, Range(1, false));
    }
    
    //////////////////////////////////////////////////
    // Phenotype simulation vars
    
    //Get the effect sizes file name.
    else if( !strcmp(argv[i],"--effect-sizes") )
    {
      this->effectsSizeFile = getFileName(argc, argv, i + 1);
    }
    //genotype file for adjusting SNP effects using their frequencies.
    else if( !strcmp(argv[i], "--adjust-bfile-list") )
    {
      this->adjustEffectsGenotypeListFile = getFileName(argc, argv, i + 1);
    }
    //Which is the h2 for simulated phenotypes?
    else if( !strcmp(argv[i],"--simu-h2") )
    {
      this->simulationh2 = getDouble(argc, argv, i + 1, Range(0., 1.));
    }
    //Prevalence in a binary trait simulation.
    else if( !strcmp(argv[i],"--prevalence") )
    {
      this->prevalence = getDouble(argc, argv, i + 1, Range(0., 1.));
    }
    //The trait is quantitative?
    else if( !strcmp(argv[i],"--simu-quantitative") )
    {
      this->simulateQuantitativeTrait = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //The trait is binary?
    else if( !strcmp(argv[i],"--simu-binary") )
    {
      this->simulateBinaryTrait = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    
    //////////////////////////////////////////////////
    // Phenotype prediction vars
    //Get the effect sizes file name.
    else if( !strcmp(argv[i],"--snp-effects") )
    {
      this->snpEffectsFile = getFileName(argc, argv, i + 1);
    }
    //file with discrete covariate effects.
    else if( !strcmp(argv[i],"--covar-effects") )
    {
      this->fileCovarEffects = getFileName(argc, argv, i + 1);
    }
    //file with quantitative covariate effects.
    else if( !strcmp(argv[i],"--qcovar-effects") )
    {
      this->fileQCovarEffects = getFileName(argc, argv, i + 1);
    }
    //Can we use covariate with keys with unestimaed effects?
    else if( !strcmp(argv[i],"--force-use-unestimated-values") )
    {
      this->forceUseNonEstimatedCovarKeys = true;
    }
    
    //////////////////////////////////////////////////
    // Regional analysis
    
    //In position regional analysis, which is the region size?
    else if( !strcmp(argv[i],"--region-size") )
    {
      this->regionSize = getInt(argc, argv, i + 1, Range(1, false));
      this->regionSize *= 1000;
      this->regionBy = byPosition;
      this->regionalAnalysis = true;
    }
    //In position regional analysis, which is the region overlapping?
    else if( !strcmp(argv[i],"--region-overlap") )
    {
      this->regionOverlapping = getInt(argc, argv, i + 1, Range(0, false));
      this->regionOverlapping *= 1000;
      this->regionBy = byPosition;
      this->regionalAnalysis = true;
    }
    //File assigning groups to SNPs
    else if( !strcmp(argv[i],"--groups") )
    {
      this->regionsFile = getFileName(argc, argv, i + 1);
      this->regionBy = byGene;
      this->regionalAnalysis = true;
    }
    //Group with all SNPs in a file
    else if( !strcmp(argv[i],"--group-all") )
    {
      appendParsedOptions(std::string(argv[i]), "");
      this->regionBy = byAll;
      this->regionalAnalysis = true;
    }
    //Which is the minimum region size?
    else if( !strcmp(argv[i],"--min-snps-region") )
    {
      this->minSNPsInRegions = getInt(argc, argv, i + 1, Range(1, false));
    }
    //All regions together?
    else if( !strcmp(argv[i],"--all-together") )
    {
      this->allRegionsTogether = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    
    ////////////////////////////
    // Filter Labeled Matrix
    else if( !strcmp(argv[i],"--imatrix") )
    {
      this->inputLabeledMatrix = getString(argc, argv, i + 1);
      misc.checkFileExists(this->inputLabeledMatrix + ".colids");
      misc.checkFileExists(this->inputLabeledMatrix + ".rowids");
      misc.checkFileExists(this->inputLabeledMatrix + ".dat");
    }
    else if( !strcmp(argv[i],"--row-labels") )
    {
      this->rowLabelsFile = getFileName(argc, argv, i + 1);;
    }
    else if( !strcmp(argv[i],"--col-labels") )
    {
      this->columnLabelsFile = getFileName(argc, argv, i + 1);;
    }
    
    //////////////////////////////////////////////////
    // Group effects analysis
    else if( !strcmp(argv[i],"--keep-groups") )
    {
      this->fileGroupsToKeep = getFileName(argc, argv, i + 1);
    }
    else if( !strcmp(argv[i],"--groups-positions") )
    {
      this->groupsPositions = getFileName(argc, argv, i + 1);;
    }
    else if( !strcmp(argv[i],"--effects-files") )
    {
      this->groupEffectsFileList = getFileNameList(argc, argv, i + 1, false);
      for(int i = 0; i<this->groupEffectsFileList.size(); i++)
      {
        misc.checkFileExists(this->groupEffectsFileList[i] + ".colids");
        misc.checkFileExists(this->groupEffectsFileList[i] + ".rowids");
        misc.checkFileExists(this->groupEffectsFileList[i] + ".dat");
      }
    }
    else if( !strcmp(argv[i],"--group-min-distance") )
    {
      this->groupDistanceForDiscarding = getInt(argc, argv, i + 1, Range(0, false));
    }
    else if( !strcmp(argv[i],"--effects-pair-files") )
    {
      this->groupEffectsPairFileList = getFileNameList(argc, argv, i + 1, false);
      this->crossedEffectCorrelations = true;
      if(this->groupEffectsPairFileList.size()%2 != 0)
      {
        misc.error("Error: The files specified with the option --effects-pair-files have to be a multiple of 2.", 0);
      }
      for(int i = 0; i<this->groupEffectsPairFileList.size(); i++)
      {
        misc.checkFileExists(this->groupEffectsPairFileList[i] + ".colids");
        misc.checkFileExists(this->groupEffectsPairFileList[i] + ".rowids");
        misc.checkFileExists(this->groupEffectsPairFileList[i] + ".dat");
      }
    }
    
    //////////////////////////////////////////////////
    // Divers
    //Analyze all phenos in phenotypes file recursively?
    else if( !strcmp(argv[i],"--all-phenos") )
    {
      this->analyzeAllPhenos = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Use GCTA approaches in REML?
    else if(!strcmp(argv[i],"--gcta-mode"))
    {
      this->remlGCTAMode = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Set the default block size for distributed matrices.
    else if( !strcmp(argv[i],"--default-block-size") )
    {
      this->defaultBlockSize = getInt(argc, argv, i + 1, Range(16, 4096));
    }
    //Allow single precision inversion?
    else if( !strcmp(argv[i],"--no-single-precision") )
    {
      this->allowSinglePrecisionInversion = false;
      appendParsedOptions(std::string(argv[i]), "");
    }
    else if( !strcmp(argv[i],"--extract") )
    {
      this->fileSNPsToKeep = getFileName(argc, argv, i + 1);
    }
    else if( !strcmp(argv[i],"--keep") )
    {
      this->fileIndividualsToKeep = getFileName(argc, argv, i + 1);
    }
    //Write a more verbose log?
    else if( !strcmp(argv[i],"--verbose") )
    {
      this->verboseLog = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Allow BGEN files with layout 1
    else if( !strcmp(argv[i],"--bgen-l1") )
    {
      this->allowBGENLayout1 = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    //Keep SNPs with zero variance
    else if( !strcmp(argv[i],"--keep-zerostd-snps") )
    {
      this->keepSNPsWithZeroVariance = true;
      appendParsedOptions(std::string(argv[i]), "");
    }
    
    
    //////////////////////////////////////////////////
    // Debugging options
    
    else if( !strcmp(argv[i],"--mpi-debug") )
    {
      this->mpiDebug = true;
    }
    else if( !strcmp(argv[i],"--debug-vars") )
    {
      this->genotypeFile = "test/test3";
      this->genotypeListFile = "";
      this->grmFile = "";
      this->phenotypesFile = "test/test3.phenos";
      this->qCovarsFile = "test/test3.qcov";
      this->qCovarsFiles.clear();
      this->qCovarsFiles.push_back("test/test3.qcov");
      this->qCovarsFiles.push_back("test/test3.qcov");
      this->covarsFile = "test/test3.cov";
      this->covarsFiles.clear();
      this->covarsFiles.push_back("test/test3.cov2");
      this->covarsFiles.push_back("test/test3.cov2");
      this->effectsSizeFile = "test/test3.effects";
    }
    else if( !strcmp(argv[i],"--debug") )
    {
      this->debug = true;
    }
    else if( !strcmp(argv[i],"--debug-default-block-size") )
    {
      this->defaultBlockSize = getInt(argc, argv, i + 1, Range(1, 4096));
    }
    
    //////////////////////////////////////////////////
    // Invalid flag?
    else
    {
      if(!strncmp(argv[i],"--",2) || !strncmp(argv[i],"-",1))
      {
        std::string temp(argv[i]);
        std::stringstream ss(temp);
        double isDoubleTest;
        if( (ss >> isDoubleTest).fail() )
        {
          //appendParsedOptions(std::string(argv[i]), " is not a valid option. It will be ignored.");
          misc.error("Error: " + std::string(argv[i]) + " is not a valid option.", 0);
        }
      }
    }
  }
  
  fillMissingOptions();
  this->restrictedCovariances = getRestrictedCovariancesList(tempRestrictedCovariances);
}

void Options::checkIncompatibilities()
{
  std::map< std::string, std::vector<std::string> > requiredOptions;
  std::map< std::string, std::vector<std::string> > incompatibleOptions;
  

  //////////////////////////////////////////////////
  // Base input files
  if( this->genotypeListFile != "" && this->genotypeFile != "" )
  {
    misc.error("Error: --bfile and --bfile-list cannot be used at same time.", 0);
  }
  if( (this->genotypeListFile != "" && this->SNPBLUPGenotypeFiles.size() != 0 ) || (this->genotypeFile != "" && this->SNPBLUPGenotypeFiles.size() != 0) )
  {
    misc.error("Error: --bfile or --bfile-list cannot be used together with --blup-bfile-list.", 0);
  }
  if( this->grmFile != "" && this->grmListFile != "" )
  {
    misc.error("Error: --grm-list and --grm cannot be used together.", 0);
  }
  if( (this->genotypeFile != "" && this->grmFile != "") && this->regionalAnalysis == false && this->computeSNPsBLUP == false)
  {
    misc.error("Error: When running a non-regional REML, --bfile and --grm can only be used at the same time in a reml regional analysis or in conjunction with option --snp-blup.", 0);
  }
  if( (this->genotypeListFile != "" || this->genotypeFile != "") && (this->grmListFile != "" && this->addAllGRMsInOne == false) )
  {
    misc.error("Error: --grm-list can only be used with --bfile or --bfile-list when combined with --add-grms.", 0);
  }
  if( this->SNPBLUPGenotypeFiles.size() != 0 && this->grmListFile != "" && this->addAllGRMsInOne == false )
  {
    misc.error("Error: --grm-list can only be used with --blup-bfile-list when combined with --add-grms.", 0);
  }
  if( this->computeSNPsBLUP == false && this->SNPBLUPGenotypeFiles.size() != 0 )
  {
    misc.error("Error: --blup-bfile-list can only be used in conjunction with option --snp-blup.", 0);
  }
  if( this->phenotypesFile == "" && ( this->analysis == bivarREMLAnalysis || this->analysis == multivarREMLAnalysis || this->analysis == GWASAnalysis || this->analysis == recursiveGWASAnalysis) )
  {
    misc.error("Error: A phenotypes file must be specified for running REML or GWAS analysis.", 0);
  }
  if( (this->phenotypesFile == "" && this->phenotypesFiles.size() == 0) && (this->analysis == REMLAnalysis) )
  {
    misc.error("Error: A phenotypes file must be specified for running REML.", 0);
  }
  if( this->parsedOptions.count("--covar") != 0 && this->parsedOptions.count("--covars") != 0 )
  {
    misc.error("Error: --covar and --covars cannot be used at same time.", 0);
  }
  if( this->parsedOptions.count("--qcovar") != 0 && this->parsedOptions.count("--qcovars") != 0 )
  {
    misc.error("Error: --qcovar and --qcovars cannot be used at same time.", 0);
  }
  if( this->parsedOptions.count("--initial-h2") != 0 && this->parsedOptions.count("--initial-h2s") != 0 )
  {
    misc.error("Error: --initial-h2 and --initial-h2s cannot be used at same time.", 0);
  }
  if( this->parsedOptions.count("--pheno-col") != 0 && this->parsedOptions.count("--pheno-cols") != 0 )
  {
    misc.error("Error: --pheno-col and --pheno-cols cannot be used at same time.", 0);
  }
  if( this->parsedOptions.count("--pheno") != 0 && this->parsedOptions.count("--phenos") != 0 )
  {
    misc.error("Error: --pheno and --phenos cannot be used at same time.", 0);
  }
  if( this->genotypeListFile == "" && this->genotypeFile == "" && this->fileSNPsToKeep != "" && this->genotypeAndGRMsListFile == "" && this->bgenGenotypeFile == "")
  {
    misc.error("Error: --extract can only be used in combination with --bfile, --bfile-list, --bgen, or --bfile-grm-list.", 0);
  }
  if( this->analysis != REMLAnalysis && this->analysis != bivarREMLAnalysis && this->analysis != multivarREMLAnalysis && this->analysis != GWASAnalysis &&
    this->analysis != recursiveGWASAnalysis && this->analysis != predictPhenotypeAnalysis && this->analysis != simulatePhenotypeAnalysis &&
    this->analysis != iGWASAnalysis  &&  this->analysis != predictCovarPhenotypeAnalysis && this->analysis != multiplePhenotypeResiduals
  )
  {
    if( this->parsedOptions.count("--covar") != 0 )
    {
      misc.error("Error: --covar can only be used with --reml, --bivar-reml, --multi-reml, --gwas, --predict, and --rgwas analyses.", 0);
    }
    if( this->parsedOptions.count("--qcovar") != 0 )
    {
      misc.error("Error: --qcovar can only be used with --reml, --bivar-reml, --multi-reml, --gwas, --predict, and --rgwas analyses.", 0);
    }
  }
  if( this->analysis != REMLAnalysis && this->analysis != GWASAnalysis && this->analysis != recursiveGWASAnalysis )
  {
    if( this->parsedOptions.count("--pheno-col") != 0 )
    {
      misc.error("Error: --pheno-col can only be used with --reml, --gwas, and --rgwas analyses.", 0);
    }
  }
  if( this->analysis != REMLAnalysis && this->analysis != bivarREMLAnalysis && this->analysis != multivarREMLAnalysis && this->analysis != GWASAnalysis  && this->analysis != iGWASAnalysis)
  {
    if( this->parsedOptions.count("--random-effects") != 0 || this->parsedOptions.count("--random-effects-cols") != 0 )
    {
      misc.error("Error: --random-effects and --random-effects-cols options can only be used with --reml, --bivar-reml, --multi-reml, or --gwas analyses.", 0);
    }
    if( this->parsedOptions.count("--initial-h2") != 0 )
    {
      misc.error("Error: --initial-h2 option can only be used with --reml, --bivar-reml, --multi-reml, or --gwas analyses.", 0);
    }
  }
  
  //////////////////////////////////////////////////
  // GRM
  if( this->parsedOptions.count("--grm-join-method") == 1 && this->genotypeListFile == "" )
  {
    misc.error("Error: --grm-join-method must be used with --bfile-list option.", 0);
  }
  if( this->parsedOptions.count("--store-both") != 0 && this->parsedOptions.count("--diagonalize") == 0)
  {
    misc.error("Error: --store-both option can only be used with --diagonalize option.", 0);
  }
  if( this->analysis != makeGRMAnalysis )
  {
    if( this->parsedOptions.count("--diagonalize") != 0 )
    {
      misc.error("Error: --diagonalize option can only be used with --make-grm option.", 0);
    }
  }

  //////////////////////////////////////////////////
  // REML
  if( ( this->analysis == bivarREMLAnalysis || this->analysis == multivarREMLAnalysis) && (this->genotypeFile == "" && this->genotypeListFile == "" && this->grmFile == "" && this->grmListFile == "") )
  {
    misc.error("Error: For running a REML analysis, at least a genotype file(s) or a grm file(s) must be specified.", 0);
  }
  if( (this->analysis == REMLAnalysis && this->regionalAnalysis == true) && (this->genotypeFile == "" && this->genotypeListFile == "") )
  {
    misc.error("Error: For running a REML regional analysis, a genotype file(s) must be specified.", 0);
  }
  if( (this->analysis == REMLAnalysis && this->regionalAnalysis == false) && (this->genotypeFile == "" && this->genotypeListFile == "" && this->grmFile == "" && this->grmListFile == "" && this->discreteRandomEffectsFile == "" && this->multiDiscreteRandomEffectsFile == "" && this->squaredExponentialKernelFiles.size() == 0 && this->GCTAGRMsFile == "") )
  {
    misc.error("Error: For running a REML analysis, at least a genotype file(s), grm file(s), or other files for defining at least one covariance matrix, must be specified.", 0);
  }
  if( this->computeSNPsBLUP && this->regionalAnalysis && (this->analysis != REMLAnalysis || this->allRegionsTogether == false) )
  {
    misc.error("Error: --snp-blup option cannot be used in a regional analysis except when all regions estimated together in a single REML.", 0);
  }
  if( this->computeSNPsBLUP && this->parsedOptions.count("--bfile") == 0 && this->parsedOptions.count("--bfile-list") == 0 && this->SNPBLUPGenotypeFiles.size() == 0 && this->grmListFile == "" )
  {
    misc.error("Error: ---bfile, --bfile-list, --blup-bfile-list, or --grm-list must be used when using --snp-blup option.", 0);
  }
  if( this->parsedOptions.count("--reml-method-em") != 0 && this->parsedOptions.count("--reml-method-ai") != 0 )
  {
    misc.error("Error: --reml-method-em and --reml-method-ai cannot be used at the same time.", 0);
  }
  if( this->parsedOptions.count("--reml-subsample") == 0 && this->parsedOptions.count("--subsample-replicates") != 0 )
  {
    misc.error("Error: Options --subsample-replicates must be used with --reml-subsample.", 0);
  }
  if( this->computeEpistasisVariance == true && isGRMDiagonal(this->grmFile) )
  {
    misc.error("Error: --epistasis-var option cannot be used with diagonal GRMs.", 0);
  }
  if( this->parsedOptions.count("--random-effects") != this->parsedOptions.count("--random-effects-cols") )
  {
    misc.error("Error: --random-effects and --random-effects-cols options have to be used together.", 0);
  }
  if( this->parsedOptions.count("--skip-test-reduced-models") != 0 && this->parsedOptions.count("--fix-correlation") )
  {
    misc.error("Error: --skip-test-reduced-models and --fix-correlation options cannot be used together.", 0);
  }
  if( this->useLogLogisticScale == true && this->useCorrelations == false && this->analysis != REMLAnalysis )
  {
    misc.error("Error: --use-log-logistic has to be used with --use-correlations option when performing a bivariate or multivariate REML analysis.", 0);
  }
  if( this->indirectEffectsCouplesFile != "" && this->phenotypeColumns.size() > 2 )
  {
    misc.error("Error: --indirect-effects-couples option can only be used with one or two phenotypes.", 0);
  }
  if( this->parsedOptions.count("--weights") == 0 && ( this->parsedOptions.count("--weights-col") != 0 || this->parsedOptions.count("--no-scale-weights") != 0 ) )
  {
    misc.error("Error: --weights-col or --no-scale-weights options can only be used with --weights option.", 0);
  }
  if( this->parsedOptions.count("--indiv-blup") == 0 && this->parsedOptions.count("--snp-blup") == 0 && this->parsedOptions.count("--indiv-blup-error") != 0)
  {
    misc.error("Error: --indiv-blup-error options can only be used with --indiv-blup and/or --indiv-snp options.", 0);
  }
  if( this->analysis != REMLAnalysis && this->analysis != bivarREMLAnalysis && this->analysis != multivarREMLAnalysis )
  {
    if( this->computeBLUE || this->computeIndividualsBLUP )
    {
      misc.error("Error: options --blue, or --indiv-blup can only be used with --reml, --bivar-reml, or --multi-reml analyses.", 0);
    }
    if( this->computeSNPsBLUP )
    {
      misc.error("Error: --snp-blup option can only be used with --reml, --bivar-reml, or --multi-reml analyses.", 0);
    }
    if( this->parsedOptions.count("--initial-variances") != 0 )
    {
      misc.error("Error: --initial-variances option can only be used with single --reml analysis.", 0);
    }
    if( this->parsedOptions.count("--skip-test-reduced-models") != 0 )
    {
      misc.error("Error: --skip-test-reduced-models option can only be used with --reml, --bivar-reml or --multi-reml analyses.", 0);
    }
    if( this->parsedOptions.count("--use-correlations") != 0 )
    {
      misc.error("Error: --use-correlations option can only be used with --reml, --bivar-reml or --multi-reml analyses.", 0);
    }
    if( this->parsedOptions.count("--use-log-logistic") != 0 )
    {
      misc.error("Error: --use-log-logistic option can only be used with --reml, --bivar-reml or --multi-reml analyses.", 0);
    }
    if( this->allRegionsTogether == true )
    {
      misc.error("Error: --all-together option can only be used with --reml, --bivar-reml or --multi-reml analyses.", 0);
    }
    if( this->indirectEffectsCouplesFile != "" )
    {
      misc.error("Error: --indirect-effects-couples option can only be used with --reml, --bivar-reml or --multi-reml analyses.", 0);
    }
    if(this->useLogLogisticScale == true)
    {
      misc.error("Error: --use-log-logistic has to be used with --reml, --bivar-reml or --multi-reml analyses.", 0);
    }
    if( this->parsedOptions.count("--pheno-cols") != 0 )
    {
      misc.error("Error: --pheno-cols option can only be used with --reml, --bivar-reml or --multi-reml analyses.", 0);
    }
  }
  if( this->analysis != bivarREMLAnalysis && this->analysis != multivarREMLAnalysis )
  {
    if( this->parsedOptions.count("--covars") != 0 )
    {
      misc.error("Error: --covars option can only be used with --bivar-reml or --multi-reml analyses.", 0);
    }
    if( this->parsedOptions.count("--qcovars") != 0 )
    {
      misc.error("Error: --qcovars option can only be used with --bivar-reml or --multi-reml analyses.", 0);
    }
    if( this->parsedOptions.count("--initial-h2s") != 0 )
    {
      misc.error("Error: --initial-h2s option can only be used with --bivar-reml or --multi-reml analyses.", 0);
    }
    if( this->parsedOptions.count("--no-environment-cov") != 0 )
    {
      misc.error("Error: --no-environment-cov option can only be used with  --bivar-reml or --multi-reml analyses.", 0);
    }
  }
  if( this->analysis != REMLAnalysis )
  {
    if( this->parsedOptions.count("--correct-ld") != 0 )
    {
      misc.error("Error: --correct-ld option can only be used when single REML analysis is performed.", 0);
    }
    if( this->computeEpistasisVariance == true )
    {
      misc.error("Error: --epistasis-var option can only be used when single REML analysis is performed.", 0);
    }
    if( this->useMLinsteadOfREML == true )
    {
      misc.error("Error: --use-ml option can only be used when single REML analysis is performed.", 0);
    }
    if( this->parsedOptions.count("--weights") != 0 )
    {
      misc.error("Error: --weights option can only be used when single REML analysis is performed.", 0);
    }
    if( this->parsedOptions.count("--phenos") != 0 )
    {
      misc.error("Error: --phenos option can only be used when single REML analysis is performed.", 0);
    }
  }
  
  //////////////////////////////////////////////////
  // GWAS
  if( (this->genotypeListFile != "" || this->genotypeFile != "") && (this->genotypeAndGRMsListFile != "" || this->genotypeBGENAndGRMsListFile != "" ))
  {
    misc.error("Error: --bfile-grm-list or --bgen-grm-list options cannot be used together with --bfile or --bfile-list options.", 0);
  }
  if( (this->genotypeAndGRMsListFile != "" && this->genotypeBGENAndGRMsListFile != "" ))
  {
    misc.error("Error: --bfile-grm-list and --bgen-grm-list options cannot be used together.", 0);
  }
  if( (this->genotypeAndResidualsListFile != "" && this->genotypeBGENAndResidualsListFile != "" ))
  {
    misc.error("Error: --bfile-residuals-list and --bgen-residuals-list options cannot be used together.", 0);
  }
  if( (this->genotypeAndResidualsListFile != "" || this->genotypeBGENAndResidualsListFile != "" )  && (this->genotypeAndGRMsListFile != "" || this->genotypeBGENAndGRMsListFile != "" ) )
  {
    misc.error("Error: --bfile-residuals-list or --bgen-residuals-list options cannot be used together with --bfile-grm-list and --bgen-grm-list.", 0);
  }
  if( this->parsedOptions.count("--gwas") != 0 )
  {
    if( this->genotypeListFile == "" && this->genotypeFile == "" && this->genotypeAndGRMsListFile == "" )
    {
      misc.error("Error: For a GWAS analysis one of the following options must be specified: --bfile, --bfile-list, or --bfile-grm-list.", 0);
    }
    if( this->parsedOptions.count("--group-effects") != 0  && this->regionBy == ungrouped )
    {
      misc.error("Error: --group-effects can only be used in a grouped GWAS analysis.", 0);
    }
  }
  if( this->parsedOptions.count("--gwas") == 0 )
  {
    if( this->parsedOptions.count("--group-all") != 0 )
    {
      misc.error("Error: --group-all can only be used in a GWAS analysis.", 0);
    }
    if( this->parsedOptions.count("--group-var") != 0 )
    {
      misc.error("Error: --group-var can only be used in a grouped GWAS analysis.", 0);
    }
    if( this->parsedOptions.count("--group-effects") != 0 )
    {
      misc.error("Error: --group-effects can only be used in a grouped GWAS analysis.", 0);
    }
  }
  if( this->parsedOptions.count("--mpgwas") == 0 )
  {
    if( this->genotypeBGENAndGRMsListFile != "" )
    {
      misc.error("Error: --bgen-grm-list options can only be used with --mpgwas analysis.", 0);
    }
  }
  
  //////////////////////////////////////////////////
  // Recursive GWAS
  if( this->parsedOptions.count("--snp-corr-threshold") != 0 && this->parsedOptions.count("--rgwas") == 0 )
  {
    misc.error("Error: --snp-corr-threshold can only be used in a Recursive GWAS analysis.", 0);
  }
  
  //////////////////////////////////////////////////
  // Multiple phenotype GWAS
  if( this->analysis == multiplePhenotypeGWASAnalysis )
  {
    if( this->genotypeAndResidualsListFile == "" && this->genotypeBGENAndResidualsListFile == "" )
    {
      misc.error("Error: --bfile-residuals-list or --bgen-residuals-list options have to be specified for --mpgwas analysis.", 0);
    }
  }
  
  //////////////////////////////////////////////////
  // Multiple phenotype Residuals
  if( this->analysis == multiplePhenotypeResiduals )
  {
    if( this->grmFile == "" )
    {
      misc.error("Error: --grm option has to be specified for --mpresiduals analysis.", 0);
    }
    if( this->phenotypesFile == "" )
    {
      misc.error("Error: --pheno option has to be specified for --mpresiduals analysis.", 0);
    }
  }
  
  //////////////////////////////////////////////////
  // Phenotype simulation vars 
  if(this->simulateBinaryTrait && this->simulateQuantitativeTrait)
  {
    misc.error("Error: flags --simu-quantitative and --simu-binary can not be used together.", 0);
  }
  if( this->parsedOptions.count("--prevalence") != 0 && this->simulateQuantitativeTrait == true )
  {
    misc.error("Error: flags --prevalence option cannot be used for simulating quantitative traits.", 0);
  }
  if( this->analysis == simulatePhenotypeAnalysis )
  {
    if(!this->simulateBinaryTrait && !this->simulateQuantitativeTrait )
    {
      misc.error("Error: For simulating phenotypes, one of these two options must be specified: --simu-quantitative or --simu-binary.", 0);
    }
    if( this->parsedOptions.count("--effect-sizes") == 0 )
    {
      misc.error("Error: For simulating phenotypes, an SNP effect sizes file must be specified with the option --effect-sizes.", 0);
    }
  }
  
  //////////////////////////////////////////////////
  // Phenotype prediction vars 
  if( this->parsedOptions.count("--predict") != 0 )
  {
    if( this->parsedOptions.count("--snp-effects") == 0 )
    {
      misc.error("Error: --snp-effects option must be specified when performing --predict analysis.", 0);
    }
    if( this->parsedOptions.count("--bfile") == 0 && this->parsedOptions.count("--bfile-list") == 0 )
    {
      misc.error("Error: Genotypes file(s) must be specified with --bfile or --bfile-list options when performing --predict analysis.", 0);
    }
  }
  
  if( this->parsedOptions.count("--predict") == 0  && this->parsedOptions.count("--accuracy-by-snp") == 0)
  {
    if( this->parsedOptions.count("--snp-effects") != 0 )
    {
      misc.error("Error: --snp-effects option must be used with --predict analysis.", 0);
    }
  }
  
  
  //////////////////////////////////////////////////
  // Regional analysis
  if( this->regionalAnalysis == true && !(this->analysis == REMLAnalysis || this->analysis == bivarREMLAnalysis || this->analysis == multivarREMLAnalysis || this->analysis == GWASAnalysis || this->analysis == makeGRMAnalysis) )
  {
    misc.error("Error: options for regional analysis can only be used with --reml, --bivar-reml, --multi-reml or --gwas.", 0);
  }
  if( (this->parsedOptions.count("--region-overlap") != 0 || this->parsedOptions.count("--region-size") != 0) && this->parsedOptions.count("--groups") != 0 ) //Only one type of regional analysis is allowed.
  {
    misc.error("Error: Only one type of regional analysis is allowed. The options --region-overlap and --region-size or --groups can not be used together.", 0);
  }
  if( this->regionSize < 0 && this->regionOverlapping > 0 )
  {
    misc.error("Error: The option --region-overlap can only be used together with --region-size.", 0);
  }
  if(this->regionalAnalysis && this->genotypeFile == "" && this->genotypeListFile == "" && this->genotypeAndGRMsListFile == "")
  {
    misc.error("Error: To perform a regional analysis, a genotype file must be specified.", 0);
  }
  if( this->allRegionsTogether == true && this->regionalAnalysis == false )
  {
    misc.error("Error: --all-together option can only be used with some of the regional options (e.g. --groups).", 0);
  }
  if( this->allRegionsTogether == true && this->grmFile != "" )
  {
    misc.error("Error: a grm cannot be specified with --all-together option.", 0);
  }
  
  ////////////////////////////
  // Filter Labeled Matrix
  if( this->analysis == filterLabeledMatrix )
  {
    if( this->inputLabeledMatrix == "" || this->rowLabelsFile == "" || this->columnLabelsFile == "" )
    {
      misc.error("Error: --imatrix, --row-labels, and --col-labels options have to be specified with --filter-matrix option.", 0);
    }
  }
  
  ////////////////////////////
  // Compute SNP stats
  if( ( this->analysis != computeSNPStatsAnalysis || this->bgenGenotypeFile == "" ) && this->keepSNPsWithZeroVariance == true )
  {
    misc.error("Error: --keep-zerostd-snps option can only be used with --snp-stats analysis and --bgen option.", 0);
  }
  if( this->analysis != computeSNPStatsAnalysis )
  {
    if( this->bgenGenotypeFile != "" )
    {
      misc.error("Error: --bgen option can only be used with --snp-stats analysis.", 0);
    }
  }
  
  //////////////////////////////////////////////////
  // Shared between analysis
  if( this->parsedOptions.count("--grm") != 0 && this->analysis != PCAAnalysis && this->analysis != REMLAnalysis && this->analysis != bivarREMLAnalysis && this->analysis != multivarREMLAnalysis && this->analysis != makeGRMAnalysis && this->analysis != GLMMAnalysis &&  this->analysis != multiplePhenotypeResiduals )
  {
    misc.error("Error: --grm option is only valid with --reml, --bivar-reml, --multi-reml, --make-grm, or --pca analyses.", 0);
  }
  if( this->parsedOptions.count("--predict") == 0 && this->parsedOptions.count("--simulate") == 0  && this->parsedOptions.count("--cov-predict") == 0)
  {
    if( this->parsedOptions.count("--covar-effects") != 0 || this->parsedOptions.count("--qcovar-effects") != 0 )
    {
      misc.error("Error: --covar-effects or --qcovar-effects options can only be used when performing --predict or --cov-predict analysis.", 0);
    }
  }
  else
  {
    if( this->parsedOptions.count("--covar-effects") != this->parsedOptions.count("--covar") )
    {
      misc.error("Error: --covar-effects and --covar options must be used together.", 0);
    }
    if( this->parsedOptions.count("--qcovar-effects") != this->parsedOptions.count("--qcovar") )
    {
      misc.error("Error: --qcovar-effects and --qcovar options must be used together.", 0);
    }
  }
  
  //////////////////////////////////////////////////
  // Others
  if( this->parsedOptions.count("--all-phenos") != 0 && ( this->analysis != REMLAnalysis || this->regionalAnalysis != false) )
  {
    misc.error("Error: --all-phenos option can only be used in nonregional --reml analysis.", 0);
  }
  for(int idx = 0; idx<this->SNPBLUPGenotypeFiles.size(); idx++)
  {
    misc.checkFileExists(this->SNPBLUPGenotypeFiles[idx] + ".bed");
    misc.checkFileExists(this->SNPBLUPGenotypeFiles[idx] + ".fam");
    misc.checkFileExists(this->SNPBLUPGenotypeFiles[idx] + ".bim");
  }
  if( this->parsedOptions.count("--out") != 0 && this->parsedOptions.count("--zout") != 0 )
  {
    misc.error("Error: --out and --zout options cannot be used together.", 0);
  }
  if(this->genotypeFile != "" && this->bgenGenotypeFile != "")
  {
    misc.error("Error: --bfile and --bgen options cannot be used together.", 0);
  }
  
  //////////////////////////////////////////////////
  // WARNINGS
  if( this->minimumNumberOverlapingSNPs == 0 )
  {
    this->warnings += "WARNING: Allowing individual relationships between individuals that share 0 non-missing SNPs may lead to unexpected results.\n";
  }
  
}

void Options::fillMissingOptions()
{
  //////////////////////////////////////////////////
  // Bivar/Multivar REML
  if( this->analysis == REMLAnalysis )
  {
    if( this->phenotypeColumns.size() == 0 )
    {
      this->phenotypeColumns.push_back(this->phenotypeColumn);
    }
    if( this->phenotypesFiles.size() == 0 )
    {
      this->phenotypesFiles.push_back(this->phenotypesFile);
    }
    if( this->covarsFiles.size() == 0 )
    {
      this->covarsFiles.push_back(this->covarsFile);
    }
    if( this->qCovarsFiles.size() == 0 )
    {
      this->qCovarsFiles.push_back(this->qCovarsFile);
    }
    if( this->initialh2Traits.size() == 0 )
    {
      this->initialh2Traits.push_back( this->initialh2Trait );
    }
  }
  
  //////////////////////////////////////////////////
  // Bivar/Multivar REML
  if( this->analysis == bivarREMLAnalysis || this->analysis == multivarREMLAnalysis )
  {
    if( this->useCorrelations == true )
    {
      this->allowSwitchFromAItoEM = false;
    }
  
    int nMaxPhenoCol = getNumberOfFileColumns(this->phenotypesFile) - 2;
    
    if( this->analysis == bivarREMLAnalysis )
    {
      if( this->phenotypeColumns.size() == 0 )
      {
        this->phenotypeColumns.push_back(1);
        this->phenotypeColumns.push_back(2);
      }
      else
      {
        if( this->phenotypeColumns.size() != 2 )
        {
          misc.error("Error: --pheno-cols have to specify two phenotype columns when performing bivariate reml analyses.", 0);
        }
      }
    }
    if( this->analysis == multivarREMLAnalysis )
    {
      if( this->phenotypeColumns.size() == 0 )
      {
        for(int i = 0; i<nMaxPhenoCol; i++)
        {
          this->phenotypeColumns.push_back(i+1);
        }
      }
    }

    std::set<int> repeatedPhenotypesColumnsCheck;
    for(int i = 0; i<this->phenotypeColumns.size(); i++)
    {
      if( this->phenotypeColumns[i] > nMaxPhenoCol )
      {
        misc.error("Error: Phenotypes file does not have the phenotype column " + i2s(this->phenotypeColumns[i]) + " it only contains information for " + i2s(nMaxPhenoCol) + " phenotypes.", 0);
      }
      if( repeatedPhenotypesColumnsCheck.find(this->phenotypeColumns[i]) != repeatedPhenotypesColumnsCheck.end() )
      {
        misc.error("Error: There cannot be repeated phenotype columns specified using --pheno-cols option.", 0);
      }
      repeatedPhenotypesColumnsCheck.insert(this->phenotypeColumns[i]);
    }
    if( this->phenotypeColumn > nMaxPhenoCol )
    {
      misc.error("Error: Phenotypes file does not have the phenotype column " + i2s(this->phenotypeColumn) + " it only contains information for " + i2s(nMaxPhenoCol) + " phenotypes.", 0);
    }
    
    if( this->analysis == bivarREMLAnalysis || this->analysis == multivarREMLAnalysis )
    {
      if( this->covarsFiles.size() != 0 )
      {
        if( this->covarsFiles.size() == 1 )
        {
          for(int i = 1; i<this->phenotypeColumns.size(); i++)
          {
            this->covarsFiles.push_back( this->covarsFiles[0] );
          }
        }
        else
        {
          if(this->covarsFiles.size() != this->phenotypeColumns.size())
          {
            misc.error("Error: Invalid number of files specified by --covars option. It have to be one file or a number equal to the number of phenotypes to analyze.", 0);
          }
        }
      }
      else
      {
        for(int i = 0; i<this->phenotypeColumns.size(); i++)
        {
          this->covarsFiles.push_back( this->covarsFile );
        }
      }
      
      if( this->qCovarsFiles.size() != 0 )
      {
        if( this->qCovarsFiles.size() == 1 )
        {
          for(int i = 1; i<this->phenotypeColumns.size(); i++)
          {
            this->qCovarsFiles.push_back( this->qCovarsFiles[0] );
          }
        }
        else
        {
          if(this->qCovarsFiles.size() != this->phenotypeColumns.size())
          {
            misc.error("Error: Invalid number of files specified by --qcovars option. It have to be one file or a number equal to the number of phenotypes to analyze.", 0);
          }
        }
      }
      else
      {
        for(int i = 0; i<this->phenotypeColumns.size(); i++)
        {
          this->qCovarsFiles.push_back( this->qCovarsFile );
        }
      }
      
      if( this->initialh2Traits.size() != 0 )
      {
        if( this->initialh2Traits.size() == 1 )
        {
          for(int i = 1; i<this->phenotypeColumns.size(); i++)
          {
            this->initialh2Traits.push_back( this->initialh2Traits[0] );
          }
        }
        else
        {
          if(this->initialh2Traits.size() != this->phenotypeColumns.size())
          {
            misc.error("Error: Invalid number of heritabilities specified by --initial-h2s option. It have to be one heritability or a number equal to the number of phenotypes to analyze.", 0);
          }
        }
      }
      else
      {
        for(int i = 0; i<this->phenotypeColumns.size(); i++)
        {
          this->initialh2Traits.push_back( this->initialh2Trait );
        }
      }
    }
  }//End bivar/Multivar analysis check if
  
  for(int i = 0; i<this->discreteRandomEffectsFileColumns.size(); i++)
  {
    this->discreteRandomEffectsFileColumns[i] = this->discreteRandomEffectsFileColumns[i] - 1;
  }
  for(int i = 0; i<this->multiDiscreteRandomEffectsFileColumns.size(); i++)
  {
    this->multiDiscreteRandomEffectsFileColumns[i] = this->multiDiscreteRandomEffectsFileColumns[i] - 1;
  }
  
  //Load individual set keys to keep
  if( this->fileIndividualsToKeep != "" )
  {
    this->keysIndividualsToKeep.clear();
    std::vector< std::vector<std::string> > tempTable;
    getTableFromFile(this->fileIndividualsToKeep, tempTable, 2);
    for(int i = 0; i < tempTable.size(); i++)
    {
      this->keysIndividualsToKeep.insert(tempTable[i][0] + "@" + tempTable[i][1]);
    }
    if( this->keysIndividualsToKeep.size() == 0 )
    {
      misc.message << "WARNING: There are no valid individuals specified in the file [ " << this->fileIndividualsToKeep << " ]. It will be ignored." << std::endl;
      this->fileIndividualsToKeep = "";
    }
  }
}

void Options::setAnalysis(AnalysisToPerform analysisToSet)
{
  if(this->analysis == withoutAnalysis)
  {
    this->analysis = analysisToSet;
  }
  else
  {
    misc.error("Error: Two analysis types cannot be performed at same time. Please, specify only one analysis.", 0);
  }
}

std::string Options::getString(int argc, char **argv, int i)
{
  std::string result;
  std::string currentFlag;
  
  if(i>0)
  {
    currentFlag = std::string(argv[i - 1]);
  }
  else
  {
    misc.error("Error: An internal error was happened while parsing the options.", 0);
  }
  
  if (i < argc)
  {
    if(!strncmp(argv[i],"--",2) || !strncmp(argv[i],"-",1))
    {
      misc.error("Error: " + std::string(argv[i]) + " is not a valid parameter for " + currentFlag + ".", 0);
    }
    result = std::string(argv[i]);
    if(result.size() == 0)
    {
      misc.error("Error: there is not a valid parameter for " + currentFlag + ".", 0);
    }
  }
  else
  {
    misc.error("Error: The option " + currentFlag + " expects a parameter.", 0);
  }
  
  appendParsedOptions(currentFlag, std::string(argv[i]));
  
  return result;
}

std::string Options::getFileName(int argc, char **argv, int i)
{
  std::string fname = getString(argc, argv, i);
  misc.checkFileExists(fname);
  return fname;
}

int Options::getInt(int argc, char **argv, int i, Range range)
{
  int result;
  std::string currentFlag;
  
  if(i>0)
  {
    currentFlag = std::string(argv[i - 1]);
  }
  else
  {
    misc.error("Error: An internal error was happened while parsing the options.", 0);
  }
  
  if (i < argc)
  {
    std::stringstream sstemp;
    sstemp << std::string(argv[i]);
    if( (sstemp >> result).fail() )
    {
      misc.error("Error: The option " + currentFlag + " expects a valid integer value.", 0);
    }
  }
  else
  {
    misc.error("Error: The option " + currentFlag + " expects a valid integer value.", 0);
  }
  
  if(!range.checkRange(result))
  {
    misc.error("Error: The value of " + currentFlag + " is not in the proper range: " + range.explainRange(result), 0);
  }
  
  appendParsedOptions(currentFlag, std::string(argv[i]));
  
  return result;
}

double Options::getDouble(int argc, char **argv, int i, Range range)
{
  double result;
  std::string currentFlag;
  
  if(i>0)
  {
    currentFlag = std::string(argv[i - 1]);
  }
  else
  {
    misc.error("Error: An internal error was happened while parsing the options.", 0);
  }
  
  if (i < argc)
  {
    std::stringstream sstemp;
    sstemp << std::string(argv[i]);
    if( (sstemp >> result).fail() )
    {
      misc.error("Error: The option " + currentFlag + " expects a number.", 0);
    }
  }
  else
  {
    misc.error("Error: The option " + currentFlag + " expects a number.", 0);
  }
  
  if(!range.checkRange(result))
  {
    misc.error("Error: The value of " + currentFlag + " is not in the proper range: " + range.explainRange(result), 0);
  }
  
  appendParsedOptions(currentFlag, std::string(argv[i]));
  
  return result;
}

std::vector<std::string> Options::getStringList(int argc, char **argv, int i)
{
  std::vector<std::string> result;
  std::string currentFlag;
  std::string parameters;
  
  if(i>0)
  {
    currentFlag = std::string(argv[i - 1]);
  }
  else
  {
    misc.error("Error: An internal error was happened while parsing the options.", 0);
  }
  
  for (int idx = i; idx < argc; idx++) 
  {
    if( !strncmp(argv[idx],"--",2) || !strncmp(argv[idx],"-",1) )
    {
      std::string temp(argv[idx]);
      std::stringstream ss(temp);
      double isDoubleTest;
      if( (ss >> isDoubleTest).fail() )
      {
        break;
      }
    }

    std::string value = std::string(argv[idx]);
    if(value.size() == 0)
    {
      misc.error("Error: there is not a valid parameter for " + currentFlag + ".", 0);
    }
    
    parameters += value + " ";
    
    result.push_back(value);
  }
  
  appendParsedOptions(currentFlag, parameters);
  
  if( result.size() == 0 )
  {
    misc.error("Error: " + currentFlag + " expects a list of parameters.", 0);
  }
  
  return result;
}

std::vector<std::string> Options::getFileNameList(int argc, char **argv, int i, bool checkExistence)
{
  std::vector<std::string> fnames = getStringList(argc, argv, i);
  if(checkExistence == true)
  {
    for(int idx = 0; idx<fnames.size(); idx++)
    {
      misc.checkFileExists(fnames[idx]);
    }
  }
  return fnames;
}

std::vector<int> Options::getIntList(int argc, char **argv, int i, Range range)
{
  std::vector<int> result;
  std::string currentFlag;
  std::string parameters;
  
  if(i>0)
  {
    currentFlag = std::string(argv[i - 1]);
  }
  else
  {
    misc.error("Error: An internal error was happened while parsing the options.", 0);
  }
  
  for (int idx = i; idx < argc; idx++) 
  {
    if(!strncmp(argv[idx],"--",2) || !strncmp(argv[idx],"-",1))
    {
      break;
    }
    
    std::stringstream sstemp;
    sstemp << std::string(argv[idx]);
    int value;
    if( (sstemp >> value).fail() )
    {
      misc.error("Error: The option " + currentFlag + " expects a valid list of integer values.", 0);
    }

    if(!range.checkRange(value))
    {
      misc.error("Error: The value of " + currentFlag + " is not in the proper range: " + range.explainRange(value), 0);
    }
    
    parameters += std::string(argv[idx]) + " ";
    
    result.push_back(value);
  }
    
  appendParsedOptions(currentFlag, parameters);
  
  if( result.size() == 0 )
  {
    misc.error("Error: " + currentFlag + " expects a list of integer values.", 0);
  }
  
  return result;
}

std::vector<double> Options::getDoubleList(int argc, char **argv, int i, Range range)
{
  std::vector<double> result;
  std::string currentFlag;
  std::string parameters;
  
  if(i>0)
  {
    currentFlag = std::string(argv[i - 1]);
  }
  else
  {
    misc.error("Error: An internal error was happened while parsing the options.", 0);
  }
  
  for (int idx = i; idx < argc; idx++) 
  {
    if(!strncmp(argv[idx],"--",2) || !strncmp(argv[idx],"-",1))
    {
      break;
    }
    
    std::stringstream sstemp;
    sstemp << std::string(argv[idx]);
    double value;
    if( (sstemp >> value).fail() )
    {
      misc.error("Error: The option " + currentFlag + " expects a valid list of numbers.", 0);
    }

    if(!range.checkRange(value))
    {
      misc.error("Error: The value of " + currentFlag + " is not in the proper range: " + range.explainRange(value), 0);
    }
    
    parameters += std::string(argv[idx]) + " ";
    
    result.push_back(value);
  }
    
  appendParsedOptions(currentFlag, parameters);
  
  if( result.size() == 0 )
  {
    misc.error("Error: " + currentFlag + " expects a list of numbers.", 0);
  }
  
  return result;
}

std::vector<REMLCovarianceRestricted> Options::getRestrictedCovariancesList( std::vector<std::string> listOfArgs)
{
  std::vector<REMLCovarianceRestricted> result;
  std::string currentFlag = "--fix-correlation";
  
  if( this->parsedOptions.count(currentFlag) == 0 )
  {
    return result;
  }
  
  std::set<std::string> validVarBaseKeys;
  std::string validVarBaseKeysStr;
  for(std::map< std::string, std::string>::iterator it = baseVarianceNames.begin(); it != baseVarianceNames.end(); ++it)
  {
    if(it->first != "cov-")
    {
      validVarBaseKeys.insert(it->second);
      validVarBaseKeysStr += " " + it->second;
    }
    else
    {
      //This is no longer required since the random effect files have header.
    }
  }
  if( this->discreteRandomEffectsFile != "" )
  {
    std::vector<std::string> header = getHeader(this->discreteRandomEffectsFile);
    for(int i =0; i < this->discreteRandomEffectsFileColumns.size(); i++)
    {
      if( (2 + this->discreteRandomEffectsFileColumns[i]) >= header.size())
      {
        misc.error("Error: The column " + i2s(this->discreteRandomEffectsFileColumns[i]) + " is larger than the number of data (non-individual FID/IID columns) columns of the file [ " + this->discreteRandomEffectsFile + " ].", 0);
      }
      std::string varName = header[ 2 + this->discreteRandomEffectsFileColumns[i] ];
      validVarBaseKeys.insert(varName);
      validVarBaseKeysStr += " " + varName;
    }
  }
  if( this->multiDiscreteRandomEffectsFile != "" )
  {
    std::vector<std::string> header = getHeader(this->multiDiscreteRandomEffectsFile);
    for(int i =0; i < this->multiDiscreteRandomEffectsFileColumns.size(); i++)
    {
      if( (2 + this->multiDiscreteRandomEffectsFileColumns[i]) >= header.size())
      {
        misc.error("Error: The column " + i2s(this->multiDiscreteRandomEffectsFileColumns[i]) + " is larger than the number of data (non-individual FID/IID columns) columns of the file [ " + this->discreteRandomEffectsFile + " ].", 0);
      }
      std::string varName = header[ 2 + this->multiDiscreteRandomEffectsFileColumns[i] ];
      validVarBaseKeys.insert(varName);
      validVarBaseKeysStr += " " + varName;
    }
  }
  validVarBaseKeys.insert("E");
  
  if( listOfArgs.size() == 0 || listOfArgs.size() % 4 != 0 )
  {
    misc.error("Error: The number of elements expected with " + currentFlag + " have to be a multiple of 4.", 0);
  }
  
  for (int idx = 0; idx < listOfArgs.size(); idx+=4) 
  {
    REMLCovarianceRestricted restricted;
    
    restricted.baseId = listOfArgs[idx];
    
    std::stringstream sstemp1;
    sstemp1 << std::string(listOfArgs[idx + 1]);
    if( (sstemp1 >> restricted.p1).fail() )
    {
      misc.error("Error: The argument " + listOfArgs[idx + 1] + " for option " + currentFlag + " is an invalid integer.", 0);
    }
    std::stringstream sstemp2;
    sstemp2 << std::string(listOfArgs[idx + 2]);
    if( (sstemp2 >> restricted.p2).fail() )
    {
      misc.error("Error: The argument " + listOfArgs[idx + 2] + " for option " + currentFlag + " is an invalid integer.", 0);
    }
    
    std::stringstream sstemp3;
    sstemp3 << std::string(listOfArgs[idx + 3]);
    if( (sstemp3 >> restricted.correlation).fail() )
    {
      misc.error("Error: The argument " + listOfArgs[idx + 3] + " for option " + currentFlag + " is an invalid number.", 0);
    }
    if( restricted.correlation < -1. || restricted.correlation > 1. )
    {
      misc.error("Error: The argument " + listOfArgs[idx + 3] + " for option " + currentFlag + " have to be between -1 and 1.", 0);
    }

    if( restricted.p1 < 1 || restricted.p1 > this->phenotypeColumns.size() || restricted.p2 < 1 || restricted.p2 > this->phenotypeColumns.size())
    {
      misc.error("Error: The list of phenotypes for option " + currentFlag + " have the be between 1 and the number of the phenotypes used in the analysis (" + i2s(this->phenotypeColumns.size()) + " in this analysis).", 0);
    }
    if( restricted.p1 == restricted.p2 )
    {
      misc.error("Error: The phenotypes pairs for option " + currentFlag + " cannot be equal.", 0);
    }
    if( restricted.p2 < restricted.p1 )
    {
      int temp = restricted.p1;
      restricted.p1 = restricted.p2;
      restricted.p2 = temp;
    }
    
    if(validVarBaseKeys.find(restricted.baseId) == validVarBaseKeys.end())
    {
      misc.error("Error: The argument " + listOfArgs[idx] + " for option " + currentFlag + " is not valid variance. Valid variances are:" + validVarBaseKeysStr, 0);
    }
      
      
    result.push_back(restricted);
  }
    
  return result;
}

void Options::appendParsedOptions(std::string option, std::string parameters)
{
  if(this->parsedOptions.count(option) != 0)
  {
    misc.error("Error: The option " + option + " is specified more than once.", 0);
  }
  this->parsedOptions[option] = parameters;
}

void Options::showParsedOptions()
{
  misc.message << "Parsed options:" << std::endl;
  for(std::map<std::string, std::string>::iterator it =  this->parsedOptions.begin(); it != this->parsedOptions.end(); ++it)
  {
    misc.message << "\t" << it->first << " " << it->second << std::endl;
  }
}

bool Options::isGRMDiagonal(std::string fn)
{
  if( fn == "" )
  {
    return false;
  }
  
  bool isDiagonal = false;
  if(communicator->mpiRoot)
  {
    //Read data
    std::ifstream file;
    std::string fname = fn + ".grm.dat";
    misc.checkFileExists(fname);
    file.open(fname.c_str(), std::ifstream::in | std::ifstream::binary);
    
    //Read header
    unsigned char *header = new unsigned char [10];
    file.read((char*)header, 10);
    if (
      header[0] != 'G' ||
      header[1] != 'R' ||
      header[2] != 'M' ||
      header[3] != '\0' ||
      header[4] != 0x5A ||
      header[5] != 0x99 ||
      header[6] != 0x2 ||
      header[7] != 0x1
    )
    {
      misc.error("Error: The file [ " + fname + " ] do not have a proper header or is not a valid file.", 0);
    }
    if( header[8] != sizeof(double) )
    {
      misc.error("Error: Data types in file [ " + fname + " ] differ with this system datatypes. Are you loading a grm computed in a different system?", 0);
    }
    if( header[9] != 0x1 && header[9] != 0x3 )
    {
      misc.error("Error: Error in file [ " + fname + " ]. Non normalized matrices still not supported.", 0);
    }
    if(header[9] == 0x3)
    {
      isDiagonal = true;
      misc.message << "Reading diagonalized GRM data from file [ " << fn << ".grm.dat/diag ] ..." << std::endl;
    }
    delete [] header;
    file.close();
  }
  communicator->broadcast(&isDiagonal);
  return isDiagonal;
}

void Options::showOptions()
{
  std::stringstream specifiedOptions;
  misc.message << "The available options are: " << std::endl << std::endl;
  misc.message << ""
  "\n"
  "    Actions\n"
  "    -------------------------------------\n"
  "\n"
  "    --check                 Just check if the options could be parsed.\n"
  "                            This could be used together with options\n"
  "                            to check.\n"
  "\n"
  "    --help, -h              Show this help\n"
  "\n"
  "    --make-grm              Compute the grm\n"
  "\n"
  "    --reml                  Perform univariate MLM analysis\n"
  "\n"
  "    --bivar-reml            Perform bivariate MLM analysis\n"
  "\n"
  "    --simulate              Simulate phenotypes\n"
  "\n"
  "    --predict               Predict phenotypes\n"
  "\n"
  "    --pca                   Perform PCA analysis\n"
  "\n"
  "    --gwas                  Perform GWAS analysis\n"
  "\n"
  "\n"
  "    Base input files\n"
  "    -------------------------------------\n"
  "\n"
  "    --bfile f               Specify the genotypes file.\n"
  "\n"
  "    --bfile-list f          Specify a file with a list of genotypes files.\n"
  "\n"
  "    --grm f                 Specify the GRM file.\n"
  "\n"
  "    --pheno f               Specify the phenotypes file.\n"
  "\n"
  "    --covar f               Specify the discrete covariance file.\n"
  "\n"
  "    --covars f1 f2          Specify the discrete covariance files\n"
  "                            for bivariate MLM, if different for both\n"
  "                            traits\n"
  "\n"
  "    --qcovar f              Specify the quantitative covariates file 1.\n"
  "\n"
  "    --qcovars f1 f2         Specify the quantitative covariates files\n"
  "                            for bivariate MLM, if different for both\n"
  "                            traits\n"
  "\n"
  "\n"
  "    Output files\n"
  "    -------------------------------------\n"
  "\n"
  "    --out                   Specify the base name for output files.\n"
  "\n"
  "\n"
  "    GRM options\n"
  "    -------------------------------------\n"
  "\n"
  "    --diagonalize           Perform the eigendecomposition of the GRM and\n"
  "                            store it.\n"
  "\n"
  "    --store-both            When using --diagonlize option, store also the.\n"
  "                            original GRM.\n"
  "\n"
  "    --grm-join-method n     Specify how grm will be constructed.\n"
  "                            n = 0 creates a GRM for each genotype file.\n"
  "                                  Then, add GRMs (default).\n"
  "                            n = 1 joins all genotype files before computing\n"
  "                                  GRM. May need much more memory.\n"
  "\n"
  "\n"
  "    MLM\n"
  "    -------------------------------------\n"
  "\n"
  "    --blue                  Compute BLUEs.\n"
  "\n"
  "    --indiv-blup            Compute BLUPs of individuals.\n"
  "\n"
  "    --snp-blup              Compute BLUPs of SNPs.\n"
  "\n"
  "    --reml-maxit n          Specify maximum REML iterations.\n"
  "\n"
  "    --variance-constrain x  The variances can be constrained.\n"
  "\n"
  //"    --fix-correlation x     Fix the correlation to x and repeat MLM.\n"
  //"\n"
  "    --pheno-col n           Specify which column use from phenotypes file.\n"
  "                            (Default n = 1)\n"
  "\n"
  "    --pheno-cols n1 n2      Specify which column use from phenotypes file\n"
  "                            for bivariate analysis.\n"
  "                            (Default n1 = 1, n2 = 2)\n"
  "\n"
  "    --initial-h2 x          Specify an initial value for h2.\n"
  "                            (Default x = 0.5)\n"
  "\n"
  "    --initial-h2s x1 x2     Specify an initial value for h2 of the traits\n"
  "                            in a MLM analysis.\n"
  "                            (Default x1 = 0.5, x2 = 0.5)\n"
  "\n"
  "    --reml-method-em        REML steps using EM method.\n"
  "\n"
  "    --reml-method-ai        REML steps using AI method.\n"
  "\n"
  "    --no-environment-cov    Exclude environment covariance in bivariate\n"
  "                            REML analysis.\n"
  "\n"
  //"    --reml-subsample        Make prior MLM computations with random\n"
  //"                            individual sampling for estimating.\n"
  //"                            initial variance values.\n"
  //"\n"
  //"    --subsample-replicates  Number of random subsample iterations\n"
  //"                            for estimating initial variances.\n"
  "\n"
  "\n"
  "    Other\n"
  "    -------------------------------------\n"
  "\n"
  "    --random-seed           Set a seed for random numbers\n"
  "                            If not specified, a random seed will be chosen.\n"
  "\n"
  "    --grm-cutoff x          Specify thresold. All GRM elements above that\n"
  "                            threshold will be eliminated.\n"
  "\n"
  //"    --min-overlap-snps n    Specify the minimum number of SNPs that have"
  //"                            to be used to compute a relationship. If lower"
  //"                            the individual pair will be removed."
  "\n"
  "\n"
  "    PCA\n"
  "    -------------------------------------\n"
  "\n"
  "    --num-eval n            Specify the number of eigenvectors/eigenvalues\n"
  "                            that will be stored.\n"
  "\n"
  "\n"
  "    GWAS\n"
  "    -------------------------------------\n"
  "\n"
  "    --group-all             Fit all SNPs in each genotype file together.\n"
  "\n"
  "\n"
  "    Phenotype simulation option\n"
  "    -------------------------------------\n"
  "\n"
  "    --effect-sizes f        Specify the effect sizes file.\n"
  "\n"
  "    --simu-h2               Specify the heritability.\n"
  "\n"
  "    --prevalence            Specify the prevalence.\n"
  "\n"
  "    --simu-quantitative     Simulate quantitative trait.\n"
  "\n"
  "    --simu-binary           Simulate binary trait.\n"
  "\n"
  "\n"
  "    Phenotype prediction options\n"
  "    -------------------------------------\n"
  "\n"
  "    --snp-effects f         Specify the file with SNPs effects.\n"
  "\n"
  "\n"
  "    Regional analysis\n"
  "    -------------------------------------\n"
  //"\n"
  //"    --region-size n         Specify region size for positional regional\n"
  //"                            analysis.\n"
  //"\n"
  //"    --region-overlap n      Specify region overlapping for positional regional\n"
  //"                            analysis.\n"
  //"\n"
  "    --groups f              Specify region groups for regional analysis.\n"
  "\n"
  "    --min-snps-region n     Specify the minimum elements a region must have.\n"
  "                            Otherwise will be ignored. (Default n = 1)\n"
  "\n" << std::endl;
}
