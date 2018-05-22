#include "accuracybysnp.h"
#include "phenotype.h"
#include "genotype.h"
#include "options.h"
#include "misc.h"
#include "predictphenotype.h"
#include "options.h"
#include "global.h"
#include "communicator.h"
#include "auxiliar.h"


#include <cmath>
#include <vector>
#include <string>

AccuracyBySNP::AccuracyBySNP()
{
  this->genotype = loadGenotypeUsingOptions();;
  
  this->realPhenotypes = new Phenotype(cyclicDistribution, options.phenotypesFile, options.phenotypeColumn);
  std::vector<std::string> tempOrderedIndividuals = orderVectorAsTemplate(this->genotype->individualIds, this->realPhenotypes->individualIds);
  if( communicator->mpiRoot == true && tempOrderedIndividuals.size() < 10 )
  {
    misc.error("Error: Less than 10 individuals shared between the genotypes file and the phenotypes file when computing snp dependent prediction.", 0);
  }
  this->realPhenotypes->filterIndividuals(tempOrderedIndividuals);
  
  Genotype *predictGenotypes = new Genotype();
  this->genotype->filterSNPsAndIndividuals(this->genotype->SNPIds, tempOrderedIndividuals, true, predictGenotypes); //This is important because in the following functions it is assumed that the individuals in the real phenotypes and the predicted ones are the same.. This could be improved for not requiring two copies of the genotypes.
  this->genotype->filterSNPsAndIndividuals(this->genotype->SNPIds, tempOrderedIndividuals); //This is important because in the following functions it is assumed that the individuals in the real phenotypes and the predicted ones are the same.
  
  
  this->predictPhenotype = new PredictPhenotype(predictGenotypes, options.snpEffectsFile);
  this->predictPhenotype->predictPhenotypes();
  predictGenotypes = NULL;
  
  std::vector<std::string> tempOrderedSNPs = orderVectorAsTemplate(this->genotype->SNPIds, this->predictPhenotype->effectSNPIds);
  if(tempOrderedSNPs != this->predictPhenotype->effectSNPIds)
  {
    misc.error("Error: There is not specified all the genotypes for all the SNPs with effects or they are not in the same order.", 0);
  }
  this->genotype->filterSNPsAndIndividuals(tempOrderedSNPs, this->genotype->individualIds);
  
  computeAccuracies();
}

AccuracyBySNP::~AccuracyBySNP()
{
  if ( this->genotype != NULL )
  {
    delete this->genotype;
    this->genotype = NULL;
  }
  if ( this->predictPhenotype != NULL )
  {
    delete this->predictPhenotype;
    this->predictPhenotype = NULL;
  }
  if ( this->realPhenotypes != NULL )
  {
    delete this->realPhenotypes;
    this->realPhenotypes = NULL;
  }
}

void AccuracyBySNP::computeAccuracies()
{
  if( this->genotype->SNPIds != this->predictPhenotype->effectSNPIds )
  {
    misc.error("Error: An internal error was happened. Discordant SNP Ids in genotypes and effects array. Accuracy by SNP efects cannot be computed.", 0);
  }
  
  Matrix * snpBLUPs = new Matrix(this->genotype->genotypes);
  
  //snpBLUPs matrix will contain the standardized blups. Individuals in columns, SNPs in rows. Each row contains the BLUP assuming that the corresponding SNP is removed.
  double * effects = new double [this->predictPhenotype->matrixEffects->nGlobRows];
  std::vector<double> gvEffects;
  this->predictPhenotype->matrixEffects->gatherMatrix(effects);
  this->predictPhenotype->matrixEffects->matrixToStandardVector(gvEffects);
  snpBLUPs->scatterVector(effects, row);
  double * vEff = snpBLUPs->v;
  snpBLUPs->v = NULL;
  delete [] effects;
  
  double * shift = new double [this->predictPhenotype->matrixShift->nGlobRows];
  this->predictPhenotype->matrixShift->gatherMatrix(shift);
  snpBLUPs->scatterVector(shift, row);
  double * vShift = snpBLUPs->v;
  snpBLUPs->v = NULL;
  delete [] shift;
  
  snpBLUPs->scatterVector(&(this->predictPhenotype->globalPhenotypes[0]), column);
  double * vPredPhenotypes = snpBLUPs->v;
  snpBLUPs->v = NULL;
  
  #pragma omp parallel for
  for(int c = 0; c<snpBLUPs->nCols; c++)
  {
    for(int r = 0; r<snpBLUPs->nRows; r++)
    {
      int temp = int(snpBLUPs->m[c*snpBLUPs->nRows + r]!=0.);
      snpBLUPs->m[c*snpBLUPs->nRows + r] = double(temp)*( ((snpBLUPs->m[c*snpBLUPs->nRows + r] - 1)*vEff[r]) + vShift[r] );
      snpBLUPs->m[c*snpBLUPs->nRows + r] = vPredPhenotypes[c] - snpBLUPs->m[c*snpBLUPs->nRows + r];
    }
  }
  delete [] vEff;
  delete [] vShift;
  delete [] vPredPhenotypes;
  
  Matrix * colOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, snpBLUPs->nGlobCols, 1);
  colOnes->fillWithConstant(1.);
  
  Matrix * mMeans = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  mMeans->multiply(snpBLUPs, 'N', colOnes, 'N', 1./double(colOnes->nGlobRows));
  double * means = new double [mMeans->nGlobRows];
  mMeans->gatherMatrix(means);
  snpBLUPs->scatterVector(means, row);
  delete mMeans;
  delete [] means;
  
  #pragma omp parallel for
  for(int c = 0; c<snpBLUPs->nCols; c++)
  {
    for(int r = 0; r<snpBLUPs->nRows; r++)
    {
      snpBLUPs->m[c*snpBLUPs->nRows + r] -= snpBLUPs->v[r];
    }
  }
  delete [] snpBLUPs->v;
  snpBLUPs->v = NULL;

  Matrix * snpBLUPsSquared = new Matrix(snpBLUPs);
  snpBLUPsSquared->elementWiseMultiplication(snpBLUPs);
  Matrix * mVars = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  mVars->multiply(snpBLUPsSquared, 'N', colOnes, 'N', 1./(double(colOnes->nGlobRows) - 1.));
  double * vars = new double [mVars->nGlobRows];
  mVars->gatherMatrix(vars);
  snpBLUPs->scatterVector(vars, row);
  delete snpBLUPsSquared;
  delete mVars;
  delete [] vars;
  
  #pragma omp parallel for
  for(int c = 0; c<snpBLUPs->nCols; c++)
  {
    for(int r = 0; r<snpBLUPs->nRows; r++)
    {
      snpBLUPs->m[c*snpBLUPs->nRows + r] /= sqrt(snpBLUPs->v[r]);  //This can be improved by precomputing the sqrts.
    }
  }
  delete [] snpBLUPs->v;  
  snpBLUPs->v = NULL;
  
  delete colOnes;
  
  
  //Standardize the phenotypes (I am not removing the mean, but I think it is not necessary).
  standardizeVector(this->realPhenotypes->phenotypes);
  
  //Compute the SNP dependent accuracies.
  Matrix * SNPDependentAccuracies = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  SNPDependentAccuracies->multiply(snpBLUPs, 'N', this->realPhenotypes->phenotypes, 'N', 1./double(this->realPhenotypes->phenotypes->nGlobRows));
  delete snpBLUPs;
  
  
  std::vector<double> gSNPDependentAccuracies;
  SNPDependentAccuracies->matrixToStandardVector(gSNPDependentAccuracies);
  delete SNPDependentAccuracies;
  
  //Compute the acccuracy using all SNPs
  double accuracy = computeAccuracy(this->predictPhenotype, this->realPhenotypes->phenotypes);
  misc.message << "The prediction accuracy using all SNP effects is: " << std::setprecision(14) <<  accuracy << std::endl;
  
  //Store results
  if(communicator->mpiRoot)
  {
    if( gSNPDependentAccuracies.size() != this->genotype->SNPs.size() )
    {
      misc.error("Error: An internal error was happened when storing SNP dependent accuracies.", 0);
    }
    
    Message message(options.outFile + ".snps.accuracies");
    message << "SNP ALLELE STDEV MEAN EFFECT CORR DELTA" << std::endl;
    for(int i = 0; i < gSNPDependentAccuracies.size(); i++)
    {
      message << this->genotype->SNPs[i].name << " " << this->genotype->SNPs[i].allele2;
      message << " " << std::setprecision(14) << this->genotype->SNPs[i].standardDev;
      message << " " << std::setprecision(14) << 2.*this->genotype->SNPs[i].p2;
      message << " " << std::setprecision(14) << gvEffects[i];
      message << " " << std::setprecision(14) << gSNPDependentAccuracies[i];
      message << " " << std::setprecision(14) << accuracy - gSNPDependentAccuracies[i];
      message << std::endl;
    }
  }
  
  bool stop = false;
  double step = 0.1;
  double currentStdScale = 3.;
  double previousAccuracy = accuracy;
  while(stop == false)
  {
    accuracy = accuracyFilteringAt(gSNPDependentAccuracies, currentStdScale);
    
    if(accuracy < previousAccuracy)
    {
      stop = true;
    }
    
    currentStdScale -= step;
    previousAccuracy = accuracy;
  }
}

double AccuracyBySNP::computeAccuracy(PredictPhenotype * pp, Matrix * rp)
{
  int temp = pp->globalPhenotypes.size();
  communicator->broadcast(&temp);
  Matrix * predictedPhenotypes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, 1, temp);
  predictedPhenotypes->scatterMatrix(&(pp->globalPhenotypes[0]));
  standardizeVector(predictedPhenotypes);
  
  Matrix * realPhenotypes = new Matrix(rp);
  standardizeVector(realPhenotypes);
  
  Matrix *totalAccuracy = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  totalAccuracy->multiply(predictedPhenotypes, 'N', this->realPhenotypes->phenotypes, 'N', 1./double(this->realPhenotypes->phenotypes->nGlobRows));
  
  std::vector<double> gTotalAccuracy;
  totalAccuracy->matrixToStandardVector(gTotalAccuracy);
  
  delete predictedPhenotypes;
  delete totalAccuracy;
  
  double accuracy;
  if(communicator->mpiRoot == true)
  {
    accuracy = gTotalAccuracy[0];
  }
  communicator->broadcast(&accuracy);
  
  return accuracy;
}

void AccuracyBySNP::standardizeVector(Matrix * vector)
{
  double vMean = untestedComputeMean(vector);
  double vStd = sqrt(computeVariance(vector));
  #pragma omp parallel for
  for(int c = 0; c<vector->nCols; c++)
  {
    for(int r = 0; r<vector->nRows; r++)
    {
      vector->m[c*vector->nRows + r] -= vMean;
      vector->m[c*vector->nRows + r] /= vStd;
    }
  }
}

double AccuracyBySNP::accuracyFilteringAt(std::vector<double> SNPAccuracies, double stdScaleThreshold)
{
  int temp = SNPAccuracies.size();
  communicator->broadcast(&temp);
  Matrix * mSNPAccuracies = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, 1, temp);
  mSNPAccuracies->scatterMatrix(&(SNPAccuracies[0]));
  //standardizeVector(predictedPhenotypes);
  
  double mean = untestedComputeMean(mSNPAccuracies);
  double std = sqrt(computeVariance(mSNPAccuracies));
  
  double threshold = mean + (std*stdScaleThreshold);
  
  std::vector<std::string> keepSNPIds;
  for(int i = 0; i<this->genotype->SNPIds.size(); i++)
  {
    if(SNPAccuracies[i] < threshold)
    {
      keepSNPIds.push_back(this->genotype->SNPIds[i]);
    }
  }
  
  if( misc.gt(keepSNPIds.size() < 1) )
  {
    return -1.;
  }
  
  Genotype * filteredGenotypes = new Genotype();
  this->genotype->filterSNPsAndIndividuals(keepSNPIds, this->genotype->individualIds, true, filteredGenotypes);
  
  PredictPhenotype * tempPredictPhenotype = new PredictPhenotype(filteredGenotypes, options.snpEffectsFile);
  tempPredictPhenotype->predictPhenotypes();
  filteredGenotypes = NULL;

  double accuracy = computeAccuracy(tempPredictPhenotype, this->realPhenotypes->phenotypes);
  
  delete tempPredictPhenotype;
  delete mSNPAccuracies;

  misc.message << int(this->genotype->SNPIds.size()) - int(keepSNPIds.size()) << " SNPs removed at threshold " << threshold << ". " << tempPredictPhenotype->effectSNPIds.size() << " SNPs are kept." << std::endl;  
  misc.message << "The prediction accuracy using filtered SNPs is " << std::setprecision(14) <<  accuracy << std::endl;
  
  return accuracy;
}