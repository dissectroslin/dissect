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

#ifndef GWAS_H
#define GWAS_H

#include "matrix.h"
#include "genotype.h"
#include "covariate.h"
#include "phenotype.h"
#include "kernel.h"
#include "labeledmatrix.h"

#include <vector>
#include <fstream>

class Matrix;
class GRM;
class Genotype;
class REML;
class CovarianceMatrix;
class Message;

enum ModelType
{
  emptyModelType,
  OLSModelType,
  REMLModelType,
  MLModelType
};

class GLMResults
{
public:
  ModelType type;
  
  bool success;
  
  std::vector<double> b;
  std::vector<double> SE;
  std::vector<double> tStatistics;              ///<Results for ordinary least squares
  std::vector<double> tStatisticPValues;        ///<Results for ordinary least squares
  
  std::vector<double> chi2Statistics;           ///<Results for REML
  std::vector<double> chi2StatisticsPValues;    ///<Results for REML
  
  double btXty;
  
  double SSE;
  double MSE;
  
  double SSR;
  double MSR;
  
  double FStatistic;
  double FStatisticPValue;                      ///<Depending of the model used, this can be a F statistic pvalue or a chi2 pvalue
  
  double groupVariance;
  
  double logLikelihood;
  
  GLMResults()
  {
    this->type = emptyModelType;
    
    this->success = false;
    
    this->b.clear();
    this->SE.clear();
    this->tStatistics.clear();
    this->tStatisticPValues.clear();
    
    this->chi2Statistics.clear();
    this->chi2StatisticsPValues.clear();
    
    this->btXty = 0.;
    
    this->SSE = 0.;
    this->MSE = 0.;
    
    this->SSR = 0.;
    this->MSR = 0.;
    
    this->FStatistic = 0.;
    this->FStatisticPValue = -1;
    
    this->groupVariance = -1.;
  }
};

class GWAS
{
public:
  
  std::vector<std::string> genotypeFiles;
  std::map<std::string, std::string> grmFiles;  ///<Map that associates a GRM file (second value in the map) to a genotype (first value/key on the map)
  std::map<std::string, std::string> residualFiles;  ///<Map that associates a Residual phenotype file (second value in the map) to a genotype (first value/key on the map)
  std::string currentFile;
  std::string currentGRMFile;
  std::string currentResidualFile;
  GenotypeAttributes::Type genotypeFilesType;
  
  Kernel * currentGRMBase;
  
  bool useCovariateMatrix;
  Kernel * currentCovariance;
  Matrix * currentCovarianceMatrix;             ///<A reference (not a copy) of the normalized covariance matrix in currentCovariance.
  //double accumulatedh2;
  std::vector<double> accumulatedVariance;      ///<Stores the accumulatedVariance of all REML fits.
  int nTests;
  
  Matrix * y;                                   ///<Phenotypes matrix.
  Matrix * X;                                   ///<Covariates matrix.
  Matrix * R;                                   ///<Covariance matrix.
  
  double significanceThreshold;                 ///<Significance over which SNPs will be stored in significantSNPIds
  std::vector<std::string> significantSNPIds;   ///<Ids of the significant SNPs
  std::set<std::string> correlatedSNPIds;       ///<Ids of the SNPs which have a correlated SNP more significant than them.
  
  GWAS();
  ~GWAS();
  
  void computeGWAS();
  void computeGroupedGWAS(Genotype * genotype, Phenotype * phenotype, Covariate * covariate);
  void computeIndividualGWAS(Genotype * genotype, Phenotype * phenotype, Covariate * covariate);
  void computeIndividualGWASGroupedCommunicator(Genotype * gcgenotype, Phenotype * gcphenotype, Covariate * gccovariate);
  
  bool computeGLM(Matrix * y, Matrix * X, Matrix * b, GLMResults & results);
  bool computeGLMWithoutCovariance(Matrix * y, Matrix * X, Matrix * b, GLMResults & results);
  bool computeGLMWithCovariance(Matrix * y, Matrix * X, Matrix * b, GLMResults & results);
  bool computeGroupSignificance(Matrix * y, GLMResults & reducedResults, GLMResults & results);
  void computeGroupVariance(GLMResults & results, Matrix * genotypes, LabeledMatrix * groupEffects, int idx);
  
  void storeResults(std::map<std::string, GLMResults > & effects, Covariate * covariate, std::map<std::string, std::vector<SNP> > & effectsSNPs);
  
  /**
   * Update to class variable correlatedSNPIds those SNPs correlated with other more significant SNP.
   * 
   * Genotypes must be normalized using the proper normalization.
   * 
   * \param threshold Threshold over which two SNPs are considered to be correlated. If threshold <= 0. do nothing.
   * \param genotypes Genotypes of the SPNs to test.
   * \param results p-values of the SNPs to test
   * \param shift Number of covariates (including mean) used in the adjustment. This is for knowing where in results there is SNP data.
   */
  void getLessSignificantCorrelatedSNPs(double threshold, Genotype* genotypes, GLMResults results, int shift);
  
  
  /**
   * Gather the results of a GWAS analysis after using a grouped communicator.
   * 
   * Assumes an individual SNP analysis.
   * 
   * \param[out] results The gathered results in the root node
   * \param[in] groupedResults The results in each communicator group
   * \param[out] resultSNPInfo The info of the SNPs in the results
   * \param[in] genotype Pointer to the Genotype class used for performing the analysis with the SNP info
   * \param[out] unfittedSNPs List with SNPs unfitted.
   */
  void gatherResults(std::map<std::string, GLMResults > & results, std::map<int, GLMResults > & groupedResults, std::map<std::string, std::vector<SNP> > & resultSNPInfo, Genotype * genotype, std::vector<std::string> & unfittedSNPs);
  
  bool computeCovariance(Genotype * genotype, Phenotype * phenotype, Covariate * covariate);
  
  void loadGenotypeResidualFiles();
  //These two pairs of two functions perform the same, but distributing the data between nodes differently.
  void computeMultiplePhenotypeGWAS();
  void computeIndividualGWASMultiplePhenotypesGroupedCommunicator(Genotype * gcgenotype, LabeledMatrix * gcphenotypes);
  //
  void computeMultiplePhenotypeGWASGroupedCommunicator();
  void computeIndividualGWASMultiplePhenotypes(Genotype * gcgenotype, LabeledMatrix * gcphenotypes, Communicator * globalCommunicator);
  /////////////////////
  bool computeGLMWithoutCovarianceMultiplePhenos(double sXtX, std::vector<double> & phenosXtys, std::vector<double> & phenoytys, int nIndividuals, GLMResults & results);
  void storeResultsMultiplePhenotype(std::map<std::string, GLMResults > & effects, std::map<std::string, std::vector<SNP> > & effectsSNPs, std::vector<std::string> & phenoLabels);
  void storeResultsMultiplePhenotype(std::map<std::string, GLMResults > & effects, std::map<std::string, std::vector<SNP> > & effectsSNPs, std::vector<std::string> & phenoLabels, std::ostream & messageSNPs, bool outputHeader);
  
  void debugWrite(Matrix * y, Matrix * X, Matrix * b);
};

#endif
