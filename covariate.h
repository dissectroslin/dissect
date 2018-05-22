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

#ifndef COVARIATE_H
#define COVARIATE_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include "matrix.h"

class Covariate
{
public:
  Matrix *covariates;
  
  int nCovariates;
  int nQuantitativeCovariates;
  int nDiscreteCovariates;
  int nIndividuals;
  
  std::vector<std::string> meanNames;                                                   ///<Name of the means.
  std::vector<std::string> quantitativeCovarNames;
  std::vector<std::string> discreteCovarNames;
  
  std::vector<std::string> individualIds;
  std::map<std::string, int> individualIdsIdx;
  
  std::set<std::string> individualIdsWithMissingData;                                   ///<Set with the individual Ids with missing data in at least one covariate.
  
  std::vector< std::vector<std::string> > rawCovars;                                    ///<Raw covariates directly copied from file
  std::vector< std::vector<std::string> > rawQCovars;                                   ///<Raw quantitative covariates directly copied from file
  std::vector< std::map<std::string, int> > covarCategories;                            ///<Discrete categories of the columns in rawCovars. This could be sincronized between different Covar objects.
  bool rawCleared;                                                                      ///<Is there data in rawCovars and rawQCovars?

  Covariate(std::string fcov, std::string fqcov, std::vector<std::string> emptyIndividualIds, bool constructMatrix = true, std::set<std::string> prefilterIndividuals = std::set<std::string>());
  ~Covariate();
  
  /**
   * Read the covariates from a file
   * 
   * Read the covariates from a file and store it in covars as appears in the file.
   * 
   * \param f File from where covariates will be read
   * \param[out] covars Return covariates in the file in a table as appear in the file.
   * \param prefilterIndividuals If it is not empty, the individuals not in this set, will be ignored when loading the data.
   */
  void readRawCovariate(std::string f, std::vector< std::vector<std::string> > & covars, std::set<std::string> prefilterIndividuals);
  
  /**
   * Parse raw covariates in rawCovars and rawQCovars and construct the covariates matrix.
   * 
   * \param emptyIndividualIds The individual ids that will be used in case no covariate files were specified.
   * \param nMeans The number of means in the covariates matrix (e.g. when constructing bivar covariate, this must be 2)
   * \param idxThisMean Each covariate for each trait will have 1's in their mean and 0's otherwise. This parameter specifies in which column it will have 1's. Must be between 0 and nMenans - 1
   */
  void parseRawCovariates(std::vector<std::string> emptyIndividualIds, int nMeans = 1, int idxThisMean = 0);
  
  /**
   * Get the categories in each column of rawCovars.
   * 
   * Get the categories in each column of rawCovars and store it in covarCategories.
   * 
   * \param rawCovars The table of discrete covariates from which categories will be created.
   */
  void getDiscreteCovariateCategories(std::vector< std::vector<std::string> > & rawCovars);
  
  /**
   * Syncronize covarCategories variable between two Covariate objects.
   * 
   * Syncronize covarCategories variable between two Covariate objects. Both objects must be created setting the parameter constructMatrix = false in the constructor.
   * After calling this function, parseRawCovariates() method must be called to construct the covariates matrix.
   * 
   * \param covariate Covariate object used for syncronize covarCategories.
   */
  void syncronizeDiscreteCovariateCategoriesWith(Covariate * covariate);
  
  /**
   * Erase data in rawCovars and rawQCovars
   */
  void clearRaw();
  
  /**
   * Process the quantitative data in rawCovars and store it in covars
   * 
   * \param[in] rawCovars Table with quantitative data for parsing
   * \param[out] covars Processsed data
   */
  void reestructureQuantitativeCovariate(std::vector< std::vector<std::string> > & rawCovars, std::vector< std::vector<double> > & covars);
  
  /**
   * Process the discrete data in rawCovars and store it in covars
   * 
   * Process the discrete data in rawCovars using information in covarCategories and store it in covars (generating a table of 0's and 1's depending whether the individual
   * is in a category or not).
   * 
   * \param[in] rawCovars Table with discrete data for parsing
   * \param[out] covars Processsed data
   */
  void reestructureDiscreteCovariate(std::vector< std::vector<std::string> > & rawCovars, std::vector< std::vector<double> > & covars);
  
  /**
   * Process the discrete data in rawCovars and store it in covars
   * 
   * Process the discrete data in rawCovars using information in covarCategories and store it in covars (generating a table of 0's and 1's depending whether the individual
   * is in a category or not). This method, as a difference of reestructureDiscreteCovariate, considers the zero category as the null cathegory and is eliminated from the matrix
   * In the results, the effects will be respect this category.
   * 
   * \param[in] rawCovars Table with discrete data for parsing
   * \param[out] covars Processsed data
   */
  void reestructureDiscreteCovariateUsingDifferences(std::vector< std::vector<std::string> > & rawCovars, std::vector< std::vector<double> > & covars);
  
  /**
   * Creates the covariates matrix from Covars and qCovars
   * 
   * \param[in] Covars Covars processed using the function reestructureDiscreteCovariate
   * \param[in] qCovars qCovars processed using the function reestructureQuantitativeCovariate
   * \param nMeans The number of means in the covariates matrix (e.g. when constructing bivar covariate, this must be 2)
   * \param idxThisMean Each covariate for each trait will have 1's in their mean and 0's otherwise. This parameter specifies in which column it will have 1's. Must be between 0 and nMenans - 1
   */
  void constructCovariateMatrix(std::vector< std::vector<double> > & Covars, std::vector< std::vector<double> > & qCovars, int nMeans, int idxThisMean);
  
  /**
   * Creates an estimation of the effect of the covariates on a list of individuals based on precomputed effects files.
   * 
   * The estimation can only be performed if the this->rawCovars and this->rawQCovars are present (this->rawCleared == false).
   * 
   * \param discretefname Name of the file with the discrete covariate effects.
   * \param quantitativefname Name of the file with the quantitative covariate effects.
   * \return A dict that associate an individual Id with the computed total effect for this individual.
   */
  std::map<std::string, double> loadEffectPrediction(std::string discretefname, std::string quantitativefname);
  
  /**
   * Helper function that loads discrete covariate effects.
   * 
   * The old version is used when old disect output files are detected.
   * 
   * \param fname File name with the effects
   */
  std::map< int, std::map<std::string, double> > loadDiscreteEffects(std::string fname);
  std::map< int, std::map<std::string, double> > loadDiscreteEffectsOld(std::string fname);
  
  /**
   * Helper function that loads quantitative covariate effects.
   * 
   * The old version is used when old disect output files are detected.
   * 
   * \param fname File name with the effects
   */
  std::vector<double> loadQuantitativeEffects(std::string fname);
  std::vector<double> loadQuantitativeEffectsOld(std::string fname);
  
  /**
   * Filter individuals from covariate matrix
   * 
   * \param keepIndividualIds Vector with the keys of individuals to keep.
   */
  void filterIndividuals(std::vector<std::string> keepIndividualIds);
  
  void printCovariate(int nSetw = 12);
  
  
  
};

#endif
