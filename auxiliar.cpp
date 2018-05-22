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

#include "auxiliar.h"
#include "global.h"
#include "misc.h"
#include "options.h"
#include "matrix.h"
#include "kernel.h"
#include "genotype.h"
#include "phenotype.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cstdarg>
#include <cmath>

#ifdef BOOSTLIB
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#endif

#if defined(BOOSTLIB) && defined(ZLIB)
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#endif


std::vector<std::string> intersectionStringVectors(int count, ...)
{
  std::vector<std::string> vResult;
  
  va_list vectors;
  va_start(vectors, count);
  
  vResult = *va_arg(vectors, std::vector<std::string>*);
  std::sort(vResult.begin(), vResult.end());
  for(int i = 1; i < count; i++)
  {
    std::vector<std::string> vTempResult;
    std::vector<std::string> vAdd = *va_arg(vectors, std::vector<std::string>*);
    std::sort(vAdd.begin(), vAdd.end());
    
    std::set_intersection(vResult.begin(), vResult.end(), vAdd.begin(), vAdd.end(), std::back_inserter(vTempResult));
    vResult = vTempResult;
  }
  
  return vResult;
}

std::vector<std::string> orderVectorAsTemplate(std::vector<std::string> & templateVector, std::vector<std::string> & sourceVector)
{
  std::set<std::string> sourceSet(sourceVector.begin(), sourceVector.end());
  std::vector<std::string> resultVector;
  for(int i = 0; i<templateVector.size(); i++)
  {
    if(sourceSet.find(templateVector[i]) != sourceSet.end())
    {
      resultVector.push_back(templateVector[i]);
    }
  }
  
  return resultVector;
}

std::vector<std::string> intersection2StringVectorsOrSecond(std::vector<std::string> base, std::vector<std::string> adding, std::string error)
{
  std::vector<std::string> result;
  if( base.size() == 0 )
  {
    result = adding;
  }
  else
  {
    result = intersectionStringVectors(2, &base, &(adding));
  }
  
  if( result.size() == 0 )
  {
    misc.error(error, 0);
  }
  
  return result;
}

std::vector<std::string> differenceBetweenTwoVectors(std::vector<std::string> & baseVector, std::vector<std::string> & vectorToSubstract)
{
  std::vector<std::string> test = orderVectorAsTemplate(vectorToSubstract, baseVector);
  if(test != vectorToSubstract)
  {
    misc.error("Error: Invalid operation. Two vectors cannot be substracted if all the elements of one are not contained into the other.", 0);
  }
  
  std::set<std::string> substractSet(vectorToSubstract.begin(), vectorToSubstract.end());
  std::set<std::string> baseSet(baseVector.begin(), baseVector.end());
  if(substractSet.size() != vectorToSubstract.size() || baseSet.size() != baseVector.size() )
  {
    misc.error("Error: Invalid operation. Two vectors cannot be substracted if one of them have repeated elements.", 0);
  }
  std::vector<std::string> resultVector;
  for(int i = 0; i<baseVector.size(); i++)
  {
    if(substractSet.find(baseVector[i]) == substractSet.end())
    {
      resultVector.push_back(baseVector[i]);
    }
  }
  
  return resultVector;
}

std::vector<std::string> getRandomSample(std::vector<std::string> & sourceVector, int nElements)
{
  if(nElements >= sourceVector.size() || nElements < 1)
  {
    misc.error("Error: An internal error was happened. Wrong number of elements for a random sample.", 0);
  }
  
  std::vector<std::string> shuffledVector = sourceVector;
  
  std::random_shuffle ( shuffledVector.begin(), shuffledVector.end() );
  
  std::vector<std::string> result(shuffledVector.begin(), shuffledVector.begin() + nElements);
  return orderVectorAsTemplate(sourceVector, result);
}

std::vector<int> extractMapValues(std::vector<std::string> & keysVector, std::map<std::string, int> & valuesMap)
{
  std::vector<int> result;
  
  for(int i = 0; i < keysVector.size(); i++)
  {
    if(valuesMap.count(keysVector[i]) == 0)
    {
      misc.error("Error: An internal error was happened. The key is not in the map.", 0);
    }
    result.push_back( valuesMap[ keysVector[i] ] );
  }
  
  return result;
}

void getListFromFile(std::string f, std::vector<std::string> & list)
{
  list.clear();
  std::ifstream file;
  misc.checkFileExists(f);
  file.open(f.c_str());
  
  std::string line;
  while(getline(file,line))
  {
    std::istringstream sstemp(line);
    std::string element;
    sstemp >> element;
    if(element != "")
    {
      list.push_back(element);
    }
  }

  file.close();
  
  if(list.size() == 0)
  {
    misc.error("Error: File [ " + f + " ] is empty.", 0);
  }
}

void getListFromFile(std::string f, std::vector<double> & list)
{
  std::vector<std::string> temp;
  
  getListFromFile(f, temp);
  
  list.clear();
  
  for(int i = 0; i<temp.size(); i++)
  {
    double value;
    std::stringstream sstemp(temp[i]);
    if( (sstemp >> value).fail() )
    {
      misc.error("Error: The value " + temp[i] + " in file [ " + f + " ] is not a valid number.", 0);
    }
    list.push_back(value);
  }
}

void getTableFromFile(std::string f, std::vector< std::vector<std::string> > & table, int nColumns)
{
  table.clear();
  std::ifstream file;
  misc.checkFileExists(f);
  file.open(f.c_str());
  
  std::string line;
  while(getline(file,line))
  {
    std::istringstream sstemp(line);
    std::string element;
    std::vector<std::string> row;
    int nElements = 0;
    while((sstemp >> element) && nElements < nColumns)
    {
      if(element != "")
      {
        row.push_back(element);
        nElements++;
      }
    }
    if(row.size() ==  0)
    {
      continue;
    }
    if(row.size() != nColumns)
    {
      misc.error(" Error: Line: \n " + line + "\n in file [ " + f + " ] does not have at least " + i2s(nColumns) + " elements.", 0);
    }
    table.push_back(row);
  }

  file.close();
  
  if(table.size() == 0)
  {
    misc.error("Error: File [ " + f + " ] is empty.", 0);
  }
}

int getNumberOfFileColumns(std::string f)
{
  std::ifstream file;
  std::string line;

  misc.checkFileExists(f);
  file.open(f.c_str());
  
  getline(file,line);
  if(!file)
  {
    misc.error("Error: File [ " + f + " ] is empty.", 0);
  }
  
  std::istringstream sstemp(line);
  std::string dummy;
  int nColumns = 0;
  while(sstemp >> dummy)
  { 
    nColumns++;
  }
    
  if(nColumns < 1)
  {
    misc.error("Error: File [ " + f + " ] is not properly formated.", 0);
  }
  
  file.close();

  return nColumns;
}

std::vector<std::string> getHeader(std::string f)
{
  std::ifstream file;
  std::string line;

  std::vector<std::string> header;
  
  misc.checkFileExists(f);
  file.open(f.c_str());
  
  getline(file,line);
  if(!file)
  {
    misc.error("Error: File [ " + f + " ] is empty.", 0);
  }
  
  std::istringstream sstemp(line);
  std::string element;
  while(sstemp >> element)
  { 
    header.push_back(element);
  }
    
  if(header.size() < 1)
  {
    misc.error("Error: First line of file [ " + f + " ] is empty.", 0);
  }
  
  file.close();

  return header;
}

void writeListToFile(std::string f, std::vector<std::string> & list)
{
  std::ofstream file;
  if(communicator->mpiRoot)
  {
    file.open(f.c_str(), std::ofstream::out);
    for(int i = 0; i<list.size(); i++)
    {
      file << list[i] << std::endl;
    }
    file.close();
  }
}

int leastCommonMultiple(int a, int b)
{
  if(a < 0 || b < 0)
  {
    misc.error("Error: An internal error was happened while computing the Least Common Multiple", 0);
  }
  return (a*b)/greatestCommonDivisor(a,b);
}

int greatestCommonDivisor(int a, int b)
{
  if(a < 0 || b < 0)
  {
    misc.error("Error: An internal error was happened while computing the Greatest Common Divisor", 0);
  }
  if(b == 0)
  {
    return a;
  }
  else
  {
    return greatestCommonDivisor(b, a % b);
  }
}

double untestedComputeMean(Matrix * m)
{
  int dimension;
  double *globalM;
  
  if(m->nGlobRows != 1 && m->nGlobCols != 1)
  {
    misc.error("Error: Unable to compute the variance of a matrix with more than 1 rows and columns.", 0);
  }
  
  if( (m->nGlobRows == 1 && m->nGlobCols <= 1) || m->nGlobRows <= 1 && m->nGlobCols == 1) 
  {
    misc.error("Error: Unable to compute the phenotypic variance of less than two elements.", 0);
  }
  
  if(m->nGlobRows > m->nGlobCols)
  {
    dimension = m->nGlobRows;
  }
  else
  {
    dimension = m->nGlobCols;
  }
  
  if (communicator->mpiRoot) {
    globalM = new double [dimension];
  }
  
  m->gatherMatrix(globalM);

  double avg = 0.;  
  if (communicator->mpiRoot)
  {
    for(int r = 0; r<dimension; r++)
    {
      avg += globalM[r];
    }
    avg /= double(dimension);
  }
  
  communicator->broadcast(&avg, 1);
  
  if (communicator->mpiRoot) {
    delete [] globalM;
  }
  
  return avg;
}

double computeVariance(Matrix * m)
{
  int dimension;
  double *globalM;
  double variance;
  
  if(m->nGlobRows != 1 && m->nGlobCols != 1)
  {
    misc.error("Error: Unable to compute the variance of a matrix with more than 1 rows and columns.", 0);
  }
  
  if( (m->nGlobRows == 1 && m->nGlobCols <= 1) || m->nGlobRows <= 1 && m->nGlobCols == 1) 
  {
    misc.error("Error: Unable to compute the phenotypic variance of less than two elements.", 0);
  }
  
  if(m->nGlobRows > m->nGlobCols)
  {
    dimension = m->nGlobRows;
  }
  else
  {
    dimension = m->nGlobCols;
  }
  
  if (communicator->mpiRoot) {
    globalM = new double [dimension];
  }
  
  m->gatherMatrix(globalM);
  
  if (communicator->mpiRoot)
  {
    double avg = 0.;
    for(int r = 0; r<dimension; r++)
    {
      avg += globalM[r];
    }
    avg /= double(dimension);
    variance = 0.;
    for(int r = 0; r<dimension; r++)
    {
      double temp = globalM[r] - avg;
      variance += temp * temp;
    }
    variance /= (double(dimension) - 1.);
  }
  
  communicator->broadcast(&variance, 1);
  
  if (communicator->mpiRoot) {
    delete [] globalM;
  }
  
  return variance;
}

template<typename T>
std::string getString(T i)
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}
template std::string getString(int i);
template std::string getString(double i);

std::string i2s(int i)
{
  return getString(i);
}

template<typename T>
T string2Number(std::string s, std::string errorMessage)
{
  T value;
  std::istringstream sstemp(s);
  if( (sstemp >> value).fail() )
  {
    misc.error(errorMessage, 0);
  }
  return value;
}
template int string2Number(std::string s, std::string errorMessage);
template double string2Number(std::string s, std::string errorMessage);
template long string2Number(std::string s, std::string errorMessage);

std::vector<std::string> splitString(std::string str, std::string delimiter)
{
  std::vector<std::string> result;
  std::string s = str;

  size_t pos = 0;
  std::string token;
  while ((pos = s.find(delimiter)) != std::string::npos) {
      token = s.substr(0, pos);
      result.push_back(token);
      s.erase(0, pos + delimiter.length());
  }
  result.push_back(s);
  
  return result;
}

std::string spacetab2underscore(std::string s)
{
  std::string text = s;
  std::replace(text.begin(), text.end(), ' ', '_');
  std::replace(text.begin(), text.end(), '\t', '_');
  return text;
}

void removeCharacters( std::string & str, std::string chars )
{
  for ( int i = 0; i < chars.size(); i++ )
  {
    str.erase( std::remove(str.begin(), str.end(), chars[i]), str.end() );
  }
}


int getMapValue(std::map<std::string, int> & smap, std::string & key)
{
  std::map<std::string, int>::iterator it = smap.find( key );
  if(it == smap.end())
  {
    misc.error("Error: An internal error has happened. Key not present in the map.", 0);
  }
  return it->second;
}


std::set<std::string> getIndividualIdsSpecifiedByOptionKeep()
{
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
    if(keepIndividualIds.size() == 0)
    {
      misc.message << "WARNING: There are no valid individuals specified in the file [ " << options.fileIndividualsToKeep << " ]. It will be ignored." << std::endl;
    }
  }
  return keepIndividualIds;
}

Genotype * loadGenotypeUsingOptions()
{
  Genotype *genotype;
  if(options.genotypeFile != "") //Read single genotype file
  {
    genotype = new Genotype(options.genotypeFile);
  }
  else if (options.genotypeListFile != "") //Read multiple genotype files
  {
    genotype = new Genotype();
    genotype->loadList(options.genotypeListFile);
  }
  else
  {
    misc.error("Error: No genotype file(s) specified.", 0);
  }
  return genotype;
}

Kernel * loadGRMUsingOptions(bool returnGenotype, Genotype **returnedGenotype)
{
  Kernel *grm = NULL;
  
  if(returnGenotype == true)
  {
    if(returnedGenotype == NULL)
    {
      misc.error("Error: An internal error was happened when loading the GRM. Wrong pointer for returning a genotype.", 0);
    }
    *returnedGenotype = NULL;
  }
  
  if( returnGenotype == true && options.GRMJoinMethod == 0 )
  {
    misc.message << "Changing GRM joining method to method 1." << std::endl;
    options.GRMJoinMethod = 1;
  }
  
  if(options.grmFile == "")
  {
    if(options.genotypeListFile == "") //Create GRM from single genotype file
    {
      Genotype *genotype = new Genotype(options.genotypeFile);
      grm = new Kernel(genotype);
      if(returnGenotype == false)
      {
        delete genotype;
      }
      else
      {
        *returnedGenotype = genotype;
      }
    }
    else //Create GRM from multiple genotype files
    {
      if(options.GRMJoinMethod == 0) //Create GRM using method 0: First compute GRMs for each genotype, then add GRMs.
      {
        misc.message << "Computing GRM using a list of genotype files. Method 0: First compute GRMs for each genotype, then add GRMs." << std::endl;
        misc.message.tab = "  ";
        misc.setGetElapsedTime("GRMPartComputation");
        std::ifstream file;
        misc.checkFileExists(options.genotypeListFile);
        file.open(options.genotypeListFile.c_str());
        
        int idx = 0;
        std::string line;
        std::vector<std::string> loadedSNPs;
        while(getline(file,line))
        {
          std::istringstream sstemp(line);
          std::string partialFileName;
          sstemp >> partialFileName;
          
          if( idx == 0 )
          {
            Genotype *genotype = new Genotype(line);
            loadedSNPs = genotype->SNPIds;
            
            grm = new Kernel(genotype, false);
            if(returnGenotype == false)
            {
              delete genotype;
            }
            else
            {
              misc.error("Error: An internal error was happened. GRM cannot be loaded using method 0 when genotypes must also be loaded.", 0);
            }
          }
          else
          {
            Genotype *genotype = new Genotype(line);
            std::vector<std::string> SNPsIntersection = intersectionStringVectors(2, &loadedSNPs, &genotype->SNPIds); //Check there are not intersections between SNPs in files.
            if(SNPsIntersection.size() != 0)
            {
              misc.error("Error: When computing a GRM from multiple genotype files. There are SNPs repeated in different files.", 0);
            }
            loadedSNPs.insert( loadedSNPs.end(), genotype->SNPIds.begin(), genotype->SNPIds.end() );
            
            Kernel *grmToAdd = new Kernel(genotype, false);
            delete genotype;
            misc.message << "Adding GRM..." << std::endl;
            grm->addKernels(1., grmToAdd);
            delete grmToAdd;
          }
          
          idx++;
        }
        file.close();
        misc.message.tab = "";
        grm->normalize();
        misc.message << "Genotype loading and GRM computation finished after " << misc.setGetElapsedTime("GRMPartComputation", true) << std::endl;
      }
      else  //Create GRM using method 1: First join genotypes, then compute GRM.
      {
        misc.message << "Computing GRM using a list of genotype files. Method 1: First join genotypes, then compute GRM." << std::endl;
        misc.message.tab = "  ";
        misc.setGetElapsedTime("GRMPartComputation");
        Genotype *genotype = new Genotype();
        genotype->loadList(options.genotypeListFile);
        grm = new Kernel(genotype);
        misc.message.tab = "";
        misc.message << "Genotype loading and GRM computation finished after " << misc.setGetElapsedTime("GRMPartComputation", true) << std::endl;
        if(returnGenotype == false)
        {
          delete genotype;
        }
        else
        {
          *returnedGenotype = genotype;
        }
      }
    }
  }
  else
  {
    grm = new Kernel(options.grmFile);
  }
  return grm;
}

void loadGRMUsingOptions(std::vector<Kernel *> & grms, bool returnGenotypes, std::map<std::string, Genotype*> & genotypesList, std::map<std::string, std::vector<std::string> > & genotypesFilesList, std::map<std::string, std::vector<std::string> > & SNPIds)
{
  if( options.grmListFile == "" && options.grmFile == "" && options.genotypeFile == "" && options.genotypeListFile == "" ) //Nothing to load, kernels should be defined with other options.
  {
    return;
  }
  
  if(options.grmListFile == "")
  {
    Genotype *genotypes = NULL;
    Kernel *grm = loadGRMUsingOptions(true, &genotypes);
    grm->name = options.baseVarianceNames["GRM"];
    grms.push_back(grm);
    
    if(options.computeSNPsBLUP)
    {
      if( genotypes == NULL && options.SNPBLUPGenotypeFiles.size() == 0 )
      {
        genotypes = loadGenotypeUsingOptions();
      }
      else if( genotypes == NULL )
      {
        if(genotypesFilesList.count(options.baseVarianceNames["GRM"]) != 0)
        {
          misc.error("Error: An internal error has happened. There is already a genotypes associated with the same GRM.", 0);
        }
        genotypesFilesList[ options.baseVarianceNames["GRM"] ] = options.SNPBLUPGenotypeFiles;
        SNPIds[ options.baseVarianceNames["GRM"] ] = grm->randomVarNames;
      }
      
      if( genotypes != NULL )
      {
        if(genotypesList.count(options.baseVarianceNames["GRM"]) != 0)
        {
          misc.error("Error: An internal error has happened. There is already a genotypes associated with the same GRM.", 0);
        }
        genotypesList[ options.baseVarianceNames["GRM"] ] = genotypes;
        SNPIds[ options.baseVarianceNames["GRM"] ] = grm->randomVarNames;
      }
    }
    else
    {
      if(genotypes != NULL)
      {
        delete genotypes;
        genotypes = NULL;
      }
    }
  }
  else
  {
    std::vector< std::vector<std::string> > listGRMFiles;
    getTableFromFile(options.grmListFile, listGRMFiles, 4);
    if(options.addAllGRMsInOne == false)
    {
      misc.message << "Loading GRMs from files specified in [ " << options.grmListFile << " ]..." << std::endl;
      misc.message.tab = "  ";
      for(int i = 0; i<listGRMFiles.size(); i++)
      {
        Kernel * grm = new Kernel(listGRMFiles[i][1]);
        grm->name = listGRMFiles[i][0];
        grms.push_back(grm);
        
        if(options.computeSNPsBLUP)
        {
          std::vector<std::string> tempList;
          if( listGRMFiles[i][2] == "N" )
          {
            misc.message << "There are not genotypes specified for the GRM " + listGRMFiles[i][0] << std::endl;
          }
          else if( listGRMFiles[i][2] == "F" )
          {
            tempList.push_back(listGRMFiles[i][3]);
            genotypesFilesList[ grm->name ] = tempList;
            SNPIds[ grm->name ] = grm->randomVarNames;
          }
          else if( listGRMFiles[i][2] == "L" )
          {
            getListFromFile(listGRMFiles[i][3], tempList);
            genotypesFilesList[ grm->name ] = tempList;
            SNPIds[ grm->name ] = grm->randomVarNames;
          }
          else
          {
            misc.error("Error: An error was happened. The file [ " + options.grmListFile + " ] contains " + listGRMFiles[i][2] + " in the third column. Valid values are N, F, or L. Please, consult the documentation for more information.", 0);
          }
        }
      }
      misc.message.tab = "";
    }
    else
    {
      misc.message << "Loading GRMs from files specified in [ " << options.grmListFile << " ] and adding all in a single GRM..." << std::endl;
      misc.message.tab = "  ";
      Kernel * grm = new Kernel(listGRMFiles[0][1]);
      for(int i = 1; i<listGRMFiles.size(); i++)
      {
        Kernel * grmToAdd = new Kernel(listGRMFiles[i][1]);
        
        std::vector<std::string> intersectionOrderedIndividuals = orderVectorAsTemplate(grm->individualIds, grmToAdd->individualIds);
        grm->filterIndividuals(intersectionOrderedIndividuals);
        grmToAdd->filterIndividuals(intersectionOrderedIndividuals);
        
        grm->addKernels(1., grmToAdd);
        delete grmToAdd;
      }
      grm->name = options.baseVarianceNames["GRM"];
      grms.push_back(grm);
      misc.message.tab = "";
      
      if( options.computeSNPsBLUP && (options.genotypeListFile != "" || options.genotypeFile != "") )
      {
        misc.message << "Genotype file(s) have been specified. Ignoring genotype files specified in [ " + options.grmListFile + " ] if any." << std::endl;
        Genotype *genotypes = loadGenotypeUsingOptions();
        genotypesList[ options.baseVarianceNames["GRM"] ] = genotypes;
        SNPIds[ options.baseVarianceNames["GRM"] ] = grm->randomVarNames;
      }
      else if( options.computeSNPsBLUP && options.SNPBLUPGenotypeFiles.size() != 0)
      {
        misc.message << "Genotype file(s) have been specified. Ignoring genotype files specified in [ " + options.grmListFile + " ] if any." << std::endl;
        genotypesFilesList[ options.baseVarianceNames["GRM"] ] = options.SNPBLUPGenotypeFiles;
        SNPIds[ options.baseVarianceNames["GRM"] ] = grm->randomVarNames;
      }
      else if( options.computeSNPsBLUP )
      {
        std::vector<std::string> tempList;
        std::set<std::string> tempSet;
        for(int i = 0; i<listGRMFiles.size(); i++)
        {
          if( listGRMFiles[i][2] == "N" )
          {
            //Do nothing
          }
          else if( listGRMFiles[i][2] == "F" )
          {
            tempSet.insert(listGRMFiles[i][3]);
          }
          else if( listGRMFiles[i][2] == "L" )
          {
            getListFromFile(listGRMFiles[i][3], tempList);
            tempSet.insert(tempList.begin(), tempList.begin());
          }
          else
          {
            misc.error("Error: An error was happened. The file [ " + options.grmListFile + " ] contains " + listGRMFiles[i][2] + " in the third column. Valid values are N, F, or L. Please, consult the documentation for more information.", 0);
          }
        }
        tempList.assign( tempSet.begin(), tempSet.end() );
        genotypesFilesList[ options.baseVarianceNames["GRM"] ] = tempList;
        SNPIds[ options.baseVarianceNames["GRM"] ] = grm->randomVarNames;
      }
    }
  }
}

void addKernelsUsingOptions(std::vector<Kernel *> & kernels, std::vector<std::string> & reducedModels, std::vector<std::string> & individualBLUPNames)
{
  //Get the intersection between the individuals of previous kernels.
  std::set<std::string> sharedIndividuals;
  std::vector<std::string> sharedIndividualsVector;
  if(communicator->mpiRoot && kernels.size() != 0)
  {
    sharedIndividualsVector = kernels[0]->individualIds;
    for(int i = 1; i<kernels.size(); i++)
    {
      sharedIndividualsVector = intersectionStringVectors(2, &sharedIndividualsVector, &(kernels[i]->individualIds));
    }
    sharedIndividuals = std::set<std::string>(sharedIndividualsVector.begin(), sharedIndividualsVector.end());
  }
  
  //Load new kernels.
  for(int icol = 0; icol < options.discreteRandomEffectsFileColumns.size(); icol++)
  {
    Kernel * covKernel = new Kernel(options.discreteRandomEffectsFile, kernelFromDiscreteCovariates, options.discreteRandomEffectsFileColumns[icol], sharedIndividuals);
    //covKernel->name = options.baseVarianceNames["cov-"] + i2s(icol);
    kernels.push_back(covKernel);
    reducedModels.push_back(covKernel->name);
    individualBLUPNames.push_back(covKernel->name);
    if(communicator->mpiRoot)
    {
      sharedIndividualsVector = intersection2StringVectorsOrSecond(sharedIndividualsVector, covKernel->individualIds, "Error: An error has happened. The intersection of the individuals in different covariance matrices (e.g. GRMs, random effects, etc) is empty.");
    }
  }
  for(int icol = 0; icol < options.multiDiscreteRandomEffectsFileColumns.size(); icol++)
  {
    Kernel * covKernel = new Kernel(options.multiDiscreteRandomEffectsFile, kernelFromMultiDiscreteCovariates, options.multiDiscreteRandomEffectsFileColumns[icol], sharedIndividuals);
    kernels.push_back(covKernel);
    reducedModels.push_back(covKernel->name);
    individualBLUPNames.push_back(covKernel->name);
    if(communicator->mpiRoot)
    {
      sharedIndividualsVector = intersection2StringVectorsOrSecond(sharedIndividualsVector, covKernel->individualIds, "Error: An error has happened. The intersection of the individuals in different covariance matrices (e.g. GRMs, random effects, etc) is empty.");
    }
  }
  for(int ifile = 0; ifile < options.squaredExponentialKernelFiles.size(); ifile++)
  {
    Kernel * sek = new Kernel(options.squaredExponentialKernelFiles[ifile], kernelSquaredExponential, -1, sharedIndividuals);
    sek->name = options.baseVarianceNames["SEK-"] + sek->name;
    kernels.push_back(sek);
    reducedModels.push_back(sek->name);
    individualBLUPNames.push_back(sek->name);
    if(communicator->mpiRoot)
    {
      sharedIndividualsVector = intersection2StringVectorsOrSecond(sharedIndividualsVector, sek->individualIds, "Error: An error has happened. The intersection of the individuals in different covariance matrices (e.g. GRMs, random effects, etc) is empty.");
    }
  }
  
  if(communicator->mpiRoot)
  {
    sharedIndividuals = std::set<std::string>(sharedIndividualsVector.begin(), sharedIndividualsVector.end());
  }
  
  //Create kernels from the previous ones, i.e. GxE kernels.
  if(options.includeGxEInteractions == true)
  {
    std::vector<Kernel *> interactionKernels;
    for(int i = 0; i < kernels.size(); i++)
    {
      if( kernels[ i ]->type == kernelGRM )
      {
        std::vector<std::string> tempIndividuals = orderVectorAsTemplate(kernels[ i ]->individualIds, sharedIndividualsVector);
        kernels[ i ]->filterIndividuals(tempIndividuals);
        for(int j = 0; j < kernels.size(); j++)
        {
          if( kernels[ j ]->type == kernelFromDiscreteCovariates && i != j)
          {
            kernels[ j ]->filterIndividuals(tempIndividuals);
            Kernel * inter = new Kernel(kernels[i], kernels[j]);
            interactionKernels.push_back(inter);
            reducedModels.push_back(inter->name);
            individualBLUPNames.push_back(inter->name);
          }
        }
      }
    }
    kernels.insert(kernels.end(), interactionKernels.begin(), interactionKernels.end());
  }
  
  //Load GCTA grms.
  if( options.GCTAGRMsFile != "" )
  {
    std::vector< std::vector<std::string> > GCTAGRMsFiles;
    getTableFromFile(options.GCTAGRMsFile, GCTAGRMsFiles, 2);
    for(int irow = 0; irow < GCTAGRMsFiles.size(); irow++)
    {
      Kernel * GCTAKernel = new Kernel(GCTAGRMsFiles[irow][1], kernelGCTAGRM);
      GCTAKernel->name = GCTAGRMsFiles[irow][0];
      kernels.push_back(GCTAKernel);
      reducedModels.push_back(GCTAKernel->name);
      individualBLUPNames.push_back(GCTAKernel->name);
      if(communicator->mpiRoot)
      {
        sharedIndividualsVector = intersection2StringVectorsOrSecond(sharedIndividualsVector, GCTAKernel->individualIds, "Error: An error has happened. The intersection of the individuals in different covariance matrices (e.g. GRMs, random effects, etc) is empty.");
      }
    }
  }
  
}

void introduceResortedGRMsByCouples(std::vector<Kernel*> & kernels, std::vector<std::string> & reducedModels, std::vector<std::string> & individualBLUPNames)
{
  if( options.indirectEffectsCouplesFile == "" )
  {
    return;
  }
  std::string couplesFile = options.indirectEffectsCouplesFile;
  
  std::map< std::string, std::string > couplesIds;
  std::map< std::string, std::pair<std::string, std::string> > relationIdsOriginalIds;

  //Read the couples
  if( communicator->mpiRoot == true )
  {
    std::vector< std::vector<std::string> > couplesTable;
    getTableFromFile(couplesFile, couplesTable, 4);
    
    std::set<std::string> testDuplicates1;
    std::set<std::string> testDuplicates2;
    
    for( int i = 0; i < couplesTable.size(); i++ )
    {
      std::string iid1 = couplesTable[i][0] + "@" + couplesTable[i][1];
      std::string iid2 = couplesTable[i][2] + "@" + couplesTable[i][3];
      
      if( testDuplicates1.find(iid1) != testDuplicates1.end() )
      {
        misc.error("Error: The individual '" + couplesTable[i][0] + " " + couplesTable[i][1]  + "' in file [ " + couplesFile + " ] is duplicated in the first column.", 0);
      }
      if( testDuplicates2.find(iid2) != testDuplicates2.end() )
      {
        misc.error("Error: The individual '" + couplesTable[i][2] + " " + couplesTable[i][3]  + "' in file [ " + couplesFile + " ] is duplicated in the second column.", 0);
      }
      
      /*if(srcKernel->individualIdsIdx.count(iid1) == 0)
      {
        misc.message << "Warning: The individual '" + couplesTable[i][0] + " " + couplesTable[i][1]  + "' in file [ " + couplesFile + " ] is not in the GRM matrix. Ignoring this couple.";
        continue;
      }
      if(srcKernel->individualIdsIdx.count(iid2) == 0)
      {
        misc.message << "Warning: The individual '" + couplesTable[i][2] + " " + couplesTable[i][3]  + "' in file [ " + couplesFile + " ] is not in the GRM matrix. Ignoring this couple.";
        continue;
      }*/
      couplesIds[iid1] = iid2;
      
      testDuplicates1.insert(iid1);
      testDuplicates2.insert(iid2);
      
      relationIdsOriginalIds[iid1] = std::pair<std::string, std::string>(couplesTable[i][0], couplesTable[i][1]);
      relationIdsOriginalIds[iid2] = std::pair<std::string, std::string>(couplesTable[i][2], couplesTable[i][3]);
    }
  }
  
  //Create indirect effect couples
  std::vector<Kernel*> newKernels = kernels;
  for(int i = 0; i < kernels.size(); i++)
  {
    Kernel * srcKernel = kernels[i];
    if( srcKernel->type != kernelGRM )
    {
      continue;
    }
    
    if(srcKernel->asymmetric == true)
    {
      misc.error("Error: An internal error was happened. This operation cannot be performed on an asymmetric GRM.", 0);
    }
    if(srcKernel->diagonalized == true)
    {
      misc.error("Error: An error was happened. Diagonalized GRM cannot be resorted using couples data.", 0);
    }
    
    std::vector<std::string> newOrderedIndividualsIds;
    std::vector<Individual> individualsKeepReplacementIds;

    //Create the vector with the reordered GRM individuals
    if( communicator->mpiRoot == true )
    {
      for(int i = 0; i < srcKernel->nIndividuals; i++)
      {
        if(couplesIds.count( srcKernel->individualIds[i] ) == 0)
        {
          misc.message << "WARNING: The individual (" + srcKernel->individuals[i].familyID + " " + srcKernel->individuals[i].individualID  + ") is in the " + srcKernel->name + " but it is not present in the first element of the couple in the file [ " + couplesFile + " ]. It will be removed from the " + srcKernel->name + "." << std::endl;
          continue;
        }
        std::string coupleId = couplesIds[srcKernel->individualIds[i]];
        if(srcKernel->individualIdsIdx.count(coupleId) == 0)
        {
          misc.message << "WARNING: The second individual of the couple (" + srcKernel->individuals[i].familyID + " " + srcKernel->individuals[i].individualID + ", " + relationIdsOriginalIds[coupleId].first + " " + relationIdsOriginalIds[coupleId].second  + ") is not in the " + srcKernel->name + ". The first individual of the couple will be removed from the " + srcKernel->name + "." << std::endl;
          continue;
        }
        newOrderedIndividualsIds.push_back( coupleId );
        individualsKeepReplacementIds.push_back( srcKernel->individuals[i] );
      }
    }
    
    if( misc.gt(newOrderedIndividualsIds.size()*4 <= srcKernel->individualIds.size()) )
    {
      misc.message << "WARNING: The " + srcKernel->name + " does not contain enough individuals defined in the file [ " + couplesFile + " ]. This GRM will be skiped for indirect effects analysis." << std::endl;
      continue;
    }
    
    Kernel * resultKernel = new Kernel(srcKernel);
    resultKernel->filterIndividuals(newOrderedIndividualsIds);
    resultKernel->replaceIndividualIds(individualsKeepReplacementIds);

    
    std::string tempName = "coup" + srcKernel->name;
    resultKernel->name = tempName;
    newKernels.push_back( resultKernel );
    reducedModels.push_back( tempName );
    individualBLUPNames.push_back( tempName );
  }
  kernels = newKernels;
}

void filterPhenotypesUsingCouples(std::vector<Phenotype*> & phenotypes)
{
  if( options.indirectEffectsCouplesFile == "" || phenotypes.size() != 2 )
  {
    return;
  }
  std::string couplesFile = options.indirectEffectsCouplesFile;
  
  std::map< std::string, std::string > couplesIds;
  std::map< std::string, std::pair<std::string, std::string> > relationIdsOriginalIds;

  //Read the couples
  if( communicator->mpiRoot == true )
  {
    std::vector< std::vector<std::string> > couplesTable;
    getTableFromFile(couplesFile, couplesTable, 4);
    
    std::set<std::string> testDuplicates1;
    std::set<std::string> testDuplicates2;
    
    for( int i = 0; i < couplesTable.size(); i++ )
    {
      std::string iid1 = couplesTable[i][0] + "@" + couplesTable[i][1];
      std::string iid2 = couplesTable[i][2] + "@" + couplesTable[i][3];
      
      if( testDuplicates1.find(iid1) != testDuplicates1.end() )
      {
        misc.error("Error: The individual '" + couplesTable[i][0] + " " + couplesTable[i][1]  + "' in file [ " + couplesFile + " ] is duplicated in the first column.", 0);
      }
      if( testDuplicates2.find(iid2) != testDuplicates2.end() )
      {
        misc.error("Error: The individual '" + couplesTable[i][2] + " " + couplesTable[i][3]  + "' in file [ " + couplesFile + " ] is duplicated in the second column.", 0);
      }
      
      couplesIds[iid1] = iid2;
      
      testDuplicates1.insert(iid1);
      testDuplicates2.insert(iid2);
      
      relationIdsOriginalIds[iid1] = std::pair<std::string, std::string>(couplesTable[i][0], couplesTable[i][1]);
      relationIdsOriginalIds[iid2] = std::pair<std::string, std::string>(couplesTable[i][2], couplesTable[i][3]);
    }
  }
  
  
  //Substitute indirect effect couples on phenotype 2
  for(int i = 1; i < 2; i++)
  {
    Phenotype * srcPhenotype = phenotypes[i];
    
    std::vector<std::string> newOrderedIndividualsIds;
    std::vector<std::string> individualsKeepReplacementIds;

    //Create the vector with the reordered Phenotype individuals
    if( communicator->mpiRoot == true )
    {
      for(int i = 0; i < srcPhenotype->nIndividuals; i++)
      {
        if(couplesIds.count( srcPhenotype->individualIds[i] ) == 0)
        {
          misc.message << "WARNING: The individual < " + srcPhenotype->individualIds[i] + " > is in the phenotypes file but it is not present in the first element of the couple in the file [ " + couplesFile + " ]. It will be removed." << std::endl;
          continue;
        }
        std::string coupleId = couplesIds[srcPhenotype->individualIds[i]];
        if(srcPhenotype->individualIdsIdx.count(coupleId) == 0)
        {
          misc.message << "WARNING: The second individual of the couple (< " + srcPhenotype->individualIds[i] + " >, < " + srcPhenotype->individualIds[i] + " >) is not in the phentypes file. The first individual of the couple will be removed." << std::endl;
          continue;
        }
        newOrderedIndividualsIds.push_back( coupleId );
        individualsKeepReplacementIds.push_back( srcPhenotype->individualIds[i] );
      }
    }
    
    phenotypes[i]->filterIndividuals(newOrderedIndividualsIds);
    phenotypes[i]->replaceIndividualIds(individualsKeepReplacementIds);
  }
}

std::vector<int> getPhenotyesForAnalysis()
{
  std::vector<int> phenotypesForAnalyze;
  if(options.analyzeAllPhenos == false)
  {
    phenotypesForAnalyze = options.phenotypeColumns;
  }
  else
  {
    int nMaxPhenoCol = getNumberOfFileColumns(options.phenotypesFile) - 2;
    for(int i = 0; i<nMaxPhenoCol; i++)
    {
      phenotypesForAnalyze.push_back(i + 1);
    }
  }
  return phenotypesForAnalyze;
}

float ran3(long *idnum)
{
  static int inext, inextp;
  static long ma[56];   /* The value 56 (range ma[1..55]) is special */
  /* and should not be modified; see Knuth.    */
  static int iff=0;
  long mj, mk;
  int i, ii, k;
  
  if (*idnum<0 || iff==0) {    /* Initialization */
    iff = 1;
    /* Initialize ma[55] using the seed idnum and the large number MSEED */
    mj = MSEED - (*idnum<0 ? -*idnum : *idnum);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    /* Now initizalize the rest of the table, in a slightly       */
    /* random order, with numbers that are not especially random. */
    for (i=1; i<=54; i++) {
      ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ) mk += MBIG;
      mj = ma[ii];
    }
    /* We randomize them by "warming up the generator." */
    for (k=1; k<=4; k++) 
      for (i=1; i<=55; i++) {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i]<MZ) ma[i] += MBIG;
      }
      inext = 0;     /* Prepare indices for our first generated number. */
      inextp = 31;   /* The constant 31 is special; see Knuth */
      *idnum = 1;
  }
  /* Here is where we start, except on initialization */
  if (++inext==56) inext = 1;   /* Initizalize inext and inextp, wrapping    */
    if (++inextp==56) inextp = 1; /* arround 56 to 1.                          */
      mj = ma[inext] - ma[inextp];  /* Generate new random number substractively */
      if (mj<MZ) mj +=MBIG;         /* Make sure that it is in range.            */
        ma[inext] = mj;               /* Store it                                  */
        
        
        return mj*FAC;                /* and output the derived uniform deviate.   */
}

double unif_rand_dbl(long *idnum)
{
  double highorder = (double) ran3(idnum);
  double loworder = (double) ran3(idnum);
  
  return highorder + loworder*FAC;
}

double box_muller(double m, double s, long *idnum)     
{                                       
  double x1, x2, w, y1;
  static double y2;
  static int use_last = 0;
  
  int nRepetitions = 0;
  if (use_last)
  {
    y1 = y2;
    use_last = 0;
  }
  else
  {
    do {
      x1 = 2.0 * unif_rand_dbl(idnum) - 1.0;
      x2 = 2.0 * unif_rand_dbl(idnum) - 1.0;
      w = x1 * x1 + x2 * x2;
      
      nRepetitions++;
      if(nRepetitions > 100000)
      {
        misc.error("Error: Sorry, an internal error was happened in the random number generator.", 0);
      }
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }
  return( m + y1 * s );
}

double chi1_CDF(int df, double x)
{
  #ifdef BOOSTLIB
    return boost::math::cdf(boost::math::complement(boost::math::chi_squared(df), x));
  #else
    return -1.;
  #endif
}

double FStatCDF(double df1, double df2, double x)
{
#ifdef BOOSTLIB
  return boost::math::cdf(boost::math::complement(boost::math::fisher_f(df1, df2), x));
#else
  return -1.;
#endif
}

double tStatCDF(double df, double x)
{
#ifdef BOOSTLIB
  return boost::math::cdf(boost::math::complement(boost::math::students_t(df), x));
#else
  return -1.;
#endif
}

std::string getFileName(std::string & path)
{
  return path.substr(path.find_last_of("/\\") + 1);
}


std::string compressData(const std::string & srcstr)
{
  std::stringstream data(srcstr);
  std::stringstream compressed;
  
#if defined(BOOSTLIB) && defined(ZLIB)
  boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
  out.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_compression)));
  out.push(data);
  boost::iostreams::copy(out, compressed);
#else
  misc.error("Error: An internal error was happened. The current version of DISSECT is compiled without zlib and boost support. These libraries are required for writting compressed output. Please, recompile it using these libraries or ask for support.", 0);
#endif
  
  return compressed.str();
}

/////////////////////////////////////////
//Just for gdb and debugging purposes:

std::string& SSS (const char* s)
{
  return *(new std::string(s));
}

int compareMatrices(Matrix * m1, Matrix * m2, double threshold)
{
  if( (m1->nRows != m2->nRows) || (m1->nCols != m2->nCols) )
  {
    misc.error("Error: Comparing matrices with different dimensions.", 0);
    return 0;
  }
  if( (m1->nBlockRows != m2->nBlockRows) || (m1->nBlockCols != m2->nBlockCols) )
  {
    misc.error("Error: Comparing matrices with different structure.", 0);
    return 0;
  }
  
  int result = 1;
  
  #pragma omp parallel for
  for(int c = 0; c<m1->nCols; c++)
  {
    for(int r = 0; r<m1->nRows; r++)
    {
      if(m1->m[c*m1->nRows + r] != m2->m[c*m1->nRows + r] && threshold < 0.)
      {
        result = 0;
        std::cout << r << " " << c << " " << m1->m[c*m1->nRows + r] << " " << m2->m[c*m1->nRows + r] << " " <<  std::endl;
      }
      if( ( fabs(m1->m[c*m1->nRows + r] - m2->m[c*m1->nRows + r])/std::min(fabs(m1->m[c*m1->nRows + r]), fabs(m2->m[c*m1->nRows + r])) ) > threshold && threshold >= 0. )
      {
        result = 0;
        std::cout << r << " " << c << " " << m1->m[c*m1->nRows + r] << " " << m2->m[c*m1->nRows + r] << " " <<  std::endl;
      }
    }
  }
  
  if(result == 0)
  {
    std::cout << "\n*****************************\nMatrices differ\n*****************************\n" << std::endl;
    return 0;
  }
  return 1;
  
//   int * results = NULL;
//   //Gather the sizes of each array to gather. Only on root. On ohter processes results == NULL.
//   results = gather(&result, 1);
//   
//   if(communicator->mpiRoot)
//   {
//     for(int i = 0; i<communicator->mpiNumTasks; i++)
//     {
//       if(results[i] == 0)
//       {
//         result = 0;
//       }
//     }
//     delete [] results;
//   }
//   communicator->broadcast(result);
//   return result;
}

int compareGlobalMatrices(Matrix * m1, Matrix * m2, double threshold)
{
  std::vector< std::vector<double> > gm1;
  std::vector< std::vector<double> > gm2;
  
  m1->matrixToStandardVector(gm1);
  m2->matrixToStandardVector(gm2);
  
  if( (m1->nGlobRows != m2->nGlobRows) || (m1->nGlobCols != m2->nGlobCols) )
  {
    misc.error("Error: Comparing matrices with different dimensions.", 0);
    return 0;
  }
  
  int result = 1;
  
  if(communicator->mpiRoot)
  {
    for(int c = 0; c<m1->nGlobCols; c++)
    {
      for(int r = 0; r<m1->nGlobRows; r++)
      {
        if(gm1[r][c] != gm2[r][c] && threshold < 0.)
        {
          result = 0;
          std::cout << r << " " << c << " " << gm1[r][c] << " " << gm2[r][c] << " " <<  std::endl;
        }
        if( ( fabs(gm1[r][c] - gm2[r][c])/std::min(fabs(gm1[r][c]), fabs(gm2[r][c])) ) > threshold && threshold >= 0. )
        {
          result = 0;
          std::cout << r << " " << c << " " << gm1[r][c] << " " << gm2[r][c] << " " <<  std::endl;
        }
      }
    }
    
    if(result == 0)
    {
      std::cout << "\n\n\n*****************************\nMatrices differ\n*****************************\n\n\n" << std::endl;
      misc.error("ERROR", 0);
      return 0;
    }
  }
  return 1;
}

/*
GRM * loadGRMUsingOptionsOld(bool returnGenotype, Genotype **returnedGenotype)
{
  GRM *grm = NULL;
  
  if(returnGenotype == true)
  {
    if(returnedGenotype == NULL)
    {
      misc.error("Error: An internal error was happened when loading the GRM. Wrong pointer for returning a genotype.", 0);
    }
    *returnedGenotype = NULL;
  }
  
  if( returnGenotype == true && options.GRMJoinMethod == 0 )
  {
    misc.message << "Changing GRM joining method to method 1." << std::endl;
    options.GRMJoinMethod = 1;
  }
  
  if(options.grmFile == "")
  {
    if(options.genotypeListFile == "") //Create GRM from single genotype file
    {
      Genotype *genotype = new Genotype(options.genotypeFile);
      grm = new GRM(genotype);
      if(returnGenotype == false)
      {
        delete genotype;
      }
      else
      {
        *returnedGenotype = genotype;
      }
    }
    else //Create GRM from multiple genotype files
    {
      if(options.GRMJoinMethod == 0) //Create GRM using method 0: First compute GRMs for each genotype, then add GRMs.
      {
        misc.message << "Computing GRM using a list of genotype files. Method 0: First compute GRMs for each genotype, then add GRMs." << std::endl;
        misc.message.tab = "  ";
        misc.setGetElapsedTime("GRMPartComputation");
        std::ifstream file;
        misc.checkFileExists(options.genotypeListFile);
        file.open(options.genotypeListFile.c_str());
        
        int idx = 0;
        std::string line;
        std::vector<std::string> loadedSNPs;
        while(getline(file,line))
        {
          std::istringstream sstemp(line);
          std::string partialFileName;
          sstemp >> partialFileName;
          
          if( idx == 0 )
          {
            Genotype *genotype = new Genotype(line);
            loadedSNPs = genotype->SNPIds;
            
            grm = new GRM(genotype, false);
            if(returnGenotype == false)
            {
              delete genotype;
            }
            else
            {
              misc.error("Error: An internal error was happened. GRM cannot be loaded using method 0 when genotypes must also be loaded.", 0);
            }
          }
          else
          {
            Genotype *genotype = new Genotype(line);
            std::vector<std::string> SNPsIntersection = intersectionStringVectors(2, &loadedSNPs, &genotype->SNPIds); //Check there are not intersections between SNPs in files.
            if(SNPsIntersection.size() != 0)
            {
              misc.error("Error: When computing a GRM from multiple genotype files. There are SNPs repeated in different files.", 0);
            }
            loadedSNPs.insert( loadedSNPs.end(), genotype->SNPIds.begin(), genotype->SNPIds.end() );
            
            GRM *grmToAdd = new GRM(genotype, false);
            delete genotype;
            misc.message << "Adding GRM..." << std::endl;
            grm->addGRMs(1., grmToAdd);
            delete grmToAdd;
          }
          
          idx++;
        }
        file.close();
        misc.message.tab = "";
        grm->normalize();
        misc.message << "GRM computation finished after " << misc.setGetElapsedTime("GRMPartComputation", true) << std::endl;
      }
      else  //Create GRM using method 1: First join genotypes, then compute GRM.
      {
        misc.message << "Computing GRM using a list of genotype files. Method 1: First join genotypes, then compute GRM." << std::endl;
        misc.message.tab = "  ";
        misc.setGetElapsedTime("GRMPartComputation");
        Genotype *genotype = new Genotype();
        genotype->loadList(options.genotypeListFile);
        grm = new GRM(genotype);
        misc.message.tab = "";
        misc.message << "GRM computation finished after " << misc.setGetElapsedTime("GRMPartComputation", true) << std::endl;
        if(returnGenotype == false)
        {
          delete genotype;
        }
        else
        {
          *returnedGenotype = genotype;
        }
      }
    }
  }
  else
  {
    grm = new GRM(options.grmFile);
  }
  return grm;
}
*/
