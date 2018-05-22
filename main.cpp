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

#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

#include <unistd.h>

#include <omp.h>

#include "main.h"
#include "global.h"
#include "matrix.h"
#include "options.h"
#include "communicator.h"
#include "misc.h"
#include "genotype.h"
#include "phenotype.h"
#include "covariate.h"
#include "reml.h"
#include "simulatephenotype.h"
#include "pca.h"
#include "auxiliar.h"
#include "results.h"
#include "analysis.h"

Communicator * communicator;
Options options;
Misc misc;

#include "test.h"

int main(int argc, char **argv)//The covariates has to be checked before running this version. This error to avoid compilation.
{
  int nMaxThreads = omp_get_max_threads();
  
  communicator = new Communicator(argc, argv);
  
  std::stringstream messageBuffer;
  messageBuffer << "@==========================================================@" << std::endl;
  messageBuffer << "+                _______________ ___ ___                   +" << std::endl;
  messageBuffer << "+               |  \\  / _/ _/ __/  _|  _|                  +" << std::endl;
  messageBuffer << "+               ||) )(__ \\__\\ _|  (_  |                    +" << std::endl;
  messageBuffer << "+               |__/__\\__/__/___\\___|_|                    +" << std::endl;
  messageBuffer << "+                                                          +" << std::endl;
  messageBuffer << "@==========================================================@" << std::endl;
  messageBuffer << "|                     version v1.15.2c                     |" << std::endl;
  messageBuffer << "|----------------------------------------------------------|" << std::endl;
  messageBuffer << "|   (C) 2014-2018 Oriol Canela-Xandri and Albert Tenesa    |" << std::endl;
  messageBuffer << "|       The Roslin Institute (University of Edinburgh)     |" << std::endl;
  messageBuffer << "|----------------------------------------------------------|" << std::endl;
  messageBuffer << "|  For documentation, citation or bug reporting see:       |" << std::endl;
  messageBuffer << "|    http://www.dissect.ed.ac.uk/                          |" << std::endl;
  messageBuffer << "@----------------------------------------------------------@" << std::endl;
  messageBuffer << std::endl;
  
  options = Options(argc, argv);
  misc.changeOutputs(std::cout, options.outFile + ".log");
  communicator->nDefaultBlockRows = options.defaultBlockSize;
  communicator->nDefaultBlockCols = options.defaultBlockSize;
  
  misc.message << messageBuffer.str();  
  misc.message << communicator->creationMessage();
  misc.message << "Number of available threads: " << nMaxThreads << std::endl; //Assuming all MPI processes have access to the same number of threads.
  options.showParsedOptions();
  misc.message << options.warnings;
  
  misc.setGetElapsedTime("Total");
  time_t now = time(NULL);
  misc.message << std::endl << "Analysis started: " << ctime(&now);
  
  
  Analysis analysis;
  
  ////////////////////////////////////////
  // Just check the options
  if(options.analysis == justCheckAnalysis)
  {
    misc.message << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    misc.message << "It seems that geneasy arguments can be parsed properly!" << std::endl;
    misc.message << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
  }
  
  ////////////////////////////////////////
  // Show help
  else if(options.analysis == showHelpAnalysis)
  {
    options.showOptions();
  }
  
  ////////////////////////////////////////
  // Compute GRM
  else if(options.analysis == makeGRMAnalysis)
  {
    analysis.makeGRM();
  }
  
  ////////////////////////////////////////
  // Compute GRM and store most related
  else if(options.analysis == makeGRMKeepMostRelated)
  {
    analysis.makeGRMAndStoreMostRelated();
  }
  
  ////////////////////////////////////////
  // REML analysis
  else if(options.analysis == REMLAnalysis)
  {
    analysis.makeREML();
  }
  
  ////////////////////////////////////////
  // bivariate REML analysis
  else if(options.analysis == bivarREMLAnalysis)
  {
    analysis.makeMultivarREML();
  }
  
  ////////////////////////////////////////
  // multivariate REML analysis
  else if(options.analysis == multivarREMLAnalysis)
  {
    analysis.makeMultivarREML();
  }
  
  ////////////////////////////////////////
  // simulate phenotypes
  else if(options.analysis == simulatePhenotypeAnalysis)
  {
    analysis.makeSimulatePhenotype();
  }
  
  ////////////////////////////////////////
  // predict phenotypes
  else if(options.analysis == predictPhenotypeAnalysis)
  {
    analysis.makePredictPhenotype();
  }
  
  ////////////////////////////////////////
  // PCA analysis
  else if(options.analysis == PCAAnalysis)
  {
    analysis.makePCA();
  }
  
  ////////////////////////////////////////
  // GWAS analysis
  else if(options.analysis == GWASAnalysis)
  {
    analysis.makeGWAS();
  }
  
  ////////////////////////////////////////
  // Recursive GWAS analysis
  else if(options.analysis == recursiveGWASAnalysis)
  {
    analysis.makeRecursiveGWAS();
  }
  
  ////////////////////////////////////////
  // Compute the accuracy as a function of SNPs included in the predictor. Experimental analysis.
  else if(options.analysis == accuracyBySNPAnalysis)
  {
    analysis.makeComputeAccuracyBySNP();
  }
  
  ////////////////////////////////////////
  // Perform analysis from gene effects matrices
  else if( options.analysis == effectsAnalysis )
  {
    analysis.makeEffectsAnalysis();
  }
  
  ////////////////////////////////////////
  // GLMM analysis (Unfinished)
  else if(options.analysis == GLMMAnalysis)
  {
    analysis.makeGLMMAnalysis();
  }
  
  else if(options.analysis == iGWASAnalysis)
  {
    analysis.makeIGWASAnalysis();
  }
  
  else if(options.analysis == predictCovarPhenotypeAnalysis )
  {
    analysis.makePredictCovarPhenotype();
  }
  //Compute multiple phenotype gwas on previously precomputed residuals
  else if(options.analysis == multiplePhenotypeGWASAnalysis )
  {
    analysis.makeMultiplePhenotypeGWAS();
  }
  //Compute residuals from multiple phenotypes for performing --mpgwas analysis
  else if(options.analysis == multiplePhenotypeResiduals)
  {
    analysis.makeMultiplePhenotypeResiduals();
  }
  //Filter labeled matrix
  else if(options.analysis == filterLabeledMatrix)
  {
    analysis.makeFilterLabeledMatrix();
  }
  //Get SNP stats (currently allele frequencies)
  else if(options.analysis == computeSNPStatsAnalysis )
  {
    analysis.makeSNPStats();
  }
  
  ////////////////////////////////////////
  // No analysis specified
  else
  {
    misc.message << "\nNot a valid analysis was specified. For a short list of options, run with --help or -h options or visit our web page.\n" << std::endl;
    //options.showOptions();
    //test15();
    //test23();
    //test26();
    //test29();
    //test29bis();
    //test20();
    //test32();
    //test39();
    //test40();
    //test41();
    //test43ARCHER();
    //test6();
    //test44();
    //test45();
    //test45bis();
    //test45bis2();
    //test45bis3();
    //test3();
    //test46();
    //test47();
    //test48();
    //test49();
    //test50();
    //test51();
    //test52();
    //test16();
    //test53();
    //test53bis();
    //test54();
    //test54bis();
    //test29();
    //test29bis2();
    //test28();
    //test55();
    //test55bis();
    //test15();
    //test56();
    //test20();
    //test20bis();
    //test19();
    //test14();
    //test57();
    //test58();
    //test59();
    //test60();
    //test61();
    //test62();
    //test63();
    //test64();
    //test65();
    //test66();
    //test67();
    //test67bis();
    //test68();
    //test9();
    //test69();
    //test70();
    //test71();
    //test72();
    //test73();
    //test74();
    //test75();
    //test76();
    //test77(argc, argv);
    //test78(argc, argv);
    //test79(argc, argv);
    //test80(argc, argv);
    //test81(argc, argv);
    //test82();
    //test83();
    //test84();
    //test85();
    //test86();
    //test87();
    //test88();
    //test89();
    //test90();
    //test91();
    //test92();
    //test93();
    //test94();
    //test95();
    //test96();
    //test97();
    //test98();
    //test99();
    //test100();
    //test101();
    //test102();
    //test103();
    //test104();
    //test105();
    //test106();
    //test107();
    //test108();
    //test108bis();
    //test109();
    //test110();
    //test111();
    //test112();
    //test113();
    //test114();
    //test114bis();
    //testDebug115();
    //test116();
    //test117();
    //test118();
    //test119();
    //test120();
  }
  
  now = time(NULL);
  misc.message << "Analysis finished: " << ctime(&now);
  misc.message << "Total time used for the analysis: " << misc.setGetElapsedTime("Total") << std::endl;
  misc.message << "Estimated memory used for the analysis: " << std::setprecision(2) << misc.maxMemory/(1024.*1024.*1024) << " GB. (" << std::setprecision(2) << misc.maxMemory/(1024.*1024.*1024*double(communicator->mpiNumTasks))  << " GB for each MPI process)" << std::endl;
  misc.message << "Thanks for using me!" << std::endl;
  
  delete communicator;
  
  return 0;
}

