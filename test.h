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

#ifndef TEST_H
#define TEST_H

#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

#include <unistd.h>

#include <omp.h>

#include "main.h"
#include "global.h"
#include "matrix.h"
#include "options.h"
#include "communicator.h"
#include "misc.h"
#include "genotype.h"
#include "kernel.h"
#include "phenotype.h"
#include "covariate.h"
#include "reml.h"
#include "simulatephenotype.h"
#include "pca.h"
#include "auxiliar.h"
#include "covariancematrix.h"
#include "analysis.h"
#include "blockmatrix.h"
#include "singlereml.h"
#include "labeledmatrix.h"
#include "glm.h"
#include "glmm.h"

Matrix * testCreateRandomMatrix(int nr, int nc, int s); //Create a random matrix;

void test1();
void test2();
void test2bis();
void test3();
void test4();
void test5();
void test6(); //Test filtering columns
void test7(); //Test filtering individuals in grm, phenotypes and covariates and genotypes
void test8(); //Test copySubMatrix() method of Matrix class.
void test9(); //Test asymmetric filtering individuals in grm
void test10(); //joinMatrices()
void test11(); //Matrix().invert()
void test12(); //Test define genotype Regions.
void test13(); //Test communicator gather and asymmetricGather
void test14(); //Test grm prunning
void test15(); //Test covariate syncronization
void test16(); //Test matrix eigenvalues and matrix bending/
void test17(); //Test covariance matrix functions
void test18(); //Test joining genotypes.
void test19(); //Test grm addition.
void test20(); //Test grm construction from different sources.
void test20bis(); //Test grm construction from different sources after switching to Kernel class.
void test21(); //Test Matrix write and read
void test22(); //Test loading matrices from different threads.
void test23(); //Test communicator asymmetricGather with doubles
void test24(); //Test parallel genotype loading
void test25(); //Test trace of product of two matrices.
void test26(); //Test new function for filter matrix rows and columns.
void test27(); //Test new matrix write/read functions.
void test28(); //Test matrix pack/unpack functions.
void test29(); //Test MPI read/write
void test29bis(); //Test MPI read/write
void test29bis2(); //Test pack matrices
void test30(); //Test pdtran/pdlaset/pdlacpy change to _
void test31(); //Test pdtran/pdlaset/pdlacpy change to _
void test32(); //Test MPI/noMPI grm outptut
void test33(); //Test symmetricInvert with float/double cases and the determinant computation.
void test34(); //Test matrix transposition.
void test35(); //Test string broadcast.
void test36(); //Test statistical functions;
void test37(); //Test diagonal() Matrix method;
void test38(); //Test covariate NA dealing;
void test39(); //Test box_muller;

void test40(); //Test generalResorting
void test41(); //Test generalResorting
void test42(double *original, int nr, int nc, int bnr, int bnc, std::map<int, int> & rowsOriginDestination, std::map<int, int> & colsOriginDestination); //Test generalResorting2
void test43ARCHER(); //Test generalResorting performance.

void test44(); //Test diagonalDistribution modifications.
void test45(); //Test multiplication after diagonalDistribution modifications.
void test45bis(); //Same as test45 but with matrices of higher dimension.
void test45bis2(); //Same as test45 but with a scaling factor.
void test45bis3(); //Testing symmetric properties after multiplication.
void test46(); //Test element wise division and multiplication.
void test47(); //Test addition with diagonal distributed matrices.
void test48(); //Test inversion with diagonal matrices.
void test49(); //Test trace, diagonal and other stuff with diagonal distributed matrices.
void test50(); //Test diagonal() and transpose()
void test51(); //joinMatrices with diagonal distributed matrices
void test52(); //Test matrixToStandardVector with diagonal distributed matrices.

void test53(); //Test GRM diagonalization.
void test53bis(); //Test GRM diagonalization after switch to Kernel class.
void test54(); //Test GRM diagonalization read/write.
void test54bis(); //Test GRM diagonalization read/write after switch to Kernel class.
void test55(); //Test GRM copying and some matrix checks with diagonal matrices.
void test55bis(); //Test GRM copiing and some matrix checks with diagonal matrices after switch to Kernel class.

void test56(); //Test Phenotype column counting and testing.

void test57(); //Test GRM substitution by Kernel. //No longer, GRM class deleted
void test57aux(Kernel *g1, Kernel *g2, Kernel *k1, Kernel *k2, bool symmetrize = true);

void test58(); //Test genotype grouping by byOrderedFixedSize

void test59(); //Test scatter int vector: scatterVectorRet() method
void test60(); //test makeIntersectionMatrix() method
void test61(); //Test diagonal() and diagonal() when the matrix is not a square matrix.

void test62(); //Test get dependent columns;

void test63(); //Test --extract option

void test64(); //Test getGlobalIndexOutsideRange() method
void test65(); //Test getGlobalIndexInsideRange() method

void test66(); //Test GRM check overlaping SNPs 

void test67(); //Test BlockMatrix;
void test67bis(); //Test BlockMatrix and error with mkl scalapack inversion
void test67bis2(int nBlocks, int nRows, int nCols, int defBlockRow, int defBlockCol); //Test BlockMatrix and error with mkl scalapack inversion

void test68(); //Test BlockMatrix multiplications with normal distributed matrices.
void test68bis(int nBlocks, int nRows, int nCols, int defBlockRow, int defBlockCol, int colDimensionTestM); //Test BlockMatrix multiplications with normal distributed matrices.

void test69(); //Test Epistasis Kernel computation.

void test70(); //Get table from file.

void test71(); //Test prepareRAWREML

void test72(); //Test filter unrelated individuals from Kernel.

void test73(); //Test loading covar effects.

void test74(); //Test creating discrete random effects kernel.

void test75(); //Test parsing covariate restricted model.

void test76(); //Test compute diagonal ABAt.

void test77(int argc, char **argv); //Test for different communicator creation.
void test78(int argc, char **argv); //Test for different communicator multiplication.
void test79(int argc, char **argv); //Test for different communicator REML.
void test80(int argc, char **argv); //Test for different communicator. Matrix redistribution.
void test81(int argc, char **argv); //Test for different communicator. Matrix copy.

void test82();  //Test introduceResortedGRMsByCouples() method in singlereml class.

void test83();  //Test OOM error.

void test84();  //Test covariance matrix with multireml.
void test85();  //Test the second derivatiaves of the covariance matrix.

void test86(); //Test trace of product of two matrices when these are or are not symmetric.

void test87(); //Test multiplication order change.

void test88(); //Test scaleBy

void test89(); //Test makeDifferenceMatrix()

void test90(); //Test quadratic exponential kernel()

void test91(); //Test matrix average
double test91aux(Matrix * m); //Auxiliar test matrix average

void test92(); //Test kernel resorting

void test93(); //Test applyExponentialOperator

void test94(); //Test gradient
double test94aux(REML & reml); //Test gradient

void test95(); //Test LabeledMatrix

void test96(); //Test Matrix method for standardizing matrices.
void test97(); //Test Matrix method for centering matrices.

void test98(); //Test glm.

void test99(); //Test box muller speed.

void test100(); //Test unif_rand_dbl.

void test101(); //Test MH sampling.

void test102(); //Test LabeledMatrix join vertically and horizontally.
LabeledMatrix * test102CreateTestLabeledMatrix(int nr, int nc, std::string rprefix, std::string cprefix, double scale); //Creates a labeled matrix;

void test103(); //Test covariate prefiltering

void test104(); //Test loading GCTA grm

void test105(); //Test random effects with multiple categories.

void test106(); //Test phenotype filtering with --keep.

void test107(); //Test vector broadcast.

void test108(); //Test genotype load with new functions.
void test108bis(); //Test genotype load with new functions.

void test109(); //Test genotype grouping by file ordered windows.

void test110(); //Test BGEN genotype loading.

void test111(); //Test LabeledMatrix loadRaw.

void test112(); //Test BGEN genotype loading 2.

void test113(); //Test BGEN genotype loading and filtering by SNP.
void test114(); //Test BED genotype loading and filtering by SNP and Inds.
void test114bis(); //Test BGEN genotype loading and filtering by SNP and Inds.

void testDebug115();

void test116(); //Test asymetricGather with strings

void test117(); //Test store different strings in different processes.

void test118(); //Test genotype loading and filtering

void test119(); //Test genotype mean center and matrix mean center.

void test120(); //Test labeled matrix filtering.

#endif
