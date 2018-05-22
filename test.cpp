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

#include "test.h"

/**
 * Test6
 */
void test6()
{
  Matrix m = Matrix(cyclicDistribution);
  Matrix mt = Matrix(cyclicDistribution);
  m.debugRead("testMatrices/matrix.dat", 8, 7);
  m.showGlobal("Matrix1");
  
  Matrix mf = Matrix(cyclicDistribution);
  /*int fr[3] = {2, 4, 6};
   * int nfr = 3;
   * int fc[6] = {0, 1, 2, 4, 5, 6};
   * int nfc = 6;
   * m.filterRowsAndColumns(&mf, &fr[0], nfr, &fc[0], nfc, false);
   * mf.showGlobal("Matrix1 filtered");*/
  
  
  
  /*m.debugRead("testMatrices/matrix11.dat", 11, 7);
   * m.showGlobal("Matrix11");
   * int fr[3] = {0, 4, 10};
   * int nfr = 3;
   * int fc[2] = {0, 6};
   * int nfc = 2;
   * m.filterRowsAndColumns(&mf, &fr[0], nfr, &fc[0], nfc, false);
   * mf.showGlobal("Matrix11 filtered");*/
  
  
  m.debugRead("testMatrices/matrix11.dat", 11, 7);
  mt.transpose(&m);
  mt.showGlobal("Matrix11 Transpose");
  int fc[3] = {0, 4, 10};
  int nfc = 3;
  int fr[2] = {0, 6};
  int nfr = 2;
  //mt.filterRowsAndColumns(&mf, &fr[0], nfr, &fc[0], nfc, false);
  mt.filterRowsAndColumns(&mf, &fr[0], nfr, &fc[0], nfc, true);
  mf.showGlobal("Matrix11 Trasnposed filtered");
}

/**
 * Test3
 */
void test3()
{
  std::string fname;
  double *mGlob;
  int ttnr, ttnc;
  
  ttnr = 5;
  ttnc = 3;
  fname = "testMatrices/matrix2.dat";
  
  Matrix m1 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  if (communicator->mpiRoot) {
    mGlob = new double[m1.nGlobRows*m1.nGlobCols];
    std::ifstream file(fname.c_str());
    for (int r = 0; r < m1.nGlobRows; ++r) {
      for (int c = 0; c < m1.nGlobCols; ++c) {
        file >> *(mGlob + m1.nGlobRows*c + r);
      }
    }
    file.close();
  }
  m1.scatterMatrix(mGlob);
  if (communicator->mpiRoot) {
    delete [] mGlob;
  }
  
  ttnr = 5;
  ttnc = 3;
  fname = "testMatrices/matrix5.dat";
  
  
  Matrix m2 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows+1, communicator->nDefaultBlockCols-1);
  if (communicator->mpiRoot) {
    mGlob = new double[m2.nGlobRows*m2.nGlobCols];
    std::ifstream file(fname.c_str());
    for (int r = 0; r < m2.nGlobRows; ++r) {
      for (int c = 0; c < m2.nGlobCols; ++c) {
        file >> *(mGlob + m2.nGlobRows*c + r);
      }
    }
    file.close();
  }
  m2.scatterMatrix(mGlob);
  if (communicator->mpiRoot) {
    delete [] mGlob;
  }
  
  m1.showGlobal("");
  m2.showGlobal("");
  
  //m2.add(&m1);
  //m2.add(&m1, 1., 1.);
  //m2.add(&m1, 2., 2.);
  //m2.add(&m1, 1., 2.);
  m2.add(&m1, 1., 1., subMatrix(1,1,2,2), subMatrix(0,0,2,2));
  //m2.add(&m1, 1., 1., subMatrix(1,1,2,2), subMatrix(1,1,2,2));
  //m2.add(&m1, 1., 1., subMatrix(1,1,5,2), subMatrix(1,1,5,2));
  //m2.add(&m1, 1., 1., subMatrix(1,1,3,2), subMatrix(1,1,4,2));
  m2.showGlobal("");
}

/**
 * Test2
 */
void test2()
{
  std::string fname;
  double *mGlob;
  int ttnr, ttnc;
  
  ttnr = 5;
  ttnc = 3;
  fname = "matrix2.dat";
  
  Matrix m1 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  if (communicator->mpiRoot) {
    mGlob = new double[m1.nGlobRows*m1.nGlobCols];
    std::ifstream file(fname.c_str());
    for (int r = 0; r < m1.nGlobRows; ++r) {
      for (int c = 0; c < m1.nGlobCols; ++c) {
        file >> *(mGlob + m1.nGlobRows*c + r);
      }
    }
    file.close();
  }
  m1.scatterMatrix(mGlob);
  if (communicator->mpiRoot) {
    delete [] mGlob;
  }
  Matrix m1copy = Matrix(&m1);
  
  //ttnr = 5;
  //ttnc = 3;
  //fname = "matrix3.dat";
  ttnr = 11;
  ttnc = 8;
  fname = "matrix9.dat";
  
  Matrix m2 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  if (communicator->mpiRoot) {
    mGlob = new double[m2.nGlobRows*m2.nGlobCols];
    std::ifstream file(fname.c_str());
    for (int r = 0; r < m2.nGlobRows; ++r) {
      for (int c = 0; c < m2.nGlobCols; ++c) {
        file >> *(mGlob + m2.nGlobRows*c + r);
      }
    }
    file.close();
  }
  m2.scatterMatrix(mGlob);
  if (communicator->mpiRoot) {
    delete [] mGlob;
  }
  Matrix m2copy = Matrix(&m2);
  
  Matrix *r = new Matrix(cyclicDistribution);
  Matrix *r2 = new Matrix(cyclicDistribution);
  Matrix *r3 = new Matrix(cyclicDistribution);
  
  m1.showGlobal("m1");
  m2.showGlobal("m2");
  
  /*r->multiply(&m1, 'n', &m2, 't');
   * r->showGlobal("m1*m2");
   * 
   * r->multiply(&m1, 'N', &m1, 'T');
   * r->showGlobal("m1*m1T");*/
  
  
  r->multiply(&m1, 't', &m1, 'n');
  r->showGlobal("m1T*m1");
  
  double determinant;
  r->multiply(&m2, 'T', &m2, 'N');
  r->showGlobal("m2T*m2");
  r->symmetricInvert(&determinant);
  std::cout << "Determinant: " << std::setprecision(9) << determinant << std::endl; std::cout.flush();
  r->showGlobal("(m2T*m2)^-1");
  //Relative error of 0.005% with respect python numpy in determinant calculation. Is this acceptable?
  
  r3->multiply(&m2copy, 'T', &m2, 'N');
  
  //r->showGlobal("test_1");
  //r3->showGlobal("m2T*m2");
  
  r2->multiply(r, 'N', r3, 'N');
  r2->showGlobal("(m2T*m2)^-1*(m2T*m2)");
  
  r2->multiply(r3, 'N', r, 'N');
  r2->showGlobal("(m2T*m2)*(m2T*m2)^-1");
  
  r->showGlobal("symmetric L", false);
  r2->transpose(r);
  r2->showGlobal("symmetric U", false);
  r3->transpose(r2);
  r3->showGlobal("symmetric L", false);
  r->symmetrizeTriangularMatrix();
  r->showGlobal("symmetric B", false);
  r3->transpose(r);
  r3->showGlobal("symmetric B", false);
  
}

/**
 * Test1
 */
void test1()
{
  //int ttnr = 8;
  //int ttnc = 7;
  //std::string fname = "matrix.dat";
  int ttnr = 5;
  int ttnc = 3;
  std::string fname = "matrix2.dat";
  
  Matrix m = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  //std::cout << "P1" << std::endl;
  //std::cout.flush();
  
  double *testvr;
  double *testvc;
  if (communicator->mpiRoot) {
    testvr = new double [ttnr];
    testvc = new double [ttnc];
    for(int i=0; i<ttnr; i++)
    {
      testvr[i] = i+1;
    }
    for(int i=0; i<ttnc; i++)
    {
      testvc[i] = i+1;
    }
  }
  
  double *mGlob, *mGlob2;
  if (communicator->mpiRoot) {
    mGlob = new double[m.nGlobRows*m.nGlobCols];
    mGlob2 = new double[m.nGlobRows*m.nGlobCols];
    std::ifstream file(fname.c_str());
    for (int r = 0; r < m.nGlobRows; ++r) {
      for (int c = 0; c < m.nGlobCols; ++c) {
        file >> *(mGlob + m.nGlobRows*c + r);
      }
    }
    
    std::cout << "Matrix A:\n";
    for (int r = 0; r < m.nGlobRows; ++r) {
      for (int c = 0; c < m.nGlobCols; ++c) {
        std::cout << std::setw(3) << *(mGlob + m.nGlobRows*c + r) << " ";
      }
      std::cout << "\n";
    }
    std::cout << std::endl;
  }
  
  //m.scatterMatrix(mGlob);
  double test[communicator->nDefaultBlockRows*communicator->nDefaultBlockCols];
  for (int r = 0; r < m.nGlobRows; r += m.nBlockRows) {
    for (int c = 0; c < m.nGlobCols; c += m.nBlockCols) {
      if(communicator->mpiRoot)
      {
        for(int i=0; i<m.nBlockRows; i++)
        {
          for(int j=0; j<m.nBlockCols; j++)
          {
            test[j*m.nBlockRows + i] = *(mGlob+(m.nGlobRows*(j+c))+(i+r));
          }
        }
      }
      m.scatterBlock(test, r, c, m.nBlockRows);
    }
  }
  
  Cblacs_barrier(communicator->context, "All");
  
  ////////////////////////////////////
  // Write Matrix to a file
  std::ofstream file;
  if(communicator->mpiRoot)
  {
    std::string fname = "test.bin";
    file.open(fname.c_str(), std::ofstream::out | std::ofstream::binary);
  }
  
  m.writeMatrixFile(file);
  
  if(communicator->mpiRoot)
  {
    file.close();
  }
  
  
  ////////////////////////////////////
  // Distribute vector between processses  
  
  misc.message << "Begin partial vectors test ***********************************" << std::endl;
  m.scatterVector(testvr, row);
  m.showPartial(row, true);
  
  m.scatterVector(testvc, column);
  m.showPartial(column, true);
  misc.message << "End partial vectors test ***********************************" << std::endl;
  
  
  m.gatherMatrix(mGlob2);
  Cblacs_barrier(communicator->context, "All");
  
  if (communicator->mpiRoot) {
    std::cout << "Matrix A2:\n";
    for (int r = 0; r < m.nGlobRows; ++r) {
      for (int c = 0; c < m.nGlobCols; ++c) {
        std::cout << std::setw(3) << *(mGlob2 + m.nGlobRows*c + r) << " ";
      }
      std::cout << "\n";
    }
    std::cout << std::endl;
  }
  
  
  Matrix mmt = Matrix(cyclicDistribution);
  mmt.multiply(&m, 'T', &m, 'N');
  
  double *mGlob3;
  if (communicator->mpiRoot) {
    //Reserve space and read matrix (with transposition!)
    mGlob3 = new double [mmt.nGlobRows*mmt.nGlobCols];
  }
  
  //std::cout << mmt.m[0] << std::endl;
  
  mmt.gatherMatrix(mGlob3);
  if (communicator->mpiRoot) {
    std::cout << "Matrix A*AT:\n";
    for (int r = 0; r < mmt.nGlobRows; ++r) {
      for (int c = 0; c < mmt.nGlobCols; ++c) {
        std::cout << std::setw(3) << *(mGlob3 + mmt.nGlobRows*c + r) << " ";
      }
      std::cout << "\n";
    }
    std::cout << std::endl;
  }
}

/**
 * Test2bis
 */
void test2bis()
{
  std::string fname;
  double *mGlob;
  int ttnr, ttnc;
  
  ttnr = 2;
  ttnc = 2;
  fname = "matrix6.dat";
  
  Matrix m1 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  if (communicator->mpiRoot) {
    mGlob = new double[m1.nGlobRows*m1.nGlobCols];
    std::ifstream file(fname.c_str());
    for (int r = 0; r < m1.nGlobRows; ++r) {
      for (int c = 0; c < m1.nGlobCols; ++c) {
        file >> *(mGlob + m1.nGlobRows*c + r);
      }
    }
    file.close();
  }
  m1.scatterMatrix(mGlob);
  if (communicator->mpiRoot) {
    delete [] mGlob;
  }
  Matrix m1copy = Matrix(&m1);
  
  ttnr = 3;
  ttnc = 2;
  fname = "matrix7.dat";
  Matrix m2 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  if (communicator->mpiRoot) {
    mGlob = new double[m2.nGlobRows*m2.nGlobCols];
    std::ifstream file(fname.c_str());
    for (int r = 0; r < m2.nGlobRows; ++r) {
      for (int c = 0; c < m2.nGlobCols; ++c) {
        file >> *(mGlob + m2.nGlobRows*c + r);
      }
    }
    file.close();
  }
  m2.scatterMatrix(mGlob);
  if (communicator->mpiRoot) {
    delete [] mGlob;
  }
  Matrix m2copy = Matrix(&m2);
  
  ttnr = 2;
  ttnc = 3;
  fname = "matrix8.dat";
  Matrix m3 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  if (communicator->mpiRoot) {
    mGlob = new double[m3.nGlobRows*m3.nGlobCols];
    std::ifstream file(fname.c_str());
    for (int r = 0; r < m3.nGlobRows; ++r) {
      for (int c = 0; c < m3.nGlobCols; ++c) {
        file >> *(mGlob + m3.nGlobRows*c + r);
      }
    }
    file.close();
  }
  m3.scatterMatrix(mGlob);
  if (communicator->mpiRoot) {
    delete [] mGlob;
  }
  Matrix m3copy = Matrix(&m3);
  
  Matrix *r = new Matrix(cyclicDistribution);
  
  m1.showGlobal("m1");
  m2.showGlobal("m2");
  m3.showGlobal("m3");
  
  m1.symmetric = true;
  m1.uplo = 'L';
  
  r->multiply(&m1, 'N', &m3, 'N');
  r->showGlobal("test1");
  
  r->multiply(&m2, 'N', &m1, 'N');
  r->showGlobal("test1");
  
  r->multiply(&m2, 'N', &m2, 'T');
  r->showGlobal("test1");
  
}

/**
 * Test4
 */
void test4()
{
  int ttnr = 8;
  int ttnc = 7;
  std::string fname = "matrix.dat";
  
  Matrix m = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  
  double *mGlob;
  if (communicator->mpiRoot) {
    mGlob = new double[m.nGlobRows*m.nGlobCols];
    std::ifstream file(fname.c_str());
    for (int r = 0; r < m.nGlobRows; ++r) {
      for (int c = 0; c < m.nGlobCols; ++c) {
        file >> *(mGlob + m.nGlobRows*c + r);
      }
    }
  }
  
  m.scatterMatrix(mGlob);
  
  m.showPartial(row, false, true);
  m.showPartial(row, false, false);
  
  misc.message << "***************************" << std::endl;
  Matrix m2(cyclicDistribution);
  m2.debugRead("matrix11.dat", 11, 7);
  m2.showPartial(row, false, false);
  
}

/**
 * Test5
 */
void test5()
{
  Genotype *genotype = new Genotype();
  genotype->load(options.genotypeFile);
  Kernel *grm = new Kernel(genotype);
  delete genotype;
  
  int dimension = grm->kernel->nGlobCols;
  
  Phenotype *temp = new Phenotype(cyclicDistribution, options.phenotypesFile, 1);
  Matrix *y = temp->phenotypes;
  Matrix *identity = new Matrix(cyclicDistribution, dimension, dimension);
  Matrix *identity2 = new Matrix(cyclicDistribution, dimension, dimension);
  
  identity->fillDiagonal(2., 0.);
  //identity->symmetric = false;
  //identity2->fillDiagonal(2., 1.);
  //identity2->symmetric = false;
  
  Matrix * subVPy = new Matrix(cyclicDistribution, dimension, 2, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  communicator->barrier();
  
  y->showGlobal("y");
  identity->showGlobal("I");
  grm->kernel->showGlobal("GRM");
  
  subVPy->fillWithConstant(0.);
  communicator->barrier();
  subVPy->multiply(grm->kernel, 'N', y, 'N', 1., subMatrix(0, 0, dimension, 1));
  subVPy->showGlobal("subVPy");
  communicator->barrier();
  subVPy->multiply(identity, 'N', y, 'N', 1., subMatrix(0, 1, dimension, 1));
  subVPy->showGlobal("subVPy");
  communicator->barrier();
  
  delete temp;
  delete grm;
}

void test7()
{
  Covariate *covariate = new Covariate(options.covarsFile, options.qCovarsFile, std::vector<std::string>());
  covariate->printCovariate();
  
  std::vector<std::string> test;
  if(communicator->mpiRoot)
  {
    test.push_back("i2@i2");
    test.push_back("i3@i3");
    test.push_back("i6@i6");
  }
  
  covariate->filterIndividuals(test);
  covariate->printCovariate();
  delete covariate;
  
  Phenotype *phenotype = new Phenotype(cyclicDistribution, options.phenotypesFile, 2);
  phenotype->printPhenotype();
  phenotype->filterIndividuals(test);
  phenotype->printPhenotype();
  delete phenotype;
  
  Genotype *genotype = new Genotype();
  std::cout << options.genotypeFile << std::endl; std::cout.flush();
  genotype->load(options.genotypeFile);
  
  Kernel grm = Kernel(genotype);
  delete genotype;
  grm.printGRM();
  grm.filterIndividuals(test, true);
  std::cout << "***************************************" << std::endl;
  std::cout << "Pointers: " << grm.kernel << " " << grm.N << std::endl;
  std::cout << "***************************************" << std::endl;
  grm.printGRM();
  std::vector<std::string> test2;
  test2.push_back("i3@i3");
  grm.filterIndividuals(test2, true);
  std::cout << "***************************************" << std::endl;
  std::cout << "Pointers: " << grm.kernel << " " << grm.N << std::endl;
  std::cout << "***************************************" << std::endl;
  grm.printGRM();
  
  std::vector<std::string> SNPs;
  std::vector<std::string> SNPs2;
  if(communicator->mpiRoot)
  {
    SNPs.push_back("rs2");
    SNPs.push_back("rs5");
    SNPs.push_back("rs6");
    SNPs.push_back("rs7");
    SNPs2.push_back("rs2");
    SNPs2.push_back("rs5");
    SNPs2.push_back("rs7");
  }
  test.push_back("i7@i7");
  genotype = new Genotype(options.genotypeFile);
  genotype->printGenotype();
  /*if(communicator->mpiRoot)
   * {
   *   SNPs.clear();
   *   test.clear();
   *   std::cout << "Unfiltered" << std::endl;
}*/
  std::cout << "Filtered" << std::endl;
  genotype->filterSNPsAndIndividuals(SNPs, test, true);
  genotype->printGenotype();
  test.pop_back();
  genotype->filterSNPsAndIndividuals(SNPs2, test, true);
  genotype->printGenotype();
  delete genotype;
  
  misc.message << "\n*******************\nFiltering out of place\n*******************\n" << std::endl;
  
  genotype = new Genotype(options.genotypeFile);
  Genotype * genotype2 = new Genotype();
  genotype->filterSNPsAndIndividuals(SNPs, test, true, genotype2);
  genotype->printGenotype();
  genotype2->printGenotype();
  genotype->printGenotype();
  genotype->filterSNPsAndIndividuals(SNPs, test, true);
  genotype->printGenotype();
}

void test8()
{
  Matrix m1(cyclicDistribution);
  m1.debugRead("matrix11.dat", 11, 7);
  Matrix m2(cyclicDistribution);
  m2.debugRead("matrix12.dat", 11, 7);
  
  m1.showGlobal("Original");
  
  int dr, dc;
  dr = 0; dc = 0;
  m1.copySubMatrix(&m2, subMatrix(0, 0, m2.nGlobRows - dr, m2.nGlobCols - dc), subMatrix(0, 0, m2.nGlobRows - dr, m2.nGlobCols - dc));
  m1.showGlobal("Matriu");
  m1.debugRead("matrix11.dat", 11, 7);
  
  //   m2.multiply(&m1, 'N', &m1, 'T');
  //   m2.symmetrizeTriangularMatrix();
  //   m2.showGlobal("symmetrized", false);
  
  /*if(communicator->mpiRank == 1){
    int i = 0;
    std::cout << "PID " << getpid() << " " << communicator->mpiRank << " on " << communicator->hostName << " ready for attach\n" << std::endl;
    fflush(stdout);
    while (0 == i)
    {
    }
  }*/
  communicator->barrier();
  
  dr = 2; dc = 2;
  m1.copySubMatrix(&m2, subMatrix(1, 1, m2.nGlobRows - dr,m2.nGlobCols - dc), subMatrix(2, 2, m2.nGlobRows - dr, m2.nGlobCols - dc));
  m1.showGlobal("Matriu");
  //m1.debugRead("matrix11.dat", 11, 7);
  
  //   dr = 2; dc = 2;
  //   m1.copySubMatrix(&m2, subMatrix(2, 2, m2.nGlobRows - dr, m2.nGlobCols - dc), subMatrix(0, 0, m2.nGlobRows - dr, m2.nGlobCols - dc));
  //   m1.showGlobal("Matriu");
  //   m1.debugRead("matrix11.dat", 11, 7);
  //   
  //   dr = 3; dc = 3;
  //   m1.copySubMatrix(&m2, subMatrix(3, 3, m2.nGlobRows - dr, m2.nGlobCols - dc), subMatrix(0, 0, m2.nGlobRows - dr, m2.nGlobCols - dc));
  //   m1.showGlobal("Matriu");
  //   m1.debugRead("matrix11.dat", 11, 7);
}

void test9()
{
  
  std::vector<std::string> test;
  std::vector<std::string> test2;
  if(communicator->mpiRoot)
  {
    test.push_back("i6@i6");
    test.push_back("i3@i3");
    test.push_back("i2@i2");
    test.push_back("i7@i7");
    
    test2.push_back("i7@i7");
    test2.push_back("i3@i3");
    test2.push_back("i4@i4");
  }
  
  /****** Expected result ******************************* 
      i7           i3          i4
      i6 i6 i6@i6 5:   -0.637531     0.0900183   0.433621
      i3 i3 i3@i3 2:   -0.523132     0.547161   -0.643315
      i2 i2 i2@i2 1:      1.2975    -0.535943   -0.255979
      i7 i7 i7@i7 6:      1.2975    -0.523132   -0.548642



      i6 i6 i6@i6 5: 6 6 7
      i3 i3 i3@i3 2: 5 6 6
      i2 i2 i2@i2 1: 6 6 7
      i7 i7 i7@i7 6: 6 5 6
   */
  
  Genotype *genotype = new Genotype();
  std::cout << options.genotypeFile << std::endl; std::cout.flush();
  genotype->load(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;

  /*
  grm->diagonalized = true;
  grm->filterIndividualsAsymmetric(grm->individualIds, test2, true); //This raises an error, as expected.
  //grm->filterIndividualsAsymmetric(test, grm->individualIds, true); //This raises an error, as expected.
  //grm->filterIndividualsAsymmetric(grm->individualIds, grm->individualIds, true); //This works fine since no filtering is needed and dos not matter if it is diagonalized.
  grm->diagonalized = false;
  */
  
  grm->printGRM();
  grm->filterIndividualsAsymmetric(test, test2, true);
  std::cout << "***************************************" << std::endl;
  std::cout << "Pointers: " << grm->kernel << " " << grm->N << std::endl;
  std::cout << "***************************************" << std::endl;
  grm->printGRM();
  
  delete grm;
  
  genotype = new Genotype();
  std::cout << options.genotypeFile << std::endl; std::cout.flush();
  genotype->load(options.genotypeFile);
  grm = new Kernel(genotype);
  delete genotype;
  
  //grm->kernel->showGlobal("Kernel", false);
  //grm->N->showGlobal("N", false);
  //grm->N->uplo='U'; //For testing exceptions
  //grm->N = NULL;    //For testing exceptions
  grm->filterIndividualsAsymmetric(test2, test, true);
  std::cout << "***************************************" << std::endl;
  std::cout << "Pointers: " << grm->kernel << " " << grm->N << std::endl;
  std::cout << "***************************************" << std::endl;
  grm->printGRM();
  //grm->kernel->showGlobal("Kernel", false);
  //std::cout << "--------------------N=" << grm->N << std::endl;
  //grm->N->showGlobal("N", false);
  delete grm;
  
  
  genotype = new Genotype();
  //std::cout << options.genotypeFile << std::endl; std::cout.flush();
  genotype->load(options.genotypeFile);
  grm = new Kernel(genotype);
  delete genotype;
  
  grm->printGRM();
  Matrix testKernel(grm->kernel);
  Matrix testN(grm->N);
  grm->filterIndividualsAsymmetric(grm->individualIds, grm->individualIds, true);
  std::cout << "***************************************" << std::endl;
  std::cout << "Pointers: " << grm->kernel << " " << grm->N << std::endl;
  std::cout << "***************************************" << std::endl;
  grm->printGRM();
  compareMatrices(&testKernel, grm->kernel);
  compareMatrices(&testN, grm->N);
  delete grm;
  
  /////////////////////////////////////////////////
  //
  //  Test normal filtering
  //
  
  /****** Expected result ******************************* 
    i5 i5 i5@i5 4:      0.3809   0.0471612   -0.318904  -0.447054
    i3 i3 i3@i3 2:   0.0471612   0.547161   -0.535943  -0.523132
    i2 i2 i2@i2 1:   -0.318904  -0.535943     1.17643     1.2975
    i7 i7 i7@i7 6:   -0.447054  -0.523132      1.2975     1.2975


    i5 i5 i5@i5 4: 7 6 7 6 
    i3 i3 i3@i3 2: 6 6 6 5 
    i2 i2 i2@i2 1: 7 6 7 6 
    i7 i7 i7@i7 6: 6 5 6 6 
   */
  
  test.clear();
  if(communicator->mpiRoot)
  {
    test.push_back("i5@i5");
    test.push_back("i3@i3");
    test.push_back("i2@i2");
    test.push_back("i7@i7");
  }
  
  genotype = new Genotype();
  //std::cout << options.genotypeFile << std::endl; std::cout.flush();
  genotype->load(options.genotypeFile);
  grm = new Kernel(genotype);
  Kernel * grm2 = new Kernel(grm);
  delete genotype;
  
  misc.message << "\n\n\n\nNormal filtering test\n\n\n\n" << std::endl;
  grm->printGRM();
  grm->filterIndividualsAsymmetric(test, test, true);
  grm2->filterIndividuals(test, true);
  std::cout << "***************************************" << std::endl;
  std::cout << "Pointers: " << grm->kernel << " " << grm->N << std::endl;
  std::cout << "***************************************" << std::endl;
  grm->printGRM();
  grm2->printGRM();
  
  Matrix testKernel1(grm->kernel);
  Matrix testN1(grm->N);
  Matrix testKernel2(grm2->kernel);
  Matrix testN2(grm2->N);
  compareMatrices(&testKernel1, &testKernel2);
  compareMatrices(&testN1, &testN2);
  delete grm;
  delete grm2;
}

void test10()
{
  Matrix m1(cyclicDistribution, 3, 4, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m1.fillWithConstant(1.1);
  Matrix m2(cyclicDistribution, 5, 2, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m2.fillWithConstant(2.2);
  
  Matrix m3(cyclicDistribution, 3, 2, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m3.fillWithConstant(-4);
  
  Matrix r(cyclicDistribution);
  
  
  m1.showGlobal("m1");
  m2.showGlobal("m2");
  m3.showGlobal("m3");
  r.joinMatrices(&m1, subMatrix(1,1,m1.nGlobRows,m1.nGlobCols), &m2, subMatrix(3+2,4+2,m2.nGlobRows,m2.nGlobCols), 3.3);
  r.showGlobal("r");
  
  r.joinMatricesVertically(&m2, &m3);
  r.showGlobal("r");
  r.joinMatricesHorizontally(&m1, &m3);
  r.showGlobal("r");
}

/**
 * Test11
 */
void test11()
{
  std::string fname;
  double *mGlob;
  int ttnr, ttnc;
  
  Matrix m1(cyclicDistribution);
  m1.debugRead("matrix2.dat", 5, 3);
  
  Matrix m1copy(&m1);
  
  Matrix m2(cyclicDistribution);
  m2.debugRead("matrix9.dat", 11, 8);
  Matrix m2copy(&m2);
  
  Genotype *genotype = new Genotype();
  std::cout << options.genotypeFile << std::endl; std::cout.flush();
  genotype->load(options.genotypeFile);
  genotype->printGenotype();
  Kernel grm = Kernel(genotype);
  delete genotype;
  
  Matrix *r = new Matrix(cyclicDistribution);
  Matrix *r2 = new Matrix(cyclicDistribution);
  Matrix *r3 = new Matrix(cyclicDistribution);
  
  std::cout << "LCM: " << leastCommonMultiple(communicator->nProcRows, communicator->nProcCols) << std::endl;
  
  double determinant;
  //r->multiply(&m2, 'T', &m2, 'N');
  //r3->multiply(&m2copy, 'T', &m2, 'N');
  r = new Matrix(grm.kernel);
  r3 = new Matrix(r);
  r->showGlobal("m2T*m2");
  r->invert(&determinant);
  //r->symmetricInvert(&determinant);
  std::cout << "Determinant: " << std::setprecision(9) << determinant << std::endl; std::cout.flush();
  r->showGlobal("(m2T*m2)^-1");
  
  
  r2->multiply(r, 'N', r3, 'N');
  r2->showGlobal("(m2T*m2)^-1*(m2T*m2)");
  
  r2->multiply(r3, 'N', r, 'N');
  r2->showGlobal("(m2T*m2)*(m2T*m2)^-1");
  
}

void test12()
{
  //testd with args:
  //--regions test/test3.groups --region-size 1
  //--regions test/test3.groups
  //--region-size 1
  //--region-size 100 --region-overlap 20
  
  Genotype *genotype = new Genotype(options.genotypeFile);
  
  if(options.regionsFile != "" && options.regionalAnalysis)
  {
    genotype->groupSNPs(byGene);
  }
  else if(options.regionalAnalysis)
  {
    genotype->groupSNPs(byPosition);
  }
  
  genotype->printGenotype(true);
  
  delete genotype;
}

void test13()
{
  int v1[] = {1+10*(communicator->mpiRank + 1), 2+10*(communicator->mpiRank + 1), 3+10*(communicator->mpiRank + 1)};
  /*for(int i = 0; i<3; i++)
   * {
   *   std::cout << v1[i] << std::endl;
}*/
  
  
  int * result;
  
  result = communicator->gather(v1, 3);
  std::cout << "=== * === " << result << std::endl;
  if(communicator->mpiRoot)
  {
    for(int i = 0; i<3*communicator->mpiNumTasks; i++)
    {
      std::cout << result[i] << std::endl;
    }
  }
  delete [] result;
  
  /*int nTotalValues;
   * result = communicator->asymmetricGather(v1, 3, &nTotalValues);
   * if(communicator->mpiRoot)
   * {
   *   for(int i = 0; i<nTotalValues; i++)
   *   {
   *     std::cout << result[i] << std::endl;
}
}
delete [] result;*/
  
  int *v2;
  int nv2;
  if(communicator->mpiRank == 0)
  {
    nv2 = 2;
    v2 = new int [nv2];
    v2[0] = 11; v2[1] = 12;
  }
  else if(communicator->mpiRank == 1)
  {
    nv2 = 7;
    v2 = new int [nv2];
    v2[0] = 21; v2[1] = 22; v2[2] = 23; v2[3] = 24; v2[4] = 25; v2[5] = 26; v2[6] = 27;
  }
  else if(communicator->mpiRank == 2)
  {
    nv2 = 3;
    v2 = new int [nv2];
    v2[0] = 31; v2[1] = 32; v2[2] = 33;
  }
  else
  {
    nv2 = 5;
    v2 = new int [nv2];
    v2[0] = -1; v2[1] = -1; v2[2] = -1; v2[3] = -1; v2[4] = -1;
  }
  int nTotalValues;
  result = communicator->asymmetricGather(v2, nv2, &nTotalValues);
  std::cout << "*** " << nTotalValues << " *** " << result << std::endl;
  if(communicator->mpiRoot)
  {
    for(int i = 0; i<nTotalValues; i++)
    {
      std::cout << result[i] << std::endl;
    }
  }
  delete [] result;
  delete [] v2;
}

void test14()
{
  //Tested with
  
  //--grm-cutoff 0.4
  //--grm-cutoff 0.5
  //--grm-cutoff 0.025
  
  //with cutoff = 0.025 it should be kept individuals: i1,i2,i4 or i4,i7,i8
  //with cutoff = 0.5 it should be filtered individuals: i2,i3 or i7,i8
  
  misc.message << "GRM cutoff: " << options.grmCutoff << std::endl;
  
  communicator->nDefaultBlockCols = 3;
  communicator->nDefaultBlockRows = 2;
  
  Genotype * genotype = new Genotype(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;
  
  std::vector<double> temp;
  
  Kernel * grmUpper = new Kernel(grm);
  grmUpper->kernel->transpose(grm->kernel);
  grmUpper->N->transpose(grm->N);
  grmUpper->kernel->matrixToStandardVector(temp);
  if(communicator->mpiRoot == true)  {temp[1] = 2.;}
  grmUpper->kernel->scatterMatrix(&(temp[0]));
  grmUpper->kernel->showGlobal("grmUpper", false);
  std::vector<std::string> individualIdsUpper = grmUpper->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  Kernel * grmLower = new Kernel(grm);
  grmLower->kernel->matrixToStandardVector(temp);
  if(communicator->mpiRoot == true)  {temp[69] = 2.;}
  grmLower->kernel->scatterMatrix(&(temp[0]));
  grmLower->kernel->showGlobal("grmLower", false);
  std::vector<std::string> individualIdsLower = grmLower->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  Kernel * grmBoth = new Kernel(grm);
  grmBoth->kernel->symmetrizeTriangularMatrix();
  grmBoth->N->symmetrizeTriangularMatrix();
  grmBoth->kernel->showGlobal("grmBoth", false);
  std::vector<std::string> individualIdsBoth = grmBoth->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  //grmUpper->filterIndividuals(individualIdsUpper);
  grmUpper->pruneKernel(options.grmCutoff);
  grmUpper->kernel->showGlobal("grmUpper", false);
  grmUpper->printGRM();
  std::vector< std::vector<double> > gUpper;
  grmUpper->kernel->matrixToStandardVector(gUpper);
  
  //grmLower->filterIndividuals(individualIdsLower);
  grmLower->pruneKernel(options.grmCutoff);
  grmLower->kernel->showGlobal("grmLower", false);
  grmLower->printGRM();
  std::vector< std::vector<double> > gLower;
  grmLower->kernel->matrixToStandardVector(gLower);
  
  //grmBoth->filterIndividuals(individualIdsBoth);
  grmBoth->pruneKernel(options.grmCutoff);
  grmBoth->printGRM();
  std::vector< std::vector<double> > gBoth;
  grmBoth->kernel->matrixToStandardVector(gBoth);
  
  if( individualIdsUpper.size() != individualIdsBoth.size() || individualIdsLower.size() != individualIdsBoth.size() || individualIdsUpper.size() != individualIdsLower.size() )
  {
    misc.error("Error: Discordant number individuals prunned.", 0);
  }
  misc.message << "**************** " << individualIdsBoth.size() << " "  << individualIdsUpper.size() << " " << individualIdsLower.size() << std::endl;
  for(int i = 0; i< individualIdsBoth.size(); i++)
  {
    misc.message << individualIdsBoth[i] << " "  << individualIdsUpper[i] << " " << individualIdsLower[i] << std::endl;
  }
  
  if(communicator->mpiRoot == true)
  {
    if(gBoth.size() != gLower.size() || gBoth.size() != gUpper.size() || gBoth.size() != individualIdsBoth.size())
    {
      misc.error("Error: Something is wrong.", 0);
    }
    for(int i = 0; i< individualIdsBoth.size(); i++)
    {
      for(int j = 0; j< individualIdsBoth.size(); j++)
      {
        if( (gBoth[i][j] > options.grmCutoff && i!=j) || (gLower[i][j] > options.grmCutoff && j<i) || (gUpper[i][j] > options.grmCutoff && i<j))
        {
          misc.error("Error: Not properly pruned.", 0);
        }
      }
    }
  }
  
  grm->printGRM();
  std::vector<std::string> individualIds = grm->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  misc.message << "**************** " << individualIds.size() << std::endl;
  for(int i = 0; i< individualIds.size(); i++)
  {
    misc.message << individualIds[i] << std::endl;
  }
  
  grm->filterIndividuals(individualIds);
  grm->printGRM();
  
  delete grm;
  
  
  //Test with a larger GRM
  
  options.grmFile = "";
  options.genotypeFile = "test/parts.snps/merged";
  options.genotypeListFile = "";
  options.GRMJoinMethod = 0;
  Kernel* grmLarge = loadGRMUsingOptions();
  grmLarge->pruneKernel(options.grmCutoff);
  std::vector< std::vector<double> > gLarge;
  grmLarge->kernel->matrixToStandardVector(gLarge);
  
  //grmLarge->kernel->showGlobal("grmLarge", false);
  
  if(communicator->mpiRoot == true)
  {
    if(gLarge.size() != grmLarge->individualIds.size() || grmLarge->kernel->uplo != 'L' )
    {
      misc.error("Error: Something is wrong.", 0);
    }
    for(int i = 0; i< grmLarge->individualIds.size(); i++)
    {
      for(int j = 0; j< grmLarge->individualIds.size(); j++)
      {
        if( (gLarge[i][j] > options.grmCutoff && j<i /*j!=i*/) )
        //if( (gLarge[i][j] > options.grmCutoff && /*j<i*/ j!=i) )
        {
          std::cout << grmLarge->kernel->uplo << " " << i << " " << j << " " << gLarge[i][j] << std::endl;
          misc.error("Error: Not properly pruned.", 0);
        }
      }
    }
  }
}

void test15()
{
  Genotype *genotype = new Genotype(options.genotypeFile);
  std::vector<std::string> individualIds = genotype->individualIds;
  delete genotype;
  
  Covariate *covariate = new Covariate(options.covarsFile, options.qCovarsFile, individualIds);
  covariate->printCovariate(6);
  delete covariate;
  
  Covariate *covariate1 = new Covariate(options.covarsFiles[0], options.qCovarsFiles[0], individualIds, false);
  Covariate *covariate2 = new Covariate(options.covarsFiles[1], options.qCovarsFiles[1], individualIds, false);
  
  covariate1->syncronizeDiscreteCovariateCategoriesWith(covariate2);
  
  covariate1->parseRawCovariates(individualIds, 2, 0);
  covariate2->parseRawCovariates(individualIds, 2, 1);
  
  covariate1->printCovariate(6);
  covariate2->printCovariate(6);
  
  if(covariate1->covarCategories == covariate2->covarCategories)
  {
    misc.message << "\n**************\nCovariate categories are concordant\n**************\n" << std::endl;
  }
  
  delete covariate1;
  delete covariate2;
  
  //   covariate1 = new Covariate(options.covarsFiles[0], options.qCovarsFiles[0], false);
  //   //covariate1->printCovariate();
  //   covariate2 = new Covariate(options.covarsFiles[1], options.qCovarsFiles[1], false);
  //   //covariate2->printCovariate();
  //   
  //   covariate2->syncronizeDiscreteCovariateCategoriesWith(covariate1);
  //   
  //   if(covariate1->covarCategories == covariate2->covarCategories)
  //   {
  //     misc.message << "\n**************\nCovariate categories are concordant\n**************\n" << std::endl;
  //   }
  //   
  //   delete covariate1;
  //   delete covariate2;
}

void test16()
{
  Genotype * genotype = new Genotype(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;
  
  Matrix * test = new Matrix(grm->kernel);
  
  test->showGlobal();
  
  Matrix * eigenValues = new Matrix();
  Matrix * eigenVectors = new Matrix();
  
  //test->eigenDecomposition(eigenValues, eigenVectors);
  
  //eigenValues->showGlobal("EigenValues");
  //eigenVectors->showGlobal("EigenVectors");
  
  test->bendMatrix();
  test->showGlobal("Bended");
  
  delete eigenValues;
  delete eigenVectors;
  delete test;
  delete grm;
}

void test17()
{
  //Commented due to a resorting of some REML functions.
  /*
  Analysis analysis;
  //CovarianceMatrix V(10);
  REML reml;
  GRM *grm = loadGRMUsingOptionsOld();
  reml.prepareSingleBivarREML(grm);
  reml.V->appendVarianceToElement("GRM1", "Var(GRM2)", variance);
  reml.V->appendVarianceToElement("GRM1", "Var(E1)", variance);
  reml.V->variances[3].variance = 0.1;
  reml.V->variances[4].variance = 0.2;
  options.fixedCorrelation = 0.5;
  
  reml.V->showCovarianceMatrix(false);
  
  misc.message << std::setw(10) << ".";
  for(int j = 0; j < reml.V->elements.size(); j++)
  {
    misc.message << std::setw(10) << reml.V->elements[j].name << " ";
  }
  misc.message << std::endl;
  for(int i = 0; i < reml.V->variances.size(); i++)
  {
    misc.message << std::setw(10) << reml.V->variances[i].name << " ";
    for(int j = 0; j < reml.V->elements.size(); j++)
    {
      if(reml.V->variances[i].inElements.find(j) != reml.V->variances[i].inElements.end())
      {
        misc.message << std::setw(10) << reml.V->computeElementVarianceDerivate(reml.V->elements[j], i) << " ";
      }
      else
      {
        misc.message << std::setw(10) << "*" << " ";
      }
    }
    misc.message << std::endl;
  }
  misc.message << std::endl;
  for(int i = 0; i < reml.V->elements.size(); i++)
  {
    misc.message << std::setw(10) << reml.V->elements[i].name << " ";
    misc.message << std::setw(10) << reml.V->computeElementVariance(reml.V->elements[i]) << std::endl;
  }
  
  
  reml.fixCovarianceMatrixVariances();
  reml.V->showCovarianceMatrix();
  
  misc.message << std::setw(10) << ".";
  for(int j = 0; j < reml.V->elements.size(); j++)
  {
    misc.message << std::setw(10) << reml.V->elements[j].name << " ";
  }
  misc.message << std::endl;
  for(int i = 0; i < reml.V->variances.size(); i++)
  {
    misc.message << std::setw(10) << reml.V->variances[i].name << " ";
    for(int j = 0; j < reml.V->elements.size(); j++)
    {
      if(reml.V->variances[i].inElements.find(j) != reml.V->variances[i].inElements.end())
      {
        misc.message << std::setw(10) << reml.V->computeElementVarianceDerivate(reml.V->elements[j], i) << " ";
      }
      else
      {
        misc.message << std::setw(10) << "*" << " ";
      }
    }
    misc.message << std::endl;
  }
  misc.message << std::endl;
  for(int i = 0; i < reml.V->elements.size(); i++)
  {
    misc.message << std::setw(10) << reml.V->elements[i].name << " ";
    misc.message << std::setw(10) << reml.V->computeElementVariance(reml.V->elements[i]) << std::endl;
  }
  
  reml.V->computeCovariance();
  reml.V->computeDerivateCovariance(0);
  reml.V->m->showGlobal("Covar");
  reml.V->mDerivate->showGlobal("DerivateCovar");
  */
}

void test18()
{
  options.genotypeListFile = "test/parts.snps/list.dat";
  
  Genotype * genotype1 = new Genotype();
  genotype1->loadList(options.genotypeListFile);
  
  Genotype * genotype2 = new Genotype("test/parts.snps/merged");
  Genotype * genotype3;
  
  compareMatrices(genotype1->genotypes, genotype2->genotypes);
  compareMatrices(genotype1->missings, genotype2->missings);
  
  //genotype1->SNPs[2].frequencies[1] = 790;
  //genotype1->individuals[2].familyID = "ehem";
  //if(communicator->mpiRoot) { std::cout << genotype1->SNPIdsIdx["rs422285"] << std::endl; genotype1->SNPIdsIdx["rs422285"] = 183; }
  if (
    genotype1->nSNPs != genotype2->nSNPs ||
    genotype1->SNPs != genotype2->SNPs ||
    genotype1->SNPIds != genotype2->SNPIds ||
    genotype1->SNPIdsIdx != genotype2->SNPIdsIdx ||
    
    genotype1->nIndividuals != genotype2->nIndividuals ||
    genotype1->individuals != genotype2->individuals ||
    genotype1->individualIds != genotype2->individualIds ||
    genotype1->individualIdsIdx != genotype2->individualIdsIdx
  )
  {
    misc.error("Error: Error when joining genotypes.", 0);
  }
  
  
  genotype1->normalizeGenotypes();
  genotype2->normalizeGenotypes();
  compareMatrices(genotype1->genotypes, genotype2->genotypes);
  
  
  
//Tested changing an individual (value and order) and repeating a SNP
//   options.genotypeListFile = "test/parts.snps/list.fake.dat";
//   genotype3 = new Genotype();
//   genotype3->loadList(options.genotypeListFile);
//   delete genotype3;

//Tested changing genotype value
//   delete genotype2;
//   genotype2 = new Genotype("test/parts.snps/merged");  
//   options.genotypeListFile = "test/parts.snps/list.fake2.dat";
//   genotype3 = new Genotype();
//   genotype3->loadList(options.genotypeListFile);
//   compareMatrices(genotype3->genotypes, genotype2->genotypes);
//   delete genotype3;
  
  delete genotype1;
  delete genotype2;
  
  
  options.genotypeListFile = "test/parts.individuals/list.dat";
  
  genotype1 = new Genotype();
  genotype1->loadList(options.genotypeListFile);
  
  genotype2 = new Genotype("test/parts.individuals/merged");
  
  compareMatrices(genotype1->genotypes, genotype2->genotypes);
  compareMatrices(genotype1->missings, genotype2->missings);
  
  //genotype1->SNPs[2].frequencies[1] = 790;
  //genotype1->individuals[2].familyID = "ehem";
  //if(communicator->mpiRoot)  genotype1->individualIdsIdx["c1_g0_0@c1_g0_0"] = 0;
  //if(communicator->mpiRoot)  genotype1->individualIdsIdx["c1_g0_20@c1_g0_20"] = 9;
  if (
    genotype1->nSNPs != genotype2->nSNPs ||
    genotype1->SNPs != genotype2->SNPs ||
    genotype1->SNPIds != genotype2->SNPIds ||
    genotype1->SNPIdsIdx != genotype2->SNPIdsIdx ||
    
    genotype1->nIndividuals != genotype2->nIndividuals ||
    genotype1->individuals != genotype2->individuals ||
    genotype1->individualIds != genotype2->individualIds ||
    genotype1->individualIdsIdx != genotype2->individualIdsIdx
  )
  {
    misc.error("Error: Error when joining genotypes.", 0);
  }
  
  genotype1->normalizeGenotypes();
  genotype2->normalizeGenotypes();
  compareMatrices(genotype1->genotypes, genotype2->genotypes);

//Tested changing a SNP (name and order) and repeating an individual
//   options.genotypeListFile = "test/parts.individuals/list.fake.dat";
//   genotype3 = new Genotype();
//   genotype3->loadList(options.genotypeListFile);
//   delete genotype3;

//Tested changing genotype value 
//     delete genotype2;
//     genotype2 = new Genotype("test/parts.individuals/merged");
//     options.genotypeListFile = "test/parts.individuals/list.fake2.dat";
//     genotype3 = new Genotype();
//     genotype3->loadList(options.genotypeListFile);
//     compareMatrices(genotype3->genotypes, genotype2->genotypes);
//     delete genotype3;
  
  delete genotype1;
  delete genotype2;
}

void test19()
{
  Genotype *genotype = new Genotype(options.genotypeFile);
  Kernel *grm1 = new Kernel(genotype);
  Kernel *grm2 = new Kernel(genotype);
  delete genotype;
  
  grm2->kernel->fillDiagonal(1.12, 0.37);
  grm2->kernel->uplo = 'L';
  Matrix *temp = new Matrix(grm2->N);
  temp->fillDiagonal(1., 3.);
  grm2->N->add(temp);
  delete temp;
  
  
  grm1->kernel->showGlobal("grm1");
  grm1->N->showGlobal("N1");
  grm2->kernel->showGlobal("grm2", false);
  grm2->N->showGlobal("N2");
  
  misc.message << "grm1:" << std::endl;
  for(int i = 0; i<grm1->randomVarNames.size(); i++)
  {
    misc.message << grm1->randomVarNames[i] << std::endl;
  }
  misc.message << "grm2:" << std::endl;
  for(int i = 0; i<grm2->randomVarNames.size(); i++)
  {
    misc.message << grm2->randomVarNames[i] << std::endl;
  }
  
  Kernel *grm3;
  Kernel *grm4;
  
  ////////////////////////////////////////////////////////////
  
  grm3 = new Kernel();
  if(communicator->mpiRoot)
  {
    grm2->randomVarNames.clear();
    grm2->randomVarNames.push_back("rs6N");
    grm2->randomVarNames.push_back("rs3N");
    //grm2->randomVarNames.push_back("rs3");
  }
  //grm2->randomVarNames.push_back("rs3");
  //grm3->addGRMs(0.5, grm1, 1., grm2);
  //grm3->addGRMs(-1., grm1, 0.1, grm2);
  //grm3->addGRMs(-1., grm1, -0.1, grm2);
  //grm3->addGRMs(-1.1, grm1, -1., grm2);
  //grm3->addGRMs(-1., grm1, 1., grm2);
  grm3->addGRMs(1., grm1, 1., grm2);
  grm3->kernel->showGlobal("g1 + g2");
  grm3->N->showGlobal("N1 + N2");
  misc.message << "grm3:" << std::endl;
  for(int i = 0; i<grm3->randomVarNames.size(); i++)
  {
    misc.message << grm3->randomVarNames[i] << std::endl;
  }
  grm4 = new Kernel(grm2);
  grm4->addGRMs(1., grm1);
  grm4->kernel->symmetrizeTriangularMatrix();
  grm4->N->symmetrizeTriangularMatrix();
  grm3->kernel->symmetrizeTriangularMatrix();
  grm3->N->symmetrizeTriangularMatrix();
  compareMatrices(grm3->kernel, grm4->kernel, 1e-10); //El normalize, denormalize changes the last digit.
  compareMatrices(grm3->N, grm4->N);
  if(grm4->randomVarNames != grm3->randomVarNames)
  {
    misc.error("Error: Discordant SNPs.", 0);
  }
  grm4 = new Kernel(grm1);
  grm4->addGRMs(1., grm2);
  grm4->kernel->symmetrizeTriangularMatrix();
  grm4->N->symmetrizeTriangularMatrix();
  grm3->kernel->symmetrizeTriangularMatrix();
  grm3->N->symmetrizeTriangularMatrix();
  compareMatrices(grm3->kernel, grm4->kernel, 1e-10);
  compareMatrices(grm3->N, grm4->N);
  
  ////////////////////////////////////////////////////////////
  
  grm3 = new Kernel();
  if(communicator->mpiRoot)
  {
    grm2->randomVarNames.clear();
    grm2->randomVarNames.push_back("rs6");
    grm2->randomVarNames.push_back("rs3");
    //grm2->randomVarNames.push_back("rs3N");
  }
  
  grm3->addGRMs(1., grm1, -1., grm2);
  grm3->kernel->showGlobal("g1 - g2");
  grm3->N->showGlobal("N1 - N2");
  misc.message << "grm3:" << std::endl;
  for(int i = 0; i<grm3->randomVarNames.size(); i++)
  {
    misc.message << grm3->randomVarNames[i] << std::endl;
  }
  grm4 = new Kernel(grm2);
  grm4->addGRMs(1., grm1, -1.);
  grm4->kernel->symmetrizeTriangularMatrix();
  grm4->N->symmetrizeTriangularMatrix();
  grm3->kernel->symmetrizeTriangularMatrix();
  grm3->N->symmetrizeTriangularMatrix();
  compareMatrices(grm3->kernel, grm4->kernel, 1e-10); //El normalize, denormalize changes the last digit.
  compareMatrices(grm3->N, grm4->N);
  if(grm4->randomVarNames != grm3->randomVarNames)
  {
    misc.error("Error: Discordant SNPs.", 0);
  }
  grm4 = new Kernel(grm1);
  grm4->addGRMs(-1., grm2);
  grm4->kernel->symmetrizeTriangularMatrix();
  grm4->N->symmetrizeTriangularMatrix();
  grm3->kernel->symmetrizeTriangularMatrix();
  grm3->N->symmetrizeTriangularMatrix();
  compareMatrices(grm3->kernel, grm4->kernel, 1e-10);
  compareMatrices(grm3->N, grm4->N);
  
  ////////////////////////////////////////////////////////////
  
  grm3 = new Kernel();
  if(communicator->mpiRoot)
  {
    grm2->randomVarNames.clear();
    grm2->randomVarNames.push_back("rs6");
    grm2->randomVarNames.push_back("rs3"); //If I comment this line, an error is raised as expected. It works properly.
    grm2->randomVarNames.push_back("rs7");
    grm2->randomVarNames.push_back("rs5");
    grm2->randomVarNames.push_back("rs4");
    grm2->randomVarNames.push_back("rs2");
    grm2->randomVarNames.push_back("rs1");
    grm2->randomVarNames.push_back("rs1N");
    grm2->randomVarNames.push_back("rs2N");
  }
  
  grm3->addGRMs(-1., grm1, 1., grm2);
  grm3->kernel->showGlobal("- g1 + g2");
  grm3->N->showGlobal("- N1 + N2");
  misc.message << "grm3:" << std::endl;
  for(int i = 0; i<grm3->randomVarNames.size(); i++)
  {
    misc.message << grm3->randomVarNames[i] << std::endl;
  }
  grm4 = new Kernel(grm2);
  grm4->addGRMs(-1., grm1);
  grm4->kernel->symmetrizeTriangularMatrix();
  grm4->N->symmetrizeTriangularMatrix();
  grm3->kernel->symmetrizeTriangularMatrix();
  grm3->N->symmetrizeTriangularMatrix();
  compareMatrices(grm3->kernel, grm4->kernel, 1e-10); //El normalize, denormalize changes the last digit.
  compareMatrices(grm3->N, grm4->N);
  if(grm4->randomVarNames != grm3->randomVarNames)
  {
    misc.error("Error: Discordant SNPs.", 0);
  }
  grm4 = new Kernel(grm1);
  grm4->addGRMs(1., grm2, -1.);
  grm4->kernel->symmetrizeTriangularMatrix();
  grm4->N->symmetrizeTriangularMatrix();
  grm3->kernel->symmetrizeTriangularMatrix();
  grm3->N->symmetrizeTriangularMatrix();
  compareMatrices(grm3->kernel, grm4->kernel, 1e-10);
  compareMatrices(grm3->N, grm4->N);
  
  ////////////////////////////////////////////////////////////
  
  grm3 = new Kernel();
  if(communicator->mpiRoot)
  {
    grm2->randomVarNames.clear();
    grm2->randomVarNames.push_back("rs6N");
    grm2->randomVarNames.push_back("rs3N");
    //grm2->randomVarNames.push_back("rs3");
  }
  
  grm3->addGRMs(-1., grm1, -1., grm2);
  grm3->kernel->showGlobal("- g1 - g2");
  grm3->N->showGlobal("- N1 - N2");
  misc.message << "grm3:" << std::endl;
  for(int i = 0; i<grm3->randomVarNames.size(); i++)
  {
    misc.message << grm3->randomVarNames[i] << std::endl;
  }
  grm4 = new Kernel(grm2);
  grm4->addGRMs(-1., grm1, -1.);
  grm4->kernel->symmetrizeTriangularMatrix();
  grm4->N->symmetrizeTriangularMatrix();
  grm3->kernel->symmetrizeTriangularMatrix();
  grm3->N->symmetrizeTriangularMatrix();
  compareMatrices(grm3->kernel, grm4->kernel, 1e-10); //El normalize, denormalize changes the last digit.
  compareMatrices(grm3->N, grm4->N);
  if(grm4->randomVarNames != grm3->randomVarNames)
  {
    misc.error("Error: Discordant SNPs.", 0);
  }
  grm4 = new Kernel(grm1);
  grm4->addGRMs(-1., grm2, -1.);
  grm4->kernel->symmetrizeTriangularMatrix();
  grm4->N->symmetrizeTriangularMatrix();
  grm3->kernel->symmetrizeTriangularMatrix();
  grm3->N->symmetrizeTriangularMatrix();
  compareMatrices(grm3->kernel, grm4->kernel, 1e-10);
  compareMatrices(grm3->N, grm4->N);
  /*if(grm4->randomVarNames != grm3->randomVarNames)
  {
    misc.message << "\n\n**" << std::endl;
    for(int i = 0; i<grm3->randomVarNames.size(); i++)
    {
      misc.message << grm3->randomVarNames[i] << std::endl;
    }
    misc.message << "**" << std::endl;
    for(int i = 0; i<grm4->randomVarNames.size(); i++)
    {
      misc.message << grm4->randomVarNames[i] << std::endl;
    }
    misc.message << "**\n\n" << std::endl;
    misc.error("Error: Discordant SNPs.", 0); //Estava testant aixo qe ha de donar error.
  }*/
  
  ////////////////////////////////////////////////////////////
  
  misc.message << "\n********************************************\n" << std::endl;
  /* No longer available factors which are not 1. or -1.
  GRM *grm3 = new GRM();
  grm3->addGRMs(0.5, grm1, 1.7, grm2);
  grm3->grm->showGlobal("grm3");
  grm3->N->showGlobal("N3");
  misc.message << "grm3:" << std::endl;
  for(int i = 0; i<grm3->SNPIds.size(); i++)
  {
    misc.message << grm3->SNPIds[i] << std::endl;
  }
  
  grm1->addGRMs(1.7, grm2, 0.5);
  grm1->grm->showGlobal("grm1r2");
  grm1->N->showGlobal("N1r2");
  misc.message << "grm1:" << std::endl;
  for(int i = 0; i<grm1->SNPIds.size(); i++)
  {
    misc.message << grm1->SNPIds[i] << std::endl;
  }

  //grm1->grm->symmetrizeTriangularMatrix();
  //grm3->grm->symmetrizeTriangularMatrix();
  //compareMatrices(grm1->grm, grm3->grm);
  
  grm1->N->fillDiagonal(0.23, 7.2);
  grm1->N->uplo = 'L';
  grm1->N->showGlobal("N1");
  grm1->addGRMs(0.5, grm2);
  grm1->grm->showGlobal("grm1r3");
  grm1->N->showGlobal("N1r3");
  
  grm1->addGRMs(0.3, grm3);
  grm1->grm->showGlobal("grm1r4", true, 20);
  grm1->N->showGlobal("N1r4");
  
  GRM *grm4 = new GRM();
  grm4->addGRMs(0.25, grm1, 0.81, grm2);
  grm4->grm->showGlobal("grm4", false, 20);
  grm4->N->showGlobal("N4");
  
  grm1->grm->showGlobal("grm1r4", true, 20);
  grm1->N->showGlobal("N1r4", true, 20);
  grm2->grm->showGlobal("grm2", true, 10);
  grm2->N->showGlobal("N2", true, 10);
  double t1 = 0.25;
  double t2 = 0.81;
  misc.message << t1 << " " << t2 << " " << ((0.25*0.41723*21.230000000000000426 + 0.81*1.1200000000000001066*7.)/(21.230000000000000426+7.)) << std::endl;
  
  GRM *grm5 = new GRM();
  grm5->addGRMs(0.25, grm4, -0.5, grm2);
  grm5->grm->showGlobal("grm5", true, 20);
  grm5->N->showGlobal("N5");
  */
}

void test20()
{
  Analysis analysis;
  
  options.outFile = "test.test20";
  options.grmFile = options.outFile;
  options.genotypeFile = "test/parts.snps/merged";
  options.genotypeListFile = "test/parts.snps/list.dat";
  options.GRMJoinMethod = 0;
  
  options.grmFile = "";
  options.genotypeFile = "test/parts.snps/merged";
  options.genotypeListFile = "";
  options.GRMJoinMethod = 0;
  Kernel* grm1 = loadGRMUsingOptions();
  grm1->kernel->symmetrizeTriangularMatrix();
  grm1->N->symmetrizeTriangularMatrix();
  grm1->writeGRM(options.outFile);
  
  options.grmFile = "";
  options.genotypeFile = "";
  options.genotypeListFile = "test/parts.snps/list.dat";
  options.GRMJoinMethod = 0;
  Kernel* grm2 = loadGRMUsingOptions();
  grm2->kernel->symmetrizeTriangularMatrix();
  grm2->N->symmetrizeTriangularMatrix();
  
  options.grmFile = "";
  options.genotypeFile = "";
  options.genotypeListFile = "test/parts.snps/list.dat";
  options.GRMJoinMethod = 1;
  Kernel* grm3 = loadGRMUsingOptions();
  grm3->kernel->symmetrizeTriangularMatrix();
  grm3->N->symmetrizeTriangularMatrix();
  
  
  options.grmFile = options.outFile;
  options.genotypeFile = "";
  options.genotypeListFile = "";
  options.GRMJoinMethod = 0;
  Kernel* grm4 = loadGRMUsingOptions();
  grm4->kernel->symmetrizeTriangularMatrix();
  grm4->N->symmetrizeTriangularMatrix();
  
  
  compareGlobalMatrices(grm1->kernel, grm2->kernel, 1e-10);
  compareGlobalMatrices(grm1->N, grm2->N);
  
  compareGlobalMatrices(grm1->kernel, grm3->kernel, 1e-14);
  compareGlobalMatrices(grm1->N, grm3->N);
  
  compareGlobalMatrices(grm1->kernel, grm4->kernel);
  compareGlobalMatrices(grm1->N, grm4->N);
  
  //grm1->randomVarNames.push_back("error");
  if(grm1->randomVarNames != grm2->randomVarNames || grm1->randomVarNames != grm3->randomVarNames)
  {
    misc.error("Error: SNP lists differ!", 0);
  }
  
//   for(int i = 0; i<grm1->randomVarNames.size(); i++)
//   {
//     misc.message << grm1->randomVarNames[i] << std::endl;
//   }
  
  //Kernel* grm4 = loadGRMUsingOptions();
  
  //Kernel* grm5 = loadGRMUsingOptions();
}

void test20bis()
{
  Analysis analysis;
  
  options.outFile = "test.test20";
  options.grmFile = options.outFile;
  options.genotypeFile = "test/parts.snps/merged";
  options.genotypeListFile = "test/parts.snps/list.dat";
  options.GRMJoinMethod = 0;
  
  options.grmFile = "";
  options.genotypeFile = "test/parts.snps/merged";
  options.genotypeListFile = "";
  options.GRMJoinMethod = 0;
  Kernel* grm1 = loadGRMUsingOptions();
  grm1->kernel->symmetrizeTriangularMatrix();
  grm1->N->symmetrizeTriangularMatrix();
  grm1->writeKernel(options.outFile);
  
  options.grmFile = "";
  options.genotypeFile = "";
  options.genotypeListFile = "test/parts.snps/list.dat";
  options.GRMJoinMethod = 0;
  Kernel* grm2 = loadGRMUsingOptions();
  grm2->kernel->symmetrizeTriangularMatrix();
  grm2->N->symmetrizeTriangularMatrix();
  
  options.grmFile = "";
  options.genotypeFile = "";
  options.genotypeListFile = "test/parts.snps/list.dat";
  options.GRMJoinMethod = 1;
  Kernel* grm3 = loadGRMUsingOptions();
  grm3->kernel->symmetrizeTriangularMatrix();
  grm3->N->symmetrizeTriangularMatrix();
  
  
  options.grmFile = options.outFile;
  options.genotypeFile = "";
  options.genotypeListFile = "";
  options.GRMJoinMethod = 0;
  Kernel* grm4 = loadGRMUsingOptions();
  grm4->kernel->symmetrizeTriangularMatrix();
  grm4->N->symmetrizeTriangularMatrix();
  
  
  compareGlobalMatrices(grm1->kernel, grm2->kernel, 1e-10);
  compareGlobalMatrices(grm1->N, grm2->N);
  
  compareGlobalMatrices(grm1->kernel, grm3->kernel, 1e-14);
  compareGlobalMatrices(grm1->N, grm3->N);
  
  compareGlobalMatrices(grm1->kernel, grm4->kernel);
  compareGlobalMatrices(grm1->N, grm4->N);
  
  //grm1->randomVarNames.push_back("error");
  if(grm1->randomVarNames != grm2->randomVarNames || grm1->randomVarNames != grm3->randomVarNames)
  {
    misc.error("Error: SNP lists differ!", 0);
  }
  
//   for(int i = 0; i<grm1->randomVarNames.size(); i++)
//   {
//     misc.message << grm1->randomVarNames[i] << std::endl;
//   }
  
  //Kernel* grm4 = loadGRMUsingOptionsOld();
  
  //Kernel* grm5 = loadGRMUsingOptionsOld();
}

void test21()
{
  Genotype *genotype = new Genotype(options.genotypeFile);
  Kernel *grm = new Kernel(genotype);
  delete genotype;
  Matrix m(grm->kernel);
  delete grm;
  
  m.fillDiagonal(2, 0.1);
  m.symmetric = true;
  m.uplo = 'U';
  m.vector = true;
  
  std::ofstream file;
  file.open("tmp/test.test21.bin", std::ofstream::out | std::ofstream::binary);
  m.writeMatrixFile(file);
  file.close();
  
  Matrix * mread = new Matrix();
  std::cout << mread->symmetric << " " << mread->uplo << " " << mread->vector << std::endl;
  
  std::ifstream file2;
  file2.open("tmp/test.test21.bin", std::ifstream::in | std::ifstream::binary);
  mread->readMatrixFile(file2, m.nGlobRows, m.nGlobCols, m.nBlockRows, m.nBlockCols);
  file2.close();
  
//   if(communicator->mpiRoot || communicator->mpiRank==2)
//   {
//     mread->symmetric = true;
//     mread->uplo = 'X';
//   }
  
  mread->showGlobal("mread");
  mread->showPartial(row, false);
}

void test22()
{
  //Tested with block sizes (2, 2), (3, 5), (5, 3), (2, 3), (9, 11) 
  Matrix m00(cyclicDistribution);
  m00.debugRead("matrix11.dat", 11, 7, 0, 0);
  
  Matrix m01(cyclicDistribution);
  m01.debugRead("matrix11.dat", 11, 7, 0, 1);
  
  Matrix m21(cyclicDistribution);
  m21.debugRead("matrix11.dat", 11, 7, 2, 1);
  
  m00.showGlobal();
  m01.showGlobal();
  m21.showGlobal();
  
  //m00.showPartial(row); //Tested with block sizes = (2, 2)
  //m01.showPartial(row);
  //m21.showPartial(row);
  
  compareMatrices(&m00, &m01);
  compareMatrices(&m00, &m21);
  
  Matrix m00v2(cyclicDistribution);
  m00v2.debugRead("matrix9.dat", 11, 8, 0, 0);
  Matrix m10v2(cyclicDistribution);
  m10v2.debugRead("matrix9.dat", 11, 8, 1, 0);
  Matrix m21v2(cyclicDistribution);
  m21v2.debugRead("matrix9.dat", 11, 8, 2, 1);
  
  //m00v2.m[0] = 0.1; //Just as a null test.
  
  compareMatrices(&m00v2, &m10v2);
  compareMatrices(&m00v2, &m21v2);

}

void test23()
{
  double v1[] = {1+10*(communicator->mpiRank + 1) + 0.1, 2+10*(communicator->mpiRank + 1) + 0.2, 3+10*(communicator->mpiRank + 1) + 0.3};
  /*for(int i = 0; i<3; i++)
  {
      std::cout << v1[i] << std::endl;
  }*/
  
  
  double * results;
  
  results = communicator->gather(v1, 3);
  std::cout << "=== * === " << results << std::endl;
  if(communicator->mpiRoot)
  {
    for(int i = 0; i<3*communicator->mpiNumTasks; i++)
    {
      std::cout << results[i] << std::endl;
    }
  }
  delete [] results;
  
  
  
  double *v2;
  int nv2;
  if(communicator->mpiRank == 0)
  {
    nv2 = 2;
    v2 = new double [nv2];
    v2[0] = 1.1; v2[1] = 1.2;
  }
  else if(communicator->mpiRank == 1)
  {
    nv2 = 7;
    v2 = new double [nv2];
    v2[0] = 2.1; v2[1] = 2.2; v2[2] = 2.3; v2[3] = 2.4; v2[4] = 2.5; v2[5] = 2.6; v2[6] = 2.7;
  }
  else if(communicator->mpiRank == 2)
  {
    nv2 = 3;
    v2 = new double [nv2];
    v2[0] = 3.1; v2[1] = 3.2; v2[2] = 3.3;
  }
  else if(communicator->mpiRank == 3)
  {
    nv2 = 0;
    v2 = new double [nv2];
  }
  else if(communicator->mpiRank == 4)
  {
    nv2 = 1;
    v2 = new double [nv2];
    v2[0] = 5.1;
  }
  else
  {
    nv2 = 5;
    v2 = new double [nv2];
    v2[0] = -(communicator->mpiRank + 1.1); v2[1] = -(communicator->mpiRank + 1.1); v2[2] = -(communicator->mpiRank + 1.1); v2[3] = -(communicator->mpiRank + 1.1); v2[4] = -(communicator->mpiRank + 1.1);
  }
  int nTotalValues;
  double *result;
  result = communicator->asymmetricGather(v2, nv2, &nTotalValues);
  std::cout << "*** " << nTotalValues << " *** " << result << std::endl;
  if(communicator->mpiRoot)
  {
    for(int i = 0; i<nTotalValues; i++)
    {
      std::cout << result[i] << std::endl;
    }
  }
  delete [] result;
  delete [] v2;
  
  
  
  std::vector<double> v22;
  if(communicator->mpiRank == 0)
  {
    v22.push_back(1.1);
    v22.push_back(1.2);
  }
  else if(communicator->mpiRank == 1)
  {
    v22.push_back(2.1);
    v22.push_back(2.2);
    v22.push_back(2.3);
    v22.push_back(2.4);
  }
  else if(communicator->mpiRank == 2)
  {
    v22.push_back(3.1);
    v22.push_back(3.2);
    v22.push_back(3.3);
  }
  else if(communicator->mpiRank == 3)
  {
  }
  else if(communicator->mpiRank == 4)
  {
    v22.push_back(5.1);
  }
  else
  {
    v22.push_back(-(communicator->mpiRank + 1.1));
  }
  int nTotalValues2;
  double *result2;
  result2 = communicator->asymmetricGather(&(v22[0]), v22.size(), &nTotalValues2);
  std::cout << "*** " << nTotalValues2 << " *** " << result2 << std::endl;
  if(communicator->mpiRoot)
  {
    for(int i = 0; i<nTotalValues2; i++)
    {
      std::cout << result2[i] << std::endl;
    }
  }
  delete [] result2;
}

void test24()
{
  Genotype *genotype = new Genotype(options.genotypeFile);
  //Kernel *grm = new Kernel(genotype);
  
  if(communicator->mpiRoot)
  {
    for(int i = 0; i < genotype->SNPs.size(); i++)
    {
      std::cout << std::setw(8) << genotype->SNPs[i].name;
      std::cout << std::setw(8) << genotype->SNPs[i].nNonMissing;
      std::cout << std::setw(8) << genotype->SNPs[i].frequencies[0];
      std::cout << std::setw(8) << genotype->SNPs[i].frequencies[1];
      std::cout << std::setw(8) << genotype->SNPs[i].frequencies[2];
      std::cout << std::setw(8) << genotype->SNPs[i].frequencies[3];
      std::cout << std::setw(10) << genotype->SNPs[i].p1;
      std::cout << std::setw(10) << genotype->SNPs[i].p2;
      std::cout << std::setw(10) << genotype->SNPs[i].standardDev << std::endl;
    }
  }
  
  genotype->printGenotype();
  genotype->genotypes->showGlobal("");
  
  delete genotype;
  //delete grm;
}

void test25()
{
  Matrix * m1 = new Matrix();
  m1->debugRead("matrix13.dat", 11, 9);
  Matrix * m2 = new Matrix();
  m2->multiply(m1, 'T', m1, 'N');
  
  Genotype *genotype = new Genotype(options.genotypeFile);
  Kernel *grm = new Kernel(genotype);
  delete genotype;
  Matrix * m3 = new Matrix(grm->kernel);
  delete grm;
  
  m2->showGlobal("m2", true, 15);
  m3->showGlobal("m3", true, 15);
  
  Matrix * temp = new Matrix();
  
  misc.message << "tr(m2*m2)" << std::endl;
  std::cout << "tr1 " << m2->traceOfMatrixProduct(m2) << std::endl;
  temp->multiply(m2, 'N', m2, 'N');
  std::cout << "tr2 " << temp->trace() << std::endl;
  misc.message << "difference: " << (m2->traceOfMatrixProduct(m2)-temp->trace())/temp->trace() << std::endl;
  communicator->barrier();
  
  misc.message << "tr(m2*m3)" << std::endl;
  std::cout << "tr1 " << m2->traceOfMatrixProduct(m3) << std::endl;
  temp->multiply(m2, 'N', m3, 'N');
  std::cout << "tr2 " << temp->trace() << std::endl;
  /*temp->multiply(m2, 'N', m3, 'T');
  std::cout << "tr22" << temp->trace();*/
  misc.message << "difference: " << (m2->traceOfMatrixProduct(m3)-temp->trace())/temp->trace() << std::endl;
  communicator->barrier();
  
  misc.message << "tr(m3*m2)" << std::endl;
  std::cout << "tr1 " << m3->traceOfMatrixProduct(m2) << std::endl;
  temp->multiply(m3, 'N', m2, 'N');
  std::cout << "tr2 " << temp->trace() << std::endl;
  misc.message << "difference: " << (m3->traceOfMatrixProduct(m2)-temp->trace())/temp->trace() << std::endl;
  communicator->barrier();
  
  misc.message << "tr(m3*m3)" << std::endl;
  std::cout << "tr1 " << m3->traceOfMatrixProduct(m3) << std::endl;
  temp->multiply(m3, 'N', m3, 'N');
  std::cout << "tr2 " << temp->trace() << std::endl;
  misc.message << "difference: " << (m3->traceOfMatrixProduct(m3)-temp->trace())/temp->trace() << std::endl;
  communicator->barrier();
  
  delete m1;
  delete m2;
  delete m3;
  
}

void test26()
{
  Matrix m = Matrix(cyclicDistribution);
  Matrix mt = Matrix(cyclicDistribution);
  Matrix mf = Matrix(cyclicDistribution);
  
  m.debugRead("matrix11.dat", 11, 7);
  mt.showGlobal("Matrix11");
  mt.transpose(&m);
  mt.showGlobal("Matrix11 Transpose");
  
  std::vector<int> fr;
  fr.push_back(0);
  fr.push_back(6);
  
  std::vector<int> fc;
  fc.push_back(0);
  fc.push_back(4);
  fc.push_back(10);
  
  mt.filterRowsAndColumns(&mf, fr, fc, false);
  mf.showGlobal("Matrix11 Trasnposed filtered");
  
  mt.filterRowsAndColumns(&mf, fr, fc, true);
  mf.showGlobal("Matrix11 Trasnposed filtered");
}

void test27()
{
  Genotype *genotype = new Genotype(options.genotypeFile);
  Kernel *grm = new Kernel(genotype);
  delete genotype;
  Matrix m(grm->kernel);
  delete grm;
  
  m.fillDiagonal(2, 0.1);
  m.symmetric = true;
  m.uplo = 'B';
  m.vector = false;
  
  std::ofstream file;
  file.open("tmp/test.test27.bin", std::ofstream::out | std::ofstream::binary);
  m.writeMatrixFilev2(file);
  file.close();
  
  Matrix * mread = new Matrix();
  std::cout << mread->symmetric << " " << mread->uplo << " " << mread->vector << std::endl;
  
  std::ifstream file2;
  file2.open("tmp/test.test27.bin", std::ifstream::in | std::ifstream::binary);
  mread->readMatrixFilev2(file2, m.nGlobRows, m.nGlobCols, m.nBlockRows, m.nBlockCols);
  file2.close();
  
  //   if(communicator->mpiRoot || communicator->mpiRank==2)
  //   {
  //     mread->symmetric = true;
  //     mread->uplo = 'X';
  //   }
  
  mread->showGlobal("mread");
  mread->showGlobal("mread", false);
  mread->showPartial(row, false);
}

void test28()
{
  //communicator->nDefaultBlockRows = 3;
  //communicator->nDefaultBlockCols = 2;
  //Tested with all combinations of nBlockRows = 2 or 3 and nBlock Cols = 2 or 3 (4 combinations). Tried also nBlockRows == nBlockCols == 34
  //Tested all this combinations multiplying m*mT and mT*m. (i.e., changing the order of 'N' and 'T' parameters on multiply method of m*m)
  
  Matrix m = Matrix(cyclicDistribution);
  Matrix ms1 = Matrix(cyclicDistribution);
  Matrix ms2 = Matrix(cyclicDistribution);
  Matrix ms1t = Matrix(cyclicDistribution);
  Matrix ms2t= Matrix(cyclicDistribution);
  
  m.debugRead("testMatrices/matrix14.dat", 4, 6);
  ms1t.multiply(&m, 'N', &m, 'T');
  ms2t.multiply(&m, 'N', &m, 'T', -1.);
  
  ms1.transpose(&ms1t);
  ms2.transpose(&ms2t);
  
  ms1.showGlobal("ms1", false);
  ms2.showGlobal("ms2", false);
  
  Matrix msPacked = Matrix(cyclicDistribution);
  
  msPacked.packMatrices(&ms1, &ms2);
  
  msPacked.showGlobal("packed");
  
  Matrix msu1 = Matrix(cyclicDistribution);
  Matrix msu2 = Matrix(cyclicDistribution);
  
  msPacked.unpackMatrices(&msu1, &msu2);
  
  msu1.showGlobal("msu1", false);
  msu2.showGlobal("msu2", false);
  
  ms1.symmetrizeTriangularMatrix();
  ms2.symmetrizeTriangularMatrix();
  msu1.symmetrizeTriangularMatrix();
  msu2.symmetrizeTriangularMatrix();
  
  msu1.showGlobal("msu1", false);
  msu2.showGlobal("msu2", false);
  
  compareMatrices(&ms1, &msu1);
  compareMatrices(&ms2, &msu2);
  //compareGlobalMatrices(&ms1, &msu1);
  //compareGlobalMatrices(&ms2, &msu2);
  
  //msu1.showPartial(row);
}

void test29()
{
  Matrix m = Matrix(cyclicDistribution);
  Matrix ms1t = Matrix(cyclicDistribution);
  Matrix ms2t= Matrix(cyclicDistribution);
  
  m.debugRead("testMatrices/matrix14.dat", 4, 6);
  ms1t.multiply(&m, 'N', &m, 'T');
  ms2t.multiply(&m, 'N', &m, 'T', -1.);
  
  ms1t.showGlobal("ms1T", false);
  ms2t.showGlobal("ms2T", false);
  
  
  Matrix *ms1;
  Matrix *ms2;
  if(true)
  {
    ms1 = new Matrix(cyclicDistribution);
    ms2 = new Matrix(cyclicDistribution);
    ms1->transpose(&ms1t);
    ms2->transpose(&ms2t);
  }
  else
  {
    ms1 = new Matrix(&ms1t);
    ms2 = new Matrix(&ms2t);
  }
  
  
  ms1->showGlobal("ms1", false);
  ms2->showGlobal("ms2", false);
  
  misc.message << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  misc.message << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  misc.message << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  
  Matrix msPacked = Matrix(cyclicDistribution);
  
  msPacked.packMatrices(ms1, ms2);
  
  msPacked.showGlobal("packed");
  //msPacked.showPartial(row);
  msPacked.writeMatrixMPI("tmp/ehem.bin");
  
  misc.message << "Reading..." << std::endl;
  
  Matrix msu1 = Matrix(cyclicDistribution);
  Matrix msu2 = Matrix(cyclicDistribution);
  //Matrix msuPacked = Matrix(cyclicDistribution, 5, 4, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  for(int i = 1; i< 5; i++)
  {
    for(int j = 1; j< 5; j++)
    {
      //communicator->nDefaultBlockRows = i;
      //communicator->nDefaultBlockCols = j;
      
      misc.message << "Testing (" << i << ", " << j << ")" << std::endl;
      
      Matrix msuPacked = Matrix(cyclicDistribution, 5, 4, i, j);
      msuPacked.readMatrixMPI("tmp/ehem.bin");
      
      msuPacked.unpackMatrices(&msu1, &msu2);
      
//       msuPacked.showGlobal("msuPacked", false);
//       msu1.showGlobal("msu1", false);
//       msu2.showGlobal("msu2", false);
      
      ms1->symmetrizeTriangularMatrix();
      ms2->symmetrizeTriangularMatrix();
      msu1.symmetrizeTriangularMatrix();
      msu2.symmetrizeTriangularMatrix();
      
      /*ms1.showGlobal("ms1", false);
      ms2.showGlobal("ms2", false);
      msu1.showGlobal("msu1", false);
      msu2.showGlobal("msu2", false);*/
      
      /*compareMatrices(&msPacked, &msuPacked);
      compareMatrices(&ms1, &msu1);
      compareMatrices(&ms2, &msu2);*/
        
      compareGlobalMatrices(&msPacked, &msuPacked);
      compareGlobalMatrices(ms1, &msu1);
      compareGlobalMatrices(ms2, &msu2);
    }
  }
  delete ms1;
  delete ms2;
}

void test29bis()
{
  //Tried with repetitions without error for all block size combinations and (i.e. by changing i and j) and using writeMatrixFile option.
  //x=0; while [ $x -lt 30 ]; do mpirun --np 9 ./geneasy >> test.log; x=$(($x + 1)); done
  //x=0; while [ $x -lt 30 ]; do mpirun --np 6 ./geneasy >> test2.log; x=$(($x + 1)); done
  //x=0; while [ $x -lt 30 ]; do mpirun --np 4 ./geneasy >> test3.log; x=$(($x + 1)); done
  
  unsigned char header[] = {'T', 'E', 'S', 'T'};
  
  for(int i = 1; i< 5; i++)
  {
    for(int j = 1; j< 5; j++)
    {
      //int i = 1; int j = 3;
      
      Matrix m = Matrix(cyclicDistribution);
      Matrix ms1 = Matrix(cyclicDistribution);
      Matrix ms2 = Matrix(cyclicDistribution);

      misc.message << "Testing (" << i << ", " << j << ")" << std::endl;      
      
      ms1.debugRead("matrix14bis2.dat", 4, 4, 0, 0, i, j);
      ms1.symmetric = true;
      ms1.uplo = 'L';
      ms2.debugRead("matrix14bis1.dat", 4, 4, 0, 0, i, j);
      ms2.symmetric = true;
      ms2.uplo = 'L';
      
//       ms1.showGlobal("ms1", false);
//       ms2.showGlobal("ms2", false);
      
//       misc.message << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
//       misc.message << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
//       misc.message << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      
      Matrix msPacked = Matrix(cyclicDistribution);
      
      msPacked.packMatrices(&ms1, &ms2);
      msPacked.fillWithConstant(0.);
      msPacked.packMatrices(&ms1, &ms2);
      
//       msPacked.showGlobal("packed", false);
//       msPacked.showPartial(row);
      int offset = 0;
      std::string fname = "ehem2.bin";
      std::vector<char> fn(fname.begin(), fname.end());
      fn.push_back('\0');
      MPI_File_delete (&fn[0], MPI_INFO_NULL);
      offset = 4;
      msPacked.writeMatrixMPI("ehem2.bin", (char*)header, 4);
      communicator->barrier();
      //---------------------
//       msPacked.writeMatrixMPI(fname);
      //---------------------
      //offset = 4;
      //std::ofstream file;
      //misc.openOutputFileInRoot(file, fname);
      //msPacked.writeMatrixFile(file);
      //if(communicator->mpiRoot)
      //{
      //  file.close();
      //}
      //---------------------
      
      
      Matrix msu1 = Matrix(cyclicDistribution);
      Matrix msu2 = Matrix(cyclicDistribution);
      
      Matrix msuPacked = Matrix(cyclicDistribution, 5, 4, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
      msuPacked.fillWithConstant(0.);
      msuPacked.readMatrixMPI("ehem2.bin", offset);
//       msuPacked.showGlobal("uPacked");
      
      msuPacked.unpackMatrices(&msu1, &msu2);
      
      ms1.symmetrizeTriangularMatrix();
      ms2.symmetrizeTriangularMatrix();
      msu1.symmetrizeTriangularMatrix();
      msu2.symmetrizeTriangularMatrix();
            
//       msu1.showGlobal("msu1", false);
      
      compareGlobalMatrices(&msPacked, &msuPacked);
      compareGlobalMatrices(&ms1, &msu1);
      compareGlobalMatrices(&ms2, &msu2);
    }
  }
}

void test29bis2()
{
  Matrix m = Matrix(cyclicDistribution);
  Matrix ms1t = Matrix(cyclicDistribution);
  Matrix ms2t= Matrix(cyclicDistribution);
  
  m.debugRead("testMatrices/matrix14.dat", 4, 6);
  ms1t.multiply(&m, 'N', &m, 'T');
  ms2t.multiply(&m, 'N', &m, 'T', -1.);
  
  ms1t.showGlobal("ms1T", false);
  ms2t.showGlobal("ms2T", false);
  
  
  Matrix *ms1u;
  Matrix *ms1l;
  Matrix *ms2u;
  Matrix *ms2l;
  ms1l = new Matrix(&ms1t);
  ms1u = new Matrix(cyclicDistribution);
  ms1u->transpose(&ms1t);
  
  ms2l = new Matrix(&ms2t);
  ms2u = new Matrix(cyclicDistribution);
  ms2u->transpose(&ms2t);
 
  
  ms1l->showGlobal("ms1l", false);
  ms1u->showGlobal("ms1u", false);
  ms2l->showGlobal("ms2l", false);
  ms2u->showGlobal("ms2u", false);
  
  misc.message << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  misc.message << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  misc.message << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  
  Matrix msPackedll = Matrix(cyclicDistribution);
  Matrix msPackeduu = Matrix(cyclicDistribution);
  Matrix msPackedul = Matrix(cyclicDistribution);
  Matrix msPackedlu = Matrix(cyclicDistribution);
  
  msPackedll.packMatrices(ms1l, ms2l);
  msPackeduu.packMatrices(ms1u, ms2u);
  msPackedul.packMatrices(ms1u, ms2l);
  msPackedlu.packMatrices(ms1l, ms2u);
  
  msPackedll.showGlobal("packed ll");
  msPackeduu.showGlobal("packed uu");
  msPackedul.showGlobal("packed ul");
  msPackedlu.showGlobal("packed lu");
  
  compareMatrices(&msPackedll, &msPackeduu);
  compareMatrices(&msPackedll, &msPackedlu);
  compareMatrices(&msPackedll, &msPackedul);
  
  delete ms1u;
  delete ms1l;
  delete ms2u;
  delete ms2l;
}

void test30()
{
  Matrix m = Matrix(cyclicDistribution);
  Matrix mt = Matrix(cyclicDistribution);
  
  Matrix ms1 = Matrix(cyclicDistribution);
  Matrix ms1t = Matrix(cyclicDistribution);
  
  m.debugRead("matrix14.dat", 4, 6);
  mt.transpose(&m);
  m.showGlobal("m", false);
  mt.showGlobal("mT", false);
  
  ms1t.multiply(&m, 'T', &m, 'N');

  ms1.transpose(&ms1t);

  ms1.showGlobal("ms1", false);
  ms1t.showGlobal("ms1T", false);
  
  ms1t.symmetrizeTriangularMatrix();
  ms1t.showGlobal("ms1T", false);
  
  ms1t.fillDiagonal(2., -1.3);
  ms1t.showGlobal("ms1T", false);
}

void test31()
{
  std::vector<std::string> sample;
  
  for(int i = 0; i<500; i++)
  {
    std::stringstream ss;
    ss << i;
    sample.push_back(ss.str());
  }

  for(int i = 0; i<5; i++)
  {
    int nElements = (i+1)*10;
    std::vector<std::string> randomSample = getRandomSample(sample, nElements);
    if(randomSample.size() != nElements)
    {
      misc.error("Error on test31", 0);
    }
    misc.message << nElements << " " << randomSample.size() << " " << sample.size() << std::endl;
    for(int j = 0; j < randomSample.size(); j++) misc.message << randomSample[j] << " ";
    misc.message << std::endl << "-------------------------------" << std::endl;
  }
}

void test32()
{
  //Tested also opening a grm with ghex with both methods. Making headers equal and then comparing with diff.
  Analysis analysis;
  
  options.outFile = "test.test20";
  options.grmFile = options.outFile;
  options.genotypeFile = "test/parts.snps/merged";
  options.genotypeListFile = "test/parts.snps/list.dat";
  options.GRMJoinMethod = 0;
  
  options.grmFile = "";
  options.genotypeFile = "test/parts.snps/merged";
  options.genotypeListFile = "";
  options.GRMJoinMethod = 0;
  Kernel* grm1 = loadGRMUsingOptions();
  grm1->kernel->symmetrizeTriangularMatrix();
  grm1->N->symmetrizeTriangularMatrix();
  grm1->writeGRM(options.outFile);
  options.useMPIForWriting = false;
  grm1->writeGRM(options.outFile + ".noMPI");

  options.grmFile = options.outFile;
  options.genotypeFile = "";
  options.genotypeListFile = "";
  options.GRMJoinMethod = 0;
  Kernel* grm4 = loadGRMUsingOptions();
  grm4->kernel->symmetrizeTriangularMatrix();
  grm4->N->symmetrizeTriangularMatrix();
  
  
  options.grmFile = options.outFile + ".noMPI";
  options.genotypeFile = "";
  options.genotypeListFile = "";
  options.GRMJoinMethod = 0;
  Kernel* grm5 = loadGRMUsingOptions();
  grm5->kernel->symmetrizeTriangularMatrix();
  grm5->N->symmetrizeTriangularMatrix();
  
  compareGlobalMatrices(grm1->kernel, grm4->kernel);
  compareGlobalMatrices(grm1->N, grm4->N);
  
  compareGlobalMatrices(grm1->kernel, grm5->kernel);
  compareGlobalMatrices(grm1->N, grm5->N);
}

/**
 * Test33
 */
void test33()
{
  Genotype *genotype = new Genotype(options.genotypeFile);
  Kernel grm = Kernel(genotype);
  delete genotype;
  
  Matrix *r = new Matrix(cyclicDistribution);
  Matrix *r2 = new Matrix(cyclicDistribution);
  Matrix *r3 = new Matrix(cyclicDistribution);
  
  double determinant;
  double determinant2;

  r = new Matrix(grm.kernel);
  
  r2->multiply(r, 'N', r, 'N');
  r2->symmetric = true;
  delete r;
  r = new Matrix(r2);
  r3 = new Matrix(r2);
  
  r2->showGlobal("grm");
  r2->symmetricInvert(&determinant);
  std::cout << "Determinant: " << std::setprecision(9) << determinant << std::endl; std::cout.flush();
  r2->showGlobal("(grm)^-1");
  
  r3->symmetricInvert(&determinant2, true);
  std::cout << "Determinant2: " << std::setprecision(9) << determinant2 << std::endl; std::cout.flush();
  r3->showGlobal("(grm)^-1 single");
  
  compareMatrices(r2, r3);
  
  Matrix * identity = new Matrix(cyclicDistribution);
  identity->multiply(r3, 'N', r, 'N');
  identity->showGlobal("identity", true, 5, 1e-10);
  identity->showGlobal("identity", true, 5, 1e-2);
  
  Matrix m =  Matrix(cyclicDistribution);
  m.debugRead("matrix14.dat", 4, 6);
  r->multiply(&m, 'T', &m, 'N');
  r->invert();
  r->multiply(&m, 'N', &m, 'T');
  r->invert();
}

void test34()
{
  for(int i = 1; i < 5; i++)
  {
    for(int j = 1; j < 5; j++)
    { 
      Matrix m = Matrix(cyclicDistribution);
      Matrix mt = Matrix(cyclicDistribution);
      Matrix mf = Matrix(cyclicDistribution);
      
      m.debugRead("matrix11.dat", 11, 7, 0, 0, i, j);
      //m.showGlobal("Matrix11");
      mt.transpose(&m);
      //mt.showGlobal("Matrix11 Transpose");
      
      std::vector< std::vector<double> > gm;
      std::vector< std::vector<double> > gmt;
      m.matrixToStandardVector(gm);
      mt.matrixToStandardVector(gmt);
      
      if(communicator->mpiRoot)
      {
        bool flag = true;
        //gm[0][0] = 0.; //Just for checking the alternative
        for(int r = 0; r < gm.size(); r++)
        {
          for(int c = 0; c < gm[0].size(); c++)
          {
            if(gm[r][c] != gmt[c][r])
            {
              flag = false;
              misc.message << "Differ: " << gm[r][c] << " " <<  gmt[c][r] << std::endl;
            }
          }
        }
        if(flag == false)
        {
          misc.error("Matrices differ!!!", 0);
        }
      }
      
    }
  }
  
  Matrix m = Matrix(cyclicDistribution);
  Matrix ms1 = Matrix(cyclicDistribution);
  Matrix ms1t = Matrix(cyclicDistribution);
  Matrix ms1tt = Matrix(cyclicDistribution);
  
  m.debugRead("matrix14.dat", 4, 6);
  
  ms1.multiply(&m, 'T', &m, 'N');

  ms1t.transpose(&ms1);
  ms1tt.transpose(&ms1t);

  ms1.showGlobal("ms1", false);
  ms1t.showGlobal("ms1T", false);
  ms1tt.showGlobal("ms1TT", false);

  ms1.symmetrizeTriangularMatrix();
  ms1t.transpose(&ms1);
  ms1.showGlobal("ms1", false);
  ms1t.showGlobal("ms1T", false);

}

void test35()
{
  std::string s;
  if(communicator->mpiRoot)
  {
    s = "helloworld";
  }
  for(int i = 0; i<communicator->mpiNumTasks; i++)
  {
    if(i == communicator->mpiRank)
    {
      std::cout << i << " nob[" << s << "]" << std::endl;
    }
    communicator->barrier();
  }
  communicator->broadcast(s);
  for(int i = 0; i<communicator->mpiNumTasks; i++)
  {
    if(i == communicator->mpiRank)
    {
      std::cout << i << " sib[" << s << "]" << std::endl;
    }
    communicator->barrier();
  }
}

void test36()
{
  int ndfvalues = 10;
  int nxvalues = 5;
  misc.message << "[" << std::endl;
  for(int i = 1; i < ndfvalues; i++)
  {
    for(int j = 1; j < ndfvalues; j++)
    {
      for(int k = 1; k < nxvalues; k++)
      {
        misc.message << "(";
        misc.message << i << "., ";
        misc.message << j << "., ";
        misc.message << k << "., ";
        misc.message << FStatCDF(double(i), double(j), double(k));
        misc.message << "),";
      }
    }
  }
  misc.message << "]" << std::endl;
}

void test37()
{
  //int nr = 5, nc = 6;
  //int nr = 5, nc = 5;
  //int nr = 5467, nc = 5467;
  int nr = 2619, nc = 2619;
  double * gm;
  if(communicator->mpiRoot)
  {
    gm = new double [nr*nc];
    for(int r = 0; r<nr; r++)
    {
      for(int c = 0; c<nc; c++)
      {
        gm[c*nr + r] = -1.;
        if(c == r)
        {
          //gm[c*nr + r] = r+1; //Check null hypothesis
          gm[c*nr + r] = r;
        }
      }
    }
  }
  
  //int nbr = 2, nbc = 3;
  //int nbr = 2, nbc = 2;
  //int nbr = 29, nbc = 29;
  int nbr = 67, nbc = 67;
  Matrix m(cyclicDistribution, nr, nc, nbr, nbc);
  m.scatterMatrix(gm);
  //m.showGlobal();
  std::vector<double> diagonal = m.diagonal();
  for(int id = 0; id<diagonal.size(); id++)
  {
    misc.message << diagonal[id] << " ";
    if(diagonal[id] != double(id))
    {
      misc.error("\n------------------------\nError: diagonal differ!\n------------------------\n", 0);
    }
  }
  
  if(diagonal.size() != nr && communicator->mpiRoot == true)
  {
    misc.error("\n------------------------\nError: diagonal dimensions differ!\n------------------------\n", 0);
  }
  
  misc.message << std::endl << diagonal.size() << " " << nr << " " << nc << std::endl;
  misc.message << m.nBlockRows << " " << m.nBlockCols << std::endl;
  
  if(communicator->mpiRoot)
  {
    delete [] gm;
  }
}

void test38()
{
  Covariate *covariate = new Covariate(options.covarsFile, options.qCovarsFile, std::vector<std::string>());
  covariate->printCovariate();
  
  std::vector<std::string> test;
  if(communicator->mpiRoot)
  {
    test.push_back("i2@i2");
    test.push_back("i4@i4");
    test.push_back("i6@i6");
  }
  
  covariate->filterIndividuals(test);
  covariate->printCovariate();
  delete covariate;
}

void test39()
{
  if(communicator->mpiRoot)
  {
    long seed = options.randomSeed;
    for(int i = 0; i < 10; i++)
    {
      misc.message << box_muller(0., 10., &seed) << std::endl;
    }
  }
}

void test40()
{
  communicator->nDefaultBlockRows = 3;
  communicator->nDefaultBlockCols = 5;
  
  Matrix m = Matrix(cyclicDistribution);
  m.debugRead("testMatrices/matrix14.dat", 4, 6);
  Matrix mf = Matrix(cyclicDistribution);
  
  
  
  std::map<int, int> rowsOriginDestination;
  rowsOriginDestination[0] = 2;
  rowsOriginDestination[2] = 0;
  rowsOriginDestination[3] = 1;
  rowsOriginDestination[1] = 3;
  
  std::map<int, int> colsOriginDestination;
  colsOriginDestination[0] = 1;
  colsOriginDestination[5] = 0;
  colsOriginDestination[3] = 2;
  colsOriginDestination[2] = 3;
  
  
  m.generalResorting(&mf, rowsOriginDestination, colsOriginDestination);
  
  m.showGlobal("original");
  mf.showGlobal("filtered");
}

void test41()
{
  std::srand(communicator->mpiRank*std::time(0));
  int nr = 1000;
  int nc = 2210;
  
  //int nr = 37;
  //int nc = 22;
  
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 3;
  //int nr = 10;
  //int nc = 22;
  
  double * original = new double [nr*nc];
  for(int i = 0; i<nr*nc; i++ ) {original[ i ] = i;}

  for(int i = 0; i<20; i++)
  {
  
    std::map<int, int> rowsOriginDestination;
    std::map<int, int> colsOriginDestination;
    
    std::vector<int> keepIndexsRows;
    std::vector<int> keepIndexsCols;
    for(int r = 0; r<nr; r++) {keepIndexsRows.push_back(r);}
    for(int c = 0; c<nc; c++) {keepIndexsCols.push_back(c);}
    std::random_shuffle ( keepIndexsRows.begin(), keepIndexsRows.end() );
    std::random_shuffle ( keepIndexsCols.begin(), keepIndexsCols.end() );
    
    communicator->broadcast(&(keepIndexsRows[0]), keepIndexsRows.size());
    communicator->broadcast(&(keepIndexsCols[0]), keepIndexsCols.size());
    
    if(keepIndexsRows == keepIndexsCols) {misc.error("Unexpected error on test41.", 0);};
    int nfr = nr/2;
    int nfc = nc/3;
    for(int r = 0; r<nfr; r++) {rowsOriginDestination[ keepIndexsRows[r] ] = r; /*std::cout << communicator->mpiRank << " " << keepIndexsRows[r] << " " << r << std::endl;*/}
    //std::cout << "-" << std::endl;
    for(int c = 0; c<nfc; c++) {colsOriginDestination[ keepIndexsCols[c] ] = c; /*std::cout << communicator->mpiRank << " " << keepIndexsCols[c] << " " << c << std::endl;*/}
    //std::cout << "-" << std::endl;
    
    communicator->nDefaultBlockRows = (rand() % 10) + 3; //7;
    communicator->nDefaultBlockCols = (rand() % 10) + 3; //3;
    communicator->broadcast(&communicator->nDefaultBlockRows);
    communicator->broadcast(&communicator->nDefaultBlockCols);
    
    int nobsr = (rand() % 10) + 3; //5;
    int nobsc = (rand() % 10) + 3; //13;
    communicator->broadcast(&nobsr);
    communicator->broadcast(&nobsc);
    
    test42(original, nr, nc, nobsr, nobsc, rowsOriginDestination, colsOriginDestination);
  
  }
  
  delete [] original;
}

void test42(double *original, int nr, int nc, int bnr, int bnc, std::map<int, int> & rowsOriginDestination, std::map<int, int> & colsOriginDestination)
{
  
  Matrix originalDistributed = Matrix(cyclicDistribution, nr, nc, bnr, bnc);
  Matrix filteredDistributed = Matrix(cyclicDistribution);
  
  originalDistributed.scatterMatrix(original);
  originalDistributed.generalResorting(&filteredDistributed, rowsOriginDestination, colsOriginDestination);
  double * filtered = new double [ rowsOriginDestination.size() * colsOriginDestination.size() ];
  if( filteredDistributed.nGlobRows != rowsOriginDestination.size() || filteredDistributed.nGlobCols != colsOriginDestination.size() ) {misc.error("Unexpected error on test42.", 0);}
  filteredDistributed.gatherMatrix(filtered);
  
  //originalDistributed.showGlobal("Original");
  //filteredDistributed.showGlobal("Filtered");
  
  if(communicator->mpiRoot == true)
  {
    //This must be only on mpiroot when RAW
    double * filteredRaw = new double [ rowsOriginDestination.size() * colsOriginDestination.size() ];
    for(std::map<int, int>::iterator itc = colsOriginDestination.begin(); itc != colsOriginDestination.end(); ++itc)
    {
      for(std::map<int, int>::iterator itr = rowsOriginDestination.begin(); itr != rowsOriginDestination.end(); ++itr)
      {
        filteredRaw[ itc->second*rowsOriginDestination.size() + itr->second ] = original[ itc->first*nr + itr->first ];
      }
    }
    
    //filteredRaw[0] = 0; //Just as a NULL case
    
    for(int i = 0; i<rowsOriginDestination.size() * colsOriginDestination.size(); i++ )
    {
      //std::cout << filteredRaw[i] << " " << filtered[i] << std::endl;
      if( filteredRaw[i] != filtered[i] )
      {
        std::cout << "\n\n Unequal matrices \n\n" << std::endl;
        misc.error("Test42 failed.", 0);
      }
    }
    delete [] filteredRaw;
    std::cout << "Comprovation finished succesfully. (" << originalDistributed.nBlockRows << ", " << originalDistributed.nBlockCols << ") (" << filteredDistributed.nBlockRows << ", " << filteredDistributed.nBlockCols << ")" << std::endl;
  }
  
  
  delete [] filtered;
  
}

void test43ARCHER()
{
  std::srand(communicator->mpiRank*std::time(0));
  int nr = 6000;
  int nc = 6000;
  
  //int nr = 6;
  //int nc = 6;
  
  double * original = new double [nr*nc];
  for(int i = 0; i<nr*nc; i++ ) {original[ i ] = i;/*NAN;*/}
  
  Matrix originalDistributed = Matrix(cyclicDistribution, nr, nc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  //originalDistributed.scatterMatrix(original);
  Matrix * filteredDistributed1 = new Matrix(cyclicDistribution);
  Matrix * filteredDistributed2 = new Matrix(cyclicDistribution);
  Matrix * filteredDistributed1Inv = new Matrix(cyclicDistribution);
  Matrix * filteredDistributed2Inv = new Matrix(cyclicDistribution);
  
  std::vector<int> randomIndexsRows;
  std::vector<int> randomIndexsCols;
  for(int r = 0; r<nr; r++) {randomIndexsRows.push_back(r);}
  for(int c = 0; c<nc; c++) {randomIndexsCols.push_back(c);}
  std::random_shuffle ( randomIndexsRows.begin(), randomIndexsRows.end() );
  std::random_shuffle ( randomIndexsCols.begin(), randomIndexsCols.end() );
  
  int nfr = nr/2;
  int nfc = nc/3;
  std::sort(randomIndexsRows.begin(), randomIndexsRows.begin() + nfr);
  std::sort(randomIndexsCols.begin(), randomIndexsCols.begin() + nfc);
  randomIndexsRows.erase(randomIndexsRows.begin() + nfr - 1,randomIndexsRows.end());
  randomIndexsCols.erase(randomIndexsCols.begin() + nfc - 1,randomIndexsCols.end());
  
  misc.setGetElapsedTime("oldMethod");
  originalDistributed.filterRowsAndColumnsOld(filteredDistributed1, randomIndexsRows, randomIndexsCols);
  misc.message << "Old method: " << misc.setGetElapsedTime("oldMethod") << " " << filteredDistributed1->nGlobRows << " " << filteredDistributed1->nGlobCols << std::endl;
  
  //delete filteredDistributed1;
  
  misc.setGetElapsedTime("newMethod");
  originalDistributed.filterRowsAndColumns(filteredDistributed2, randomIndexsRows, randomIndexsCols);
  misc.message << "New method: " << misc.setGetElapsedTime("newMethod") << " " << filteredDistributed2->nGlobRows << " " << filteredDistributed2->nGlobCols << std::endl;
  
  compareMatrices(filteredDistributed1, filteredDistributed2);
  //compareMatrices(&originalDistributed, filteredDistributed2);
  
  
  //Repeat when keep = False
  
  misc.setGetElapsedTime("oldMethodInv");
  originalDistributed.filterRowsAndColumnsOld(filteredDistributed1Inv, randomIndexsRows, randomIndexsCols, false);
  misc.message << "Old methodInv: " << misc.setGetElapsedTime("oldMethodInv") << " " << filteredDistributed1Inv->nGlobRows << " " << filteredDistributed1Inv->nGlobCols << std::endl;
  
  misc.setGetElapsedTime("newMethodInv");
  originalDistributed.filterRowsAndColumns(filteredDistributed2Inv, randomIndexsRows, randomIndexsCols, false);
  misc.message << "New methodInv: " << misc.setGetElapsedTime("newMethodInv") << " " << filteredDistributed2Inv->nGlobRows << " " << filteredDistributed2Inv->nGlobCols << std::endl;
  
  compareMatrices(filteredDistributed1Inv, filteredDistributed2Inv);
  
  //filteredDistributed1Inv->showGlobal("1Inv");
  //filteredDistributed2Inv->showGlobal("2Inv");
  
  delete [] original;
}

void test44()
{
  int randombr = 3;
  int randombc = 7;
  Matrix * m = new Matrix(diagonalDistribution, 17, 17, randombr, randombc);
  misc.message << m->vector << std::endl;
  //Matrix * m2 = new Matrix(cyclicDistribution, 40, 1);
  Matrix * m2 = new Matrix(cyclicDistribution, 10, 10, randombr, randombc);
  misc.message << m2->vector << std::endl;
  //Matrix * m3 = new Matrix(diagonalDistribution, 40, 1);
  Matrix * m3 = new Matrix(cyclicDistribution, 5, 7);
  
  m->fillWithConstant(0.);
  //m->fillDiagonal(1.1, 0.2);
  m->fillDiagonal(1.1, 0.);
  
  m2->fillWithConstant(0.);
  m2->fillDiagonal(1.1, 0.2);
  
  m->showGlobal("m");
  m2->showGlobal("m2");
  
  double * diagonal1;
  double * diagonal2;
  int dimension1;
  int dimension2;
  if(communicator->mpiRoot == true)
  {
    dimension1 = m->nGlobRows;
    dimension2 = m2->nGlobRows;
    diagonal1 = new double [m->nGlobRows];
    diagonal2 = new double [m2->nGlobRows];
    for(int i = 0; i < m->nGlobRows; i++)
    {
      diagonal1[i] = i + .1;
    }
    for(int i = 0; i < m2->nGlobRows; i++)
    {
      diagonal2[i] = i + .1;
    }
  }
  
  m->setDiagonal(diagonal1, dimension1);
  m2->setDiagonal(diagonal2, dimension2);
  
  
  m->showGlobal("m");
  m2->showGlobal("m2");
  
  Matrix * mCopy = new Matrix(m);
  Matrix * m2Copy = new Matrix(m2);
  
  mCopy->showGlobal("mCopy");
  m2Copy->showGlobal("m2Copy");
  
  compareMatrices(m, mCopy);
  compareMatrices(m2, m2Copy);
}

void test45()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 3;
  
  Matrix * m1 = new Matrix();
  m1->debugRead("testMatrices2/data/random_1_11x7.dat", 11, 7);
  Matrix * m2 = new Matrix();
  m2->debugRead("testMatrices2/data/random_2_7x11.dat", 7, 11);
  Matrix * m3D = new Matrix(diagonalDistribution);
  m3D->debugRead("testMatrices2/data/random_3_11x11_diagonal.dat", 11, 11);
  Matrix * m3 = new Matrix();
  m3->debugRead("testMatrices2/data/random_3_11x11.dat", 11, 11);
  Matrix * m4D = new Matrix(diagonalDistribution);
  m4D->debugRead("testMatrices2/data/random_4_11x11_diagonal.dat", 11, 11);
  
  //m3D->showGlobal("m3D");
  
  //Test diagonal by non-diagonal
  
  Matrix * mR1P = new Matrix();
  Matrix * mR2P = new Matrix();
  Matrix * mR3P = new Matrix();
  Matrix * mR4P = new Matrix();
  mR1P->debugRead("testMatrices2/data/result_1tx3.dat", 7, 11);
  mR2P->debugRead("testMatrices2/data/result_2x3.dat", 7, 11);
  mR3P->debugRead("testMatrices2/data/result_3x1.dat", 11, 7);
  mR4P->debugRead("testMatrices2/data/result_3x2t.dat", 11, 7);
  
  Matrix * mR1 = new Matrix();
  Matrix * mR2 = new Matrix();
  Matrix * mR3 = new Matrix();
  Matrix * mR4 = new Matrix();
  
  mR1->multiply(m1, 'T', m3D, 'N');
  mR2->multiply(m2, 'N', m3D, 'N');
  mR3->multiply(m3D, 'N', m1, 'N');
  mR4->multiply(m3D, 'N', m2, 'T');

  //mR1P->showGlobal("mR1P");
  //mR1->showGlobal("mR1", true, 10);
  
  compareMatrices(mR1, mR1P);
  compareMatrices(mR2, mR2P);
  compareMatrices(mR3, mR3P);
  compareMatrices(mR4, mR4P);
  
  //Test diagonal by diagonal
  
  Matrix * mR5P = new Matrix(diagonalDistribution);
  mR5P->debugRead("testMatrices2/data/result_3x4.dat", 11, 11);
  
  Matrix * mR5 = new Matrix();
  mR5->multiply(m3D, 'T', m4D, 'N');
  
  //mR5P->showGlobal("mR5P");
  //mR5->showGlobal("mR5", true, 10);
  
  compareMatrices(mR5, mR5P);
  
  //Test diagonal multiplications when result is on a submatrix
  
  Matrix * mR6 = new Matrix(cyclicDistribution, 55, 47);
  mR6->fillWithConstant(0.);
  mR6->multiply(m1, 'T', m3D, 'N', 1., subMatrix(3,5,7,11) );
  
  Matrix * temp1 =  new Matrix(cyclicDistribution, 55, 47);
  temp1->fillWithConstant(0.);
  Matrix * mR6T = new Matrix();
  mR6T->joinMatrices(temp1, subMatrix(0, 0, 55, 47), mR1, subMatrix(3,5,7,11) );
  
  //mR6->showGlobal("mR6");
  compareMatrices(mR6, mR6T);
  
  
  Matrix * mR7 = new Matrix(cyclicDistribution, 55, 47);
  mR7->fillWithConstant(0.);
  mR7->multiply(m3D, 'T', m4D, 'N', 1., subMatrix(3,5,11,11) );
  
  Matrix * temp2 =  new Matrix(cyclicDistribution, 55, 47);
  temp2->fillWithConstant(0.);
  Matrix * temp3 =  new Matrix(cyclicDistribution, 11, 11);
  temp3->fillWithConstant(0.);
  temp3->setDiagonal(mR5->m, mR5->nGlobRows);
  Matrix * mR7T = new Matrix();
  mR7T->joinMatrices(temp2, subMatrix(0, 0, 55, 47), temp3, subMatrix(3,5,11,11) );
  
  //mR7->showGlobal("mR7");
  compareMatrices(mR7, mR7T);
}

void test45bis2()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 3;
  
  Matrix * m1 = new Matrix();
  m1->debugRead("testMatrices2/data/scale_2_random_1_11x7.dat", 11, 7);
  Matrix * m2 = new Matrix();
  m2->debugRead("testMatrices2/data/scale_2_random_2_7x11.dat", 7, 11);
  Matrix * m3D = new Matrix(diagonalDistribution);
  m3D->debugRead("testMatrices2/data/scale_2_random_3_11x11_diagonal.dat", 11, 11);
  Matrix * m3 = new Matrix();
  m3->debugRead("testMatrices2/data/scale_2_random_3_11x11.dat", 11, 11);
  Matrix * m4D = new Matrix(diagonalDistribution);
  m4D->debugRead("testMatrices2/data/scale_2_random_4_11x11_diagonal.dat", 11, 11);
  
  //m3D->showGlobal("m3D");
  
  //Test diagonal by non-diagonal
  
  Matrix * mR1P = new Matrix();
  Matrix * mR2P = new Matrix();
  Matrix * mR3P = new Matrix();
  Matrix * mR4P = new Matrix();
  mR1P->debugRead("testMatrices2/data/scale_2_result_1tx3.dat", 7, 11);
  mR2P->debugRead("testMatrices2/data/scale_2_result_2x3.dat", 7, 11);
  mR3P->debugRead("testMatrices2/data/scale_2_result_3x1.dat", 11, 7);
  mR4P->debugRead("testMatrices2/data/scale_2_result_3x2t.dat", 11, 7);
  
  Matrix * mR1 = new Matrix();
  Matrix * mR2 = new Matrix();
  Matrix * mR3 = new Matrix();
  Matrix * mR4 = new Matrix();
  
  mR1->multiply(m1, 'T', m3D, 'N', 2.);
  mR2->multiply(m2, 'N', m3D, 'N', 2.);
  mR3->multiply(m3D, 'N', m1, 'N', 2.);
  mR4->multiply(m3D, 'N', m2, 'T', 2.);

  //mR1P->showGlobal("mR1P");
  //mR1->showGlobal("mR1", true, 10);
  
  compareMatrices(mR1, mR1P);
  compareMatrices(mR2, mR2P);
  compareMatrices(mR3, mR3P);
  compareMatrices(mR4, mR4P);
  
  //Test diagonal by diagonal
  
  Matrix * mR5P = new Matrix(diagonalDistribution);
  mR5P->debugRead("testMatrices2/data/scale_2_result_3x4.dat", 11, 11);
  
  Matrix * mR5 = new Matrix();
  mR5->multiply(m3D, 'T', m4D, 'N', 2.);
  
  //mR5P->showGlobal("mR5P");
  //mR5->showGlobal("mR5", true, 10);
  
  compareMatrices(mR5, mR5P);
  
  //Test diagonal multiplications when result is on a submatrix
  
  Matrix * mR6 = new Matrix(cyclicDistribution, 55, 47);
  mR6->fillWithConstant(0.);
  mR6->multiply(m1, 'T', m3D, 'N', 2., subMatrix(3,5,7,11) );
  
  Matrix * temp1 =  new Matrix(cyclicDistribution, 55, 47);
  temp1->fillWithConstant(0.);
  Matrix * mR6T = new Matrix();
  mR6T->joinMatrices(temp1, subMatrix(0, 0, 55, 47), mR1, subMatrix(3,5,7,11) );
  
  //mR6->showGlobal("mR6");
  compareMatrices(mR6, mR6T);
  
  
  Matrix * mR7 = new Matrix(cyclicDistribution, 55, 47);
  mR7->fillWithConstant(0.);
  mR7->multiply(m3D, 'T', m4D, 'N', 2., subMatrix(3,5,11,11) );
  
  Matrix * temp2 =  new Matrix(cyclicDistribution, 55, 47);
  temp2->fillWithConstant(0.);
  Matrix * temp3 =  new Matrix(cyclicDistribution, 11, 11);
  temp3->fillWithConstant(0.);
  temp3->setDiagonal(mR5->m, mR5->nGlobRows);
  Matrix * mR7T = new Matrix();
  mR7T->joinMatrices(temp2, subMatrix(0, 0, 55, 47), temp3, subMatrix(3,5,11,11) );
  
  //mR7->showGlobal("mR7");
  compareMatrices(mR7, mR7T);
}

void test45bis()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 3;
  
  Matrix * m1 = new Matrix();
  m1->debugRead("testMatrices2/data/random_5_1231x977.dat", 1231, 977);
  Matrix * m2 = new Matrix();
  m2->debugRead("testMatrices2/data/random_6_977x1231.dat", 977, 1231);
  Matrix * m3D = new Matrix(diagonalDistribution);
  m3D->debugRead("testMatrices2/data/random_7_1231x1231_diagonal.dat", 1231, 1231);
  Matrix * m3 = new Matrix();
  m3->debugRead("testMatrices2/data/random_7_1231x1231.dat", 1231, 1231);
  Matrix * m4D = new Matrix(diagonalDistribution);
  m4D->debugRead("testMatrices2/data/random_8_1231x1231_diagonal.dat", 1231, 1231);
  
  //m3D->showGlobal("m3D");
  
  //Test diagonal by non-diagonal
  
  Matrix * mR1P = new Matrix();
  Matrix * mR2P = new Matrix();
  Matrix * mR3P = new Matrix();
  Matrix * mR4P = new Matrix();
  mR1P->debugRead("testMatrices2/data/result_5tx7.dat", 977, 1231);
  mR2P->debugRead("testMatrices2/data/result_6x7.dat", 977, 1231);
  mR3P->debugRead("testMatrices2/data/result_7x5.dat", 1231, 977);
  mR4P->debugRead("testMatrices2/data/result_7x6t.dat", 1231, 977);
  
  Matrix * mR1 = new Matrix();
  Matrix * mR2 = new Matrix();
  Matrix * mR3 = new Matrix();
  Matrix * mR4 = new Matrix();
  
  mR1->multiply(m1, 'T', m3D, 'N');
  mR2->multiply(m2, 'N', m3D, 'N');
  mR3->multiply(m3D, 'N', m1, 'N');
  mR4->multiply(m3D, 'N', m2, 'T');

  //mR1P->showGlobal("mR1P");
  //mR1->showGlobal("mR1", true, 10);
  
  compareMatrices(mR1, mR1P);
  compareMatrices(mR2, mR2P);
  compareMatrices(mR3, mR3P);
  compareMatrices(mR4, mR4P);
  
  //Test diagonal by diagonal
  
  Matrix * mR5P = new Matrix(diagonalDistribution);
  mR5P->debugRead("testMatrices2/data/result_7x8.dat", 1231, 1231);
  
  Matrix * mR5 = new Matrix();
  mR5->multiply(m3D, 'T', m4D, 'N');
  
  //mR5P->showGlobal("mR5P");
  //mR5->showGlobal("mR5", true, 10);
  
  compareMatrices(mR5, mR5P);
  
  //Test diagonal multiplications when result is on a submatrix
  
  Matrix * mR6 = new Matrix(cyclicDistribution, 5555, 4777);
  mR6->fillWithConstant(0.);
  mR6->multiply(m1, 'T', m3D, 'N', 1., subMatrix(3,5,977,1231) );
  
  Matrix * temp1 =  new Matrix(cyclicDistribution, 5555, 4777);
  temp1->fillWithConstant(0.);
  Matrix * mR6T = new Matrix();
  mR6T->joinMatrices(temp1, subMatrix(0, 0, 5555, 4777), mR1, subMatrix(3,5,977,1231) );
  
  //mR6->showGlobal("mR6");
  compareMatrices(mR6, mR6T);
  
  
  Matrix * mR7 = new Matrix(cyclicDistribution, 5555, 4777);
  mR7->fillWithConstant(0.);
  mR7->multiply(m3D, 'T', m4D, 'N', 1., subMatrix(3,5,1231,1231) );
  
  Matrix * temp2 =  new Matrix(cyclicDistribution, 5555, 4777);
  temp2->fillWithConstant(0.);
  Matrix * temp3 =  new Matrix(cyclicDistribution, 1231, 1231);
  temp3->fillWithConstant(0.);
  temp3->setDiagonal(mR5->m, mR5->nGlobRows);
  Matrix * mR7T = new Matrix();
  mR7T->joinMatrices(temp2, subMatrix(0, 0, 5555, 4777), temp3, subMatrix(3,5,1231,1231) );
  
  //mR7->showGlobal("mR7");
  compareMatrices(mR7, mR7T);
}

void test45bis3()
{
  double global[4] = {2,1,0,3};
  double globalD[2] = {0.5,1.5};

  Matrix * m1 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, 2,2);
  m1->scatterMatrix(global);
  Matrix * m2D = new Matrix(diagonalDistribution, 2, 2);
  m2D->setDiagonal(globalD, 2);
  
  Matrix * mR = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, 2,2);
  Matrix * mR2 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, 3,3);

  m1->showGlobal("m1");
  m2D->showGlobal("m2D");
  
  //Test diagonal by non-diagonal
  
  mR->fillWithConstant(-1);
  mR->multiply(m1, 'T', m2D, 'N');
  mR->showGlobal("m1txm2D");
  
  mR->fillWithConstant(-1);
  mR->multiply(m1, 'N', m2D, 'N');
  mR->showGlobal("m1xm2D");
  
  mR->fillWithConstant(-1);
  mR->multiply(m2D, 'N', m1, 'N');
  mR->showGlobal("m2Dxm1");
  
  mR->fillWithConstant(-1);
  mR->multiply(m2D, 'N', m1, 'T');
  mR->showGlobal("m2Dxm1t");
  
  //Test diagonal by non-diagonal symmetric
  
  misc.message << std::endl;
  
  m1->symmetric = true;
  m1->uplo = 'L';
  
  mR->fillWithConstant(-1);
  mR->multiply(m1, 'N', m2D, 'N');
  mR->showGlobal("m1s x m2D");
  
  mR->fillWithConstant(-1);
  mR->multiply(m2D, 'N', m1, 'N');
  mR->showGlobal("m2D x m1s");
  
  m1->symmetric = true;
  m1->uplo = 'U';
  
  mR->fillWithConstant(-1);
  mR->multiply(m1, 'N', m2D, 'N');
  mR->showGlobal("m1s x m2D");
  
  mR->fillWithConstant(-1);
  mR->multiply(m2D, 'N', m1, 'N');
  mR->showGlobal("m2D x m1s");
  
  
  //Diagonal by diagonal symmetry
  
  misc.message << std::endl;
  
  mR->fillWithConstant(-1);
  mR->multiply(m2D, 'N', m2D, 'N');
  mR->showGlobal("m2D x m2D");
  
  //Submatrix symmetry
  
  misc.message << std::endl;
  
  mR2->fillWithConstant(-1);
  mR2->multiply(m2D, 'N', m2D, 'N', 1., subMatrix(1,1,2,2) );
  mR2->showGlobal("m2Dxm2D sub");

  //Matriux by matrix when they are not symmetrix in a submatrix
  
  m1->symmetric = false;
  m1->uplo = 'B';
  mR2->fillWithConstant(-1);
  mR2->multiply(m1, 'N', m1, 'T', 1., subMatrix(1,1,2,2) );
  mR2->showGlobal("m1xm1 sub");
}

void test46()
{
  int ttnr = 5;
  int ttnc = 3;
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 3;
  
  Matrix m1 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m1.fillWithConstant(1.);
  Matrix m2 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m2.fillWithConstant(2.);
  
  m1.elementWiseMultiplication(&m2);
  m1.showGlobal("m1*m2");
  
  m1.fillWithConstant(1.);
  m1.elementWiseMultiplication(&m2, 1.5);
  m1.showGlobal("m1*m2*1.5");
  
  m1.fillWithConstant(1.);
  m1.elementWiseDivision(&m2);
  m1.showGlobal("m1/m2");
  
  m1.fillWithConstant(1.);
  m1.elementWiseDivision(&m2, 4.);
  m1.showGlobal("4.*m1/m2");
  
  m1.fillWithRandom();
  Matrix m1Test = Matrix(m1);
  m2.fillWithRandom();
  m1.showGlobal("randomM1");
  m2.showGlobal("randomM2");
  m1.elementWiseMultiplication(&m2, 3.4);
  m1.showGlobal("randomM1*randomM2*3.4");
  m1.elementWiseDivision(&m2, 1./3.4);
  m1.showGlobal("randomM1/(randomM2*3.4)");
  
  compareMatrices(&m1, &m1Test);
}

void test47()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 2;
  
  std::string fname;
  double *mGlob;
  int ttnr, ttnc;
  
  ttnr = 5;
  ttnc = 3;
  fname = "testMatrices/matrix2.dat";
  
  Matrix m1 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m1.debugRead(fname, ttnr, ttnc);
  
  ttnr = 5;
  ttnc = 3;
  fname = "testMatrices/matrix5.dat";
  
  Matrix m2 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows+1, communicator->nDefaultBlockCols-1);
  m2.debugRead(fname, ttnr, ttnc);
  
  Matrix m3D = Matrix(diagonalDistribution, ttnr, ttnr, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  Matrix m4D = Matrix(diagonalDistribution, ttnr, ttnr, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  if(communicator->mpiRoot)
  {
    for(int i = 0; i < m3D.nGlobRows; i++)
    {
      m3D.m[i] = (i + 1)*10 + (i + 1);
    }
    for(int i = 0; i < m4D.nGlobRows; i++)
    {
      m4D.m[i] = (i + 1)*1000 + (i + 1)*100;
    }
  }
  
  m1.showGlobal("m1");
  m2.showGlobal("m2");
  m3D.showGlobal("m3D");
  m4D.showGlobal("m4D");
  
    //m2.add(&m1);
  //m2.add(&m1, 1., 1.);
  //m2.add(&m1, 2., 2.);
  //m2.add(&m1, 1., 2.);
  //m2.add(&m1, 1., 1., subMatrix(1,1,2,2), subMatrix(0,0,2,2));
  //m2.add(&m1, 1., 1., subMatrix(1,1,2,2), subMatrix(1,1,2,2));
  //m2.add(&m1, 1., 1., subMatrix(0,0,2,2), subMatrix(1,1,2,2));
  //m2.add(&m1, 2., 0.5, subMatrix(1,1,2,2), subMatrix(1,1,2,2));
  //m2.add(&m1, 1., 1., subMatrix(1,1,5,2), subMatrix(1,1,5,2));
  //m2.add(&m1, 1., 1., subMatrix(1,1,3,2), subMatrix(1,1,4,2));
  
  //m2.add(&m3D, 1., 1., subMatrix(1,1,3,2), subMatrix(0,0,3,2));
  //m2.add(&m3D, 1., 1., subMatrix(0,0,1,2), subMatrix(0,0,1,2));
  //m2.add(&m3D, 1., 1., subMatrix(1,1,3,2), subMatrix(1,0,3,2));
  //m2.add(&m3D, 2., 0.5, subMatrix(1,1,3,2), subMatrix(1,0,3,2));
  
  m3D.add(&m2, 1., 1., subMatrix(1,3,3,2), subMatrix(1,0,3,2));
  //m3D.add(&m2, 0.5, 2., subMatrix(1,3,3,2), subMatrix(1,0,3,2));
  
  //m3D.add(&m4D, 1., 1., subMatrix(1,1,2,2), subMatrix(3,3,2,2));
  //m3D.add(&m4D, 0.5, 2., subMatrix(1,1,2,2), subMatrix(3,3,2,2));
  
  //m3D.add(&m4D, 1., 1.);
  //m3D.add(&m4D, 0.5, 2.);
  
  m2.showGlobal("m2");
  m3D.showGlobal("m3D", true, 10);
}

void test48()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 2;
  
  int ttnr, ttnc;
 
  ttnr = 1000;
  ttnc = 1000;
  
  std::string det1 = "600.03485";
  Matrix m1 = Matrix(diagonalDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m1.debugRead("testMatrices2/data/invert.1_1000x1000_logdet_" + det1 + ".dat", ttnr, ttnc);
  
  Matrix m1Test = Matrix(diagonalDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m1Test.debugRead("testMatrices2/data/invert.1_1000x1000_logdet_" + det1 + ".dat", ttnr, ttnc);
  
  Matrix m1i = Matrix(diagonalDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m1i.debugRead("testMatrices2/data/invert.inverted.1_1000x1000.dat", ttnr, ttnc);
  
  double logdet1 = 0.;
  if(communicator->mpiRoot == true)
  {
    //m1.m[0] = -m1.m[0];
    //m1.m[1] = 0.;
  }
  //bool test = m1.symmetricInvert(&logdet1);
  //bool test = m1.symmetricInvert();
  //bool test = m1.invert();
  bool test = m1.invert(&logdet1);
  if(test == false)
  {
    misc.error("SINGULAR MATRIX", 0);
  }
  compareMatrices(&m1, &m1i, 1e-5);
  std::cout << det1 << " " << logdet1 << std::endl;
  
  //Matrix mTest = Matrix(cyclicDistribution);
  //mTest.multiply(&m1, 'T', &m1Test, 'N');
  //mTest.showGlobal("mTest");
  
  std::string det2 = "8036.45007";
  Matrix m2 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m2.debugRead("testMatrices2/data/invert.1000x1000_logdet_" + det2 + ".dat", ttnr, ttnc);
  m2.symmetric = true;
  
  Matrix m2i = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m2i.debugRead("testMatrices2/data/invert.inverted.1000x1000.dat", ttnr, ttnc);
  //m2i.symmetric = true;
  
  //m2i.showGlobal("m2i");
  
  double logdet2 = 0.;
  //m2.symmetricInvert(&logdet2);
  m2.invert(&logdet2);
  //m2.invert();
  m2.symmetrizeTriangularMatrix();
  compareMatrices(&m2, &m2i, 1e-3);
  misc.message << det2 << " " << logdet2 << std::endl;
  
}

void test49()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 2;
  
  int ttnr, ttnc;
 
  //ttnr = 100;
  //ttnc = 100;
  ttnr = 2;
  ttnc = 2;
  
  //std::string trace1 = "242.73270";
  std::string trace1 = "6.82240";
  Matrix m1 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  //m1.debugRead("testMatrices2/data/trace_1.100x100_trace_" + trace1 + ".dat", ttnr, ttnc);
  m1.debugRead("testMatrices2/data/trace_1.2x2_trace_" + trace1 + ".dat", ttnr, ttnc);
  m1.symmetric = true;
  //m1.m[2] = 23.; //Just for testing that traceofMatrixProduc symmetrizes well the matrix if only have data on on triangular part. Test for 2x2 matrices. Test passed.
  //m1.uplo = 'L';
  
  Matrix m1d = Matrix(diagonalDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  //m1d.debugRead("testMatrices2/data/trace_1d.100x100_trace_" + trace1 + ".dat", ttnr, ttnc);
  m1d.debugRead("testMatrices2/data/trace_1d.2x2_trace_" + trace1 + ".dat", ttnr, ttnc);
  
  //std::string trace2 = "235.63770";
  std::string trace2 = "1.07710";
  Matrix m2 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  //m2.debugRead("testMatrices2/data/trace_2.100x100_trace_" + trace2 + ".dat", ttnr, ttnc);
  m2.debugRead("testMatrices2/data/trace_2.2x2_trace_" + trace2 + ".dat", ttnr, ttnc);
  m2.symmetric = true;
  //m2.m[2] = 23.; //Just for testing that traceofMatrixProduc symmetrizes well the matrix if only have data on on triangular part. Test for 2x2 matrices. Test passed.
  //m2.uplo = 'L';
  
  Matrix m2d = Matrix(diagonalDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  //m2d.debugRead("testMatrices2/data/trace_2d.100x100_trace_" + trace2 + ".dat", ttnr, ttnc);
  m2d.debugRead("testMatrices2/data/trace_2d.2x2_trace_" + trace2 + ".dat", ttnr, ttnc);
  
  std::cout << trace1 << " " << m1.trace() << " " << m1d.trace() << std::endl;
  std::cout << trace2 << " " << m2.trace() << " " << m2d.trace() << std::endl;

  communicator->barrier();
  
  //std::string trace1x2 = "63446.31982";
  //std::string trace1dx2d = "613.83754";
  std::string trace1x2 = "14.90212";
  std::string trace1dx2d = "3.22304";
  
  Matrix product12 = Matrix();
  Matrix product21 = Matrix();
  Matrix product12d = Matrix();
  Matrix product21d = Matrix();
  Matrix product1d2 = Matrix();
  Matrix product2d1 = Matrix();
  Matrix product1d2d = Matrix();
  Matrix product2d1d = Matrix();
  
  product12.multiply(&m1, 'N', &m2, 'N');
  product21.multiply(&m2, 'N', &m1, 'N');
  product12d.multiply(&m1, 'N', &m2d, 'N');
  product21d.multiply(&m2, 'N', &m1d, 'N');
  product1d2.multiply(&m1d, 'N', &m2, 'N');
  product2d1.multiply(&m2d, 'N', &m1, 'N');
  product1d2d.multiply(&m1d, 'N', &m2d, 'N');
  product2d1d.multiply(&m2d, 'N', &m1d, 'N');
  
  
  misc.message << std::endl;
  misc.message << trace1x2 << " " << m1.traceOfMatrixProduct(&m2) << " " << m2.traceOfMatrixProduct(&m1) << std::endl;
  misc.message << trace1dx2d << " " << m1d.traceOfMatrixProduct(&m2) << " " << m1d.traceOfMatrixProduct(&m2d) << " " << m1.traceOfMatrixProduct(&m2d) << std::endl;
  misc.message << trace1dx2d << " " << m2d.traceOfMatrixProduct(&m1) << " " << m2d.traceOfMatrixProduct(&m1d) << " " << m2.traceOfMatrixProduct(&m1d) << std::endl;
  misc.message << trace1x2 << " " << product12.trace() << " " << product21.trace() << std::endl;
  misc.message << trace1dx2d << " " << product12d.trace() << " " << product21d.trace() << " " << product1d2.trace() << " " << product2d1.trace() << std::endl;
  misc.message << trace1dx2d << " " << product1d2d.trace() << " " << product2d1d.trace() << std::endl;
  
  /*for(int i = 0; i < communicator->mpiNumTasks; i++) //Just for testing broadasting values works properly.
  {
    std::stringstream temp;
    temp << std::endl;
    temp << trace1x2 << " " << m1.traceOfMatrixProduct(&m2) << " " << m2.traceOfMatrixProduct(&m1) << std::endl;
    temp << trace1dx2d << " " << m1d.traceOfMatrixProduct(&m2) << " " << m1d.traceOfMatrixProduct(&m2d) << " " << m1.traceOfMatrixProduct(&m2d) << std::endl;
    temp << trace1dx2d << " " << m2d.traceOfMatrixProduct(&m1) << " " << m2d.traceOfMatrixProduct(&m1d) << " " << m2.traceOfMatrixProduct(&m1d) << std::endl;
    temp << trace1x2 << " " << product12.trace() << " " << product21.trace() << std::endl;
    temp << trace1dx2d << " " << product12d.trace() << " " << product21d.trace() << " " << product1d2.trace() << " " << product2d1.trace() << std::endl;
    temp << trace1dx2d << " " << product1d2d.trace() << " " << product2d1d.trace() << std::endl;
    if(communicator->mpiRank == i)
    {
      std::cout << temp.str();
    }
    communicator->barrier();
  }*/
  
}

void test50()
{
  //communicator->nDefaultBlockRows = 3;
  //communicator->nDefaultBlockCols = 3;
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 2;
  
  int nr = 10, nc = 10;
  //int nr = 11, nc = 11;
  double * gm;
  if(communicator->mpiRoot)
  {
    gm = new double [nr*nc];
    for(int r = 0; r<nr; r++)
    {
      for(int c = 0; c<nc; c++)
      {
        if(c>r)
        {
          gm[c*nr + r] = -1.;
        }
        else
        {
          gm[c*nr + r] = 1.;
        }
        if(c == r)
        {
          //gm[c*nr + r] = r+1; //Check null hypothesis
          gm[c*nr + r] = r*100 + c;
        }
      }
    }
  }
  
  Matrix m(cyclicDistribution, nr, nc);
  m.scatterMatrix(gm);
  m.showGlobal("m");
  std::vector<double> diagonal = m.diagonal();
  
  Matrix md(diagonalDistribution, nr, nc);
  md.setDiagonal(&(diagonal[0]), diagonal.size());
  md.showGlobal("md");
  std::vector<double> diagonalD = md.diagonal();
  
  misc.message << diagonal.size() << " " << diagonalD.size() << std::endl << "--" << std::endl;
  
  for(int id = 0; id<diagonal.size(); id++)
  {
    misc.message << "d: " <<  diagonal[id] << " " << diagonalD[id] << std::endl;
    if(diagonal[id] != double(id*100. + id) || diagonalD[id] != double(id*100. + id))
    {
      misc.error("\n------------------------\nError: diagonal differ!\n------------------------\n", 0);
    }
  }
  
  if(diagonal.size() != nr && communicator->mpiRoot == true)
  {
    misc.error("\n------------------------\nError: diagonal dimensions differ!\n------------------------\n", 0);
  }
  
  Matrix mt(diagonalDistribution, 1, 1);
  Matrix mtt(cyclicDistribution, 876, 272);
  mt.transpose(&m);
  mt.showGlobal("mt");
  mtt.transpose(&mt);
  mtt.showGlobal("mtt");
  compareMatrices(&m, &mtt);
  
  Matrix mdt(cyclicDistribution, 876, 12);
  mdt.uplo = 'L'; //Just for testing with gdb it properly changes to 'B'
  Matrix mdtt(cyclicDistribution, 1, 2);
  mdt.transpose(&md);
  mdt.showGlobal("mdt");
  mdtt.transpose(&mdt);
  mdtt.showGlobal("mdtt");
  compareMatrices(&md, &mdtt);
  
  mdt.symmetrizeTriangularMatrix();
  
  if(communicator->mpiRoot)
  {
    delete [] gm;
  }
}

void test51()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 3;
  
  Matrix m1(cyclicDistribution, 3, 4, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m1.fillWithConstant(1.1);
  Matrix m2(cyclicDistribution, 5, 2, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  m2.fillWithConstant(2.2);
  
  Matrix m3(diagonalDistribution, 3, 3, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  if(communicator->mpiRoot)
  {
    for(int i = 0; i<3; i++)
    {
      m3.m[i] = i + 0.4;
    }
  }
  Matrix m4(diagonalDistribution, 5, 5, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  if(communicator->mpiRoot)
  {
    for(int i = 0; i<5; i++)
    {
      m4.m[i] = i + 0.5;
    }
  }
  
  Matrix r(cyclicDistribution);
  
  
  m1.showGlobal("m1");
  m2.showGlobal("m2");
  m3.showGlobal("m3");
  m4.showGlobal("m4");
  
  /*r.joinMatrices(&m1, subMatrix(1,1,m1.nGlobRows,m1.nGlobCols), &m2, subMatrix(3+2,4+2,m2.nGlobRows,m2.nGlobCols), 3.3);
  r.showGlobal("r");
  r.joinMatrices(&m1, subMatrix(1,1,m1.nGlobRows,m1.nGlobCols), &m2, subMatrix(3,4,m2.nGlobRows,m2.nGlobCols), 3.3);
  r.showGlobal("r");
  r.joinMatrices(&m1, subMatrix(3,4,m1.nGlobRows,m1.nGlobCols), &m2, subMatrix(2,1,m2.nGlobRows,m2.nGlobCols), 3.3);
  r.showGlobal("r");*/
  
  /*r.joinMatrices(&m4, subMatrix(3,3,m4.nGlobRows,m4.nGlobCols), &m3, subMatrix(4,4,m3.nGlobRows,m3.nGlobCols), 3.3);
  r.showGlobal("r");
  
  r.joinMatrices(&m4, subMatrix(3,3,m4.nGlobRows,m4.nGlobCols), &m3, subMatrix(4,4,m3.nGlobRows,m3.nGlobCols), 0.0);
  r.showGlobal("r");
  
  r.joinMatrices(&m4, subMatrix(3,6,m4.nGlobRows,m4.nGlobCols), &m3, subMatrix(4,4,m3.nGlobRows,m3.nGlobCols), 0.0);
  r.showGlobal("r");
  
  r.joinMatrices(&m4, subMatrix(3,3,m4.nGlobRows,m4.nGlobCols), &m3, subMatrix(4,1,m3.nGlobRows,m3.nGlobCols), 0.0);
  r.showGlobal("r");
  
  r.joinMatrices(&m4, subMatrix(3,6,m4.nGlobRows,m4.nGlobCols), &m3, subMatrix(4,4,m3.nGlobRows,m3.nGlobCols), 1.1);
  r.showGlobal("r");
  
  r.joinMatrices(&m4, subMatrix(3,3,m4.nGlobRows,m4.nGlobCols), &m3, subMatrix(4,1,m3.nGlobRows,m3.nGlobCols), 1.1);
  r.showGlobal("r");*/
  
  r.joinMatrices(&m1, subMatrix(3,3,m1.nGlobRows,m1.nGlobCols), &m3, subMatrix(4,4,m3.nGlobRows,m3.nGlobCols), 0.0);
  r.showGlobal("r");
  
  r.joinMatrices(&m3, subMatrix(3,3,m3.nGlobRows,m3.nGlobCols), &m1, subMatrix(2,2,m1.nGlobRows,m1.nGlobCols), 0.0);
  r.showGlobal("r");
  
  r.joinMatrices(&m1, subMatrix(2,2,m1.nGlobRows,m1.nGlobCols), &m3, subMatrix(4,4,m3.nGlobRows,m3.nGlobCols), 3.3);
  r.showGlobal("r");
  
  r.joinMatrices(&m4, subMatrix(1,1,m4.nGlobRows,m4.nGlobCols), &m3, subMatrix(9,9,m3.nGlobRows,m3.nGlobCols), 0.0);
  r.showGlobal("r");
  
  r.joinMatrices(&m4, subMatrix(3,3,m4.nGlobRows,m4.nGlobCols), &m3, subMatrix(1,1,m3.nGlobRows,m3.nGlobCols), 0.0);
  r.showGlobal("r");
}

void test52()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 3;
  
  Matrix m = Matrix(cyclicDistribution);
  m.debugRead("testMatrices/matrix11.dat", 11, 7);
  
  Matrix mD = Matrix(diagonalDistribution, 10, 10, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  if(communicator->mpiRoot)
  {
    for(int i = 0; i < mD.nGlobRows; i++)
    {
      mD.m[i] = i;
    }
  }
  
  std::vector< std::vector<double> > gm;
  std::vector< double > gmflat;
  std::vector< std::vector<double> > gmd;
  std::vector< double > gmdflat;
  m.matrixToStandardVector(gm);
  m.matrixToStandardVector(gmflat);
  mD.matrixToStandardVector(gmd);
  mD.matrixToStandardVector(gmdflat);

  m.showGlobal("m");
  
  if(communicator->mpiRoot)
  {
    for(int r = 0; r < gm.size(); r++)
    {
      for(int c = 0; c < gm[0].size(); c++)
      {
        misc.message << gm[r][c] << " ";
      }
      misc.message << std::endl;
    }
    misc.message << std::endl;
    if(gmflat.size() != m.nGlobRows*m.nGlobCols)
    {
      misc.error("ERROR, dimensions do not match!", 0);
    }
    for(int r = 0; r < gmflat.size(); r++)
    {
      misc.message << gmflat[r] << " ";
    }
    misc.message << std::endl;
    misc.message << std::endl;
  }
  
  mD.showGlobal("mD");
  
  if(communicator->mpiRoot)
  {
    for(int r = 0; r < gmd.size(); r++)
    {
      for(int c = 0; c < gmd[0].size(); c++)
      {
        misc.message << gmd[r][c] << " ";
      }
      misc.message << std::endl;
    }
    misc.message << std::endl;
    if(gmdflat.size() != mD.nGlobRows*mD.nGlobCols)
    {
      misc.error("ERROR, dimensions do not match!", 0);
    }
    for(int r = 0; r < gmdflat.size(); r++)
    {
      misc.message << gmdflat[r] << " ";
    }
    misc.message << std::endl;
    misc.message << std::endl;
  }
      
}

void test53()
{
  Genotype * genotype = new Genotype(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;
  
  Matrix * originalGRM = new Matrix(grm->kernel);
  originalGRM->symmetrizeTriangularMatrix();
  
  //grm->printGRM();
  grm->kernel->showGlobal("grm", false);
  grm->diagonalizeKernel();
  //compareMatrices(originalGRM, grm->eigenVectors, 1e-2);
  grm->printGRM();
  grm->recoverKernelFromEigenDecomposition();
  //grm->printGRM();
  grm->kernel->showGlobal("grm", false);
  
  compareMatrices(originalGRM, grm->kernel, 1e-2);
}

void test53bis()
{
  Genotype * genotype = new Genotype(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;
  
  Matrix * originalGRM = new Matrix(grm->kernel);
  originalGRM->symmetrizeTriangularMatrix();
  
  //grm->printGRM();
  grm->kernel->showGlobal("grm", false);
  grm->diagonalizeKernel();
  //compareMatrices(originalGRM, grm->eigenVectors, 1e-2);
  grm->printGRM();
  grm->recoverKernelFromEigenDecomposition();
  //grm->printGRM();
  grm->kernel->showGlobal("grm", false);
  
  compareMatrices(originalGRM, grm->kernel, 1e-2);
}

void test54()
{
  communicator->nDefaultBlockRows = 3;
  communicator->nDefaultBlockCols = 3;
  
  Genotype * genotype = new Genotype(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;
  
  Kernel * originalGRM = new Kernel(grm);
  
  grm->writeGRM("tmp/original");
  grm->diagonalizeKernel();
  grm->writeGRM("tmp/diagonalized");
  
  Kernel * normalGRM = new Kernel("tmp/original");
  Kernel * diagonalGRM = new Kernel("tmp/diagonalized");
  
  originalGRM->kernel->showGlobal("original GRM", false);
  normalGRM->kernel->showGlobal("read GRM", false);
  
  originalGRM->N->showGlobal("original N", false);
  normalGRM->N->showGlobal("read N", false);
  
  grm->printGRM();
  diagonalGRM->printGRM();
  
  misc.message << "Checking read grm" << std::endl;
  originalGRM->kernel->symmetrizeTriangularMatrix();
  normalGRM->kernel->symmetrizeTriangularMatrix();
  compareMatrices(originalGRM->kernel, normalGRM->kernel);
  misc.message << "Checking read N" << std::endl;
  originalGRM->N->symmetrizeTriangularMatrix();
  normalGRM->N->symmetrizeTriangularMatrix();
  compareMatrices(originalGRM->N, normalGRM->N);
  
  misc.message << "Checking read eigenVectors" << std::endl;
  compareMatrices(grm->eigenVectors, diagonalGRM->eigenVectors);
  misc.message << "Checking read eigenValues" << std::endl;
  compareMatrices(grm->eigenValues, diagonalGRM->eigenValues);
  
  misc.message << "Checking original grm recovery from eigenValues/eigenVectors." << std::endl;
  diagonalGRM->recoverKernelFromEigenDecomposition();
  compareMatrices(originalGRM->kernel, diagonalGRM->kernel, 1e-2);
}

void test54bis()
{
  communicator->nDefaultBlockRows = 3;
  communicator->nDefaultBlockCols = 3;
  
  Genotype * genotype = new Genotype(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;
  
  Kernel * originalGRM = new Kernel(grm);
  
  grm->writeKernel("tmp/original");
  grm->diagonalizeKernel();
  grm->writeKernel("tmp/diagonalized");
  
  Kernel * normalGRM = new Kernel("tmp/original");
  Kernel * diagonalGRM = new Kernel("tmp/diagonalized");
  
  originalGRM->kernel->showGlobal("original GRM", false);
  normalGRM->kernel->showGlobal("read GRM", false);
  
  originalGRM->N->showGlobal("original N", false);
  normalGRM->N->showGlobal("read N", false);
  
  grm->printGRM();
  diagonalGRM->printGRM();
  
  misc.message << "Checking read grm" << std::endl;
  originalGRM->kernel->symmetrizeTriangularMatrix();
  normalGRM->kernel->symmetrizeTriangularMatrix();
  compareMatrices(originalGRM->kernel, normalGRM->kernel);
  misc.message << "Checking read N" << std::endl;
  originalGRM->N->symmetrizeTriangularMatrix();
  normalGRM->N->symmetrizeTriangularMatrix();
  compareMatrices(originalGRM->N, normalGRM->N);
  
  misc.message << "Checking read eigenVectors" << std::endl;
  compareMatrices(grm->eigenVectors, diagonalGRM->eigenVectors);
  misc.message << "Checking read eigenValues" << std::endl;
  compareMatrices(grm->eigenValues, diagonalGRM->eigenValues);
  
  misc.message << "Checking original grm recovery from eigenValues/eigenVectors." << std::endl;
  diagonalGRM->recoverKernelFromEigenDecomposition();
  compareMatrices(originalGRM->kernel, diagonalGRM->kernel, 1e-2);
}

void test55()
{
  Genotype * genotype = new Genotype(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;
  
  //grm->denormalize();
  
  grm->kernel->symmetrizeTriangularMatrix();
  grm->N->symmetrizeTriangularMatrix();
  
  Kernel * originalGRM = new Kernel(grm);

  compareMatrices(originalGRM->kernel, grm->kernel);
  compareMatrices(originalGRM->N, grm->N);
  
  grm->diagonalizeKernel();
  
  Kernel * diagonalizedGRM = new Kernel(grm);
  
  compareMatrices(diagonalizedGRM->eigenValues, grm->eigenValues);
  compareMatrices(diagonalizedGRM->eigenVectors, grm->eigenVectors);
  
  grm->recoverKernelFromEigenDecomposition();
  compareMatrices(originalGRM->kernel, grm->kernel, 1e-5);
  misc.message << "GRM current pointers: grm = " << grm->kernel << " N = " << grm->N << std::endl;
  
  Matrix test = Matrix();
  std::vector<int> vtest;
  //All the following functions should raise an error: (Now, not allshould raise an error)
  //diagonalizedGRM->eigenValues->filterRowsAndColumns(&test, vtest, vtest);
  //diagonalizedGRM->eigenValues->eigenDecomposition(&test, &test);
  //diagonalizedGRM->normalize();
  //diagonalizedGRM->denormalize();
  //grm->denormalize();
  //grm->normalize(); //This is the only which should not raise an error because is already normalized;
  //grm->normalized = false; //If we change the flag, then it raises an error.
  //grm->normalize();
}

void test55bis()
{
  Genotype * genotype = new Genotype(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;
  
  //grm->denormalize();
  
  grm->kernel->symmetrizeTriangularMatrix();
  grm->N->symmetrizeTriangularMatrix();
  
  Kernel * originalGRM = new Kernel(grm);

  compareMatrices(originalGRM->kernel, grm->kernel);
  compareMatrices(originalGRM->N, grm->N);
  
  grm->diagonalizeKernel();
  
  Kernel * diagonalizedGRM = new Kernel(grm);
  
  compareMatrices(diagonalizedGRM->eigenValues, grm->eigenValues);
  compareMatrices(diagonalizedGRM->eigenVectors, grm->eigenVectors);
  
  grm->recoverKernelFromEigenDecomposition();
  compareMatrices(originalGRM->kernel, grm->kernel, 1e-5);
  misc.message << "GRM current pointers: grm = " << grm->kernel << " N = " << grm->N << std::endl;
  
  Matrix test = Matrix();
  std::vector<int> vtest;
  //All the following functions should raise an error: (Now, not allshould raise an error)
  //diagonalizedGRM->eigenValues->filterRowsAndColumns(&test, vtest, vtest);
  //diagonalizedGRM->eigenValues->eigenDecomposition(&test, &test);
  //diagonalizedGRM->normalize();
  //diagonalizedGRM->denormalize();
  //grm->denormalize();
  //grm->normalize(); //This is the only which should not raise an error because is already normalized;
  //grm->normalized = false; //If we change the flag, then it raises an error.
  //grm->normalize();
}

void test56()
{
  Phenotype phenotype = Phenotype(MATRIX_DEFAULT_DISTRIBUTION, "./test/testColumns.phenos", 2);
  
  std::cout << phenotype.nPhenotypesInFile << std::endl;
}

void test57()
{
  Genotype *genotype = new Genotype(options.genotypeFile);
  Kernel *grm1 = new Kernel(genotype);
  Kernel *grm2 = new Kernel(genotype);
  Kernel *k1 = new Kernel(genotype);
  Kernel *k2 = new Kernel(genotype);
  delete genotype;
  
  grm2->kernel->fillDiagonal(1.12, 0.37);
  grm2->kernel->uplo = 'L';
  Matrix *temp = new Matrix(grm2->N);
  temp->fillDiagonal(1., 3.);
  grm2->N->add(temp);
  delete temp;
  
  k2->kernel->fillDiagonal(1.12, 0.37);
  k2->kernel->uplo = 'L';
  temp = new Matrix(k2->N);
  temp->fillDiagonal(1., 3.);
  k2->N->add(temp);
  delete temp;
  
  
  /*grm1->kernel->showGlobal("grm1");
  grm1->N->showGlobal("N1");
  grm2->kernel->showGlobal("grm2", false);
  grm2->N->showGlobal("N2");*/
  
  misc.message << "grm1:" << std::endl;
  for(int i = 0; i<grm1->randomVarNames.size(); i++)
  {
    misc.message << grm1->randomVarNames[i] << std::endl;
  }
  misc.message << "grm2:" << std::endl;
  for(int i = 0; i<grm2->randomVarNames.size(); i++)
  {
    misc.message << grm2->randomVarNames[i] << std::endl;
  }
  
  Kernel *grm3;
  Kernel *grm4;
  Kernel *k3;
  Kernel *k4;
  
  ////////////////////////////////////////////////////////////
  
  grm3 = new Kernel();
  k3 = new Kernel();
  if(communicator->mpiRoot)
  {
    grm2->randomVarNames.clear();
    grm2->randomVarNames.push_back("rs6N");
    grm2->randomVarNames.push_back("rs3N");
    //grm2->randomVarNames.push_back("rs3");
    k2->randomVarNames.clear();
    k2->randomVarNames.push_back("rs6N");
    k2->randomVarNames.push_back("rs3N");
  }
  //grm3->addGRMs(-1., grm1, 1., grm2);
  grm3->addGRMs(1., grm1, 1., grm2);
  k3->addKernels(1., k1, 1., k2);

  grm4 = new Kernel(grm2);
  grm4->addGRMs(1., grm1);
  k4 = new Kernel(k2);
  k4->addKernels(1., k1);
  test57aux(grm3, grm4, k3, k4);
  
  grm4 = new Kernel(grm1);
  grm4->addGRMs(1., grm2);
  k4 = new Kernel(k1);
  k4->addKernels(1., k2);
  test57aux(grm3, grm4, k3, k4);
  
  ////////////////////////////////////////////////////////////
  
  grm3 = new Kernel();
  k3 = new Kernel();
  if(communicator->mpiRoot)
  {
    grm2->randomVarNames.clear();
    grm2->randomVarNames.push_back("rs6");
    grm2->randomVarNames.push_back("rs3");
    //grm2->randomVarNames.push_back("rs3N");
    k2->randomVarNames.clear();
    k2->randomVarNames.push_back("rs6");
    k2->randomVarNames.push_back("rs3");
  }
  
  grm3->addGRMs(1., grm1, -1., grm2);
  k3->addKernels(1., k1, -1., k2);
  
  grm4 = new Kernel(grm2);
  grm4->addGRMs(1., grm1, -1.);
  k4 = new Kernel(k2);
  k4->addKernels(1., k1, -1.);
  test57aux(grm3, grm4, k3, k4);

  grm4 = new Kernel(grm1);
  grm4->addGRMs(-1., grm2);
  k4 = new Kernel(k1);
  k4->addKernels(-1., k2);
  test57aux(grm3, grm4, k3, k4);
  
  
  ////////////////////////////////////////////////////////////
  
  grm3 = new Kernel();
  k3 = new Kernel();
  if(communicator->mpiRoot)
  {
    grm2->randomVarNames.clear();
    grm2->randomVarNames.push_back("rs6");
    grm2->randomVarNames.push_back("rs3"); //If I comment this line, an error is raised as expected. It works properly.
    grm2->randomVarNames.push_back("rs7");
    grm2->randomVarNames.push_back("rs5");
    grm2->randomVarNames.push_back("rs4");
    grm2->randomVarNames.push_back("rs2");
    grm2->randomVarNames.push_back("rs1");
    grm2->randomVarNames.push_back("rs1N");
    grm2->randomVarNames.push_back("rs2N");
    k2->randomVarNames.clear();
    k2->randomVarNames.push_back("rs6");
    k2->randomVarNames.push_back("rs3"); //If I comment this line, an error is raised as expected. It works properly.
    k2->randomVarNames.push_back("rs7");
    k2->randomVarNames.push_back("rs5");
    k2->randomVarNames.push_back("rs4");
    k2->randomVarNames.push_back("rs2");
    k2->randomVarNames.push_back("rs1");
    k2->randomVarNames.push_back("rs1N");
    k2->randomVarNames.push_back("rs2N");
  }
  
  grm3->addGRMs(-1., grm1, 1., grm2);
  k3->addKernels(-1., k1, 1., k2);
  
  grm4 = new Kernel(grm2);
  grm4->addGRMs(-1., grm1);
  k4 = new Kernel(k2);
  k4->addKernels(-1., k1);
  test57aux(grm3, grm4, k3, k4);

  grm4 = new Kernel(grm1);
  grm4->addGRMs(1., grm2, -1.);
  k4 = new Kernel(k1);
  k4->addKernels(1., k2, -1.);
  test57aux(grm3, grm4, k3, k4);
  
  ////////////////////////////////////////////////////////////
  
  grm3 = new Kernel();
  k3 = new Kernel();
  if(communicator->mpiRoot)
  {
    grm2->randomVarNames.clear();
    grm2->randomVarNames.push_back("rs6N");
    grm2->randomVarNames.push_back("rs3N");
    //grm2->randomVarNames.push_back("rs3");
    k2->randomVarNames.clear();
    k2->randomVarNames.push_back("rs6N");
    k2->randomVarNames.push_back("rs3N");
  }
  
  grm3->addGRMs(-1., grm1, -1., grm2);
  k3->addKernels(-1., k1, -1., k2);
  
  grm4 = new Kernel(grm2);
  grm4->addGRMs(-1., grm1, -1.);
  k4 = new Kernel(k2);
  k4->addKernels(-1., k1, -1.);
  test57aux(grm3, grm4, k3, k4);

  grm4 = new Kernel(grm1);
  grm4->addGRMs(-1., grm2, -1.);
  k4 = new Kernel(k1);
  k4->addKernels(-1., k2, -1.);
  test57aux(grm3, grm4, k3, k4);

  
  ////////////////////////////////////////////////////////////
  
  
  std::vector<std::string> fIds;
  if(communicator->mpiRoot)
  {
    fIds.push_back("i2@i2");
    fIds.push_back("i3@i3");
    fIds.push_back("i6@i6");
  }
  
  std::vector<std::string> fIds2;
  if(communicator->mpiRoot)
  {
    fIds2.push_back("i1@i1");
    fIds2.push_back("i5@i5");
  }
  grm4 = new Kernel(grm1);
  grm4->filterIndividuals(fIds);
  k4 = new Kernel(k1);
  k4->filterIndividuals(fIds);
  
  test57aux(grm4, grm4, k4, k4);
  
  grm4 = new Kernel(grm1);
  grm4->filterIndividualsAsymmetric(fIds, fIds2);
  k4 = new Kernel(k1);
  k4->filterIndividualsAsymmetric(fIds, fIds2);
  
  test57aux(grm4, grm4, k4, k4, false);
  
  grm4 = new Kernel(grm1);
  grm4->diagonalizeKernel();
  k4 = new Kernel(k1);
  k4->diagonalizeKernel();
  
  compareMatrices(grm4->eigenValues, k4->eigenValues);
  compareMatrices(grm4->eigenVectors, k4->eigenVectors);
  
  misc.message << "\n********************************************\n" << std::endl;
}

void test57aux(Kernel *g1, Kernel *g2, Kernel *k1, Kernel *k2, bool symmetrize)
{
  if(symmetrize == true)
  {
    g2->kernel->symmetrizeTriangularMatrix();
    g2->N->symmetrizeTriangularMatrix();
    g1->kernel->symmetrizeTriangularMatrix();
    g1->N->symmetrizeTriangularMatrix();
    
    k2->kernel->symmetrizeTriangularMatrix();
    k2->N->symmetrizeTriangularMatrix();
    k1->kernel->symmetrizeTriangularMatrix();
    k1->N->symmetrizeTriangularMatrix();
  }
  
  compareMatrices(g1->kernel, g2->kernel, 1e-10);
  compareMatrices(g1->N, g2->N, 1e-10);
  
  compareMatrices(g1->kernel, k1->kernel);
  compareMatrices(g1->N, k1->N);
  compareMatrices(g2->kernel, k2->kernel);
  compareMatrices(g2->N, k2->N);
  compareMatrices(g1->kernel, k2->kernel, 1e-10);
  compareMatrices(g1->N, k2->N, 1e-10);
  
  std::set<std::string> sSNP1(g1->randomVarNames.begin(), g1->randomVarNames.end());
  std::set<std::string> sSNP2(g2->randomVarNames.begin(), g2->randomVarNames.end());
  
  if(g1->randomVarNames.size() != g2->randomVarNames.size() || sSNP1 != sSNP2 || g1->randomVarNames != k1->randomVarNames || g2->randomVarNames != k2->randomVarNames)
  {
    misc.error("Error: Discordant SNPs.", 0);
  }
  
  /*g1->kernel->showGlobal("grm3");
  g1->N->showGlobal("N grm3");
  
  k1->kernel->showGlobal("k1");
  k1->N->showGlobal("N k1");
  
  k2->kernel->showGlobal("k2");
  k2->N->showGlobal("N k2");
  
  misc.message << g1->randomVarNames.size() << " " << g2->randomVarNames.size() << " " << k1->randomVarNames.size() << " " << k2->randomVarNames.size() << std::endl;
  for(int i = 0; i<g1->randomVarNames.size(); i++)
  {
    misc.message << g1->randomVarNames[i] << " " << k1->randomVarNames[i] << " " << k2->randomVarNames[i] << std::endl;
  }*/
  
  misc.message << "GRMs and SNPs equal to Kernels and randomVarNames" << std::endl;
}

void test58()
{
  //mpirun --np 6 ./dissect --debug-vars
  //mpirun --np 6 ./dissect --bfile ./test/test3_withcoords
  
  Genotype *genotype = new Genotype(options.genotypeFile);
  
  std::vector<std::string> SNPIdssubset;
  //SNPIdssubset.push_back("rs7");
  //SNPIdssubset.push_back("rs5");
  //SNPIdssubset.push_back("rs1");
  //SNPIdssubset.push_back("rs4");
  
  //options.fixedGroupSize = 1;
  //options.fixedGroupSize = 2;
  options.fixedGroupSize = 3;
  genotype->groupSNPs(byOrderedFixedSize, SNPIdssubset);
  
  genotype->printGenotype(true);
  
  delete genotype;
}

void test59()
{
  communicator->nDefaultBlockRows = 3;
  communicator->nDefaultBlockCols = 5;
  
  for(int rep = 0; rep<10; rep++)
  {
    int dimension = 17;
    std::vector<int> categories;
    std::vector<double> categories2;
    for(int i = 0; i<dimension; i++)
    {
      if(rand() % 1000 < 500)
      {
        categories.push_back(1);
        categories2.push_back(1.); //This must be comented for the NULL test
      }
      else
      {
        categories.push_back(0);
        categories2.push_back(0); //This must be comented for the NULL test
      }
      //As a NULL test:
      /*if(rand() % 1000 < 500)
      {
        categories2.push_back(1.);
      }
      else
      {
        categories2.push_back(0.);
      }*/
    }
    
    Matrix * m = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, dimension, dimension);
    misc.message << "Blocks: " << m->nBlockRows << " " << m->nBlockCols << std::endl;
    m->scatterVector(&(categories2[0]), row);
    int *vt = m->scatterVectorRet(&(categories[0]), row);
    for(int i = 0; i<m->nRows; i++)
    {
      if(int(m->v[i]) != vt[i])
      {
        misc.error("Error: Vectors differ", 0);
      }
      //std::cout << m->v[i] << " " << vt[i] << std::endl;
    }
    delete [] vt;
    
    m->scatterVector(&(categories2[0]), column);
    vt = m->scatterVectorRet(&(categories[0]), column);
    for(int i = 0; i<m->nCols; i++)
    {
      if(int(m->v[i]) != vt[i])
      {
        misc.error("Error: Vectors differ", 0);
      }
      //std::cout << m->v[i] << " " << vt[i] << std::endl;
    }
    delete [] vt;
    
    delete m;
    
    communicator->nDefaultBlockRows = rand() % 4 + 1;
    communicator->nDefaultBlockCols = rand() % 4 + 1;
  }
  
}

void test60()
{
  communicator->nDefaultBlockRows = 3;
  communicator->nDefaultBlockCols = 5;
  
  //std::vector<double> diagonal;
  //std::vector<double> diagonalTest;
  
  for(int rep = 0; rep<10; rep++)
  {
    double intersection = 5. - float(rand() % 1000)/100.;
    double background = 5. - float(rand() % 1000)/100.;
    
    int dimension = 7;
    std::vector<int> categories;
    std::vector<double> categories2;
    //diagonal.clear();
    //diagonalTest.clear();
    misc.message << "vector: ";
    for(int i = 0; i<dimension; i++)
    {
      //diagonal.push_back(2);
      int var = rand() % 1000;
      if(var < 500)
      {
        categories.push_back(1);
        categories2.push_back(0);
        categories2.push_back(1);
        categories2.push_back(0);
        //diagonalTest.push_back(1);
        misc.message << "1 ";
      }
      else if (var > 500 && var < 800)
      {
        categories.push_back(0);
        categories2.push_back(1);
        categories2.push_back(0);
        categories2.push_back(0);
        //diagonalTest.push_back(1);
        misc.message << "0 ";
      }
      else
      {
        categories.push_back(2);
        categories2.push_back(0);
        categories2.push_back(0);
        categories2.push_back(1);
        //diagonalTest.push_back(2); //This is the value used for missings.
        misc.message << "2 ";
      }
    }
    misc.message << std::endl;
    
    Matrix * m = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, dimension, dimension);
    misc.message << "Blocks: " << m->nBlockRows << " " << m->nBlockCols << std::endl;
    //makeIntersectionMatrix(categories, double onIntersection = 1., double onNoIntersection = 0.);
    //m->makeIntersectionMatrix(categories, -2.3, 0.5);
    m->makeIntersectionMatrix(categories, intersection, background);
    //m->makeIntersectionMatrix(categories, 1., 0.);
    
    /*
    m->fillWithConstant(0.);
    m->setDiagonal(&(diagonal[0]), diagonal.size());
    //m->makeIntersectionMatrix(categories);
    m->makeIntersectionMatrix(categories, 1., true, 2);
    */ //deprecated
    
   
        
    
    Matrix * gm = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, 3, dimension);
    Matrix * test = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    gm->scatterMatrix(&(categories2[0]));
    test->multiply(gm, 'T', gm, 'N');
    test->symmetrizeTriangularMatrix();
    
    Matrix * temp = new Matrix(test);
    temp->fillWithConstant(1.);
    temp->add(test, 1., -1.);
    test->add(temp, intersection, background);
    
    
    //test->setDiagonal(&(diagonal[0]), diagonal.size()); //As expected this fails //deprecated
    //test->setDiagonal(&(diagonalTest[0]), diagonalTest.size()); //As expected this fails //deprecated

    gm->showGlobal("Test", false);
    
    m->showGlobal("Intersection");
    test->showGlobal("Test", false);
    
    compareGlobalMatrices(test, m);
    
    
    delete m;
    delete test;
    delete gm;
    delete temp;
    
    communicator->nDefaultBlockRows = rand() % 4 + 1;
    communicator->nDefaultBlockCols = rand() % 4 + 1;
  }
  
}

void test61()
{
  for(int idb = 1; idb < 13; idb++ )
  {
    for(int nr = 1; nr<50; nr++)
    {
      for(int nc = 1; nc<50; nc++)
      {
        communicator->nDefaultBlockRows = 2;
        communicator->nDefaultBlockCols = 2;
        
        //int nr = 5, nc = 3;
        int minimumDimension = (nc<nr)?nc:nr;
        double * gm;
        if(communicator->mpiRoot)
        {
          gm = new double [nr*nc];
          for(int r = 0; r<nr; r++)
          {
            for(int c = 0; c<nc; c++)
            {
              if(c>r)
              {
                gm[c*nr + r] = -1.;
              }
              else
              {
                gm[c*nr + r] = 1.;
              }
              if(c == r)
              {
                //gm[c*nr + r] = r+1; //Check null hypothesis
                gm[c*nr + r] = r*100 + c;
              }
            }
          }
        }
        
        Matrix m(cyclicDistribution, nr, nc);
        m.scatterMatrix(gm);
        //m.showGlobal("m");
        std::vector<double> diagonal = m.diagonal(true);
        
        Matrix md(diagonalDistribution, minimumDimension, minimumDimension);
        md.setDiagonal(&(diagonal[0]), diagonal.size());
        //md.showGlobal("md");
        std::vector<double> diagonalD = md.diagonal();
        
        //misc.message << diagonal.size() << " " << diagonalD.size() << std::endl << "--" << std::endl;
        
        //diagonal[0] = 3.3; //Null test
        //diagonal.push_back(4.4); //Null test
        //diagonal.pop_back();
        //diagonal.erase(diagonal.begin()+5);
        //diagonal.erase(diagonal.begin()+1);
        //diagonal.erase(diagonal.begin());
        
        for(int id = 0; id<diagonal.size(); id++)
        {
          //misc.message << "d: " <<  diagonal[id] << " " << diagonalD[id] << std::endl;
          if(diagonal[id] != double(id*100. + id) || diagonalD[id] != double(id*100. + id))
          {
            misc.error("\n------------------------\nError: diagonal differ!\n------------------------\n", 0);
          }
        }
        
        
        if(diagonal.size() != minimumDimension && communicator->mpiRoot == true)
        {
          misc.error("\n------------------------\nError: diagonal dimensions differ!\n------------------------\n", 0);
        }
        
        if(communicator->mpiRoot)
        {
          delete [] gm;
        }
      }
    }
  }
}

void test62()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 2;
  
  //int nr = 7, nc = 13; //Tested and does not work as expected.
  int nr = 9, nc = 7;
  std::vector<int> dependent;
  dependent.push_back(3);
  dependent.push_back(6);
  double * gm;
  if(communicator->mpiRoot)
  {
    gm = new double [nr*nc];
    for(int r = 0; r<nr; r++)
    {
      for(int c = 0; c<nc; c++)
      {
        if(c>r)
        {
          gm[c*nr + r] = -1.;
          gm[c*nr + r] = rand() % 1000;
        }
        else
        {
          gm[c*nr + r] = 1.;
          gm[c*nr + r] = rand() % 1000;
        }
        if(c == r)
        {
          //gm[c*nr + r] = r+1; //Check null hypothesis
          gm[c*nr + r] = r*100 + c;
          gm[c*nr + r] = rand() % 1000;
        }
      }
    }
    
    for(int sc = 0; sc<dependent.size(); sc++)
    {
      int dc = dependent[sc];
      double factor = (rand() % 1000);
      for(int r = 0; r<nr; r++)
      {
        gm[dc*nr + r] = gm[sc*nr + r]*factor;
        //gm[dc*nr + r] = gm[sc*nr + r]*factor + (rand() % 1); //For adding noise
        //gm[dc*nr + r] = gm[sc*nr + r]*factor + (rand() % 2); //For adding noise
      }
    }
  }
  
  Matrix m(cyclicDistribution, nr, nc);
  m.scatterMatrix(gm);
  //m.symmetric = true; //Just for testing
  //m.uplo = 'L'; //Just for testing
  m.showGlobal("m");
  std::vector<int> dependentColumns = m.getDependentColumns();
  m.showGlobal("QR");
  
  for(int id = 0; id<dependentColumns.size(); id++)
  {
    misc.message << "dependent: " <<  dependentColumns[id] << std::endl;
  }
  
  if(communicator->mpiRoot)
  {
    delete [] gm;
  } 
}

void test63()
{
  //mpirun --np 6 ./dissect --debug-vars
  //mpirun --np 6 ./dissect --debug-vars --extract  ./test/test3.snps
  
  Genotype *genotype = new Genotype(options.genotypeFile);
  
  std::vector<std::string> SNPIdssubset;
  options.fixedGroupSize = 3;
  genotype->groupSNPs(byOrderedFixedSize, SNPIdssubset);
  
  genotype->printGenotype(true);
  
  delete genotype;
  
  options.genotypeFile = "test/parts.snps/merged";
  options.genotypeListFile = "test/parts.snps/list.dat";

  Genotype * genotype1 = loadGenotypeUsingOptions();
  
  options.genotypeFile = "";
  options.genotypeListFile = "test/parts.snps/list.dat";
  
  Genotype * genotype2 = loadGenotypeUsingOptions();

  compareMatrices(genotype1->genotypes, genotype2->genotypes);
  compareMatrices(genotype1->missings, genotype2->missings);
  
  genotype1->normalizeGenotypes();
  genotype2->normalizeGenotypes();
  
  compareMatrices(genotype1->genotypes, genotype2->genotypes);
  compareMatrices(genotype1->missings, genotype2->missings);
}

void test64()
{
  communicator->nDefaultBlockCols = 3;
  communicator->nDefaultBlockRows = 2;
  
  Genotype * genotype = new Genotype(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;
  
  std::vector<double> temp;
  
  Kernel * grmUpper = new Kernel(grm);
  grmUpper->kernel->transpose(grm->kernel);
  grmUpper->N->transpose(grm->N);
  grmUpper->kernel->matrixToStandardVector(temp);
  if(communicator->mpiRoot == true)  {temp[1] = 2.;}
  grmUpper->kernel->scatterMatrix(&(temp[0]));
  std::vector<std::string> individualIdsUpper = grmUpper->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  Kernel * grmLower = new Kernel(grm);
  grmLower->kernel->matrixToStandardVector(temp);
  if(communicator->mpiRoot == true)  {temp[69] = 2.;}
  grmLower->kernel->scatterMatrix(&(temp[0]));
  std::vector<std::string> individualIdsLower = grmLower->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  Kernel * grmBoth = new Kernel(grm);
  grmBoth->kernel->symmetrizeTriangularMatrix();
  grmBoth->N->symmetrizeTriangularMatrix();
  std::vector<std::string> individualIdsBoth = grmBoth->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  //double lowerT = -0.8;
  //double upperT = 0.44;
  double lowerT = -0.87;
  double upperT = 1.3;
  
  std::vector<int> gri1;
  std::vector<int> gci1;
  //grmUpper->filterIndividuals(individualIdsUpper);
  grmUpper->kernel->getGlobalIndexOutsideRange(lowerT, upperT, gri1, gci1);
  grmUpper->kernel->showGlobal("grmUpper", false);
  std::vector< std::vector<double> > gUpper;
  grmUpper->kernel->matrixToStandardVector(gUpper);
  
  std::vector<int> gri2;
  std::vector<int> gci2;
  //grmLower->filterIndividuals(individualIdsLower);
  grmLower->kernel->getGlobalIndexOutsideRange(lowerT, upperT, gri2, gci2);
  grmLower->kernel->showGlobal("grmLower", false);
  std::vector< std::vector<double> > gLower;
  grmLower->kernel->matrixToStandardVector(gLower);
  
  std::vector<int> gri3;
  std::vector<int> gci3;
  //grmBoth->filterIndividuals(individualIdsBoth);
  grmBoth->kernel->getGlobalIndexOutsideRange(lowerT, upperT, gri3, gci3);
  std::vector< std::vector<double> > gBoth;
  grmBoth->kernel->matrixToStandardVector(gBoth);
  
  for(int i = 0; i<communicator->mpiNumTasks; i++)
  {
    if(communicator->mpiRank == i)
    {
      std::cout << i << " " << gri1.size() << " " << gci1.size() << std::endl;
      std::cout << i << " " << gri2.size() << " " << gci2.size() << std::endl;
      std::cout << i << " " << gri3.size() << " " << gci3.size() << std::endl;
    }
    communicator->barrier();
  }
  
  if(gri1.size() != gci1.size() || gri1.size() != gri2.size() || gri1.size() != gri3.size() || gri1.size() != gci2.size() || gri1.size() != gci3.size())
  {
    misc.error("Error: Error in indices.", 0);
  }
  
  std::set< std::pair<int, int> > s1;
  std::set< std::pair<int, int> > s2;
  std::set< std::pair<int, int> > s3;
  
  if(communicator->mpiRoot == true)
  {
    for(int i = 0; i<gri1.size(); i++)
    {
      std::cout << gri1[i] << " " << gci1[i];
      std::cout << " " << gri2[i] << " " << gci2[i];
      std::cout << " " << gri3[i] << " " << gci3[i] << std::endl;
      
      std::pair<int, int> p1(gri1[i], gci1[i]);
      std::pair<int, int> p2(gri2[i], gci2[i]);
      std::pair<int, int> p3(gri3[i], gci3[i]);
      
      if(s1.find(p1) != s1.end() || s2.find(p2) != s2.end() || s3.find(p3) != s3.end())
      {
        misc.error("Error: Repeated indices.", 0);
      }
      
      s1.insert(p1);
      s2.insert(p2);
      s3.insert(p3);
    }
    
    if(s1!=s2 || s1 != s3)
    {
      misc.error("Error: Difering indices.", 0);
    }
  }

  if(communicator->mpiRoot == true)
  {
    if(gBoth.size() != gLower.size() || gBoth.size() != gUpper.size())
    {
      misc.error("Error: Something is wrong.", 0);
    }
    std::cout << "s1 size: " << s1.size() << std::endl;
    for(int i = 0; i< gBoth.size(); i++)
    {
      for(int j = 0; j< gBoth.size(); j++)
      {
        if( ( gBoth[i][j] > upperT || gBoth[i][j] < lowerT ) ) /*|| (gLower[i][j] > options.grmCutoff && j<i) || (gUpper[i][j] > options.grmCutoff && i<j)*/
        {
          std::pair<int, int> p1(i, j);
          if(s1.find(p1) != s1.end())
          { 
            std::cout << gBoth[i][j] << " properly removed" << std::endl;
            s1.erase(p1);
          }
          else
          {
            misc.error("Error: Something is wrong. Index undetected", 0);
          }
        }
      }
    }
    std::cout << "s1 size: " << s1.size() << std::endl;
    if(s1.size() != 0)
    {
      misc.error("Error: Something is wrong. Detected more indexs than expected.", 0);
    }
  }
  
}

void test65()
{
  communicator->nDefaultBlockCols = 3;
  communicator->nDefaultBlockRows = 2;
  
  Genotype * genotype = new Genotype(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;
  
  std::vector<double> temp;
  
  Kernel * grmUpper = new Kernel(grm);
  grmUpper->kernel->transpose(grm->kernel);
  grmUpper->N->transpose(grm->N);
  grmUpper->kernel->matrixToStandardVector(temp);
  if(communicator->mpiRoot == true)  {temp[1] = 2.;}
  grmUpper->kernel->scatterMatrix(&(temp[0]));
  std::vector<std::string> individualIdsUpper = grmUpper->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  Kernel * grmLower = new Kernel(grm);
  grmLower->kernel->matrixToStandardVector(temp);
  if(communicator->mpiRoot == true)  {temp[69] = 2.;}
  grmLower->kernel->scatterMatrix(&(temp[0]));
  std::vector<std::string> individualIdsLower = grmLower->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  Kernel * grmBoth = new Kernel(grm);
  grmBoth->kernel->symmetrizeTriangularMatrix();
  grmBoth->N->symmetrizeTriangularMatrix();
  std::vector<std::string> individualIdsBoth = grmBoth->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  //double lowerT = -0.87;
  //double upperT = 1.3;
  double lowerT = 1.2;
  double upperT = 2.9;
  //double upperT = 1.8;

  
  std::vector<int> gri1;
  std::vector<int> gci1;
  //grmUpper->filterIndividuals(individualIdsUpper);
  grmUpper->kernel->getGlobalIndexInsideRange(lowerT, upperT, gri1, gci1);
  grmUpper->kernel->showGlobal("grmUpper", false);
  std::vector< std::vector<double> > gUpper;
  grmUpper->kernel->matrixToStandardVector(gUpper);
  
  std::vector<int> gri2;
  std::vector<int> gci2;
  //grmLower->filterIndividuals(individualIdsLower);
  grmLower->kernel->getGlobalIndexInsideRange(lowerT, upperT, gri2, gci2);
  grmLower->kernel->showGlobal("grmLower", false);
  std::vector< std::vector<double> > gLower;
  grmLower->kernel->matrixToStandardVector(gLower);
  
  std::vector<int> gri3;
  std::vector<int> gci3;
  //grmBoth->filterIndividuals(individualIdsBoth);
  grmBoth->kernel->getGlobalIndexInsideRange(lowerT, upperT, gri3, gci3);
  std::vector< std::vector<double> > gBoth;
  grmBoth->kernel->matrixToStandardVector(gBoth);
  
  for(int i = 0; i<communicator->mpiNumTasks; i++)
  {
    if(communicator->mpiRank == i)
    {
      std::cout << i << " " << gri1.size() << " " << gci1.size() << std::endl;
      std::cout << i << " " << gri2.size() << " " << gci2.size() << std::endl;
      std::cout << i << " " << gri3.size() << " " << gci3.size() << std::endl;
    }
    communicator->barrier();
  }
  
  if(gri1.size() != gci1.size() || gri1.size() != gri2.size() || gri1.size() != gri3.size() || gri1.size() != gci2.size() || gri1.size() != gci3.size())
  {
    misc.error("Error: Error in indices.", 0);
  }
  
  std::set< std::pair<int, int> > s1;
  std::set< std::pair<int, int> > s2;
  std::set< std::pair<int, int> > s3;
  
  if(communicator->mpiRoot == true)
  {
    for(int i = 0; i<gri1.size(); i++)
    {
      std::cout << gri1[i] << " " << gci1[i];
      std::cout << " " << gri2[i] << " " << gci2[i];
      std::cout << " " << gri3[i] << " " << gci3[i] << std::endl;
      
      std::pair<int, int> p1(gri1[i], gci1[i]);
      std::pair<int, int> p2(gri2[i], gci2[i]);
      std::pair<int, int> p3(gri3[i], gci3[i]);
      
      if(s1.find(p1) != s1.end() || s2.find(p2) != s2.end() || s3.find(p3) != s3.end())
      {
        misc.error("Error: Repeated indices.", 0);
      }
      
      s1.insert(p1);
      s2.insert(p2);
      s3.insert(p3);
    }
    
    if(s1!=s2 || s1 != s3)
    {
      misc.error("Error: Difering indices.", 0);
    }
  }

  if(communicator->mpiRoot == true)
  {
    if(gBoth.size() != gLower.size() || gBoth.size() != gUpper.size())
    {
      misc.error("Error: Something is wrong.", 0);
    }
    std::cout << "s1 size: " << s1.size() << std::endl;
    for(int i = 0; i< gBoth.size(); i++)
    {
      for(int j = 0; j< gBoth.size(); j++)
      {
        if( ( gBoth[i][j] < upperT && gBoth[i][j] > lowerT ) ) /*|| (gLower[i][j] > options.grmCutoff && j<i) || (gUpper[i][j] > options.grmCutoff && i<j)*/
        {
          std::pair<int, int> p1(i, j);
          if(s1.find(p1) != s1.end())
          { 
            std::cout << gBoth[i][j] << " properly removed" << std::endl;
            s1.erase(p1);
          }
          else
          {
            misc.error("Error: Something is wrong. Index undetected", 0);
          }
        }
      }
    }
    std::cout << "s1 size: " << s1.size() << std::endl;
    if(s1.size() != 0)
    {
      misc.error("Error: Something is wrong. Detected more indexs than expected.", 0);
    }
  }
  
}


void test66()
{
  //Tested with
  
  //--grm-cutoff 0.4
  //--grm-cutoff 0.5
  //--grm-cutoff 0.025
  
  options.maximumProportionOfElimitaedIndividuals = 0.8;
  //options.maximumProportionOfElimitaedIndividuals = 0.01;
  
  misc.message << "Prune: " << options.pruneGRM << " GRM cutoff: " << options.grmCutoff << std::endl;
  misc.message << "GRM min overlaping snps: " << options.minimumNumberOverlapingSNPs << std::endl;
  misc.message << "Allowed proportion of eliminated individuals: " << options.maximumProportionOfElimitaedIndividuals << std::endl;
  
  
  
  communicator->nDefaultBlockCols = 3;
  communicator->nDefaultBlockRows = 2;
  
  Genotype * genotype = new Genotype(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;
  
  std::vector<double> temp;
  
  Kernel * grmUpper = new Kernel(grm);
  grmUpper->kernel->transpose(grm->kernel);
  grmUpper->N->transpose(grm->N);
  grmUpper->kernel->matrixToStandardVector(temp);
  if(communicator->mpiRoot == true)  {temp[1] = 2.;}
  grmUpper->kernel->scatterMatrix(&(temp[0]));
  grmUpper->kernel->showGlobal("grmUpper", false);
  //std::vector<std::string> individualIdsUpper = grmUpper->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  Kernel * grmLower = new Kernel(grm);
  grmLower->kernel->matrixToStandardVector(temp);
  if(communicator->mpiRoot == true)  {temp[69] = 2.;}
  grmLower->kernel->scatterMatrix(&(temp[0]));
  grmLower->kernel->showGlobal("grmLower", false);
  //std::vector<std::string> individualIdsLower = grmLower->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  Kernel * grmBoth = new Kernel(grm);
  grmBoth->kernel->symmetrizeTriangularMatrix();
  grmBoth->N->symmetrizeTriangularMatrix();
  grmBoth->kernel->showGlobal("grmBoth", false);
  //std::vector<std::string> individualIdsBoth = grmBoth->searchNoHighRelatedIndividuals(options.grmCutoff);
  
  //grmUpper->filterIndividuals(individualIdsUpper);
  misc.message << grmUpper->sanitizeKernel() << std::endl;
  //grmUpper->kernel->showGlobal("grmUpper", false);
  grmUpper->printGRM();
  std::vector< std::vector<double> > gUpper;
  grmUpper->kernel->matrixToStandardVector(gUpper);
  
  //grmLower->filterIndividuals(individualIdsLower);
  misc.message << grmLower->sanitizeKernel() << std::endl;
  //grmLower->kernel->showGlobal("grmLower", false);
  grmLower->printGRM();
  std::vector< std::vector<double> > gLower;
  grmLower->kernel->matrixToStandardVector(gLower);
  
  //grmBoth->filterIndividuals(individualIdsBoth);
  misc.message << grmBoth->sanitizeKernel() << std::endl;
  grmBoth->printGRM();
  std::vector< std::vector<double> > gBoth;
  grmBoth->kernel->matrixToStandardVector(gBoth);
  
    
//   if(communicator->mpiRoot == true)
//   {
//     if(gBoth.size() != gLower.size() || gBoth.size() != gUpper.size() || gBoth.size() != individualIdsBoth.size())
//     {
//       misc.error("Error: Something is wrong.", 0);
//     }
//     for(int i = 0; i< individualIdsBoth.size(); i++)
//     {
//       for(int j = 0; j< individualIdsBoth.size(); j++)
//       {
//         if( (gBoth[i][j] > options.grmCutoff && i!=j) || (gLower[i][j] > options.grmCutoff && j<i) || (gUpper[i][j] > options.grmCutoff && i<j))
//         {
//           misc.error("Error: Not properly pruned.", 0);
//         }
//       }
//     }
//   }
  
  
  delete grm;
  
  
//   //Test with a larger GRM
//   
//   options.grmFile = "";
//   options.genotypeFile = "test/parts.snps/merged";
//   options.genotypeListFile = "";
//   options.GRMJoinMethod = 0;
//   GRM* grmLarge = loadGRMUsingOptionsOld();
//   grmLarge->pruneKernel(options.grmCutoff);
//   std::vector< std::vector<double> > gLarge;
//   grmLarge->grm->matrixToStandardVector(gLarge);
//   
//   //grmLarge->grm->showGlobal("grmLarge", false);
//   
//   if(communicator->mpiRoot == true)
//   {
//     if(gLarge.size() != grmLarge->individualIds.size() || grmLarge->grm->uplo != 'L' )
//     {
//       misc.error("Error: Something is wrong.", 0);
//     }
//     for(int i = 0; i< grmLarge->individualIds.size(); i++)
//     {
//       for(int j = 0; j< grmLarge->individualIds.size(); j++)
//       {
//         if( (gLarge[i][j] > options.grmCutoff && j<i /*j!=i*/) )
//         //if( (gLarge[i][j] > options.grmCutoff && /*j<i*/ j!=i) )
//         {
//           std::cout << grmLarge->grm->uplo << " " << i << " " << j << " " << gLarge[i][j] << std::endl;
//           misc.error("Error: Not properly pruned.", 0);
//         }
//       }
//     }
//   }
}

void test67()
{
  std::srand(communicator->mpiRank*std::time(0));

  communicator->nDefaultBlockRows = 3; //With mpirun --np 4 and this block size, this is buggy.
  communicator->nDefaultBlockCols = 3;
  
  //int dimensionBlockMatrix = 2; //Number of blocks on each axis.
  //int rDimensionBlocks = 1;
  //int cDimensionBlocks = 1;
  
  //int dimensionBlockMatrix = 5; 
  //int rDimensionBlocks = 4;
  //int cDimensionBlocks = 4;
  
  int dimensionBlockMatrix = 5; 
  int rDimensionBlocks = 4;
  int cDimensionBlocks = 4;
  
  std::vector< std::vector<Matrix*> > blockm(dimensionBlockMatrix, std::vector<Matrix*>(dimensionBlockMatrix) );
  std::vector< std::vector<Matrix*> > blockmd(dimensionBlockMatrix, std::vector<Matrix*>(dimensionBlockMatrix) );
  
  double * tgm;
  double * tgmd;
  if(communicator->mpiRoot)
  {
    tgm = new double [rDimensionBlocks*cDimensionBlocks*dimensionBlockMatrix*dimensionBlockMatrix];
    tgmd = new double [rDimensionBlocks*rDimensionBlocks*dimensionBlockMatrix*dimensionBlockMatrix];
  }

  for(int br = 0; br < dimensionBlockMatrix; br++ )
  {
    for(int bc = 0; bc < dimensionBlockMatrix; bc++ )
    {
      double * gm;
      double * gmd;
      if(communicator->mpiRoot)
      {
        gm = new double [rDimensionBlocks*cDimensionBlocks];
        gmd = new double [rDimensionBlocks];
        for(int r = 0; r<rDimensionBlocks; r++)
        {
          for(int c = 0; c<cDimensionBlocks; c++)
          {
            if(c>=r)
            {
              gm[c*rDimensionBlocks + r] = (rand() % 20) + 3;
            }
            else
            {
              gm[c*rDimensionBlocks + r] = gm[r*rDimensionBlocks + c];
            }
            
            if(bc>=br)
            {
              tgm[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = gm[c*rDimensionBlocks + r];
            }
            else
            {
              tgm[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = tgm[ (br*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (bc*rDimensionBlocks + r) ];
            }
            
            tgmd[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = 0.;
            if(c == r)
            {
              gmd[r] = (rand() % 20) + 3;
              if(bc>=br)
              {
                tgmd[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = gmd[r];
              }
              else
              {
                tgmd[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = tgmd[ (br*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (bc*rDimensionBlocks + r) ];
              }
            }
          }
        }
      }
      
      Matrix * m = new Matrix(cyclicDistribution, rDimensionBlocks, cDimensionBlocks);
      m->scatterMatrix(gm);
      m->symmetric = true;
      m->uplo='B';
      if(bc>=br)
      {
        blockm[br][bc] = m;
      }
      else
      {
        blockm[br][bc] = blockm[bc][br];
      }
      
      Matrix * md = new Matrix(diagonalDistribution, rDimensionBlocks, rDimensionBlocks);
      md->setDiagonal(gmd, rDimensionBlocks);
      if(bc>=br)
      {
        blockmd[br][bc] = md;
      }
      else
      {
        blockmd[br][bc] = blockmd[bc][br];
      }
      
      //m->showGlobal(i2s(br) + "-" + i2s(bc));
      //md->showGlobal();
    }
  }
    
  Matrix * mTotal = new Matrix(cyclicDistribution, rDimensionBlocks*dimensionBlockMatrix, cDimensionBlocks*dimensionBlockMatrix);
  mTotal->scatterMatrix(tgm);
  mTotal->showGlobal();
  
  Matrix * mdTotal = new Matrix(cyclicDistribution, rDimensionBlocks*dimensionBlockMatrix, rDimensionBlocks*dimensionBlockMatrix);
  mdTotal->scatterMatrix(tgmd);
  mdTotal->showGlobal();
  Matrix * mdTotalInv = new Matrix(cyclicDistribution, rDimensionBlocks*dimensionBlockMatrix, rDimensionBlocks*dimensionBlockMatrix);
  mdTotalInv->scatterMatrix(tgmd);
  mdTotalInv->showGlobal();

  BlockMatrix mBlocks(blockm);
  BlockMatrix mdBlocks(blockmd);
  BlockMatrix mdBlocksInv(blockmd);
        
  //mBlocks.invert();
  mdBlocks.showGlobal();
  Matrix * test1 = mdBlocks.block2distributed();
  compareMatrices(test1, mdTotal, 1e-8);
  
  mdTotalInv->invert();  
  mdTotalInv->showGlobal("Global");
  mdBlocksInv.invert();
  mdBlocksInv.showGlobal();
  
  Matrix * test2 = mdBlocksInv.block2distributed();
  compareMatrices(test2, mdTotalInv, 1e-8);
  
  Matrix * identity1 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  identity1->multiply(mdTotal, 'N', mdTotalInv, 'N');
  identity1->showGlobal("Identity1", false, 2, 1e-5);
  
  Matrix * identity2 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  identity2->multiply(test1, 'N', test2, 'N');
  identity2->showGlobal("Identity2", false, 2, 1e-5);
  
  delete test1;
  delete test2;
  
  if(communicator->mpiRoot)
  {
    delete [] tgm;
    delete [] tgmd;
  }

}

void test67bis()
{
  for(int b = 1; b < 15; b++ )
  {
    for(int dim = 1; dim < 15; dim++ )
    {
      for(int defblock = 1; defblock < 11; defblock++ )
      {
        test67bis2(b, dim, dim, defblock, defblock);
      } 
    }
  }
}

void test67bis2(int nBlocks, int nRows, int nCols, int defBlockRow, int defBlockCol)
{
  std::srand(communicator->mpiRank*std::time(0));

  //communicator->nDefaultBlockRows = 3;
  //communicator->nDefaultBlockCols = 3;
  
  //int dimensionBlockMatrix = 2; //Number of blocks on each axis.
  //int rDimensionBlocks = 1;
  //int cDimensionBlocks = 1;
  
  //int dimensionBlockMatrix = 5; 
  //int rDimensionBlocks = 4;
  //int cDimensionBlocks = 4;
  
  //int dimensionBlockMatrix = 7; 
  //int rDimensionBlocks = 4;
  //int cDimensionBlocks = 4;
  
  communicator->nDefaultBlockRows = defBlockRow;
  communicator->nDefaultBlockCols = defBlockCol;
  
  int dimensionBlockMatrix = nBlocks; 
  int rDimensionBlocks = nRows;
  int cDimensionBlocks = nCols;
  
  std::vector< std::vector<Matrix*> > blockm(dimensionBlockMatrix, std::vector<Matrix*>(dimensionBlockMatrix) );
  std::vector< std::vector<Matrix*> > blockmd(dimensionBlockMatrix, std::vector<Matrix*>(dimensionBlockMatrix) );
  
  double * basicBlock;
  double * tgm;
  double * tgmd;
  if(communicator->mpiRoot)
  {
    tgm = new double [rDimensionBlocks*cDimensionBlocks*dimensionBlockMatrix*dimensionBlockMatrix];
    tgmd = new double [rDimensionBlocks*rDimensionBlocks*dimensionBlockMatrix*dimensionBlockMatrix];
    basicBlock = new double [rDimensionBlocks*rDimensionBlocks];
  }

  if(communicator->mpiRoot)
  {
    for(int r = 0; r<rDimensionBlocks; r++)
    {
      for(int c = 0; c<cDimensionBlocks; c++)
      {
        if(c==r)
        {
          basicBlock[c*rDimensionBlocks + r] = (rand() % 600) + 3;
        }
        else
        {
          basicBlock[c*rDimensionBlocks + r] = 0.;
        }
      }
    }
  }
  
  double factor = 1.;
  for(int br = 0; br < dimensionBlockMatrix; br++ )
  {
    for(int bc = 0; bc < dimensionBlockMatrix; bc++ )
    {
      double * gm;
      double * gmd;
      if(bc != br){factor = 0.1;} else {factor = 1.;}
      if(communicator->mpiRoot)
      {
        gm = new double [rDimensionBlocks*cDimensionBlocks];
        gmd = new double [rDimensionBlocks];
        for(int r = 0; r<rDimensionBlocks; r++)
        {
          for(int c = 0; c<cDimensionBlocks; c++)
          {
            if(c>=r)
            {
              gm[c*rDimensionBlocks + r] = basicBlock[c*rDimensionBlocks + r]*factor;
            }
            else
            {
              gm[c*rDimensionBlocks + r] = gm[r*rDimensionBlocks + c];
            }
            
            if(bc>=br)
            {
              tgm[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = gm[c*rDimensionBlocks + r];
            }
            else
            {
              tgm[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = tgm[ (br*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (bc*rDimensionBlocks + r) ];
            }
            
            tgmd[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = 0.;
            if(c == r)
            {
              gmd[r] = basicBlock[c*rDimensionBlocks + r]*factor;
              if(bc>=br)
              {
                tgmd[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = gmd[r];
              }
              else
              {
                tgmd[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = tgmd[ (br*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (bc*rDimensionBlocks + r) ];
              }
            }
          }
        }
      }
      
      Matrix * m = new Matrix(cyclicDistribution, rDimensionBlocks, cDimensionBlocks);
      m->scatterMatrix(gm);
      m->symmetric = true;
      m->uplo='B';
      if(bc>=br)
      {
        blockm[br][bc] = m;
      }
      else
      {
        blockm[br][bc] = blockm[bc][br];
      }
      
      Matrix * md = new Matrix(diagonalDistribution, rDimensionBlocks, rDimensionBlocks);
      md->setDiagonal(gmd, rDimensionBlocks);
      if(bc>=br)
      {
        blockmd[br][bc] = md;
      }
      else
      {
        blockmd[br][bc] = blockmd[bc][br];
      }
      
      //m->showGlobal(i2s(br) + "-" + i2s(bc));
      //md->showGlobal();
    }
  }
    
  double logDet1 = 0;
  double logDet2 = 0;
  double logDet3 = 0;
    
  Matrix * mTotal = new Matrix(cyclicDistribution, rDimensionBlocks*dimensionBlockMatrix, cDimensionBlocks*dimensionBlockMatrix);
  mTotal->scatterMatrix(tgm);
  //mTotal->showGlobal();
  
  Matrix * mdTotal = new Matrix(cyclicDistribution, rDimensionBlocks*dimensionBlockMatrix, rDimensionBlocks*dimensionBlockMatrix);
  mdTotal->scatterMatrix(tgmd);
  //mdTotal->showGlobal();
  Matrix * mdTotalInv = new Matrix(cyclicDistribution, rDimensionBlocks*dimensionBlockMatrix, rDimensionBlocks*dimensionBlockMatrix);
  mdTotalInv->scatterMatrix(tgmd);
  //mdTotalInv->showGlobal();
  Matrix * mdTotalInv2 = new Matrix(cyclicDistribution, rDimensionBlocks*dimensionBlockMatrix, rDimensionBlocks*dimensionBlockMatrix);
  mdTotalInv2->scatterMatrix(tgmd);
  mdTotalInv2->symmetric = true;
  mdTotalInv2->uplo = 'B';
  //mdTotalInv2->showGlobal();

  BlockMatrix mBlocks(blockm);
  BlockMatrix mdBlocks(blockmd);
  BlockMatrix mdBlocksInv(blockmd);
        
  //mBlocks.invert();
  //mdBlocks.showGlobal();
  Matrix * test1 = mdBlocks.block2distributed();
  compareMatrices(test1, mdTotal, 1e-8);
  
  mdTotalInv->invert(&logDet1);  
  //mdTotalInv->showGlobal("Global");
  mdTotalInv2->symmetricInvert(&logDet2);  
  mdTotalInv2->symmetrizeTriangularMatrix();
  //mdTotalInv2->showGlobal("Global");
  mdBlocksInv.invert(&logDet3);
  //mdBlocksInv.showGlobal();
  
  Matrix * test2 = mdBlocksInv.block2distributed();
  compareMatrices(test2, mdTotalInv, 1e-8);
  compareMatrices(test2, mdTotalInv2, 1e-8);
  
  Matrix * identity1 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  identity1->multiply(mdTotal, 'N', mdTotalInv, 'N');
  //identity1->showGlobal("Identity1", false, 2, 1e-5);
  
  Matrix * identity2 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  identity2->multiply(test1, 'N', test2, 'N');
  //identity2->showGlobal("Identity2", false, 2, 1e-5);
  
  //misc.message << std::setprecision(20) << logDet1 << " " << logDet2 << " " << logDet3 << std::endl;
  
  if( (fabs(logDet1-logDet2)/logDet1) > 1e-13 || (fabs(logDet1-logDet3)/logDet1) > 1e-13 )
  {
    misc.message << "Warning, determinants differ: " << std::setprecision(20) << logDet1 << " " << logDet2 << " " << logDet3 << " " << (fabs(logDet1-logDet2)/logDet1) << " " << (fabs(logDet1-logDet3)/logDet1) << std::endl;
  }
  
  delete test1;
  delete test2;
  
  if(communicator->mpiRoot)
  {
    delete [] tgm;
    delete [] tgmd;
    delete [] basicBlock;
  }

}

void test68()
{
  test68bis(3, 2, 2, 3, 3, 4);
  for(int b = 1; b < 15; b++ )
  {
    for(int dim = 1; dim < 15; dim++ )
    {
      for(int defblockr = 1; defblockr < 11; defblockr++ )
      {
        for(int defblockc = 1; defblockc < 11; defblockc++ )
        {
          for(int coldimtest = 1; coldimtest < 7; coldimtest++ )
          {
            test68bis(b, dim, dim, defblockr, defblockc, coldimtest);
          }
        } 
      }
    }
  }
}

void test68bis(int nBlocks, int nRows, int nCols, int defBlockRow, int defBlockCol, int colDimensionTestM)
{
  std::srand(communicator->mpiRank*std::time(0));

  //communicator->nDefaultBlockRows = 3;
  //communicator->nDefaultBlockCols = 3;
  
  //int dimensionBlockMatrix = 2; //Number of blocks on each axis.
  //int rDimensionBlocks = 1;
  //int cDimensionBlocks = 1;
  
  //int dimensionBlockMatrix = 5; 
  //int rDimensionBlocks = 4;
  //int cDimensionBlocks = 4;
  
  //int dimensionBlockMatrix = 7; 
  //int rDimensionBlocks = 4;
  //int cDimensionBlocks = 4;
  
  communicator->nDefaultBlockRows = defBlockRow;
  communicator->nDefaultBlockCols = defBlockCol;
  
  int dimensionBlockMatrix = nBlocks; 
  int rDimensionBlocks = nRows;
  int cDimensionBlocks = nCols;
  
  std::vector< std::vector<Matrix*> > blockm(dimensionBlockMatrix, std::vector<Matrix*>(dimensionBlockMatrix) );
  std::vector< std::vector<Matrix*> > blockmd(dimensionBlockMatrix, std::vector<Matrix*>(dimensionBlockMatrix) );
  
  double * columnMatrix;
  double * columnMatrix2;
  double * tgm;
  double * tgmd;
  if(communicator->mpiRoot)
  {
    tgm = new double [rDimensionBlocks*cDimensionBlocks*dimensionBlockMatrix*dimensionBlockMatrix];
    tgmd = new double [rDimensionBlocks*rDimensionBlocks*dimensionBlockMatrix*dimensionBlockMatrix];
    columnMatrix = new double [rDimensionBlocks*dimensionBlockMatrix];
    columnMatrix2 = new double [rDimensionBlocks*dimensionBlockMatrix*colDimensionTestM];
  }

  if(communicator->mpiRoot)
  {
    for(int r = 0; r<rDimensionBlocks*dimensionBlockMatrix; r++)
    {
      columnMatrix[r] = (rand() % 600) + 3;
      for(int c = 0; c<colDimensionTestM; c++)
      {
        columnMatrix2[c*(rDimensionBlocks*dimensionBlockMatrix) + r] = (rand() % 600) + 3;
      }
    }
  }
  Matrix * mColumn = new Matrix(cyclicDistribution, rDimensionBlocks*dimensionBlockMatrix, 1);
  mColumn->scatterMatrix(columnMatrix);
  Matrix * mColumn2 = new Matrix(cyclicDistribution, rDimensionBlocks*dimensionBlockMatrix, colDimensionTestM);
  mColumn2->scatterMatrix(columnMatrix2);
  
  
  
  for(int br = 0; br < dimensionBlockMatrix; br++ )
  {
    for(int bc = 0; bc < dimensionBlockMatrix; bc++ )
    {
      double * gm;
      double * gmd;
      if(communicator->mpiRoot)
      {
        gm = new double [rDimensionBlocks*cDimensionBlocks];
        gmd = new double [rDimensionBlocks];
        for(int r = 0; r<rDimensionBlocks; r++)
        {
          for(int c = 0; c<cDimensionBlocks; c++)
          {
            if(c>=r)
            {
              gm[c*rDimensionBlocks + r] = (rand() % 600) + 3 + sqrt(rand() % 600);
            }
            else
            {
              gm[c*rDimensionBlocks + r] = gm[r*rDimensionBlocks + c];
            }
            
            if(bc>=br)
            {
              tgm[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = gm[c*rDimensionBlocks + r];
            }
            else
            {
              tgm[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = tgm[ (br*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (bc*rDimensionBlocks + r) ];
            }
            
            tgmd[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = 0.;
            if(c == r)
            {
              gmd[r] = (rand() % 600) + 3 + sqrt(rand() % 600);
              if(bc>=br)
              {
                tgmd[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = gmd[r];
              }
              else
              {
                tgmd[ (bc*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (br*rDimensionBlocks + r) ] = tgmd[ (br*cDimensionBlocks + c)*(rDimensionBlocks*dimensionBlockMatrix) + (bc*rDimensionBlocks + r) ];
              }
            }
          }
        }
      }
      
      Matrix * m = new Matrix(cyclicDistribution, rDimensionBlocks, cDimensionBlocks);
      m->scatterMatrix(gm);
      m->symmetric = true;
      m->uplo='B';
      if(bc>=br)
      {
        blockm[br][bc] = m;
      }
      else
      {
        blockm[br][bc] = blockm[bc][br];
      }
      
      Matrix * md = new Matrix(diagonalDistribution, rDimensionBlocks, rDimensionBlocks);
      md->setDiagonal(gmd, rDimensionBlocks);
      if(bc>=br)
      {
        blockmd[br][bc] = md;
      }
      else
      {
        blockmd[br][bc] = blockmd[bc][br];
      }
      
      //m->showGlobal(i2s(br) + "-" + i2s(bc));
      //md->showGlobal();
    }
  }
    
   
  Matrix * mTotal = new Matrix(cyclicDistribution, rDimensionBlocks*dimensionBlockMatrix, cDimensionBlocks*dimensionBlockMatrix);
  mTotal->scatterMatrix(tgm);
  //mTotal->showGlobal();
  
  Matrix * mdTotal = new Matrix(cyclicDistribution, rDimensionBlocks*dimensionBlockMatrix, rDimensionBlocks*dimensionBlockMatrix);
  mdTotal->scatterMatrix(tgmd);
  //mdTotal->showGlobal();
  
  BlockMatrix mBlocks(blockm);
  BlockMatrix mdBlocks(blockmd);
  
  //mBlocks.showGlobal();
  //mdBlocks.showGlobal();
  
  Matrix * test1 = new Matrix();
  test1->multiply(mTotal, 'N', mColumn, 'N');
  //test1->showGlobal("test1");
  
  Matrix * test2 = new Matrix();
  test2->multiply(mTotal, 'N', mColumn2, 'N');
  //test2->showGlobal("test2");
  
  BlockMatrix btest1;
  btest1.multiply(mBlocks, mColumn);
  //btest1.showGlobal();
  Matrix * btest1g = btest1.block2distributed();
  BlockMatrix btest2;
  btest2.multiply(mBlocks, mColumn2);
  //btest2.showGlobal();
  Matrix * btest2g = btest2.block2distributed();
  
  compareMatrices(test1, btest1g, 1e-8);
  compareMatrices(test2, btest2g, 1e-8);    
  
  //Same with diagonal matrices
  Matrix * dtest1 = new Matrix();
  dtest1->multiply(mdTotal, 'N', mColumn, 'N');
  //dtest1->showGlobal("dtest1");
  
  Matrix * dtest2 = new Matrix();
  dtest2->multiply(mdTotal, 'N', mColumn2, 'N');
  //dtest2->showGlobal("dtest2");
  
  BlockMatrix dbtest1;
  dbtest1.multiply(mdBlocks, mColumn);
  //dbtest1.showGlobal();
  Matrix * dbtest1g = dbtest1.block2distributed();
  BlockMatrix dbtest2;
  dbtest2.multiply(mdBlocks, mColumn2);
  //dbtest2.showGlobal();
  Matrix * dbtest2g = dbtest2.block2distributed();
  
  compareMatrices(dtest1, dbtest1g, 1e-8);
  compareMatrices(dtest2, dbtest2g, 1e-8); 
  
  
  delete test1;
  delete test2;
  
  if(communicator->mpiRoot)
  {
    delete [] tgm;
    delete [] tgmd;
    delete [] columnMatrix;
    delete [] columnMatrix2;
  }

}

void test69()
{
  Genotype * genotype = new Genotype();
  std::cout << options.genotypeFile << std::endl; std::cout.flush();
  genotype->load(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  
  genotype->genotypes->showGlobal("Genotypes");
  
  grm->printGRM();
  
  Kernel * epi = new Kernel(grm, kernelEpistaticGRM, genotype);
  
  epi->printGRM();
  
  grm->kernel->showGlobal();
  epi->kernel->showGlobal();
  
  delete epi;
  delete grm;
  delete genotype;
}

void test70()
{
  std::vector< std::vector<std::string> > table;
  getTableFromFile("test/testColumns.dat", table, 2);
  for(int ir = 0; ir < communicator->mpiNumTasks; ir++)
  {
    if( ir == communicator->mpiRank )
    {
      std::cout << ir << " -----------------" << std::endl;
      for(int i = 0; i<table.size(); i++)
      {
        for(int j = 0; j<table[i].size(); j++)
        {
          std::cout << table[i][j] << " ";
        }
        std::cout << std::endl;
      }
    }
    communicator->barrier();
  }
}

void test71()
{
  //--reml --grm testREML/autosome.geneasy --pheno testREML/MERGED.phenos --covar testREML/MERGED.covar --qcovar testREML/MERGED.qcovar
  
  REML reml1;
  
  Kernel * srcGRM = loadGRMUsingOptions();
  
  std::vector<Kernel*> grms;
  srcGRM->name = "GRM";
  grms.push_back(srcGRM);
  
  std::vector<double> weights;
  weights.push_back(1.);
  
  std::vector<int> phenotypeColumns;
  phenotypeColumns.push_back(1);

  std::vector<double> heritabilities;
  heritabilities.push_back(0.5);
  
  std::vector<std::pair<std::string, std::string> > covariateFiles;
  covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFile, options.qCovarsFile));
  
  bool prepared = reml1.prepare(singleREMLType, grms, weights, phenotypeColumns, heritabilities, covariateFiles);
  if( prepared == true )
  {
    reml1.computeREML();
  }
  else
  {
    misc.message << "Sorry, a problem was happened while preparing data for performing REML. The MLM cannot be fitted. Please, check the logs." << std::endl;
  }
  
  ///////////////////////////////////////////////////////
  
  REML reml2;
  
  Kernel * srcGRM2 = loadGRMUsingOptions();
  Phenotype * phenotype = new Phenotype(cyclicDistribution, options.phenotypesFile, options.phenotypeColumn);
  Covariate * covariate = new Covariate(options.covarsFile, options.qCovarsFile, phenotype->individualIds);

  std::vector<std::string> commonIndividuals, commonIndividualsInGRMOrder;
  if(communicator->mpiRoot)
  {
    commonIndividuals = intersectionStringVectors(3, &srcGRM2->individualIds, &phenotype->individualIds, &covariate->individualIds);
    commonIndividualsInGRMOrder = orderVectorAsTemplate(srcGRM2->individualIds, commonIndividuals);
  }
  
  phenotype->filterIndividuals(commonIndividualsInGRMOrder);
  covariate->filterIndividuals(commonIndividualsInGRMOrder);
  srcGRM2->filterIndividuals(commonIndividualsInGRMOrder, false);
  std::vector<Matrix*> grms2;
  srcGRM2->name = "GRM";
  grms2.push_back(srcGRM2->kernel);
  
  std::vector<KernelType> kernelTypes;
  kernelTypes.push_back(kernelCovarianceMatrix);
  
  prepared = reml2.prepare(phenotype->phenotypes, covariate->covariates, grms2, kernelTypes, 0.5,  weights);
  if( prepared == true )
  {
    reml2.computeREML();
  }
  else
  {
    misc.message << "Sorry, a problem was happened while preparing data for performing REML. The MLM cannot be fitted. Please, check the logs." << std::endl;
  }
  
  misc.message << std::setprecision(9) << reml1.V->variances.size() << " " << reml2.V->variances.size() << std::endl;
  misc.message << std::setprecision(9) << reml1.V->variances[0].variance << " " << reml2.V->variances[0].variance << " " << (reml1.V->variances[0].variance - reml2.V->variances[0].variance) << " " << (reml1.V->variances[0].variance - reml2.V->variances[0].variance)/reml1.V->variances[0].variance << std::endl;
  misc.message << std::setprecision(9) << reml1.V->variances[1].variance << " " << reml2.V->variances[1].variance << " " << (reml1.V->variances[1].variance - reml2.V->variances[1].variance) << " " << (reml1.V->variances[1].variance - reml2.V->variances[1].variance)/reml1.V->variances[1].variance << std::endl;
}

void test72()
{
  std::vector< std::pair<double, double> > thresholds;

  thresholds.push_back( std::pair<double, double>(0.0, 0.1) );
  //thresholds.push_back( std::pair<double, double>(-110.0, 1.7) ); //This should filter all.
  thresholds.push_back( std::pair<double, double>(-110.0, 1.2) );
  thresholds.push_back( std::pair<double, double>(-0.82, 1.2) );
  thresholds.push_back( std::pair<double, double>(-0.8, 20.2) );
  thresholds.push_back( std::pair<double, double>(-110., 0.5) );
  thresholds.push_back( std::pair<double, double>(-0.8, 0.5) );
  thresholds.push_back( std::pair<double, double>(-0.82, 0.5) );
  
  for(int rep = 0; rep<thresholds.size(); rep++)
  {
    double lowerThreshold = thresholds[rep].first;
    double upperThreshold = thresholds[rep].second;
    misc.message << lowerThreshold << " " << upperThreshold << std::endl;
    
    communicator->nDefaultBlockCols = 3;
    communicator->nDefaultBlockRows = 2;
    
    Genotype * genotype = new Genotype(options.genotypeFile);
    Kernel * grm = new Kernel(genotype);
    delete genotype;
    
    std::vector<double> temp;
    
    Kernel * grmUpper = new Kernel(grm);
    grmUpper->kernel->transpose(grm->kernel);
    grmUpper->N->transpose(grm->N);
    grmUpper->kernel->matrixToStandardVector(temp);
    if(communicator->mpiRoot == true)  {temp[1] = 2.;}
    grmUpper->kernel->scatterMatrix(&(temp[0]));
    grmUpper->keepWithRelatednessOutside(lowerThreshold, upperThreshold);
    
    Kernel * grmLower = new Kernel(grm);
    grmLower->kernel->matrixToStandardVector(temp);
    if(communicator->mpiRoot == true)  {temp[69] = 2.;}
    grmLower->kernel->scatterMatrix(&(temp[0]));
    grmLower->keepWithRelatednessOutside(lowerThreshold, upperThreshold);
    
    Kernel * grmBoth = new Kernel(grm);
    grmBoth->kernel->symmetrizeTriangularMatrix();
    grmBoth->N->symmetrizeTriangularMatrix();
    grmBoth->kernel->showGlobal("Original", false);
    grmBoth->keepWithRelatednessOutside(lowerThreshold, upperThreshold);
    
    
    grmUpper->kernel->showGlobal("grmUpper", false);
    std::vector< std::vector<double> > gUpper;
    grmUpper->kernel->symmetrizeTriangularMatrix();
    grmUpper->kernel->matrixToStandardVector(gUpper);
    
    grmLower->kernel->showGlobal("grmLower", false);
    std::vector< std::vector<double> > gLower;
    grmLower->kernel->symmetrizeTriangularMatrix();
    grmLower->kernel->matrixToStandardVector(gLower);
    
    std::vector< std::vector<double> > gBoth;
    grmBoth->kernel->matrixToStandardVector(gBoth);
    
    std::vector< std::vector<double> > gAll;
    grm->kernel->symmetrizeTriangularMatrix();
    grm->kernel->matrixToStandardVector(gAll);
    
    compareMatrices(grmUpper->kernel, grmBoth->kernel);
    compareMatrices(grmLower->kernel, grmBoth->kernel);
    
    if(communicator->mpiRoot == true)
    {
      //gAll[0][1] = 10; //NULL hypothesis
      
      if(gBoth.size() != gLower.size() || gBoth.size() != gUpper.size())
      {
        misc.error("Error: Something is wrong.", 0);
      }

      for(int i = 0; i< gBoth.size(); i++) //Check if there are individuals that had to be filtered
      {
        bool indToKeep = false;
        for(int j = 0; j< gBoth.size(); j++)
        {
          
          if( ( gBoth[i][j] > upperThreshold || gBoth[i][j] < lowerThreshold ) ) /*|| (gLower[i][j] > options.grmCutoff && j<i) || (gUpper[i][j] > options.grmCutoff && i<j)*/
          {
            if( i != j)
            {
              indToKeep = true;
            }
          }
        }
        if( indToKeep == false )
        {
          misc.error("Error: keept wrong individual", 0);
        }
      }
      
      std::vector< std::vector<double> > gFiltered;
      std::vector<int> toKeep;
      
      for(int i = 0; i< gAll.size(); i++) //Look at individuals that have to be kept in the global original matrix
      {
        bool indToKeep = false;
        for(int j = 0; j< gAll.size(); j++)
        {
          
          if( ( gAll[i][j] > upperThreshold || gAll[i][j] < lowerThreshold ) ) /*|| (gLower[i][j] > options.grmCutoff && j<i) || (gUpper[i][j] > options.grmCutoff && i<j)*/
          {
            if( i != j)
            {
              indToKeep = true;
            }
          }
        }
        if( indToKeep == true )
        {
          toKeep.push_back(i);
        }
      }
      for(int i = 0; i< toKeep.size(); i++) //Remove the individuals that have to be removed in the global matrix
      {
        std::vector<double> temp;
        for(int j = 0; j< toKeep.size(); j++)
        {
          temp.push_back( gAll[ toKeep[i] ][ toKeep[j] ] );
        }
        gFiltered.push_back(temp);
      }
      
      if( gFiltered.size() != gBoth.size() || gFiltered[0].size() != gBoth[0].size() )
      {
        misc.error("Error: Something is wrong. gFiltered!=gBoth", 0);
      }
      
      for(int i = 0; i< gBoth.size(); i++) //Compare the globally filtered matrix with the distributed one
      {
        for(int j = 0; j< gBoth.size(); j++)
        {
          
          if( gFiltered[i][j] != gBoth[i][j] )
          {
            misc.message << i << " " << j << " " << gFiltered[i][j] << " " << gBoth[i][j] << std::endl;
            misc.error("Error: Something is wrong. gFiltered[i][j]!=gBoth[i][j]", 0);
          }
        }
      }
    }
  }
}

void test73()
{
  Covariate covariate("test.covareffects/covar.txt", "test.covareffects/qcovar.txt", std::vector<std::string>(), false);
  std::map<std::string, double> covariateEffects = covariate.loadEffectPrediction("test.covareffects/covareffects.txt", "test.covareffects/qcovareffects.txt");
  
  //Checked this both in loadEffectPrediction with gdb
  //std::map< int, std::map<std::string, double> > covarEffects
  //std::vector<double> qCovarEffects
  //Tested also different bad formnating files
  //Tested also the requirement of --force-use-unestimated-values
  //Tested also the requirement of --covar and --covar-effects to go together
  //Tested also the requirement of --qcovar and --qcovar-effects to go together
}

void test74()
{
  options.discreteRandomEffectsFile = "test.covareffects/covar2.txt";
  options.discreteRandomEffectsFileColumns.push_back(0);
  //options.discreteRandomEffectsFileColumns.push_back(1);
  
  Kernel * covKernel = new Kernel(options.discreteRandomEffectsFile, kernelFromDiscreteCovariates, options.discreteRandomEffectsFileColumns[0]);
  covKernel->printGRM(true, 2);
  
  
}

void test75()
{
  
  //./dissect --fix-correlation cov-1 2 1 -1.0 GRM 1 3 -0.5 --pheno-cols 3 1 2 --multi-reml --pheno test.covareffects/qcovar.txt --bfile test/test3 --random-effects test.covareffects/qcovar.txt --random-effects-cols 1 3
  //options.analysis = withoutAnalysis;
  
  for(std::vector<REMLCovarianceRestricted>::iterator it = options.restrictedCovariances.begin(); it != options.restrictedCovariances.end(); ++it)
  {
    std::cout << it->baseId << " " << it->p1 << " " << it->p2 << " " << it->correlation << std::endl;
  }
}

void test76()
{
  Matrix m = Matrix(cyclicDistribution);
  m.debugRead("matrix1.dat", 3, 2);
  m.showGlobal("Matrix1");
  
  Matrix m2 = Matrix(cyclicDistribution);
  m2.debugRead("matrix2.dat", 2, 2);
  m2.showGlobal("Matrix2");
  
  Matrix m3 = Matrix(cyclicDistribution);
  m3.diagonalOfABAt(&m, &m2);
  m3.showGlobal("Diagonal M1*M2*M1t");
  
  Matrix temp = Matrix(cyclicDistribution);
  Matrix m4 = Matrix(cyclicDistribution);
  temp.multiply(&m, 'N', &m2, 'N');
  m4.multiply(&temp, 'N', &m, 'T');
  m4.showGlobal("M1*M2*M1t");
}

void test77(int argc, char **argv)
{
  options.numberOfProcessesDefaultMPIGroup = 6; //Tested with mpirun --np 14 dissect
  //options.numberOfProcessesDefaultMPIGroup = 1; //Tested with mpirun --np 14 dissect
  
  Communicator * communicator1 = new Communicator(communicator, basicGroupedCommunicator);
  
  for(int i = 0; i<communicator->mpiNumTasks; i++)
  {
    if(i == communicator->mpiRank)
    {
      std::cout << i << " (" << communicator1->myRow << ", " << communicator1->myCol << ") " << communicator1->mpiRank << " " << communicator1->group << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  //Communicator * communicator2 = new Communicator(communicator, basicGroupedCommunicator);
  
  delete communicator1;
}

void test78(int argc, char **argv)
{
  //options.numberOfProcessesDefaultMPIGroup = 2;
  options.numberOfProcessesDefaultMPIGroup = 3;
  int nr = 2;
  int nc = 3;
  
  Communicator * communicator1 = new Communicator(communicator, basicGroupedCommunicator);
  Communicator * communicatorGlobal = communicator;
  
  double * mg1 = NULL;
  double * mg2 = NULL;
  if(communicator1->mpiRoot == true)
  {
    mg1 = new double [nr*nc];
    mg2 = new double [nr*nc];
    for(int i =0; i<nr*nc; i++)
    {
      if(communicatorGlobal->mpiRoot == true)
      {
        mg1[i] = i;
        mg2[i] = 2*i;
      }
      else
      {
        mg1[i] = i;
        mg2[i] = -2*i;
      }
    }
  }
  
  communicator = communicator1;
  
  Matrix * m1 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nr, nc);
  m1->scatterMatrix(mg1);
  Matrix * m2 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nr, nc);
  m2->scatterMatrix(mg2);
  
  m1->showGlobal();
  m2->showGlobal();
  
  Matrix * mr = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nr, nc);
  mr->multiply(m1, 'N', m2, 'T');
  
  mr->showGlobal("MR");
  
  delete m1;
  delete m2;
  delete mr;
  
  communicator = communicatorGlobal;
  
  if(communicator1->mpiRoot == true)
  {
    delete [] mg1;
    delete [] mg2;
  }
  delete communicator1;
}

void test79(int argc, char **argv)
{
  options.grmFile = "testREML/autosome.geneasy";
  options.phenotypesFile = "testREML/MERGED.phenos";
  options.covarsFile = "testREML/MERGED.covar";
  options.qCovarsFile = "testREML/MERGED.qcovar";
  options.outFile = "test/test";
  
  options.numberOfProcessesDefaultMPIGroup = 2;
  //options.numberOfProcessesDefaultMPIGroup = 3;
  //options.numberOfProcessesDefaultMPIGroup = 1;
  
  Communicator * communicator1 = new Communicator(communicator, basicGroupedCommunicator);
  Communicator * communicatorGlobal = communicator;
  communicator = communicator1;
  
  options.outFile += "." + i2s(communicatorGlobal->group);
  
  //--reml --grm testREML/autosome.geneasy --pheno testREML/MERGED.phenos --covar testREML/MERGED.covar --qcovar testREML/MERGED.qcovar
  
  REML reml1;
  
  Kernel * srcGRM = loadGRMUsingOptions();
  
  std::vector<Kernel*> grms;
  srcGRM->name = "GRM";
  grms.push_back(srcGRM);
  
  std::vector<double> weights;
  weights.push_back(1.);
  
  std::vector<int> phenotypeColumns;
  if(communicator->group == 0)
  {
    phenotypeColumns.push_back(1);
  }
  else
  {
    phenotypeColumns.push_back(2);
  }

  std::vector<double> heritabilities;
  heritabilities.push_back(0.5);
  
  std::vector<std::pair<std::string, std::string> > covariateFiles;
  covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFile, options.qCovarsFile));
  
  bool prepared = reml1.prepare(singleREMLType, grms, weights, phenotypeColumns, heritabilities, covariateFiles);
  if( prepared == true )
  {
    reml1.computeREML();
  }
  else
  {
    misc.message << "Sorry, a problem was happened while preparing data for performing REML. The MLM cannot be fitted. Please, check the logs." << std::endl;
  }
  
  communicator = communicatorGlobal;
  
  delete communicator1;
}

void test80(int argc, char **argv)
{
  options.numberOfProcessesDefaultMPIGroup = 1; //Works with mpirun --np 2
  //options.numberOfProcessesDefaultMPIGroup = 2; //Works with mpirun --np 3
  int nr = 2;
  int nc = 3;
  
  Communicator * communicator1 = new Communicator(communicator, basicGroupedCommunicator);
  Communicator * communicatorGlobal = communicator;
  
  double * mg1 = NULL;
  if(communicator1->mpiRoot == true)
  {
    mg1 = new double [nr*nc];
    for(int i =0; i<nr*nc; i++)
    {
      if(communicatorGlobal->mpiRoot == true)
      {
        mg1[i] = i;
      }
      else
      {
        mg1[i] = i;
      }
    }
  }
  
  /*std::map<int, std::pair<int, int> > rowsOriginDestination;
  rowsOriginDestination[0] = std::pair<int, int>(0,0);
  rowsOriginDestination[1] = std::pair<int, int>(0,1);
  
  std::map<int, std::pair<int, int> > colsOriginDestination;
  colsOriginDestination[0] = std::pair<int, int>(0,0);
  colsOriginDestination[1] = std::pair<int, int>(1,0);
  colsOriginDestination[2] = std::pair<int, int>(0,1);*/
  
  std::vector<int> nGlobalRowsInGroup;
  nGlobalRowsInGroup.push_back(2);
  nGlobalRowsInGroup.push_back(2);
  std::vector<int> nGlobalColsInGroup;
  nGlobalColsInGroup.push_back(2);
  nGlobalColsInGroup.push_back(1);
  
  std::map<int, std::pair<int, int> > rowsOriginDestination;
  rowsOriginDestination[0] = std::pair<int, int>(0,-1);
  rowsOriginDestination[1] = std::pair<int, int>(1,-1);
  
  std::map<int, std::pair<int, int> > colsOriginDestination;
  colsOriginDestination[0] = std::pair<int, int>(0,0);
  colsOriginDestination[1] = std::pair<int, int>(1,0);
  colsOriginDestination[2] = std::pair<int, int>(0,1);
  
  
  
  
  Matrix * m1 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nr, nc);
  m1->scatterMatrix(mg1);
  m1->showGlobal();

  
  Matrix * ms = m1->redistributionToGroupedCommunicatorMatrices(communicator1, nGlobalRowsInGroup, nGlobalColsInGroup, rowsOriginDestination, colsOriginDestination, false);
  delete m1;
  
  communicator = communicator1;
  
  ms->showGlobal("MS");
  
  delete ms;
  
  communicator = communicatorGlobal;
  
  
  if(communicator1->mpiRoot == true)
  {
    delete [] mg1;
  }
  delete communicator1;
}

void test81(int argc, char **argv)
{
  options.numberOfProcessesDefaultMPIGroup = 1; //Works with mpirun --np 2
  //options.numberOfProcessesDefaultMPIGroup = 2; //Works with mpirun --np 3
  int nr = 2;
  int nc = 3;
  
  Communicator * communicator1 = new Communicator(communicator, basicGroupedCommunicator);
  Communicator * communicatorGlobal = communicator;
  
  double * mg1 = NULL;
  if(communicator1->mpiRoot == true)
  {
    mg1 = new double [nr*nc];
    for(int i =0; i<nr*nc; i++)
    {
      if(communicatorGlobal->mpiRoot == true)
      {
        mg1[i] = i;
      }
      else
      {
        mg1[i] = i;
      }
    }
  }
  
  Matrix * m1 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nr, nc);
  m1->scatterMatrix(mg1);
  m1->showGlobal();

  
  Matrix * ms = m1->copyToGroupedCommunicator(communicator1);
  delete m1;
  
  communicator = communicator1;
  
  ms->showGlobal("MS");
  
  delete ms;
  
  communicator = communicatorGlobal;
  
  
  if(communicator1->mpiRoot == true)
  {
    delete [] mg1;
  }
  delete communicator1;
}

void test82()
{
  //Tested also when an individual in the couples is not in the GRM, or when an individual is in the grm but not in the first column of the couples.
  
  options.indirectEffectsCouplesFile = "test/couples.dat";
  
  std::vector<std::string> test;
  std::vector<std::string> test2;
  if(communicator->mpiRoot)
  {
    test.push_back("i1@i1");
    test.push_back("i2@i2");
    test.push_back("i3@i3");
    test.push_back("i4@i4");
  }
  
  Genotype *genotype = new Genotype();
  std::cout << options.genotypeFile << std::endl; std::cout.flush();
  genotype->load(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  delete genotype;

  grm->filterIndividuals(test, true);
  grm->printGRM();
  
  std::vector<Kernel*> kernels;
  std::vector<std::string> reducedModels;
  std::vector<std::string> individualBLUPNames;
  kernels.push_back(grm);
  reducedModels.push_back(grm->name);
  individualBLUPNames.push_back(grm->name);
  
  SingleREML sr = SingleREML();
  introduceResortedGRMsByCouples(kernels, reducedModels, individualBLUPNames);
  
  for(int i = 0; i<kernels.size(); i++)
  {
    misc.message << reducedModels[i] << " " << kernels[i]->name << std::endl;
    kernels[i]->printGRM();
  }
}

void test83()
{
  std::srand(communicator->mpiRank*std::time(0));
  //Matrix * m = new Matrix(cyclicDistribution, 4, 20000);
  //Matrix * m2 = new Matrix(cyclicDistribution, 40000, 20000);
  Matrix * m = new Matrix(cyclicDistribution, 40000, 2);
  Matrix * m2 = new Matrix(cyclicDistribution, 40000, 2);
  Matrix * mr = new Matrix();
  for(int i = 0; i<m->nRows*m->nCols; i++)
  {
    m->m[i] = double(rand());
    m2->m[i] = double(rand());
  }
  
  mr->multiply(m, 'N', m2, 'T');
  
  std::cout << mr->m[3];
  //delete m;
}

void test84()
{
  Genotype *genotype = new Genotype();
  genotype->load(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  grm->kernel->fillWithConstant(1.);
  
  
  options.useCorrelations = true;
  //options.useCorrelations = false;
  //options.useLogLogisticScale = true;
  //options.firstStepEM = false;
  
  options.phenotypeColumns = std::vector<int>();
  options.phenotypeColumns.push_back(1);
  options.phenotypeColumns.push_back(2);
  options.phenotypeColumns.push_back(3);
  options.phenotypeColumns.push_back(4);
  options.phenotypeColumns.push_back(5);
  options.initialh2Traits = std::vector<double>();
  options.initialh2Traits.push_back(0.2);
  options.initialh2Traits.push_back(0.2);
  options.initialh2Traits.push_back(0.2);
  options.initialh2Traits.push_back(0.2);
  options.initialh2Traits.push_back(0.2);
  options.qCovarsFiles.clear();
  options.covarsFiles.clear();
  
  options.pruneGRM = false;
  options.minimumNumberOverlapingSNPs = 0;
  /*
  if(options.useCorrelations == true)
  {
    options.initialVariancesFile = "test/test3.initialVars";
    //options.initialVariancesFile = "test/test3.v2.initialVars";
  }
  else
  {
    options.initialVariancesFile = "test/test3.covars.initialVars";
  }
  */
  
  std::vector<std::string> keepIndividualIds;
  keepIndividualIds.push_back("i1@i1");
  keepIndividualIds.push_back("i2@i2");
  keepIndividualIds.push_back("i3@i3");
  grm->filterIndividuals(keepIndividualIds);
  
  std::vector<Kernel*> kernels;
  grm->name = options.baseVarianceNames["GRM"];
  kernels.push_back(grm);
  
  std::vector<double> weights;
  weights.assign( kernels.size() , 1./double(kernels.size()) );
  
  std::vector<std::pair<std::string, std::string> > covariateFiles;
  for(int i = 0; i<options.phenotypeColumns.size(); i++)
  {
    covariateFiles.push_back(std::pair<std::string, std::string>("", ""));
  }
  
  REML *reml = new REML();
  
  bool prepared = reml->prepare(bivariateREMLType, kernels, weights, options.phenotypeColumns, options.initialh2Traits, covariateFiles);
  
  reml->V->reinitializeVariances();
  reml->V->showCovarianceMatrix(true, true);
  
  //reml->V->computeDerivateCovariance(0);
  //reml->V->mDerivate->showGlobal("Derivate 0", false);
  //reml->V->computeDerivateCovariance(1);
  //reml->V->mDerivate->showGlobal("Derivate 1", false);
  
  //reml->V->computeDerivateCovariance(5);
  //reml->V->mDerivate->showGlobal("Derivate 5", false);
  //reml->V->computeDerivateCovariance(6);
  //reml->V->mDerivate->showGlobal("Derivate 6", false);
  //reml->V->computeDerivateCovariance(7);
  //reml->V->mDerivate->showGlobal("Derivate 7", false);
  //reml->V->computeDerivateCovariance(8);
  //reml->V->mDerivate->showGlobal("Derivate 8", false);
  
  //reml->V->computeDerivateCovariance(9);
  //reml->V->mDerivate->showGlobal("Derivate 9", false);
  //reml->V->computeDerivateCovariance(10);
  //reml->V->mDerivate->showGlobal("Derivate 10", false);
  //reml->V->computeDerivateCovariance(11);
  //reml->V->mDerivate->showGlobal("Derivate 11", false);
  
  //reml->V->computeDerivateCovariance(12);
  //reml->V->mDerivate->showGlobal("Derivate 12", false);
  //reml->V->computeDerivateCovariance(13);
  //reml->V->mDerivate->showGlobal("Derivate 13", false);
  //reml->V->computeDerivateCovariance(14);
  //reml->V->mDerivate->showGlobal("Derivate 14", false);
  
  
  //reml->V->computeDerivateCovariance(15);
  //reml->V->mDerivate->showGlobal("Derivate 15", false);
  //reml->V->computeDerivateCovariance(16);
  //reml->V->mDerivate->showGlobal("Derivate 16", false);
  //reml->V->computeDerivateCovariance(17);
  //reml->V->mDerivate->showGlobal("Derivate 17", false);
  
  //reml->V->computeDerivateCovariance(24);
  //reml->V->mDerivate->showGlobal("Derivate 24", false);
  //reml->V->computeDerivateCovariance(25);
  //reml->V->mDerivate->showGlobal("Derivate 25", false);
  //reml->V->computeDerivateCovariance(26);
  //reml->V->mDerivate->showGlobal("Derivate 26", false);
  
  //reml->V->computeDerivateCovariance(27);
  //reml->V->mDerivate->showGlobal("Derivate 27", false);
  //reml->V->computeDerivateCovariance(28);
  //reml->V->mDerivate->showGlobal("Derivate 28", false);
  //reml->V->computeDerivateCovariance(29);
  //reml->V->mDerivate->showGlobal("Derivate 29", false);
  
  reml->computeREML();
  
  delete reml;
}

void test85()
{
  Genotype *genotype = new Genotype();
  genotype->load(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  grm->kernel->fillWithConstant(1.);
  
  
  options.useCorrelations = true;
  //options.useCorrelations = false;
  //options.useLogLogisticScale = true;
  //options.firstStepEM = false;
  
  options.phenotypeColumns = std::vector<int>();
  options.phenotypeColumns.push_back(1);
  options.phenotypeColumns.push_back(2);
  options.phenotypeColumns.push_back(3);
  options.initialh2Traits = std::vector<double>();
  options.initialh2Traits.push_back(0.2);
  options.initialh2Traits.push_back(0.2);
  options.initialh2Traits.push_back(0.2);
  options.qCovarsFiles.clear();
  options.covarsFiles.clear();
  
  options.pruneGRM = false;
  options.minimumNumberOverlapingSNPs = 0;
  /*
  if(options.useCorrelations == true)
  {
    options.initialVariancesFile = "test/test3.initialVars";
    //options.initialVariancesFile = "test/test3.v2.initialVars";
  }
  else
  {
    options.initialVariancesFile = "test/test3.covars.initialVars";
  }
  */
  
  std::vector<std::string> keepIndividualIds;
  keepIndividualIds.push_back("i1@i1");
  keepIndividualIds.push_back("i2@i2");
  keepIndividualIds.push_back("i3@i3");
  grm->filterIndividuals(keepIndividualIds);
  
  std::vector<Kernel*> kernels;
  grm->name = options.baseVarianceNames["GRM"];
  kernels.push_back(grm);
  
  std::vector<double> weights;
  weights.assign( kernels.size() , 1./double(kernels.size()) );
  
  std::vector<std::pair<std::string, std::string> > covariateFiles;
  for(int i = 0; i<options.phenotypeColumns.size(); i++)
  {
    covariateFiles.push_back(std::pair<std::string, std::string>("", ""));
  }
  
  REML *reml = new REML();
  
  bool prepared = reml->prepare(bivariateREMLType, kernels, weights, options.phenotypeColumns, options.initialh2Traits, covariateFiles);
  
  reml->V->reinitializeVariances();
  reml->V->variances[1].variance = 0.23;
  reml->V->variances[2].variance = 0.53;
  reml->V->variances[4].variance = 0.07;
  reml->V->variances[8].variance = 0.4;
  reml->V->elements[1].constantFactor = 1.23;
  reml->V->showCovarianceMatrix(true, true);
  
  //just to test;
  //reml->V->elements[1].iVariances.push_back(0);
  //reml->V->elements[1].iVariances.push_back(1);
  
  Matrix * t;
  int i;
  int j;
  //GRMS tests
  i = 0;
  j = 1;
  t = reml->V->computeDerivateCovariance(i, j);
  misc.message << "-------------------------------------" << std::endl << std::endl;
  if( t == NULL )
  {
    misc.message << reml->V->variances[i].name << " " << reml->V->variances[j].name << " is NULL" << std::endl;
  }
  else
  {
    t->showGlobal("Derivate " + reml->V->variances[i].name + " " + reml->V->variances[j].name, false);
  }
  delete t;
  misc.message << "-------------------------------------" << std::endl;
  
  i = 3;
  j = 0;
  t = reml->V->computeDerivateCovariance(i, j);
  misc.message << "-------------------------------------" << std::endl << std::endl;
  if( t == NULL )
  {
    misc.message << reml->V->variances[i].name << " " << reml->V->variances[j].name << " is NULL" << std::endl;
  }
  else
  {
    t->showGlobal("Derivate " + reml->V->variances[i].name + " " + reml->V->variances[j].name, false);
  }
  delete t;
  misc.message << "-------------------------------------" << std::endl;
  
  i = 3;
  j = 3;
  t = reml->V->computeDerivateCovariance(i, j);
  misc.message << "-------------------------------------" << std::endl << std::endl;
  if( t == NULL )
  {
    misc.message << reml->V->variances[i].name << " " << reml->V->variances[j].name << " is NULL" << std::endl;
  }
  else
  {
    t->showGlobal("Derivate " + reml->V->variances[i].name + " " + reml->V->variances[j].name, false);
  }
  delete t;
  misc.message << "-------------------------------------" << std::endl;
  
  i = 1;
  j = 1;
  t = reml->V->computeDerivateCovariance(i, j);
  misc.message << "-------------------------------------" << std::endl << std::endl;
  if( t == NULL )
  {
    misc.message << reml->V->variances[i].name << " " << reml->V->variances[j].name << " is NULL" << std::endl;
  }
  else
  {
    t->showGlobal("Derivate " + reml->V->variances[i].name + " " + reml->V->variances[j].name, false);
  }
  delete t;
  misc.message << "-------------------------------------" << std::endl;
  
  i = 1;
  j = 2;
  t = reml->V->computeDerivateCovariance(i, j);
  misc.message << "-------------------------------------" << std::endl << std::endl;
  if( t == NULL )
  {
    misc.message << reml->V->variances[i].name << " " << reml->V->variances[j].name << " is NULL" << std::endl;
  }
  else
  {
    t->showGlobal("Derivate " + reml->V->variances[i].name + " " + reml->V->variances[j].name, false);
  }
  delete t;
  misc.message << "-------------------------------------" << std::endl;
  
  //E tests
  i = 6;
  j = 7;
  //reml->V->variances[6].variance = 0.2667;
  t = reml->V->computeDerivateCovariance(i, j);
  misc.message << "-------------------------------------" << std::endl << std::endl;
  if( t == NULL )
  {
    misc.message << reml->V->variances[i].name << " " << reml->V->variances[j].name << " is NULL" << std::endl;
  }
  else
  {
    t->showGlobal("Derivate " + reml->V->variances[i].name + " " + reml->V->variances[j].name, false);
  }
  delete t;
  misc.message << "-------------------------------------" << std::endl;
  
  i = 6;
  j = 6;
  t = reml->V->computeDerivateCovariance(i, j);
  misc.message << "-------------------------------------" << std::endl << std::endl;
  if( t == NULL )
  {
    misc.message << reml->V->variances[i].name << " " << reml->V->variances[j].name << " is NULL" << std::endl;
  }
  else
  {
    t->showGlobal("Derivate " + reml->V->variances[i].name + " " + reml->V->variances[j].name, false);
  }
  delete t;
  misc.message << "-------------------------------------" << std::endl;
  
  i = 7;
  j = 7;
  t = reml->V->computeDerivateCovariance(i, j);
  misc.message << "-------------------------------------" << std::endl << std::endl;
  if( t == NULL )
  {
    misc.message << reml->V->variances[i].name << " " << reml->V->variances[j].name << " is NULL" << std::endl;
  }
  else
  {
    t->showGlobal("Derivate " + reml->V->variances[i].name + " " + reml->V->variances[j].name, false);
  }
  delete t;
  misc.message << "-------------------------------------" << std::endl;
  
  //Single Derivate tests
  i = 1;
  t = reml->V->computeDerivateCovariance(i);
  misc.message << "-------------------------------------" << std::endl << std::endl;
  if( t == NULL )
  {
    misc.message << reml->V->variances[i].name << " is NULL" << std::endl;
  }
  else
  {
    t->showGlobal("Derivate " + reml->V->variances[i].name, false);
  }
  //delete t;
  misc.message << "-------------------------------------" << std::endl;
  
  i = 3;
  t = reml->V->computeDerivateCovariance(i);
  misc.message << "-------------------------------------" << std::endl << std::endl;
  if( t == NULL )
  {
    misc.message << reml->V->variances[i].name << " is NULL" << std::endl;
  }
  else
  {
    t->showGlobal("Derivate " + reml->V->variances[i].name, false);
  }
  //delete t;
  misc.message << "-------------------------------------" << std::endl;
  
  
  //reml->computeREML();
  
  delete reml;
}

Matrix * testCreateRandomMatrix(int nr, int nc, int s)
{
  double * tgm;
  
  if(communicator->mpiRoot)
  {
    tgm = new double [nr*nc];
    
    for(int r = 0; r<nr; r++)
    {
      for(int c = 0; c<nc; c++)
      {
        int temp = rand() % s;
        tgm[c*(nr) + r] = double(s)/2. - double(temp);
      }
    }
  }
  Matrix * result = new Matrix(cyclicDistribution, nr, nc);
  result->scatterMatrix(tgm);
  
  return result;
}

void test86()
{
  std::srand(communicator->mpiRank*std::time(0));
  
  Matrix * m1 = testCreateRandomMatrix(20, 20, 100);
  m1->showGlobal("m1", false);
  
  Matrix * m2 = testCreateRandomMatrix(20, 20, 100);
  m2->showGlobal("m2", false);
  
  Matrix * mult = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  mult->multiply(m1, 'N', m2, 'N');
  
  Matrix * mult2 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  mult2->multiply(m2, 'N', m1, 'N');
  
  double b1 = mult->trace();
  double b2 = mult2->trace();
  
  double t1 = m1->traceOfMatrixProduct(m2);
  double t2 = m2->traceOfMatrixProduct(m1);
  misc.message << b1 << " " << b2 << " " << t1 << " " << t2 << std::endl;
  if( b1 == b2 && t1 == t2 && t1 == b1)
  {
    //changed to do nothing
  }
  else
  {
    misc.error("Error: Traces are different.", 0);
  }
  
  delete m1;
  delete m2;
  
  //Repeat when one matrix is symmetric
  
  //m1 = testCreateRandomMatrix(20, 20, 100);
  //delete m1;
  m1 = testCreateRandomMatrix(20, 20, 100);
  m1->symmetric = true;
  //m1->uplo = 'U';
  m1->uplo = 'L';
  m1->showGlobal("m1", true);
  
  m2 = testCreateRandomMatrix(20, 20, 100);
  m2->showGlobal("m2", false);
  
  mult->multiply(m1, 'N', m2, 'N');
  
  mult2->multiply(m2, 'N', m1, 'N');
  
  b1 = mult->trace();
  b2 = mult2->trace();
  
  t1 = m1->traceOfMatrixProduct(m2);
  t2 = m2->traceOfMatrixProduct(m1);
  misc.message << b1 << " " << b2 << " " << t1 << " " << t2 << std::endl;
  if( b1 == b2 && t1 == t2 && t1 == b1)
  {
    //changed to do nothing
  }
  else
  {
    misc.error("Error: Traces are different.", 0);
  }
  
  //repeat when both matrices are symmetric
  
  m1 = testCreateRandomMatrix(20, 20, 100);
  delete m1;
  m1 = testCreateRandomMatrix(20, 20, 100);
  m1->symmetric = true;
  m1->uplo = 'U';
  m1->showGlobal("m1", false);
  
  m2 = testCreateRandomMatrix(20, 20, 100);
  m2->symmetric = true;
  m2->uplo = 'U';
  m2->showGlobal("m2", false);
  
  m1->uplo = 'U'; m2->uplo = 'U';
  //m1->uplo = 'U'; m2->uplo = 'L';
  //m1->uplo = 'L'; m2->uplo = 'U';
  //m1->uplo = 'L'; m2->uplo = 'L';
  
  mult = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  mult->multiply(m1, 'N', m2, 'N');
  
  mult2 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  mult2->multiply(m2, 'N', m1, 'N');
  
  b1 = mult->trace();
  b2 = mult2->trace();
  
  t1 = m1->traceOfMatrixProduct(m2);
  t2 = m2->traceOfMatrixProduct(m1);
  misc.message << b1 << " " << b2 << " " << t1 << " " << t2 << std::endl;
  if( b1 == b2 && t1 == t2 && t1 == b1)
  {
    //changed to do nothing
  }
  else
  {
    misc.error("Error: Traces are different.", 0);
  }
  

}

void test87()
{
  std::srand(communicator->mpiRank*std::time(0));
  
  Matrix * m1 = testCreateRandomMatrix(4, 1, 100);
  m1->showGlobal("m1", false);
  
  Matrix * m2 = testCreateRandomMatrix(4, 3, 100);
  m2->showGlobal("m2", false);
  
  Matrix * m3 = testCreateRandomMatrix(3, 4, 100);
  m3->showGlobal("m3", false);
  
  Matrix * m4 = testCreateRandomMatrix(1, 4, 100);
  m4->showGlobal("m4", false);
  
  misc.message << "-----------------------" << std::endl;
  
  Matrix * mm = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  
  mm->multiply(m1, 'T', m2, 'N');
  mm->showGlobal("mm", false);
  mm->multiply(m2, 'T', m1, 'N');
  mm->showGlobal("mm", false);
  
  misc.message << "-----------------------" << std::endl;
  
  mm->multiply(m4, 'N', m3, 'T');
  mm->showGlobal("mm", false);
  mm->multiply(m3, 'N', m4, 'T');
  mm->showGlobal("mm", false);
  
}

void test88()
{
  int ttnr = 7;
  int ttnc = 5;
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 3;
  
  Matrix md1 = Matrix(diagonalDistribution, ttnc, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  md1.fillWithConstant(0.3);
  md1.showGlobal("md1");
  Matrix md2 = Matrix(diagonalDistribution, ttnc, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  md2.fillWithConstant(0.2);
  md2.showGlobal("md2");
  
  md1.elementWiseMultiplication(&md2);
  md1.showGlobal("md1*md2");
  
  md1.fillWithConstant(0.3);
  md1.elementWiseMultiplication(&md2, 1.5);
  md1.showGlobal("md1*md2*1.5");

  Matrix m1 = Matrix(cyclicDistribution, ttnr, ttnc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  m1.fillWithRandom();
  m1.showGlobal("randomM1");
  m1.scaleBy(0.5);
  m1.showGlobal("randomM1*0.5");
  
  md1.fillWithRandom();
  md1.showGlobal("randomMd1");
  md1.scaleBy(0.5);
  md1.showGlobal("randomMd1*0.5");
}

void test89()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 3;
  
  std::vector<double> v1;
  std::vector<double> v2;
  
  v1.push_back(10);
  v1.push_back(20);
  v1.push_back(30);
  v1.push_back(40);
  v1.push_back(50);
  
  v2.push_back(1);
  v2.push_back(2);
  v2.push_back(3);
  v2.push_back(4);
  v2.push_back(5);
  v2.push_back(6);
  v2.push_back(7);
  
  Matrix * m = new Matrix();
  m->makeDifferenceMatrix(v1, v2);
  m->showGlobal("v1, v2");
  
  m->makeDifferenceMatrix(v2, v1);
  m->showGlobal("v2, v1");
}

void test90()
{
  communicator->nDefaultBlockCols = 2;
  communicator->nDefaultBlockRows = 2;
  
  Kernel test("./test/testExponentialKernel.dat", kernelSquaredExponential);
  test.printGRM();
  
  /*
  >>> i1 = [-1, 1.2, 123.2, -3.2]
  >>> i4 = [3, 2.7, -739.1, 1]
  >>> i9 = [3, -3.3, -190.3, 2.]
  >>> import numpy as np
  >>> i1 = np.array(i1)
  >>> i4 = np.array(i4)
  >>> i9 = np.array(i9)
  >>> (i1-i4)*(i1-i4)
  array([  1.60000000e+01,   2.25000000e+00,   7.43561290e+05,          1.76400000e+01])
  >>> (i1-i4)*(i1-i4).sum()
  array([  -3410.4 ,   -1278.9 ,  735196.98,   -3580.92])
  >>> ((i1-i4)*(i1-i4)).sum()
  743597.18000000017
  >>> ((i1-i9)*(i1-i9)).sum()
  98345.539999999994
  >>> ((i1-i4)*(i1-i4)).sum() -  743597
  0.18000000016763806
  >>> ((i1-i9)*(i1-i9)).sum() -  98345.5
  0.039999999993597157
  >>> ((i9-i1)*(i9-i1)).sum() -  98345.5
  0.039999999993597157
  >>> ((i9-i1)*(i9-i1)).sum()
  98345.539999999994
  >>> ((i9-i4)*(i9-i4)).sum()
  301218.43999999994
  >>> ((i9-i4)*(i9-i4)).sum() - 301218
  0.43999999994412065
  >>> ((i9-i4)*(i9-i4)).sum()
  301218.43999999994
  >>> (-1. - 3.)*(-1. - 3.) + (1.2-2.7)*(1.2-2.7) + (123.2 - (-739.1))*(123.2 - (-739.1)) + (-3.2 - 1.)*(-3.2 - 1.)
  743597.1800000002
  >>> (-1. - 3.)*(-1. - 3.) + (1.2-2.7)*(1.2-2.7) + (123.2 - (-739.1))*(123.2 - (-739.1)) + (-3.2 - 1.)*(-3.2 - 1.) - ((i1-i4)*(i1-i4)).sum()
  0.0
  >>> (-1. - 3.)*(-1. - 3.) + (1.2-2.7)*(1.2-2.7) + (123.2 - (-739.1))*(123.2 - (-739.1)) + (-3.2 - 1.)*(-3.2 - 1.)
  743597.1800000002
  */
}

void test91()
{
  communicator->nDefaultBlockCols = 2;
  communicator->nDefaultBlockRows = 2;
  
  //std::srand(communicator->mpiRank*std::time(0));
  std::srand(45034);
  
  Matrix * m1 = new Matrix(diagonalDistribution, 20, 20);
  m1->fillWithRandom();
  m1->showGlobal("m1");
  
  Matrix * m2 = new Matrix(cyclicDistribution, 20, 20);
  m2->fillWithRandom();
  m2->showGlobal("m2");
  
  Matrix * m3 = new Matrix(cyclicDistribution, 17, 29);
  m3->fillWithRandom();
  m3->showGlobal("m3");
  
  Matrix * m2u = new Matrix(m2);
  Matrix * m2l = new Matrix(m2);
  
  m2u->symmetric = true;
  m2u->uplo = 'U';
  m2l->symmetric = true;
  m2l->uplo = 'L';
  
  /*
  std::cout << "m1 " << m1->elementsAverage() << " " << test91aux(m1) << " " << m1->elementsAverage() - test91aux(m1) << std::endl;
  std::cout << "m2 " << m2->elementsAverage() << " " << test91aux(m2) << " " << m2->elementsAverage() - test91aux(m2) << std::endl;
  std::cout << "m2u " << m2u->elementsAverage() << " " << test91aux(m2u) << " " << m2u->elementsAverage() - test91aux(m2u) << std::endl;
  std::cout << "m2l " << m2l->elementsAverage() << " " << test91aux(m2l) << " " << m2l->elementsAverage() - test91aux(m2l) << std::endl;
  std::cout << "m3 " << m3->elementsAverage() << " " << test91aux(m3) << " " << m3->elementsAverage() - test91aux(m3) << std::endl;
  */
  misc.message << "m1 " << m1->elementsAverage() << " " << test91aux(m1) << " " << m1->elementsAverage() - test91aux(m1) << std::endl;
  misc.message << "m2 " << m2->elementsAverage() << " " << test91aux(m2) << " " << m2->elementsAverage() - test91aux(m2) << std::endl;
  misc.message << "m2u " << m2u->elementsAverage() << " " << test91aux(m2u) << " " << m2u->elementsAverage() - test91aux(m2u) << std::endl;
  misc.message << "m2l " << m2l->elementsAverage() << " " << test91aux(m2l) << " " << m2l->elementsAverage() - test91aux(m2l) << std::endl;
  misc.message << "m3 " << m3->elementsAverage() << " " << test91aux(m3) << " " << m3->elementsAverage() - test91aux(m3) << std::endl;
 
  
  //m1 = np.diag(m1)
  //m2u = np.triu(m2) + np.triu(m2, 1).T
  //m2l = np.tril(m2) + np.tril(m2, -1).T

  // >>> m1.sum()/float(m1.size)
  // 0.108755
  // >>> m2.sum()/float(m2.size)
  // -0.30570000000000019
  // >>> m2u.sum()/float(m2u.size)
  // -0.50612499999999916
  // >>> m2l.sum()/float(m2l.size)
  // -0.10527500000000015
  // >>> m3.sum()/float(m3.size)
  // 0.12518052738336735

}


double test91aux(Matrix * m)
{
  std::vector< std::vector<double> > globalVector;
  m->matrixToStandardVector(globalVector);
  
  double result = 0;
  if(communicator->mpiRoot == true)
  {
    for( int r = 0; r< globalVector.size(); r++)
    {
      for( int c = 0; c< globalVector[r].size(); c++)
      {
        if(m->symmetric == true && m->uplo != 'B')
        {
          if(m->uplo == 'U' && r<c)
          {
            result += 2.*globalVector[r][c];
          }
          else if (m->uplo == 'L' && r > c)
          {
            result += 2.*globalVector[r][c];
          }
          else if(r==c)
          {
            result += globalVector[r][c];
          }
        }
        else
        {
          result += globalVector[r][c];
        }
      }
    }
  }
  result = result/(double(m->nGlobRows)*double(m->nGlobCols));
  communicator->broadcast(&result);
  
  return result;
}

void test92()
{
  communicator->nDefaultBlockCols = 2;
  communicator->nDefaultBlockRows = 2;
  
  Genotype *genotype = new Genotype();
  genotype->load(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  
  std::vector<double> valuesToChange;
  for(int i = 0; i < grm->nIndividuals; i++)
  {
    for(int j = 0; j < grm->nIndividuals; j++)
    {
      valuesToChange.push_back((i+1) + (j+1)*10);
    }
  }
  grm->kernel->scatterMatrix(&(valuesToChange[0]));
  Kernel * grm2 = new Kernel(grm);
  grm->kernel->showGlobal("mod. grm", false);
  
  
  grm->printGRM();
  
  std::vector<std::string> keepIndividualIds;
  keepIndividualIds.push_back("i1@i1");
  keepIndividualIds.push_back("i7@i7");
  keepIndividualIds.push_back("i3@i3");
  keepIndividualIds.push_back("i5@i5");
  keepIndividualIds.push_back("i2@i2");
  
  
  grm->filterIndividuals(keepIndividualIds);
  grm->printGRM();
  
  std::vector<std::string> keepIndividualIdsRows;
  keepIndividualIdsRows.push_back("i1@i1");
  keepIndividualIdsRows.push_back("i7@i7");
  keepIndividualIdsRows.push_back("i3@i3");
  keepIndividualIdsRows.push_back("i5@i5");
  keepIndividualIdsRows.push_back("i2@i2");
  
  std::vector<std::string> keepIndividualIdsCols;
  keepIndividualIdsCols.push_back("i9@i9");
  keepIndividualIdsCols.push_back("i6@i6");
  keepIndividualIdsCols.push_back("i7@i7");
  
  grm2->filterIndividualsAsymmetric(keepIndividualIdsRows, keepIndividualIdsCols);
  grm2->printGRM();
}

void test93()
{
  communicator->nDefaultBlockCols = 2;
  communicator->nDefaultBlockRows = 2;
  
  //std::srand(communicator->mpiRank*std::time(0));
  std::srand(45034);
  
  Matrix * m1 = new Matrix(diagonalDistribution, 3, 3);
  m1->fillWithRandom();
  m1->showGlobal("m1");
  m1->applyExponentialOperator(0.5);
  m1->showGlobal("exp(m1*0.5)");
  
  Matrix * m2 = new Matrix(cyclicDistribution, 3, 3);
  m2->fillWithRandom();
  m2->showGlobal("m2");
  m2->applyExponentialOperator(0.5);
  m2->showGlobal("exp(m2*0.5)");
  
  

}


void test94()
{
  //mpirun --np 2 ./dissect --debug-vars --sqrt-exp-coord-files test/testExponentialKernel.dat
  
  REML reml;
  
  std::vector<std::string> individualBLUPNames;
  std::vector<std::string> reducedModels;
  std::vector<Kernel *> kernels;  

  
  Genotype *genotype = new Genotype();
  genotype->load(options.genotypeFile);
  Kernel * grm = new Kernel(genotype);
  grm->name = options.baseVarianceNames["GRM"];
  kernels.push_back(grm);
  addKernelsUsingOptions(kernels, reducedModels, individualBLUPNames);
  
  options.useCorrelations = true;
  //options.useCorrelations = false;
  //options.useLogLogisticScale = true;
  //options.firstStepEM = false;
  
  options.phenotypeColumns = std::vector<int>();
  options.phenotypeColumns.push_back(1);
  options.phenotypeColumns.push_back(2);
  /*options.phenotypeColumns.push_back(3);
  options.initialh2Traits = std::vector<double>();
  options.initialh2Traits.push_back(0.2);
  options.initialh2Traits.push_back(0.2);
  options.initialh2Traits.push_back(0.2);
  options.qCovarsFiles.clear();
  options.covarsFiles.clear();
  */
  
  options.maxREMLIterations = 2;
  
  options.pruneGRM = false;
  options.minimumNumberOverlapingSNPs = 0;
  
  std::vector<double> weights;
  weights.assign( kernels.size() , 1./double(kernels.size()) );
  
  std::vector<std::pair<std::string, std::string> > covariateFiles;
  for(int i = 0; i<options.phenotypeColumns.size(); i++)
  {
    covariateFiles.push_back(std::pair<std::string, std::string>("", ""));
  }
  
  bool prepared = reml.prepare(singleREMLType, kernels, weights, options.phenotypeColumns, options.initialh2Traits, covariateFiles);
  
  ////////////////////////
  // Test gradient
  double step = 0.001;
  
  reml.type = rawREMLType;
  reml.computeREML();
  
  std::vector<double> baseVariances = reml.V->getVectorVariances();
  std::vector<double> modifiedVariances = reml.V->getVectorVariances();
  
  communicator->broadcast(&baseVariances[0], baseVariances.size());
  communicator->broadcast(&modifiedVariances[0], modifiedVariances.size());
  
  std::vector<double> gradient;
  std::vector<double> testingGradient;
  double baseLogLikelihood = test94aux(reml);
  reml.yPsubVPy_trPsubV->matrixToStandardVector(gradient);
  if(communicator->mpiRoot == false)
  {
    gradient = std::vector<double>(baseVariances.size());
  }
  communicator->broadcast(&gradient[0], gradient.size()); 
  
  misc.message << std::endl;
  for(int i = 0; i < baseVariances.size(); i++)
  {
    modifiedVariances[i] = baseVariances[i] + baseVariances[i]*step;
    //modifiedVariances[i] = baseVariances[i] + step;
    communicator->broadcast(&modifiedVariances[0], modifiedVariances.size()); 

    reml.V->storeVectorVariances(modifiedVariances);
    testingGradient.push_back( (test94aux(reml) - baseLogLikelihood)/(modifiedVariances[i]-baseVariances[i]) );
    
    modifiedVariances = baseVariances;
    reml.V->storeVectorVariances(baseVariances);
    
    for(int p = 0; p<communicator->mpiNumTasks; p++)
    {
      if(communicator->mpiRank == p)
      {
        std::cout << std::setw(20) << reml.V->variances[i].name << std::setw(20) << gradient[i] <<  std::setw(20) << testingGradient[i] << std::setw(20) <<  gradient[i]/testingGradient[i] << std::endl;
      }
      communicator->barrier();
    }
    
  }
  misc.message << std::endl;
  
}

double test94aux(REML & reml)
{
  reml.V->computeCovariance();
  reml.success = reml.V->invertCovariance(&reml.logDetV, reml.singlePrecisionInversion);
  if(reml.success == false)
  {
    misc.error("Error 1", 0);
  }
  
  reml.computePMatrix();
  if(reml.success == false)
  {
    misc.error("Error 2", 0);
  }
  
  reml.computeSubVPyMatrix();
  
  reml.computeAIMatrix();
  
  reml.computeyPsubVPy_trPsubVVector(0.5);
  
  reml.computeLogLikelihood();
  
  double result = reml.logLikelihood;
  communicator->broadcast(&result);
  return result;
}

void test95()
{
  communicator->nDefaultBlockRows = 3;
  communicator->nDefaultBlockCols = 2;
  
  int nRows = 5;
  int nCols = 7;
  
  LabeledMatrix * effects = new LabeledMatrix(nRows, nCols);
  
  std::vector<std::string> temp;
  for(int r = 0; r < nRows; r++)
  {
    temp.push_back("r" + i2s(r));
  }
  effects->setRowLabels(temp);
  temp.clear();
  for(int c = 0; c < nCols; c++)
  {
    temp.push_back("c" + i2s(c));
  }
  effects->setColLabels(temp);
  effects->show();
  
  //////////////////////////////////////
  // Insert rows
  for(int r = 0; r < nRows; r++)
  {
    Matrix * row = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, 1, nCols);
    std::vector<double> g;
    for(int c = 0; c < nCols; c++)
    {
      g.push_back(double(r) + double(c)/10.);
    }
    row->scatterMatrix(&(g[0]));
    
    effects->insertRow(row, r);
    
    delete row;
  }
  effects->show();
  
  for(int r = 0; r < nRows; r++) //Check that replaces properly previous values
  {
    if(r != 1 && r != 3)
    {
      continue;
    }
    
    Matrix * row = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, 1, nCols);
    std::vector<double> g;
    for(int c = 0; c < nCols; c++)
    {
      //g.push_back(0.);
      g.push_back(1.);
      //g.push_back(100.);
    }
    row->scatterMatrix(&(g[0]));
    
    effects->insertRow(row, r);
    
    delete row;
  }
  effects->show();
  
  //////////////////////////////////////
  // Insert cols
  for(int c = 0; c < nCols; c++) //Check that replaces properly previous values
  {
    if(c != 1 && c != 4)
    {
      continue;
    }
    Matrix * col = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nRows, 1);
    std::vector<double> g;
    for(int r = 0; r < nRows; r++)
    {
      g.push_back(0.);
    }
    col->scatterMatrix(&(g[0]));
    
    effects->insertCol(col, c);
    
    delete col;
  }
  effects->show();
  
  for(int c = 0; c < nCols; c++)
  {
    Matrix * col = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nRows, 1);
    std::vector<double> g;
    for(int r = 0; r < nRows; r++)
    {
      g.push_back(double(c) + double(r)/10.);
    }
    col->scatterMatrix(&(g[0]));
    
    effects->insertCol(col, c);
    
    delete col;
  }
  effects->show();
  
  //////////////////////////////////////
  // Save/load
  effects->save("./tmp/effects");
  LabeledMatrix * loadEffects = new LabeledMatrix("./tmp/effects");
  loadEffects->show();
  
  compareMatrices(effects->getMatrix(), loadEffects->getMatrix());
  if(effects->getRowLabels() != loadEffects->getRowLabels() || effects->getColLabels() != loadEffects->getColLabels() )
  {
    misc.error("Error: Different labels.", 0);
  }
  
  //////////////////////////////////////
  // Filter
  std::vector<std::string> keepRows;
  keepRows.push_back("r3");
  keepRows.push_back("r0");
  keepRows.push_back("r2");
  std::vector<std::string> keepCols;
  keepCols.push_back("c5");
  keepCols.push_back("c2");
  keepCols.push_back("c4");
  keepCols.push_back("c1");
  
  effects->filterRowsAndCols(keepRows, keepCols);
  effects->show();
  
  
  delete effects;
}

void test96()
{
  //python:
  //from sklearn import preprocessing as pp
  //np.abs((pp.scale(m2, axis = 0) - m2s)).max()
  
  communicator->nDefaultBlockRows = 3;
  communicator->nDefaultBlockCols = 2;
  
  int nRows = 27;
  int nColumns = 13;
  
  Matrix * m = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nRows, nColumns);
  Matrix * mt = new Matrix();
  m->fillWithRandom();
  mt->transpose(m);
  
  m->showGlobal("m");
  //mt->showGlobal("mt");
  m->standardizeMatrix(row);
  mt->standardizeMatrix(column);
  m->showGlobal("m row standardized");
  //mt->showGlobal("m row standardized");
  
  Matrix * mtt = new Matrix();
  mtt->transpose(mt);
  compareMatrices(m, mtt, 1e-8);
  
  Matrix * m2 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nRows, nColumns);
  m2->fillWithRandom();
  m2->showGlobal("m2");
  m2->standardizeMatrix(column);
  m2->showGlobal("m2 col standardized");

}

void test97()
{
  //python:
  //import numpy as np
  //np.abs((m - m.mean(axis = 1)[:,np.newaxis]) - mc).max()
  
  communicator->nDefaultBlockRows = 3;
  communicator->nDefaultBlockCols = 2;
  
  int nRows = 17;
  int nColumns = 13;
  /*int nRows = 3;
  int nColumns = 3;*/
  
  Matrix * m = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nRows, nColumns);
  Matrix * mt = new Matrix();
  m->fillWithRandom();
  mt->transpose(m);
  
  m->showGlobal("m");
  //mt->showGlobal("mt");
  m->centerMatrixRowsColumns(row);
  mt->centerMatrixRowsColumns(column);
  m->showGlobal("m row centered");
  //mt->showGlobal("m row standardized");
  
  Matrix * mtt = new Matrix();
  mtt->transpose(mt);
  compareMatrices(m, mtt, 1e-8);
  
  Matrix * m2 = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nRows, nColumns);
  m2->fillWithRandom();
  m2->showGlobal("m2");
  m2->centerMatrixRowsColumns(column);
  m2->showGlobal("m2 col centered");

  //Test means
  Matrix * mROnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nRows, 1);
  mROnes->fillWithConstant(1.);
  Matrix * mCOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nColumns, 1);
  mCOnes->fillWithConstant(1.);
  
  Matrix * mMeans = new Matrix();
  Matrix * mMeans2 = new Matrix();
  
  mMeans->multiply(m, 'N', mCOnes, 'N', 1./double(nColumns));
  mMeans->showGlobal("Means row centered");
  mMeans2->multiply(m2, 'T', mROnes, 'N', 1./double(nRows));
  mMeans2->showGlobal("Means column centered");
  
  
  
}

void test98()
{
  /*python
  X = np.random.normal(0, 1, (100, 10))
  betas = np.random.normal(5, 20, (10, 1))
  blups = np.dot(X, betas)
  p = blups + np.random.norm(0, np.std(blups)*0.2, (100, 1) )
  y = p.copy()
  y[y<50] = 0
  y[y>0] = 1
  np.savetxt("./glm/y.dat", y)
  np.savetxt("./glm/X.dat", X)
  np.savetxt("./glm/blups.dat", blups)
  np.savetxt("./glm/p.dat", p)
  np.savetxt("./glm/betas.dat", betas)
  
  randomEffects = np.random.normal(0, 0.1, (100, 5))
  np.savetxt("./glm/randomEffects.dat", randomEffects)
  
  
  import statsmodels.api as sm
  result = sm.GLM(y, X, family = sm.families.Binomial()).fit()
  result.summary()
  */
  
  Matrix * y = new Matrix();
  y->debugRead("./test/glm/y.dat", 100, 1);
  
  Matrix * X = new Matrix();
  X->debugRead("./test/glm/X.dat", 100, 10);
  
  Matrix * randomEffects = new Matrix();
  randomEffects->debugRead("./test/glm/randomEffects.dat", 100, 10);
  
  GLM glm(&y, &X);
  //GLM glm(&y, &X, NULL, &randomEffects);
  
  glm.fit();
  
  if(glm.success == true)
  {
    misc.message << "GLM fit after " << glm.nIterations << " iterations." << std::endl;
  }
  else
  {
    misc.message << "GLM failed to fit." << std::endl;
  }
  
  glm.betas->showGlobal("betas");
}

void test99()
{
  /*
   from scipy import stats
   100*2*stats.norm.sf(4, 0, 1)
   */
  
  long seed = options.randomSeed;
  int nNumbers = 100000000;
  double * numbers = new double[nNumbers];
  
  misc.setGetElapsedTime("Random");
  for(int i = 0; i< nNumbers; i++)
  {
    numbers[i] = box_muller(0., 1., &seed);
  }
  misc.message << nNumbers << " random numbers generated after " << misc.setGetElapsedTime("Random", true) << "." << std::endl;
  
  long nByStd[10];
  for(int j = 1; j< 10; j++)
  {
    nByStd[j] = 0;
  }
  for(int i = 0; i< nNumbers; i++)
  {
    for(int j = 1; j< 10; j++)
    {
      if(fabs(numbers[i]) > double(j))
      {
        nByStd[j] += 1;
      }
    }
  }
  for(int j = 1; j< 10; j++)
  {
    std::cout << j << ": " << nByStd[j] << " " << std::setprecision(3) << 100*double(nByStd[j])/double(nNumbers) << "%" << std::endl;
  }
  
  delete [] numbers;
}

void test100()
{
  long seed = options.randomSeed;
  int nNumbers = 1000000;
  double * numbers = new double[nNumbers];
  
  misc.setGetElapsedTime("Random");
  Message message("./test/uniformRandom.dat");
  for(int i = 0; i< nNumbers; i++)
  {
    numbers[i] = unif_rand_dbl(&seed);
    message << numbers[i] << std::endl;
    if(numbers[i] < 0. || numbers[i] > 1.)
    {
      misc.error("Error: Outside bonds.", 0);
    }
  }
  misc.message << nNumbers << " random numbers generated after " << misc.setGetElapsedTime("Random", true) << "." << std::endl;
  
  int nBins = 50;
  long * bins = new long [nBins];
  for(int j = 0; j< nBins; j++)
  {
    bins[j] = 0;
  }
  for(int i = 0; i< nNumbers; i++)
  {
    for(int j = 0; j< nBins; j++)
    {
      if((numbers[i] > double(j)/double(nBins)) && (numbers[i] < double(j+1)/double(nBins)))
      {
        bins[j] += 1;
      }
    }
  }
  for(int j = 0; j< nBins; j++)
  {
    std::cout << j << ": " << bins[j] << " " << std::setprecision(3) << double(nBins)*double(bins[j])/double(nNumbers) << std::endl;
  }
  
  delete [] numbers;
  delete [] bins;
}

void test101()
{
  //long seed = options.randomSeed;
  long seed = 27;
  int nIndividuals = 100.;
  int nSNPs = 1000;
  
  double * mGlob;
  double * yGlob;

  Matrix * gen = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nIndividuals, nSNPs);
  gen->fillWithRandom(-10., 10., &seed);
  gen->standardizeMatrix(column);
  Matrix * m = new Matrix();
  m->multiply(gen, 'N', gen, 'T', 1./double(nSNPs));
  m->symmetrizeTriangularMatrix();
  //m->debugWrite("./test/m.dat");
  
  Matrix * diagonalAdd = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nIndividuals, nIndividuals);
  diagonalAdd->fillDiagonal(1., 0.);
  m->add(diagonalAdd);
  
  Matrix * y = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nIndividuals, 1);
  
  if (communicator->mpiRoot) {
    mGlob = new double[m->nGlobRows*m->nGlobCols];
    yGlob = new double[y->nGlobRows];
    for (int r = 0; r < m->nGlobRows; ++r) {
      if(unif_rand_dbl(&seed) > 0.7)
      {
        *(yGlob + r) = 1.;
      }
      else
      {
        *(yGlob + r) = 0.;
      }
      for (int c = 0; c < m->nGlobCols; ++c) {
        if(r == c)
        {
          *(mGlob + m->nGlobRows*c + r) = 1.;
        }
        else
        {
          *(mGlob + m->nGlobRows*c + r) = unif_rand_dbl(&seed);
        }
      }
    }
  }
  //m->scatterMatrix(mGlob);
  m->symmetric = true;
  m->uplo = 'U';
  y->scatterMatrix(yGlob);
  if (communicator->mpiRoot) {
    delete [] mGlob;
    delete [] yGlob;
  }
  
  Matrix * mtest = new Matrix(m);
  
  REML * reml = new REML;
  reml->V = new CovarianceMatrix(nIndividuals);
  reml->V->m = m;
  reml->V->mStatus = inverseCovarianceMatrix;
  m->debugWrite("./test/m.dat");
  bool success = reml->V->m->symmetricInvert(NULL);
  //bool success = reml->V->m->invert(NULL);
  reml->V->m->debugWrite("./test/V.dat");
  if(success == false)
  {
    misc.error("Error: matrix not inverted.", 0);
  }
  reml->V->m->symmetric = true;
  reml->V->m->uplo = 'U';
  reml->V->m->symmetrizeTriangularMatrix();
  
  Matrix * itest = new Matrix();
  itest->multiply(mtest, 'N', reml->V->m, 'N');
  //itest->showGlobal("I");
  //m->showGlobal("m");
  
  
  Matrix * X = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nIndividuals, 1);  
  X->fillWithConstant(1.);
  //Matrix * betas = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, 1, 1);  
  //betas->fillWithConstant(1.);
  
  reml->y = y;
  reml->X = X;
  GLMM * glmm = new GLMM(reml);
  //glmm->betas = betas;
  
  int nSamples = 1000;
  Matrix * Us = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nIndividuals, nSamples);  
  Us->fillWithConstant(0.);
  Matrix * previousU = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, nIndividuals, 1);  
  Matrix * currentU;
  previousU->fillWithRandom(-10., 10., &seed);
  for(int i = 0; i < nSamples; i++)
  {
    misc.message << i << " ";
    currentU = glmm->MHSampling(previousU);
    delete previousU;
    previousU = currentU;
    Us->add(currentU, 1., 1., subMatrix(0, i, nIndividuals, 1), subMatrix(currentU) );
  }
  delete previousU;
  
  Us->debugWrite("./test/Us.dat");
  //reml->V->m->debugWrite("./test/V.dat");
  
}

void test102()
{
  LabeledMatrix * lm11 = test102CreateTestLabeledMatrix(5, 3, "1", "1", 1.);
  LabeledMatrix * lm12 = test102CreateTestLabeledMatrix(3, 3, "2", "1", 10.);
  lm11->show();
  lm12->show();
  
  lm11->appendVertically(lm12);
  lm11->show();
  
  LabeledMatrix * lm21 = test102CreateTestLabeledMatrix(3, 3, "1", "1", 1.);
  LabeledMatrix * lm22 = test102CreateTestLabeledMatrix(3, 5, "1", "2", 10.);
  lm21->show();
  lm22->show();
  
  lm21->appendHorizontally(lm22);
  lm21->show();
  
  
  
  //effects->show();
}

LabeledMatrix * test102CreateTestLabeledMatrix(int nr, int nc, std::string rprefix, std::string cprefix, double scale)
{
  LabeledMatrix * lm = new LabeledMatrix(nr, nc);
  
  Matrix m = Matrix(cyclicDistribution, nr, nc, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  double * mGlob;
  if (communicator->mpiRoot) {
    mGlob = new double[m.nGlobRows*m.nGlobCols];
    for (int r = 0; r < m.nGlobRows; ++r) {
      for (int c = 0; c < m.nGlobCols; ++c) {
        *(mGlob + m.nGlobRows*c + r) = ((c + 1) + 10*(r + 1))*scale;
      }
    }
  }
  lm->getMatrix()->scatterMatrix(mGlob);
  if (communicator->mpiRoot) {
    delete [] mGlob;
  }
  
  std::vector<std::string> temp;
  for(int r = 0; r < nr; r++)
  {
    temp.push_back(rprefix + "r" + i2s(r));
  }
  lm->setRowLabels(temp);
  temp.clear();
  for(int c = 0; c < nc; c++)
  {
    temp.push_back(cprefix + "c" + i2s(c));
  }
  lm->setColLabels(temp);
  temp.clear();
  
  
  return lm;
}

void test103()
{
  options.maximumProportionOfCategories = 1.;
  
  //options.qCovarsFile = "test/test3.qcov";
  options.qCovarsFile = "test/test3.prefilter.qcov";
  options.covarsFile = "test/test3.cov";
  
  //Filtering results
  Covariate *covariate = new Covariate(options.covarsFile, options.qCovarsFile, std::vector<std::string>());
  misc.message << "\n\n---------------------------\n\n" << "BASE" << "\n\n---------------------------\n\n" << std::endl;
  covariate->printCovariate(14);
  
  std::vector<std::string> test;
  if(communicator->mpiRoot)
  {
    test.push_back("i1@i1");
    test.push_back("i3@i3");
    test.push_back("i4@i4");
    test.push_back("i5@i5");
    test.push_back("i6@i6");
    test.push_back("i7@i7");
    test.push_back("i8@i8");
    
    
    //test.push_back("i2@i2");
    //test.push_back("i9@i9");
  }
  
  misc.message << "\n\n---------------------------\n\n" << "FILTERED" << "\n\n---------------------------\n\n" << std::endl;
  covariate->filterIndividuals(test);
  covariate->printCovariate(14);
  //delete covariate;
  
  //Prefiltering results
  Covariate *covariate2 = new Covariate(options.covarsFile, options.qCovarsFile, std::vector<std::string>(), true, std::set<std::string>(test.begin(), test.end()));
  misc.message << "\n\n---------------------------\n\n" << "PREFILTERED" << "\n\n---------------------------\n\n" << std::endl;
  covariate2->printCovariate(14);
  
  //Comparing matrices after removing and resorting columns.
  std::vector<int> rowsToFilter;
  std::vector<int> colsToFilter;
  colsToFilter.push_back(1);
  colsToFilter.push_back(6);
  Matrix * m = new Matrix();
  covariate->covariates->filterRowsAndColumns(m, rowsToFilter, colsToFilter, false);
  
  std::map<int, int> rowsOriginDestination;
  std::map<int, int> colsOriginDestination;
  for(int i = 0; i < covariate->nIndividuals; i++)
  {
    rowsOriginDestination[ i ] = i;
  }
  for(int i = 0; i < m->nCols; i++)
  {
    if( i == 3 )
    {
      colsOriginDestination[ 3 ] = 4;
    }
    else if( i == 4)
    {
      colsOriginDestination[ 4 ] = 3;
    }
    else
    {
      colsOriginDestination[ i ] = i;
    }
  }
  Matrix * m2 = new Matrix();
  m->generalResorting(m2, rowsOriginDestination, colsOriginDestination, false);
  communicator->barrier();
  
  
  
  compareGlobalMatrices(m2, covariate2->covariates);
  compareMatrices(m2, covariate2->covariates);
  //m2->showGlobal();
  //covariate2->covariates->showGlobal();
  
  //Deleting
  delete covariate;
  delete covariate2;
  delete m;
  delete m2;
}

void test104()
{
  Kernel * kdissect = new Kernel("./test/test3.dissect");
  Kernel * kgcta = new Kernel("./test/test3", kernelGCTAGRM);

  if(kdissect->individualIds != kgcta->individualIds)
  {
    misc.error("Error: Different ids.", 0);
  }
  if(kdissect->individualIdsIdx != kgcta->individualIdsIdx)
  {
    misc.error("Error: Different idIdxs.", 0);
  }
  if(kdissect->individuals != kgcta->individuals)
  {
    misc.error("Error: Different individuals.", 0);
  }
  
  kdissect->kernel->symmetrizeTriangularMatrix();
  compareMatrices(kdissect->kernel, kgcta->kernel, 1e-5);
  
  kdissect->N->symmetrizeTriangularMatrix();
  compareMatrices(kdissect->N, kgcta->N);
  
  //kdissect->printGRM(true);
  //kgcta->printGRM(true);
}

void test105()
{
  //also check with gdb the communicator->broadcast with c++ vectors.
  communicator->nDefaultBlockRows = 3;
  communicator->nDefaultBlockCols = 2;
  
  Kernel * multiCovKernel = new Kernel("./test/test.multirandom.dat", kernelFromMultiDiscreteCovariates, 0);
  multiCovKernel->printGRM(true, 15);
}

void test106()
{
  Phenotype phenotype = Phenotype(MATRIX_DEFAULT_DISTRIBUTION, "./test/testColumns.phenos", 2);
  
  //std::cout << phenotype.nPhenotypesInFile << std::endl;
  
  phenotype.printPhenotype();
}

void test107()
{
  std::vector<int> test;
  if(communicator->mpiRoot)
  {
    test.push_back(2);
    test.push_back(4);
    test.push_back(6);
    test.push_back(1);
    test.push_back(3);
    test.push_back(5);
    test.push_back(7);
  }
  if(communicator->mpiRank == 1)
  {
    test.push_back(-1);
    test.push_back(-2);
    test.push_back(-3);
  }
  if(communicator->mpiRank == 2)
  {
    test.clear();
  }
  if(communicator->mpiRank == 3)
  {
    test.push_back(-10);
    test.push_back(-20);
    test.push_back(-30);
    test.push_back(-12);
    test.push_back(-22);
    test.push_back(-32);
    test.push_back(-13);
    test.push_back(-23);
    test.push_back(-33);
    test.push_back(-14);
    test.push_back(-24);
    test.push_back(-34);
    test.push_back(-15);
    test.push_back(-25);
    test.push_back(-35);
    test.push_back(-16);
    test.push_back(-26);
    test.push_back(-36);
  }
  
  for(int i = 0; i<communicator->mpiNumTasks; i++)
  {
    if(i == communicator->mpiRank)
    {
      std::cout << i << "\n--------" << std::endl;
      for(int iv = 0; iv<test.size(); iv++)
      {
        std::cout << test[iv] << " ";
      }
      std::cout << std::endl;
    }
    std::cout.flush();
    communicator->barrier();
  }
  
  communicator->broadcast(test);
  
  misc.message << "\n**************************\n" << std::endl;
  std::cout.flush();
  
  for(int i = 0; i<communicator->mpiNumTasks; i++)
  {
    if(i == communicator->mpiRank)
    {
      std::cout << i << "\n--------" << std::endl;
      for(int iv = 0; iv<test.size(); iv++)
      {
        std::cout << test[iv] << " ";
      }
      std::cout << std::endl;
    }
    std::cout.flush();
    communicator->barrier();
  }
}

void test108()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 2;
  
  Genotype gbase("./test/test3", false);
  Genotype gfull("./test/test3");
  
  gfull.printGenotype(false);
  
  std::set<std::string> snps;
  snps.insert("rs7");
  snps.insert("rs3");
  snps.insert("rs6");
  std::set<std::string> inds;
  inds.insert("i9@i9");
  inds.insert("i2@i2");
  inds.insert("i7@i7");
  inds.insert("i1@i1");
  Genotype gpart("./test/test3", &gbase, snps, inds);
  gpart.printGenotype(false);

}

void test108bis()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 2;
  
  Genotype gfull("./test/test3");
  
  //gfull.normalizeGenotypes();
  
  /*
  std::set<std::string> snps;//(gfull.individualIds.begin(), gfull.individualIds.end());
  snps.insert("rs1");
  snps.insert("rs2");
  snps.insert("rs3");
  snps.insert("rs4");
  snps.insert("rs5");
  snps.insert("rs6");
  snps.insert("rs7");
  std::set<std::string> inds;//(gfull.SNPIds.begin(), gfull.SNPIds.end());
  inds.insert("i1@i1");
  inds.insert("i2@i2");
  inds.insert("i3@i3");
  inds.insert("i4@i4");
  inds.insert("i5@i5");
  inds.insert("i6@i6");
  inds.insert("i7@i7");
  inds.insert("i8@i8");
  inds.insert("i9@i9");
  */
  std::set<std::string> snps(gfull.SNPIds.begin(), gfull.SNPIds.end());
  std::set<std::string> inds(gfull.individualIds.begin(), gfull.individualIds.end());
  
  Genotype gfullr("./test/test3", snps, inds);
  //gfullr.normalizeGenotypes();
  
  compareMatrices(gfullr.genotypes, gfull.genotypes);
}

void test109()
{
  communicator->nDefaultBlockRows = 2;
  communicator->nDefaultBlockCols = 2;
  
  Genotype gfull("./test/test3");
  Genotype gfull2("./test/test3", false);
  gfull.groupSNPsByFileOrderedWindows(5);
  gfull2.groupSNPsByFileOrderedWindows(2);
  gfull.printGenotype(true);
  
  if(communicator->mpiRoot)
  {
    misc.message << "********" << std::endl;
    misc.message << "Genotype grouped by: " << gfull2.groupedBy << std::endl;
    misc.message << gfull2.groupedSNPs.size() << " groups." << std::endl;
    misc.message << "********" << std::endl;
    
    std::map<std::string, std::set<std::string> >::iterator it;
    for(it = gfull2.groupedSNPs.begin(); it != gfull2.groupedSNPs.end(); ++it)
    {
      misc.message << it->first << std::endl;
      std::set<std::string>::iterator it2;
      for(it2 = it->second.begin(); it2 != it->second.end(); ++it2)
      {
        misc.message << *it2 << std::endl;
      }
    }
  }
 
}

void test110()
{
  std::set<std::string> keepSNPs;
  std::set<std::string> keepIndividualIds;
  
  Genotype genotype("../bgentest/data/example.8bits.bgen", keepSNPs, keepIndividualIds, GenotypeAttributes::probabilities );
  
  //genotype.printGenotype(false, true);
  //genotype.printGenotype(false, false);
  
  genotype.printGenotypeT(true);
  
  Matrix * m = new Matrix();
  m->multiply(genotype.genotypes, 'N', genotype.genotypes, 'T');
  m->showGlobal("");
  
  Matrix * mOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, genotype.nIndividuals, 1);
  mOnes->fillWithConstant(1.);
  Matrix * nmissings = new Matrix();
  nmissings->multiply(genotype.missings, 'N', mOnes, 'N');
  nmissings->showGlobal("");  
  
  std::vector<double> diagonal = m->diagonal();
  if(communicator->mpiRoot == true)
  {
    misc.message << "\n\nN individuals: " << genotype.nIndividuals << std::endl;
    for(int i = 0; i<diagonal.size(); i++)
    {
      misc.message << diagonal[i] << std::endl;
    }
  }
  
  delete m;
  delete mOnes;
  delete nmissings;
}

void test111()
{
  //Checked errors by changing in the file: repeated header elements, not enough or valid elements in the values. Starting header values check.
  std::vector<std::string> lCols;
  lCols.push_back("FID");
  lCols.push_back("IID");
  
  LabeledMatrix *lm = new LabeledMatrix("./test/mphenos.phenos", lCols);
  lm->show();
  
  //////////////////////////////////////
  // Filter
  std::vector<std::string> keepRows;
  keepRows.push_back("i4@i4");
  keepRows.push_back("i1@i1");
  
  std::vector<std::string> keepCols;
  keepCols.push_back("P4");
  keepCols.push_back("P2");
  keepCols.push_back("P5");
  
  lm->filterRowsAndCols(keepRows, keepCols);
  lm->show();
}

void test112()
{
  options.allowBGENLayout1 = true;
  
  std::set<std::string> keepSNPs;
  std::set<std::string> keepIndividualIds;
  
  Genotype genotype1("./test/bgen/t", keepSNPs, keepIndividualIds, GenotypeAttributes::probabilities );
  Genotype genotype2("./test/bgen/t" );
  genotype2.normalizeGenotypes();
 
  genotype1.genotypes->showGlobal("");
  genotype2.genotypes->showGlobal("");
  
  for(int i = 0; i < genotype1.SNPs.size(); i++)
  {
    std::cout << genotype1.SNPs[i].name << " " << genotype1.SNPs[i].standardDev << " " << genotype2.SNPs[i].standardDev << std::endl;
  }
  
  //genotype1.printGenotype(false);
  //genotype2.printGenotype(false);
  
  Matrix * m1 = new Matrix();
  m1->multiply(genotype1.genotypes, 'N', genotype1.genotypes, 'T');
  m1->showGlobal("c1");
  
  Matrix * m2 = new Matrix();
  m2->multiply(genotype2.genotypes, 'N', genotype2.genotypes, 'T');
  m2->showGlobal("c2");
}

void test113()
{
  //Tested with
  //mpirun --np 6 ./dissect
  //mpirun --np 9 ./dissect
  //mpirun --np 4 ./dissect
  //mpirun --np 8 ./dissect
  //mpirun --np 16 ./dissect
  
  options.allowBGENLayout1 = true;
  
  std::vector<int> subsets;
  subsets.push_back(100);
  subsets.push_back(101);
  subsets.push_back(102);
  subsets.push_back(103);
  subsets.push_back(104);
  subsets.push_back(250);
  subsets.push_back(500);
  subsets.push_back(750);
  subsets.push_back(900);
  
  for(int nr = 2; nr <= 5; nr++)
  {
    for(int nc = 2; nc <= 5; nc++)
    {
      communicator->nDefaultBlockRows = nr;
      communicator->nDefaultBlockCols = nc;
      
      misc.message << "\n-----------------------\n" << communicator->nDefaultBlockRows << " " << communicator->nDefaultBlockCols << "\n-----------------------\n" << std::endl;
      
      for(int i = 0; i<subsets.size(); i++)
      {
        std::set<std::string> keepIndividualIdsEmpty;
        std::set<std::string> keepSNPsSetEmpty;
        
        std::string subsetFile = "../testScripts/inputs/randomData/subsets/randomdata.randomsnps." + i2s(subsets[i]);

        std::vector<std::string> keepSNPs;
#if defined(BGEN) && defined(ZLIB)
        Genotype genotypeSNPs;
        keepSNPs = genotypeSNPs.getBGENFileSNPs(subsetFile);
#else
        misc.error("Error", 0);
#endif
        std::set<std::string> keepSNPsSet(keepSNPs.begin(), keepSNPs.end());
        
        Genotype genotypeFull("../testScripts/inputs/randomData/subsets/randomdata", keepSNPsSet, keepIndividualIdsEmpty, GenotypeAttributes::probabilities );
        Genotype genotypeSubset(subsetFile, keepSNPsSetEmpty, keepIndividualIdsEmpty, GenotypeAttributes::probabilities );
        
        //Genotype genotypeSubset("../testScripts/inputs/randomData/subsets/randomdata.randomsnps." + i2s(subsets[i + 1]), keepSNPsSetEmpty, keepIndividualIdsEmpty, GenotypeAttributes::probabilities ); //null test.
        //genotypeFull.genotypes->m[0] = 0.; //null test.
        
        int r = compareMatrices(genotypeFull.genotypes, genotypeSubset.genotypes);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        
        if(genotypeFull.SNPs != genotypeSubset.SNPs)
        {
          misc.error("Discordant SNPs.", 0);
        }
        if(genotypeFull.individuals != genotypeSubset.individuals)
        {
          misc.error("Discordant Individuals.", 0);
        }
        
        if(genotypeFull.SNPIds != genotypeSubset.SNPIds)
        {
          misc.error("Discordant SNP Ids.", 0);
        }
        if(genotypeFull.individualIds != genotypeSubset.individualIds)
        {
          misc.error("Discordant Individual Ids.", 0);
        }
        
        if(genotypeFull.SNPIdsIdx != genotypeSubset.SNPIdsIdx)
        {
          misc.error("Discordant SNP Idxs.", 0);
        }
        if(genotypeFull.individualIdsIdx != genotypeSubset.individualIdsIdx)
        {
          misc.error("Discordant Individual Idxs.", 0);
        }
        
        if( i == 0 && false)
        {
          genotypeFull.printGenotype();
        }
        misc.message << genotypeFull.genotypes->nBlockRows << " " << genotypeFull.genotypes->nBlockCols << std::endl;
      }
    }
  }
  
  //testar amd keepSNPs buit
  
}

void test114()
{
  //Tested with
  //To test: mpirun --np 6 ./dissect, mpirun --np 8 ./dissect, mpirun --np 16 ./dissect
  //mpirun --np 9 ./dissect
  //mpirun --np 4 ./dissect
  //Checked if detects a manual change in SNPs allele frequencies raises error.
  
  std::vector<int> SNPsubsets;
  SNPsubsets.push_back(100);
  SNPsubsets.push_back(101);
  SNPsubsets.push_back(102);
  SNPsubsets.push_back(103);
  SNPsubsets.push_back(104);
  SNPsubsets.push_back(250);
  SNPsubsets.push_back(500);
  SNPsubsets.push_back(750);
  SNPsubsets.push_back(900);
  
  std::vector<int> Indsubsets;
  Indsubsets.push_back(100);
  Indsubsets.push_back(500);
  Indsubsets.push_back(501);
  Indsubsets.push_back(502);
  Indsubsets.push_back(503);
  Indsubsets.push_back(504);
  Indsubsets.push_back(700);
  Indsubsets.push_back(900);
  
  for(int nr = 2; nr <= 5; nr++)
  {
    for(int nc = 2; nc <= 5; nc++)
    {
      communicator->nDefaultBlockRows = nr;
      communicator->nDefaultBlockCols = nc;
      
      misc.message << "\n-----------------------\n" << communicator->nDefaultBlockRows << " " << communicator->nDefaultBlockCols << "\n-----------------------\n" << std::endl;
      
      for(int i = 0; i<SNPsubsets.size(); i++)
      {
      for(int j = 0; j<Indsubsets.size(); j++)
      {
        std::set<std::string> keepIndividualIdsEmpty;
        std::set<std::string> keepSNPsSetEmpty;
        
        std::string subsetFile = "../testScripts/inputs/randomData/subsets2/randomdata.randomsnps.s" + i2s(SNPsubsets[i]) + ".i" + i2s(Indsubsets[j]);
        
        /////////////
        //Load subset SNPs and Inds
        std::vector< std::vector<std::string> > temptable;
        
        getTableFromFile(subsetFile + ".fam", temptable, 6);
        std::set<std::string> keepInds;
        for(int tidx = 0; tidx < temptable.size(); tidx++)
        {
          keepInds.insert(temptable[tidx][0] + "@" + temptable[tidx][1]);
        }
        getTableFromFile(subsetFile + ".bim", temptable, 6);
        std::set<std::string> keepSNPs;
        for(int tidx = 0; tidx < temptable.size(); tidx++)
        {
          keepSNPs.insert(temptable[tidx][1]);
        }
        
        std::set<std::string> keepSNPsEmpty;
        std::set<std::string> keepIndsEmpty;
        Genotype genotypeFull("../testScripts/inputs/randomData/randomdata", keepSNPs, keepInds);
        Genotype genotypeSubset(subsetFile, keepSNPsEmpty, keepIndsEmpty);
        Genotype genotypeSubset2(subsetFile);
        
        //Genotype genotypeSubset("../testScripts/inputs/randomData/subsets/randomdata.randomsnps." + i2s(subsets[i + 1]), keepSNPsSetEmpty, keepIndividualIdsEmpty, GenotypeAttributes::probabilities ); //null test.
        //genotypeFull.genotypes->m[0] = 0.; //null test.
        
        //////////////////////
        // Check matrices
        int r;
        r = compareMatrices(genotypeFull.genotypes, genotypeSubset.genotypes);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        r = compareMatrices(genotypeFull.missings, genotypeSubset.missings);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        
        r = compareMatrices(genotypeFull.genotypes, genotypeSubset2.genotypes);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        r = compareMatrices(genotypeFull.missings, genotypeSubset2.missings);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        
        genotypeFull.normalizeGenotypes();
        genotypeSubset.normalizeGenotypes();
        genotypeSubset2.normalizeGenotypes();
        
        r = compareMatrices(genotypeFull.genotypes, genotypeSubset.genotypes);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        
        r = compareMatrices(genotypeFull.genotypes, genotypeSubset2.genotypes);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        
        //////////////////////
        // Add check it captures some test.
        //genotypeFull.SNPs[0].frequencies[1] = 23;
        //genotypeFull.SNPs[0].p1 = 0.1;
        
        //////////////////////
        // Check other data
        if(genotypeFull.SNPs != genotypeSubset.SNPs || genotypeFull.SNPs != genotypeSubset2.SNPs)
        {
          misc.error("Discordant SNPs.", 0);
        }
        if(genotypeFull.individuals != genotypeSubset.individuals || genotypeFull.individuals != genotypeSubset2.individuals)
        {
          misc.error("Discordant Individuals.", 0);
        }
        
        if(genotypeFull.SNPIds != genotypeSubset.SNPIds || genotypeFull.SNPIds != genotypeSubset2.SNPIds)
        {
          misc.error("Discordant SNP Ids.", 0);
        }
        if(genotypeFull.individualIds != genotypeSubset.individualIds || genotypeFull.individualIds != genotypeSubset2.individualIds)
        {
          misc.error("Discordant Individual Ids.", 0);
        }
        
        if(genotypeFull.SNPIdsIdx != genotypeSubset.SNPIdsIdx || genotypeFull.SNPIdsIdx != genotypeSubset2.SNPIdsIdx)
        {
          misc.error("Discordant SNP Idxs.", 0);
        }
        if(genotypeFull.individualIdsIdx != genotypeSubset.individualIdsIdx || genotypeFull.individualIdsIdx != genotypeSubset2.individualIdsIdx)
        {
          misc.error("Discordant Individual Idxs.", 0);
        }
        
        if( i == 0 && false)
        {
          genotypeFull.printGenotype();
        }
        misc.message << genotypeFull.genotypes->nBlockRows << " " << genotypeFull.genotypes->nBlockCols << std::endl;
      }
      }
    }
  }
  
}

void test114bis()
{
  //Tested with
  //To test: mpirun --np 6 ./dissect, mpirun --np 8 ./dissect, mpirun --np 16 ./dissect
  //mpirun --np 9 ./dissect
  //mpirun --np 4 ./dissect
  //Checked if detects a manual change in SNPs allele frequencies raises error.
  
  std::vector<int> SNPsubsets;
  SNPsubsets.push_back(100);
  SNPsubsets.push_back(101);
  SNPsubsets.push_back(102);
  SNPsubsets.push_back(103);
  SNPsubsets.push_back(104);
  SNPsubsets.push_back(250);
  SNPsubsets.push_back(500);
  SNPsubsets.push_back(750);
  SNPsubsets.push_back(900);
  
  std::vector<int> Indsubsets;
  Indsubsets.push_back(100);
  Indsubsets.push_back(500);
  Indsubsets.push_back(501);
  Indsubsets.push_back(502);
  Indsubsets.push_back(503);
  Indsubsets.push_back(504);
  Indsubsets.push_back(700);
  Indsubsets.push_back(900);
  
  for(int nr = 2; nr <= 5; nr++)
  {
    for(int nc = 2; nc <= 5; nc++)
    {
      communicator->nDefaultBlockRows = nr;
      communicator->nDefaultBlockCols = nc;
      
      misc.message << "\n-----------------------\n" << communicator->nDefaultBlockRows << " " << communicator->nDefaultBlockCols << "\n-----------------------\n" << std::endl;
      
      for(int i = 0; i<SNPsubsets.size(); i++)
      {
      for(int j = 0; j<Indsubsets.size(); j++)
      {
        std::set<std::string> keepIndividualIdsEmpty;
        std::set<std::string> keepSNPsSetEmpty;
        
        std::string subsetFile = "../testScripts/inputs/randomData/subsets2/randomdata.randomsnps.s" + i2s(SNPsubsets[i]) + ".i" + i2s(Indsubsets[j]);
        
        /////////////
        //Load subset SNPs and Inds
        std::vector< std::vector<std::string> > temptable;
        
        getTableFromFile(subsetFile + ".fam", temptable, 6);
        std::set<std::string> keepInds;
        for(int tidx = 0; tidx < temptable.size(); tidx++)
        {
          keepInds.insert(temptable[tidx][0] + "@" + temptable[tidx][1]);
        }
        getTableFromFile(subsetFile + ".bim", temptable, 6);
        std::set<std::string> keepSNPs;
        for(int tidx = 0; tidx < temptable.size(); tidx++)
        {
          keepSNPs.insert(temptable[tidx][1]);
        }
        
        std::set<std::string> keepSNPsEmpty;
        std::set<std::string> keepIndsEmpty;
        Genotype genotypeFull("../testScripts/inputs/randomData/subsets2/randomdata", keepSNPs, keepInds, GenotypeAttributes::probabilities);
        Genotype genotypeSubset(subsetFile, keepSNPsEmpty, keepIndsEmpty, GenotypeAttributes::probabilities);
        Genotype genotypeSubset2(subsetFile, keepSNPsEmpty, keepIndsEmpty, GenotypeAttributes::probabilities);
        //Genotype genotypeSubset2(subsetFile);
        genotypeSubset2.normalizeGenotypes();
        
        //Genotype genotypeSubset("../testScripts/inputs/randomData/subsets/randomdata.randomsnps." + i2s(subsets[i + 1]), keepSNPsSetEmpty, keepIndividualIdsEmpty, GenotypeAttributes::probabilities ); //null test.
        //genotypeFull.genotypes->m[0] = 0.; //null test.
        
        //////////////////////
        // Check matrices
        int r;
        r = compareMatrices(genotypeFull.genotypes, genotypeSubset.genotypes);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        r = compareMatrices(genotypeFull.missings, genotypeSubset.missings);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        
        r = compareMatrices(genotypeFull.genotypes, genotypeSubset2.genotypes);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        r = compareMatrices(genotypeFull.missings, genotypeSubset2.missings);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        
        genotypeFull.normalizeGenotypes();
        genotypeSubset.normalizeGenotypes();
        genotypeSubset2.normalizeGenotypes();
        
        r = compareMatrices(genotypeFull.genotypes, genotypeSubset.genotypes);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        
        r = compareMatrices(genotypeFull.genotypes, genotypeSubset2.genotypes);
        if(r == 0)
        {
          misc.error("Different genotype data.", 0);
        }
        
        //////////////////////
        // Add check it captures some test.
        //genotypeFull.SNPs[0].frequencies[1] = 23;
        //genotypeFull.SNPs[0].p1 = 0.1;
        
        //////////////////////
        // Check other data
        if(genotypeFull.SNPs != genotypeSubset.SNPs ) //|| genotypeFull.SNPs != genotypeSubset2.SNPs #The frequencies are not the same between bgen and bed loading files
        {
          misc.error("Discordant SNPs.", 0);
        }
        if(genotypeFull.individuals != genotypeSubset.individuals ) //|| genotypeFull.individuals != genotypeSubset2.individuals
        {
          misc.error("Discordant Individuals.", 0);
        }
        
        if(genotypeFull.SNPIds != genotypeSubset.SNPIds || genotypeFull.SNPIds != genotypeSubset2.SNPIds)
        {
          misc.error("Discordant SNP Ids.", 0);
        }
        if(genotypeFull.individualIds != genotypeSubset.individualIds || genotypeFull.individualIds != genotypeSubset2.individualIds)
        {
          misc.error("Discordant Individual Ids.", 0);
        }
        
        if(genotypeFull.SNPIdsIdx != genotypeSubset.SNPIdsIdx || genotypeFull.SNPIdsIdx != genotypeSubset2.SNPIdsIdx)
        {
          misc.error("Discordant SNP Idxs.", 0);
        }
        if(genotypeFull.individualIdsIdx != genotypeSubset.individualIdsIdx || genotypeFull.individualIdsIdx != genotypeSubset2.individualIdsIdx)
        {
          misc.error("Discordant Individual Idxs.", 0);
        }
        
        if( i == 0 && false)
        {
          genotypeFull.printGenotype();
        }
        misc.message << genotypeFull.genotypes->nBlockRows << " " << genotypeFull.genotypes->nBlockCols << std::endl;
      }
      }
    }
  }
  
}


void testDebug115()
{
  std::vector<int> SNPsubsets;
  SNPsubsets.push_back(100);
  SNPsubsets.push_back(101);
  SNPsubsets.push_back(102);
  SNPsubsets.push_back(103);
  SNPsubsets.push_back(104);
  SNPsubsets.push_back(250);
  SNPsubsets.push_back(500);
  SNPsubsets.push_back(750);
  SNPsubsets.push_back(900);
  
  std::vector<int> Indsubsets;
  Indsubsets.push_back(100);
  Indsubsets.push_back(500);
  Indsubsets.push_back(501);
  Indsubsets.push_back(502);
  Indsubsets.push_back(503);
  Indsubsets.push_back(504);
  Indsubsets.push_back(700);
  Indsubsets.push_back(900);
  
  
  GenotypeAttributes::Type genotypeFilesType = GenotypeAttributes::calls;
  
  Communicator * globalCommunicator = communicator;  
  Communicator * groupedCommunicator = new Communicator(globalCommunicator, basicGroupedCommunicator);
  communicator = groupedCommunicator;
  
  LabeledMatrix * residuals = NULL;
  Matrix * groupedCommunicatorDistributedResiduals = NULL;
  
  
  
  for(int nr = 2; nr <= 5; nr++)
  {
    for(int nc = 2; nc <= 5; nc++)
    {
      communicator->nDefaultBlockRows = nr;
      communicator->nDefaultBlockCols = nc;
      for(int i = 0; i<SNPsubsets.size(); i++)
      {
      for(int j = 0; j<Indsubsets.size(); j++)
      {
        //Get SNPs to keep for this group.
        std::string currentFile = "../testScripts/inputs/randomData/subsets2/randomdata.randomsnps.s" + i2s(SNPsubsets[i]) + ".i" + i2s(Indsubsets[j]);;
        std::vector<std::string> allSNPIds;
        if(genotypeFilesType == GenotypeAttributes::calls)
        {
          Genotype loadSNPs;
          std::set<std::string> keepSNPs;
          std::vector<int> SNPsBEDIdxsToKeep;
          loadSNPs.readBIMFile(currentFile, keepSNPs, SNPsBEDIdxsToKeep);
          allSNPIds = loadSNPs.SNPIds;
        }
        else
        {
#if defined(BGEN) && defined(ZLIB)
          Genotype loadSNPs;
          allSNPIds = loadSNPs.getBGENFileSNPs(currentFile);
#else
          misc.error("Error", 0);
#endif
        }
        int allNSNPs = allSNPIds.size();
        int nSNPsBaseInGroup = allNSNPs/groupedCommunicator->nGroups;
        int remainingSNPs = allNSNPs - (nSNPsBaseInGroup*groupedCommunicator->nGroups);
        if(remainingSNPs < 0 || remainingSNPs >= groupedCommunicator->nGroups)
        {
          misc.error("Error: An internal error was happened. Unexpected remainder when distributing SNPs between grouped processes.", 0);
        }
        std::set<std::string> currentGroupSNPIds;
        int isShift = 0;
        for(int ig = 0; ig < groupedCommunicator->nGroups; ig++)
        {
          int nSNPsInGroup = nSNPsBaseInGroup;
          if( ig < remainingSNPs )
          {
            nSNPsInGroup++;
          }
          if(ig == communicator->group)
          {
            for(int is = 0; is < nSNPsInGroup; is++)
            {
              int gsi = is + isShift; //Global SNP index.
              currentGroupSNPIds.insert(allSNPIds[gsi]);
            }
          }
          isShift += nSNPsInGroup;
        }
        if(isShift != allNSNPs)
        {
          misc.error("Error: An internal error has happened. Not all SNPs have been redistributed between grouped processes.", 0);
        }
        
        //Load genotypes
        Genotype * genotype = NULL;
        std::set<std::string> keepIndividualIds = getIndividualIdsSpecifiedByOptionKeep();
        if(genotypeFilesType == GenotypeAttributes::calls)
        {
          genotype = new Genotype(currentFile, currentGroupSNPIds, keepIndividualIds);
          genotype->normalizeGenotypes();
        }
        else
        {
          genotype = new Genotype(currentFile, currentGroupSNPIds, keepIndividualIds, GenotypeAttributes::probabilities);
          genotype->normalizeGenotypes();
        }
      
      delete genotype;
      MPI_Barrier(MPI_COMM_WORLD);
      }
      }

    }
  }
  delete communicator;
  communicator = globalCommunicator;
}

void test116()
{
  std::srand(communicator->mpiRank*std::time(0));
  
  int nelements = (rand() % 10) + 1;
  std::stringstream sdata;
  std::string data;
  
  sdata << " r" << communicator->mpiRank << ":" << nelements << "r";
  for(int i = 0; i<nelements; i++)
  {
    sdata << "-" << communicator->mpiRank << "*" << i;
  }
  data = sdata.str();
  
  std::string result = communicator->asymmetricGather(data);
  if(communicator->mpiRoot)
  {
    std::cout << result;
    std::cout << "END" << std::endl << std::endl;
  }
}

void test117()
{
  std::srand(communicator->mpiRank*std::time(0));
  
  int nelements = (rand() % 10) + 1;
  std::stringstream sdata;
  std::string data;
  
  sdata << " r" << communicator->mpiRank << ":" << nelements << "r";
  for(int i = 0; i<nelements; i++)
  {
    sdata << "-" << communicator->mpiRank << "*" << i;
  }
  sdata << std::endl;
  data = sdata.str();
  
  communicator->storeArraysMPI("test/test.dat", data);
}

void test118()
{
  std::set<std::string> keepSNPs;
  keepSNPs.insert("rs2");
  keepSNPs.insert("rs5");
  keepSNPs.insert("rs6");
  
  std::set<std::string> keepInds;
  keepInds.insert("i1@i1");
  keepInds.insert("i3@i3");
  keepInds.insert("i4@i4");
  keepInds.insert("i8@i8");
  
  for(int nr = 2; nr <= 5; nr++)
  {
    for(int nc = 2; nc <= 5; nc++)
    {
      communicator->nDefaultBlockRows = nr;
      communicator->nDefaultBlockCols = nc;
      
      misc.message << "\n-----------------------\n" << communicator->nDefaultBlockRows << " " << communicator->nDefaultBlockCols << "\n-----------------------\n" << std::endl;
      
      
      Genotype genotypeFull("./test/test3", keepSNPs, keepInds);
      genotypeFull.printGenotype();
        
    }
  }
  
  /* The genotype output should be:
   * rs2 rs5 rs6
   * i1 i1 0 0 2 -9 G G G A T T
   * i3 i3 0 0 2 -9 G G A A T T
   * i4 i4 0 0 2 -9 G G A A A T
   * i8 i8 0 0 2 -9 G G A A T T
   */
}

void test119()
{
  std::set<std::string> keepSNPs;
  /*keepSNPs.insert("rs2");
  keepSNPs.insert("rs5");
  keepSNPs.insert("rs6");*/
  
  std::set<std::string> keepInds;
  /*keepInds.insert("i1@i1");
  keepInds.insert("i3@i3");
  keepInds.insert("i4@i4");
  keepInds.insert("i8@i8");*/
  
  Genotype genotypeFull("./test/test3", keepSNPs, keepInds);
  if(communicator->mpiRoot == true)
  {
    for(int is = 0; is<genotypeFull.nSNPs; is++)
    {
      std::string a1 = genotypeFull.SNPs[is].allele1;
      std::string a2 = genotypeFull.SNPs[is].allele2;
      
      misc.message << genotypeFull.SNPs[is].name << " ";
      misc.message << genotypeFull.SNPs[is].frequencies[0] << "00" << " ";
      misc.message << genotypeFull.SNPs[is].frequencies[1] << a1 + a1 << " ";
      misc.message << genotypeFull.SNPs[is].frequencies[2] << a1 + a2 << " ";
      misc.message << genotypeFull.SNPs[is].frequencies[3] << a2 + a2 << std::endl;
    }
  }
  genotypeFull.printGenotype();
  genotypeFull.normalizeGenotypes();
  genotypeFull.genotypes->showGlobal();
  
  Matrix * m1 = new Matrix(genotypeFull.genotypes);
  m1->centerMatrixRowsColumns(column);
  m1->showGlobal("Column centered");
  
  m1->transpose(genotypeFull.genotypes);
  m1->centerMatrixRowsColumns(column);
  m1->showGlobal("Transposed and Column centered");
  
  
  Matrix * m2 = new Matrix(genotypeFull.genotypes);
  m2->centerMatrixRowsColumns(row);
  m2->showGlobal("Row centered");
}

void test120()
{
  communicator->nDefaultBlockRows = 3;
  communicator->nDefaultBlockCols = 5;
  
  bool makeRandom = true;
  
  //int nRows = 11;
  //int nCols = 7;
  int nRows = 167;
  int nCols = 193;
  
  LabeledMatrix * lm = test102CreateTestLabeledMatrix(nRows, nCols, "", "", 1.);
  if(makeRandom == true)
  {
    lm->matrix->fillWithRandom();
  }
  lm->show();
  lm->matrix->showGlobal();
  lm->save("../testScripts/results2/labeledMatrix/filter/originalMatrix/lm_r" + i2s(nRows) + "_c" + i2s(nCols) + (makeRandom?"_random":""));
}