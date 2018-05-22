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

#ifndef GLOBAL_H
#define GLOBAL_H

extern "C" {
  void Cblacs_pinfo(int*, int*);
  void Cblacs_get(int, int, int*);
  void Cblacs_gridinit(int*, const char*, int, int);
  void Cblacs_gridmap(int*, int*, int, int, int);
  void Cblacs_pcoord(int, int, int*, int*);
  void Cblacs_gridinfo(int, int*, int*, int*, int*);
  void Cblacs_gridexit(int);
  void Cblacs_barrier(int, const char*);
  void Cdgerv2d(int, int, int, double*, int, int, int);
  void Cdgesd2d(int, int, int, double*, int, int, int);
  void Cigerv2d(int, int, int, int*, int, int, int);
  void Cigesd2d(int, int, int, int*, int, int, int);
  
  int numroc_(int*, int*, int*, int*, int*);
  void descinit_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
  
  void pdgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int* ia, int* ja, int* desca, double* b, int* ib, int* jb, int* descb, double* beta, double* c, int* ic, int* jc, int* descc);
  void pdsymm_(char* side, char* uplo, int* m, int* n, double* alpha, double* a, int* ia, int* ja, int* desca, double* b, int* ib, int* jb, int* descb, double* beta, double* c, int* ic, int* jc, int* descc);
  void pdsyrk_(char* uplo, char* trans, int* n, int* k, double* alpha, double* a, int* ia, int* ja, int* desca, double* beta, double* c, int* ic, int* jc, int* descc);
  
  void pdgeadd_(char* trans, int* m, int* n, double* alpha, double* a, int* ia, int* ja, int* desca, double* beta, double* c, int* ic, int* jc, int* descc);
  
  void pdpotrf_(char* uplo, int* n, double* m, int* ia, int* ja, int* desca, int* info);
  void pdpotri_(char* uplo, int* n, double* m, int* ia, int* ja, int* desca, int* info);
  void pspotrf_(char* uplo, int* n, float* m, int* ia, int* ja, int* desca, int* info);
  void pspotri_(char* uplo, int* n, float* m, int* ia, int* ja, int* desca, int* info);
  void pdgetrf_(int* m, int* n, double* a, int* ia, int* ja, int* desca, int* ipiv, int* info);
  void pdgetri_(int* n, double* a, int* ia, int* ja, int* desca, int* ipiv, double* work, int* lwork, int* iwork, int* liwork, int* info);
  
  void pdgeqrf_(int* m, int* n, double* a, int* ia, int* ja, int* desca, double* tau, double* work, int* lwork, int* info);
  
  void pdlaset_(char* uplo, int* m, int* n, double* alpha, double* beta, double* a, int* ia, int* ja, int* desca);
  void pdlacpy_(char* uplo, int* m, int* n, double* a, int* ia, int* ja, int* desca, double* b, int* ib, int* jb, int* descb);
  void pdtran_(int* m, int* n, double* alpha, double* a, int* ia, int* ja, int* desca, double* beta, double* c, int* ic, int* jc, int* descc);
  
  void pdsyev_(char* jobz, char* uplo, int* n, double* a, int* ia, int* ja, int* desca, double* w, double* z, int* iz, int* jz, int* descz, double* work, int* lwork, int* info);
}

class Communicator;
class Options;
class Misc;

extern Communicator * communicator;
extern Options options;
extern Misc misc;

#endif
