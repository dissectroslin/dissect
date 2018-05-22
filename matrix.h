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

#ifndef MATRIX_H
#define MATRIX_H

#include "communicator.h"
#include "global.h"

#include <fstream>
#include <vector>
#include <map>

#define DESCRIPTOR_SIZE 10
#define MATRIX_DEFAULT_DISTRIBUTION cyclicDistribution

class Matrix;

enum DistributionType
{
  global,				///<The matrix is not splitted between different nodes. The m pointer in class Matrix points to the global matrix. ATTENTION: This distribution is not fully implemented and should not be used.
  cyclicDistribution,			///<The matrix is distributed between different nodes. The m pointer in class Matrix points to each local matrix.
  diagonalDistribution                  ///<The matrix is diagonal. Only the diagonal is stored in the root node. The m pointer in class Matrix is != NULL only in root node and stores the diagonal.
};

enum RowColumn
{
  row,
  column
};

struct LocalPosition
{
  int proc; //The row/col proc
  int position; //The local row/col within current proc
};

/**
 * Class for defining a submatrix
 */
class subMatrix
{
public:
  bool active;				///<When value is false, the submatrix is not active. In this case, the functions that uses them, should use the full matrix.
  
  int ir;				///<First row of the submatrix
  int ic;				///<First column of the submatrix
  int nr;				///<Number of rows of the submatrix
  int nc;				///<Number of columns of the submatrix
  
  subMatrix(int ir, int ic, int nr, int nc);
  subMatrix(Matrix * m);
  subMatrix();
};

/**
 * Class that is the basic representation of a matrix
 * This class contains the matrix data as well as functions for distributing and gather the
 * data between nodes, and performing basic matrix operations.
 */
class Matrix
{
public:
  DistributionType distribution;	///<The matrix is distributed or is in a single node? See DistributionType enum for current possibilities
  double *m;				///<A pointer to the distributed matrix. When distribution == diagonalDistribution it is NULL on all nodes except on root node where only stores the diagonal.
  double *v;				///<A distributed vector. It can be used for applying different normalizations in the distributed matrices.
  float  *mSinglePrecision;             ///<A single precission copy of m or NULL if unused
  
  int descriptor[DESCRIPTOR_SIZE];	///<A matrix descriptor according to scaLAPACK documentation.
  
  int nGlobRows;			///<The number of rows of the global matrix
  int nGlobCols;			///<The number of columns of the global matrix
  int nRows;				///<The number of rows of the local matrix. If distribution == diagonalDistribution, this is nGlobRows (and nGlobCols which are equal) on root node, 0 otherwise.
  int nCols;				///<The number of columns of the local matrix. If distribution == diagonalDistribution, this is 1 on root node, 0 otherwise.
  int nBlockRows;			///<The number of rows of the distributed matrix blocks
  int nBlockCols;			///<The number of columns of the distributed matrix blocks
  
  bool vector;				///<Is this vector or matrix? true if is vector. Not fully implemented/tested. This flag can not be always trusted.
  bool symmetric;			///<The matrix is symmetric?
  char uplo;				///<It is stored in the upper part 'U' or in the lower part 'L' or both 'B'. If symmetric is false uplo must be 'B'.
  
  /**
   * Matrix constructor that allocates memory for it
   * \param dist Define how the matrix will be distributed
   * \param ngr The number of rows of the global matrix
   * \param ngc The number of columns of the global matrix
   * \param nbr The number of rows of a matrix block
   * \param nbc The number of columns of a matrix block
   */
  Matrix(DistributionType dist, int ngr, int ngc, int nbr = communicator->nDefaultBlockRows, int nbc = communicator->nDefaultBlockCols);
  
  /**
   * Matrix constructor without allocating memory
   * \param dist Define how the matrix will be distributed
   */
  Matrix(DistributionType dist);
  
  /**
   * Matrix constructor without allocating memory
   * 
   * The matrix will be distributed according to the defined constant MATRIX_DEFAULT_DISTRIBUTION
   */
  Matrix();
  
  /**
   * Matrix constructor from an existing matrix
   * 
   * Creates a new matrix which is a copy of a previous one. The data in vector v will not be copied.
   * 
   * \param srcMatrix src matrix from which data will be copied. 
   */
  Matrix(Matrix *srcMatrix);
  
  /**
   * Matrix destructor
   */
  ~Matrix();
  
  /**
   * Sets the base default parameters for a new matrix.
   */
  void defaultMatrix(DistributionType dist);

  /**
   * Copy a matrix to this matrix
   * 
   * \param srcMatrix Pointer to the matrix which will be copied to this.
   */
  void duplicateMatrix(Matrix *srcMatrix);
  
  /**
   * Check matrix structure between this and m1
   * 
   * If nGlobRows, nGlobCols, nRows, nCols, nBlockCols, and nBlockCols differ between both matrices, an internal error is raised.
   * 
   *\param m1 matrix to compare with this matrix
   */
  void checkMatrixStructure(Matrix * m1);
  
  /**
   * Distributes a global (column-major) matrix between nodes
   * \param mGlobal a pointer to the global matrix (only used by the process defined by sourceProcessRow and sourceProcessCol)
   * \param sourceProcessRow The row of the process that have the data in mGlobal. The root process id the default.
   * \param sourceProcessCol The column of the process that have the data in mGlobal. The root process id the default.
   */
  void scatterMatrix(double *mGlobal, int sourceProcessRow = 0, int sourceProcessCol = 0);
  
  /**
   * Distributes a block of a global (column-major) matrix between nodes
   * 
   * Distributes a block of size nBlockRows*nBlockCols to the corresponding process.
   * The pointer mBlockGlobal can point to a part of a bigger size block (or the full matrix). The mumber of rows
   * of the bigger block must be indicated in the blockRowLength parameter.
   * 
   * \param mBlockGlobal a pointer to the begining of the block in the process defined by sourceProcessRow and sourceProcessCol
   * \param r row position of the block in the global matrix. It must be a multiple of nBlockRows
   * \param c column position of the block in the global matrix. It must be a multiple of nBlockCols
   * \param blockRowLength The row leading dimension of mBlockGlobal. i.e. the number of rows of this block.
   * \param sourceProcessRow The row of the process that have the data in mBlockProcessRow. The root process id the default.
   * \param sourceProcessCol The column of the process that have the data in mBlockProcessRow. The root process id the default.
   */
  void scatterBlock(double *mBlockGlobal, int r, int c, int blockRowLength, int sourceProcessRow = 0, int sourceProcessCol = 0);
  
  /**
   * Gathers a block from the distributed matrix
   * 
   * Gathers a block of size nBlockRows*nBlockCols from the corresponding porcess and storers in 
   * mBlockGlobal. The pointer can point to a bigger size block (or the full matrix). The mumber of rows
   * of the bigger block must be indicated in the blockRowLength parameter.
   * 
   * \param mBlockGlobal a pointer to the begining of the block in the root node
   * \param r row position of the block in the global matrix. It must be a multiple of nBlockRows
   * \param c column position of the block in the global matrix. It must be a multiple of nBlockCols
   * \param blockRowLength The row leading dimension of mBlockGlobal. i.e. the number of rows of this block.
   */
  void gatherBlock(double *mBlockGlobal, int r, int c, int blockRowLength);
  
  /**
   * Gathers a global (column-major) matrix between nodes
   * \param mGlobal a pointer to the global matrix (only used by the root node)
   */
  void gatherMatrix(double *mGlobal);
  
  /**
   * Distributes a vector between nodes
   * 
   * Distributes a vector between nodes. Each element of the distributed vector corresponds with a row/column of the distributed matrix.
   * This will allow to perform operations on each row or column of the matrix.
   * 
   * \param vGlobal a pointer to the vector in the root node
   * \param rowcolumn can have the values row/column depending whether the length of the vector is equal to the number of matrix columns or rows.
   */
  void scatterVector(double *vGlobal, RowColumn rowcolumn);
  
  /**
   * Distributes a vector between nodes
   * 
   * Distributes a vector between nodes. Each element of the distributed vector corresponds with a row/column of the distributed matrix.
   * This will allow to perform operations on each row or column of the matrix. This function returns the scattered vector.
   * 
   * \param vGlobal a pointer to the vector in the root node
   * \param rowcolumn can have the values row/column depending whether the length of the vector is equal to the number of matrix columns or rows.
   * \return Returns a pointer to the scattered vector. It must be freed after use.
   */
  int * scatterVectorRet(int *vGlobal, RowColumn rowcolumn);
  
  /**
   * Returns the global coordinate of a local element
   * 
   * return the global row/column from the local row/column.
   * 
   * \param rc specifies if we want convert rows or columns.
   * \param n the local row/column
   * \return the global row/column
   */
  int local2global(RowColumn rc, int n);
  
  /**
   * Returns the local coordinate position based on global position and matrix distribution
   */
  LocalPosition global2local(int globalPosition, int sizeBlock, int nProcs);
  
  
  /**
   * This function filters and resorts some matrix elements.
   * 
   * 
   * \param resultantMatrix The resultant matrix with the columns and rows filtered resorted
   * \param rowsOriginDestination A map tha maps the global row index of the original matrix to the global row index on the resultant matrix. It must be defined and equal on all nodes.
   * \param colsOriginDestination A map tha maps the global column index of the original matrix to the global column index on the resultant matrix. It must be defined and equal on all nodes.
   * \param unallocateThisMatrix When true, this matrix data will be deleted after filtering. This helps using less memory. (the matrix instance still has to be deleted).
   */
  void generalResorting(Matrix *resultantMatrix, std::map<int, int> & rowsOriginDestination, std::map<int, int> & colsOriginDestination, bool unallocateThisMatrix = false);
  
  /**
   * Filter rows and columns of current matrix
   * 
   * Filter the indicated rows and columns of the current matrix and store the resultant matrix in the matrix indicated by resultantMatrix.
   * It is responsibility of the calling funtion of setting the new symmetric, and uplo flags properly.
   * 
   * \param resultantMatrix The resultant matrix with the columns and rows filtered
   * \param rows Pointer to an ordered array which stores the rows that will be keep/filtered. Must be defined in all nodes.
   * \param nElemRows Number of elements in the rowsKeep array. Must be defined in all nodes.
   * \param cols Pointer to an ordered array which stores the columns that will be keep/filtered. Must be defined in all nodes.
   * \param nElemCols Number of elements in the colsKeep array. Must be defined in all nodes.
   * \param keep if true the rows/columns in rows/cols indicate rows/columns to be kept, otherwise they indicate rows to be filtered. Must be defined in all nodes.
   * \param unallocateThisMatrix When true, this matrix data will be deleted after filtering. This helps using less memory.
   */
  void filterRowsAndColumns(Matrix *resultantMatrix, int *rows, int nElemRows, int *cols, int nElemCols, bool keep = true, bool unallocateThisMatrix = false);
  
  /**
   * Filter rows and columns of current matrix.
   * 
   * Filter the indicated rows and columns of the current matrix and store the resultant matrix in the matrix indicated by resultantMatrix.
   * It is responsibility of the calling funtion of setting the new symmetric, and uplo flags properly. This function calls previous filterRosAndColumns (that
   * uses raw arrays). Is a helping function that deals with vectors only in root node.
   * 
   * \param resultantMatrix The resultant matrix with the columns and rows filtered
   * \param rows vector with the row indices to keep/filter (only needed in the root node).
   * \param cols vector with the column indices to keep/filter (only needed in the root node).
   * \param keep if true the rows/columns in rows/cols indicate rows/columns to be kept, otherwise they indicate rows to be filtered.
   * \param unallocateThisMatrix When true, this matrix data will be deleted after filtering. This helps using less memory.
   */
  void filterRowsAndColumns(Matrix *resultantMatrix, std::vector<int> & rows, std::vector<int> & cols, bool keep = true, bool unallocateThisMatrix = false);
  
  /**
   * This function splits and redistributes this matrix to grouped matrices defined with a new communicator.
   * 
   * It creates the matrices with dimensions based on those defined in nGlobalRowsInGroup and nGlobalColsInGroup. If there are elements on the destination matrices that remain empty,
   * (i.e. they are not receiving any element from the origin matrices) those will be filled with zero.
   * 
   * \param newCommunicator The new grouped communicator which will be used for redistributing the matrix. Its type has to be basicGroupedCommunicator.
   * \param nGlobalRowsInGroup Vector with the number of rows for the matrix in each communicator group.
   * \param nGlobalColsInGroup Vector with the number of columns for the matrix in each communicator group.
   * \param rowsOriginDestination A map tha maps the global row index of the original matrix to the global row index on the destination matrix and the group: std::pair<int, int>(row, group). It must be defined and equal on all nodes.
   * \param colsOriginDestination A map tha maps the global column index of the original matrix to the global column index on the destination matrix and the group: std::pair<int, int>(col, group). It must be defined and equal on all nodes.
   * \param unallocateThisMatrix When true, this matrix data will be deleted after filtering. This helps using less memory.
   * \return The redistributed matrices.
   */
  Matrix* redistributionToGroupedCommunicatorMatrices(Communicator * newCommunicator, std::vector<int> nGlobalRowsInGroup, std::vector<int> nGlobalColsInGroup, std::map<int, std::pair<int, int> > & rowsOriginDestination, std::map<int, std::pair<int, int> > & colsOriginDestination, bool unallocateThisMatrix);
  
  /**
   * This function distributes this matrix to different groups of a distributed new communicator.
   * 
   * \param groupedCommunicator The new grouped communicator which will be used for redistributing the matrix. Its type has to be basicGroupedCommunicator.
   * \return The redistributed matrix.
   */
  Matrix* copyToGroupedCommunicator(Communicator * groupedCommunicator);
  
  /**
   * Stores the matrix in a file (column major)
   * 
   * Stores the matrix data in a file. The function can be called many times from different matrices the the same file
   * for storing all the matrices in one file.
   * 
   * \param file the already opened file object where the matrix will be stored
   */
  void writeMatrixFile(std::ofstream & file);
  
  /**
   * ATTENTION this function is not tested.
   */
  void writeMatrixFilev2(std::ofstream & file);
  
  /**
   * Loads the matrix from a file (column major)
   * 
   * Loads the matrix data from a file. The function can be called many times from different matrices with the same file
   * for loading all the matrices from one file.
   * 
   * \param file the already opened file object from where the matrix will be loaded
   * \param ngr The number of rows of the global matrix
   * \param ngc The number of solumns of the global matrix
   * \param nbr The number of rows of a matrix block
   * \param nbc The number of columns of a matrix block
   */
  void readMatrixFile(std::ifstream & file, int ngr, int ngc, int nbr, int nbc);
  
  /**
   * ATTENTION this function is not tested.
   */
  void readMatrixFilev2(std::ifstream & file, int ngr, int ngc, int nbr, int nbc);
  
  void readMatrixMPI(std::string fname, int offset = 0);
  void writeMatrixMPI(std::string fname, char * header = NULL, int offset = 0);  //ATTENTION This function can lead to wrong results in some file systems (e.g. sshfs or NFS?).
  
  /**
   * This function packs two symmetric matrices into a new one
   * 
   * The function gets matrix m1 and stores it in the lower part of this and m2 and stores it in the upper part of this.
   * This matrix have the same number of columns of m1 and m2, but one more row.
   * m1 and m2 must be square, symmetric and they must have the same dimensions.
   * 
   * \param m1 first matrix to pack
   * \param m2 second matrix to pack
   */
  void packMatrices(Matrix * m1, Matrix * m2);
  
  /**
   * This function unpacks two symmetric matrices from this
   * 
   * The function gets matrix m1 from the lower part of this and m2 from the upper part of this.
   * This matrix have the same number of columns of m1 and m2, but one more row (i.e. this->nGlobRows - 1 == m1->nGlobRows == m2->nGlobRows)
   * m1 and m2 will be square, symmetric and will have the same dimensions. Both will be stored on the lower part.
   * 
   * \param[out] m1 first unpacked matrix
   * \param[out] m2 second unpacked matrix
   * \param unallocateThisMatrix When true, this matrix data will be deleted after filtering. This helps using less memory.
   */
  void unpackMatrices(Matrix * m1, Matrix * m2, bool unallocateThisMatrix = false);
 
  /**
   * Allocate memory for matrix data
   */
  void allocateMemory();
  
  /**
   * Unallocate memory for matrix data
   */
  void unallocateMemory();
  
  /**
   * Initiate matrix parameters and allocate memory
   * 
   * \param ngr The number of rows of the global matrix
   * \param ngc The number of columns of the global matrix
   * \param nbr The number of rows of a matrix block
   * \param nbc The number of columns of a matrix block
   * \param newDist The new distribution type of the matrix.
   * \param allocateMemoryNow Memory for this matrix will be allocated when true. 
   */
  void initParameters(int ngr, int ngc, int nbr, int nbc, DistributionType newDist = cyclicDistribution, bool allocateMemoryNow = true);
  
  /**
   * Fill the matrix with a constant value
   * 
   * \param value The value used for filling the matrix
   */
  void fillWithConstant(double value);
  
  /**
   * Initializes the offdiagonal elements of a matrix to one value and the diagonal elements to another.
   * 
   * \param diagonal The value used for filling the diagonal
   * \param background The value used for filling the offdiagonal elements. If not specified: default value 0.
   */
  void fillDiagonal(double diagonal, double background = 0.);
  
  /**
   * Returns a local matrix element
   * 
   * 
   * \param r The row of the element
   * \param c The column of the element
   * \param defaultValue The default value returned if local==false
   * \param[out] local Is set to true if the element is in the local matrix, false otherwise.
   * \param useSinglePrecision if true, get value from this->mSinglePrecision instead of this->m.
   * \return The value in row r and column c. Returns 1. if this element is not in the local matrix (local=false).
   */
  double matrixElement(int r, int c, double defaultValue, bool *local, bool useSinglePrecision = false);
  
  /**
   * Gathers a row block from the distributed matrix ATTENTION Needs testing
   * 
   * Gathers a block of size nBlockRows x nGlobCols from the distributed matrix.
   * 
   * \param r row position of the block in the global matrix. It must be a multiple of nBlockRows
   * \param blockRow a pointer to the begining of the block in the root node where block will be stored. The block allocated memory must be of size nGlobRows x nBlockCols
   * \param nr returns the number of rows in this block.
   */
  void gatherRowBlock(int r, double * blockRow, int * nr);
  
  /**
   * Gathers a column block from the distributed matrix ATTENTION Needs testing
   * 
   * Gathers a block of size nGlobRows x nBlockCols from the distributed matrix.
   * 
   * \param c col position of the block in the global matrix. It must be a multiple of nBlockCols
   * \param blockCol a pointer to the begining of the block in the root node where block will be stored. The block allocated memory must be of size nBlockRows x nGlobCols
   * \param nc returns the number of columns in this block.
   */
  void gatherColBlock(int c, double * blockCol, int * nc);
  
  /**
   * Multiply two matrices/vectors and store the result in this matrix
   * 
   * Performs the matrix/vector product op(m1)*op(m2). op() transposes the matrix/vector depending on t1 and t2 values.
   * 
   * \param m1 Matrix 1
   * \param t1 'T' or 'N'. 'T' indicates the transpose of matrix m1 is used.
   * \param m2 Matrix 2
   * \param t2 'T' or 'N'. 'T' indicates the transpose of matrix m2 is used.
   * \param scale Constant which will multiply the product.
   */
  void multiply(Matrix * m1, char t1, Matrix * m2, char t2, double scale = 1., subMatrix smr = subMatrix());
  
  /**
   * Multiply two matrices and store the result in this matrix
   * 
   * Performs the matrix product op(m1)*op(m2). op() transposes the matrix depending on t1 and t2 values.
   * 
   * \param m1 Matrix 1
   * \param t1 'T' or 'N'. 'T' indicates the transpose of matrix m1 is used.
   * \param m2 Matrix 2
   * \param t2 'T' or 'N'. 'T' indicates the transpose of matrix m2 is used.
   */
  void multiplyMatrixMatrix(Matrix * m1, char t1, Matrix * m2, char t2, double scale, subMatrix smr);
  
  /**
   * Multiply A matrix with a vector (must be implemented)
   * 
   * Performs the matrix product op(m1)*op(m2). op() transposes the matrix depending on t1 and t2 values.
   * 
   * \param m1 Matrix/Vector
   * \param t1 'T' or 'N'. 'T' indicates the transpose of matrix m1 is used.
   * \param m2 Matrix/Vector
   * \param t2 'T' or 'N'. 'T' indicates the transpose of matrix m2 is used.
   */
  void multiplyMatrixVector(Matrix * m1, char t1, Matrix * m2, char t2, double scale, subMatrix smr);
  
  /**
   * Multiply two matrices and store the result in this matrix. On or both multiplied matrices must be diagonally distributed.
   * 
   * Performs the matrix product op(m1)*op(m2). op() transposes the matrix depending on t1 and t2 values.
   * 
   * \param m1 Matrix 1
   * \param t1 'T' or 'N'. 'T' indicates the transpose of matrix m1 is used.
   * \param m2 Matrix 2
   * \param t2 'T' or 'N'. 'T' indicates the transpose of matrix m2 is used.
   */
  void multiplyDiagonalMatrixMatrix(Matrix * m1, char t1, Matrix * m2, char t2, double scale, subMatrix smr);
  
  /**
   * Adds this and m1 matrices and store the result in this matrix.
   * 
   * Performs the matrix sum sub(this->m) = beta*sub(this->m) + alpha*op(sub(m1)). The summation can be performed to all the matrix or using submatrices.
   * If the summation must be performed using submatrices, both smt and sm1 must be specified. If this is diagonal distributed, but m1 is not diagonal, then this
   * will be redistributed as cyclic distributed.
   * 
   * \param m1 Matrix 1
   * \param beta scaling factor for this matrix
   * \param alpha scaling factor for the m1 matrix
   * \param smt Indicates which submatrix of this matrix must be used
   * \param sm1 Indicates which submatrix of m1 matrix must be used
   */
  void add(Matrix * m1, double beta = 1., double alpha = 1., subMatrix smt = subMatrix(), subMatrix sm1 = subMatrix());
  
  /**
   * Adds this and m1 matrices and store the result in this matrix. Helping function of add() which sums diagonal distributed matrices
   * 
   * Performs the matrix sum sub(this->m) = beta*sub(this->m) + alpha*op(sub(m1)). The summation can be performed to all the matrix or using submatrices.
   * If the summation must be performed using submatrices, both smt and sm1 must be specified. If this is diagonal distributed, but m1 is not diagonal, then this
   * will be redistributed as cyclic distributed.
   * 
   * \param m1 Matrix 1
   * \param beta scaling factor for this matrix
   * \param alpha scaling factor for the m1 matrix
   * \param smt Indicates which submatrix of this matrix must be used
   * \param sm1 Indicates which submatrix of m1 matrix must be used
   */
  void addDiagonalMatrix(Matrix * m1, double beta, double alpha, subMatrix smt, subMatrix sm1);
  
  /**
   * Adds this and m1 matrices and store the result in this matrix. Helping function of add() which sums non-diagonal distributed matrices.
   * 
   * Performs the matrix sum sub(this->m) = beta*sub(this->m) + alpha*op(sub(m1)). The summation can be performed to all the matrix or using submatrices.
   * If the summation must be performed using submatrices, both smt and sm1 must be specified.
   * 
   * \param m1 Matrix 1
   * \param beta scaling factor for this matrix
   * \param alpha scaling factor for the m1 matrix
   * \param smt Indicates which submatrix of this matrix must be used
   * \param sm1 Indicates which submatrix of m1 matrix must be used
   */
  void addMatrix(Matrix * m1, double beta, double alpha, subMatrix smt, subMatrix sm1);
  
  /**
   * Function that inverts the current matrix
   * 
   * The matrix to be inverted must be symmetric. The function checks the class flag symmetric before inversion.
   * The function also returns the natural logarithm of the determinant if logDeterminant != NULL.
   * 
   * \param logDeterminant If != NULL, the natural logarithm of the determinant is computed and returned to the var pointed by the pointer.
   * \param useSinglePrecision if true, single precision inversion is performed.
   * \return true if success false otherwise
   */
  bool symmetricInvert(double * logDeterminant = NULL, bool useSinglePrecision = false);
  
  /**
   * Function that inverts the current matrix
   * 
   * The function also returns the natural logarithm of the determinant if logDeterminant != NULL. ATTENTION The determinant computation is not exact.
   * Probably it is only correct when the matrix is 'nearly' positive definite. It is using the GCTA approach.
   * 
   * \param logDeterminant If != NULL, the natural logarithm of the determinant is computed and returned to the var pointed by the pointer.
   * \return true if success false otherwise
   */
  bool invert(double * logDeterminant = NULL);
  
  /**
   * Compute the eigendecomposition of this matrix
   * 
   * Compute the eigedecomposition of this matrix. The matrix must be a square symmetric matrix. The values on the matrix will be deleted.
   * 
   * \param eigenValues pointer to a matrix where the eigenvalues will be stored. The matrix is a diagonal distributed matrix, and the eigenvalues are in the diagonal.
   * \param eigenVectors pointer to a matrix where the eigenvectors will be stored.
   */
  void eigenDecomposition(Matrix * eigenValues, Matrix * eigenVectors);
  
  /**
   * Bend this matrix to be positive definite ATTENTION This function is not tested extensively. Maybe need some corrections.
   */
  void bendMatrix();
  
  /**
   * Performs a QR decomposition in place.
   * 
   * Performs a QR decomposition in place. The current matrix will be substitued by the decomposition. See scaLAPACK documentation for details.
   */
  void QRDecomposition();
  
  /**
   * Returns linearly dependent columns of the matrix.
   * 
   * The function requires that the number of columns is less or equal than the number of rows. The function performs a QR decomposition in place.
   * Thus the current matrix will be deleted. The independent columns are found by looking at the diagonal elements of the R matrix.
   * A column it is considered dependent if the corresponding diagonal element i smaller than threshold parameter.
   * 
   * \param threshold A column is considered dependent if the corresponding diagonal element is below the threshold.
   * \return Returns a vector with the indices of the matrix linearly dependent columns.
   */
  std::vector<int> getDependentColumns(double threshold = 1e-5);
  
  /**
   * Function that returns the trace of the matrix
   * 
   * This function returns the trace of the matrix. The matrix must be a square matrix.
   * 
   * \return trace of the matrix
   */
  double trace();
  
  /**
   * Function that returns the average of the sum of all matrix elements
   * 
   * \return average of all element
   */
  double elementsAverage();
  
  /**
   * Function that computes the exponential of the matrix elements multiplied by a factor
   * 
   * Computes this->m(i,j) = exp(this->m(i,j)*alpha)
   * 
   * \param alpha The factor used in the exponential that multiplies the matrix
   */
  void applyExponentialOperator(double alpha);
  
  /**
   * Function that computes the inverse logistic function on the matrix elements
   * 
   * ATTENTION: This function has not been checked.
   * 
   * this->m(i,j) = 1./(1 + exp( -this->m(i,j) ))
   */
  void applyInverseLogistic();
  
  /**
   * Function that returns the diagonal elements of the matrix
   * 
   * This function returns the diagonal elements of the matrix. It can only be used in square matrices (if allowNonSquareMatrices is set to false). I can also only
   * be used on matrices with symmetric block sizes. A copy is returned on root process.
   * 
   * \param allowNonSquareMatrices If set to true, non square matrices are allowed.
   * \return diagonal of the matrix in the root process
   */
  std::vector<double> diagonal(bool allowNonSquareMatrices = false);
  
  /**
   * Set the diagonal on a matrix
   * 
   * The matrix must be a square matrix.
   * 
   * \param diagonal Pointer to an array with the diagonal elements. Only needs to be set on root process.
   * \param nElements The number of elements on diagonal parameter. These must be equal to the idmensions of the destination matrix.
   */
  void setDiagonal(double * diagonal, int nElements);
  
  /**
   * Function that returns the trace of the matrix pruduct this*m1;
   * 
   * This function returns the trace of the product tr(this*m1). It is much more efficient that computing the product and then use the function trace().
   * 
   * \param m1 The matrix that will be multiplied with this before computing the trace.
   * \return trace of the matrix
   */
  double traceOfMatrixProduct(Matrix * m1);
  
  /**
   * Function that computes diagonal of the matrix product A*B*At;
   * 
   * The result is a diagonal matrix which will be stored in this. B has to be a square matrix. This is much more efficient than performing the full
   * product if The number of columns of A >> than the number of rows.
   * 
   * \param A The A matrix
   * \param B The B matrix. B has to be a square matrix.
   * \return trace of the matrix
   */
  double diagonalOfABAt(Matrix * A, Matrix * B);
  
  /**
   * Function for transposing a matrix
   * 
   * This function gets a matrix as an argument and stores its transpose in the current matrix.
   * 
   * \param m1 Matrix to transpose
   */
  void transpose(Matrix * m1);
  
  /**
   * Copies the lower or upper triangular part of a matrix to the other triangular part.
   */
  void symmetrizeTriangularMatrix();
  
  /**
   * Performs an element wise multiplication
   * 
   * this(i,j) = this(i,j)*m1(i,j). this and m1 must have the same dimensions and block sizes.
   * 
   * \param m1 matrix used for multiplication;
   * \param scale Constant which will multiply the product.
   */
  void elementWiseMultiplication(Matrix * m1, double scale = 1.);
  
  /**
   * Performs an element wise division
   * 
   * this(i,j) = this(i,j)/m1(i,j). this and m1 must have the same dimensions and block sizes.
   * 
   * \param m1 matrix used for division;
   * \param scale Constant which will multiply the product.
   */
  void elementWiseDivision(Matrix * m1, double scale = 1.);
  
  /**
   * Multiply all elements of the matrix by a constant factor.
   * 
   * \param scale factor used to multiply all matrix elements.
   */
  void scaleBy(double scale);
  
  /**
   * Create an intersection matrix from a vector.
   * 
   * Create an intersection matrix from a generating vector. If the vector element v(i) is equal to vector element v(j),
   * then the matrix element m(i,j) is set to onIntersection, onNoIntersection value will be used otherwise.
   * 
   * \param categories The vector for generating the intersection matrix. Only needs to be defined in root.
   * \param onIntersection The value put on element ij of the matrix when vector element i equals vector element j
   * \param onNoIntersection The value put on element ij of the matrix when vector element i differs from vector element j
   */
  void makeIntersectionMatrix(std::vector<int> categories, double onIntersection = 1., double onNoIntersection = 0.);
  
  /**
   * Create a matrix of differences between the different elements of two vectors.
   * 
   * Create a matrix of differences between the different elements of a vector. Matrix element m(i ,j) = v1(i) - v2(j).
   * 
   * \param v1 The first vector for generating the matrix. Only needs to be defined in root.
   * \param v2 The second vector for generating the matrix. Only needs to be defined in root.
   */
  void makeDifferenceMatrix(std::vector<double> v1, std::vector<double> v2);
  
  /**
   * Copy the data in this->m to this->mSinglePrecision
   * 
   * Copy the data in this->m to this->mSinglePrecision. If it is necessary (this->mSinglePrecission == NULL),
   * the function will allocate the needed memory.
   */
  void copyDoubleToSingle();
  
  /**
   * Copy the data in this->mSinglePrecision to this->m
   * 
   * \param unallocateSinglePrecisionMatrix if true, the memory of this->mSinglePrecission will be unallocated.
   */
  void copySingleToDouble(bool unallocateSinglePrecisionMatrix = false);
  
  /**
   * Return the global indices of the elements in local matrix greather than a particular threshold
   * 
   * Returns indices only of the local part of the matrix in the current process. Note that since an element has a row and a column,
   * the number of indices in rows and cols must be the same.
   * 
   * \param threshold The elements must be greather than this threshold
   * \param rows Indices of the rows of the elements greather than threshold
   * \param cols Indices of the columns of the elements greather than threshold
   */
  void getGlobalIndexElementsGreatherThan(double threshold, std::vector<int> & rows, std::vector<int> & cols);
  
  /**
   * Returns the global indices of the matrix elements outside a defined range in the root process.
   * 
   * Note that since an element has a row and a column, the number of indices in rows and cols must be the same.
   * 
   * \param lowerThreshold Lower threshold defining the range. Indices of all elements below this threshold will be returned.
   * \param upperThreshold Upper threshold defining the range. Indices of all elements above this threshold will be returned.
   * \param rows Row indices of the elements.
   * \param cols Column indices of the elements.
   */
  void getGlobalIndexOutsideRange(double lowerThreshold, double upperThreshold, std::vector<int> & idxGlobalRows, std::vector<int> & idxGlobalCols);
  
  /**
   * Returns the global indices of the matrix elements inside a defined range in the root process.
   * 
   * Note that since an element has a row and a column, the number of indices in rows and cols must be the same.
   * 
   * \param lowerThreshold Lower threshold defining the range. Indices of all elements below this threshold will be returned.
   * \param upperThreshold Upper threshold defining the range. Indices of all elements above this threshold will be returned.
   * \param rows Row indices of the elements.
   * \param cols Column indices of the elements.
   */
  void getGlobalIndexInsideRange(double lowerThreshold, double upperThreshold, std::vector<int> & idxGlobalRows, std::vector<int> & idxGlobalCols);
  
  /**
   * Joins two matrices in a new one.
   * 
   * This function joins two matrices in a new one. To this end, this matrix is resized to the minimum dimensions needed for fitting the two matrices.
   * Each matrix will be placed in a position defined by sm1 and sm2. Overlapping of both matrices is allowed. In this case, the value of the overlaping
   * region will be the values of the second matrix. As a function of source matrix properties, the resultant matrix will be distributued diagonal or cyclic.
   * If source matrices are bot diagonal, they must be placed on the diagonal and tha backgroundValue is 0, then the resultant matrix is diagonal distributed. 
   * It will be cyclic distributed otherwise.
   * 
   * \param m1 Matrix 1
   * \param sm1 region where Matrix 1 will be placed. nr and nc of sm1 must be equal to the dimensions of m1;
   * \param m2 Matrix 2
   * \param sm2 region where Matrix 2 will be placed. nr and nc of sm2 must be equal to the dimensions of m2;
   * \param backgroundValue This matrix will have this value if there are regions filled by m1 or m2.
   */
  void joinMatrices(Matrix * m1, subMatrix sm1, Matrix * m2, subMatrix sm2, double backgroundValue = 0.);
  
  /**
   * Joins two matrices in a new one. Only for diagonal distributed matrices.
   * 
   * Only for diagonal distributed matrices. The place of the matrices must be also in the diagonal. This function joins two matrices in a new one. To this end, this matrix is resized to the minimum dimensions needed for fitting the two matrices.
   * Each matrix will be placed in a position defined by sm1 and sm2. Overlapping of both matrices is allowed. In this case, the value of the overlaping
   * region will be the values of the second matrix.
   * 
   * \param m1 Matrix 1
   * \param sm1 region where Matrix 1 will be placed. nr and nc of sm1 must be equal to the dimensions of m1;
   * \param m2 Matrix 2
   * \param sm2 region where Matrix 2 will be placed. nr and nc of sm2 must be equal to the dimensions of m2;
   * \param backgroundValue This matrix will have this value if there are regions filled by m1 or m2. In this functions it must be 0.
   */
  void joinDiagonalMatrices(Matrix * m1, subMatrix sm1, Matrix * m2, subMatrix sm2, double backgroundValue = 0.);
  
  /**
   * Joins two matrices in a new one. Only for cyclic distributed matrices.
   * 
   * Only for cyclic distributed matrices. This function joins two matrices in a new one. To this end, this matrix is resized to the minimum dimensions needed for fitting the two matrices.
   * Each matrix will be placed in a position defined by sm1 and sm2. Overlapping of both matrices is allowed. In this case, the value of the overlaping
   * region will be the values of the second matrix.
   * 
   * \param m1 Matrix 1
   * \param sm1 region where Matrix 1 will be placed. nr and nc of sm1 must be equal to the dimensions of m1;
   * \param m2 Matrix 2
   * \param sm2 region where Matrix 2 will be placed. nr and nc of sm2 must be equal to the dimensions of m2;
   * \param backgroundValue This matrix will have this value if there are regions filled by m1 or m2.
   */
  void joinGeneralMatrices(Matrix * m1, subMatrix sm1, Matrix * m2, subMatrix sm2, double backgroundValue = 0.);
  
  /**
   * Joins two matrices vertically.
   * 
   * This function joins two matrices vertically. To this end, this matrix is resized to the minimum dimensions needed for fitting the two matrices.
   * m1 and m2 must have the same number of columns.
   * 
   * \param m1 Matrix 1
   * \param m2 Matrix 2
   */
  void joinMatricesVertically(Matrix * m1, Matrix * m2);
  
  /**
   * Joins two matrices horizontally.
   * 
   * This function joins two matrices horizontally. To this end, this matrix is resized to the minimum dimensions needed for fitting the two matrices.
   * m1 and m2 must have the same number of rows.
   * 
   * \param m1 Matrix 1
   * \param m2 Matrix 2
   */
  void joinMatricesHorizontally(Matrix * m1, Matrix * m2);
  
  /**
   * Returns this matrix as a standard vector in root process
   * 
   * \param globalVector On exit, stores this matrix as a std::vector in root process. empty vector on other processes.
   */
  void matrixToStandardVector(std::vector< std::vector<double> > & globalVector);
  
  /**
   * Returns this matrix as a standard "flat" vector in root process
   * 
   * \param globalVector On exit, stores this matrix as a std::vector in root process. empty vector on other processes.
   */
  void matrixToStandardVector(std::vector<double> & globalVector);
  
  /**
   * Standardizes the columns/rows of the matrix. ATTENTION: This function has not been tested, yet. Especially the ddof parameter != 0.
   * 
   * \param rowcolumn Indicates whether the rows or columns will be standardized (e.g. column standardizes the columns).
   * \param ddof Divisor correction used for computing the std. i.e. sqrt(sum( (v - mean(v))^2 )/(N-ddof))
   */
  void standardizeMatrix(RowColumn rowcolumn, int ddof = 0);
  
  /**
   * Remove the colum/row mean to the matrix columns/rows.
   * 
   * \param rowcolumn Indicates whether the rows or columns will be standardized (e.g. column standardizes the columns).
   */
  void centerMatrixRowsColumns(RowColumn rowcolumn);
  
  /**
   * ATTENTION this function do not work as expected. pdlacpy not communicates???
   * 
   * (ATTENTION this function do not work as expected. pdlacpy not communicates???) Performs the matrix copy sub(this->m) = sub(m1).
   * 
   * \param m1 Matrix 1
   * \param smt Indicates which submatrix of this matrix must be used
   * \param sm1 Indicates which submatrix of m1 matrix must be used. If not specified, all m1 matrix will be used.
   */
  void copySubMatrix(Matrix * m1, subMatrix smt, subMatrix sm1 = subMatrix());
  
  /**
   * Function for debugging pourposes that prints the local matrix
   */
  void showPartial(RowColumn rowcolumn, bool showv = false, bool showContext = false);
  
  /**
   * Function for debugging pourposes that gathers and prints the global matrix
   */
  void showGlobal(std::string name = "", bool symmetrize = true, int precision = 5, double zeroThreshold = 0.);
  
  /**
   * Function for debugging pourposes that reads and scatters a matrix
   */
  void debugRead(std::string fname, int ttnr, int ttnc,  int sourceProcessRow = 0, int sourceProcessCol = 0, int ttnblockr = communicator->nDefaultBlockRows, int ttnblockc = communicator->nDefaultBlockCols);
  
  /**
   * Function for debugging pourposes that gathers and writes a matrix
   */
  void debugWrite(std::string fname);
  
  /**
   * Function for debugging pourposes that fills the matrix with random numbers.
   */
  void fillWithRandom(double min = -10., double max = 10., long * seed = NULL);

  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Old deprecated functions
  
  
  /**
   * Filter rows and columns of current matrix
   * 
   * Filter the indicated rows and columns of the current matrix and store the resultant matrix in the matrix indicated by resultantMatrix.
   * It is responsibility of the calling funtion of setting the new symmetric, and uplo flags properly.
   * 
   * \param resultantMatrix The resultant matrix with the columns and rows filtered
   * \param rows Pointer to an ordered array which stores the rows that will be keep/filtered. Must be defined in all nodes.
   * \param nElemRows Number of elements in the rowsKeep array. Must be defined in all nodes.
   * \param cols Pointer to an ordered array which stores the columns that will be keep/filtered. Must be defined in all nodes.
   * \param nElemCols Number of elements in the colsKeep array. Must be defined in all nodes.
   * \param keep if true the rows/columns in rows/cols indicate rows/columns to be kept, otherwise they indicate rows to be filtered. Must be defined in all nodes.
   */
  void filterRowsAndColumnsOld(Matrix *resultantMatrix, int *rows, int nElemRows, int *cols, int nElemCols, bool keep = true);
  
  /**
   * Filter rows and columns of current matrix.
   * 
   * Filter the indicated rows and columns of the current matrix and store the resultant matrix in the matrix indicated by resultantMatrix.
   * It is responsibility of the calling funtion of setting the new symmetric, and uplo flags properly. This function calls previous filterRosAndColumns (that
   * uses raw arrays). Is a helping function that deals with vectors only in root node.
   * 
   * \param resultantMatrix The resultant matrix with the columns and rows filtered
   * \param rows vector with the row indices to keep/filter (only needed in the root node).
   * \param cols vector with the column indices to keep/filter (only needed in the root node).
   * \param keep if true the rows/columns in rows/cols indicate rows/columns to be kept, otherwise they indicate rows to be filtered.
   */
  void filterRowsAndColumnsOld(Matrix *resultantMatrix, std::vector<int> & rows, std::vector<int> & cols, bool keep = true);
  
  /**
   * This function unpacks two symmetric matrices from this
   * 
   * The function gets matrix m1 from the lower part of this and m2 from the upper part of this.
   * This matrix have the same number of columns of m1 and m2, but one more row (i.e. this->nGlobRows - 1 == m1->nGlobRows == m2->nGlobRows)
   * m1 and m2 will be square, symmetric and will have the same dimensions. Both will be stored on the lower part.
   * 
   * \param[out] m1 first unpacked matrix
   * \param[out] m2 second unpacked matrix
   */
  void unpackMatricesOld(Matrix * m1, Matrix * m2);

};




#endif
