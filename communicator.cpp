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

#include "communicator.h"
#include "global.h"
#include "misc.h"
#include "options.h"

#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <stdlib.h>

Communicator::Communicator(int argc, char **argv, CommunicatorType communicatorType)
{
  this->localargc = argc;
  this->localargv = argv;
  
  // Initiate MPI
  if(communicatorType == globalCommunicator)
  {
    int tmp;
    tmp = MPI_Init(&argc, &argv);
    if (tmp != MPI_SUCCESS) {
      error("Error: MPI can not be started. Terminating.");
    }
    this->type = globalCommunicator;
    this->mpiCommunicator = MPI_COMM_WORLD; 
  
    MPI_Comm_rank(MPI_COMM_WORLD, &this->mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &this->mpiNumTasks);
    MPI_Get_processor_name(this->hostName, &this->lenHostName);
    this->setMPIRoot();
    
    // Initiate Cblas context 
    // For the moment I assume that the number of processes have a integer square root 
    //int temp = ceil(sqrt(float(this->mpiNumTasks)));
    //this->nProcCols = temp;
    //this->nProcRows = temp;

    // allow a non-square number of processes
    //this->nProcCols=sqrt(this->mpiNumTasks);
    //this->nProcRows=(this->mpiNumTasks+this->nProcCols-1)/this->nProcCols; //This solution have problemes for some numbers: e.g. 60
    
    this->nProcCols = floor(sqrt(this->mpiNumTasks));
    while(this->mpiNumTasks % this->nProcCols != 0)
    {
      this->nProcCols--;
      if(this->nProcCols <= 1)
      {
        error("Error: The number of MPI tasks can not be a prime number.");
      }
    }
    this->nProcRows = this->mpiNumTasks/this->nProcCols;
    if(this->nProcCols*this->nProcRows != this->mpiNumTasks)
    {
      error("Error: The number of MPI tasks can not be a prime number.");
    }

    //this->nDefaultBlockRows = 3;
    this->nDefaultBlockRows = 64;
    this->nDefaultBlockCols = 64;


    //allow block size to be overwritten by an environment variable
    char* tmpBS;
    tmpBS = getenv ("BLACS_BLOCKSIZE");
    if (tmpBS!=NULL){
      this->nDefaultBlockRows = atoi(tmpBS);
      this->nDefaultBlockCols = atoi(tmpBS);
      if(this->nDefaultBlockRows <= 0 || this->nDefaultBlockCols <= 0)
      {
        error("Error: The environment variable BLACS_BLOCKSIZE is not a valid integer value.");
      }
    }



    Cblacs_pinfo(&this->myId, &this->nProc);
    Cblacs_get(0, 0, &this->context);
    Cblacs_gridinit(&this->context, "Row-major", this->nProcRows, this->nProcCols); //Attention: This must be Row-major for compatibility with MPI IO
    Cblacs_pcoord(this->context, this->myId, &this->myRow, &this->myCol);
    
    this->group = 0;
    this->nGroups = 1;
    this->groupNProcRows.clear();
    this->groupNProcCols.clear();
    this->groupNProcRows.push_back(this->nProcRows);
    this->groupNProcCols.push_back(this->nProcCols);
    this->basicGroupSize = this->mpiNumTasks;
  }
  else if(communicatorType == basicGroupedCommunicator)
  {
    error("Error: An internal error was happened. This communicator type has to be initialized from another global communicator.");
  }
  
  if(this->mpiRoot == true && (this->myRow != 0 || this->myCol != 0))
  {
    error("Error: An internal error happened after initializing BLACS grid.");
  }
}

Communicator::Communicator(Communicator * baseCommunicator, CommunicatorType communicatorType)
{
  this->localargc = baseCommunicator->localargc;
  this->localargv = baseCommunicator->localargv;
  
  if(communicatorType == globalCommunicator)
  {
    error("Error: An internal error was happened. Please, initiate this kind of communicator using other constructors.");
  }
  else if(communicatorType == basicGroupedCommunicator)
  {
    this->type = basicGroupedCommunicator;
    
    //MPI communicator
    
    MPI_Group globalGroup, localGroup;
    
    std::vector<int> localRanks;
    std::vector< std::vector<int> > allRanks;
    this->basicGroupSize = options.numberOfProcessesDefaultMPIGroup;
    for(int i = 0; i < baseCommunicator->mpiNumTasks; i++)
    {
      if( (i % options.numberOfProcessesDefaultMPIGroup) == 0 && i != 0 )
      {
        allRanks.push_back(localRanks);
        localRanks.clear();
      }
      localRanks.push_back(i);
      if(baseCommunicator->mpiRank == i)
      {
        this->group = allRanks.size();
      }
    }
    allRanks.push_back(localRanks);
    this->nGroups = allRanks.size();
    
    if( this->group >= this->nGroups )
    {
      error("Error: An internal error was happened while creating MPI group communicators.");
    }
    this->mpiNumTasks = allRanks[this->group].size();

    MPI_Comm_group(baseCommunicator->mpiCommunicator, &globalGroup);
    MPI_Group_incl(globalGroup, this->mpiNumTasks, &(allRanks[this->group][0]), &localGroup);
    MPI_Comm_create(baseCommunicator->mpiCommunicator, localGroup, &this->mpiCommunicator);
    MPI_Group_rank (localGroup, &this->mpiRank);
    
    this->setMPIRoot();
    
    
    //BLACs context
    
    std::vector<int> allNProcRows;
    std::vector<int> allNProcCols;
    for(int i = 0; i < allRanks.size(); i++)
    {
      int tempCols = floor(sqrt(allRanks[i].size()));
      while(int(allRanks[i].size()) % tempCols != 0)
      {
        tempCols--;
        if(tempCols < 1)
        {
          error("Error: An internal error was happened while creating a BLACs grid.");
        }
      }
      int tempRows = int(allRanks[i].size())/tempCols;
      if(tempCols*tempRows != allRanks[i].size())
      {
        error("Error: An internal error was happened while creating a BLACs grid.");
      }
      allNProcRows.push_back(tempRows);
      allNProcCols.push_back(tempCols);
    }
    this->groupNProcRows = allNProcRows;
    this->groupNProcCols = allNProcCols;
    
    Cblacs_pinfo(&this->myId, &this->nProc);
    std::vector<int> BLACSContexts;
    this->groupNProcRows.clear();
    this->groupNProcCols.clear();
    for(int i = 0; i < allRanks.size(); i++)
    {
      int temp;
      Cblacs_get(0, 0, &temp);
      BLACSContexts.push_back(temp);
    }
    for(int i = 0; i < allRanks.size(); i++)
    {
      Cblacs_gridmap (&BLACSContexts[i], &(allRanks[i][0]), allNProcRows[i], allNProcRows[i], allNProcCols[i] );
    }
    for(int i = 0; i < allRanks.size(); i++)
    {
      int tempRow;
      int tempCol;
      //Cblacs_pcoord(BLACSContexts[i], this->myId, &tempRow, &tempCol);
      int tempNProcRows;
      int tempNProcCols;
      Cblacs_gridinfo(BLACSContexts[i], &tempNProcRows, &tempNProcCols, &tempRow, &tempCol );
      
      if(i == this->group)
      {
        this->myRow = tempRow;
        this->myCol = tempCol;
        
        this->nProcRows = tempNProcRows;
        this->nProcCols = tempNProcCols;
        
        if( tempNProcRows != allNProcRows[i] || tempNProcCols != allNProcCols[i] )
        {
          error("Error: An internal error happened while initializing a grouped BLACS grid.");
        }
      }
    }
    this->context = BLACSContexts[this->group];
    
    this->nDefaultBlockRows = options.defaultBlockSize;
    this->nDefaultBlockCols = options.defaultBlockSize;
  }
  if(this->mpiRoot == true && (this->myRow != 0 || this->myCol != 0))
  {
    error("Error: An internal error happened after initializing BLACS grid.");
  }
}

Communicator::~Communicator()
{
  Cblacs_gridexit(this->context);
  if(this->type == globalCommunicator)
  {
    MPI_Finalize();
  }
}

std::string Communicator::creationMessage()
{
  std::stringstream messageBuffer;
  messageBuffer << this->nProcRows << "x" << this->nProcCols << " BLACS grid successfully created." << std::endl;
  messageBuffer << "Using " << this->nProcRows*this->nProcCols << " MPI process" << ((this->nProcRows*this->nProcCols>1)?"es":"") << "." << std::endl;
  messageBuffer << "BLACS block size: "<< this->nDefaultBlockRows << "x" << this->nDefaultBlockCols << std::endl;
  return messageBuffer.str();
}

void Communicator::error(std::string e)
{
  std::cerr << "**********************************" << std::endl;
  std::cerr << e << std::endl;
  std::cerr << "**********************************" << std::endl;
  MPI_Abort(MPI_COMM_WORLD, -1);
}

void Communicator::setMPIRoot()
{
  this->mpiRoot = this->mpiRank == 0;
}

void Communicator::barrier()
{
  Cblacs_barrier(communicator->context, "All");
}

void Communicator::broadcast(int * values, int nValues)
{
  MPI_Bcast(values, nValues, MPI_INT, 0, this->mpiCommunicator);
}

void Communicator::broadcast(double * values, int nValues)
{
  MPI_Bcast(values, nValues, MPI_DOUBLE, 0, this->mpiCommunicator);
}

void Communicator::broadcast(char * values, int nValues)
{
  MPI_Bcast(values, nValues, MPI_CHAR, 0, this->mpiCommunicator);
}

void Communicator::broadcast(std::string & s)
{
  int size = s.size();
  broadcast(&size);
  
  char * cs = new char [size + 1];
  if(this->mpiRoot)
  {
    std::strcpy(cs, s.c_str());
  }
  
  broadcast(cs, size + 1);
  
  s = std::string(cs, size);
  
  delete [] cs;
}

void Communicator::broadcast(bool * value)
{
  int sendValue;
  if(this->mpiRoot)
  {
    sendValue = ((*value==true)?1:0);
  }
  MPI_Bcast(&sendValue, 1, MPI_INT, 0, this->mpiCommunicator);
  *value = ((sendValue==1)?true:false);
}

void Communicator::broadcast(std::vector<int> & v)
{
  int size = 0;
  if(this->mpiRoot)
  {
    size = v.size();
  }
  broadcast(&size);
  if(this->mpiRoot == false)
  {
    v.clear();
    v = std::vector<int>(size, 0);
  }
  broadcast(&(v[0]), size);
}

int * Communicator::gather(int * values, int nValues)
{
  int success;
  int * result = NULL;
  if(this->mpiRoot)
  {
    result = new int [this->mpiNumTasks*nValues];
  }
  success = MPI_Gather( values, nValues, MPI_INT, result, nValues, MPI_INT, 0, this->mpiCommunicator);
  if(success != MPI_SUCCESS)
  {
    misc.error("Error: Error on MPI gather communication. Terminating.", 0);
  }
  return result;
}

double * Communicator::gather(double * values, int nValues)
{
  int success;
  double * result = NULL;
  if(this->mpiRoot)
  {
    result = new double [this->mpiNumTasks*nValues];
  }
  success = MPI_Gather( values, nValues, MPI_DOUBLE, result, nValues, MPI_DOUBLE, 0, this->mpiCommunicator);
  if(success != MPI_SUCCESS)
  {
    misc.error("Error: Error on MPI gather communication. Terminating.", 0);
  }
  return result;
}

char * Communicator::gather(char * values, int nValues)
{
  int success;
  char * result = NULL;
  if(this->mpiRoot)
  {
    result = new char [this->mpiNumTasks*nValues];
  }
  success = MPI_Gather( values, nValues, MPI_CHAR, result, nValues, MPI_DOUBLE, 0, this->mpiCommunicator);
  if(success != MPI_SUCCESS)
  {
    misc.error("Error: Error on MPI gather communication. Terminating.", 0);
  }
  return result;
}

int * Communicator::asymmetricGather(int * values, int nValues, int * totalSize)
{
  int success;
  int * sizes = NULL;
  int * displacements = NULL;
  
  if(this->mpiRoot)
  {
    displacements = new int [this->mpiNumTasks];
  }
  
  //Gather the sizes of each array to gather. Only on root. On ohter processes sizes == NULL.
  sizes = gather(&nValues, 1);

  //Compute the total resulting size and displacements. Only in root.
  *totalSize = 0;  
  if(this->mpiRoot)
  {
    for(int i = 0; i<this->mpiNumTasks; i++)
    {
      displacements[i] = *totalSize;
      *totalSize += sizes[i];
    }
  }
  
  int * result = NULL;
  if(this->mpiRoot)
  {
    result = new int [*totalSize];
  }

  success = MPI_Gatherv(values, nValues, MPI_INT, result, sizes, displacements, MPI_INT, 0, this->mpiCommunicator);
  if(success != MPI_SUCCESS)
  {
    //char error_string[BUFSIZ];
    //int length_of_error_string, error_class;
    
    //MPI_Error_class(success, &error_class);
    //MPI_Error_string(error_class, error_string, &length_of_error_string);
    std::cout << success << std::endl;
    misc.error("Error: Error on MPI gather communication. Terminating.", 0);
  }
  
  if(this->mpiRoot)
  {
    delete [] sizes;
    delete [] displacements;
  }
  
  return result;
}

double * Communicator::asymmetricGather(double * values, int nValues, int * totalSize)
{
  int success;
  int * sizes = NULL;
  int * displacements = NULL;
  
  if(this->mpiRoot)
  {
    displacements = new int [this->mpiNumTasks];
  }
  
  //Gather the sizes of each array to gather. Only on root. On ohter processes sizes == NULL.
  sizes = gather(&nValues, 1);
  
  //Compute the total resulting size and displacements. Only in root.
  *totalSize = 0;  
  if(this->mpiRoot)
  {
    for(int i = 0; i<this->mpiNumTasks; i++)
    {
      displacements[i] = *totalSize;
      *totalSize += sizes[i];
    }
  }
  
  double * result = NULL;
  if(this->mpiRoot)
  {
    result = new double [*totalSize];
  }
  
  success = MPI_Gatherv(values, nValues, MPI_DOUBLE, result, sizes, displacements, MPI_DOUBLE, 0, this->mpiCommunicator);
  if(success != MPI_SUCCESS)
  {
    //char error_string[BUFSIZ];
    //int length_of_error_string, error_class;
    
    //MPI_Error_class(success, &error_class);
    //MPI_Error_string(error_class, error_string, &length_of_error_string);
    std::cout << success << std::endl;
    misc.error("Error: Error on MPI gather communication. Terminating.", 0);
  }
  
  if(this->mpiRoot)
  {
    delete [] sizes;
    delete [] displacements;
  }
  
  return result;
}

char * Communicator::asymmetricGather(const char * values, int nValues, int * totalSize)
{
  int success;
  int * sizes = NULL;
  int * displacements = NULL;
  
  if(this->mpiRoot)
  {
    displacements = new int [this->mpiNumTasks];
  }
  
  //Gather the sizes of each array to gather. Only on root. On ohter processes sizes == NULL.
  sizes = gather(&nValues, 1);
  
  //Compute the total resulting size and displacements. Only in root.
  *totalSize = 0;  
  if(this->mpiRoot)
  {
    for(int i = 0; i<this->mpiNumTasks; i++)
    {
      displacements[i] = *totalSize;
      *totalSize += sizes[i];
    }
  }
  
  char * result = NULL;
  if(this->mpiRoot)
  {
    result = new char [*totalSize];
  }
  
  success = MPI_Gatherv(values, nValues, MPI_CHAR, result, sizes, displacements, MPI_CHAR, 0, this->mpiCommunicator);
  if(success != MPI_SUCCESS)
  {
    //char error_string[BUFSIZ];
    //int length_of_error_string, error_class;
    
    //MPI_Error_class(success, &error_class);
    //MPI_Error_string(error_class, error_string, &length_of_error_string);
    std::cout << success << std::endl;
    misc.error("Error: Error on MPI gather communication. Terminating.", 0);
  }
  
  if(this->mpiRoot)
  {
    delete [] sizes;
    delete [] displacements;
  }
  
  return result;
}

std::string Communicator::asymmetricGather(std::string & values)
{
  int totalSize;
  char * charresult = communicator->asymmetricGather(values.data(), values.size(), &totalSize);
  
  if(this->mpiRoot)
  {
    std::string result(charresult, totalSize);
    delete charresult;
    return result;
  }
  else
  {
    return std::string();
  }
}

void Communicator::storeArraysMPI(std::string fn, std::string data)
{
  int rc;
  MPI_Status status;
  MPI_File fh;
  
  if(this->mpiCommunicator != MPI_COMM_WORLD  || this->type != globalCommunicator)
  {
    misc.error("Error: An internal error has happened. storeArrayMPI can only be used with global communicators.", 0);
  }
  
  rc = MPI_File_open(this->mpiCommunicator, fn.c_str(), MPI_MODE_EXCL | MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  if(rc != MPI_SUCCESS)
  {
    misc.error("Error: I have been unable to open the file [ " + fn + " ] for writting. Is it possible that the file already exists?", 0);
  }

  long int maxIntTest = std::numeric_limits<int>::max();
  long long maxLongTest = std::numeric_limits<long>::max();
  if(data.size() > maxIntTest)
  {
    misc.error("Error: An internal error has happened. MPI data size for storage is larger than expected for a process. Please, increase the numnber of processes, or contact us.", 0);
  }
  int size = data.size();
  
  int * sizes = gather(&size, 1);
  long * offsets = new long [this->mpiNumTasks];
  
  if( this->mpiRoot == true )
  {
    offsets[0] = 0;
    for(int i = 1; i < this->mpiNumTasks; i++)
    {
      offsets[i] = offsets[ i - 1 ] + sizes[ i - 1 ];
      
      long long t1 = offsets[ i - 1 ];
      long long t2 = sizes[ i - 1 ];
      if(t1 + t2 > maxLongTest)
      {
        misc.error("Error: An internal error has happened. MPI data size for storage is larger than expected for a process. Please, contact us.", 0);
      }
    }
  }
  MPI_Bcast(offsets, this->mpiNumTasks, MPI_LONG, 0, this->mpiCommunicator);
  
  rc = MPI_File_write_at_all(fh, offsets[this->mpiRank], data.data(), size, MPI_CHAR, &status);
  if(rc != MPI_SUCCESS)
  {
    misc.error("Error: I have been unable to write on file [ " + fn + " ] for writting. Is there enough disk space?", 0);
  }
  
  rc = MPI_File_close(&fh);
  if(rc != MPI_SUCCESS)
  {
    misc.error("Error: I have been unable to close the file [ " + fn + " ].", 0);
  }

  
  delete [] offsets;
  if(this->mpiRoot == true)
  {
    delete [] sizes;
  }
}

void Communicator::mpiDebug(int pid)
{
  //Code for attaching gdb to a particular process
  if((communicator->mpiRank == pid || pid < 0) /*&& options.mpiDebug*/){
   int i = 0;
   std::cout << "PID " << communicator->mpiRank << " on " << communicator->hostName << " ready for attach\n" << std::endl;
   fflush(stdout);
   while (0 == i)
   {
    }
  }
}
