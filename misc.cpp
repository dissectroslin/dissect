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

#include "misc.h"
#include "global.h"
#include "options.h"
#include "communicator.h"
#include "message.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>

#include <mpi.h>

Misc::Misc()
{
  this->logFile = NULL;
  
  this->localMemory = 0.;
  this->currentMemory = 0.;
  this->maxMemory = 0.;
}

Misc::~Misc()
{
  if(this->logFile != NULL)
  {
    if(this->logFile->is_open())
    {
      this->logFile->close();
    }
    delete this->logFile;
    this->logFile = NULL;
  }
}

void Misc::changeOutputs(std::ostream & out)
{
  this->message = Message(&out);
  if(this->logFile != NULL)
  {
    if(this->logFile->is_open())
    {
      this->logFile->close();
    }
    delete this->logFile;
    this->logFile = NULL;
  }
}

void Misc::changeOutputs(std::ostream & out, std::string fname)
{
  if(this->logFile != NULL)
  {
    if(this->logFile->is_open())
    {
      this->logFile->close();
    }
    delete this->logFile;
    this->logFile = NULL;
  }
  this->logFile = new std::ofstream();
  openOutputFileInRoot(*this->logFile, fname);
  this->message = Message(&out, this->logFile);
}

void Misc::checkFileExists(std::string f)
{
  if(doesFileExists(f) == false)
  {
    std::string e = "Error: The file [ " + f + " ] does not exist.";
    error(e, 0);
  }
}

bool Misc::doesFileExists(std::string f)
{
  std::ifstream test;
 
  test.open(f.c_str(), std::ifstream::in);
  if(test.fail())
  {
    test.clear(std::ios::failbit);
    test.close();
    return false;
  }
  test.close();
  return true;
}

void Misc::checkFileNoExists(std::string f)
{
  std::ifstream test;
  
  test.open(f.c_str(), std::ifstream::in);
  if(!test.fail())
  {
    test.clear(std::ios::failbit);
    test.close();
    std::string e = "Error: The file [ " + f + " ] exists. We detected problems in this cases for some MPI implementations. Please, delete it manually before running the analysis.";
    error(e, 0);
  }
  test.close();
}

void Misc::openOutputFileInRoot(std::ofstream &file, std::string fname)
{
  if(communicator->mpiRoot)
  {
    file.open(fname.c_str(), std::ofstream::out);
    if(!file.is_open())
    {
      error("Error: The file [ " + fname +  " ] can not be opened for writting.", 0);
    }
  }
}

void Misc::openBinaryOutputFileInRoot(std::ofstream &file, std::string fname)
{
  if(communicator->mpiRoot)
  {
    file.open(fname.c_str(), std::ofstream::out | std::ios_base::binary);
    if(!file.is_open())
    {
      error("Error: The file [ " + fname +  " ] can not be opened for writting.", 0);
    }
  }
}

void Misc::error(std::string e, int type)
{
  //this->message.flush();
  std::cerr << "\n**********************************\n";
  std::cerr << "Process " << communicator->mpiRank << " failed.\n";
  std::cerr << e << std::endl;
  std::cerr << "**********************************\n" << std::endl;
  std::cerr.flush();
  
  std::cout << "\n**********************************\n";
  std::cout << "Process " << communicator->mpiRank << " failed.\n";
  std::cout << e << std::endl;
  std::cout << "**********************************\n" << std::endl;
  std::cout.flush();
  
  #ifdef CREATETRACEBACK
  abort();//Just For creating backtrace
  #endif
  MPI_Abort(MPI_COMM_WORLD, -1);
}

void Misc::write(std::string s, int type)
{
  std::cerr << s << std::endl;
}

bool Misc::checkFileHeader(std::string header, std::string templateHeader)
{
  std::istringstream ssTestHeader(header);
  std::istringstream ssTemplateHeader(templateHeader);
  
  std::string templateElement;
  std::string testElement;
  while( !(ssTemplateHeader >> templateElement).fail() )
  {
    if( (ssTestHeader >> testElement).fail() )
    {
      return false;
    }
    if( testElement != templateElement )
    {
      return false;
    }
  }
  if( !(ssTestHeader >> testElement).fail() )
  {
    return false;
  }
  
  return true;
}

bool Misc::gt(bool test)
{
  bool result = test;
  communicator->broadcast(&result);
  return result;
}

std::string Misc::setGetElapsedTime(std::string name, bool clear)
{
  double result;
  std::stringstream ssresult;
  
  time_t now = time(NULL);
  
  if(this->elapsedTime.count(name) == 1)
  {
    result = difftime(now, this->elapsedTime[name]);
  }
  else
  {
    result = 0.;
  }
  if(clear == true)
  {
    if(this->elapsedTime.count(name) == 1)
    {
      this->elapsedTime.erase(name);
    }
  }
  else
  {
    this->elapsedTime[name] = now;
  }
  
  if(result < 60.)
  {
    ssresult << result << "s";
    return ssresult.str();
  }
  if(result > 60. && result < 3600.)
  {
    ssresult << std::setprecision(4) << result/60. << "m";
    return ssresult.str();
  }
  else
  {
    ssresult << std::setprecision(4) << result/3600. << "h";
    return ssresult.str();
  }
}

void Misc::estimateMaxMemory(double changeInMemory)
{
  /*this->localMemory += changeInMemory;
  if(this->localMemory < 0.)
  {
    this->localMemory = 0.;
  }
  
  double * arrayMemoryThreads = communicator->gather(&this->localMemory, 1);
  
  double newMaxMemory = 0.;
  if(communicator->mpiRoot)
  {
    for(int i = 0; i<communicator->mpiNumTasks; i++)
    {
      newMaxMemory += arrayMemoryThreads[i];
    }
  }
   
  if(newMaxMemory > this->maxMemory)
  {
    this->maxMemory = newMaxMemory;
  }
  
  delete [] arrayMemoryThreads;*/
  
  this->currentMemory += changeInMemory;
  if(this->currentMemory < 0.)
  {
    this->currentMemory = 0.;
  }
  
  if(this->currentMemory > this->maxMemory)
  {
    this->maxMemory = this->currentMemory;
  }
}

