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

#include "message.h"
#include "communicator.h"
#include "global.h"
#include "misc.h"
#include "options.h"

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

#if defined(BOOSTLIB) && defined(ZLIB)
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#endif

Message::Message()
{
  this->active = true;
  
  this->output = &std::cout;
  this->log = NULL;
  this->logActive = false;
  this->outputIsFile = false;
  
  this->compressedFile = NULL;
  
  this->tab = "";
  this->newLine = false;
}

Message::Message(std::ostream * out)
{
  this->active = true;
  
  this->output = out;
  this->log = NULL;
  this->logActive = false;
  this->outputIsFile = false;
  
  this->compressedFile = NULL;
  
  this->tab = "";
  this->newLine = false;
}

Message::Message(std::ostream * out, std::ostream * logstream)
{
  this->active = true;
  
  this->output = out;
  this->log = logstream;
  this->logActive = true;
  this->outputIsFile = false;
  
  this->compressedFile = NULL;
  
  this->tab = "";
  this->newLine = false;
}

Message::Message(std::string fname)
{
  this->active = true;
  
  std::ofstream * file = new std::ofstream();
  if(options.compressedOutput == false)
  {
    misc.openOutputFileInRoot(*file, fname);
    this->output = file;
    this->compressedFile = NULL;
  }
  else
  {
#if defined(BOOSTLIB) && defined(ZLIB)
    
    misc.openBinaryOutputFileInRoot(*file, fname + ".gz");
    boost::iostreams::filtering_streambuf<boost::iostreams::output> * gzout = new boost::iostreams::filtering_streambuf<boost::iostreams::output>();
    gzout->push(boost::iostreams::gzip_compressor());
    gzout->push(*file);
    
    this->output = new std::ostream(gzout);
    this->compressedFile = file;
    
    misc.error("Error: Sorry, this function has to be further tested before use.", 0); //It has to be tested what happens with gzout, which has to be freed at the end. It also has to be tested that this->output has only to be created at root, otherwise it will not be properly freed on all processes.
#else
  misc.error("Error: The current version of DISSECT is compiled without zlib and boost support. These libraries are required for writting compressed output. Please, recompile it using these libraries or ask for support.", 0);
#endif    
  }
  
  this->log = NULL;
  this->logActive = false;
  this->outputIsFile = true;
  
  this->tab = "";
  this->newLine = false;
}

Message::~Message()
{
  if(this->outputIsFile == true)
  {
    if(communicator->mpiRoot)
    {
      if(this->compressedFile != NULL)
      {
        ((std::ofstream*)this->compressedFile)->close();
      }
      ((std::ofstream*)this->output)->close();
      delete this->output;
      if(this->compressedFile != NULL)
      {
        delete this->compressedFile;
      }
    }
  }
}

void Message::redirect(std::string fname)
{
  if(this->outputIsFile == true)
  {
    if(communicator->mpiRoot)
    {
      ((std::ofstream*)this->output)->close();
    }
    delete this->output;
    std::ofstream * file = new std::ofstream();
    misc.openOutputFileInRoot(*file, fname);
    this->output = file;
  }
  else
  {
    misc.error("Error: An internal error was happened. This message can not be reconverted to a file type.", 0);
  }
}

void Message::setWidth(int width)
{
  if(communicator->mpiRoot && this->active == true)
  {
    this->output->width(width);
    if(this->logActive == true)
    {
      this->log->width(width);
    }
  }
}

void Message::flush()
{
  if(communicator->mpiRoot && this->active == true)
  {
    this->output->flush();
    if(this->logActive)
    {
      this->log->flush();
    }
  }
}

Message& Message::operator<<(std::ostream& (*fn)(std::ostream&))
{
  if(communicator->mpiRoot && this->active == true)
  {
    this->newLine = true;
    
    *(this->output) << fn;
    if(this->logActive)
    {
      *(this->log) << fn;
    }
  }
  return *this;
}

Message& Message::operator<<(std::_Setw manip)
{
  if(communicator->mpiRoot && this->active == true)
  {
    *(this->output) << manip;
    if(this->logActive)
    {
      *(this->log) << manip;
    }
  }
  return *this;
}

Message& Message::operator<<(std::_Setprecision manip)
{
  if(communicator->mpiRoot && this->active == true)
  {
    *(this->output) << manip;
    if(this->logActive)
    {
      *(this->log) << manip;
    }
  }
  return *this;
}
