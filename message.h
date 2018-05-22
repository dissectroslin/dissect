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

#ifndef MESSAGE_H
#define MESSAGE_H

#include "communicator.h"
#include "global.h"

#include <string>
#include <fstream>
#include <iomanip>

#if defined(BOOSTLIB) && defined(ZLIB)
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#endif

class Message
{
public:
  bool active;
  
  std::ostream * output;
  bool outputIsFile;
  
  bool logActive;
  std::ostream * log;
  
  bool newLine;
  std::string tab;
  
  std::ostream * compressedFile;
  boost::iostreams::filtering_streambuf<boost::iostreams::output> * gzout;
  
  Message();
  Message(std::ostream * out);
  Message(std::ostream * out, std::ostream * logstream);
  Message(std::string fname);
  
  ~Message();

  void redirect(std::string fname);
  void setWidth(int width);
  void flush();
  
  template<typename T>
  Message& operator<< (T x)
  {
    int width = this->output->width();
    this->output->width(0);
    std::string tab = this->newLine?this->tab:"";
    
    if(communicator->mpiRoot && this->active == true)
    {
      this->newLine = false;
      *(this->output) << tab << std::setw(width) << x;
      if(this->logActive)
      {
        this->log->width(0);
        *(this->log) << tab << std::setw(width) << x;
      }
    }
    return *this;
  }
  
 
  Message& operator<<(std::ostream& (*fn)(std::ostream&));
  Message& operator<<(std::_Setw manip);
  Message& operator<<(std::_Setprecision manip);

};

#endif
