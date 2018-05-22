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

#ifndef MISC_H
#define MISC_H

#include "communicator.h"
#include "global.h"
#include "message.h"

#include <string>
#include <fstream>
#include <iomanip>
#include <map>
#include <ctime>

class Misc
{
public:

  Message message;
  
  std::ofstream * logFile;
  
  std::map<std::string, time_t> elapsedTime;
  
  double localMemory;
  double currentMemory;
  double maxMemory;
  
  Misc();
  ~Misc();
  
  void changeOutputs(std::ostream & out);
  void changeOutputs(std::ostream & out, std::string fname);
  void checkFileExists(std::string f);
  bool doesFileExists(std::string f);
  void checkFileNoExists(std::string f);
  void openOutputFileInRoot(std::ofstream &file, std::string fname);
  void openBinaryOutputFileInRoot(std::ofstream &file, std::string fname);
  void error(std::string e, int type);
  void write(std::string s, int type);
  
  bool checkFileHeader(std::string header, std::string templateHeader);
  
  /**
   * Evaluate an boolean expression in root and spreads the result in all processes
   * 
   * Evaluate an boolean expression in root and spreads the result in all nodes. The result is returned by the function in each process
   * 
   * \param test expression to test on root
   * \return test result scattered on all nodes
   */
  bool gt(bool test);
  
  std::string setGetElapsedTime(std::string name, bool clear = false);
  void estimateMaxMemory(double changeInMemory);
};

#endif
