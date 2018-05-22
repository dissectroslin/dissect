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

#include "memusage.h"
#include "misc.h"

#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

MemUsage::MemUsage()
{
  this->lastResidentSet = 0.;
  this->lastVM = 0.;
}

void MemUsage::memUsage(double& vm_usage, double& resident_set)
{
  using std::ios_base;
  using std::ifstream;
  using std::string;
  
  vm_usage     = 0.0;
  resident_set = 0.0;
  
  // 'file' stat seems to give the most reliable results
  ifstream stat_stream("/proc/self/stat",ios_base::in);
  
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;
  
  // the two relevant fields for us
  unsigned long vsize;
  long rss;
  
  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
  >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
  >> utime >> stime >> cutime >> cstime >> priority >> nice
  >> O >> itrealvalue >> starttime >> vsize >> rss;
  
  stat_stream.close();
  
  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  vm_usage     = vsize / 1024.0;
  resident_set = rss * page_size_kb;
}

void MemUsage::showMemUsage(std::string title)
{
  double vm, rss;
  memUsage(vm, rss);
  
  double dvm = vm - this->lastVM;
  double drss = rss - this->lastResidentSet;
  
  this->lastResidentSet = rss;
  this->lastVM = vm;
  
  misc.message << "**  Memory usage  *******************************" << std::endl;
  misc.message << "-- " + title + " --" << std::endl;
  misc.message << "Current: " << vm/1024. << "MB " << rss/1024. << "MB" << std::endl;
  misc.message << "Delta: " << dvm/1024. << "MB " << drss/1024. << "MB" << std::endl;
  misc.message << "*************************************************" << std::endl;
  
}
