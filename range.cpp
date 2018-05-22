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

#include "range.h"
#include "misc.h"

#include <string>
#include <sstream>

Range::Range()
{
  this->downThresholdActive = false;
  this->upThresholdActive = false;
};

Range::Range(bool flag, int upThreshold)
{
  if(flag)
  {
    misc.error("Error: An internal error was happened. Invalid flag when defining range.", 0);
  }
  this->downThresholdActive = false;
  this->upThresholdActive = true;
  this->iUpThreshold = upThreshold;
  
};

Range::Range(int downThreshold, bool flag)
{
  if(flag)
  {
    misc.error("Error: An internal error was happened. Invalid flag when defining range.", 0);
  }
  this->downThresholdActive = true;
  this->upThresholdActive = false;
  this->iDownThreshold = downThreshold;
};

Range::Range(int downThreshold, int upThreshold)
{
  this->downThresholdActive = true;
  this->upThresholdActive = true;
  this->iDownThreshold = downThreshold;
  this->iUpThreshold = upThreshold;
};

Range::Range(bool flag, double upThreshold)
{
  if(flag)
  {
    misc.error("Error: An internal error was happened. Invalid flag when defining range.", 0);
  }
  this->downThresholdActive = false;
  this->upThresholdActive = true;
  this->dUpThreshold = upThreshold;
  
};

Range::Range(double downThreshold, bool flag)
{
  if(flag)
  {
    misc.error("Error: An internal error was happened. Invalid flag when defining range.", 0);
  }
  this->downThresholdActive = true;
  this->upThresholdActive = false;
  this->dDownThreshold = downThreshold;
};

Range::Range(double downThreshold, double upThreshold)
{
  this->downThresholdActive = true;
  this->upThresholdActive = true;
  this->dDownThreshold = downThreshold;
  this->dUpThreshold = upThreshold;
};

Range::~Range() {};

bool Range::checkRange(int x)
{
  if( (this->downThresholdActive == false) && (this->upThresholdActive == false) )
  {
    return true;
  }
  else if( (this->downThresholdActive == true) && (this->upThresholdActive == false) )
  {
    if(x >= this->iDownThreshold)
    {
      return true;
    }
    return false;
  }
  else if( (this->downThresholdActive == false) && (this->upThresholdActive == true) )
  {
    if(x <= this->iUpThreshold)
    {
      return true;
    }
    return false;
  }
  else
  {
    if( (x >= this->iDownThreshold) && (x <= this->iUpThreshold) )
    {
      return true;
    }
    return false;
  }
  return false;
}

bool Range::checkRange(double x)
{
  if( (this->downThresholdActive == false) && (this->upThresholdActive == false) )
  {
    return true;
  }
  else if( (this->downThresholdActive == true) && (this->upThresholdActive == false) )
  {
    if(x >= this->dDownThreshold)
    {
      return true;
    }
    return false;
  }
  else if( (this->downThresholdActive == false) && (this->upThresholdActive == true) )
  {
    if(x <= this->dUpThreshold)
    {
      return true;
    }
    return false;
  }
  else
  {
    if( (x >= this->dDownThreshold) && (x <= this->dUpThreshold) )
    {
      return true;
    }
    return false;
  }
  return false;
}

std::string Range::explainRange(int x)
{
  std::ostringstream ss;
  
  if( (this->downThresholdActive == false) && (this->upThresholdActive == false) )
  {
    return "";
  }
  else if( (this->downThresholdActive == true) && (this->upThresholdActive == false) )
  {
    ss << "It has to be >= " << this->iDownThreshold;
    return ss.str();
  }
  else if( (this->downThresholdActive == false) && (this->upThresholdActive == true) )
  {
    ss << "It has to be <= " << this->iUpThreshold;
    return ss.str();
  }
  else
  {
    ss << "It has to be >= " << this->iDownThreshold << " and <= " << this->iUpThreshold;
    return ss.str();
  }
  return ""; 
}

std::string Range::explainRange(double x)
{
  std::ostringstream ss;
  
  if( (this->downThresholdActive == false) && (this->upThresholdActive == false) )
  {
    return "";
  }
  else if( (this->downThresholdActive == true) && (this->upThresholdActive == false) )
  {
    ss << "It has to be >= " << this->dDownThreshold;
    return ss.str();
  }
  else if( (this->downThresholdActive == false) && (this->upThresholdActive == true) )
  {
    ss << "It has to be <= " << this->dUpThreshold;
    return ss.str();
  }
  else
  {
    ss << "It has to be >= " << this->dDownThreshold << " and <= " << this->dUpThreshold;
    return ss.str();
  }
  return ""; 
}
