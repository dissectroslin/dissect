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

#include "results.h"
#include "reml.h"
#include "covariancematrix.h"
#include "global.h"
#include "options.h"
#include "message.h"
#include "auxiliar.h"

#include <iomanip>

Results::Results(REML * reml)
{
  this->logLikelihood = reml->logLikelihood;
  this->nVariances = reml->V->variances.size();
}

void Results::compareREMLs(REML * reml)
{
  double LogRatio = 2.0*(this->logLikelihood - reml->logLikelihood);
  if(LogRatio < 0.)
  {
    LogRatio = 0.;
  }
  int df = this->nVariances - reml->V->variances.size();
  
  Message message(options.outFile + ".reml.fixed");
  
  message << "LRT\t" << std::setprecision(3) << LogRatio << std::endl;
  message << "df\t" << std::setprecision(1) << df << std::endl;
  message << "Pval\t" << std::setprecision(4) << 0.5*chi1_CDF(df, LogRatio) << std::endl;
}

