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

#include <vector>
#include <map>

#include "groupeffects.h"
#include "labeledmatrix.h"
#include "matrix.h"
#include "misc.h"
#include "communicator.h"
#include "options.h"
#include "auxiliar.h"


GroupEffects::GroupEffects(std::string fn)
{
  this->effects = new LabeledMatrix(fn);
}

GroupEffects::GroupEffects(std::vector<std::string> fns, RowColumn rowcolumn)
{
  if(fns.size() < 1)
  {
    misc.error("Error: An internal error has happened. There file list is empty when loading effects.", 0);
  }
  this->effects = new LabeledMatrix(fns[0]);
  for(int i = 1; i<fns.size(); i++)
  {
    LabeledMatrix * lmt = new LabeledMatrix(fns[i]);
    if(rowcolumn == column)
    {
      this->effects->appendVertically(lmt);
    }
    else
    {
      this->effects->appendHorizontally(lmt);
    }
    delete lmt;
  }
}

GroupEffects::~GroupEffects()
{
  if(this->effects != NULL)
  {
    this->effects = NULL;
    delete this->effects;
  }
}

LabeledMatrix * GroupEffects::computeCorrelations(RowColumn rowcolumn)
{
  LabeledMatrix * correlations = new LabeledMatrix();
  
  Matrix * standardizedEffects = new Matrix( this->effects->getMatrix() );
  standardizedEffects->standardizeMatrix(rowcolumn);
  
  if(rowcolumn == column)
  {
    correlations->getMatrix()->multiply(standardizedEffects, 'T', standardizedEffects, 'N', 1./double(standardizedEffects->nGlobRows));
    
    correlations->setRowLabels( this->effects->getColLabels() );
    correlations->setColLabels( this->effects->getColLabels() );
  }
  else
  {
    correlations->getMatrix()->multiply(standardizedEffects, 'N', standardizedEffects, 'T', 1./double(standardizedEffects->nGlobCols));
    
    correlations->setRowLabels( this->effects->getRowLabels() );
    correlations->setColLabels( this->effects->getRowLabels() );
  }
  
  correlations->getMatrix()->symmetrizeTriangularMatrix();
  
  delete standardizedEffects;

  return correlations;
}

LabeledMatrix * GroupEffects::computeCovariances(RowColumn rowcolumn)
{
  LabeledMatrix * covariances = new LabeledMatrix();
  
  Matrix * centeredEffects = new Matrix( this->effects->getMatrix() );
  centeredEffects->centerMatrixRowsColumns(rowcolumn);
  
  if(rowcolumn == column)
  {
    covariances->getMatrix()->multiply(centeredEffects, 'T', centeredEffects, 'N', 1./double(centeredEffects->nGlobRows));
    
    covariances->setRowLabels( this->effects->getColLabels() );
    covariances->setColLabels( this->effects->getColLabels() );
  }
  else
  {
    covariances->getMatrix()->multiply(centeredEffects, 'N', centeredEffects, 'T', 1./double(centeredEffects->nGlobCols));
    
    covariances->setRowLabels( this->effects->getRowLabels() );
    covariances->setColLabels( this->effects->getRowLabels() );
  }
  
  covariances->getMatrix()->symmetrizeTriangularMatrix();
  
  delete centeredEffects;

  return covariances;
}

void GroupEffects::getPairsHighlyCorrelated(LabeledMatrix * correlations, double threshold, std::vector<int> & rowIdxs, std::vector<int> & colIdxs)
{
  correlations->getMatrix()->getGlobalIndexOutsideRange(-threshold, threshold, rowIdxs, colIdxs);
}

std::map<std::string, GroupAttributes> GroupEffects::getGroupPositions(std::string fn)
{
  std::vector< std::vector<std::string> > table;
  getTableFromFile(fn, table, 4);
  
  std::map<std::string, GroupAttributes> results;
  
  for( int i = 0; i<table.size(); i++ )
  {
    GroupAttributes grow;
    
    grow.name = table[i][0];
    grow.chrom = table[i][1];
    grow.minpos = string2Number<long>(table[i][2], "Error: There is a value in the third row of file [ " + fn + " ] which is not a number.");
    grow.maxpos = string2Number<long>(table[i][3], "Error: There is a value in the fourth row of file [ " + fn + " ] which is not a number.");
    
    if( results.count(grow.name) != 0 )
    {
      misc.error("Error: The identificator " + grow.name + " appears at least two times in the file [ " + fn + " ]", 0);
    }
    
    results[grow.name] = grow;
  }
  
  misc.message << "Positions for " << results.size() << " groups has been loaded from [ " << fn << " ]." << std::endl;
  
  return results;
}

void GroupEffects::filterCorrelatedGroups(RowColumn rowcolumn, double threshold, std::string fn)
{
  std::vector<int> rowIdxs;
  std::vector<int> colIdxs;
  
  misc.message << "Filtering correlated groups..." << std::endl;
  LabeledMatrix * correlations = computeCorrelations(rowcolumn);
  getPairsHighlyCorrelated(correlations, threshold, rowIdxs, colIdxs);
  misc.message << "There are " << rowIdxs.size() << " group pairs which have a correlation larger than " << threshold << std::endl;
  
  std::vector<std::string> labelsToKeep;
  if( communicator->mpiRoot == true )
  {
    std::map<std::string, GroupAttributes> positions;
    if( fn != "" )
    {
      positions = getGroupPositions(fn);
    }
    
    if( rowIdxs.size() != colIdxs.size() )
    {
      misc.error("Error: An internal error was happened. Mismatching number of rows and columns when filtering groups by their correlation.", 0);
    }
    
    std::vector<std::string> labels = correlations->getRowLabels();
    
    std::vector<int> filteredRowIdxs;
    std::vector<int> filteredColIdxs;
    std::map<int, int> indexFrequency;
    
    //From the indexs of related labels, remove those which are far enough and store the number of times they appear.
    for( int i = 0; i < rowIdxs.size(); i++ )
    {
      if( rowIdxs[i] == colIdxs[i] )
      {
        continue;
      }
      std::string group1 = labels[ rowIdxs[i] ];
      std::string group2 = labels[ colIdxs[i] ];
      
      bool overlap = false;
      if( fn != "" )
      {
        if( positions.count( group1 ) == 0 || positions.count( group2 ) == 0 )
        {
          misc.error("Error: The set " + group1 + " and/or the set " + group2 + " do not have an entry on the file [" + fn + "] defining their position.", 0);
        }
        std::pair<bool, double> distance = positions[ group1 ].getDistance(positions[ group2 ]);
        if( distance.first == false )
        {
          continue;
        }
        if( distance.first == true && distance.second > options.groupDistanceForDiscarding )
        {
          continue;
        }
      }
      
      filteredRowIdxs.push_back( rowIdxs[i] );
      filteredColIdxs.push_back( colIdxs[i] );
      if( indexFrequency.count( rowIdxs[i] ) == 0 )
      {
        indexFrequency[ rowIdxs[i] ] = 0;
      }
      if( indexFrequency.count( colIdxs[i] ) == 0 )
      {
        indexFrequency[ colIdxs[i] ] = 0;
      }
      indexFrequency[ rowIdxs[i] ] += 1;
      indexFrequency[ colIdxs[i] ] += 1;
    }
    if( fn != "" )
    {
      misc.message << "From " << rowIdxs.size() << " group pairs, only " << filteredRowIdxs.size() << " are close enough to be filtered based on the distances defined in [ " + fn +  " ]." << std::endl;
    }
    
    //Select indices to delete
    std::set<int> idxsToDelete;
    for(int i = 0; i < filteredRowIdxs.size(); i++)
    {
      if(filteredRowIdxs[i] == filteredColIdxs[i])
      {
        continue;
      }
      
      if(indexFrequency[ filteredRowIdxs[i] ] < indexFrequency[ filteredColIdxs[i] ])
      {
        if( idxsToDelete.find(filteredRowIdxs[i]) == idxsToDelete.end() )
        {
          idxsToDelete.insert(filteredColIdxs[i]);
        }
      }
      else
      {
        if( idxsToDelete.find(filteredColIdxs[i]) == idxsToDelete.end() )
        {
          idxsToDelete.insert(filteredRowIdxs[i]);
        }
      }
    }
    misc.message << idxsToDelete.size() << " groups will be filtered." << std::endl;
    
    //Which individuals must be kept?
    labelsToKeep.clear();
    for(int i = 0; i < labels.size(); i++)
    {
      if(idxsToDelete.find(i) == idxsToDelete.end())
      {
        labelsToKeep.push_back(labels[i]);
      }
    }
  }
  
  delete correlations;
  
  if( rowcolumn == column )
  {
    std::vector<std::string> unfiltered = this->effects->getRowLabels();
    this->effects->filterRowsAndCols(unfiltered, labelsToKeep);
  }
  else
  {
    std::vector<std::string> unfiltered = this->effects->getColLabels();
    this->effects->filterRowsAndCols(labelsToKeep, unfiltered);
  }
}
