/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        RefinablePointSet
//- Description:  Implementation code for RefinablePointSet class
//-               
//- Owner:        Christopher Miller, University of Maryland
//- Contact:      cmiller@math.umd.edu

#include "RefinablePointSet.hpp"

namespace Pecos{

RefinablePointSet::RefinablePointSet(size_t map_starting_capacity):
  current_level(1),
  current_level_size(1),
  num_interp_points(1),
  interpPointsValid(false),
  highestLevelPointsValid(false),
  interpPts(1),
  highestLevelPoints(1),
  int_index_to_level_index_map(1)
{
  int_index_to_level_index_map.reserve(map_starting_capacity);
}

RefinablePointSet::~RefinablePointSet()
{
}

void RefinablePointSet::refine_all()
{
  const BoolDeque bitset(current_level_size,true);
  this->refine(bitset);
}

const Real RefinablePointSet::get_interp_point(unsigned int point) const
{
  return get_interp_point( int_index_to_level_index_map.at(point) );
}
   

const RealArray& RefinablePointSet::get_interp_points()
{
  update_interpPts();
  return interpPts;
}

const RealArray& RefinablePointSet::get_highest_level_points()
{
  update_highestLevelPoints();
  return highestLevelPoints;
}

const unsigned int RefinablePointSet::get_current_level() const
{
  return current_level;
}

const unsigned int RefinablePointSet::get_current_level_size() const
{
  return current_level_size;
}

const unsigned int RefinablePointSet::get_num_interp_points() const
{
  return num_interp_points;
}

const IntArray& RefinablePointSet::
get_level_index_pair(unsigned int i) const
{
  return int_index_to_level_index_map.at(i);
}

void RefinablePointSet::update_interpPts()
{
  if ( !interpPointsValid ) {
    interpPts.resize(num_interp_points);
    for ( unsigned int idx = 0; idx < num_interp_points; idx++ ) {
      interpPts[idx] = get_interp_point(int_index_to_level_index_map.at(idx));
    }
    std::sort( interpPts.begin(),interpPts.end() );
  }

  interpPointsValid = true;
}

void RefinablePointSet::update_highestLevelPoints()
{
  if ( !highestLevelPointsValid ) {
    highestLevelPoints.resize(current_level_size);
    IntArray levelIndex(2); 
    for ( unsigned int idx = num_interp_points-current_level_size;
          idx < num_interp_points; idx++ ) {
      levelIndex = int_index_to_level_index_map.at(idx);
      highestLevelPoints[idx - ( num_interp_points-current_level_size ) ] = 
                  get_interp_point( levelIndex );
    }
  }

  highestLevelPointsValid = true;
}

} // End namespace Pecos
