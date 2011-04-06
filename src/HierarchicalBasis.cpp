/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchicalBasis
//- Description:  Implementation code for HierarchicalBasis class
//-               
//- Owner:        Christopher Miller, University of Maryland
//- Contact:      cmiller@math.umd.edu

#include "HierarchicalBasis.hpp"

namespace Pecos{

/// Default constructor throws an exception if left_end_ >= right_end_.
HierarchicalBasis::HierarchicalBasis(RefinablePointSet& pointSet_):
  pointSet(pointSet_)
{}

HierarchicalBasis::~HierarchicalBasis()
{}

/// Point evaluation of a basis element.
const Real& HierarchicalBasis::
get_type1_value(const Real& x, const unsigned int i) 
{

  const Real nodalPoint = pointSet.get_interp_point(i);
  const Real leftEndPoint = pointSet.get_left_neighbor(i);
  const Real rightEndPoint = pointSet.get_right_neighbor(i); 
  if ( i == 0 ) {
    if ( ( x < leftEndPoint ) || ( x > rightEndPoint ) ) basisPolyValue = 0;
    else basisPolyValue = 1;
  } else {
    if ( ( x < nodalPoint ) && ( x >= leftEndPoint ) ) {
      basisPolyValue = (x - leftEndPoint)/(nodalPoint - leftEndPoint);
    } else if ( ( x >= nodalPoint ) && ( x <= rightEndPoint ) ) {
      if ( rightEndPoint == nodalPoint ) basisPolyValue = 1;
      else basisPolyValue = (rightEndPoint - x)/(rightEndPoint - nodalPoint);
    } else basisPolyValue = 0;
  }
  return basisPolyValue;
}

/// Point evaluation of the gradient of a basis element.
const Real& HierarchicalBasis::
get_type1_gradient(const Real& x, const unsigned int i) 
{
  const Real nodalPoint = pointSet.get_interp_point(i);
  const Real leftEndPoint = pointSet.get_left_neighbor(i);
  const Real rightEndPoint = pointSet.get_right_neighbor(i); 
  if ( i == 0 ) {
    basisPolyGradient = 0;
  } else {
    if ( ( x < nodalPoint ) && ( x > leftEndPoint ) ) {
      basisPolyGradient = 1/(nodalPoint - leftEndPoint);
    } else if ( ( x > nodalPoint ) && ( x < rightEndPoint ) ) {
      basisPolyGradient = -1/(rightEndPoint - nodalPoint);
    } else basisPolyGradient = 0;
  }
  return basisPolyGradient;
}

const Real& HierarchicalBasis::
get_type2_value(const Real& x, const unsigned int i)
{
  std::cout << "Function not yer implemented." << std::endl;
  basisPolyValue = 1.0;
  return basisPolyValue;
}

const Real& HierarchicalBasis::
get_type2_gradient(const Real& x, const unsigned int i)
{
  std::cout << "Function not yer implemented." << std::endl;
  basisPolyGradient = 1.0;
  return basisPolyGradient;
}

} // End namespace Pecos
