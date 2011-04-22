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
  HierarchicalBasis::HierarchicalBasis(RefinablePointSet& pointSet_, short interpType_):
    PiecewiseInterpPolynomial(interpType_,CLENSHAW_CURTIS),
    pointSet(pointSet_),
    interpType(interpType_)
{
  interpPts = (pointSet.get_interp_points());
  if ( interpType == PIECEWISE_QUADRATIC_INTERP ) {
    PCerr << "Quadratic interpolation not currently implemented in"
	  << " HierarchicalBasis. Defaulting to linear interpolation.";
    interpType = PIECEWISE_LINEAR_INTERP;
  }
}

HierarchicalBasis::~HierarchicalBasis()
{}

/// Point evaluation of a basis element.
const Real& HierarchicalBasis::
get_type1_value(const Real& x, const unsigned int i) 
{

  const Real nodalPoint = pointSet.get_interp_point(i);
  const Real leftEndPoint = pointSet.get_left_neighbor(i);
  const Real rightEndPoint = pointSet.get_right_neighbor(i);

  switch (interpType) {
    case PIECEWISE_LINEAR_INTERP: {
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
      break;
    }
    case PIECEWISE_CUBIC_INTERP: {
      if ( i == 0 ) {
        if ( ( x < leftEndPoint ) || ( x > rightEndPoint ) ) basisPolyValue = 0;
        else basisPolyValue = 1;
      } else {
        if ( ( x < nodalPoint ) && ( x >= leftEndPoint ) ) {
	  Real t = ( x - leftEndPoint )/( nodalPoint - leftEndPoint );
  	  basisPolyValue = t*t*(3.-2.*t);
        }
        else if ( ( x >= nodalPoint ) && ( x < rightEndPoint ) ) {
	  Real t = ( x - nodalPoint )/( rightEndPoint - nodalPoint );
	  Real tm1 = t - 1.;
	  basisPolyValue = tm1*tm1*(1.+2.*t);
        } else if ( ( x == nodalPoint ) && ( x == rightEndPoint ) ) {
	  basisPolyValue = 1;
	}
        else basisPolyValue = 0.;
      }
      break;
    }
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
  switch (interpType) {
    case PIECEWISE_LINEAR_INTERP: {
      if ( i == 0 ) {
	basisPolyGradient = 0;
      } else {
	if ( ( x < nodalPoint ) && ( x > leftEndPoint ) ) {
	  basisPolyGradient = 1/(nodalPoint - leftEndPoint);
	} else if ( ( x > nodalPoint ) && ( x < rightEndPoint ) ) {
	  basisPolyGradient = -1/(rightEndPoint - nodalPoint);
	} else basisPolyGradient = 0;
      }
      break;
    }
    case PIECEWISE_CUBIC_INTERP: {
      if ( i == 0 ) {
	basisPolyGradient = 0;
      } else { 
	if ( ( x < nodalPoint ) && ( x > leftEndPoint ) ) {
	  Real interval = nodalPoint - leftEndPoint;
	  Real t = ( x - leftEndPoint )/interval;
	  Real dt_dx = 1./interval;
	  basisPolyGradient = 6.*t*(1.-t)*dt_dx; // dh01/dt * dt/dx
	}
	else if ( ( x > nodalPoint ) && ( x < rightEndPoint ) ) {
	  Real interval = rightEndPoint - nodalPoint;
	  Real t = ( x - nodalPoint )/interval;
	  Real dt_dx = 1./interval;
	  basisPolyGradient = 6.*t*(t-1.)*dt_dx; // dh00/dt * dt/dx
	}
	else basisPolyGradient = 0;
      
      }
      break;
    }
  }
  return basisPolyGradient;
}

const Real& HierarchicalBasis::
get_type2_value(const Real& x, const unsigned int i)
{
  const Real nodalPoint = pointSet.get_interp_point(i);
  const Real leftEndPoint = pointSet.get_left_neighbor(i);
  const Real rightEndPoint = pointSet.get_right_neighbor(i);
  
  switch (interpType) {
    case PIECEWISE_CUBIC_INTERP: {
      if ( i == 0 ) {
	basisPolyValue = 0;
      }
      else if (x <= nodalPoint && x > leftEndPoint) {
        // left half interval: m_k+1=1, m_k=p_k=p_k+1=0
	Real interval = nodalPoint - leftEndPoint;
	Real t = (x-leftEndPoint)/interval;
	basisPolyValue = interval*t*t*(t-1.); // interval*h11(t) -> h11(\xi)
      }
      // final condition here prevents double evaluations in the general case
      else if (x >= nodalPoint && x < rightEndPoint ) { 
	// right half interval: m_k=1, m_k+1=p_k=p_k+1=0
	Real interval = rightEndPoint-nodalPoint;
	Real t = (x-nodalPoint)/interval;
	Real tm1 = t-1.;
	basisPolyValue = interval*tm1*tm1*t;  // interval*h10(t) -> h10(\xi)
      }
      else
	basisPolyValue = 0.;
      break;
    }
    default: {
      basisPolyValue = 0;
      break;
    }
  }   
  return basisPolyValue;
}

const Real& HierarchicalBasis::
get_type2_gradient(const Real& x, const unsigned int i)
{
  const Real nodalPoint = pointSet.get_interp_point(i);
  const Real leftEndPoint = pointSet.get_left_neighbor(i);
  const Real rightEndPoint = pointSet.get_right_neighbor(i);
  
  switch (interpType) {
    case PIECEWISE_CUBIC_INTERP: {
      if ( i == 0 ) {
	basisPolyGradient = 0;
      } else {
        if (x < nodalPoint && x > leftEndPoint) {
	  // left half interval: m_k+1=1, m_k=p_k=p_k+1=0
	  Real interval = nodalPoint - leftEndPoint;
	  Real t = (x-leftEndPoint)/interval;
	  basisPolyGradient = t*(3.*t-2.); 
        }
        else if (x > nodalPoint && x < rightEndPoint) {
	  // right half interval: m_k=1, m_k+1=p_k=p_k+1=0
	  Real interval = rightEndPoint-nodalPoint;
	  Real t = (x-nodalPoint)/interval;
	  basisPolyGradient = t*(3.*t-4.)+1;
        }
	else if ( x == nodalPoint ) {
	  basisPolyGradient = 1;
	}
        else
	  basisPolyGradient = 0.;
        break;
      }
    }
    default: {
      basisPolyGradient = 0;
      break;
    }
  }
  return basisPolyGradient;
}

} // End namespace Pecos
