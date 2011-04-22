/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        PiecewiseInterpPolynomial
//- Description:  Implementation code for PiecewiseInterpPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "PiecewiseInterpPolynomial.hpp"

namespace Pecos {


void PiecewiseInterpPolynomial::precompute_data()
{
  //interpIntervals.size();
  //for (size_t i=0; i<; ++i)
  //  interpIntervals[i] = x[k+1] - x[k];

  // For numInterpPts == 1, could unconditionally return value = 1, grad = 0.
  // For now, enforce numInterpPts to be at least 2.
  if (numInterpPts < 2) {
    PCerr << "Error: PiecewiseInterpPolynomial requires at least two points."
	  << std::endl;
    abort_handler(-1);
  }

  // for equally spaced pts (requires at least 2 points at interval bounds):
  if (interpMode == NEWTON_COTES) {
    size_t num_intervals = numInterpPts - 1;
    interpInterval = (interpPts[num_intervals] - interpPts[0])/num_intervals;
  }
}


/** Compute value of the piecewise interpolation polynomial corresponding
    to interpolation point i. */
const Real& PiecewiseInterpPolynomial::
get_type1_value(const Real& x, unsigned short i)
{
  // does x lie within interval corresponding to interpolation point i
  const Real& pt_i = interpPts[i];
  switch (interpType) {
  case PIECEWISE_LINEAR_INTERP:
    switch (interpMode) {
    case NEWTON_COTES: {
      // linear spline interpolant with equidistant pts on [a,b]
      Real abs_dist = std::abs(x - pt_i);
      basisPolyValue = (abs_dist < interpInterval) ?
	1. - abs_dist/interpInterval : 0.;
      break;
    }
    default: {
      // linear spline interpolant with general point spacing on [a,b];
      // forward/backward looking indices protected by x{<,>}pt_i (closed rules)
      const Real& pt_im1 = interpPts[i-1]; const Real& pt_ip1 = interpPts[i+1];
      if (x <= pt_i && x > pt_im1) // for x=pt_i: val=1, not 0
	basisPolyValue = 1. - (x - pt_i)/(pt_im1 - pt_i);
      else if (x > pt_i && x < pt_ip1)
	basisPolyValue = 1. - (x - pt_i)/(pt_ip1 - pt_i);
      else
	basisPolyValue = 0.;
      break;
    }
    }
    break;
  case PIECEWISE_QUADRATIC_INTERP:
    switch (interpMode) {
    case NEWTON_COTES: {
      // quadratic spline interpolant with equidistant pts on [a,b]
      Real dist = x - pt_i;
      if (std::abs(dist) < interpInterval) {
	Real ratio = dist/interpInterval;
	basisPolyValue = 1.-ratio*ratio;
      }
      else
	basisPolyValue = 0.;
      break;
    }
    default: {
      // quadratic spline interpolant with general point spacing on [a,b]
      const Real& pt_im1 = interpPts[i-1]; const Real& pt_ip1 = interpPts[i+1];
      if (i == 0 && x < pt_ip1) {
	Real ratio = (x - pt_i)/(pt_ip1 - pt_i);
	basisPolyValue = 1.-ratio*ratio; // 1-sided as equidistant 2-sided
      }
      else if (i == numInterpPts - 1 && x > pt_im1) {
	Real ratio = (x - pt_i)/(pt_i - pt_im1);
	basisPolyValue = 1.-ratio*ratio; // 1-sided as equidistant 2-sided
      }
      else if (x > pt_im1 && x < pt_ip1)// 2-sided non-equidistant
	// Note: this interpolant does not have zero derivative at x=pt_i
	basisPolyValue = (x-pt_im1)*(pt_ip1-x)/(pt_i-pt_im1)/(pt_ip1-pt_i);
      else
	basisPolyValue = 0.;
      break;
    }
    }
    break;
  case PIECEWISE_CUBIC_INTERP: {
    // cubic Hermite spline interpolant with general point spacing
    // glue together shape fn h01 from [i-1,i] with h00 from [i,i+1]
    const Real& pt_im1 = interpPts[i-1]; const Real& pt_ip1 = interpPts[i+1];
    if (x <= pt_i && x > pt_im1) { // for x=pt_i: val=1, not 0
      // left half interval: p_k+1=1, p_k=m_k=m_k+1=0
      Real t = (x-pt_im1)/(pt_i-pt_im1); 
      basisPolyValue = t*t*(3.-2.*t);     // h01(t)
    }
    else if (x > pt_i && x < pt_ip1) {
      // right half interval: p_k=1, p_k+1=m_k=m_k+1=0
      Real t = (x - pt_i)/(pt_ip1-pt_i), tm1 = t-1.;
      basisPolyValue = tm1*tm1*(1.+2.*t); // h00(t)
    }
    else
      basisPolyValue = 0.;
    break;
  }
  }
  return basisPolyValue;
}


const Real& PiecewiseInterpPolynomial::
get_type2_value(const Real& x, unsigned short i)
{
  switch (interpType) {
  case PIECEWISE_LINEAR_INTERP: case PIECEWISE_QUADRATIC_INTERP:
    basisPolyValue = 0.;
    break;
  case PIECEWISE_CUBIC_INTERP: {
    // cubic Hermite spline interpolant with equidistant pts
    // glue together shape fn h11 from [i-1,i] with h10 from [i,i+1]
    const Real& pt_i   = interpPts[i];
    const Real& pt_im1 = interpPts[i-1]; const Real& pt_ip1 = interpPts[i+1];
    if (x < pt_i && x > pt_im1) {
      // left half interval: m_k+1=1, m_k=p_k=p_k+1=0
      Real interval = pt_i-pt_im1, t = (x-pt_im1)/interval;
      basisPolyValue = interval*t*t*(t-1.); // interval*h11(t) -> h11(\xi)
    }
    else if (x > pt_i && x < pt_ip1) {
      // right half interval: m_k=1, m_k+1=p_k=p_k+1=0
      Real interval = pt_ip1-pt_i, t = (x-pt_i)/interval, tm1 = t-1.;
      basisPolyValue = interval*tm1*tm1*t;  // interval*h10(t) -> h10(\xi)
    }
    else // includes x == pt_i
      basisPolyValue = 0.;
    break;
  }
  }
  return basisPolyValue;
}


/** Compute derivative with respect to x of the piecewise interpolation
    polynomial corresponding to interpolation point i. */
const Real& PiecewiseInterpPolynomial::
get_type1_gradient(const Real& x, unsigned short i)
{ 
  // does x lie within interval corresponding to interpolation point i
  const Real& pt_i = interpPts[i];
  switch (interpType) {
  case PIECEWISE_LINEAR_INTERP:
    switch (interpMode) {
    case NEWTON_COTES: {
      // linear spline interpolant with equidistant pts on [a,b]
      Real dist = x - pt_i, abs_dist = std::abs(dist);
      // zero gradient at transition points and outside local support
      if (abs_dist == 0. || abs_dist >= interpInterval)
	basisPolyGradient = 0.;
      else
	basisPolyGradient = (dist < 0.) ?
	  1./interpInterval : -1./interpInterval;
      break;
    }
    default: {
      const Real& pt_im1 = interpPts[i-1]; const Real& pt_ip1 = interpPts[i+1];
      if      (x < pt_i && x > pt_im1)
	basisPolyGradient =  1./(pt_i - pt_im1);
      else if (x > pt_i && x < pt_ip1)
	basisPolyGradient = -1./(pt_ip1 - pt_i);
      else // includes x == pt_i
	basisPolyGradient = 0.;
      break;
    }
    }
    break;
  case PIECEWISE_QUADRATIC_INTERP:
    switch (interpMode) {
    case NEWTON_COTES: {
      Real dist = x - pt_i;
      basisPolyGradient = (std::abs(dist) < interpInterval) ?
	-2.*dist/(interpInterval*interpInterval) : 0.; // includes dist=0 case
      break;
    }
    default: {
      const Real& pt_im1 = interpPts[i-1]; const Real& pt_ip1 = interpPts[i+1];
      if (i == 0 && x < pt_ip1) {
	Real interval = pt_ip1 - pt_i; // 1-sided as equidistant 2-sided
	basisPolyGradient = -2.*(x - pt_i)/(interval*interval);
      }
      else if (i == numInterpPts - 1 && x > pt_im1) {
	Real interval = pt_i - pt_im1; // 1-sided as equidistant 2-sided
	basisPolyGradient = -2.*(x - pt_i)/(interval*interval);
      }
      else if (x > pt_im1 && x < pt_ip1) {
	// 2-sided non-equidistant
	// Note: this interpolant does not have zero derivative at x=pt_i
	Real interval1 = pt_i - pt_im1,       interval2 = pt_ip1 - pt_i,
	  ratio1 = (x - pt_im1)/interval1, ratio2 = (pt_ip1 - x)/interval2;
	// ratio1 * dratio2/dx + ratio2 * dratio1/dx
	basisPolyGradient = ratio2/interval1 - ratio1/interval2;
      }
      else
	basisPolyGradient = 0.;
      break;
    }
    }
    break;
  case PIECEWISE_CUBIC_INTERP: {
    // cubic Hermite spline interpolant with equidistant pts
    // glue together shape fn h01 from [i-1,i] with h00 from [i,i+1]
    const Real& pt_im1 = interpPts[i-1]; const Real& pt_ip1 = interpPts[i+1];
    if (x < pt_i && x > pt_im1) {
      // left half interval: m_k+1=1, m_k=p_k=p_k+1=0
      Real interval = pt_i-pt_im1, t = (x-pt_im1)/interval, dt_dx = 1./interval;
      basisPolyGradient = 6.*t*(1.-t)*dt_dx; // dh01/dt * dt/dx
    }
    else if (x > pt_i && x < pt_ip1) {
      // right half interval: m_k=1, m_k+1=p_k=p_k+1=0
      Real interval = pt_ip1-pt_i, t = (x-pt_i)/interval, dt_dx = 1./interval;
      basisPolyGradient = 6.*t*(t-1.)*dt_dx; // dh00/dt * dt/dx
    }
    else // includes x == pt_i
      basisPolyGradient = 0.;
    break;
  }
  }
  return basisPolyGradient;
}


const Real& PiecewiseInterpPolynomial::
get_type2_gradient(const Real& x, unsigned short i)
{
  switch (interpType) {
  case PIECEWISE_LINEAR_INTERP: case PIECEWISE_QUADRATIC_INTERP:
    basisPolyGradient = 0.;
    break;
  case PIECEWISE_CUBIC_INTERP: {
    // cubic Hermite spline interpolant with equidistant pts
    // glue together shape fn h11 from [i-1,i] with h10 from [i,i+1]
    const Real& pt_i   = interpPts[i];
    const Real& pt_im1 = interpPts[i-1]; const Real& pt_ip1 = interpPts[i+1];
    if (x <= pt_i && x > pt_im1) { // for x=pt_i: grad=1, not 0
      // left half interval: m_k+1=1, m_k=p_k=p_k+1=0
      Real interval = pt_i-pt_im1, t = (x-pt_im1)/interval;//,dt_dx=1./interval;
      basisPolyGradient = t*(3.*t-2.);//*interval*dt_dx (terms cancel)
    }
    else if (x > pt_i && x < pt_ip1) {
      // right half interval: m_k=1, m_k+1=p_k=p_k+1=0
      Real interval = pt_ip1-pt_i, t = (x-pt_i)/interval;//,dt_dx = 1./interval;
      basisPolyGradient = t*(3.*t-4.)+1.;//*interval*dt_dx (terms cancel)
    }
    else
      basisPolyGradient = 0.;
    break;
  }
  }
  return basisPolyGradient;
}

} // namespace Pecos
