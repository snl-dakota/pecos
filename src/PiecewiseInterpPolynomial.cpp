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
  if (interpType == LINEAR_EQUIDISTANT || interpType == QUADRATIC_EQUIDISTANT ||
      interpType ==  CUBIC_EQUIDISTANT) {
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
  Real dist = x - pt_i;
  switch (interpType) {
  case LINEAR_EQUIDISTANT: {
    // linear spline interpolant with equidistant pts on [a,b]
    Real abs_dist = std::abs(dist);
    basisPolyValue = (abs_dist < interpInterval) ?
      1. - abs_dist/interpInterval : 0.;
    break;
  }
  case LINEAR:
    // linear spline interpolant with general point spacing on [a,b];
    // forward/backward looking indices protected by x{<,>}pt_i (closed rules)
    if      (x < pt_i && x > interpPts[i-1])
      basisPolyValue = 1. - dist/(interpPts[i-1] - pt_i);
    else if (x > pt_i && x < interpPts[i+1])
      basisPolyValue = 1. - dist/(interpPts[i+1] - pt_i);
    else
      basisPolyValue = 0.;
    break;
  case QUADRATIC_EQUIDISTANT: {
    // quadratic spline interpolant with equidistant pts on [a,b]
    if (std::abs(dist) < interpInterval) {
      Real ratio = dist/interpInterval;
      basisPolyValue = 1.-ratio*ratio;
    }
    else
      basisPolyValue = 0.;
    break;
  }
  case QUADRATIC:
    // quadratic spline interpolant with general point spacing on [a,b]
    if (i == 0 && x < interpPts[i+1]) {
      Real ratio = dist/(interpPts[i+1]-pt_i);
      basisPolyValue = 1.-ratio*ratio; // 1-sided as equidistant 2-sided
    }
    else if (i == numInterpPts - 1 && x > interpPts[i-1]) {
      Real ratio = dist/(pt_i-interpPts[i-1]);
      basisPolyValue = 1.-ratio*ratio; // 1-sided as equidistant 2-sided
    }
    else if (x > interpPts[i-1] && x < interpPts[i+1])// 2-sided non-equidistant
      // Note: this interpolant does not have zero derivative at x=pt_i
      basisPolyValue = (x    - interpPts[i-1])*(interpPts[i+1] - x)
	             / (pt_i - interpPts[i-1])/(interpPts[i+1] - pt_i);
    else
      basisPolyValue = 0.;
    break;
  case CUBIC: case CUBIC_EQUIDISTANT:
    // cubic Hermite spline interpolant with equidistant pts
    // glue together shape fn h01 from [i-1,i] with h00 from [i,i+1]
    if      (x < pt_i && x > interpPts[i-1]) { // p_{k+1}=1, p_k=m_k=m_{k+1}=0
      Real t = (x-interpPts[i-1])/(pt_i-interpPts[i-1]); // left half interval
      basisPolyValue = t*t*(3.-2.*t);     // h01(t)
    }
    else if (x > pt_i && x < interpPts[i+1]) { // p_k=1, p_{k+1}=m_k=m_{k+1}=0
      Real t = (x-pt_i)/(interpPts[i+1]-pt_i), // right half interval
	 tm1 = t-1.;
      basisPolyValue = tm1*tm1*(1.+2.*t); // h00(t)
    }
    else
      basisPolyValue = 0.;

    /* Previous code:
    Real abs_dist = std::abs(dist);
    if (abs_dist < interpInterval) {
      Real t = (x-pt_i)/interpInterval, t_sq = t*t, tm1 = t-1.,
	tm1_sq = tm1*tm1, t2 = 2.*t, h00_t = tm1_sq*(1.+t2),
	h10_t = t*tm1_sq, h01_t = t_sq*(3.-t2), h11_t = t_sq*tm1;
      Real p0 = 1., p1 = 0., m0 = 0., m1 = 0.; --> bPV = h00(t)
      basisPolyValue = h00_t*p0 + h01_t*p1;// + h10_t*m0 + h11_t*m1;
    }
    else
      basisPolyValue = 0.;
    */
    break;
  }
  return basisPolyValue;
}


const Real& PiecewiseInterpPolynomial::
get_type2_value(const Real& x, unsigned short i)
{
  switch (interpType) {
  case LINEAR_EQUIDISTANT:    case LINEAR:
  case QUADRATIC_EQUIDISTANT: case QUADRATIC:
    basisPolyValue = 0.;
    break;
  case CUBIC: case CUBIC_EQUIDISTANT: {
    // cubic Hermite spline interpolant with equidistant pts
    // glue together shape fn h11 from [i-1,i] with h10 from [i,i+1]
    const Real& pt_i = interpPts[i];
    Real dist = x - pt_i;
    if      (x < pt_i && x > interpPts[i-1]) { // m_{k+1}=1, m_k=p_k=p_{k+1}=0
      Real t = (x-interpPts[i-1])/(pt_i-interpPts[i-1]); // left half interval
      basisPolyValue = t*t*(t-1.); // h11(t)
    }
    else if (x > pt_i && x < interpPts[i+1]) { // m_k=1, m_{k+1}=p_k=p_{k+1}=0
      Real t = (x-pt_i)/(interpPts[i+1]-pt_i), // right half interval
	 tm1 = t-1.;
      basisPolyValue = tm1*tm1*t;  // h10(t)
    }
    else
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
  Real dist = x - pt_i;
  switch (interpType) {
  case LINEAR_EQUIDISTANT: {
    // linear spline interpolant with equidistant pts on [a,b]
    Real abs_dist = std::abs(dist);
    // zero gradient at transition points and outside local support
    if (abs_dist == 0. || abs_dist >= interpInterval)
      basisPolyGradient = 0.;
    else
      basisPolyGradient = (dist < 0.) ? 1./interpInterval : -1./interpInterval;
    break;
  }
  case LINEAR:
    if      (x < pt_i && x > interpPts[i-1])
      basisPolyGradient =  1./(pt_i - interpPts[i-1]);
    else if (x > pt_i && x < interpPts[i+1])
      basisPolyGradient = -1./(interpPts[i+1] - pt_i);
    else
      basisPolyGradient = 0.;
    break;
  case QUADRATIC_EQUIDISTANT:
    basisPolyGradient = (std::abs(dist) < interpInterval) ?
      -2.*dist/(interpInterval*interpInterval) : 0.;
    break;
  case QUADRATIC:
    if (i == 0 && x < interpPts[i+1]) {
      Real interval = interpPts[i+1]-pt_i; // 1-sided as equidistant 2-sided
      basisPolyGradient = -2.*dist/(interval*interval);
    }
    else if (i == numInterpPts - 1 && x > interpPts[i-1]) {
      Real interval = pt_i-interpPts[i-1]; // 1-sided as equidistant 2-sided
      basisPolyGradient = -2.*dist/(interval*interval);
    }
    else if (x > interpPts[i-1] && x < interpPts[i+1]) {
      // 2-sided non-equidistant
      // Note: this interpolant does not have zero derivative at x=pt_i
      Real interval1 = pt_i-interpPts[i-1], interval2 = interpPts[i+1]-pt_i,
	ratio1 = (x-interpPts[i-1])/interval1,
	ratio2 = (interpPts[i+1]-x)/interval2;
      // ratio1 * dratio2/dx + ratio2 * dratio1/dx
      basisPolyGradient = ratio2/interval1 - ratio1/interval2;
    }
    else
      basisPolyGradient = 0.;
    break;
  case CUBIC: case CUBIC_EQUIDISTANT:
    // cubic Hermite spline interpolant with equidistant pts
    // glue together shape fn h01 from [i-1,i] with h00 from [i,i+1]
    if (x < pt_i && x > interpPts[i-1]) {        // m_{k+1}=1, m_k=p_k=p_{k+1}=0
      const Real& pt_im1 = interpPts[i-1];
      Real interval = pt_i-pt_im1, t = (x-pt_im1)/interval,// left half interval
	dt_dx = 1./interval;
      basisPolyGradient = 6.*t*(1.-t)*dt_dx;                  // dh01/dt * dt/dx
    }
    else if (x > pt_i && x < interpPts[i+1]) {   // m_k=1, m_{k+1}=p_k=p_{k+1}=0
      const Real& pt_ip1 = interpPts[i+1];
      Real interval = pt_ip1-pt_i, t = (x-pt_i)/interval, // right half interval
	 dt_dx = 1./interval;
      basisPolyGradient = 6.*t*(t-1.)*dt_dx;                  // dh00/dt * dt/dx
    }
    else
      basisPolyGradient = 0.;
    break;
  }
  return basisPolyGradient;
}


const Real& PiecewiseInterpPolynomial::
get_type2_gradient(const Real& x, unsigned short i)
{
  switch (interpType) {
  case LINEAR_EQUIDISTANT:    case LINEAR:
  case QUADRATIC_EQUIDISTANT: case QUADRATIC:
    basisPolyGradient = 0.;
    break;
  case CUBIC: case CUBIC_EQUIDISTANT: {
    // cubic Hermite spline interpolant with equidistant pts
    // glue together shape fn h11 from [i-1,i] with h10 from [i,i+1]
    const Real& pt_i = interpPts[i];
    if (x < pt_i && x > interpPts[i-1]) {        // m_{k+1}=1, m_k=p_k=p_{k+1}=0
      const Real& pt_im1 = interpPts[i-1];
      Real interval = pt_i-pt_im1, t = (x-pt_im1)/interval,// left half interval
	dt_dx = 1./interval;
      basisPolyGradient = t*(3.*t-2.)*dt_dx;                  // dh11/dt * dt/dx
    }
    else if (x > pt_i && x < interpPts[i+1]) {   // m_k=1, m_{k+1}=p_k=p_{k+1}=0
      const Real& pt_ip1 = interpPts[i+1];
      Real interval = pt_ip1-pt_i, t = (x-pt_i)/interval, // right half interval
	 dt_dx = 1./interval;
      basisPolyGradient = (t*(3.*t-4.)+1.)*dt_dx;             // dh10/dt * dt/dx
    }
    else
      basisPolyGradient = 0.; // TO DO: x == pt_i -> grad = 1
    break;
  }
  }
  return basisPolyGradient;
}

} // namespace Pecos
