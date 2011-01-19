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
get_value(const Real& x, unsigned short i)
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
    Real abs_dist = std::abs(dist);
    basisPolyValue = (abs_dist < interpInterval) ? 
      (dist/interpInterval - 1.)*(-dist/interpInterval - 1.) : 0.;
    break;
  }
  case QUADRATIC:
    // quadratic spline interpolant with general point spacing on [a,b]
    if (i == 0 && x < interpPts[i+1]) {
      Real interval = interpPts[i+1] - pt_i; // 1-sided as equidistant 2-sided
      basisPolyValue = (dist/interval - 1.)*(-dist/interval - 1.);
    }
    else if (i == numInterpPts - 1 && x > interpPts[i-1]) {
      Real interval = pt_i - interpPts[i-1]; // 1-sided as equidistant 2-sided
      basisPolyValue = (dist/interval-1.)*(-dist/interval-1.);
    }
    else if (x > interpPts[i-1] && x < interpPts[i+1])// 2-sided non-equidistant
      // Note: this interpolant does not have zero derivative at x=pt_i
      basisPolyValue = (x    - interpPts[i-1])*(interpPts[i+1] - x)
	             / (pt_i - interpPts[i-1])/(interpPts[i+1] - pt_i);
    else
      basisPolyValue = 0.;
    break;
  case CUBIC_EQUIDISTANT: {
    // cubic Hermite spline interpolant with equidistant pts
    Real abs_dist = std::abs(dist);
    if (abs_dist < interpInterval) { // TO DO: verify
      Real t = (x-pt_i)/interpInterval, t_sq = t*t, tm1 = t-1.,
	tm1_sq = tm1*tm1, t2 = 2.*t, h00_t = tm1_sq*(1.+t2),
	h10_t = t*tm1_sq, h01_t = t_sq*(3.-t2), h11_t = t_sq*tm1;
      Real p0 = 1., m0 = 0., p1 = 0., m1 = 0.; // TO DO
      basisPolyValue = h00_t*p0 + h10_t*m0 + h01_t*p1 + h11_t*m1;
    }
    else
      basisPolyValue = 0.;
    break;
  }
  case CUBIC:
    // cubic Hermite spline interpolant with general point spacing
    // TO DO
    break;
  }
  return basisPolyValue;
}


/** Compute derivative with respect to x of the piecewise interpolation
    polynomial corresponding to interpolation point i. */
const Real& PiecewiseInterpPolynomial::
get_gradient(const Real& x, unsigned short i)
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
  //case LINEAR:
  //case QUADRATIC_EQUIDISTANT:
  //case QUADRATIC:
  //case CUBIC_EQUIDISTANT:
  //case CUBIC:
  }
  }
  return basisPolyGradient;
}

} // namespace Pecos
