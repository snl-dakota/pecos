/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        ChebyshevOrthogPolynomial
//- Description:  Implementation code for ChebyshevOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "ChebyshevOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.H"
#endif

//#define DEBUG


namespace Pecos {

const Real& ChebyshevOrthogPolynomial::
get_value(const Real& x, unsigned short order)
{
  switch (order) {
  case 0:
    basisPolyValue = 1.;
    break;
  case 1:
    basisPolyValue = x;
    break;
  case 2:
    basisPolyValue = 2.*x*x - 1.;
    break;
  case 3:
    basisPolyValue = x*(4.*x*x - 3.);
    break;
  case 4: {
    Real x2 = x*x;
    basisPolyValue = 8.*x2*(x2 - 1.) + 1.;
    break;
  }
  case 5: {
    Real x2 = x*x;
    basisPolyValue = x*(16.*x2*x2 - 20.*x2 + 5.);
    break;
  }
  case 6: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue = 32.*x2*x4 - 48.*x4 + 18.*x2 - 1.;
    break;
  }
  case 7: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue = x*(64.*x2*x4 - 112.*x4 + 56.*x2 - 7.);
    break;
  }
  case 8: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue = 128.*x4*x4 - 256.*x2*x4 + 160.*x4 - 32.*x2 + 1.;
    break;
  }
  case 9: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyValue = x*(256.*x4*x4 - 576.*x2*x4 + 432.*x4 - 120.*x2 + 9.);
    break;
  }
  default:
    // Support higher order polynomials using the 3 point recursion formula:
    Real x2 = x*x, x4 = x2*x2, x6 = x2*x4, x8 = x4*x4,
      T_n       = x*(256.*x8 - 576.*x6 + 432.*x4 - 120.*x2 + 9.), // T_9
      T_nminus1 = 128.*x8 - 256.*x6 + 160.*x4 - 32.*x2 + 1.;      // T_8
    for (size_t i=9; i<order; i++) {
      basisPolyValue = 2.*x*T_n - T_nminus1; // T_nplus1
      if (i != order-1) {
	T_nminus1 = T_n;
	T_n       = basisPolyValue;
      }
    }
    break;
  }

  return basisPolyValue;
}


const Real& ChebyshevOrthogPolynomial::
get_gradient(const Real& x, unsigned short order)
{
  switch (order) {
  case 0:
    basisPolyGradient = 0.;
    break;
  case 1:
    basisPolyGradient = 1;
    break;
  case 2:
    basisPolyGradient = 4.*x;
    break;
  case 3:
    basisPolyGradient = 12.*x*x - 3.;
    break;
  case 4:
    basisPolyGradient = x*(32.*x*x - 16.);
    break;
  case 5: {
    Real x2 = x*x;
    basisPolyGradient = 80.*x2*x2 - 60.*x2 + 5.;
    break;
  }
  case 6: {
    Real x2 = x*x;
    basisPolyGradient = x*(192.*x2*x2 - 192.*x2 + 36.);
    break;
  }
  case 7: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyGradient = 448.*x2*x4 - 560.*x4 + 168.*x2 - 7.;
    break;
  }
  case 8: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyGradient = x*(1024.*x2*x4 - 1536.*x4 + 640.*x2 - 64.);
    break;
  }
  case 9: {
    Real x2 = x*x, x4 = x2*x2;
    basisPolyGradient = 2304.*x4*x4 - 4032.*x4*x2 + 2160.*x4 - 360.*x2 + 9.;
    break;
  }
  default:
    // Support higher order polynomials using a 3 point recursion formula:
    Real x2 = x*x, x4 = x2*x2,
      dTdx_n = 2304.*x4*x4 - 4032.*x4*x2 + 2160.*x4 - 360.*x2 + 9., // P'_9
      dTdx_nminus1 = x*(1024.*x2*x4 - 1536.*x4 + 640.*x2 - 64.);    // P'_8
    for (size_t i=9; i<order; i++) {
      // dTdx_nplus1:
      basisPolyGradient	= 2.*x*dTdx_n + 2.*get_value(x,i) - dTdx_nminus1;
      if (i != order-1) {
	dTdx_nminus1 = dTdx_n;
	dTdx_n       = basisPolyGradient;
      }
    }
    break;
  }

  return basisPolyGradient;
}


const Real& ChebyshevOrthogPolynomial::norm_squared(unsigned short order)
{
  orthogPolyNormSq = (order) ? PI/2. : PI;
  return orthogPolyNormSq;
}


const RealArray& ChebyshevOrthogPolynomial::gauss_points(unsigned short order)
{
  // pull this out from default below since order=0 is initial gauss pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "ChebyshevOrthogPolynomial::gauss_points()." << std::endl;
    abort_handler(-1);
  }

  if (gaussPoints.size() != order) { // if not already computed
    gaussPoints.resize(order);

#ifdef HAVE_SPARSE_GRID
    if (gaussWeights.size() != order)
      gaussWeights.resize(order);
    if (gaussMode == CLENSHAW_CURTIS ||
	gaussMode == CLENSHAW_CURTIS_MODERATE ||
	gaussMode == CLENSHAW_CURTIS_SLOW)
      webbur::clenshaw_curtis_compute(order, &gaussPoints[0], &gaussWeights[0]);
    else if (gaussMode == FEJER2 ||
	     gaussMode == FEJER2_MODERATE ||
	     gaussMode == FEJER2_SLOW)
      webbur::fejer2_compute(order, &gaussPoints[0], &gaussWeights[0]);
    else {
      PCerr << "Error: unsupported Gauss point type in "
	    << "ChebyshevOrthogPolynomial::gauss_points()." << std::endl;
      abort_handler(-1);
    }
    for (size_t i=0; i<order; i++)
      gaussWeights[i] *= wtFactor;
#else
    PCerr << "Error: configuration with VPISparseGrid package required in "
	  << "ChebyshevOrthogPolynomial::gauss_points()." << std::endl;
    abort_handler(-1);
#endif
  }

  return gaussPoints;
}


const RealArray& ChebyshevOrthogPolynomial::gauss_weights(unsigned short order)
{
  // The sums of the weights = 1, which is the integral of the density
  // function 1/2 over the support range of [-1,+1].  These differ from
  // VPISparseGrid by a constant factor of 1/2.

  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "ChebyshevOrthogPolynomial::gauss_weights()." << std::endl;
    abort_handler(-1);
  }

  if (gaussWeights.size() != order) { // if not already computed
    gaussWeights.resize(order);

#ifdef HAVE_SPARSE_GRID
    if (gaussPoints.size() != order)
      gaussPoints.resize(order);
    if (gaussMode == CLENSHAW_CURTIS ||
	gaussMode == CLENSHAW_CURTIS_MODERATE ||
	gaussMode == CLENSHAW_CURTIS_SLOW)
      webbur::clenshaw_curtis_compute(order, &gaussPoints[0], &gaussWeights[0]);
    else if (gaussMode == FEJER2 ||
	     gaussMode == FEJER2_MODERATE ||
	     gaussMode == FEJER2_SLOW)
      webbur::fejer2_compute(order, &gaussPoints[0], &gaussWeights[0]);
    else {
      PCerr << "Error: unsupported Gauss weight type in "
	    << "ChebyshevOrthogPolynomial::gauss_weights()." << std::endl;
      abort_handler(-1);
    }
    for (size_t i=0; i<order; i++)
      gaussWeights[i] *= wtFactor;
#else
    PCerr << "Error: configuration with VPISparseGrid package required in "
	  << "ChebyshevOrthogPolynomial::gauss_weights()." << std::endl;
    abort_handler(-1);
#endif
  }

  return gaussWeights;
}

} // namespace Pecos
