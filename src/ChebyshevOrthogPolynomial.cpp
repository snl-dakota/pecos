/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        ChebyshevOrthogPolynomial
//- Description:  Implementation code for ChebyshevOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "ChebyshevOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.hpp"
#endif

//#define DEBUG


namespace Pecos {

Real ChebyshevOrthogPolynomial::type1_value(const Real& x, unsigned short order)
{
  Real t1_val;
  switch (order) {
  case 0:
    t1_val = 1.;                                                          break;
  case 1:
    t1_val = x;                                                           break;
  case 2:
    t1_val = 2.*x*x - 1.;                                                 break;
  case 3:
    t1_val = x*(4.*x*x - 3.);                                             break;
  case 4: {
    Real x2 = x*x;
    t1_val = 8.*x2*(x2 - 1.) + 1.;                                        break;
  }
  case 5: {
    Real x2 = x*x;
    t1_val = x*(16.*x2*x2 - 20.*x2 + 5.);                                 break;
  }
  case 6: {
    Real x2 = x*x, x4 = x2*x2;
    t1_val = 32.*x2*x4 - 48.*x4 + 18.*x2 - 1.;                            break;
  }
  case 7: {
    Real x2 = x*x, x4 = x2*x2;
    t1_val = x*(64.*x2*x4 - 112.*x4 + 56.*x2 - 7.);                       break;
  }
  case 8: {
    Real x2 = x*x, x4 = x2*x2;
    t1_val = 128.*x4*x4 - 256.*x2*x4 + 160.*x4 - 32.*x2 + 1.;             break;
  }
  case 9: {
    Real x2 = x*x, x4 = x2*x2;
    t1_val = x*(256.*x4*x4 - 576.*x2*x4 + 432.*x4 - 120.*x2 + 9.);        break;
  }
  default:
    // Support higher order polynomials using the 3 point recursion formula:
    Real x2 = x*x, x4 = x2*x2, x6 = x2*x4, x8 = x4*x4,
      T_n       = x*(256.*x8 - 576.*x6 + 432.*x4 - 120.*x2 + 9.), // T_9
      T_nminus1 = 128.*x8 - 256.*x6 + 160.*x4 - 32.*x2 + 1.;      // T_8
    for (size_t i=9; i<order; i++) {
      t1_val = 2.*x*T_n - T_nminus1; // T_nplus1
      if (i != order-1) {
	T_nminus1 = T_n;
	T_n       = t1_val;
      }
    }
    break;
  }

  return t1_val;
}


Real ChebyshevOrthogPolynomial::
type1_gradient(const Real& x, unsigned short order)
{
  Real t1_grad;
  switch (order) {
  case 0:
    t1_grad = 0.;                                                         break;
  case 1:
    t1_grad = 1;                                                          break;
  case 2:
    t1_grad = 4.*x;                                                       break;
  case 3:
    t1_grad = 12.*x*x - 3.;                                               break;
  case 4:
    t1_grad = x*(32.*x*x - 16.);                                          break;
  case 5: {
    Real x2 = x*x;
    t1_grad = 80.*x2*x2 - 60.*x2 + 5.;                                    break;
  }
  case 6: {
    Real x2 = x*x;
    t1_grad = x*(192.*x2*x2 - 192.*x2 + 36.);                             break;
  }
  case 7: {
    Real x2 = x*x, x4 = x2*x2;
    t1_grad = 448.*x2*x4 - 560.*x4 + 168.*x2 - 7.;                        break;
  }
  case 8: {
    Real x2 = x*x, x4 = x2*x2;
    t1_grad = x*(1024.*x2*x4 - 1536.*x4 + 640.*x2 - 64.);                 break;
  }
  case 9: {
    Real x2 = x*x, x4 = x2*x2;
    t1_grad = 2304.*x4*x4 - 4032.*x4*x2 + 2160.*x4 - 360.*x2 + 9.;        break;
  }
  default:
    // Support higher order polynomials using a 3 point recursion formula:
    Real x2 = x*x, x4 = x2*x2,
      dTdx_n = 2304.*x4*x4 - 4032.*x4*x2 + 2160.*x4 - 360.*x2 + 9., // P'_9
      dTdx_nminus1 = x*(1024.*x2*x4 - 1536.*x4 + 640.*x2 - 64.);    // P'_8
    for (size_t i=9; i<order; i++) {
      // dTdx_nplus1:
      t1_grad	= 2.*x*dTdx_n + 2.*type1_value(x,i) - dTdx_nminus1;
      if (i != order-1) {
	dTdx_nminus1 = dTdx_n;
	dTdx_n       = t1_grad;
      }
    }
    break;
  }

  return t1_grad;
}


Real ChebyshevOrthogPolynomial::norm_squared(unsigned short order)
{ return (order) ? PI/2. : PI; }


const RealArray& ChebyshevOrthogPolynomial::
collocation_points(unsigned short order)
{
  // pull this out from default below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in Chebyshev"
	  << "OrthogPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  if (collocPoints.size() != order) { // if not already computed
    collocPoints.resize(order);

#ifdef HAVE_SPARSE_GRID
    // separable calculation of points/weights in sandia_rules.C
    if (collocRule == CLENSHAW_CURTIS)
      webbur::clenshaw_curtis_compute_points(order, &collocPoints[0]);
    else if (collocRule == FEJER2)
      webbur::fejer2_compute_points(order, &collocPoints[0]);
    else {
      PCerr << "Error: unsupported collocation point type in ChebyshevOrthog"
	    << "Polynomial::collocation_points()." << std::endl;
      abort_handler(-1);
    }
#else
    PCerr << "Error: configuration with VPISparseGrid package required in "
	  << "ChebyshevOrthogPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
#endif
  }

  return collocPoints;
}


const RealArray& ChebyshevOrthogPolynomial::
type1_collocation_weights(unsigned short order)
{
  // The sums of the weights = 1, which is the integral of the density
  // function 1/2 over the support range of [-1,+1].  These differ from
  // VPISparseGrid by a constant factor of 1/2.

  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in Chebyshev"
	  << "OrthogPolynomial::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  if (collocWeights.size() != order) { // if not already computed
    collocWeights.resize(order);

#ifdef HAVE_SPARSE_GRID
    // separable calculation of points/weights in sandia_rules.C
    if (collocRule == CLENSHAW_CURTIS)
      webbur::clenshaw_curtis_compute_weights(order, &collocWeights[0]);
    else if (collocRule == FEJER2)
      webbur::fejer2_compute_weights(order, &collocWeights[0]);
    else {
      PCerr << "Error: unsupported collocation weight type in ChebyshevOrthog"
	    << "Polynomial::type1_collocation_weights()." << std::endl;
      abort_handler(-1);
    }
    for (size_t i=0; i<order; i++)
      collocWeights[i] *= wtFactor;
#else
    PCerr << "Error: configuration with VPISparseGrid package required in "
	  << "ChebyshevOrthogPolynomial::type1_collocation_weights()."
	  << std::endl;
    abort_handler(-1);
#endif
  }

  return collocWeights;
}

} // namespace Pecos
