/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        GenLaguerreOrthogPolynomial
//- Description:  Implementation code for GenLaguerreOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "GenLaguerreOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.H"
#endif
#include "pecos_stat_util.hpp"

//#define DEBUG


namespace Pecos {

const Real& GenLaguerreOrthogPolynomial::
get_value(const Real& x, unsigned short order)
{
  switch (order) {
  case 0:
    basisPolyValue = 1.;
    break;
  case 1:
    basisPolyValue = -x + alphaPoly + 1.;
    break;
  case 2:
    basisPolyValue = (x*x - 2.*(alphaPoly + 2.)*x +
		      (alphaPoly+1)*(alphaPoly+2.))/2.;
    break;
  case 3: {
    Real x2 = x*x;
    basisPolyValue = (-x*x2 + 3.*(alphaPoly+3.)*x2 -
		      3.*(alphaPoly+2.)*(alphaPoly+3.)*x +
		      (alphaPoly+1.)*(alphaPoly+2.)*(alphaPoly+3.))/6.;
    break;
  }
  default: {
    // Support higher order polynomials using the 3 point recursion formula
    Real x2 = x*x,
      La_n = (-x*x2 + 3.*(alphaPoly+3.)*x2 - 3.*(alphaPoly+2.)*(alphaPoly+3.)*x
	      + (alphaPoly+1.)*(alphaPoly+2.)*(alphaPoly+3.))/6., // La_3
      La_nm1 = (x*x - 2.*(alphaPoly + 2.)*x +
		(alphaPoly+1)*(alphaPoly+2.))/2.; // La_2
    for (size_t i=3; i<order; i++) {
      basisPolyValue = ( (2.*i+1.+alphaPoly-x)*La_n - (i+alphaPoly)*La_nm1 )
	/ (i+1.); // La_np1
      if (i != order-1) {
	La_nm1 = La_n;
	La_n   = basisPolyValue;
      }
    }
    break;
  }
  }

  return basisPolyValue;
}


const Real& GenLaguerreOrthogPolynomial::
get_gradient(const Real& x, unsigned short order)
{
#ifdef DEBUG
  // See Abramowitz & Stegun, Section 22.8, p.783
  //basisPolyGradient = (order) ? (order*get_value(x, order)
  //                  - (order+alphaPoly)*get_value(x, order-1))/x : 0.;
  if (order) { // be careful with reference to changing basisPolyValue
    Real La_n = get_value(x, order), La_nminus1 = get_value(x, order-1);
    basisPolyGradient = (order*La_n - (order+alphaPoly)*La_nminus1)/x;
  }
  else
    basisPolyGradient = 0.;
  PCout << "Gen Laguerre gradient approach 1: " << basisPolyGradient << '\n';
#endif // DEBUG

  // The previous approach, while very compact, produces 0/0 = NaN at x = 0.
  // To avoid NaN issue at lower bound, differentiate the 3 pt value recursion
  // to get a 3 point gradient recursion
  switch (order) {
  case 0:
    basisPolyGradient = 0.;
    break;
  case 1:
    basisPolyGradient = -1.;
    break;
  case 2:
    basisPolyGradient = x - (alphaPoly + 2.);
    break;
  case 3:
    basisPolyGradient
      = (-x*x + 2.*(alphaPoly+3.)*x - (alphaPoly+2.)*(alphaPoly+3.) )/2.;
    break;
  default: {
    // Support higher order polynomials using the 3 point recursion formula
    Real x2 = x*x, dLadx_n
      = (-x*x + 2.*(alphaPoly+3.)*x - (alphaPoly+2.)*(alphaPoly+3.) )/2.,//L'a_3
      dLadx_nm1 = x - (alphaPoly + 2.);                                 // L'a_2
    for (size_t i=3; i<order; i++) {
      basisPolyGradient = ( (2.*i+1.+alphaPoly-x)*dLadx_n - get_value(x,i) -
			    (i+alphaPoly)*dLadx_nm1 ) / (i+1.); // dLadx_np1
      if (i != order-1) {
	dLadx_nm1 = dLadx_n;
	dLadx_n   = basisPolyGradient;
      }
    }
    break;
  }
  }
#ifdef DEBUG
  PCout << "Gen Laguerre gradient approach 2: " << basisPolyGradient << '\n';
#endif // DEBUG

  return basisPolyGradient;
}


const Real& GenLaguerreOrthogPolynomial::norm_squared(unsigned short order)
{
  // For integer alphaPoly, Gamma(alphaPoly+n+1)/n!/Gamma(alphaPoly+1)
  // = (alphaPoly+n)!/n!/alphaPoly!
  //orthogPolyNormSq = n_choose_k(alphaPoly+n,n);

  // For real alphaPoly: Gamma(alphaPoly+n+1)/Gamma(alphaPoly+1)/n!
  // = (alphaPoly+1)(alphaPoly+2)...(alphaPoly+1+(n-1))/n!
  // = pochhammer(alphaPoly+1,n)/n!
  orthogPolyNormSq = pochhammer(alphaPoly+1,order) / factorial(order);

  return orthogPolyNormSq;
}


const RealArray& GenLaguerreOrthogPolynomial::
gauss_points(unsigned short order)
{
  // pull this out from default below since order=0 is initial gauss pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "GenLaguerreOrthogPolynomial::gauss_points()." << std::endl;
    abort_handler(-1);
  }

  if (gaussPoints.size() != order) { // if not already computed
    gaussPoints.resize(order);
    switch (order) {
    case 1: // zeros of L^(alphaPoly)_1(x) for one gen Gauss-Laguerre pt:
      gaussPoints[0] =  1. + alphaPoly;
      break;
    case 2: { // zeros of L^(alphaPoly)_2(x) for two gen Gauss-Laguerre pts:
      Real srap2 = sqrt(alphaPoly + 2.);
      gaussPoints[0] = alphaPoly + 2. - srap2;
      gaussPoints[1] = alphaPoly + 2. + srap2;
      break;
    }
    default:
#ifdef HAVE_SPARSE_GRID
      if (gaussWeights.size() != order)
	gaussWeights.resize(order);
      webbur::gen_laguerre_compute(order, alphaPoly, &gaussPoints[0],
				   &gaussWeights[0]);
      const Real& wt_factor = weight_factor();
      for (size_t i=0; i<order; i++)
	gaussWeights[i] *= wt_factor; // polynomial weight fn -> PDF
#else
      PCerr << "Error: overflow in maximum quadrature order limit (2) in "
	    << "GenLaguerreOrthogPolynomial::gauss_points().  Configure with "
	    << "VPISparseGrid to extend range." << std::endl;
      abort_handler(-1);
#endif
      break;
    }
  }

  return gaussPoints;
}


const RealArray& GenLaguerreOrthogPolynomial::
gauss_weights(unsigned short order)
{
  // Derived from -(A_{n+1} gamma_n)/(A_n Phi_n'(x_i) Phi_{n+1}(x_i)),
  // which for L^(alphaPoly)(x), is Gamma(n+alphaPoly) x_i /
  // (n! (n+alphaPoly) Gamma(1+alphaPoly) (L^(alphaPoly)_{n-1}(x_i))^2).

  // The sums of the weights = 1, which is the integral of the density function
  // x^alphaPoly exp(-x)/Gamma(alpha+1) over the support range of [0,+infinity].

  if (gaussWeights.size() != order) { // if not already computed
    gaussWeights.resize(order);
    switch (order) {
    case 1: // weight for one generalized Gauss-Laguerre point:
      gaussWeights[0] = 1.;
      break;
    default:
#ifdef HAVE_SPARSE_GRID
      if (gaussPoints.size() != order)
	gaussPoints.resize(order);
      webbur::gen_laguerre_compute(order, alphaPoly, &gaussPoints[0],
				   &gaussWeights[0]);
      const Real& wt_factor = weight_factor();
      for (size_t i=0; i<order; i++)
	gaussWeights[i] *= wt_factor; // polynomial weight fn -> PDF
#else
      // define Gauss wts from Gauss pts using formula above
      const RealArray& gauss_pts = gauss_points(order);
      for (size_t i=0; i<order; i++) {
	const Real& x_i = gauss_pts[i];
	// For integer alphaPoly:
	//gaussWeights[i] = factorial_ratio(order+alphaPoly-1, order) * x_i /
	//  (order+alphaPoly) / factorial(alphaPoly) /
	//  std::pow(get_value(x_i,order-1),2);
	// For real alphaPoly:
	gaussWeights[i] = pochhammer(alphaPoly+1., order) * x_i /
	  factorial(order)/std::pow((order+alphaPoly)*get_value(x_i,order-1),2);
      }
#endif
      break;
    }
  }

  return gaussWeights;
}


const Real& GenLaguerreOrthogPolynomial::weight_factor()
{
//#ifdef HAVE_BOOST
  wtFactor = 1./bmth::tgamma(alphaPoly + 1.);
/*
#elif HAVE_GSL
  wtFactor = 1./gsl_sf_gamma(alphaPoly + 1.);
#else
  PCerr << "Error: BOOST or GSL required in GenLaguerreOrthogPolynomial::"
        << "weight_factor()." << std::endl;
  abort_handler(-1);
#endif
*/
  return wtFactor;
}

} // namespace Pecos
