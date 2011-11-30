/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        JacobiOrthogPolynomial
//- Description:  Implementation code for JacobiOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "JacobiOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.H"
#endif
#include "pecos_stat_util.hpp"

//#define DEBUG


namespace Pecos {

Real JacobiOrthogPolynomial::type1_value(const Real& x, unsigned short order)
{
  Real t1_val;
  switch (order) {
  case 0:
    t1_val = 1.;                                                          break;
  case 1:
    t1_val = (alphaPoly + betaPoly + 2.)*(x-1.)/2. + alphaPoly + 1.;      break;
  case 2: {
    Real xm1 = x - 1.;
    t1_val = ( (alphaPoly + betaPoly + 3.)*(alphaPoly + betaPoly + 4.)*xm1*xm1 +
	       4.*(alphaPoly + betaPoly + 3.)*(alphaPoly + 2.)*xm1 +
	       4.*(alphaPoly + 1.)*(alphaPoly + 2.) ) / 8.;
    break;
  }
  default: {
    // Support higher order polynomials using the 3 point recursion formula:
    Real xm1 = x - 1.,
      Pab_n = ( (alphaPoly + betaPoly + 3.)*(alphaPoly + betaPoly + 4.)*xm1*xm1
		+ 4*(alphaPoly + betaPoly + 3.)*(alphaPoly + 2.)*xm1
		+ 4.*(alphaPoly + 1.)*(alphaPoly + 2.) ) / 8.,       // Pab_2
      Pab_nm1 = (alphaPoly + betaPoly + 2.)*xm1/2. + alphaPoly + 1.; // Pab_1
    for (size_t i=2; i<order; i++) {
      Real ab2i = alphaPoly + betaPoly + 2.*i;
      t1_val // Pab_np1
	= ( ( (ab2i+1.)*(alphaPoly*alphaPoly-betaPoly*betaPoly)
	      + x*pochhammer(ab2i,3) )*Pab_n
	    - 2.*(i+alphaPoly)*(i+betaPoly)*(ab2i+2.)*Pab_nm1 )
	/ ( 2.*(i+1.)*(i+alphaPoly+betaPoly+1.)*ab2i );
      if (i != order-1) {
	Pab_nm1 = Pab_n;
	Pab_n   = t1_val;
      }
    }
    break;
  }
  }

  return t1_val;
}


Real JacobiOrthogPolynomial::type1_gradient(const Real& x, unsigned short order)
{
  Real t1_grad;
#ifdef DEBUG
  // See Abramowitz & Stegun, Section 22.8, p.783
  //t1_grad = (order) ?
  //  (g1*type1_value(x,order) + g0*type1_value(x,order-1))/g2 :0.;
  if (order) {
    Real ab2n = 2.*order+alphaPoly+betaPoly,
      g0 = 2.*(order+alphaPoly)*(order+betaPoly),
      g1 = order*(alphaPoly-betaPoly-ab2n*x), g2 = ab2n*(1.-x*x),
      Pab_n = type1_value(x, order), Pab_nminus1 = type1_value(x, order-1);
    t1_grad = (g1*Pab_n + g0*Pab_nminus1)/g2;
  }
  else
    t1_grad = 0.;
  PCout << "Jacobi gradient approach 1: " << t1_grad << '\n';
#endif // DEBUG

  // The previous approach, while very compact, produces 0/0 = NaN at x = +/-1.
  // To avoid NaN issues at bounds, differentiate the 3 pt value recursion
  // to get a 3 point gradient recursion
  switch (order) {
  case 0:
    t1_grad = 0.;                                                         break;
  case 1:
    t1_grad = (alphaPoly + betaPoly + 2.)/2.;                             break;
  case 2: {
    Real xm1 = x - 1.;
    t1_grad = ( (alphaPoly + betaPoly + 3.)*(alphaPoly + betaPoly + 4.)*xm1 +
		2.*(alphaPoly + betaPoly + 3.)*(alphaPoly + 2.) ) / 4.;
    break;
  }
  default: {
    // Support higher order polynomials using the 3 point recursion formula:
    Real xm1 = x - 1.,
      dPabdx_n = ( (alphaPoly + betaPoly + 3.)*(alphaPoly + betaPoly + 4.)*xm1 +
		   2.*(alphaPoly + betaPoly + 3.)*(alphaPoly + 2.) ) / 4.,//P'_2
      dPabdx_nm1 = (alphaPoly + betaPoly + 2.)/2.;                       // P'_1
    for (size_t i=2; i<order; i++) {
      Real ab2i = alphaPoly + betaPoly + 2.*i, pab2i3 = pochhammer(ab2i, 3);
      t1_grad // dPabdx_np1
	= ( ( (ab2i+1.)*(alphaPoly*alphaPoly-betaPoly*betaPoly) + x*pab2i3 )
	    * dPabdx_n + pab2i3*type1_value(x,i)
	    - 2.*(i+alphaPoly)*(i+betaPoly)*(ab2i+2.)*dPabdx_nm1 )
	/ ( 2.*(i+1.)*(i+alphaPoly+betaPoly+1.)*ab2i );
      if (i != order-1) {
	dPabdx_nm1 = dPabdx_n;
	dPabdx_n   = t1_grad;
      }
    }
    break;
  }
  }
#ifdef DEBUG
  PCout << "Jacobi gradient approach 2: " << t1_grad << '\n';
#endif // DEBUG

  return t1_grad;
}


Real JacobiOrthogPolynomial::norm_squared(unsigned short order)
{
  return (alphaPoly+betaPoly+1.) / (2.*order+alphaPoly+betaPoly+1.)
    * pochhammer(alphaPoly+1.,order) * pochhammer(betaPoly+1.,order)
    / pochhammer(alphaPoly+betaPoly+1.,order) / factorial(order);
}


const RealArray& JacobiOrthogPolynomial::
collocation_points(unsigned short order)
{
  // pull this out from default below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "JacobiOrthogPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  if (collocPoints.size() != order) { // if not already computed
    collocPoints.resize(order);
    switch (order) {
    case 1: // zeros of Pab_1(x) for one Gauss-Jacobi point:
      collocPoints[0] = (betaPoly - alphaPoly) / (alphaPoly + betaPoly + 2.);
      break;
    case 2: { // zeros of Pab_2(x) for two Gauss-Jacobi points:
      Real a = (alphaPoly+betaPoly+3.)*(alphaPoly+betaPoly+4.),
	   b = 4.*(alphaPoly+betaPoly+3.)*(alphaPoly+2.),
	   c = 4.*(alphaPoly+1.)*(alphaPoly+2.),
	   srdiscrim = std::sqrt(b*b-4.*a*c), a2 = 2.*a;
      collocPoints[0] = 1. - (b+srdiscrim)/a2;
      collocPoints[1] = 1. - (b-srdiscrim)/a2;
      break;
    }
    default:
#ifdef HAVE_SPARSE_GRID
      // sandia_rules.C calculates points/weights together
      if (collocWeights.size() != order)
	collocWeights.resize(order);
      webbur::jacobi_compute(order, alphaPoly, betaPoly, &collocPoints[0],
			     &collocWeights[0]);
      const Real& wt_factor = weight_factor();
      for (size_t i=0; i<order; i++)
	collocWeights[i] *= wt_factor; // polynomial weight fn -> PDF
#else
      PCerr << "Error: overflow in maximum quadrature order limit (2) in "
	    << "JacobiOrthogPolynomial::collocation_points().  Configure with "
	    << "VPISparseGrid to extend range." << std::endl;
      abort_handler(-1);
#endif
      break;
    }
  }

  return collocPoints;
}


const RealArray& JacobiOrthogPolynomial::
type1_collocation_weights(unsigned short order)
{
  // Derived from (A_n gamma_{n-1})/(A_{n-1} Phi_n'(x_i) Phi_{n-1}(x_i))

  // The sums of the weights = 1, which is the integral of the density
  // function (1-x)^alpha (1+x)^beta / (2^(alpha+beta+1) B(alpha+1,beta+1))
  // over the support range of [-1,+1].

  if (collocWeights.size() != order) { // if not already computed
    collocWeights.resize(order);
    switch (order) {
    case 1: // weights for one Gauss-Jacobi point:
      collocWeights[0] = 1.0;
      break;
    default:
#ifdef HAVE_SPARSE_GRID
      // sandia_rules.C calculates points/weights together
      if (collocPoints.size() != order)
	collocPoints.resize(order);
      webbur::jacobi_compute(order, alphaPoly, betaPoly, &collocPoints[0],
			     &collocWeights[0]);
      const Real& wt_factor = weight_factor();
      for (size_t i=0; i<order; i++)
	collocWeights[i] *= wt_factor; // polynomial weight fn -> PDF
#else
      // define Gauss wts from Gauss pts using formula above
      const RealArray& colloc_pts = collocation_points(order);
      for (size_t i=0; i<order; i++) {
	const Real& x_i = colloc_pts[i];
	Real AnoAnm1
	  = (2.*order+alphaPoly+betaPoly) * (2.*order+alphaPoly+betaPoly-1.)
	  / (2.*order) / (order+alphaPoly+betaPoly);
	collocWeights[i]
	  = AnoAnm1 * norm_squared(order-1) / type1_value(x_i, order-1)
	  / type1_gradient(x_i, order);
      }
#endif
      break;
    }
  }

  return collocWeights;
}


const Real& JacobiOrthogPolynomial::weight_factor()
{
//#ifdef HAVE_BOOST
  wtFactor = 1. / std::pow(2., alphaPoly + betaPoly + 1.) /
    bmth::beta(alphaPoly + 1., betaPoly + 1.);
/*
#elif HAVE_GSL
  wtFactor = 1. / std::pow(2., alphaPoly + betaPoly + 1.) /
    gsl_sf_beta(alphaPoly + 1., betaPoly + 1.);
#else
  PCerr << "Error: GSL required in JacobiOrthogPolynomial::weight_factor()."
        << std::endl;
  abort_handler(-1);
#endif
*/
  return wtFactor;
}

} // namespace Pecos
