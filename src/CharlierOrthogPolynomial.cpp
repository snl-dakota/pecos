/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "CharlierOrthogPolynomial.hpp"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace Pecos {

Real CharlierOrthogPolynomial::type1_value( Real x, unsigned short order )
{
  Real result = 1.;
  switch ( order ) {
  case 0:
    result = 1.; break;
  case 1:
    result = (lambdaStat-x)/lambdaStat; break;
  case 2:{
    Real alpha2 = lambdaStat*lambdaStat;
    result = (alpha2+x*(-2.*lambdaStat+x-1.))/alpha2;
    break;
  }
  case 3:{
    Real alpha2 = lambdaStat*lambdaStat, alpha3 = lambdaStat*alpha2;
    result=(alpha3+(-3.*alpha2+(2.+3.*lambdaStat-x)*(-1.+x))*x)/alpha3;
    break;
  }
  case 4:{
    Real alpha2 = lambdaStat*lambdaStat, alpha3 = lambdaStat*alpha2, x2 = x*x;
    result = (alpha3*(2.+lambdaStat)-2.*(3.+2.*lambdaStat)*(1.+lambdaStat+alpha2)*x+(11.+2.*lambdaStat*(7.+3.*lambdaStat))*x2-2.*(3.+2.*lambdaStat)*x2*x + x2*x2)/(alpha2*alpha2);
    break;
  }
  default:{
    // Support higher order polynomials using the 3 point recursion formula:
    Real alpha2 = lambdaStat*lambdaStat, alpha3 = lambdaStat*alpha2,x2 = x*x,
      Ch_nm1 = (alpha3+(-3.*alpha2+(2.+3.*lambdaStat-x)*(-1.+x))*x)/alpha3, //3
      Ch_n = (alpha3*(2.+lambdaStat)-2.*(3.+2.*lambdaStat)*
	      (1.+lambdaStat+alpha2)*x+(11.+2.*lambdaStat*(7.+3.*lambdaStat))
	      *x2-2.*(3.+2.*lambdaStat)*x2*x + x2*x2)/(alpha2*alpha2); //4
    for ( size_t i=4; i<order; i++ ) {
      result = ((i+lambdaStat-x)*Ch_n-i*Ch_nm1)/lambdaStat; // Ch_nplus1
      if (i != order-1) {
	Ch_nm1 = Ch_n;
	Ch_n   = result;
      }
    }
    break;
    // The recusion above produces a different result to the following recursion
    // Ch_n(x,a)=1/a*x*Ch_n(x-1,a) - Ch_n(x,a)
    // Specifically every odd polynomial is different by a factor of -1.
  }
  }
  return result;
}

Real CharlierOrthogPolynomial::type1_gradient( Real x, unsigned short order )
{
  Real result = 0.;
  switch ( order ) {
  case 0:
    result = 0.; break;
  case 1:
    result = -1/lambdaStat; break;
  case 2:{
    Real alpha2 = lambdaStat*lambdaStat;
    result = (2.*(-lambdaStat+x)-1.)/alpha2;
    break;
  }
  case 3:{
    Real alpha3 = lambdaStat*lambdaStat*lambdaStat;
    result=(-2.+(6.-3.*x)*x+lambdaStat*(-3.-3.*lambdaStat+6.*x))/alpha3;
    break;
  }
  case 4:{
    Real alpha2 = lambdaStat*lambdaStat;
    result = (-6.+lambdaStat*(-10.+(-10.-4.*lambdaStat)*lambdaStat)+x*(22.+lambdaStat*(28.+12.*lambdaStat)+x*(-18.-12.*lambdaStat+4.*x)))/(alpha2*alpha2);
    break;
  }
  default:{
    // Support higher order polynomials using the 3 point recursion formula:
    Real dChdx_nm1 = type1_gradient(x,3), dChdx_n =  type1_gradient(x,4);
    for ( size_t i=4; i<order; i++ ) {
      result = ((i+lambdaStat-x)*dChdx_n-type1_value(x,order)-i*dChdx_nm1)/lambdaStat;
      if (i != order-1) {
	dChdx_nm1 = dChdx_n;
	dChdx_n   = result;
      }
    }
    break;
  }
  }
  return result;
};

Real CharlierOrthogPolynomial::type1_hessian( Real x, unsigned short order )
{
  Real result = 0.;
  switch ( order ) {
  case 0:{
    result = 0.;
    break;
  }
  case 1:{
    result = 0.;
    break;
  }
  case 2:{
    Real alpha2 = lambdaStat*lambdaStat;
    result = 2./alpha2;
    break;
  }
  case 3:{
    Real alpha3 = lambdaStat*lambdaStat*lambdaStat;
    result=6.*(lambdaStat-x+1.)/alpha3;
    break;
  }
  case 4:{
    Real alpha2 = lambdaStat*lambdaStat;
    result = (2.*(11.+6.*alpha2+2.*lambdaStat*(7.-6.*x)+6.*(-3.+x)*x))/alpha2*alpha2;
    break;
  }
  default:{
    // Support higher order polynomials using the 3 point recursion formula:
    Real d2Chdx2_nm1 = type1_hessian(x,3), d2Chdx2_n = type1_hessian(x,4);
    for ( size_t i=4; i<order; i++ ) {
      result = ((i+lambdaStat-x)*d2Chdx2_n-2.*type1_gradient(x,order)-i*d2Chdx2_nm1)/lambdaStat;
      if (i != order-1) {
	d2Chdx2_nm1 = d2Chdx2_n;
	d2Chdx2_n   = result;
      }
    }
    break;
  }
  }
  return result;
};


Real CharlierOrthogPolynomial::norm_squared( unsigned short order )
{
  return std::pow( lambdaStat, order ) *
    boost::math::factorial<Real>( (Real)order );
}

} // namespace Pecos
