/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HahnOrthogPolynomial
//- Description:  Implementation code for HahnOrthogPolynomial class
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#include "HahnOrthogPolynomial.hpp"
#include "pecos_stat_util.hpp"


namespace Pecos {

Real HahnOrthogPolynomial::type1_value(Real x, unsigned short order)
{
  Real t1_val;

  switch (order) {
  case 0:
    t1_val = 1.;
    break;

  case 1:
    t1_val = 1.0 + (2.0+totalPop+selectPop)/((-numDrawn)*(totalPop+1.0))*x;
    break;

  case 2:
    t1_val = 1.0 - 2.0*(3.0+totalPop+selectPop)*x/(numDrawn*(totalPop+1.0)) + (3.0+totalPop+selectPop)*(4.0+totalPop+selectPop)/((totalPop+1.0)*(totalPop+2.0)*numDrawn*(numDrawn-1.0))*x*(x-1.0);
    break;

  default: {
    // Support higher order polynomials using the 3 point recursion formula:
    Real omN = 1.0 - selectPop, om1 = Real(order) - 1.0;
    // TO DO: unroll this call stack with a loop (as for other Polynomials)
    Real fm2 = type1_value(x, order-2), fm1 = type1_value(x, order-1);
    Real A = (om1+totalPop+selectPop+1.0)*(om1+totalPop+1.0)*(numDrawn-om1)/((2.0*om1+totalPop+selectPop+1.0)*(2.0*om1+totalPop+selectPop+2.0));
    Real C = om1*(om1+totalPop+selectPop+numDrawn+1.0)*(om1+selectPop)/((2.0*om1+totalPop+selectPop)*(2.0*om1+totalPop+selectPop+1.0));
    t1_val = ((A+C-x)*fm1 - C*fm2)/A;
    break;
  }
  }

  return t1_val;
}


} // namespace Pecos
