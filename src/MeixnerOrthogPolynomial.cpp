/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        MeixnerOrthogPolynomial
//- Description:  Implementation code for MeixnerOrthogPolynomial class
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#include "MeixnerOrthogPolynomial.hpp"
#include "pecos_stat_util.hpp"


namespace Pecos {

Real MeixnerOrthogPolynomial::type1_value(Real x, unsigned short order)
{
  Real t1_val;

  switch (order) {
  case 0:
    t1_val = 1.;
    break;

  case 1:
    t1_val = (probPerTrial*numTrials + (probPerTrial-1.0)*x)/(probPerTrial*numTrials);
    break;

  case 2:
    t1_val = (probPerTrial*probPerTrial*numTrials*(numTrials+1.0) + (2.0*probPerTrial*(numTrials+1.0)-probPerTrial+1.0)*(probPerTrial-1.0)*x + (probPerTrial-1)*(probPerTrial-1)*x*x)/(probPerTrial*probPerTrial*numTrials*(numTrials+1.0));
    break;

  default: {
    // Support higher order polynomials using the 3 point recursion formula:
    Real om1 = Real(order-1);
    // TO DO: unroll this call stack with a loop (as for other Polynomials)
    Real fm2 = type1_value(x, order-2), fm1 = type1_value(x, order-1);
    t1_val = ((om1+(om1+numTrials)*probPerTrial+(probPerTrial-1.0)*x)*fm1 - om1*fm2)/(probPerTrial*(om1+numTrials));
    break;
  }
  }

  return t1_val;
}


} // namespace Pecos
