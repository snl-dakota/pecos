/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        KrawtchoukOrthogPolynomial
//- Description:  Implementation code for KrawtchoukOrthogPolynomial class
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#include "KrawtchoukOrthogPolynomial.hpp"
#include "pecos_stat_util.hpp"


namespace Pecos {

Real KrawtchoukOrthogPolynomial::type1_value(Real x, unsigned short order)
{
  Real t1_val;

  switch (order) {
  case 0:
    t1_val = 1.;
    break;

  case 1:
    t1_val = (probPerTrial*Real(numTrials) - x)/(probPerTrial*Real(numTrials));
    break;

  case 2: {
    Real omN = 1. - numTrials;
    t1_val = (probPerTrial*probPerTrial*numTrials*omN + (1.0-2.0*probPerTrial*omN)*x - x*x)/(probPerTrial*probPerTrial*numTrials*omN);
    break;
  }
  default: {
    // Support higher order polynomials using the 3 point recursion formula:
    Real om1 = Real(order) - 1.;
    // TO DO: unroll this call stack with a loop (as for other Polynomials)
    Real fm2 = type1_value(x, order-2), fm1 = type1_value(x, order-1);
    t1_val = ((probPerTrial*(numTrials-om1)+om1*(1.0-probPerTrial)-x)*fm1 - om1*(1.0-probPerTrial)*fm2)/(probPerTrial*(numTrials-om1));
    break;
  }
  }

  return t1_val;
}


} // namespace Pecos
