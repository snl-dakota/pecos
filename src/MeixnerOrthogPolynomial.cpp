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
  Real t1_val, nt = (Real)numTrials;

  switch (order) {
  case 0:
    t1_val = 1.;
    break;

  case 1:
    t1_val = (probPerTrial*nt + (probPerTrial-1.)*x)/(probPerTrial*nt);
    break;

  case 2: {
    Real ppt2 = probPerTrial*probPerTrial, pm1 = probPerTrial-1., ntp1 = nt+1.;
    t1_val = (ppt2*nt*ntp1 + (2.*probPerTrial*ntp1-pm1)*pm1*x + pm1*pm1*x*x)/
      (ppt2*nt*ntp1);
    break;
  }

  default: {
    // Support higher order polynomials using the 3 point recursion formula:
    Real om1 = (Real)order-1., ppt2 = probPerTrial*probPerTrial,
      pm1 = probPerTrial-1., pnt = probPerTrial*nt, ntp1 = nt+1.,
      ppt2ntntp1 = ppt2*nt*ntp1,  Kc_nm1 = (pnt + pm1*x) / pnt, //1
      Kc_n = (ppt2ntntp1 + x*pm1*((2.*probPerTrial*ntp1-pm1) + pm1*x))
           /  ppt2ntntp1,//2
      A = probPerTrial*(om1+nt);
    for (size_t i=3; i<order; i++) {
      t1_val = ((om1 + A + pm1*x)*Kc_n - om1*Kc_nm1) / A; // Kc_nplus1
      if (i != order-1) {
	Kc_nm1 = Kc_n;
	Kc_n   = t1_val;
      }
    }
    break;
  }
  }

  return t1_val;
}


} // namespace Pecos
