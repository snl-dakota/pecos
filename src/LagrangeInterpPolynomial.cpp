/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LagrangeInterpPolynomial
//- Description:  Implementation code for LagrangeInterpPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "LagrangeInterpPolynomial.hpp"


namespace Pecos {

/** Pre-compute denominator products that are only a function of the
    interpolationPts. */
void LagrangeInterpPolynomial::precompute_data()
{
  if (lagDenominators.empty())
    lagDenominators.resize(numInterpPts);
  lagDenominators = 1.;
  size_t i, j;
  for (i=0; i<numInterpPts; i++) {
    const Real& interp_pt_i = interpolationPts[i];
    Real&       lag_denom_i = lagDenominators[i];
    for (j=0; j<numInterpPts; j++)
      if (i != j)
	lag_denom_i *= interp_pt_i - interpolationPts[j];
  }
}


/** Compute value of Lagrange polynomial for interpolation point i. */
const Real& LagrangeInterpPolynomial::get_value(const Real& x, unsigned short i)
{
  basisPolyValue = 1. / lagDenominators[i];
  for (size_t j=0; j<numInterpPts; j++)
    if (i != j)
      basisPolyValue *= x - interpolationPts[j];
  return basisPolyValue;
}


/** Compute derivative with respect to x of Lagrange polynomial for
    interpolation point i. */
const Real& LagrangeInterpPolynomial::
get_gradient(const Real& x, unsigned short i)
{ 
  size_t j, k;
  Real numer = 0.;
  for (j=0; j<numInterpPts; j++) {
    if (j != i) {
      Real prod = 1.;
      for (k=0; k<numInterpPts; k++)
	if (k != j && k != i)
	  prod *= x - interpolationPts[k];
      numer += prod;
    }
  }
  basisPolyGradient = numer / lagDenominators[i];
  return basisPolyGradient;
}

} // namespace Pecos
