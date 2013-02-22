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
    interpolation points. */
void LagrangeInterpPolynomial::precompute_data()
{
  // precompute w_j for all forms of Lagrange interpolation
  size_t i, j, num_interp_pts = interpPts.size();
  if (bcWeights.empty())
    bcWeights.resize(num_interp_pts);
  bcWeights = 1.;
  for (i=0; i<num_interp_pts; i++) {
    Real interp_pt_i = interpPts[i];
    Real&    bc_wt_i = bcWeights[i];
    for (j=0; j<num_interp_pts; j++)
      if (i != j)
	bc_wt_i /= interp_pt_i - interpPts[j];
  }
}


/** Define the bcWeightFactors (and exactIndex if needed) corresponding to x. */
void LagrangeInterpPolynomial::set_new_point(Real x)
{
  if (x == newPoint)
    return;

  newPoint = x; exactIndex = _NPOS;
  size_t j, num_interp_pts = interpPts.size(); Real diff;

  /* precompute l(x) for 1st form of barycentric interpolation formula
  xProduct = 1.;
  for (j=0; j<num_interp_pts; ++j) {
    diff = newPoint - interpPts[j];
    if (diff == 0.)
      { exactIndex = j; break; }
    else
      xProduct *= diff;
  }
  */

  // precompute w_j/(x-x_j) for 2nd form of barycentric interpolation formula
  if (bcWeightFactors.length() != num_interp_pts)
    bcWeightFactors.resize(num_interp_pts);
  if (bcWeights.length() != num_interp_pts) {
    PCerr << "Error: length of precomputed bcWeights (" << bcWeights.length()
	  << ") is inconsistent with number of collocation points ("
	  << num_interp_pts << ")." << std::endl;
    abort_handler(-1);
  }
  bcWeightFactorSum = 0.;
  for (j=0; j<num_interp_pts; j++) {
    diff = newPoint - interpPts[j];
    if (diff == 0.) // no tolerance needed due to favorable stability analysis
      { exactIndex = j; break; }
    else
      bcWeightFactorSum += bcWeightFactors[j] = bcWeights[j] / diff;
  }
}


/** Compute value of the Lagrange polynomial (1st barycentric form)
    corresponding to interpolation point i using data from previous
    call to set_new_point().
Real LagrangeInterpPolynomial::type1_value(unsigned short i)
{
  // first form of the barycentric interpolation formula
  if (exactIndex == _NPOS)
    return bcWeights[i] * xProduct / (newPoint - interpPts[i]);
  else
    return (exactIndex == i) ? 1. : 0.;
}
*/


/** Compute derivative with respect to x of the Lagrange polynomial
    (1st barycentric form) corresponding to interpolation point i
    using data from previous call to set_new_point().
Real LagrangeInterpPolynomial::type1_gradient(unsigned short i)
{ 
  // first form of the barycentric interpolation formula
  if (exactIndex == _NPOS) {
    Real sum = 0.,
      t1_i = bcWeights[i] * xProduct / (newPoint - interpPts[i]);
    size_t j, num_interp_pts = interpPts.size();
    for (j=0; j<num_interp_pts; j++)
      if (j != i)
	sum += t1_i / (newPoint - interpPts[j]);
    return sum;
  }
  else // double loop fallback
    return type1_gradient(newPoint, i);
}
*/


/** Compute value of the Lagrange polynomial (traditional characteristic
    polynomial form) corresponding to interpolation point i. */
Real LagrangeInterpPolynomial::type1_value(Real x, unsigned short i)
{
  size_t j, num_interp_pts = interpPts.size();
  Real t1_val = bcWeights[i];
  for (j=0; j<num_interp_pts; j++)
    if (i != j)
      t1_val *= x - interpPts[j];
  return t1_val;
}


/** Compute derivative with respect to x of the Lagrange polynomial (traditional
    characteristic polynomial form) corresponding to interpolation point i. */
Real LagrangeInterpPolynomial::type1_gradient(Real x, unsigned short i)
{ 
  size_t j, k, num_interp_pts = interpPts.size();
  Real sum = 0., prod;
  for (j=0; j<num_interp_pts; j++) {
    if (j != i) {
      prod = 1.;
      for (k=0; k<num_interp_pts; k++)
	if (k != j && k != i)
	  prod *= x - interpPts[k];
      sum += prod;
    }
  }
  return sum * bcWeights[i];
}

} // namespace Pecos
