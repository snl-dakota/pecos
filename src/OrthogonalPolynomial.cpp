/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        OrthogonalPolynomial
//- Description:  Class implementation of base class for orthogonal polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

void OrthogonalPolynomial::gauss_check(unsigned short order)
{
  // Check that Gauss points are roots of corresponding polynomial and
  // that Gauss weights sum to inner product of weight function
  PCout << "\nUnit test for Gauss points/weights for order " << order << '\n';
  const RealArray& x = collocation_points(order);
  const RealArray& w = type1_collocation_weights(order);
  Real sum = 0.;
  for (size_t i=0; i<order; i++) {
    PCout << "Root x = " << x[i] << " Poly(x) = " << type1_value(x[i], order)
	  << '\n';
    sum += w[i];
  }
  PCout << "Weights sum to " << sum << "\n\n";
}


/** There are a number of ways to do this precomputation.  The PECOS
    approach favors memory over flops by storing nonzero Cijk only for
    unique index sets.  An alternative approach (used by Stokhos) that
    favors flops would store all non-zeros and return iterators to
    allow efficient iteration over these non-zeros. */
void OrthogonalPolynomial::
precompute_triple_products(unsigned short max_basis_order)
{
  // Since orthogonal polynomial instances may be shared among
  // multiple dimensions, check to see if this precomputation has
  // already been performed to sufficient order.
  if (tripleProductOrder && max_basis_order <= tripleProductOrder)
    return;

  // Could tailor quad rule to each ijk order: OK if lookup, but too expensive
  // if numerically generated.  Instead, retrieve a single rule of max order.
  size_t i, j, k, l, max_int_order = 3*max_basis_order,
    max_quad_order = max_int_order/2 + 1; // rounds up using truncation
  // Override any nested rule setting to ensure i=2m-1.
  short orig_rule = NO_RULE;
  if (collocRule == GENZ_KEISTER)
    { orig_rule = collocRule; collocRule = GAUSS_HERMITE; }
  else if (collocRule == GAUSS_PATTERSON)
    { orig_rule = collocRule; collocRule = GAUSS_LEGENDRE; }
  const RealArray& pts = collocation_points(max_quad_order);
  const RealArray& wts = type1_collocation_weights(max_quad_order);
  if (orig_rule) // restore
    collocRule = orig_rule;

  UShortMultiSet ijk_triple;
  Real c_ijk, norm_sq_i, norm_sq_j, tol = 1.e-12; // tol sync'd with Stokhos
  for (i=tripleProductOrder; i<max_basis_order; ++i) {
    norm_sq_i = norm_squared(i);
    for (j=0; j<=i; ++j) {
      norm_sq_j = norm_squared(j);
      for (k=0; k<=j; ++k) {
	c_ijk = 0.;
	for (l=0; l<max_quad_order; ++l) {
	  const Real& pt = pts[l];
	  c_ijk += wts[l] * type1_value(pt, i) * type1_value(pt, j)
	                  * type1_value(pt, k);
	}
	if (std::abs(c_ijk) / std::sqrt(norm_sq_i*norm_sq_j*norm_squared(k))
	    > tol) {
	  ijk_triple.clear();
	  ijk_triple.insert(i); ijk_triple.insert(j); ijk_triple.insert(k);
	  tripleProductMap[ijk_triple] = c_ijk;
	}
      }
    }
  }
  tripleProductOrder = max_basis_order;
}

} // namespace Pecos
