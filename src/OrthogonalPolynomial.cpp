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


void OrthogonalPolynomial::precompute_triple_products(unsigned short max_order)
{
  Real c_ijk, tol = 1.e-10;
  size_t i, j, k, l, max_int_order = 3*max_order,
    max_quad_order = max_int_order/2 + 1; // rounds up using truncation
  UShortMultiSet ijk_triple;
  // Could tailor quad rule to each ijk order: OK if lookup, but too expensive
  // if numerically generated.  Instead, evaluate a single rule of max order.
  const RealArray& pts = collocation_points(max_quad_order);
  const RealArray& wts = type1_collocation_weights(max_quad_order);
  for (i=0; i<max_order; ++i)
    for (j=0; j<=i; ++j)
      for (k=0; k<=j; ++k) {
	c_ijk = 0.;
	for (l=0; l<max_quad_order; ++l) {
	  const Real& pt = pts[l];
	  c_ijk += wts[l] * type1_value(pt, i) * type1_value(pt, j)
	                  * type1_value(pt, k);
	}
	if (std::abs(c_ijk) > tol) {
	  ijk_triple.clear();
	  ijk_triple.insert(i); ijk_triple.insert(j); ijk_triple.insert(k);
	  tripleProductMap[ijk_triple] = c_ijk;
	}
      }
}


Real OrthogonalPolynomial::triple_product(const UShortMultiSet& ijk_key) const
{
  UShortMultiSetRealMap::const_iterator cit = tripleProductMap.find(ijk_key);
  return (cit == tripleProductMap.end()) ? 0. : cit->second;
}

} // namespace Pecos
