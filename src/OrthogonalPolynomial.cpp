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
  const RealArray& x = gauss_points(order);
  const RealArray& w = gauss_weights(order);
  Real sum = 0.;
  for (size_t i=0; i<order; i++) {
    PCout << "Root x = " << x[i] << " Poly(x) = " << get_value(x[i], order)
	  << '\n';
    sum += w[i];
  }
  PCout << "Weights sum to " << sum << "\n\n";
}

} // namespace Pecos
