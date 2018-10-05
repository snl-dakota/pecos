/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef MONOMIAL_HPP
#define MONOMIAL_HPP
#include "PolyApproximation.hpp"

namespace Surrogates {
/**
\class Monomial
\brief A multivariate monomial approximation.

This class was generated mainly for unit-testing purposes. Much greater
functionality can be reached by including the PECOS library and utilizing
the polynomial chaos wrappers.
*/
class Monomial : public PolyApproximation {
public:
  Monomial();

  ~Monomial();

  void set_options(const OptionsList &opts);

  /** \copydoc PolyApproximation::generate_canonical_basis_matrix() */
  void generate_canonical_basis_matrix(const RealMatrix &samples, RealMatrix &result_0);

}; // class Monomial

} // namespace Surrogates
#endif // MONOMIAL_HPP
