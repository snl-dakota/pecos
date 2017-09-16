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

  /** \copydoc PolyApproximation::generate_basis_matrix() */
  void generate_basis_matrix(const RealMatrix &samples, RealMatrix &result_0);

}; // class Monomial

} // namespace Surrogates
#endif // MONOMIAL_HPP
