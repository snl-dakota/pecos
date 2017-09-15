#ifndef POLYNOMIAL_CHAOS_EXPANSION_WRAP_HPP
#define POLYNOMIAL_CHAOS_EXPANSION_WRAP_HPP

#include "PolynomialApproximation.hpp"
#include "BasisApproximation.hpp"
#include "OptionsList.hpp"
#include "teuchos_data_types.hpp"
#include "pecos_data_types.hpp"
#include "SharedBasisApproxData.hpp"

namespace Surrogates{

/**
\class PolynomialChaosExpansionWrap

\brief Polynomial chaos expansion of a \f$L^2\f$ function, i.e a function with finite variance.
*/
class  PolynomialChaosExpansionWrap: public PolynomialApproximation{
protected:
/// TODO Do I need this???
Pecos::SharedBasisApproxData sharedData_;

/// The Pecos basis approximation which is being wrapped
Pecos::BasisApproximation poly_;

/// Model specific options
OptionsList opts_;

/// The number of QoI (outputs) of the vector-valued function
int numQOI_;

public:
  /// Default constructtor
  PolynomialChaosExpansionWrap();

  /// Destructor
  virtual ~PolynomialChaosExpansionWrap();


  /** \copydoc PolynomialApproximation::generate_basis_matrix() */
  void generate_basis_matrix(const RealMatrix &samples, RealMatrix &result_0);

}; // class PolynomialChaosExpansionWrap

} // namespace Surrogates

#endif // POLYNOMIAL_CHAOS_EXPANSION_WRAP_HPP
