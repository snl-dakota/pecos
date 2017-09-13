#ifndef POLYNOMIAL_CHAOS_EXPANSION_HPP
#define POLYNOMIAL_CHAOS_EXPANSION_HPP

#include <PolynomialApproximation.hpp>

namespace Surrogates{

/**
\class PolynomialChaosExpansion

\brief Polynomial chaos expansion of a \f$L^2\f$ function, i.e a function with finite variance.
*/
  class PolynomialChaosExpansion : public PolynomialApproximation{
protected:
  /// Model specific options
  Teuchos::ParameterList opts_;

  /// The number of QoI (outputs) of the vector-valued function
  int numQOI_;

public:
  /// Default constructtor
  PolynomialChaosExpansion();

  /// Destructor
  virtual ~PolynomialChaosExpansion();

  /** \copydoc Function::value() */
  virtual void value(const RealMatrix &samples, RealMatrix &values) = 0;

  /** \copydoc Function::gradient() */
  virtual void gradient(const RealMatrix &samples, int qoi, RealMatrix &gradients);

  /** \copydoc Function::jacobian() */
  virtual void jacobian(const RealVector &sample, RealMatrix &jacobian);

  /** \copydoc Function::hessian() */
  virtual void hessian(const RealMatrix &samples, int qoi, RealMatrixList &hessians);

  /**\brief set options specific to the model
   */
  virtual void set_options(const Teuchos::ParameterList &opts);

  /**\brief get the model specific options
   */
  virtual void get_options(Teuchos::ParameterList &opts);

}; // class PolynomialChaosExpansion

}; // namespace Surrogates

#endif // POLYNOMIAL_CHAOS_EXPANSION_HPP
