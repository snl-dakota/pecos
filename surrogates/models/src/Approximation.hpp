#ifndef APPROXIMATION_HPP
#define APPROXIMATION_HPP

#include "Function.hpp"
#include "VariableTransformation.hpp"
#include "OptionsList.hpp"
#include "teuchos_data_types.hpp"

namespace Surrogates {

/**
\class Approximation

\brief The abstract base class of an approximation \f$p(x)\f$ of a function \f$f(x)\f$.
 */
class Approximation : public Function {
protected:
  boost::shared_ptr<VariableTransformation> varTransform_;

public:

  /// Default constructor
  Approximation(){};

  /// Destructor
  virtual ~Approximation(){};

  /** \copydoc Variables::num_vars() */
  int num_vars();

  /** \brief Set the variable transformation
  \param[in] var_transform
      The definition of the function variables x of the approximation
      and how to map them to a standardized space.

  \todo Perhaps this should be hidden from the user and instead we
  provide an interface that takes a variables object and a set of
  transformation types. Then the transformation can be constructed
  internally.
   */
  void set_variable_transformation(const boost::shared_ptr<VariableTransformation> &var_transform);


  /**
  \brief Get the variable transformation.
  \param[out] var_transform
      The definition of the function variables x of the approximation
      and how to map them to a standardized space.
   */
  boost::shared_ptr<VariableTransformation> get_variable_transformation();

}; // class Approximation

}; // namespace Surrogates

#endif // APPROXIMATION_HPP