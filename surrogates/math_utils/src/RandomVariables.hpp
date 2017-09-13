/**
 * \file RandomVariables.hpp
 * \brief 
 * \author Russell Hooper
 */


#ifndef RANDOMVARIABLES_HPP
#define RANDOMVARIABLES_HPP

#include <Teuchos_ParameterList.hpp> 

#include <teuchos_data_types.hpp>
#include <RandomNumberGenerator.hpp>

namespace Surrogates {

/**
 * \class RandomVariables
 * \brief 
 */
class RandomVariables 
{
 public:
  /**
   * Default constructor
   */
  RandomVariables() { }

  /**
   * Destructor
   */
  ~RandomVariables() { }
  
  static void generate_samples(int num_vars, 
                               int num_samples,
      		               const Teuchos::ParameterList & params, 
                               IntMatrix &samples);

  static void generate_samples(int num_vars, 
                               int num_samples,
      		               const Teuchos::ParameterList & params, 
                               RealMatrix &samples);
};
}; // namespace Surrogates

#endif
