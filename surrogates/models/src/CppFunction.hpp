#ifndef CPPFUNCTION_HPP
#define CPPFUNCTION_HPP
#include <Function.hpp>

namespace Surrogates {

/**
\class CppFunction

\brief \f$C_0\f$ function that evaluates wraps a function pointer.

This class is very useful for unit testing.
*/
class CppFunction : public Function {
public:

  CppFunction();

  ~CppFunction();

  void set_function(
     void (*function)(const Real* sample, Real* func_vals,
		      const OptionsList &opts));

  /** \copydoc Function::value() */
  void value(const RealMatrix &samples, RealMatrix &values_out);

  /** \brief evaluate the function \f$f(x)\f$ at single sample x

  \param[in] sample (num_vars x 1) vector
       The coordindates of the sample x.

   \param[out] result (num_qoi x 1) vector
       The values of the vectored function at the sample x.
   */
  void value(const RealVector &sample, RealVector &result);

protected:
  /// The function \f$f(x)\f$
  void (*targetFunction_)(const Real*, Real*,
		       const OptionsList &opts);

}; // class CppFunction

}; // namespace Surrogates

#endif // CPPFUNCTION_HPP
