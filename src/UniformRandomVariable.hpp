/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 UniformRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef UNIFORM_RANDOM_VARIABLE_HPP
#define UNIFORM_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for uniform random variables.

/** Manages lower and upper bounds. */

class UniformRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  UniformRandomVariable();                   ///< default constructor
  UniformRandomVariable(Real lwr, Real upr); ///< alternate constructor
  ~UniformRandomVariable();                  ///< destructor

  //
  //- Heading: Member functions
  //

  Real cdf(Real x) const;
  Real cdf_inverse(Real p) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;

protected:

  //
  //- Heading: Data
  //

  /// lower bound of uniform random variable
  Real lowerBnd;
  /// upper bound of uniform random variable
  Real upperBnd;
};


inline UniformRandomVariable::UniformRandomVariable():
  RandomVariable(BaseConstructor()), lowerBnd(-1.), upperBnd(1.)
{ }


inline UniformRandomVariable::UniformRandomVariable(Real lwr, Real upr):
  RandomVariable(BaseConstructor()), lowerBnd(lwr), upperBnd(upr)
{ }


inline UniformRandomVariable::~UniformRandomVariable()
{ }


inline Real UniformRandomVariable::cdf(Real x) const
{ return uniform_cdf(x, lowerBnd, upperBnd); }


inline Real UniformRandomVariable::cdf_inverse(Real p) const
{ return lowerBnd + (upperBnd - lowerBnd) * p; }


inline Real UniformRandomVariable::pdf(Real x) const
{ return uniform_pdf(lowerBnd, upperBnd); }


inline Real UniformRandomVariable::pdf_gradient(Real x) const
{ return 0.; }


inline Real UniformRandomVariable::pdf_hessian(Real x) const
{ return 0.; }

} // namespace Pecos

#endif
