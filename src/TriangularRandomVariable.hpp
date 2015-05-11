/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 TriangularRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef TRIANGULAR_RANDOM_VARIABLE_HPP
#define TRIANGULAR_RANDOM_VARIABLE_HPP

#include "UniformRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for triangular random variables.

/** Manages mode and inherits bounds. */

class TriangularRandomVariable: public UniformRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  TriangularRandomVariable();
  /// alternate constructor
  TriangularRandomVariable(Real mode, Real lwr, Real upr);
  /// destructor
  ~TriangularRandomVariable();

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

  /// mode of triangular random variable
  Real triangularMode;
};


inline TriangularRandomVariable::TriangularRandomVariable():
  UniformRandomVariable(), triangularMode(0)
{ }


inline TriangularRandomVariable::
TriangularRandomVariable(Real mode, Real lwr, Real upr):
  UniformRandomVariable(lwr, upr), triangularMode(mode)
{ }


inline TriangularRandomVariable::~TriangularRandomVariable()
{ }


inline Real TriangularRandomVariable::cdf(Real x) const
{ return triangular_cdf(x, triangularMode, lowerBnd, upperBnd); }


inline Real TriangularRandomVariable::cdf_inverse(Real p) const
{
  // assume x < mode and then check
  Real range = upperBnd - lowerBnd,
       x = lowerBnd + std::sqrt(p*range*(triangularMode-lowerBnd));
  Real x_pdf = 2.*(x-lowerBnd)/range/(triangularMode-lowerBnd),
       m_pdf = 2./range;
  // check pdf value to ensure that proper equation used
  if ( x_pdf > m_pdf )
    x = upperBnd - std::sqrt(p*range*(upperBnd-triangularMode));
  return x;
}


inline Real TriangularRandomVariable::pdf(Real x) const
{ return triangular_pdf(x, triangularMode, lowerBnd, upperBnd); }


inline Real TriangularRandomVariable::pdf_gradient(Real x) const
{
  Real range = upperBnd - lowerBnd;
  if (x < triangularMode)
    return  2. / ( range * (triangularMode - lowerBnd) );
  else if (x > triangularMode)
    return -2. / ( range * (upperBnd - triangularMode) );
  else // x == triangularMode
    return  0.; // f'(x) is undefined: use 0.
}


inline Real TriangularRandomVariable::pdf_hessian(Real x) const
{ return 0.; }

} // namespace Pecos

#endif
