/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 LognormalRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef LOGNORMAL_RANDOM_VARIABLE_HPP
#define LOGNORMAL_RANDOM_VARIABLE_HPP

#include "NormalRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for Gaussian random variables.

/** Manages mean and standard deviation. */

class LognormalRandomVariable: public NormalRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  LognormalRandomVariable();                      ///< constructor
  LognormalRandomVariable(Real mean, Real stdev); ///< alternate constructor
  ~LognormalRandomVariable();                     ///< destructor

  //
  //- Heading: Member functions
  //

  Real cdf(Real x) const;
  Real cdf_inverse(Real p) const;

  Real pdf(Real x) const;
  //Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

protected:

  //
  //- Heading: Data
  //

  /// lambda parameter for lognormal random variable
  Real lambda;
  /// zeta parameter for lognormal random variable
  Real zeta;
};


inline LognormalRandomVariable::LognormalRandomVariable():
  NormalRandomVariable()
{ lognormal_params_from_moments(gaussMean, gaussStdDev, lambda, zeta); }


inline LognormalRandomVariable::LognormalRandomVariable(Real mean, Real stdev):
  NormalRandomVariable(mean, stdev)
{ lognormal_params_from_moments(mean, stdev, lambda, zeta); }


inline LognormalRandomVariable::~LognormalRandomVariable()
{ }


inline Real LognormalRandomVariable::cdf(Real x) const
{ return Phi((std::log(x) - lambda)/zeta); }


inline Real LognormalRandomVariable::cdf_inverse(Real p) const
{
  // p = Phi((std::log(x) - lambda)/zeta)
  return std::exp(Phi_inverse(p) * zeta + lambda);
}


inline Real LognormalRandomVariable::pdf(Real x) const
{
  lognormal_dist logn1(lambda, zeta);
  return bmth::pdf(logn1, x);
}


//inline Real LognormalRandomVariable::pdf_gradient(Real x) const
//{
//  return pdf(x);// * ...; // TO DO
//}


//inline Real LognormalRandomVariable::pdf_hessian(Real x) const
//{
//  return pdf(x);// * ...; // TO DO
//}

} // namespace Pecos

#endif
