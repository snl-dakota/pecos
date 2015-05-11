/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 GammaRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef GAMMA_RANDOM_VARIABLE_HPP
#define GAMMA_RANDOM_VARIABLE_HPP

#include "ExponentialRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for gamma random variables.

/** Manages alpha and inherits beta. */

class GammaRandomVariable: public ExponentialRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  GammaRandomVariable();
  /// alternate constructor
  GammaRandomVariable(Real alpha, Real beta);
  /// destructor
  ~GammaRandomVariable();

  //
  //- Heading: Member functions
  //

  Real cdf(Real x) const;
  Real cdf_inverse(Real p) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

protected:

  //
  //- Heading: Data
  //

  /// alpha parameter of gamma random variable (statistical PDF
  /// convention; differs from GenLaguerre polynomial convention)
  Real alphaStat;
};


inline GammaRandomVariable::GammaRandomVariable():
  ExponentialRandomVariable(), alphaStat(0.) // TO DO: same as Exponential?
{ }


inline GammaRandomVariable::GammaRandomVariable(Real alpha, Real beta):
  ExponentialRandomVariable(beta), alphaStat(alpha)
{ }


inline GammaRandomVariable::~GammaRandomVariable()
{ }


inline Real GammaRandomVariable::cdf(Real x) const
{ return gamma_cdf(x, alphaStat, betaStat); }


inline Real GammaRandomVariable::cdf_inverse(Real p) const
{ return gamma_cdf_inverse(p, alphaStat, betaStat); }


inline Real GammaRandomVariable::pdf(Real x) const
{ return gamma_pdf(x, alphaStat, betaStat); }


inline Real GammaRandomVariable::pdf_gradient(Real x) const
{
  return std::pow(betaStat,-alphaStat) / gamma_function(alphaStat) *
    (std::exp(-x/betaStat)*(alphaStat-1.) * std::pow(x,alphaStat-2.) -
     std::pow(x,alphaStat-1.) * std::exp(-x/betaStat)/betaStat);
}


//inline Real GammaRandomVariable::pdf_hessian(Real x) const
//{
//  return gamma_pdf(x, alphaStat, betaStat); // * ...; // TO DO
//}

} // namespace Pecos

#endif
