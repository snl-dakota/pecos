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

/** Manages alpha and inherits beta.  This follows the
    Gamma(alpha,beta) definition in Law & Kelton, which differs from
    the LHS definition (beta_LK = 1/beta_LHS). */

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
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

  // inherited from ExponentialRandomVariable
  //Real to_std(Real x) const;
  //Real from_std(Real z) const;

  void update(Real alpha, Real beta);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real alpha, Real beta);
  static Real cdf(Real x, Real alpha, Real beta);
  static Real inverse_cdf(Real cdf, Real alpha, Real beta);

  static void moments_from__params(Real alpha, Real beta,
				   Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// alpha parameter of gamma random variable (statistical PDF
  /// convention; differs from GenLaguerre polynomial convention)
  Real alphaStat;

  /// pointer to the Boost gamma_distribution instance
  gamma_dist* gammaDist;
};


inline GammaRandomVariable::GammaRandomVariable():
  ExponentialRandomVariable(), alphaStat(1.), gammaDist(NULL)
{ }


inline GammaRandomVariable::GammaRandomVariable(Real alpha, Real beta):
  ExponentialRandomVariable(beta), alphaStat(alpha),
  gammaDist(new gamma_dist(alphaStat, betaStat))
{ }


inline GammaRandomVariable::~GammaRandomVariable()
{ if (gammaDist) delete gammaDist; }


inline Real GammaRandomVariable::cdf(Real x) const
{ return bmth::cdf(*gammaDist, x); }


inline Real GammaRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*gammaDist, x)); }


inline Real GammaRandomVariable::inverse_cdf(Real p_cdf) const
{ return bmth::quantile(*gammaDist, p_cdf); }


inline Real GammaRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*gammaDist, p_ccdf)); }


inline Real GammaRandomVariable::pdf(Real x) const
{ return bmth::pdf(*gammaDist, x); }


inline Real GammaRandomVariable::pdf_gradient(Real x) const
{
  return std::pow(betaStat,-alphaStat) / bmth::tgamma(alphaStat) *
    (std::exp(-x/betaStat)*(alphaStat-1.) * std::pow(x,alphaStat-2.) -
     std::pow(x,alphaStat-1.) * std::exp(-x/betaStat)/betaStat);
}


//inline Real GammaRandomVariable::pdf_hessian(Real x) const
//{
//  return gamma_pdf(x, alphaStat, betaStat); // * ...; // TO DO
//}


inline void GammaRandomVariable::update(Real alpha, Real beta)
{
  if (!gammaDist || alphaStat != alpha || betaStat != beta) {
    alphaStat = alpha; betaStat = beta;
    if (gammaDist) delete gammaDist;
    gammaDist = new gamma_dist(alphaStat, betaStat);
  }
}


inline Real GammaRandomVariable::pdf(Real x, Real alpha, Real beta)
{
  gamma_dist gamma1(alpha, beta);
  return bmth::pdf(gamma1, x);
}


inline Real GammaRandomVariable::cdf(Real x, Real alpha, Real beta)
{
  gamma_dist gamma1(alpha, beta);
  return bmth::cdf(gamma1, x);
}


inline Real GammaRandomVariable::inverse_cdf(Real cdf, Real alpha, Real beta)
{
  gamma_dist gamma1(alpha, beta);
  return bmth::quantile(gamma1, cdf);
}


inline void GammaRandomVariable::
moments_from__params(Real alpha, Real beta, Real& mean, Real& std_dev)
{ mean = alpha*beta; std_dev = std::sqrt(alpha)*beta; }

} // namespace Pecos

#endif
