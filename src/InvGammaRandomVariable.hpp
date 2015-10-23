/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 InvGammaRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef INV_GAMMA_RANDOM_VARIABLE_HPP
#define INV_GAMMA_RANDOM_VARIABLE_HPP

#include "ExponentialRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for gamma random variables.

/** Manages alpha and inherits beta.  This follows the
    Gamma(alpha,beta) definition in Law & Kelton, which differs from
    the LHS definition (beta_LK = 1/beta_LHS). */

class InvGammaRandomVariable: public ExponentialRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  InvGammaRandomVariable();
  /// alternate constructor
  InvGammaRandomVariable(Real alpha, Real beta);
  /// destructor
  ~InvGammaRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  /*
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;
  Real log_pdf(Real x) const;
  Real log_pdf_gradient(Real x) const;
  Real log_pdf_hessian(Real x) const;

  Real standard_pdf(Real z) const;
  Real log_standard_pdf(Real z) const;
  Real log_standard_pdf_gradient(Real z) const;
  Real log_standard_pdf_hessian(Real z) const;
  */

  // inherited from ExponentialRandomVariable
  //Real to_standard(Real x) const;
  //Real from_standard(Real z) const;

  Real parameter(short dist_param) const;
  void parameter(short dist_param, Real val);

  RealRealPair moments() const;

  Real coefficient_of_variation() const;

  //
  //- Heading: Member functions
  //

  void update(Real alpha, Real beta);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real alpha, Real beta);
  static Real cdf(Real x, Real alpha, Real beta);
  static Real inverse_cdf(Real cdf, Real alpha, Real beta);

  static void moments_from_params(Real alpha, Real beta,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Member functions
  //

  /// create a new invGammaDist instance
  void update_boost();

  //
  //- Heading: Data
  //

  /// alpha shape parameter of inverse gamma random variable (statistical
  /// PDF convention; differs from GenLaguerre polynomial convention)
  Real alphaShape;

  /// pointer to the Boost inv_gamma_dist instance
  inv_gamma_dist* invGammaDist;
};


inline InvGammaRandomVariable::InvGammaRandomVariable():
  ExponentialRandomVariable(), alphaShape(3.), invGammaDist(NULL)
{ ranVarType = INV_GAMMA; }


inline InvGammaRandomVariable::InvGammaRandomVariable(Real alpha, Real beta):
  ExponentialRandomVariable(beta), alphaShape(alpha),
  invGammaDist(new inv_gamma_dist(alphaShape, betaScale))
{ ranVarType = INV_GAMMA; }


inline InvGammaRandomVariable::~InvGammaRandomVariable()
{ if (invGammaDist) delete invGammaDist; }


inline Real InvGammaRandomVariable::cdf(Real x) const
{ return bmth::cdf(*invGammaDist, x); }


inline Real InvGammaRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*invGammaDist, x)); }


inline Real InvGammaRandomVariable::inverse_cdf(Real p_cdf) const
{ return bmth::quantile(*invGammaDist, p_cdf); }


inline Real InvGammaRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*invGammaDist, p_ccdf)); }


inline Real InvGammaRandomVariable::pdf(Real x) const
{ return bmth::pdf(*invGammaDist, x); }


/*
inline Real InvGammaRandomVariable::pdf_gradient(Real x) const
{
  if (x <= 0.) {
    if      (alphaShape < 1.) return -std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return  std::numeric_limits<Real>::quiet_NaN();
    else return -ExponentialRandomVariable::pdf(x) / betaScale;
  }
  else
    return pdf(x) * ((alphaShape-1.)/x - 1./betaScale);
}


inline Real InvGammaRandomVariable::pdf_hessian(Real x) const
{
  if (x <= 0.) {
    if      (alphaShape < 1.) return std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return std::numeric_limits<Real>::quiet_NaN();
    else return ExponentialRandomVariable::pdf(x) / (betaScale*betaScale);
  }
  else {
    Real am1 = alphaShape - 1., term = am1 / x - 1. / betaScale;
    return pdf(x) * (term*term - am1 / (x*x));
  }
}


inline Real InvGammaRandomVariable::log_pdf(Real x) const
{
  if (x <= 0.) {
    if      (alphaShape < 1.) return  std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return -std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_pdf(x);
  }
  else
    return (alphaShape-1.)*std::log(x) - x/betaScale
      - std::log(bmth::tgamma(alphaShape)) - alphaShape*std::log(betaScale);
}


inline Real InvGammaRandomVariable::log_pdf_gradient(Real x) const
{
  if (x <= 0.) {
    if      (alphaShape < 1.) return -std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return  std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_pdf_gradient(x);
  }
  else
    return (alphaShape-1.)/x - 1./betaScale;
}


inline Real InvGammaRandomVariable::log_pdf_hessian(Real x) const
{
  if (x <= 0.) {
    if      (alphaShape < 1.) return  std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return -std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_pdf_hessian(x);
  }
  else
    return (1.-alphaShape)/(x*x);
}


inline Real InvGammaRandomVariable::standard_pdf(Real z) const
{
  inv_gamma_dist inv_gamma1(alphaShape, 1.);
  return bmth::pdf(inv_gamma1, z);
}


inline Real InvGammaRandomVariable::log_standard_pdf(Real z) const
{
  if (z <= 0.) {
    if      (alphaShape < 1.) return  std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return -std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_standard_pdf(z);
  }
  else
    return (alphaShape-1.)*std::log(z) - z - std::log(bmth::tgamma(alphaShape)); 
}


inline Real InvGammaRandomVariable::log_standard_pdf_gradient(Real z) const
{
  if (z <= 0.) {
    if      (alphaShape < 1.) return -std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return  std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_standard_pdf_gradient(z);
  }
  else
    return (alphaShape-1.)/z - 1.;
}


inline Real InvGammaRandomVariable::log_standard_pdf_hessian(Real z) const
{
  if (z <= 0.) {
    if      (alphaShape < 1.) return  std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return -std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_standard_pdf_hessian(z);
  }
  else
    return (1.-alphaShape)/(z*z);
}
*/

inline Real InvGammaRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case IGA_ALPHA: return alphaShape; break;
  case IGA_BETA:  return betaScale;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in InvGammaRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void InvGammaRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case IGA_ALPHA: alphaShape = val; break;
  case IGA_BETA:  betaScale  = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in InvGammaRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
  update_boost(); // create a new invGammaDist instance
}


inline RealRealPair InvGammaRandomVariable::moments() const
{
  Real mean, std_dev;
  moments_from_params(alphaShape, betaScale, mean, std_dev);
  return RealRealPair(mean, std_dev);
}


inline Real InvGammaRandomVariable::coefficient_of_variation() const
{
  Real mean, std_dev;
  moments_from_params(alphaShape, betaScale, mean, std_dev);
  return std_dev / mean;
}


inline void InvGammaRandomVariable::update_boost()
{
  if (invGammaDist) delete invGammaDist;
  invGammaDist = new inv_gamma_dist(alphaShape, betaScale);
}


inline void InvGammaRandomVariable::update(Real alpha, Real beta)
{
  if (!invGammaDist || alphaShape != alpha || betaScale != beta)
    { alphaShape = alpha; betaScale = beta; update_boost(); }
}


inline Real InvGammaRandomVariable::pdf(Real x, Real alpha, Real beta)
{
  inv_gamma_dist inv_gamma1(alpha, beta);
  return bmth::pdf(inv_gamma1, x);
}


inline Real InvGammaRandomVariable::cdf(Real x, Real alpha, Real beta)
{
  inv_gamma_dist inv_gamma1(alpha, beta);
  return bmth::cdf(inv_gamma1, x);
}


inline Real InvGammaRandomVariable::inverse_cdf(Real cdf, Real alpha, Real beta)
{
  inv_gamma_dist inv_gamma1(alpha, beta);
  return bmth::quantile(inv_gamma1, cdf);
}


inline void InvGammaRandomVariable::
moments_from_params(Real alpha, Real beta, Real& mean, Real& std_dev)
{
  inv_gamma_dist inv_gamma1(alpha, beta);
  mean    = bmth::mean(inv_gamma1);                // domain_error if alpha <= 1
  std_dev = std::sqrt(bmth::variance(inv_gamma1)); // domain_error if alpha <= 2
}

} // namespace Pecos

#endif
