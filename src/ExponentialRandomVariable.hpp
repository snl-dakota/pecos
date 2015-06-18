/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 ExponentialRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef EXPONENTIAL_RANDOM_VARIABLE_HPP
#define EXPONENTIAL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"
#include "NormalRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for exponential random variables.

/** Manages beta parameter.  Pecos employs the 1/beta exp(-x/beta)
    definition, which differs from the lambda exp(-lambda x) LHS
    definition (lambda_LHS = 1/beta_PECOS). */

class ExponentialRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  ExponentialRandomVariable();
  /// alternate constructor
  ExponentialRandomVariable(Real beta);
  /// destructor
  ~ExponentialRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;

  Real inverse_log_ccdf(Real log_p_ccdf) const;
  Real log_pdf(Real x) const;

  Real to_std(Real x) const;
  Real from_std(Real z) const;

  Real parameter(short dist_param) const;
  void parameter(short dist_param, Real val);

  RealRealPair moments() const;
  RealRealPair bounds() const;

  Real coefficient_of_variation() const;
  Real correlation_warping_factor(const RandomVariable& rv, Real corr) const;

  Real dx_ds(short dist_param, short u_type, Real x, Real z) const;
  Real dz_ds_factor(short u_type, Real x, Real z) const;

  //
  //- Heading: Member functions
  //

  void update(Real beta);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real std_pdf(Real x);
  static Real std_cdf(Real x);

  static Real pdf(Real x, Real beta);
  static Real cdf(Real x, Real beta);

  static void moments_from_params(Real beta, Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// beta parameter of exponential random variable
  Real betaStat;
};


inline ExponentialRandomVariable::ExponentialRandomVariable():
  RandomVariable(BaseConstructor()), betaStat(1.)
{ }


inline ExponentialRandomVariable::ExponentialRandomVariable(Real beta):
  RandomVariable(BaseConstructor()), betaStat(beta)
{ }


inline ExponentialRandomVariable::~ExponentialRandomVariable()
{ }


inline Real ExponentialRandomVariable::cdf(Real x) const
{ return -bmth::expm1(-x/betaStat); }


inline Real ExponentialRandomVariable::ccdf(Real x) const
{ return std::exp(-x/betaStat); }


inline Real ExponentialRandomVariable::inverse_cdf(Real p_cdf) const
{
  // p_cdf = 1 - exp(-x/beta)  -->  -x/beta = log(1-p_cdf)
  return -betaStat * bmth::log1p(-p_cdf);
}


inline Real ExponentialRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return -betaStat * std::log(p_ccdf); }


inline Real ExponentialRandomVariable::inverse_log_ccdf(Real log_p_ccdf) const
{ return -betaStat * log_p_ccdf; }


//  F(x) = 1. - e^(-x/beta)
//  f(x) = e^(-x/beta) / beta
// f'(x) = - e^(-x/beta) / beta^2
inline Real ExponentialRandomVariable::pdf(Real x) const
{ return std::exp(-x/betaStat)/betaStat; }


inline Real ExponentialRandomVariable::pdf_gradient(Real x) const
{ return -pdf(x, betaStat) / betaStat; }


inline Real ExponentialRandomVariable::pdf_hessian(Real x) const
{ return pdf(x, betaStat) / (betaStat * betaStat); }


inline Real ExponentialRandomVariable::log_pdf(Real x) const
{ return -x / betaStat - std::log(betaStat); }


inline Real ExponentialRandomVariable::to_std(Real x) const
{ return x / betaStat; }


inline Real ExponentialRandomVariable::from_std(Real z) const
{ return z * betaStat; }


inline Real ExponentialRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case E_BETA: return betaStat; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in ExponentialRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void ExponentialRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case E_BETA: betaStat = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in ExponentialRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
}


inline RealRealPair ExponentialRandomVariable::moments() const
{ return RealRealPair(betaStat, betaStat); }


inline RealRealPair ExponentialRandomVariable::bounds() const
{ return RealRealPair(0., std::numeric_limits<Real>::infinity()); }


inline Real ExponentialRandomVariable::coefficient_of_variation() const
{ return 1.; }


inline Real ExponentialRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // correlation warping factor for transformations to STD_NORMAL space
  // Der Kiureghian and Liu, ASCE JEM 112:1, 1986
  Real COV;
  switch (rv.type()) { // x-space types mapped to STD_NORMAL u-space

  // Der Kiureghian & Liu: Table 4 (quadratic approximations in corr)
  case EXPONENTIAL:
    return 1.229 + (-0.367 + 0.153*corr)*corr;      break; // Max Error 1.5%
  case GUMBEL:
    return 1.142 + (-0.154*corr + 0.031*corr)*corr; break; // Max Error 0.2%

  // Der Kiureghian & Liu: Table 5 (quadratic approximations in corr,COV)
  case GAMMA:
    COV = rv.coefficient_of_variation();
    return 1.104 + (0.003 + 0.014*corr)*corr
      + (-0.008 + 0.173*COV - 0.296*corr)*COV; break; // Max Error 0.9%
  case FRECHET:
    COV = rv.coefficient_of_variation();
    return 1.109 + (-0.152 + 0.130*corr)*corr
      + ( 0.361 + 0.455*COV - 0.728*corr)*COV; break; // Max Error 4.5%
  case WEIBULL:
    COV = rv.coefficient_of_variation();
    return 1.147 + (0.145 + 0.010*corr)*corr
      + (-0.271 + 0.459*COV - 0.467*corr)*COV; break; // Max Error 0.4%

  // warping factors are defined once for lower triangle based on uv order
  case NORMAL: case LOGNORMAL: case UNIFORM:
    return rv.correlation_warping_factor(*this, corr); break;

  default: // Unsupported warping (should be prevented upsteam)
    PCerr << "Error: unsupported correlation warping for ExponentialRV."
	  << std::endl;
    abort_handler(-1); return 1.; break;
  }
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  dz/ds is zero if uncorrelated, 
    while dz_ds_factor() manages contributions in the correlated case. */
inline Real ExponentialRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  // to STD_EXPONENTIAL: x = beta*z
  // to STD_NORMAL:      Phi(z) = 1. - exp(-x/beta)
  //                     x = -beta ln(1. - Phi(z))
  bool u_type_err = false, dist_err = false;
  switch (dist_param) {
  case E_BETA: // Deriv of exponential w.r.t. beta
    switch (u_type) {
    case STD_NORMAL:      return x / betaStat; break;
    //case STD_UNIFORM:   TO DO;               break;
    case STD_EXPONENTIAL: return z;            break;
    default:              u_type_err = true;   break;
    }
    break;
  default:                dist_err = true;     break;
  }

  if (u_type_err)
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in ExponentialRandomVariable::dx_ds()." << std::endl;
  if (dist_err)
    PCerr << "Error: mapping failure for distribution parameter " << dist_param
	  << " in ExponentialRandomVariable::dx_ds()." << std::endl;
  if (u_type_err || dist_err)
    abort_handler(-1);
  return 0.;
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  For the uncorrelated case,
    u and z are constants.  For the correlated case, u is a constant, but 
    z(s) = L(s) u due to Nataf dependence on s and dz/ds = dL/ds u. */
inline Real ExponentialRandomVariable::
dz_ds_factor(short u_type, Real x, Real z) const
{
  switch (u_type) {
  case STD_NORMAL:
    return betaStat * NormalRandomVariable::std_pdf(z) /
      NormalRandomVariable::std_ccdf(z);                    break;
  //case STD_UNIFORM:   TO DO;                              break;
  case STD_EXPONENTIAL: return betaStat;                    break;
  default:
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in ExponentialRandomVariable::dz_ds_factor()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void ExponentialRandomVariable::update(Real beta)
{ betaStat = beta; }

// static functions:

inline Real ExponentialRandomVariable::std_pdf(Real x)
{ return std::exp(-x); }


inline Real ExponentialRandomVariable::std_cdf(Real x)
{
  // as with log1p(), avoid numerical probs when exp(~0) is ~ 1
  return -bmth::expm1(-x); //1. - std::exp(-x);
}


inline Real ExponentialRandomVariable::pdf(Real x, Real beta)
{ return std::exp(-x/beta)/beta; }


inline Real ExponentialRandomVariable::cdf(Real x, Real beta)
{
  // as with log1p(), avoid numerical probs when exp(~0) is ~ 1
  return -bmth::expm1(-x/beta);
}


inline void ExponentialRandomVariable::
moments_from_params(Real beta, Real& mean, Real& std_dev)
{ mean = std_dev = beta; }

} // namespace Pecos

#endif