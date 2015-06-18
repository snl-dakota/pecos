/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 GumbelRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef GUMBEL_RANDOM_VARIABLE_HPP
#define GUMBEL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for gumbel random variables.

/** Manages alpha and beta parameters.  See Haldar and Mahadevan, p. 90. */

class GumbelRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Virtual function redefinitions
  //

  /// default constructor
  GumbelRandomVariable();
  /// alternate constructor
  GumbelRandomVariable(Real alpha, Real beta);
  /// destructor
  ~GumbelRandomVariable();

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

  Real inverse_log_cdf(Real log_p) const;
  Real log_pdf(Real x) const;

  Real parameter(short dist_param) const;
  void parameter(short dist_param, Real val);

  RealRealPair moments() const;
  RealRealPair bounds() const;

  Real correlation_warping_factor(const RandomVariable& rv, Real corr) const;

  Real dx_ds(short dist_param, short u_type, Real x, Real z) const;
  Real dz_ds_factor(short u_type, Real x, Real z) const;

  //
  //- Heading: Member functions
  //

  void update(Real alpha, Real beta);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real alpha, Real beta);
  static Real cdf(Real x, Real alpha, Real beta);

  static void moments_from_params(Real alpha, Real beta,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// alpha parameter of gumbel random variable (inverse of scale)
  Real alphaStat;
  /// beta parameter of gumbel random variable (location)
  Real betaStat;

  /// overflow and underflow values for argument of std::exp
  RealRealPair expLimits;

  /// pointer to the Boost weibull_distribution instance
  extreme_value_dist* gumbelDist;
};


inline GumbelRandomVariable::GumbelRandomVariable():
  RandomVariable(BaseConstructor()), alphaStat(1.), betaStat(0.),
  expLimits(std::log(DBL_MAX),std::log(DBL_MIN)), gumbelDist(NULL)
{ }


inline GumbelRandomVariable::GumbelRandomVariable(Real alpha, Real beta):
  RandomVariable(BaseConstructor()), alphaStat(alpha), betaStat(beta),
  expLimits(std::log(DBL_MAX),std::log(DBL_MIN)),
  gumbelDist(new extreme_value_dist(beta, 1./alpha)) // location, scale
{ }


inline GumbelRandomVariable::~GumbelRandomVariable()
{ }


inline Real GumbelRandomVariable::cdf(Real x) const
{
  return bmth::cdf(*gumbelDist, x);
  //return std::exp(-std::exp(-alphaStat * (x - betaStat)));
}


inline Real GumbelRandomVariable::ccdf(Real x) const
{
  return bmth::cdf(complement(*gumbelDist, x));
  //return -bmth::expm1(-std::exp(-alphaStat * (x - betaStat)));
}


inline Real GumbelRandomVariable::inverse_cdf(Real p_cdf) const
{
  return bmth::quantile(*gumbelDist, p_cdf);

  // p = std::exp(-std::exp(-alpha * (x - beta)))
  //return betaStat - std::log(-std::log(p_cdf)) / alphaStat;
}


inline Real GumbelRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  return bmth::quantile(complement(*gumbelDist, p_ccdf));
  //return betaStat - std::log(-bmth::log1p(-p_ccdf)) / alphaStat;
}


inline Real GumbelRandomVariable::inverse_log_cdf(Real log_p) const
{ return betaStat - std::log(-log_p) / alphaStat; }


//  F(x) = e^(-e^(-alpha*(x-u)))
//  f(x) = alpha e^(-alpha*(x-u)) F(x)
// f'(x) = alpha (e^(-alpha*(x-u)) f(x) - alpha F(x) e^(-alpha*(x-u)))
inline Real GumbelRandomVariable::pdf(Real x) const
{ return bmth::pdf(*gumbelDist, x); }


inline Real GumbelRandomVariable::pdf_gradient(Real x) const
{
  Real num = std::exp(-alphaStat*(x-betaStat)),
    cdf = std::exp(-num), pdf = alphaStat * num * cdf;
  return alphaStat * num * (pdf - alphaStat * cdf);
}


//inline Real GumbelRandomVariable::pdf_hessian(Real x) const
//{
//  return pdf(x, alphaStat, betaStat) * ...; // TO DO
//}


inline Real GumbelRandomVariable::log_pdf(Real x) const
{
  Real num = -alphaStat*(x-betaStat);
  return std::log(alphaStat) + num - std::exp(num);
}


inline Real GumbelRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case GU_ALPHA: return alphaStat; break;
  case GU_BETA:  return betaStat;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in GumbelRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void GumbelRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case GU_ALPHA: alphaStat = val; break;
  case GU_BETA:  betaStat  = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in GumbelRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
}


inline RealRealPair GumbelRandomVariable::moments() const
{
  Real mean, std_dev;
  moments_from_params(alphaStat, betaStat, mean, std_dev);
  return RealRealPair(mean, std_dev);
}


inline RealRealPair GumbelRandomVariable::bounds() const
{
  Real dbl_inf = std::numeric_limits<Real>::infinity();
  return RealRealPair(-dbl_inf, dbl_inf);
}


inline Real GumbelRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // correlation warping factor for transformations to STD_NORMAL space
  // Der Kiureghian and Liu, ASCE JEM 112:1, 1986
  Real COV;
  switch (rv.type()) { // x-space types mapped to STD_NORMAL u-space

  // Der Kiureghian & Liu: Table 4 (quadratic approx in corr)
  case GUMBEL:      // Max Error 0.0%
    return 1.064 + (-0.069 + 0.005*corr)*corr; break;

  // Der Kiureghian & Liu: Table 5 (quadratic approx in COV,corr)
  case FRECHET: // Max Error 1.0%
    COV = rv.coefficient_of_variation();
    return 1.056 + (-0.06 + 0.02*corr)*corr
      + (0.263 + 0.383*COV - 0.332*corr)*COV; break;
  case WEIBULL: // Max Error 0.2%
    COV = rv.coefficient_of_variation();
    return 1.064 + (0.065 +  0.003*corr)*corr
      + (-0.21 + 0.356*COV - 0.211*corr)*COV; break;

  // warping factors are defined once for lower triangle based on uv order
  case NORMAL: case LOGNORMAL: case UNIFORM: case EXPONENTIAL: case GAMMA:
    return rv.correlation_warping_factor(*this, corr); break;

  default: // Unsupported warping (should be prevented upsteam)
    PCerr << "Error: unsupported correlation warping for GumbelRV."<< std::endl;
    abort_handler(-1); return 1.; break;
  }
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  dz/ds is zero if uncorrelated, 
    while dz_ds_factor() manages contributions in the correlated case. */
inline Real GumbelRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  // to STD_NORMAL: x = beta - ln(-ln(Phi(z)))/alpha
  bool u_type_err = false, dist_err = false;
  switch (u_type) {
  case STD_NORMAL:
    switch (dist_param) {
    case GU_ALPHA: return (betaStat - x) / alphaStat; break; // ln(-ln Phi)/a^2
    case GU_BETA:  return 1.;                         break;
    // Gumbel Mean:
    //   num = -alpha*(x-z); return -alpha*exp(num-exp(num))/phi(z); break;
    // Gumbel Standard Deviation:
    //   num = -alpha*(x-z); return num*exp(num-exp(num))/stdev/phi(z); break;
    default:      dist_err = true; break;
    }
    break;
  //case GUMBEL:  TO DO;                                        break;
  default:        u_type_err = true;                            break;
  }

  if (u_type_err)
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in GumbelRandomVariable::dx_ds()." << std::endl;
  if (dist_err)
    PCerr << "Error: mapping failure for distribution parameter " << dist_param
	  << " in GumbelRandomVariable::dx_ds()." << std::endl;
  if (u_type_err || dist_err)
    abort_handler(-1);
  return 0.;
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  For the uncorrelated case,
    u and z are constants.  For the correlated case, u is a constant, but 
    z(s) = L(s) u due to Nataf dependence on s and dz/ds = dL/ds u. */
inline Real GumbelRandomVariable::
dz_ds_factor(short u_type, Real x, Real z) const
{
  switch (u_type) {
  case STD_NORMAL:
    return -NormalRandomVariable::std_pdf(z) / (alphaStat *
      NormalRandomVariable::std_cdf(z) * NormalRandomVariable::log_std_cdf(z));
    break;
  //case GUMBEL: TO DO; break;
  default:
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in GumbelRandomVariable::dz_ds_factor()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void GumbelRandomVariable::update(Real alpha, Real beta)
{
  if (!gumbelDist || alphaStat != alpha || betaStat != beta) {
    alphaStat = alpha; betaStat = beta;
    if (gumbelDist) delete gumbelDist;
    gumbelDist = new extreme_value_dist(beta, 1./alpha); // location, scale
  }
}

// static functions:

inline Real GumbelRandomVariable::pdf(Real x, Real alpha, Real beta)
{
  /*
  // Unprotected implementation observed to generate nans
  extreme_value_dist gumbel1(beta, 1./alpha); // location, scale
  Real pdf1 = bmth::pdf(gumbel1, x);

  Real num = std::exp(alpha*(beta-x));
  // if num overflows, apply l'Hopital's rule
  Real pdf2 = (num > DBL_MAX) ? 0. : alpha*num*std::exp(-num);

  //  std::log(DBL_MAX) = 709.783 std::log(DBL_MIN) = -708.396

  PCout << "x = " << x << " num = " << num << " Boost pdf = " << pdf1
	<< " Pecos pdf = " << pdf2 << " abs diff = " << std::abs(pdf1-pdf2)
	<< '\n';
    // yields "num = inf Boost EVD pdf = -nan explicit pdf = 0 abs diff = nan"
  // for small x

  return pdf2;
  */
  ///////////////////

  Real num = alpha*(beta-x);
  if (num > 700.) // 1st exp overflows; log(DBL_MAX) = 709.783
    return 0.; // Boost generates nan; use l'Hopital's rule: denominator exp(-num) grows faster
  //else if (num < -700.)) // underflows; log(DBL_MIN) = -708.396
  //  return 0.; // Boost handles this case OK
  else {
    extreme_value_dist gumbel1(beta, 1./alpha); // location, scale
    return bmth::pdf(gumbel1, x);
  }
}


inline Real GumbelRandomVariable::cdf(Real x, Real alpha, Real beta)
{
  extreme_value_dist gumbel1(beta, 1./alpha); // location, scale
  return bmth::cdf(gumbel1, x);
  //return std::exp(-std::exp(alpha*(beta-x)));
}


inline void GumbelRandomVariable::
moments_from_params(Real alpha, Real beta, Real& mean, Real& std_dev)
{ mean = beta + 0.57721566490153286/alpha; std_dev = PI/std::sqrt(6.)/alpha; }
/* Euler-Mascheroni constant is 0.5772... */

} // namespace Pecos

#endif
