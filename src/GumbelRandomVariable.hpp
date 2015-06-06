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

  Real correlation_warping_factor(const RandomVariable& rv, Real corr) const;

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

  /// pointer to the Boost weibull_distribution instance
  extreme_value_dist* gumbelDist;
};


inline GumbelRandomVariable::GumbelRandomVariable():
  RandomVariable(BaseConstructor()), alphaStat(1.), betaStat(0.),
  gumbelDist(NULL)
{ }


inline GumbelRandomVariable::GumbelRandomVariable(Real alpha, Real beta):
  RandomVariable(BaseConstructor()), alphaStat(alpha), betaStat(beta),
  gumbelDist(new extreme_value_dist(1./alpha, beta)) // note alpha inversion
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


inline void GumbelRandomVariable::update(Real alpha, Real beta)
{
  if (!gumbelDist || alphaStat != alpha || betaStat != beta) {
    alphaStat = alpha; betaStat = beta;
    if (gumbelDist) delete gumbelDist;
    gumbelDist = new extreme_value_dist(1./alpha, beta); // note alpha inversion
  }
}

// static functions:

inline Real GumbelRandomVariable::pdf(Real x, Real alpha, Real beta)
{
  extreme_value_dist gumbel1(1./alpha, beta);
  return bmth::pdf(gumbel1, x);

  //Real num = std::exp(-alpha*(x-beta));
  // if num overflows, apply l'Hopital's rule
  //return (num > DBL_MAX) ? 0. : alpha*num*std::exp(-num);
}


inline Real GumbelRandomVariable::cdf(Real x, Real alpha, Real beta)
{
  extreme_value_dist gumbel1(1./alpha, beta);
  return bmth::cdf(gumbel1, x);
  //return std::exp(-std::exp(-alpha * (x - beta)));
}


inline void GumbelRandomVariable::
moments_from_params(Real alpha, Real beta, Real& mean, Real& std_dev)
{ mean = beta + 0.57721566490153286/alpha; std_dev = PI/std::sqrt(6.)/alpha; }
/* Euler-Mascheroni constant is 0.5772... */

} // namespace Pecos

#endif
