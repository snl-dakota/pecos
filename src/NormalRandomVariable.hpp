/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 NormalRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef NORMAL_RANDOM_VARIABLE_HPP
#define NORMAL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for Gaussian random variables.

/** Manages mean and standard deviation. */

class NormalRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  NormalRandomVariable();                      ///< constructor
  NormalRandomVariable(Real mean, Real stdev); ///< alternate constructor
  ~NormalRandomVariable();                     ///< destructor

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

  Real log_pdf(Real x) const;

  Real to_std(Real x) const;
  Real from_std(Real z) const;

  Real correlation_warping_factor(const RandomVariable& rv, Real corr) const;

  //
  //- Heading: Member functions
  //

  void update(Real mean, Real stdev);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real mean, Real std_dev);
  static Real cdf(Real x, Real mean, Real std_dev);

  // standard normal distributions

  static Real std_pdf(Real beta);

  static Real std_cdf(Real beta);
  static Real std_ccdf(Real beta);
  static Real inverse_std_cdf(Real p_cdf);
  static Real inverse_std_ccdf(Real p_ccdf);

  //static Real log_std_cdf(Real beta);
  //static Real log_std_ccdf(Real beta);
  //static Real inverse_log_std_cdf(Real p_cdf);
  //static Real inverse_log_std_ccdf(Real p_ccdf);

  static Real phi(Real beta);
  static Real phi(Real beta, size_t n);
  static Real phi(const RealVector& u);

  static Real Phi(Real beta);
  static Real Phi_inverse(Real p_cdf);

protected:

  //
  //- Heading: Data
  //

  /// mean of Gaussian random variable
  Real gaussMean;
  /// standard deviation of Gaussian random variable
  Real gaussStdDev;

  // normal distribution instance from Boost math
  //normal_dist* normDist;
};


inline NormalRandomVariable::NormalRandomVariable():
  RandomVariable(BaseConstructor()), gaussMean(0), gaussStdDev(1.)
{ }


inline NormalRandomVariable::NormalRandomVariable(Real mean, Real stdev):
  RandomVariable(BaseConstructor()), gaussMean(mean), gaussStdDev(stdev)
{ }


inline NormalRandomVariable::~NormalRandomVariable()
{ }


inline Real NormalRandomVariable::cdf(Real x) const
{
  normal_dist norm(gaussMean, gaussStdDev);
  return bmth::cdf(norm, x);
}


inline Real NormalRandomVariable::ccdf(Real x) const
{
  normal_dist norm(gaussMean, gaussStdDev);
  return bmth::cdf(complement(norm, x)); //bmth::cdf(norm, -x);
}


inline Real NormalRandomVariable::inverse_cdf(Real p_cdf) const
{
  normal_dist norm(gaussMean, gaussStdDev);
  return bmth::quantile(norm, p_cdf); 
}


inline Real NormalRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  normal_dist norm(gaussMean, gaussStdDev);
  return bmth::quantile(complement(norm, p_ccdf));
}


inline Real NormalRandomVariable::pdf(Real x) const
{
  normal_dist norm(gaussMean, gaussStdDev);
  return bmth::pdf(norm, x);
}


inline Real NormalRandomVariable::pdf_gradient(Real x) const
{ return pdf(x) * (gaussMean - x) / (gaussStdDev * gaussStdDev); }


inline Real NormalRandomVariable::pdf_hessian(Real x) const
{
  Real var = gaussStdDev * gaussStdDev, mu_minus_x = gaussMean - x;
  return pdf(x) * ( mu_minus_x * mu_minus_x / var - 1. ) / var;
}


inline Real NormalRandomVariable::log_pdf(Real x) const
{
  Real xms = (x - gaussMean) / gaussStdDev;
  return -xms*xms / 2. - std::log(gaussStdDev * std::sqrt(2.*PI));
}


inline Real NormalRandomVariable::to_std(Real x) const
{ return (x - gaussMean) / gaussStdDev; }


inline Real NormalRandomVariable::from_std(Real z) const
{ return z * gaussStdDev + gaussMean; }


inline Real NormalRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // correlation warping factor for transformations to STD_NORMAL space
  // Der Kiureghian and Liu, ASCE JEM 112:1, 1986
  Real COV;
  switch (rv.type()) { // x-space types mapped to STD_NORMAL u-space

  case NORMAL:      return 1.; break; // No warping

  // Der Kiureghian & Liu: Table 2 (constants)
  case UNIFORM:     return 1.023326707946488488; break; // Max Error 0.0%
  case EXPONENTIAL: return 1.107; break;                // Max Error 0.0%
  case GUMBEL:      return 1.031; break;                // Max Error 0.0%

  // Der Kiureghian & Liu: Table 3 (quadratic approximations in COV)
  case LOGNORMAL:
    COV = rv.coefficient_of_variation();
    return COV/std::sqrt(bmth::log1p(COV*COV)); break; // Exact
  case GAMMA:
    COV = rv.coefficient_of_variation();
    return 1.001 + (-0.007 + 0.118*COV)*COV; break; // Max Error 0.0%
  case FRECHET:
    COV = rv.coefficient_of_variation();
    return 1.03  + ( 0.238 + 0.364*COV)*COV; break; // Max Error 0.1%
  case WEIBULL:
    COV = rv.coefficient_of_variation();
    return 1.031 + (-0.195 + 0.328*COV)*COV; break; // Max Error 0.1%

  default: // Unsupported warping (should be prevented upsteam)
    PCerr << "Error: unsupported correlation warping for NormalRV."<< std::endl;
    abort_handler(-1); return 1.; break;
  }
}


inline void NormalRandomVariable::update(Real mean, Real stdev)
{ gaussMean = mean; gaussStdDev = stdev; }

// Static member functions:

inline Real NormalRandomVariable::pdf(Real x, Real mean, Real std_dev)
{
  normal_dist norm(mean, std_dev);
  return bmth::pdf(norm, x);
  //return phi((x-mean)/std_dev)/std_dev;
}


inline Real NormalRandomVariable::cdf(Real x, Real mean, Real std_dev)
{
  normal_dist norm(mean, std_dev);
  return bmth::cdf(norm, x);
  //return Phi((x-mean)/std_dev);
}


/** returns a cumulative probability < 0.5 for negative beta and a
    cumulative probability > 0.5 for positive beta. */
inline Real NormalRandomVariable::std_pdf(Real beta)
{
  normal_dist norm(0., 1.);
  return bmth::pdf(norm, beta);
  //return std::exp(-beta*beta/2.)/std::sqrt(2.*PI);
}


/** returns a cumulative probability < 0.5 for negative beta and a
    cumulative probability > 0.5 for positive beta. */
inline Real NormalRandomVariable::std_cdf(Real beta)
{
  normal_dist norm(0., 1.);
  return bmth::cdf(norm, beta);
  //return .5 + .5*erf(beta/std::sqrt(2.));
}


inline Real NormalRandomVariable::inverse_std_cdf(Real p_cdf)
{
  normal_dist norm(0., 1.);
  return bmth::quantile(norm, p_cdf); 
  //return std::sqrt(2.)*erf_inverse(2.*p_cdf - 1.);
}


inline Real NormalRandomVariable::std_ccdf(Real beta)
{
  normal_dist norm(0., 1.);
  return bmth::cdf(complement(norm, beta)); //bmth::cdf(norm, -beta);
}


inline Real NormalRandomVariable::inverse_std_ccdf(Real p_ccdf)
{
  normal_dist norm(0., 1.);
  return bmth::quantile(complement(norm, p_ccdf));
}


/** log pdf would be useful for NormalRandomVariable, but log cdf offers
    no discernible benefit due to the presence of the error function. 

// avoid precision loss for large beta > 0 (cdf indistinguishable from 1)
inline Real NormalRandomVariable::log_std_cdf(Real beta)
{ return (beta > 0.) ? bmth::log1p(-std_cdf(-beta)) : std::log(std_cdf(beta)); }


inline Real NormalRandomVariable::inverse_log_std_cdf(Real log_p_cdf)
{
  normal_dist norm(0., 1.);
  // check p > .5 --> ln p > -.7
  return (log_p_cdf > -.7) ?
    bmth::quantile(complement(norm, -bmth::expm1(log_p_cdf))) :
    bmth::quantile(norm, std::exp(log_p_cdf));
}


// avoid precision loss for large beta < 0 (ccdf indistinguishable from 1)
inline Real NormalRandomVariable::log_std_ccdf(Real beta)
{
  return (beta < 0.) ?
    bmth::log1p(-std_ccdf(-beta)) : std::log(std_ccdf(beta));
}


inline Real NormalRandomVariable::inverse_log_std_ccdf(Real log_p_ccdf)
{
  normal_dist norm(0., 1.);
  return (log_p_ccdf > -.7) ?
    bmth::quantile(norm, -bmth::expm1(log_p_ccdf)) :
    bmth::quantile(complement(norm, std::exp(log_p_ccdf)));
}
*/


/** Univariate standard normal density function. */
inline Real NormalRandomVariable::phi(Real beta)
{ return std_pdf(beta); }


/** Multivariate standard normal density function with aggregate distance. */
inline Real NormalRandomVariable::phi(Real beta, size_t n)
{
  // need n instances of 1/sqrt(2Pi), but 1D pdf only includes 1:
  return (n > 1) ?
    std_pdf(beta) * std::pow(2.*PI, -((Real)(n-1))/2.) :// correct 1D pdf for nD
    std_pdf(beta);
}


/** Multivariate standard normal density function from vector. */
inline Real NormalRandomVariable::phi(const RealVector& u)
{
  return phi(u.normFrobenius(), u.length());

  // Alternate implementation invokes exp() repeatedly:
  //normal_dist norm(0., 1.);
  //size_t i, n = u.length(); Real pdf = 1.;
  //for (i=0; i<n; ++i)
  //  pdf *= bmth::pdf(norm, u[i]);
}


/** returns a cumulative probability < 0.5 for negative beta and a
    cumulative probability > 0.5 for positive beta. */
inline Real NormalRandomVariable::Phi(Real beta)
{ return std_cdf(beta); }

/** returns a negative beta for cumulative probability < 0.5 and a
    positive beta for cumulative probability > 0.5. */
inline Real NormalRandomVariable::Phi_inverse(Real p_cdf)
{ return inverse_std_cdf(p_cdf); }

} // namespace Pecos

#endif
