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
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;

  Real to_std(Real x) const;
  Real from_std(Real z) const;

  Real correlation_warping_factor(const RandomVariable& rv, Real corr) const;

  //
  //- Heading: Member functions
  //

  void update(Real lwr, Real upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real lwr, Real upr);
  static Real cdf(Real x, Real lwr, Real upr);
  static Real inverse_cdf(Real p_cdf, Real lwr, Real upr);

  static Real std_pdf();

  static Real std_cdf(Real beta);
  static Real std_ccdf(Real beta);
  static Real inverse_std_cdf(Real p_cdf);
  static Real inverse_std_ccdf(Real p_ccdf);

  static void moments_from_params(Real lwr, Real upr, Real& mean,
				  Real& std_dev);

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


inline Real UniformRandomVariable::std_pdf()
{ return 0.5; } // equal probability on [-1,1]


inline Real UniformRandomVariable::std_cdf(Real x)
{
  if      (x >=  1.) return 1.;
  else if (x <= -1.) return 0.;
  else               return (x + 1.)/2.; // linear x \in [-1,1] -> p \in [0,1]
}


inline Real UniformRandomVariable::std_ccdf(Real x)
{
  if      (x >=  1.) return 0.;
  else if (x <= -1.) return 1.;
  else               return (1. - x)/2.; // linear x \in [-1,1] -> p \in [1,0]
}


inline Real UniformRandomVariable::inverse_std_cdf(Real p_cdf)
{
  if      (p_cdf >= 1.) return  1.;
  else if (p_cdf <= 0.) return -1.;
  else return 2.*p_cdf - 1.; // linear p \in [0,1] -> x \in [-1,1]
}


inline Real UniformRandomVariable::inverse_std_ccdf(Real p_ccdf)
{
  if      (p_ccdf >= 1.) return -1.;
  else if (p_ccdf <= 0.) return  1.;
  else return 1. - 2.*p_ccdf; // linear p \in [1,0] -> x \in [-1,1]
}


inline Real UniformRandomVariable::pdf(Real lwr, Real upr)
{ return 1./(upr - lwr); } // equal probability on [lwr,upr]


inline Real UniformRandomVariable::cdf(Real x, Real lwr, Real upr)
{
  if      (x >= upr) return 1.;
  else if (x <= lwr) return 0.;
  else               return (x - lwr)/(upr - lwr); // linear [l,u] -> [0,1]
}


inline Real UniformRandomVariable::inverse_cdf(Real p_cdf, Real lwr, Real upr)
{
  if      (p_cdf >= 1.) return upr;
  else if (p_cdf <= 0.) return lwr;
  else                  return lwr + (upr - lwr) * p_cdf;
}


inline Real UniformRandomVariable::cdf(Real x) const
{ return cdf(x, lowerBnd, upperBnd); }


inline Real UniformRandomVariable::ccdf(Real x) const
{
  if      (x >= upperBnd) return 0.;
  else if (x <= lowerBnd) return 1.;
  else                    return (upperBnd - x)/(upperBnd - lowerBnd);
}


inline Real UniformRandomVariable::inverse_cdf(Real p_cdf) const
{ return inverse_cdf(p_cdf, lowerBnd, upperBnd); }


inline Real UniformRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  if      (p_ccdf >= 1.) return lowerBnd;
  else if (p_ccdf <= 0.) return upperBnd;
  else                   return upperBnd - (upperBnd - lowerBnd) * p_ccdf;
}


inline Real UniformRandomVariable::pdf(Real x) const
{ return pdf(lowerBnd, upperBnd); }


inline Real UniformRandomVariable::pdf_gradient(Real x) const
{ return 0.; }


inline Real UniformRandomVariable::pdf_hessian(Real x) const
{ return 0.; }


inline Real UniformRandomVariable::to_std(Real x) const
{
  // [L,U] -> [-1,1]
  if      (x >= upperBnd) return  1.;
  else if (x <= lowerBnd) return -1.;
  else                    return 2.*(x - lowerBnd)/(upperBnd - lowerBnd) - 1.;
}


inline Real UniformRandomVariable::from_std(Real z) const
{
  // [-1,1] -> [L,U]
  if      (z >=  1.) return upperBnd;
  else if (z <= -1.) return lowerBnd;
  else               return lowerBnd + (upperBnd - lowerBnd) * (z + 1.) / 2.;
}


inline Real UniformRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // correlation warping factor for transformations to STD_NORMAL space
  // Der Kiureghian and Liu, ASCE JEM 112:1, 1986
  Real COV;
  switch (rv.type()) { // x-space types mapped to STD_NORMAL u-space

  // Der Kiureghian & Liu: Table 4
  case UNIFORM:     return 1.047 - 0.047*corr*corr; break; // Max Error 0.0%
  case EXPONENTIAL: return 1.133 + 0.029*corr*corr; break; // Max Error 0.0%
  case GUMBEL:      return 1.055 + 0.015*corr*corr; break; // Max Error 0.0%

  // Der Kiureghian & Liu: Table 5 (quadratic approximations in COV,corr)
  case GAMMA:     // Max Error 0.1%
    COV = rv.coefficient_of_variation();
    return 1.023 + (-0.007 + 0.127*COV)*COV + 0.002*corr*corr; break;
  case FRECHET:   // Max Error 2.1%
    COV = rv.coefficient_of_variation();
    return 1.033 + ( 0.305 + 0.405*COV)*COV + 0.074*corr*corr; break;
  case WEIBULL:   // Max Error 0.5%
    COV = rv.coefficient_of_variation();
    return 1.061 + (-0.237 + 0.379*COV)*COV - 0.005*corr*corr; break;

  // warping factors are defined once for lower triangle based on uv order
  case NORMAL: case LOGNORMAL:
    return rv.correlation_warping_factor(*this, corr); break;

  default: // Unsupported warping (should be prevented upsteam)
    PCerr << "Error: unsupported correlation warping for UniformRV."<<std::endl;
    abort_handler(-1); return 1.; break;
  }
}


inline void UniformRandomVariable::update(Real lwr, Real upr)
{ lowerBnd = lwr; upperBnd = upr; }


inline void UniformRandomVariable::
moments_from_params(Real lwr, Real upr, Real& mean, Real& std_dev)
{ mean = (lwr + upr)/2.; std_dev = (upr - lwr)/std::sqrt(12.); }

} // namespace Pecos

#endif
