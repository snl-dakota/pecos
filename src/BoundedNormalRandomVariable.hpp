/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 BoundedNormalRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef BOUNDED_NORMAL_RANDOM_VARIABLE_HPP
#define BOUNDED_NORMAL_RANDOM_VARIABLE_HPP

#include "NormalRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for bounded normal random variables.

/** Manages lower and upper bounds. */

class BoundedNormalRandomVariable: public NormalRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  BoundedNormalRandomVariable();
  /// alternate constructor
  BoundedNormalRandomVariable(Real mean, Real stdev, Real lwr, Real upr);
  /// destructor
  ~BoundedNormalRandomVariable();

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

  //Real log_pdf(Real x) const;

  Real parameter(short dist_param) const;
  void parameter(short dist_param, Real val);

  RealRealPair bounds() const;

  Real dx_ds(short dist_param, short u_type, Real x, Real z) const;
  Real dz_ds_factor(short u_type, Real x, Real z) const;

  void update(Real mean, Real stdev, Real lwr, Real upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real mean, Real std_dev, Real lwr, Real upr);
  static Real cdf(Real x, Real mean, Real std_dev, Real lwr, Real upr);

protected:

  //
  //- Heading: Data
  //

  /// lower bound of bounded_normal random variable
  Real lowerBnd;
  /// upper bound of bounded_normal random variable
  Real upperBnd;
};


inline BoundedNormalRandomVariable::BoundedNormalRandomVariable():
  NormalRandomVariable(), lowerBnd(-std::numeric_limits<Real>::infinity()),
  upperBnd(std::numeric_limits<Real>::infinity())
{ }


inline BoundedNormalRandomVariable::
BoundedNormalRandomVariable(Real mean, Real stdev, Real lwr, Real upr):
  NormalRandomVariable(mean, stdev), lowerBnd(lwr), upperBnd(upr)
{ }


inline BoundedNormalRandomVariable::~BoundedNormalRandomVariable()
{ }


inline Real BoundedNormalRandomVariable::
pdf(Real x, Real mean, Real std_dev, Real lwr, Real upr)
{
  if (x < lwr || x > upr) return 0.;
  else {
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    Real Phi_lms = (lwr > -dbl_inf) ? std_cdf((lwr-mean)/std_dev) : 0.;
    Real Phi_ums = (upr <  dbl_inf) ? std_cdf((upr-mean)/std_dev) : 1.;
    return std_pdf((x-mean)/std_dev)/(Phi_ums - Phi_lms)/std_dev;
  }
}


inline Real BoundedNormalRandomVariable::
cdf(Real x, Real mean, Real std_dev, Real lwr, Real upr)
{
  if      (x < lwr) return 0.;
  else if (x > upr) return 1.;
  else {
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    Real Phi_lms = (lwr > -dbl_inf) ? std_cdf((lwr-mean)/std_dev) : 0.;
    Real Phi_ums = (upr <  dbl_inf) ? std_cdf((upr-mean)/std_dev) : 1.;
    return (std_cdf((x-mean)/std_dev) - Phi_lms) / (Phi_ums - Phi_lms);
  }
}


inline Real BoundedNormalRandomVariable::cdf(Real x) const
{ return cdf(x, gaussMean, gaussStdDev, lowerBnd, upperBnd); }


inline Real BoundedNormalRandomVariable::ccdf(Real x) const
{
  if      (x < lowerBnd) return 1.;
  else if (x > upperBnd) return 0.;
  else {
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    Real Phi_lms = (lowerBnd > -dbl_inf) ?
      std_cdf((lowerBnd-gaussMean)/gaussStdDev) : 0.;
    Real Phi_ums = (upperBnd <  dbl_inf) ?
      std_cdf((upperBnd-gaussMean)/gaussStdDev) : 1.;
    return (Phi_ums - std_cdf((x-gaussMean)/gaussStdDev)) / (Phi_ums - Phi_lms);
  }
}


inline Real BoundedNormalRandomVariable::inverse_cdf(Real p_cdf) const
{
  if      (p_cdf <= 0.) return lowerBnd;
  else if (p_cdf >= 1.) return upperBnd;
  else {
    // p = (Phi((x-mean)/std_dev) - Phi_lms)/(Phi_ums - Phi_lms)
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    Real Phi_lms = (lowerBnd > -dbl_inf) ?
      std_cdf((lowerBnd-gaussMean)/gaussStdDev) : 0.;
    Real Phi_ums = (upperBnd <  dbl_inf) ?
      std_cdf((upperBnd-gaussMean)/gaussStdDev) : 1.;
    return gaussMean + gaussStdDev *
      inverse_std_cdf(p_cdf * (Phi_ums - Phi_lms) + Phi_lms);
  }
}


inline Real BoundedNormalRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  if      (p_ccdf >= 1.) return lowerBnd;
  else if (p_ccdf <= 0.) return upperBnd;
  else {
    // p = (Phi_ums - Phi((x-mean)/std_dev))/(Phi_ums - Phi_lms)
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    Real Phi_lms = (lowerBnd > -dbl_inf) ?
      std_cdf((lowerBnd-gaussMean)/gaussStdDev) : 0.;
    Real Phi_ums = (upperBnd <  dbl_inf) ?
      std_cdf((upperBnd-gaussMean)/gaussStdDev) : 1.;
    return gaussMean + gaussStdDev *
      inverse_std_cdf(Phi_ums - p_ccdf * (Phi_ums - Phi_lms));
  }
}


inline Real BoundedNormalRandomVariable::pdf(Real x) const
{ return pdf(x, gaussMean, gaussStdDev, lowerBnd, upperBnd); }


inline Real BoundedNormalRandomVariable::pdf_gradient(Real x) const
{ return pdf(x) * (gaussMean - x) / (gaussStdDev * gaussStdDev); }


//inline Real BoundedNormalRandomVariable::pdf_hessian(Real x) const
//{
//  return pdf(x) * ...; // TO DO
//}


inline Real BoundedNormalRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case N_MEAN:    case N_LOCATION: return gaussMean;   break;
  case N_STD_DEV: case N_SCALE:    return gaussStdDev; break;
  case N_LWR_BND:                  return lowerBnd; break;
  case N_UPR_BND:                  return upperBnd; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BoundedNormalRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void BoundedNormalRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case N_MEAN:    gaussMean   = val; break;
  case N_STD_DEV: gaussStdDev = val; break;
  // Bounded normal case must translate/scale bounds for N_LOCATION,N_SCALE
  // (see NestedModel::real_variable_mapping())
  //case N_LOCATION: TO DO; break;
  //case N_SCALE:    TO DO; break;
  case N_LWR_BND: lowerBnd = val;    break;
  case N_UPR_BND: upperBnd = val;    break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BoundedNormalRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
}


inline RealRealPair BoundedNormalRandomVariable::bounds() const
{ return RealRealPair(lowerBnd, upperBnd); }


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  dz/ds is zero if uncorrelated, 
    while dz_ds_factor() manages contributions in the correlated case. */
inline Real BoundedNormalRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  bool u_type_err = false, dist_err = false;
  switch (u_type) {
  case STD_NORMAL: {
    Real num, lms, ums, xms = (x-gaussMean)/gaussStdDev, phi_x = std_pdf(xms),
      dbl_inf = std::numeric_limits<Real>::infinity();
    switch (dist_param) {
    case N_MEAN:
      num = 0.;
      if (lowerBnd > -dbl_inf) {
	lms = (lowerBnd-gaussMean)/gaussStdDev;
	num += std_ccdf(z) * std_pdf(lms);
      }
      if (upperBnd <  dbl_inf) {
	ums = (upperBnd-gaussMean)/gaussStdDev;
	num += std_cdf(z)  * std_pdf(ums);
      }
      return 1. - num / phi_x;                                           break;
    case N_STD_DEV:
      num = 0.;
      if (lowerBnd > -dbl_inf) {
        lms = (lowerBnd-gaussMean)/gaussStdDev;
	num += std_ccdf(z) * std_pdf(lms) * lms;
      }
      if (upperBnd <  dbl_inf) {
        ums = (upperBnd-gaussMean)/gaussStdDev;
	num += std_cdf(z)  * std_pdf(ums) * ums;
      }
      return xms - num / phi_x;                                          break;
    case N_LWR_BND:
      lms = (lowerBnd-gaussMean)/gaussStdDev;
      return std_ccdf(z) * std_pdf(lms) / phi_x;                         break;
    case N_UPR_BND:
      ums = (upperBnd-gaussMean)/gaussStdDev;
      return std_cdf(z)  * std_pdf(ums) / phi_x;                         break;
    //case N_LOCATION:
    //case N_SCALE:
    default: dist_err = true;                                            break;
    }
    break;
  }
  //case BOUNDED_NORMAL: TO DO;                                        break;
  default: u_type_err = true;                                          break;
  }

  if (u_type_err)
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in BoundedNormalRandomVariable::dx_ds()." << std::endl;
  if (dist_err)
    PCerr << "Error: mapping failure for distribution parameter " << dist_param
	  << " in BoundedNormalRandomVariable::dx_ds()." << std::endl;
  if (u_type_err || dist_err)
    abort_handler(-1);
  return 0.;
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  For the uncorrelated case,
    u and z are constants.  For the correlated case, u is a constant, but 
    z(s) = L(s) u due to Nataf dependence on s and dz/ds = dL/ds u. */
inline Real BoundedNormalRandomVariable::
dz_ds_factor(short u_type, Real x, Real z) const
{
  Real num = 0., lms, ums, xms = (x-gaussMean)/gaussStdDev,
    dbl_inf = std::numeric_limits<Real>::infinity();
  switch (u_type) {
  case STD_NORMAL:
    if (upperBnd <  dbl_inf)
      { ums = (upperBnd-gaussMean)/gaussStdDev; num += std_cdf(ums); }
    else num += 1.;
    if (lowerBnd > -dbl_inf)
      { lms = (lowerBnd-gaussMean)/gaussStdDev; num -= std_cdf(lms); }
    return gaussStdDev * std_pdf(z) * num / std_pdf(xms);
    break;
  //case BOUNDED_NORMAL: TO DO; break;
  default:
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in BoundedNormalRandomVariable::dz_ds_factor()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void BoundedNormalRandomVariable::
update(Real mean, Real stdev, Real lwr, Real upr)
{ gaussMean = mean; gaussStdDev = stdev; lowerBnd = lwr; upperBnd = upr; }

} // namespace Pecos

#endif
