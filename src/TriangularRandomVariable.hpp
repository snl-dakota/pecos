/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 TriangularRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef TRIANGULAR_RANDOM_VARIABLE_HPP
#define TRIANGULAR_RANDOM_VARIABLE_HPP

#include "UniformRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for triangular random variables.

/** Manages mode and inherits bounds.  See Haldar and Mahadevan, p. 99. */

class TriangularRandomVariable: public UniformRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  TriangularRandomVariable();
  /// alternate constructor
  TriangularRandomVariable(Real lwr, Real mode, Real upr);
  /// destructor
  ~TriangularRandomVariable();

  //
  //- Heading: Member functions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;

  void update(Real lwr, Real mode, Real upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real lwr, Real mode, Real upr);
  static Real cdf(Real x, Real lwr, Real mode, Real upr);

  static void moments_from_params(Real lwr, Real upr, Real mode, Real& mean,
				  Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// mode of triangular random variable
  Real triangularMode;

  /// pointer to the Boost gamma_distribution instance
  triangular_dist* triangDist;
};


inline TriangularRandomVariable::TriangularRandomVariable():
  UniformRandomVariable(), triangularMode(0), triangDist(NULL)
{ }


inline TriangularRandomVariable::
TriangularRandomVariable(Real lwr, Real mode, Real upr):
  UniformRandomVariable(lwr, upr), triangularMode(mode),
  triangDist(new triangular_dist(lwr, mode, upr))
{ }


inline TriangularRandomVariable::~TriangularRandomVariable()
{ if (triangDist) delete triangDist; }


inline Real TriangularRandomVariable::pdf(Real x, Real lwr, Real mode, Real upr)
{
  triangular_dist tri1(lwr, mode, upr);
  return bmth::pdf(tri1, x);

  //return (x < mode) ? 2.*(x-lwr)/(upr-lwr)/(mode-lwr) :
  //                    2.*(upr-x)/(upr-lwr)/(upr-mode);
}


inline Real TriangularRandomVariable::cdf(Real x, Real lwr, Real mode, Real upr)
{
  triangular_dist tri1(lwr, mode, upr);
  return bmth::cdf(tri1, x);

  //return (x < mode) ? std::pow(x-lwr,2.)/(upr-lwr)/(mode-lwr) :
  //  ((mode-lwr) - (x+mode-2*upr)*(x-mode)/(upr-mode))/(upr-lwr);
}


inline Real TriangularRandomVariable::cdf(Real x) const
{ return bmth::cdf(*triangDist, x); }


inline Real TriangularRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*triangDist, x)); }


inline Real TriangularRandomVariable::inverse_cdf(Real p_cdf) const
{
  return bmth::quantile(*triangDist, p_cdf);

  /*
  // assume x < mode and then check
  Real range = upperBnd - lowerBnd,
       x = lowerBnd + std::sqrt(p*range*(triangularMode-lowerBnd));
  Real x_pdf = 2.*(x-lowerBnd)/range/(triangularMode-lowerBnd),
       m_pdf = 2./range;
  // check pdf value to ensure that proper equation used
  if ( x_pdf > m_pdf )
    x = upperBnd - std::sqrt((1.-p)*range*(upperBnd-triangularMode));
  return x;
  */
}


inline Real TriangularRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*triangDist, p_ccdf)); }


inline Real TriangularRandomVariable::pdf(Real x) const
{ return bmth::pdf(*triangDist, x); }


inline Real TriangularRandomVariable::pdf_gradient(Real x) const
{
  Real range = upperBnd - lowerBnd;
  if (x < triangularMode)
    return  2. / ( range * (triangularMode - lowerBnd) );
  else if (x > triangularMode)
    return -2. / ( range * (upperBnd - triangularMode) );
  else // x == triangularMode
    return  0.; // f'(x) is undefined: use 0.
}


inline Real TriangularRandomVariable::pdf_hessian(Real x) const
{ return 0.; }


inline void TriangularRandomVariable::update(Real lwr, Real mode, Real upr)
{
  if (!triangDist ||
      lowerBnd != lwr || triangularMode != mode || upperBnd != upr) {
    lowerBnd = lwr; triangularMode = mode; upperBnd = upr;
    if (triangDist) delete triangDist;
    triangDist = new triangular_dist(lwr, mode, upr);
  }
}


inline void TriangularRandomVariable::
moments_from_params(Real lwr, Real mode, Real upr, Real& mean, Real& std_dev)
{
  mean = (lwr + mode + upr)/3.;
  std_dev
    = std::sqrt((lwr*(lwr - mode) + mode*(mode - upr) + upr*(upr - lwr))/18.);
}

} // namespace Pecos

#endif
