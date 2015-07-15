/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 PoissonRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef POISSON_RANDOM_VARIABLE_HPP
#define POISSON_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for poisson random variables.

/** Manages alpha and beta parameters. */

class PoissonRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  PoissonRandomVariable();
  /// alternate constructor
  PoissonRandomVariable(Real lambda);
  /// destructor
  ~PoissonRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;

  Real parameter(short dist_param) const;
  void parameter(short dist_param, Real val);

  RealRealPair moments() const;
  RealRealPair bounds() const;

  //
  //- Heading: Member functions
  //

  void update(Real lambda);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real lambda);
  static Real cdf(Real x, Real lambda);

  static void moments_from_params(Real lambda, Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// lambda parameter of poisson random variable
  Real poissonLambda;

  /// pointer to the Boost poisson_distribution instance
  poisson_dist* poissonDist;
};


inline PoissonRandomVariable::PoissonRandomVariable():
  RandomVariable(BaseConstructor()), poissonLambda(1.), poissonDist(NULL)
{ }


inline PoissonRandomVariable::PoissonRandomVariable(Real lambda):
  RandomVariable(BaseConstructor()), poissonLambda(lambda),
  poissonDist(new poisson_dist(lambda))
{ }


inline PoissonRandomVariable::~PoissonRandomVariable()
{ if (poissonDist) delete poissonDist; }


inline Real PoissonRandomVariable::cdf(Real x) const
{ return bmth::cdf(*poissonDist, x); }


inline Real PoissonRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*poissonDist, x)); }


inline Real PoissonRandomVariable::inverse_cdf(Real p_cdf) const
{ return bmth::quantile(*poissonDist, p_cdf); }


inline Real PoissonRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*poissonDist, p_ccdf)); }


inline Real PoissonRandomVariable::pdf(Real x) const
{ return bmth::pdf(*poissonDist, x); }


inline Real PoissonRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case P_LAMBDA: return poissonLambda; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in PoissonRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void PoissonRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case P_LAMBDA: poissonLambda = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in PoissonRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
}


inline RealRealPair PoissonRandomVariable::moments() const
{
  Real mean, std_dev;
  moments_from_params(poissonLambda, mean, std_dev);
  return RealRealPair(mean, std_dev);
}


inline RealRealPair PoissonRandomVariable::bounds() const
{ return RealRealPair(0., std::numeric_limits<Real>::infinity()); }


inline void PoissonRandomVariable::update(Real lambda)
{
  if (!poissonDist || poissonLambda != lambda) {
    poissonLambda = lambda;
    if (poissonDist) delete poissonDist;
    poissonDist = new poisson_dist(poissonLambda);
  }
}

// Static member functions:

inline Real PoissonRandomVariable::pdf(Real x, Real lambda)
{ poisson_dist poisson1(lambda); return bmth::pdf(poisson1, x); }


inline Real PoissonRandomVariable::cdf(Real x, Real lambda)
{ poisson_dist poisson1(lambda); return bmth::cdf(poisson1, x); }


inline void PoissonRandomVariable::
moments_from_params(Real lambda, Real& mean, Real& std_dev)
{ mean = lambda; std_dev = std::sqrt(lambda); }

} // namespace Pecos

#endif
