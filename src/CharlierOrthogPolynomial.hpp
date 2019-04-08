/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        CharlierOrthogPolynomial
//- Description:  Class for Charlier Orthogonal Polynomial
//-               
//- Owner:        John Jakeman, Sandia National Laboratories

#ifndef CHARLIER_ORTHOG_POLYNOMIAL_HPP
#define CHARLIER_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"

namespace Pecos {

/**
 * \class CharlierOrthogPolynomial
 * \brief One-dimensional Charlier polynomial
 */
class CharlierOrthogPolynomial : public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //
  
  /// default constructor
  CharlierOrthogPolynomial();
  /// destructor
  ~CharlierOrthogPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

protected:

  Real type1_value( Real x, unsigned short order );
  Real type1_gradient( Real x, unsigned short order );
  Real type1_hessian( Real x, unsigned short order );
  Real norm_squared( unsigned short order );

  /// set betaPoly
  void parameter(short dist_param, Real param);
  /// get betaPoly
  Real parameter(short dist_param);

private: 
  
  /// Poisson distribution is the probability that a realization of a random
  /// variable X with mean lambda occurring k times in a fixed interval.

  /// expected value of the random variable X
  Real lambdaParam;
};


inline CharlierOrthogPolynomial::CharlierOrthogPolynomial():
  lambdaParam(-1.) // dummy value prior to update
{ }


inline CharlierOrthogPolynomial::~CharlierOrthogPolynomial()
{ }


inline Real CharlierOrthogPolynomial::parameter(short dist_param)
{
  switch (dist_param) {
  case P_LAMBDA: return lambdaParam; break;
  default:
    PCerr << "Error: unsupported distribution parameter in CharlierOrthog"
	  << "Polynomial::parameter()." << std::endl;
    abort_handler(-1);
    return 0.;
  }
}


inline void CharlierOrthogPolynomial::parameter(short dist_param, Real param)
{
  if (dist_param != P_LAMBDA) {
    PCerr << "Error: unsupported distribution parameter in CharlierOrthog"
	  << "Polynomial::parameter()." << std::endl;
    abort_handler(-1);
  }

  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Therefore, set parametricUpdate to false unless an actual parameter change.
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPoints.empty() || collocWeights.empty()) { // first pass
    parametricUpdate = true; // prevent false if default value assigned
    lambdaParam = param;
  }
  else {
    parametricUpdate = false;
    if (!real_compare(lambdaParam, param))
      { lambdaParam = param; parametricUpdate = true; reset_gauss(); }
  }
}

} // namespace Pecos

#endif // CHARLIER_ORTHOG_POLYNOMIAL_HPP
