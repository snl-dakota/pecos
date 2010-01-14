/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HermiteOrthogPolynomial
//- Description:  Class for Hermite Orthogonal Polynomial
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef HERMITE_ORTHOG_POLYNOMIAL_HPP
#define HERMITE_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Hermite polynomials

/** The HermiteOrthogPolynomial class evaluates a univariate Hermite
    polynomial of a particular order.  It uses the "probabilist's"
    formulation for which the polynomials are orthogonal with respect
    to the weight function 1/std::sqrt(2*Pi) exp(-x^2/2) when integrated
    over the support range of [-infinity,+infinity].  It enables
    (mixed) multidimensional orthogonal polynomial basis functions
    within OrthogPolyApproximation. */

class HermiteOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  HermiteOrthogPolynomial();  ///< default constructor
  ~HermiteOrthogPolynomial(); ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// retrieve the Hermite polynomial value for a given parameter x 
  const Real& get_value(const Real& x, unsigned short order);
  /// retrieve the Hermite polynomial gradient for a given parameter x 
  const Real& get_gradient(const Real& x, unsigned short order);

  /// return the inner product <He_n,He_n> = ||He_n||^2
  const Real& norm_squared(unsigned short order);

  /// return the Gauss-Hermite quadrature points corresponding to
  /// polynomial order
  const RealArray& gauss_points(unsigned short order);
  /// return the Gauss-Hermite quadrature weights corresponding to
  /// polynomial order
  const RealArray& gauss_weights(unsigned short order);

private:

  //
  //- Heading: Data
  //

  /// numerical value of Pi
  static const Real Pi;
};


inline HermiteOrthogPolynomial::HermiteOrthogPolynomial()
{ ptFactor = std::sqrt(2.); wtFactor = 1./std::sqrt(Pi); }


inline HermiteOrthogPolynomial::~HermiteOrthogPolynomial()
{ }

} // namespace Pecos

#endif