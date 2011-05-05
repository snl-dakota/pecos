/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        OrthogonalPolynomial
//- Description:  Abstract base class for orthogonal polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef ORTHOGONAL_POLYNOMIAL_HPP
#define ORTHOGONAL_POLYNOMIAL_HPP

#include "BasisPolynomial.hpp"
#include "pecos_data_types.hpp"

namespace Pecos {


/// Base class for the orthogonal polynomial class hierarchy.

/** The OrthogonalPolynomial class is the base class for the
    univariate orthogonal polynomial class hierarchy in PECOS.  One
    instance of an OrthogonalPolynomial is created for each variable
    within a multidimensional orthogonal polynomial basis function (a
    vector of OrthogonalPolynomials is contained in
    OrthogPolyApproximation, which may be mixed and matched in, e.g.,
    the Wiener-Askey scheme for polynomial chaos). */

class OrthogonalPolynomial: public BasisPolynomial
{
public:

  //
  //- Heading: Constructors, destructor, assignment operator
  //

  OrthogonalPolynomial();  /// default constructor
  ~OrthogonalPolynomial(); /// destructor

  //
  //- Heading: Virtual function redefinitions
  //

  /// destroy history of Gauss pts/wts due to change in distribution parameters
  void reset_gauss();

  //
  //- Heading: Member functions
  //

  /// perform unit testing on Gauss points/weights
  void gauss_check(unsigned short order);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// set collocRule
  void collocation_rule(short rule);
  /// get collocRule
  short collocation_rule() const;

  //
  //- Heading: Data
  //

  /// norm-squared of the n_th order polynomial defined by the inner product
  /// <Poly_n, Poly_n> = ||Poly_n||^2 (returned by norm_squared())
  Real orthogPolyNormSq;

  /// collocation points for one-dimensional quadrature
  /// (x parameter values for which Poly_n(x) = 0)
  RealArray collocPoints;
  /// collocation weights for one-dimensional quadrature
  RealArray collocWeights;

  /// the type of integration rule associated with the orthogonal polynomial
  /** In most cases, this is just the corresponding Gauss quadrature
      rule.  However, for Legendre, collocRule manages the option of
      GAUSS_LEGENDRE or GAUSS_PATTERSON, for Chebyshev, it manages the
      option of CLENSHAW_CURTIS or FEJER2, and for Hermite, it manages
      the option of GAUSS_HERMITE or GENZ_KEISTER. */
  short collocRule;

private:

  //
  //- Heading: Data
  //
};


inline OrthogonalPolynomial::OrthogonalPolynomial():
  BasisPolynomial(BaseConstructor())
{ }


inline OrthogonalPolynomial::~OrthogonalPolynomial()
{ }


inline void OrthogonalPolynomial::reset_gauss()
{ collocPoints.clear(); collocWeights.clear(); }


inline void OrthogonalPolynomial::collocation_rule(short rule)
{ collocRule = rule; }


inline short OrthogonalPolynomial::collocation_rule() const
{ return collocRule; }

} // namespace Pecos

#endif
