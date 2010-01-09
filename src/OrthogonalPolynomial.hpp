/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
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

  /// destroy history of Gauss pts/wts due to change in alpha/beta stats
  void reset_gauss();

  //
  //- Heading: Member functions
  //

  /// perform unit testing on the Gauss points/weights
  void gauss_check(unsigned short order);

protected:

  //
  //- Heading: Data
  //

  /// norm-squared of the n_th order polynomial defined by the inner product
  /// <Poly_n, Poly_n> = ||Poly_n||^2 (returned by norm_squared())
  Real orthogPolyNormSq;

  /// Gauss points for one-dimensional Gaussian quadrature
  /// (x parameter values for which Poly_n(x) = 0)
  RealArray gaussPoints;
  /// Gauss weights for one-dimensional Gaussian quadrature
  RealArray gaussWeights;

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
{ gaussPoints.clear(); gaussWeights.clear(); }

} // namespace Pecos

#endif
