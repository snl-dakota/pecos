/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LagrangeInterpPolynomial
//- Description:  Class for 1-D Lagrange Interpolation Polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef LAGRANGE_INTERP_POLYNOMIAL_HPP
#define LAGRANGE_INTERP_POLYNOMIAL_HPP

#include "BasisPolynomial.hpp"
#include "pecos_data_types.hpp"


namespace Pecos {

/// Derived basis polynomial class for 1-D Lagrange interpolation polynomials

/** The LagrangeInterpPolynomial class evaluates a univariate Lagrange
    interpolation polynomial.  The order of the polynomial is dictated
    by the number of interpolation points (order = N_p - 1).  It enables
    multidimensional interpolants within InterpPolyApproximation. */

class LagrangeInterpPolynomial: public BasisPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  LagrangeInterpPolynomial();
  /// standard constructor
  LagrangeInterpPolynomial(const RealArray& interpolation_pts);
  /// destructor
  ~LagrangeInterpPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  /// retrieve the value of the i_th Lagrange polynomial for a given
  /// parameter x
  const Real& get_value(const Real& x, unsigned short i);
  /// retrieve the gradient of the i_th Lagrange polynomial for a
  /// given parameter x
  const Real& get_gradient(const Real& x, unsigned short i);

  //
  //- Heading: Set/get functions
  //

  /// set interpolationPts
  void interpolation_points(const RealArray& interpolation_pts);

private:

  //
  //- Heading: Convenience functions
  //

  /// precompute data that is reused repeatedly within Lagrange interpolation
  void precompute_data();

  //
  //- Heading: Data
  //

  /// set of 1-D interpolation points: the i_th interpolation polynomial
  /// evaluated at the j_th interpolation point produces Kronecker delta_ij
  RealArray interpolationPts;
  /// number of 1-D interpolation points
  size_t numInterpPts;

  /// set of denominator products calculated from interpolationPts in
  /// precompute_data()
  RealVector lagDenominators;
};


inline LagrangeInterpPolynomial::LagrangeInterpPolynomial():
  BasisPolynomial(BaseConstructor())
{ }


inline LagrangeInterpPolynomial::
LagrangeInterpPolynomial(const RealArray& interpolation_pts):
  BasisPolynomial(BaseConstructor()), interpolationPts(interpolation_pts),
  numInterpPts(interpolationPts.size())
{ precompute_data(); }


inline LagrangeInterpPolynomial::~LagrangeInterpPolynomial()
{ }


inline void LagrangeInterpPolynomial::
interpolation_points(const RealArray& interpolation_pts)
{
  interpolationPts = interpolation_pts;
  numInterpPts = interpolationPts.size();
  precompute_data();
}

} // namespace Pecos

#endif
