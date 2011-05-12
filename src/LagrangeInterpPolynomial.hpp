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

#include "InterpolationPolynomial.hpp"
#include "pecos_data_types.hpp"


namespace Pecos {

/// Derived basis polynomial class for 1-D Lagrange interpolation polynomials

/** The LagrangeInterpPolynomial class evaluates a univariate Lagrange
    interpolation polynomial.  The order of the polynomial is dictated
    by the number of interpolation points (order = N_p - 1).  It enables
    multidimensional interpolants within InterpPolyApproximation. */

class LagrangeInterpPolynomial: public InterpolationPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  LagrangeInterpPolynomial();
  /// standard constructor
  LagrangeInterpPolynomial(const RealArray& interp_pts);
  /// destructor
  ~LagrangeInterpPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  /// retrieve the value of the i_th Lagrange polynomial for a given
  /// parameter x
  const Real& type1_value(const Real& x, unsigned short i);
  /// retrieve the gradient of the i_th Lagrange polynomial for a given
  /// parameter x
  const Real& type1_gradient(const Real& x, unsigned short i);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void precompute_data();

private:

  //
  //- Heading: Data
  //

  /// set of denominator products calculated from interpPts in
  /// precompute_data()
  RealVector lagDenominators;
};


inline LagrangeInterpPolynomial::LagrangeInterpPolynomial():
  InterpolationPolynomial()
{ }


inline LagrangeInterpPolynomial::
LagrangeInterpPolynomial(const RealArray& interp_pts):
  InterpolationPolynomial(interp_pts)
{ }


inline LagrangeInterpPolynomial::~LagrangeInterpPolynomial()
{ }

} // namespace Pecos

#endif
