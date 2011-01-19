/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        PiecewiseInterpPolynomial
//- Description:  Class for 1-D piecewise Interpolation Polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef PIECEWISE_INTERP_POLYNOMIAL_HPP
#define PIECEWISE_INTERP_POLYNOMIAL_HPP

#include "InterpolationPolynomial.hpp"
#include "pecos_data_types.hpp"

namespace Pecos {


/// Derived basis polynomial class for 1-D piecewise interpolation polynomials

/** The PiecewiseInterpPolynomial class evaluates a univariate
    interpolation polynomial with local support.  The order of the
    polynomial may be linear, based only on interpolated values, or
    cubic, based on interpolated values and gradients.  It enables
    multidimensional interpolants within InterpPolyApproximation. */

class PiecewiseInterpPolynomial: public InterpolationPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  PiecewiseInterpPolynomial();
  /// standard constructor
  PiecewiseInterpPolynomial(const RealArray& interp_pts);
  /// destructor
  ~PiecewiseInterpPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  /// retrieve the value of the i_th piecewise polynomial for a given
  /// parameter x
  const Real& get_value(const Real& x, unsigned short i);
  /// retrieve the gradient of the i_th piecewise polynomial for a
  /// given parameter x
  const Real& get_gradient(const Real& x, unsigned short i);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void precompute_data();

private:

  //
  //- Heading: Data
  //

  /// type of polynomial interpolant: LINEAR, QUADRATIC, CUBIC,
  /// {LINEAR,QUADRATIC,CUBIC}_EQUIDISTANT
  short interpType;

  /// the constant interval between points for *_EQUIDISTANT interpTypes
  Real interpInterval;
};


inline PiecewiseInterpPolynomial::PiecewiseInterpPolynomial():
  InterpolationPolynomial()
{ }


inline PiecewiseInterpPolynomial::
PiecewiseInterpPolynomial(const RealArray& interp_pts):
  InterpolationPolynomial(interp_pts)
{ }


inline PiecewiseInterpPolynomial::~PiecewiseInterpPolynomial()
{ }

} // namespace Pecos

#endif
