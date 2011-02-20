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
  /// constructor with mode argument
  PiecewiseInterpPolynomial(short interp_type);
  /// constructor with mode and set of points to interpolate
  PiecewiseInterpPolynomial(const RealArray& interp_pts, short interp_type);
  /// destructor
  ~PiecewiseInterpPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  const Real& get_type1_value(const Real& x, unsigned short i);
  const Real& get_type2_value(const Real& x, unsigned short i);

  const Real& get_type1_gradient(const Real& x, unsigned short i);
  const Real& get_type2_gradient(const Real& x, unsigned short i);

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


inline PiecewiseInterpPolynomial::PiecewiseInterpPolynomial(short interp_type):
  InterpolationPolynomial(), interpType(interp_type)
{ }


inline PiecewiseInterpPolynomial::
PiecewiseInterpPolynomial(const RealArray& interp_pts, short interp_type):
  InterpolationPolynomial(interp_pts), interpType(interp_type)
{ }


inline PiecewiseInterpPolynomial::~PiecewiseInterpPolynomial()
{ }

} // namespace Pecos

#endif
