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
  PiecewiseInterpPolynomial(short poly_type, short mode = NEWTON_COTES);
  /// constructor with mode and set of points to interpolate
  PiecewiseInterpPolynomial(const RealArray& interp_pts, short poly_type,
			    short mode = NEWTON_COTES);
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

  /// return the interpolation points corresponding to a point set of size order
  const RealArray& interpolation_points(unsigned short order);
  /// return the type 1 interpolation weights corresponding to a point
  /// set of size order
  const RealArray& type1_interpolation_weights(unsigned short order);
  /// return the type 2 interpolation weights corresponding to a point
  /// set of size order
  const RealArray& type2_interpolation_weights(unsigned short order);

private:

  //
  //- Heading: Data
  //

  /// type of polynomial interpolant: PIECEWISE_LINEAR_INTERP,
  /// PIECEWISE_QUADRATIC_INTERP, or PIECEWISE_CUBIC_INTERP
  short interpType;
  /// name of closed nested rule: NEWTON_COTES (equidistant) or
  /// CLENSHAW_CURTIS (non-equidistant)
  short interpMode;

  /// the constant interval between points for an equidistant interpMode
  Real interpInterval;

  /// set of 1-D weights for interpolation of values
  RealArray type1InterpWts;
  /// set of 1-D] weights for interpolation of gradients
  RealArray type2InterpWts;
};


inline PiecewiseInterpPolynomial::PiecewiseInterpPolynomial():
  InterpolationPolynomial()
{ }


inline PiecewiseInterpPolynomial::
PiecewiseInterpPolynomial(short poly_type, short mode):
  InterpolationPolynomial(), interpType(poly_type), interpMode(mode)
{ }


inline PiecewiseInterpPolynomial::
PiecewiseInterpPolynomial(const RealArray& interp_pts, short poly_type,
			  short mode):
  InterpolationPolynomial(interp_pts), interpType(poly_type), interpMode(mode)
{ }


inline PiecewiseInterpPolynomial::~PiecewiseInterpPolynomial()
{ }

} // namespace Pecos

#endif
