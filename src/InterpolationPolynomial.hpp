/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        InterpolationPolynomial
//- Description:  Class for 1-D Interpolation Polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef INTERPOLATION_POLYNOMIAL_HPP
#define INTERPOLATION_POLYNOMIAL_HPP

#include "BasisPolynomial.hpp"
#include "pecos_data_types.hpp"

namespace Pecos {


/// Derived basis polynomial class for 1-D Lagrange interpolation polynomials

/** The InterpolationPolynomial class evaluates a univariate
    interpolation polynomial.  It enables multidimensional
    interpolants within InterpPolyApproximation. */

class InterpolationPolynomial: public BasisPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  InterpolationPolynomial();
  /// standard constructor
  InterpolationPolynomial(const RealArray& interp_pts);
  /// destructor
  ~InterpolationPolynomial();

  //
  //- Heading: Virtual functions
  //

  /// precompute data that is reused repeatedly within interpolation
  virtual void precompute_data() = 0;

  //
  //- Heading: Set/get functions
  //

  /// set interpPts
  void interpolation_points(const RealArray& interp_pts);

protected:

  //
  //- Heading: Data
  //

  /// set of 1-D interpolation points: the i_th interpolation polynomial
  /// evaluated at the j_th interpolation point produces Kronecker delta_ij
  RealArray interpPts;
  /// number of 1-D interpolation points
  size_t numInterpPts;

private:

  //
  //- Heading: Data
  //

};


inline InterpolationPolynomial::InterpolationPolynomial():
  BasisPolynomial(BaseConstructor())
{ }


inline InterpolationPolynomial::
InterpolationPolynomial(const RealArray& interp_pts):
  BasisPolynomial(BaseConstructor()), interpPts(interp_pts),
  numInterpPts(interpPts.size())
{ precompute_data(); }


inline InterpolationPolynomial::~InterpolationPolynomial()
{ }


inline void InterpolationPolynomial::
interpolation_points(const RealArray& interp_pts)
{
  interpPts    = interp_pts;
  numInterpPts = interpPts.size();
  precompute_data();
}

} // namespace Pecos

#endif
