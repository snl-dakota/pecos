/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef BASIS_APPROXIMATION_HPP
#define BASIS_APPROXIMATION_HPP

#include "pecos_data_types.hpp"


namespace Pecos {

/// Base class for multivariate basis approximations used for
/// projection of random variables through time or space

/** The base class for basis approximations defined from Fourier
    functions, eigenfunctions, or polynomial functions. */

class BasisApproximation
{
public:

  //
  //- Heading: Constructors, destructor, operator=
  //

  /// default constructor
  BasisApproximation();
  /// standard constructor for envelope
  BasisApproximation(const String& basis_approx_type, size_t num_vars);
  /// copy constructor
  BasisApproximation(const BasisApproximation& basis_approx);

  /// destructor
  virtual ~BasisApproximation();

  /// assignment operator
  BasisApproximation operator=(const BasisApproximation& basis_approx);

  //
  //- Heading: Virtual functions
  //

  /// retrieve the approximate function value for a given parameter vector
  virtual const Real& get_value(const RealVector& x);
  /// retrieve the approximate function gradient for a given parameter vector
  virtual const RealVector& get_gradient(const RealVector& x);
  /// retrieve the approximate function Hessian for a given parameter vector
  virtual const RealSymMatrix& get_hessian(const RealVector& x);

  /// return the minimum number of samples (unknowns) required to
  /// build the derived class approximation type in numVars dimensions
  virtual int min_coefficients() const;
  /// calculate the data fit coefficients using currentPoints and anchorPoint
  virtual void find_coefficients();
  /// print the coefficient array computed in find_coefficients()
  virtual void print_coefficients(std::ostream& s) const;

  /// return the coefficient array computed by find_coefficients()
  virtual const RealVector& approximation_coefficients() const;
  /// set the coefficient array from external sources, rather than
  /// computing with find_coefficients()
  virtual void approximation_coefficients(const RealVector& approx_coeffs);

  //
  //- Heading: Member functions
  //

  /// returns approxRep for access to derived class member functions
  /// that are not mapped to the top Approximation level
  BasisApproximation* approx_rep() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  BasisApproximation(BaseConstructor, size_t num_vars);

  //
  //- Heading: Data members
  //

  /// number of variables used in the approximation
  size_t numVars;

  /// value of the approximation returned by get_value()
  Real approxValue;
  /// gradient of the approximation returned by get_gradient()
  RealVector approxGradient;

private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// basisApproxRep to the appropriate derived type.
  BasisApproximation* get_basis_approx(const String& basis_approx_type,
				       size_t num_vars);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  BasisApproximation* basisApproxRep;
  /// number of objects sharing basisApproxRep
  int referenceCount;
};


inline BasisApproximation* BasisApproximation::approx_rep() const
{ return basisApproxRep; }

} // namespace Pecos

#endif