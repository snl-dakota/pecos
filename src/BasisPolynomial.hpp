/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        BasisPolynomial
//- Description:  Abstract base class for basis polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef BASIS_POLYNOMIAL_HPP
#define BASIS_POLYNOMIAL_HPP

#include "pecos_data_types.hpp"

namespace Pecos {


/// Base class for the basis polynomial class hierarchy.

/** The BasisPolynomial class is the base class for the univariate
    basis polynomial class hierarchy in PECOS.  One instance of an
    BasisPolynomial is created for each variable within a
    multidimensional polynomial basis function (a vector of
    BasisPolynomials is contained in BasisPolyApproximation, which may
    be mixed and matched in, e.g., the Wiener-Askey scheme for
    polynomial chaos).  For memory efficiency and enhanced
    polymorphism, the basis polynomial hierarchy employs the
    "letter/envelope idiom" (see Coplien "Advanced C++", p. 133), for
    which the base class (BasisPolynomial) serves as the envelope and
    one of the derived classes (selected in
    BasisPolynomial::get_polynomial()) serves as the letter. */

class BasisPolynomial
{
public:

  //
  //- Heading: Constructors, destructor, assignment operator
  //

  /// default constructor
  BasisPolynomial();
   /// alternate constructor
  BasisPolynomial(short poly_type, short rule = 0);
  /// copy constructor
  BasisPolynomial(const BasisPolynomial& polynomial);

  /// destructor
  virtual ~BasisPolynomial();

  /// assignment operator
  BasisPolynomial operator=(const BasisPolynomial& polynomial);

  //
  //- Heading: Virtual functions
  //

  /// retrieve the value of the n_th type 1 polynomial for a given parameter x
  /** For orthogonal polynomials, n specifies the order of the polynomial,
      whereas for interpolation polynomials, it identifies the interpolant
      for the n-th point. */
  virtual Real type1_value(const Real& x, unsigned short n);
  /// retrieve the value of the n_th type 2 polynomial for a given parameter x
  /** For orthogonal polynomials, n specifies the order of the polynomial,
      whereas for interpolation polynomials, it identifies the interpolant
      for the n-th point. */
  virtual Real type2_value(const Real& x, unsigned short n);

  /// retrieve the gradient of the n_th type 1 polynomial for a given
  /// parameter x
  /** For orthogonal polynomials, n specifies the order of the polynomial,
      whereas for interpolation polynomials, it identifies the interpolant
      for the n-th point. */
  virtual Real type1_gradient(const Real& x, unsigned short n);
  /// retrieve the gradient of the n_th type 2 polynomial for a given
  /// parameter x
  /** For orthogonal polynomials, n specifies the order of the polynomial,
      whereas for interpolation polynomials, it identifies the interpolant
      for the n-th point. */
  virtual Real type2_gradient(const Real& x, unsigned short n);

  /// returns the norm-squared of the n_th order polynomial defined by the
  /// inner product <Poly_n, Poly_n> = ||Poly_n||^2
  /** This is defined only for orthogonal polynomials. */
  virtual Real norm_squared(unsigned short n);

  /// return collocation points corresponding to orthogonal polynomial order n
  /** This is defined for orthogonal and piecewise interpolation polynomials. */
  virtual const RealArray& collocation_points(unsigned short n);
  /// return the type 1 collocation weights corresponding to a point
  /// set of size order
  /** This is defined for orthogonal and piecewise interpolation polynomials. */
  virtual const RealArray& type1_collocation_weights(unsigned short order);
  /// return the type 2 collocation weights corresponding to a point
  /// set of size order
  /** This is defined for piecewise interpolation polynomials. */
  virtual const RealArray& type2_collocation_weights(unsigned short order);

  /// destroy history of Gauss pts/wts (due to distribution parameter changes)
  /** This is defined only for orthogonal polynomials. */
  virtual void reset_gauss();

  /// (calculate and) return ptFactor
  virtual const Real& point_factor();
  /// (calculate and) return wtFactor
  virtual const Real& weight_factor();

  /// return {Jacobi,GenLaguerre}OrthogPolynomial::alphaPoly
  /** This is defined only for parameterized orthogonal polynomials. */
  virtual const Real& alpha_polynomial() const;
  /// return JacobiOrthogPolynomial::betaPoly
  /** This is defined only for parameterized orthogonal polynomials. */
  virtual const Real& beta_polynomial() const;
  /// set JacobiOrthogPolynomial::betaPoly or
  /// GenLaguerreOrthogPolynomial::alphaPoly from statistical defn of alpha
  /** This is defined only for parameterized orthogonal polynomials. */
  virtual void alpha_stat(const Real& alpha);
  /// set JacobiOrthogPolynomial::alphaPoly from statistical defn of beta
  /** This is defined only for parameterized orthogonal polynomials. */
  virtual void beta_stat(const Real& beta);

  /// set OrthogonalPolynomial::collocRule
  virtual void collocation_rule(short rule);
  /// get OrthogonalPolynomial::collocRule
  virtual short collocation_rule() const;

  /// set InterpolationPolynomial::interpolationPts
  /** This is defined only for interpolation polynomials. */
  virtual void interpolation_points(const RealArray& interpolation_pts);
  /// set InterpolationPolynomial::interpolationPts
  /** This is defined only for interpolation polynomials. */
  virtual const RealArray& interpolation_points() const;

  //
  //- Heading: Member functions
  //

  /// compute n!
  static Real factorial(unsigned short n);
  /// compute num!/den!
  static Real factorial_ratio(unsigned short num, unsigned short den);
  /// compute n!/(k!(n-k)!)
  /** Note: n+p choose n or p = (n+p)!/n!/p! */
  static Real n_choose_k(unsigned short n, unsigned short k);
  /// compute the Pochhammer symbol (m)_n = m*(m+1)...*(m+n-1)
  static Real pochhammer(const Real& m, unsigned short n);

  /// return basisPolyType
  short basis_type() const;

  /// returns polyRep for access to derived class member functions
  /// that are not mapped to the top BasisPolynomial level
  BasisPolynomial* polynomial_rep() const;
  /// function to check polyRep (does this handle contain a body)
  bool is_null() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  BasisPolynomial(BaseConstructor);

  //
  //- Heading: Data
  //

  /// basis polynomial type:
  /// {HERMITE,LEGENDRE,LAGUERRE,JACOBI,GEN_LAGUERRE,NUM_GEN_ORTHOG}_ORTHOG,
  /// {LAGRANGE,HERMITE}_INTERP, or PIECEWISE_{LINEAR,QUADRATIC,CUBIC}_INTERP
  short basisPolyType;

  /// weight discrepancy factor between Abramowitz-Stegun and PDF orthogonality
  Real wtFactor;
  /// point discrepancy factor between Abramowitz-Stegun and PDF orthogonality
  Real ptFactor;

private:

  //
  //- Heading: Member functions
  //

  /// Used by the envelope constructor to initialize polyRep to the
  /// appropriate derived type.
  BasisPolynomial* get_polynomial(short poly_type, short rule);

  //
  //- Heading: Data
  //

  /// pointer to the letter (initialized only for the envelope)
  BasisPolynomial* polyRep;
  /// number of objects sharing polyRep
  int referenceCount;
};


inline short BasisPolynomial::basis_type() const
{ return (polyRep) ? polyRep->basisPolyType : basisPolyType; }


inline BasisPolynomial* BasisPolynomial::polynomial_rep() const
{ return polyRep; }


inline bool BasisPolynomial::is_null() const
{ return (polyRep) ? false : true; }


/** This implementation is unprotected from overflow, but this should
    be fine for the polynomial orders that we would expect to
    encounter.  Whenever possible, orthogonal polynomial
    implementations should use factorial_ratio() or n_choose_k()
    instead of factorial() to avoid overflow. */
inline Real BasisPolynomial::factorial(unsigned short n)
{
  Real fact = 1.;
  for (unsigned short i=2; i<=n; i++)
    fact *= i;
  return fact;
}


/** This implementation sequences products in order to minimize the
    chances of overflow, and its use should be preferred to
    factorial() whenever possible. */
inline Real BasisPolynomial::
factorial_ratio(unsigned short num, unsigned short den)
{
  Real fact_ratio = 1.;
  if (num > den)
    for (unsigned short i=den+1; i<=num; i++)
      fact_ratio *= i;
  else if (den > num)
    for (unsigned short i=num+1; i<=den; i++)
      fact_ratio /= i;
  return fact_ratio;
}


/** This implementation sequences products in order to minimize the
    chances of overflow, and its use should be preferred to
    factorial() whenever possible. */
inline Real BasisPolynomial::n_choose_k(unsigned short n, unsigned short k)
{
  if (n < k) {
    PCerr << "Error: bad inputs resulting in negative factorial in "
	  << "BasisPolynomial::n_choose_k()." << std::endl;
    abort_handler(-1);
  }
  Real combination = 1.;
  unsigned short nmk = n - k;

  // For integer n,k
  if (k <= nmk) // loop from 0 to k-1
    for (unsigned short i=0; i<k; i++)
      combination *= (Real)(n-i)/(Real)(k-i);
  else          // loop from 0 to nmk-1
    for (unsigned short i=0; i<nmk; i++)
      combination *= (Real)(n-i)/(Real)(nmk-i);

  // For real n,k
  //combination = gsl_sf_gamma(n+1.)/gsl_sf_gamma(k+1.)/gsl_sf_gamma(nmk+1.);

  return combination;
}


/** This is the rising/upper factorial formulation of the Pochhammer
    symbol (m)_n. */
inline Real BasisPolynomial::pochhammer(const Real& m, unsigned short n)
{
  // For integer n
  Real poch = (n == 0) ? 1. : m;
  for (unsigned short i=1; i<n; i++)
    poch *= m+i;

  // For real n:
  //Real poch = gsl_sf_gamma(m+n)/gsl_sf_gamma(m);

  return poch;
}

} // namespace Pecos

#endif
