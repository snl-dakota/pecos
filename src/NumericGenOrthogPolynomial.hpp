/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        NumericGenOrthogPolynomial
//- Description:  Class for numerically generated orthogonal polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef NUMERIC_GEN_ORTHOG_POLYNOMIAL_HPP
#define NUMERIC_GEN_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"
#include "pecos_global_defs.hpp"
#include "pecos_stat_util.hpp"


namespace Pecos {

/// pointer to a PDF evaluation function used within integral evaluators
typedef Real ( *NGFPType ) (const Real& x, const RealVector& params);
//typedef Real ( *NGFPType ) (const Real& x, const Real& p1, const Real& p2);


/// Derived orthogonal polynomial class for numerically-generated
/// orthogonal polynomials

/** The NumericGenOrthogPolynomial class numerically generates a
    univariate orthogonal polynomial of a particular order, along with
    its Gauss points, Gauss weights, and norms.  It uses a variety of
    algorithms due to Chebyshev and Stieltjes as reported by Golub and
    Welsch (Mathematics of Computation, Vol. 23, No. 106, 1969) and
    Gautschi (SIAM J. Sci. Stat. Comput., Vol. 3, No. 3, 1982).  It
    enables (mixed) multidimensional orthogonal polynomial basis
    functions within OrthogPolyApproximation. */

class NumericGenOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  NumericGenOrthogPolynomial();  ///< default constructor
  ~NumericGenOrthogPolynomial(); ///< destructor

  //
  //- Heading: Public functions
  //

  /// calculate and return alpha3TR[order]
  const Real& alpha_recursion(unsigned short order);
  /// calculate and return beta3TR[order]
  const Real& beta_recursion(unsigned short order);

  /// set distribution type and parameters for a BOUNDED_NORMAL distribution
  void bounded_normal_distribution(const Real& mean,  const Real& std_dev,
				   const Real& l_bnd, const Real& u_bnd);
  /// set distribution type and parameters for a LOGNORMAL distribution
  void lognormal_distribution(const Real& mean, const Real& std_dev);
  /// set distribution type and parameters for a BOUNDED_LOGNORMAL distribution
  void bounded_lognormal_distribution(const Real& mean,  const Real& std_dev,
				      const Real& l_bnd, const Real& u_bnd);
  /// set distribution type and parameters for a LOGUNIFORM distribution
  void loguniform_distribution(const Real& l_bnd, const Real& u_bnd);
  /// set distribution type and parameters for a TRIANGULAR distribution
  void triangular_distribution(const Real& mode, const Real& l_bnd,
			       const Real& u_bnd);
  /// set distribution type and parameters for a GUMBEL distribution
  void gumbel_distribution(const Real& alpha, const Real& beta);
  /// set distribution type and parameters for a FRECHET distribution
  void frechet_distribution(const Real& alpha, const Real& beta);
  /// set distribution type and parameters for a WEIBULL distribution
  void weibull_distribution(const Real& alpha, const Real& beta);
  /// set distribution type and parameters for a WEIBULL distribution
  void histogram_bin_distribution(const RealVector& bin_pairs);

  /// set coeffsNormsFlag
  void coefficients_norms_flag(bool flag);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  Real type1_value(const Real& x, unsigned short order);
  Real type1_gradient(const Real& x, unsigned short order);
  Real norm_squared(unsigned short order);

  const RealArray& collocation_points(unsigned short order);
  const RealArray& type1_collocation_weights(unsigned short order);

private:

  //
  //- Heading: Convenience functions
  //

  /// thin wrapper for bounded_normal_pdf for NGFPType API
  static Real bounded_normal_pdf(const Real& x, const RealVector& params);
  /// thin wrapper for lognormal_pdf for NGFPType API
  static Real lognormal_pdf(const Real& x, const RealVector& params);
  /// thin wrapper for bounded_lognormal_pdf for NGFPType API
  static Real bounded_lognormal_pdf(const Real& x, const RealVector& params);
  /// thin wrapper for loguniform_pdf for NGFPType API
  static Real loguniform_pdf(const Real& x, const RealVector& params);
  /// thin wrapper for triangular_pdf for NGFPType API
  static Real triangular_pdf(const Real& x, const RealVector& params);
  /// thin wrapper for gumbel_pdf for NGFPType API
  static Real gumbel_pdf(const Real& x, const RealVector& params);
  /// thin wrapper for frechet_pdf for NGFPType API
  static Real frechet_pdf(const Real& x, const RealVector& params);
  /// thin wrapper for weibull_pdf for NGFPType API
  static Real weibull_pdf(const Real& x, const RealVector& params);

  /// solve a symmetric tridiagonal eigenvalue problem for the Gauss
  /// points and weights for an orthogonal polynomial of order m
  void solve_eigenproblem(unsigned short m);

  /// compute three point recursion for polyCoeffs[i+1]
  void polynomial_recursion(RealVector& poly_coeffs_ip1, const Real& alpha_i,
			    const RealVector& poly_coeffs_i, const Real& beta_i,
			    const RealVector& poly_coeffs_im1);
  /// compute truncated three point recursion for polyCoeffs[i+1]
  void polynomial_recursion(RealVector& poly_coeffs_ip1, const Real& alpha_i,
			    const RealVector& poly_coeffs_i);

  /// compute inner product of specified polynomial orders
  Real inner_product(const RealVector& poly_coeffs1,
		     const RealVector& poly_coeffs2);

  /// compute an unbounded integral using Gauss-Hermite integration
  Real hermite_unbounded_integral(const RealVector& poly_coeffs1,
				  const RealVector& poly_coeffs2,
				  NGFPType weight_fn);
  /// compute an unbounded integral using Fejer integration and a
  /// change of variables
  Real fejer_unbounded_integral(const RealVector& poly_coeffs1,
				const RealVector& poly_coeffs2,
				NGFPType weight_fn, unsigned short quad_order);
  /// compute a semibounded integral using Gauss-Laguerre integration
  Real laguerre_semibounded_integral(const RealVector& poly_coeffs1,
				     const RealVector& poly_coeffs2,
				     NGFPType weight_fn);
  /// compute a semibounded integral using Fejer integration and a
  /// change of variables
  Real fejer_semibounded_integral(const RealVector& poly_coeffs1,
				  const RealVector& poly_coeffs2,
				  NGFPType weight_fn,
				  unsigned short quad_order);
  /// compute a bounded integral over the specified range using
  /// Gauss-Legendre integration
  Real legendre_bounded_integral(const RealVector& poly_coeffs1,
				 const RealVector& poly_coeffs2,
				 NGFPType weight_fn, Real start, Real end);
  /// compute a bounded integral over the specified range using
  /// Clenshaw-Curtis integration
  Real cc_bounded_integral(const RealVector& poly_coeffs1,
			   const RealVector& poly_coeffs2, NGFPType weight_fn,
			   Real start, Real end, unsigned short quad_order);
  /// compute a bounded integral over the specified range using Riemann sums
  Real riemann_bounded_integral(const RealVector& poly_coeffs1,
				const RealVector& poly_coeffs2,
				NGFPType weight_fn, Real start, Real end);
  /// compute an integral using the native Gaussian quadrature rule
  /// (up to order 2m-1 based on collocPoints and collocWeights of order m)
  Real native_quadrature_integral(const RealVector& poly_coeffs1,
				  const RealVector& poly_coeffs2);

  /// retrieve the value of the 1-D generated polynomial (of given
  /// coefficients) for a given parameter value
  Real type1_value(const Real& x, const RealVector& poly_coeffs);
  /// retrieve the gradient of the 1-D generated polynomial (of given
  /// coefficients) with respect to its dimension for a given parameter value
  Real type1_gradient(const Real& x, const RealVector& poly_coeffs);

  //
  //- Heading: Data
  //

  /// the type of non-Askey distribution: BOUNDED_NORMAL, LOGNORMAL,
  /// BOUNDED_LOGNORMAL, LOGUNIFORM, TRIANGULAR, GUMBEL, FRECHET,
  /// WEIBULL, or STOCHASTIC_EXPANSION
  short distributionType;

  /// distribution parameters (e.g., mean, std_dev, alpha, beta)
  RealVector distParams;

  /// flag identifying the need to compute polyCoeffs and orthogPolyNormsSq
  /// (if false, only collocPoints and collocWeights are computed)
  bool coeffsNormsFlag;

  /// coefficients of the orthogonal polynomials, from order 0 to m
  RealVectorArray polyCoeffs;

  /// alpha three-term recurrence parameters: alpha3TR[i] multiplied
  /// by polyCoeffs[i] contributes to polyCoeffs[i+1]
  RealVector alpha3TR;
  /// beta three-term recurrence parameters: beta3TR[i] multiplied
  /// by polyCoeffs[i-1] contributes to polyCoeffs[i+1]
  RealVector beta3TR;

  /// norm-squared of all orthogonal polynomials, from order 0 to m,
  /// as defined by the inner product <Poly_i, Poly_i> = ||Poly_i||^2
  RealVector orthogPolyNormsSq;
};


inline NumericGenOrthogPolynomial::NumericGenOrthogPolynomial() :
  coeffsNormsFlag(false)
{ collocRule = GOLUB_WELSCH; ptFactor = wtFactor = 1.; }


inline NumericGenOrthogPolynomial::~NumericGenOrthogPolynomial()
{ }


inline void NumericGenOrthogPolynomial::
polynomial_recursion(RealVector& poly_coeffs_ip1, const Real& alpha_i,
		     const RealVector& poly_coeffs_i)
{
  // compute alpha[i] recursion contribution to polyCoeffs[i+1]:
  int i_len = poly_coeffs_i.length();
  poly_coeffs_ip1.size(i_len+1); // initialize to zero
  for (int j=0; j<i_len; ++j) {
    poly_coeffs_ip1[j]   -= alpha_i*poly_coeffs_i[j]; // -alpha_i * poly_i
    poly_coeffs_ip1[j+1] += poly_coeffs_i[j];         // xi * poly_i
  }
}


inline void NumericGenOrthogPolynomial::
polynomial_recursion(RealVector& poly_coeffs_ip1, const Real& alpha_i,
		     const RealVector& poly_coeffs_i, const Real& beta_i,
		     const RealVector& poly_coeffs_im1)
{
  // compute alpha[i] recursion contribution to polyCoeffs[i+1]:
  polynomial_recursion(poly_coeffs_ip1, alpha_i, poly_coeffs_i);
  // compute beta[i] recursion contribution to polyCoeffs[i+1]:
  int im1_len = poly_coeffs_im1.length();
  for (int j=0; j<im1_len; ++j)
    poly_coeffs_ip1[j] -= beta_i*poly_coeffs_im1[j]; // -beta_i * poly_{i-1}
}


inline void NumericGenOrthogPolynomial::
bounded_normal_distribution(const Real& mean,  const Real& std_dev,
			    const Real& l_bnd, const Real& u_bnd)
{
  distributionType = BOUNDED_NORMAL;
  distParams.sizeUninitialized(4);
  distParams[0] = mean;  distParams[1] = std_dev;
  distParams[2] = l_bnd; distParams[3] = u_bnd;
  reset_gauss();
}


inline void NumericGenOrthogPolynomial::
lognormal_distribution(const Real& mean, const Real& std_dev)
{
  distributionType = LOGNORMAL;
  distParams.sizeUninitialized(2);
  distParams[0] = mean; distParams[1] = std_dev;
  reset_gauss();
}


inline void NumericGenOrthogPolynomial::
bounded_lognormal_distribution(const Real& mean,  const Real& std_dev,
			       const Real& l_bnd, const Real& u_bnd)
{
  distributionType = BOUNDED_LOGNORMAL;
  distParams.sizeUninitialized(4);
  distParams[0] = mean;  distParams[1] = std_dev;
  distParams[2] = l_bnd; distParams[3] = u_bnd;
  reset_gauss();
}


inline void NumericGenOrthogPolynomial::
loguniform_distribution(const Real& l_bnd, const Real& u_bnd)
{
  distributionType = LOGUNIFORM;
  distParams.sizeUninitialized(2);
  distParams[0] = l_bnd; distParams[1] = u_bnd;
  reset_gauss();
}


inline void NumericGenOrthogPolynomial::
triangular_distribution(const Real& mode, const Real& l_bnd, const Real& u_bnd)
{
  distributionType = TRIANGULAR;
  distParams.sizeUninitialized(3);
  distParams[0] = mode; distParams[1] = l_bnd; distParams[2] = u_bnd;
  reset_gauss();
}


inline void NumericGenOrthogPolynomial::
gumbel_distribution(const Real& alpha, const Real& beta)
{
  distributionType = GUMBEL;
  distParams.sizeUninitialized(2);
  distParams[0] = alpha; distParams[1] = beta;
  reset_gauss();
}


inline void NumericGenOrthogPolynomial::
frechet_distribution(const Real& alpha, const Real& beta)
{
  distributionType = FRECHET;
  distParams.sizeUninitialized(2);
  distParams[0] = alpha; distParams[1] = beta;
  reset_gauss();
}


inline void NumericGenOrthogPolynomial::
weibull_distribution(const Real& alpha, const Real& beta)
{
  distributionType = WEIBULL;
  distParams.sizeUninitialized(2);
  distParams[0] = alpha; distParams[1] = beta;
  reset_gauss();
}


inline void NumericGenOrthogPolynomial::
histogram_bin_distribution(const RealVector& bin_pairs)
{
  distributionType = HISTOGRAM_BIN;
  distParams = bin_pairs;
  reset_gauss();
}


inline Real NumericGenOrthogPolynomial::
bounded_normal_pdf(const Real& x, const RealVector& params)
{
  return Pecos::bounded_normal_pdf(x, params[0], params[1],
				   params[2], params[3]);
}


inline Real NumericGenOrthogPolynomial::
lognormal_pdf(const Real& x, const RealVector& params)
{ return Pecos::lognormal_pdf(x, params[0], params[1]); }


inline Real NumericGenOrthogPolynomial::
bounded_lognormal_pdf(const Real& x, const RealVector& params)
{
  return Pecos::bounded_lognormal_pdf(x, params[0], params[1],
				      params[2], params[3]);
}


inline Real NumericGenOrthogPolynomial::
loguniform_pdf(const Real& x, const RealVector& params)
{ return Pecos::loguniform_pdf(x, params[0], params[1]); }


inline Real NumericGenOrthogPolynomial::
triangular_pdf(const Real& x, const RealVector& params)
{ return Pecos::triangular_pdf(x, params[0], params[1], params[2]); }


inline Real NumericGenOrthogPolynomial::
gumbel_pdf(const Real& x, const RealVector& params)
{ return Pecos::gumbel_pdf(x, params[0], params[1]); }


inline Real NumericGenOrthogPolynomial::
frechet_pdf(const Real& x, const RealVector& params)
{ return Pecos::frechet_pdf(x, params[0], params[1]); }


inline Real NumericGenOrthogPolynomial::
weibull_pdf(const Real& x, const RealVector& params)
{ return Pecos::weibull_pdf(x, params[0], params[1]); }


inline void NumericGenOrthogPolynomial::coefficients_norms_flag(bool flag)
{ coeffsNormsFlag = flag; }

} // namespace Pecos

#endif
