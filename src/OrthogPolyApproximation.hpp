/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        OrthogPolyApproximation
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        Mike Eldred

#ifndef ORTHOG_POLY_APPROXIMATION_HPP
#define ORTHOG_POLY_APPROXIMATION_HPP

#include "PolynomialApproximation.hpp"
#include "NumericGenOrthogPolynomial.hpp"
#include "SharedOrthogPolyApproxData.hpp"

namespace Pecos {


/// Derived approximation class for orthogonal polynomials (global
/// approximation).

/** The OrthogPolyApproximation class provides a global approximation
    based on orthogonal polynomials.  It is used primarily for polynomial
    chaos expansions (for stochastic finite element approaches to
    uncertainty quantification). */

class OrthogPolyApproximation: public PolynomialApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  OrthogPolyApproximation(const SharedBasisApproxData& shared_data);
  /// destructor
  ~OrthogPolyApproximation();

  //
  //- Heading: Virtual functions
  //

  /// retrieve number of terms in the orthogonal polynomial expansion
  virtual size_t expansion_terms() const;

  /// estimate chaos expansion coefficient decay rates for each random
  /// variable dimension using linear least squares in semilog space
  virtual const RealVector& dimension_decay_rates();

  /// retrieve or form a set of dense coefficients that correspond to
  /// SharedOrthogPolyApproxData::multiIndex
  virtual RealVector dense_coefficients() const;

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  int min_coefficients() const;

  void store_coefficients();
  void combine_coefficients(short combine_type);

  void print_coefficients(std::ostream& s, bool normalized = false);

  const RealVector& approximation_coefficients() const;
  void approximation_coefficients(const RealVector& approx_coeffs);

  void coefficient_labels(std::vector<std::string>& all_coeff_tags) const;

  /// initialize multiIndex, expansionCoeffs, et al.
  void allocate_arrays();

  /// Performs global sensitivity analysis via variance-based decomposition;
  /// computes component (main and interaction) Sobol' indices
  void compute_component_sobol();
  /// Performs global sensitivity analysis via variance-based decomposition;
  /// computes total Sobol' indices
  void compute_total_sobol();

  /// retrieve the response PCE value for a given parameter vector
  Real value(const RealVector& x);
  /// retrieve the gradient of the response PCE with respect to all
  /// variables included in the polynomial basis (e.g., probabilistic
  /// variables) for a given parameter vector
  const RealVector& gradient_basis_variables(const RealVector& x);
  /// retrieve the gradient of the response PCE with respect to variables
  /// included in the polynomial basis (e.g., probabilistic or "all"
  /// variables) for a given parameter vector and a given DVV subset
  const RealVector& gradient_basis_variables(const RealVector& x,
					     const SizetArray& dvv);
  /// retrieve the gradient of the response PCE with respect to variables not
  /// included in the polynomial basis (nonprobabilistic variables such as
  /// design or epistemic when not in "all" mode) for a given parameter vector
  const RealVector& gradient_nonbasis_variables(const RealVector& x);

  Real stored_value(const RealVector& x);
  const RealVector& stored_gradient_basis_variables(const RealVector& x);
  const RealVector& stored_gradient_nonbasis_variables(const RealVector& x);

  Real mean();
  Real mean(const RealVector& x);
  const RealVector& mean_gradient();
  const RealVector& mean_gradient(const RealVector& x, const SizetArray& dvv);

  Real variance();
  Real variance(const RealVector& x);
  const RealVector& variance_gradient();
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

  Real covariance(PolynomialApproximation* poly_approx_2);
  Real covariance(const RealVector& x, PolynomialApproximation* poly_approx_2);

  /// compute expansion moments to order 2
  void compute_moments();
  /// compute expansion moments in all-variables mode to order 2
  void compute_moments(const RealVector& x);

  /// return expansionMoments
  const RealVector& moments() const;

  //
  //- Heading: Member functions
  //

  /// size expansion{Coeffs,CoeffGrads} based on multiIndex
  void size_expansion();
  /// size expansion{Coeffs,CoeffGrads} based on multiIndex
  void size_expansion(size_t num_exp_terms);
  /// synchronize expansion{Coeffs,CoeffGrads} with an updated multiIndex
  void resize_expansion();
  /// synchronize expansion{Coeffs,CoeffGrads} with an updated multiIndex
  void resize_expansion(size_t num_exp_terms);

  /// overlay the passed expansion with the aggregate
  /// expansion{Coeffs,CoeffGrads} as managed by the multi_index_map
  void overlay_expansion(const SizetArray& multi_index_map,
			 const RealVector& exp_coeffs,
			 const RealMatrix& exp_grads, int coeff);
  /// multiply current expansion ("a") with incoming expansion ("b")
  /// and store in product expansion ("c")
  void multiply_expansion(const UShort2DArray& multi_index_b,
			  const RealVector& exp_coeffs_b,
			  const RealMatrix& exp_grads_b,
			  const UShort2DArray& multi_index_c);

  /// update add_val and add_gradient based on surrData's failure map
  void fail_booleans(SizetShortMap::const_iterator& fit, size_t j,
		     bool& add_val, bool& add_grad);

  /// perform sanity checks prior to numerical integration
  void integration_checks();

  /// utility function for solving the least squares estimation of decay rates
  void solve_decay_rates(RealVectorArray& A_vectors, RealVectorArray& b_vectors,
			 UShortArray& max_orders);

  //
  //- Heading: Data
  //

  /// the coefficients of the expansion
  RealVector expansionCoeffs;
  /// the gradients of the expansion coefficients
  /** may be interpreted as either the gradients of the expansion
      coefficients or the coefficients of expansions for the response
      gradients.  This array is used when sensitivities of moments are
      needed with respect to variables that do not appear in the
      expansion (e.g., with respect to design or epistemic variables
      for an expansion only over probabilistic variables). */
  RealMatrix expansionCoeffGrads;

  /// copy of expansionCoeffs (aggregated expansion) stored in
  /// store_coefficients() for use in combine_coefficients()
  RealVector storedExpCoeffs;
  /// copy of expansionCoeffGrads (aggregated expansion) stored in
  /// store_coefficients() for use in combine_coefficients()
  RealMatrix storedExpCoeffGrads;

  /// spectral coefficient decay rates estimated by LLS on log of
  /// univariate expansion coefficients
  RealVector decayRates;

private:

  //
  //- Heading: Data
  //

};


inline OrthogPolyApproximation::
OrthogPolyApproximation(const SharedBasisApproxData& shared_data):
  PolynomialApproximation(shared_data)
{ }


inline OrthogPolyApproximation::~OrthogPolyApproximation()
{ }


inline void OrthogPolyApproximation::compute_moments()
{
  // standard variables mode supports two expansion and four numerical moments
  mean(); variance(); // updates expansionMoments[0] and [1]
  //standardize_moments(expansionMoments);
}


inline void OrthogPolyApproximation::compute_moments(const RealVector& x)
{
  // all variables mode only supports first two moments
  mean(x); variance(x); // updates expansionMoments[0] and [1]
  //standardize_moments(expansionMoments);
}


inline const RealVector& OrthogPolyApproximation::moments() const
{ return expansionMoments; }


/** default implementation if no sparsity (overridden in
    RegressOrthogPolyApproximation for CS) */
inline size_t OrthogPolyApproximation::expansion_terms() const
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  return data_rep->multiIndex.size();
}


/** In this case, all expansion variables are random variables and the
    variance of the expansion is the sum over all but the first term
    of the coefficients squared times the polynomial norms squared. */
inline Real OrthogPolyApproximation::variance()
{ return covariance(this); }


/** In this case, a subset of the expansion variables are random variables
    and the variance of the expansion involves summations over this subset. */
inline Real OrthogPolyApproximation::variance(const RealVector& x)
{ return covariance(x, this); }


inline const RealVector& OrthogPolyApproximation::
approximation_coefficients() const
{ return expansionCoeffs; }


inline void OrthogPolyApproximation::
approximation_coefficients(const RealVector& approx_coeffs)
{
  expansionCoeffs = approx_coeffs;

  // allocate arrays in support of external coefficient import (mirrors
  // allocate_arrays() except for redundant size_expansion())
  allocate_total_sobol();
  allocate_component_sobol();
  if (expansionMoments.empty())
    expansionMoments.sizeUninitialized(2);
}


inline RealVector OrthogPolyApproximation::dense_coefficients() const
{
  // default implementation
  return RealVector(Teuchos::View, expansionCoeffs.values(),
		    expansionCoeffs.length());
}


inline void OrthogPolyApproximation::size_expansion()
{ size_expansion(expansion_terms()); }


inline void OrthogPolyApproximation::size_expansion(size_t num_exp_terms)
{
  if (expansionCoeffFlag && expansionCoeffs.length() != num_exp_terms)
    expansionCoeffs.sizeUninitialized(num_exp_terms);
  if (expansionCoeffGradFlag) {
    size_t num_deriv_vars = surrData.num_derivative_variables();
    if (expansionCoeffGrads.numRows() != num_deriv_vars ||
	expansionCoeffGrads.numCols() != num_exp_terms)
      expansionCoeffGrads.shapeUninitialized(num_deriv_vars, num_exp_terms);
  }
}


inline void OrthogPolyApproximation::resize_expansion()
{ resize_expansion(expansion_terms()); }


inline void OrthogPolyApproximation::resize_expansion(size_t num_exp_terms)
{
  if (expansionCoeffFlag)
    expansionCoeffs.resize(num_exp_terms); // new terms initialized to 0
  if (expansionCoeffGradFlag) {
    size_t num_deriv_vars = expansionCoeffGrads.numRows();
    expansionCoeffGrads.reshape(num_deriv_vars, num_exp_terms);
    // new terms initialized to 0
  }
}


inline void OrthogPolyApproximation::
fail_booleans(SizetShortMap::const_iterator& fit, size_t j,
	      bool& add_val, bool& add_grad)
{
  if (fit != surrData.failed_response_data().end() && fit->first == j) {
    short fail_bits = fit->second;
    if (fail_bits & 1) add_val  = false;
    if (fail_bits & 2) add_grad = false;
    ++fit;
  }
}

} // namespace Pecos

#endif
