/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        NodalInterpPolyApproximation
//- Description:  Class for Nodal Interpolation Polynomial Approximation
//-               
//- Owner:        Mike Eldred

#ifndef NODAL_INTERP_POLY_APPROXIMATION_HPP
#define NODAL_INTERP_POLY_APPROXIMATION_HPP

#include "InterpPolyApproximation.hpp"
#include "InterpolationPolynomial.hpp"


namespace Pecos {


/// Derived approximation class for nodal interpolation polynomials
/// (global approximation interpolating function values and
/// potentially gradients at collocation points).

/** The NodalInterpPolyApproximation class provides a global polynomial
    approximation based on either Lagrange or Hermite interpolation
    polynomials using a nodal basis approach.  It is used primarily
    for stochastic collocation approaches to uncertainty quantification. */

class NodalInterpPolyApproximation: public InterpPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  NodalInterpPolyApproximation(short basis_type, size_t num_vars,
			       bool use_derivs, short output_level);
  /// destructor
  ~NodalInterpPolyApproximation();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_expansion_coefficients();
  void compute_expansion_coefficients();
  void increment_expansion_coefficients();
  void decrement_expansion_coefficients();
  void restore_expansion_coefficients();
  void finalize_expansion_coefficients();

  /// store current state within storedExpType{1Coeffs,2Coeffs,1CoeffGrads},
  /// storedColloc{Key,Indices}, and storedLev{MultiIndex,Coeffs}
  void store_coefficients();
  /// augment current interpolant using
  /// storedExpType{1Coeffs,2Coeffs,1CoeffGrads}, storedColloc{Key,Indices},
  /// and storedLev{MultiIndex,Coeffs}
  void combine_coefficients(short combine_type);

  void compute_numerical_response_moments(size_t num_moments);
  void compute_numerical_expansion_moments(size_t num_moments);

  Real value(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x,
					     const SizetArray& dvv);
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

  void compute_total_sobol_indices();
  void compute_partial_variance(const BitArray& set_value);

  const RealVector& approximation_coefficients() const;
  void approximation_coefficients(const RealVector& approx_coeffs);

  void set_new_point(const RealVector& x, const UShortArray& basis_index,
		     short order);
  void set_new_point(const RealVector& x, const UShortArray& basis_index,
		     const SizetList& subset_indices, short order);

  size_t barycentric_exact_index(const UShortArray& basis_index);
  size_t barycentric_exact_index(const UShortArray& basis_index,
				 const SizetList& subset_indices);

private:

  //
  //- Heading: Convenience functions
  //

  /// compute the mean of a tensor interpolant on a tensor grid;
  /// contributes to mean(x)
  Real tensor_product_mean(const RealVector& x, const UShortArray& lev_index,
    const UShort2DArray& key, const SizetArray& colloc_index);

  /// compute the gradient of the mean of a tensor interpolant on a
  /// tensor grid; contributes to mean_gradient(x)
  const RealVector& tensor_product_mean_gradient(const RealVector& x,
    const UShortArray& lev_index,   const UShort2DArray& key,
    const SizetArray& colloc_index, const SizetArray& dvv);

  /// compute the covariance of two tensor interpolants on the same
  /// tensor grid using an interpolation of products or product of
  /// interpolants approach; contributes to covariance(x, poly_approx_2)
  Real tensor_product_covariance(const RealVector& x,
    const UShortArray& lev_index,   const UShort2DArray& key,
    const SizetArray& colloc_index, NodalInterpPolyApproximation* nip_approx_2);
  /// compute the covariance of two tensor interpolants on different
  /// tensor grids using a product of interpolants approach;
  /// contributes to covariance(x, poly_approx_2)
  Real tensor_product_covariance(const RealVector& x,
    const UShortArray& lev_index_1, const UShort2DArray& key_1,
    const SizetArray& colloc_index_1, const UShortArray& lev_index_2,
    const UShort2DArray& key_2, const SizetArray& colloc_index_2,
    NodalInterpPolyApproximation* nip_approx_2);

  /// compute the gradient of the variance of a tensor interpolant on
  /// a tensor grid using an interpolation of products or product of
  /// interpolants approach; contributes to variance_gradient(x)
  const RealVector& tensor_product_variance_gradient(const RealVector& x,
    const UShortArray& lev_index, const UShort2DArray& key,
    const SizetArray& colloc_index, const SizetArray& dvv);

  /// update accumulators for barycentric type1 contributions to moment value
  void accumulate_barycentric(RealVector& t1_accumulator,
			      const UShortArray& lev_index,
			      const UShortArray& key_p);
  /// update accumulators for type1 contributions to moment value
  void accumulate_horners(RealVector& t1_accumulator,
			  const UShortArray& lev_index,
			  const UShortArray& key_p, const RealVector& x);
  /// update accumulators for type1 + type2 contributions to moment value
  void accumulate_horners(RealVector& t1_accumulator,
			  RealMatrix& t2_accumulator,
			  const UShortArray& lev_index,
			  const UShortArray& key_p, const RealVector& x);

  /// update accumulators for barycentric type1 contributions to moment gradient
  void accumulate_barycentric_gradient(RealMatrix& t1_accumulator,
				       const UShortArray& lev_index,
				       const UShortArray& key_p,
				       const SizetArray& dvv);
  /// update accumulators for type1 contributions to moment gradient
  void accumulate_horners_gradient(RealMatrix& t1_accumulator,
				   const UShortArray& lev_index,
				   const UShortArray& key_p,
				   const SizetArray& dvv, const RealVector& x);
  /// update accumulators for type1 + type2 contributions to moment gradient
  void accumulate_horners_gradient(RealMatrix& t1_accumulator,
				   RealMatrixArray& t2_accumulators,
				   const UShortArray& lev_index,
				   const UShortArray& key_p,
				   const SizetArray& dvv, const RealVector& x);

  /// update precomputation of nonzero multidimensional integrals of
  /// products of interpolation polynomials
  void update_nonzero_basis_products(const UShort2DArray& sm_multi_index);

  /// evaluate 1D integral of product of interpolation polynomials
  bool basis_product_1d(InterpolationPolynomial* poly_rep_1,
			InterpolationPolynomial* poly_rep_2,
			unsigned short key_1, unsigned short key_2,
			const RealArray& pts, const RealArray& wts, Real& prod);
  /// lookup multidimensional integral of products of interpolation polynomials
  bool basis_product(const UShortArray& lev_index_1, const UShortArray& key_1,
		     const UShortArray& lev_index_2, const UShortArray& key_2,
		     Real& prod);

  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using weights from the CombinedSparseGridDriver
  Real expectation(const RealVector& t1_coeffs, const RealMatrix& t2_coeffs);
  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using t{1,2}_wts
  Real expectation(const RealVector& t1_coeffs, const RealVector& t1_wts,
		   const RealMatrix& t2_coeffs, const RealMatrix& t2_wts);

  /// compute value of reduced-dimension interpolant
  Real value(const RealVector& x, const RealVectorArray& t1_coeffs,
	     const RealMatrixArray& t2_coeffs, const UShort3DArray& colloc_key,
	     const SizetList& subset_indices);
  /// compute gradient of reduced-dimension interpolant with respect
  /// to basis variables
  const RealVector& gradient_basis_variables(const RealVector& x,
					     const RealVectorArray& t1_coeffs,
					     const RealMatrixArray& t2_coeffs,
					     const UShort3DArray& colloc_key,
					     const SizetList& subset_indices);

  /// computes higher-order grid for tensor reinterpolation of the
  /// covariance fn for non-integrated dimensions in all_variables mode
  void reinterpolated_level(const UShortArray& lev_index);

  /// compute integral for total Sobol' index for variables in a set
  Real member_integral(const BitArray& member_bits, Real mean);
  /// defines member_coeffs and member_wts for a particular membership set
  void member_coefficients_weights(const BitArray& member_bits,
    const UShortArray& quad_order,    const UShortArray& lev_index,
    const UShort2DArray& colloc_key,  const SizetArray& colloc_index,
    RealVector& member_t1_coeffs,     RealVector& member_t1_wts,
    RealMatrix& member_t2_coeffs,     RealMatrix& member_t2_wts,
    UShort2DArray& member_colloc_key, SizetArray& member_colloc_index);
  /// create a unique map key for value() and gradient() calculation reuse
  void update_member_key(const UShortArray& data,
			 const SizetList& member_indices,
			 UShortArray& member_map_key, size_t cntr);

  //
  //- Heading: Data
  //

  /// type of interpolation for all-variables covariance and variance gradient
  short momentInterpType;

  /// the type1 coefficients of the expansion for interpolating values
  RealVector expansionType1Coeffs;
  /// the type2 coefficients of the expansion for interpolating gradients
  RealMatrix expansionType2Coeffs;
  /// the gradients of the type1 expansion coefficients
  /** may be interpreted as either the gradients of the expansion
      coefficients or the coefficients of expansions for the response
      gradients.  This array is used when sensitivities of moments are
      needed with respect to variables that do not appear in the
      expansion (e.g., with respect to design variables for an
      expansion only over the random variables). */
  RealMatrix expansionType1CoeffGrads;

  /// map from random index to unique nonZerosMapArray
  SizetArray nonZerosMapIndices;
  /// tracks level maxima already populated within nonZerosMap
  UShortArray nonZerosMapMaxLevels;
  /// expectations of products of interpolation polynomials,
  /// precomputed in update_nonzero_basis_products() for efficiency
  std::vector<UShort2DMultiSetRealMap> nonZerosMapArray;

  /// storage of expansionType1Coeffs state for subsequent restoration
  RealVector storedExpType1Coeffs;
  /// storage of expansionType2Coeffs state for subsequent restoration
  RealMatrix storedExpType2Coeffs;
  /// storage of expansionType1CoeffGrads state for subsequent restoration
  RealMatrix storedExpType1CoeffGrads;
  /// storage of level multi-index (levels for tensor or sparse grids)
  /// for subsequent restoration
  UShort2DArray storedLevMultiIndex;
  /// storage of IntegrationDriver combinatorial coefficients state
  /// for subsequent restoration
  IntArray storedLevCoeffs;
  /// storage of IntegrationDriver collocation key state for
  /// subsequent restoration
  UShort3DArray storedCollocKey;
  /// storage of IntegrationDriver collocation indices state for
  /// subsequent restoration
  Sizet2DArray storedCollocIndices;
};


inline NodalInterpPolyApproximation::
NodalInterpPolyApproximation(short basis_type, size_t num_vars,
			     bool use_derivs, short output_level):
  InterpPolyApproximation(basis_type, num_vars, use_derivs, output_level),
  // These 3 compile-time options are relevant for all-variables covariance
  // involving expectations over variable subsets.  Covariance for hierarchical
  // interpolants, nodal covariance in the standard view mode, uses of
  // PolynomialApproximation::compute_numerical_moments(), and Sobol' index
  // calculations all employ an INTERPOLATION_OF_PRODUCTS approach, so that
  // setting is the most self-consistent.  Gradient enhancement is also not
  // currently supported for PRODUCT_OF_INTERPOLANTS approaches.
  momentInterpType(INTERPOLATION_OF_PRODUCTS)
  //momentInterpType(REINTERPOLATION_OF_PRODUCTS)
  //momentInterpType(PRODUCT_OF_INTERPOLANTS_FULL)
  //momentInterpType(PRODUCT_OF_INTERPOLANTS_FAST)
{
  // We are migrating towards consistent usage of INTERPOLATION_OF_PRODUCTS,
  // but its usage of higher-order reinterpolation of covariance is currently
  // too slow for production usage.  Thus, we only activate it when needed to
  // support new capability, such as gradient-enhanced interpolation.
  //momentInterpType = (use_derivs) ?
  //  REINTERPOLATION_OF_PRODUCTS : PRODUCT_OF_INTERPOLANTS_FAST;
}


inline NodalInterpPolyApproximation::~NodalInterpPolyApproximation()
{ }


inline void NodalInterpPolyApproximation::increment_expansion_coefficients()
{ restore_expansion_coefficients(); }


inline void NodalInterpPolyApproximation::decrement_expansion_coefficients()
{
  // not necessary to prune; next increment/restore/finalize takes care of this
  //if (expConfigOptions.expansionCoeffFlag) {
  //  expansionType1Coeffs.resize(numCollocPts);
  //  if (basisConfigOptions.useDerivs) {
  //    size_t num_deriv_vars = expansionType2Coeffs.numRows();
  //    expansionType2Coeffs.reshape(num_deriv_vars, numCollocPts);
  //  }
  //}
  //if (expConfigOptions.expansionCoeffGradFlag) {
  //  size_t num_deriv_vars = expansionType1CoeffGrads.numRows();
  //  expansionType1CoeffGrads.reshape(num_deriv_vars, numCollocPts);
  //}
}


inline void NodalInterpPolyApproximation::finalize_expansion_coefficients()
{ restore_expansion_coefficients(); }


inline const RealVector& NodalInterpPolyApproximation::
approximation_coefficients() const
{
  if (basisConfigOptions.useDerivs) {
    PCerr << "Error: approximation_coefficients() not supported in "
	  << "InterpPolyApproximation for type2 coefficients." << std::endl;
    return abort_handler_t<const RealVector&>(-1);
  }
  else
    return expansionType1Coeffs;
}


inline void NodalInterpPolyApproximation::
approximation_coefficients(const RealVector& approx_coeffs)
{
  if (basisConfigOptions.useDerivs) {
    PCerr << "Error: approximation_coefficients() not supported in "
	  << "InterpPolyApproximation for type2 coefficients." << std::endl;
    abort_handler(-1);
  }
  else
    expansionType1Coeffs = approx_coeffs;
}


inline bool NodalInterpPolyApproximation::
basis_product_1d(InterpolationPolynomial* poly_rep_1,
		 InterpolationPolynomial* poly_rep_2,
		 unsigned short key_1, unsigned short key_2,
		 const RealArray& pts, const RealArray& wts, Real& prod)
{
  Real tol = 1.e-12; // consistent with OrthogonalPolynomial triple product tol
  prod = 0.;
  size_t i, num_pts = pts.size();
  for (i=0; i<num_pts; ++i)
    prod += wts[i] * poly_rep_1->type1_value(pts[i], key_1)
                   * poly_rep_2->type1_value(pts[i], key_2);
  return (std::abs(prod) > tol) ? true : false;
}


inline Real NodalInterpPolyApproximation::
expectation(const RealVector& t1_coeffs, const RealMatrix& t2_coeffs)
{
  // This version defaults to type1/2 weights from CombinedSparseGridDriver
  return expectation(t1_coeffs, driverRep->type1_weight_sets(),
		     t2_coeffs, driverRep->type2_weight_sets());
}


inline void NodalInterpPolyApproximation::
update_member_key(const UShortArray& data,
		  const SizetList&   member_indices,
		  UShortArray& member_map_key, size_t cntr)
{
  for (SizetList::const_iterator cit=member_indices.begin();
       cit!=member_indices.end(); ++cit, ++cntr)
    member_map_key[cntr] = data[*cit];
}


inline void NodalInterpPolyApproximation::
set_new_point(const RealVector& x, const UShortArray& basis_index, short order)
{
  unsigned short bi_j;
  for (size_t j=0; j<numVars; ++j) {
    bi_j = basis_index[j];
    if (bi_j) // exclusion of pt must be sync'd w/ factors/scalings
      polynomialBasis[bi_j][j].set_new_point(x[j], order);
  }
}


inline void NodalInterpPolyApproximation::
set_new_point(const RealVector& x, const UShortArray& basis_index,
	      const SizetList& subset_indices, short order)
{
  SizetList::const_iterator cit; size_t j; unsigned short bi_j;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j = *cit; bi_j = basis_index[j];
    if (bi_j) // exclusion of pt must be sync'd w/ factors/scalings
      polynomialBasis[bi_j][j].set_new_point(x[j], order);
  }
}


inline size_t NodalInterpPolyApproximation::
barycentric_exact_index(const UShortArray& basis_index)
{
  size_t j, pt_index = 0, prod = 1; unsigned short bi_j;
  for (j=0; j<numVars; ++j) {
    bi_j = basis_index[j];
    // Note: if bi_j == 0, then constant interp with 1 point: we can replace
    // this constant interpolation with the value at the 1 colloc index (ei=0)
    if (bi_j) {
      BasisPolynomial& poly_i = polynomialBasis[bi_j][j];
      pt_index += poly_i.exact_index() * prod;
      prod     *= poly_i.interpolation_size();
    }
  }
  return pt_index;
}


inline size_t NodalInterpPolyApproximation::
barycentric_exact_index(const UShortArray& basis_index,
			const SizetList& subset_indices)
{
  size_t j, pt_index = 0, prod = 1; unsigned short bi_j;
  SizetList::const_iterator cit;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j = *cit; bi_j = basis_index[j];
    // Note: if bi_j == 0, then constant interp with 1 point: we can replace
    // this constant interpolation with the value at the 1 colloc index (ei=0)
    if (bi_j) {
      BasisPolynomial& poly_j = polynomialBasis[bi_j][j];
      pt_index += poly_j.exact_index() * prod;
      prod     *= poly_j.interpolation_size();
    }
  }
  return pt_index;
}

} // namespace Pecos

#endif
