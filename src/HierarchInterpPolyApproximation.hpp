/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchInterpPolyApproximation
//- Description:  Class for Nodal Interpolation Polynomial Approximation
//-               
//- Owner:        Mike Eldred

#ifndef HIERARCH_INTERP_POLY_APPROXIMATION_HPP
#define HIERARCH_INTERP_POLY_APPROXIMATION_HPP

#include "InterpPolyApproximation.hpp"
#include "SharedHierarchInterpPolyApproxData.hpp"
#include "HierarchSparseGridDriver.hpp"

namespace Pecos {


/// Derived approximation class for hierarchical interpolation polynomials
/// (interpolating values and potentially gradients at collocation points).

/** The HierarchInterpPolyApproximation class provides a polynomial
    approximation based on hierarchical interpolation.  Both local and
    global hierarchical basis functions are available.  It is used
    primarily for stochastic collocation approaches to uncertainty
    quantification. */

class HierarchInterpPolyApproximation: public InterpPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// Default constructor
  HierarchInterpPolyApproximation(const SharedBasisApproxData& shared_data);
  /// destructor
  ~HierarchInterpPolyApproximation();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// size expansionType{1,2}Coeffs and expansionType1CoeffGrads
  void allocate_arrays();

  void compute_coefficients();

  /// update the coefficients for the expansion of interpolation polynomials:
  /// increment expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void increment_coefficients();
  /// restore the coefficients to their previous state prior to last increment:
  /// decrement expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void decrement_coefficients(bool save_data);
  /// restore the coefficients to a previously incremented state as
  /// identified by the current increment to the Smolyak multi index:
  /// push expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void push_coefficients();
  /// finalize the coefficients by applying all previously evaluated increments:
  /// finalize expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void finalize_coefficients();

  /// clear inactive expansions from expansionType{1Coeffs,2Coeffs,1CoeffGrads}
  void clear_inactive();

  /// update combinedExpT{1Coeffs,2Coeffs,1CoeffGrads}
  void combine_coefficients();

  /// replace active expansions with combinedExpT{1Coeffs,2Coeffs,1CoeffGrads}
  void combined_to_active();

  void integrate_response_moments(size_t num_moments);
  void integrate_expansion_moments(size_t num_moments);

  Real value(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x,
					     const SizetArray& dvv);
  const RealVector& gradient_nonbasis_variables(const RealVector& x);
  const RealSymMatrix& hessian_basis_variables(const RealVector& x);

  Real stored_value(const RealVector& x, const UShortArray& key);
  const RealVector& stored_gradient_basis_variables(const RealVector& x,
						    const UShortArray& key);
  const RealVector& stored_gradient_basis_variables(const RealVector& x,
						    const SizetArray& dvv,
						    const UShortArray& key);
  const RealVector& stored_gradient_nonbasis_variables(const RealVector& x,
						       const UShortArray& key);
  const RealSymMatrix& stored_hessian_basis_variables(const RealVector& x,
						      const UShortArray& key);

  Real mean();
  Real mean(const RealVector& x);
  
  const RealVector& mean_gradient();
  const RealVector& mean_gradient(const RealVector& x,
				  const SizetArray& dvv);

  Real variance();
  Real variance(const RealVector& x);
  const RealVector& variance_gradient();
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

  Real covariance(PolynomialApproximation* poly_approx_2);
  Real covariance(const RealVector& x,
		  PolynomialApproximation* poly_approx_2);
  Real combined_covariance(PolynomialApproximation* poly_approx_2);
  Real combined_covariance(const RealVector& x,
			   PolynomialApproximation* poly_approx_2);

  Real delta_covariance(PolynomialApproximation* poly_approx_2);
  Real delta_covariance(const RealVector& x,
			PolynomialApproximation* poly_approx_2);
  Real delta_combined_covariance(PolynomialApproximation* poly_approx_2);
  Real delta_combined_covariance(const RealVector& x,
				 PolynomialApproximation* poly_approx_2);

  Real delta_mean();
  Real delta_mean(const RealVector& x);
  Real delta_std_deviation();
  Real delta_std_deviation(const RealVector& x);
  Real delta_beta(bool cdf_flag, Real z_bar);
  Real delta_beta(const RealVector& x, bool cdf_flag, Real z_bar);
  Real delta_z(bool cdf_flag, Real beta_bar);
  Real delta_z(const RealVector& x, bool cdf_flag, Real beta_bar);

  void compute_total_sobol_indices();
  void compute_partial_variance(const BitArray& set_value);

private:

  //
  //- Heading: Convenience functions
  //

  /// update {expT1Coeffs,expT2Coeffs,expT1CoeffGrads}Iter for new
  /// activeKey from sharedDataRep
  void update_active_iterators();

  /// compute the value at a point for a particular interpolation level
  Real value(const RealVector& x, const UShort3DArray& sm_mi,
	     const UShort4DArray& key, const RealVector2DArray& t1_coeffs,
	     const RealMatrix2DArray& t2_coeffs, unsigned short level);
  /// compute the value at a point for a particular interpolation
  /// level and for a specified subset of the variables
  Real value(const RealVector& x, const UShort3DArray& sm_mi,
	     const UShort4DArray& key, const RealVector2DArray& t1_coeffs,
	     const RealMatrix2DArray& t2_coeffs, unsigned short level,
	     const SizetList& subset_indices);

  /// compute the approximate gradient with respect to the basis variables
  /// at a particular point for a particular interpolation level
  const RealVector& gradient_basis_variables(const RealVector& x,
    const UShort3DArray& sm_mi, const UShort4DArray& key,
    const RealVector2DArray& t1_coeffs, const RealMatrix2DArray& t2_coeffs,
    unsigned short level);
  /// compute the approximate gradient with respect to the basis variables
  /// at a particular point for a particular interpolation level
  const RealVector& gradient_basis_variables(const RealVector& x,
    const UShort3DArray& sm_mi, const UShort4DArray& key,
    const RealVector2DArray& t1_coeffs, const RealMatrix2DArray& t2_coeffs,
    unsigned short level, const SizetList& subset_indices);
  /// compute the approximate gradient with respect to the basis variables
  /// for a particular point, interpolation level, and DVV
  const RealVector& gradient_basis_variables(const RealVector& x,
    const UShort3DArray& sm_mi, const UShort4DArray& key,
    const RealVector2DArray& t1_coeffs, const RealMatrix2DArray& t2_coeffs,
    const SizetArray& dvv, unsigned short level);
  /// compute the approximate gradient with respect to the nonbasis
  /// variables at a particular point for a particular interpolation level
  const RealVector& gradient_nonbasis_variables(const RealVector& x,
    const UShort3DArray& sm_mi, const UShort4DArray& key,
    const RealMatrix2DArray& t1_coeff_grads, unsigned short level);
  /// compute the approximate Hessian with respect to the basis variables
  /// at a particular point for a particular interpolation level
  const RealSymMatrix& hessian_basis_variables(const RealVector& x,
    const UShort3DArray& sm_mi,	const UShort4DArray& colloc_key,
    const RealVector2DArray& t1_coeffs, unsigned short level);

  /// update bookkeeping when adding a grid increment relative to the
  /// grid reference
  void increment_current_from_reference();
  /// update bookkeeping when removing a grid increment and returning
  /// to the grid reference
  void decrement_current_to_reference();

  /// compute the reference mean, excluding the current grid
  /// increment, using ref_key
  Real reference_mean(const UShort2DArray& ref_key);
  /// compute the reference mean, excluding the current grid
  /// increment, using ref_key
  Real reference_mean(const RealVector& x, const UShort2DArray& ref_key);
  /// compute the reference variance, excluding the current grid
  /// increment, using ref_key
  Real reference_variance(const UShort2DArray& ref_key);
  /// compute the reference variance, excluding the current grid
  /// increment, using ref_key
  Real reference_variance(const RealVector& x, const UShort2DArray& ref_key);

  /// compute the covariance increment due to the current grid increment
  Real delta_covariance(const RealVector2DArray& r1_t1_coeffs,
			const RealMatrix2DArray& r1_t2_coeffs,
			const RealVector2DArray& r2_t1_coeffs,
			const RealMatrix2DArray& r2_t2_coeffs, bool same,
			const RealVector2DArray& r1r2_t1_coeffs,
			const RealMatrix2DArray& r1r2_t2_coeffs,
			const RealVector2DArray& t1_wts,
			const RealMatrix2DArray& t2_wts,
			const UShort2DArray& ref_key,
			const UShort2DArray& incr_key);
  /// compute the covariance increment at x due to the current grid increment
  Real delta_covariance(const RealVector& x,
			const RealVector2DArray& r1_t1_coeffs,
			const RealMatrix2DArray& r1_t2_coeffs,
			const RealVector2DArray& r2_t1_coeffs,
			const RealMatrix2DArray& r2_t2_coeffs, bool same,
			const RealVector2DArray& r1r2_t1_coeffs,
			const RealMatrix2DArray& r1r2_t2_coeffs,
			const UShort3DArray& sm_mi,
			const UShort4DArray& colloc_key,
			const UShort2DArray& ref_key,
			const UShort2DArray& incr_key);

  /// compute the mean increment due to the current grid increment
  Real delta_mean(const UShort2DArray& incr_key);
  /// compute the mean increment due to the current grid increment
  Real delta_mean(const RealVector& x, const UShort2DArray& incr_key);
  /// compute the variance increment due to the current grid increment
  Real delta_variance(const UShort2DArray& ref_key,
		      const UShort2DArray& incr_key);
  /// compute the variance increment due to the current grid increment
  Real delta_variance(const RealVector& x, const UShort2DArray& ref_key,
		      const UShort2DArray& incr_key);
  /// compute the standard deviation increment due to the current grid increment
  Real delta_std_deviation(const UShort2DArray& ref_key,
			   const UShort2DArray& incr_key);
  /// compute the standard deviation increment due to the current grid increment
  Real delta_std_deviation(const RealVector& x, const UShort2DArray& ref_key,
			   const UShort2DArray& incr_key);
  /// compute the reliability index increment due to the current grid increment
  Real delta_beta(bool cdf_flag, Real z_bar, const UShort2DArray& ref_key,
		  const UShort2DArray& incr_key);
  /// compute the reliability index increment due to the current grid increment
  Real delta_beta(const RealVector& x, bool cdf_flag, Real z_bar,
		  const UShort2DArray& ref_key, const UShort2DArray& incr_key);
  /// shared logic for handling exceptional cases
  Real delta_beta_map(Real mu0, Real delta_mu, Real var0, Real delta_sigma,
		      bool cdf_flag, Real z_bar);
  /// compute the response level increment due to the current grid increment
  Real delta_z(bool cdf_flag, Real beta_bar, const UShort2DArray& ref_key,
	       const UShort2DArray& incr_key);
  /// compute the response level increment due to the current grid increment
  Real delta_z(const RealVector& x, bool cdf_flag, Real beta_bar,
	       const UShort2DArray& ref_key, const UShort2DArray& incr_key);

  /// form type 1/2 coefficients for interpolation of R_1 R_2
  void product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
    RealVector2DArray& r1r2_t1_coeffs, RealMatrix2DArray& r1r2_t2_coeffs,
    const UShort2DArray& reference_key = UShort2DArray());
  /// form type 1/2 coefficients for interpolation of R_1 R_2
  void product_interpolant(const RealMatrix2DArray& var_sets,
    const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
    const RealVector2DArray& r1_t1_coeffs,
    const RealMatrix2DArray& r1_t2_coeffs,
    const RealVector2DArray& r2_t1_coeffs,
    const RealMatrix2DArray& r2_t2_coeffs, bool same,
    RealVector2DArray& r1r2_t1_coeffs, RealMatrix2DArray& r1r2_t2_coeffs,
    const UShort2DArray& reference_key = UShort2DArray());

  /// form type 1/2 coefficients for interpolation of (R_1 - mu_1)(R_2 - mu_2)
  void central_product_interpolant(
    HierarchInterpPolyApproximation* hip_approx_2, Real mean_1, Real mean_2,
    RealVector2DArray& cov_t1_coeffs, RealMatrix2DArray& cov_t2_coeffs,
    const UShort2DArray& reference_key = UShort2DArray());
  /// form type 1/2 coefficients for interpolation of (R_1 - mu_1)(R_2 - mu_2)
  void central_product_interpolant(const RealMatrix2DArray& var_sets,
    const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
    const RealVector2DArray& r1_t1_coeffs,
    const RealMatrix2DArray& r1_t2_coeffs,
    const RealVector2DArray& r2_t1_coeffs,
    const RealMatrix2DArray& r2_t2_coeffs, bool same, Real mean_1, Real mean_2,
    RealVector2DArray& cov_t1_coeffs, RealMatrix2DArray& cov_t2_coeffs,
    const UShort2DArray& reference_key = UShort2DArray());

  /// form type1 coefficient gradients for interpolation of 
  /// d/ds [(R_1 - mu_1)(R_2 - mu_2)]
  void central_product_gradient_interpolant(
    HierarchInterpPolyApproximation* hip_approx_2, Real mean_1, Real mean_2,
    const RealVector& mean1_grad, const RealVector& mean2_grad, 
    RealMatrix2DArray& cov_t1_coeff_grads,
    const UShort2DArray& reference_key = UShort2DArray());
  /// form type1 coefficient gradients for interpolation of 
  /// d/ds [(R_1 - mu_1)(R_2 - mu_2)]
  void central_product_gradient_interpolant(const RealMatrix2DArray& var_sets,
    const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
    const RealVector2DArray& r1_t1_coeffs,
    const RealMatrix2DArray& r1_t2_coeffs,
    const RealMatrix2DArray& r1_t1_coeff_grads,
    const RealVector2DArray& r2_t1_coeffs,
    const RealMatrix2DArray& r2_t2_coeffs,
    const RealMatrix2DArray& r2_t1_coeff_grads, bool same, Real mean_1,
    Real mean_2, const RealVector& mean1_grad, const RealVector& mean2_grad, 
    RealMatrix2DArray& cov_t1_coeff_grads,
    const UShort2DArray& reference_key = UShort2DArray());

  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using active weights from the HierarchSparseGridDriver
  Real expectation(const RealVector2DArray& t1_coeffs,
		   const RealMatrix2DArray& t2_coeffs,
		   const UShort2DArray& set_partition = UShort2DArray());
  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using t{1,2}_wts
  Real expectation(const RealVector2DArray& t1_coeffs,
		   const RealVector2DArray& t1_wts,
		   const RealMatrix2DArray& t2_coeffs,
		   const RealMatrix2DArray& t2_wts,
		   const UShort2DArray& set_partition = UShort2DArray());
  // compute the expected value of the interpolant given by t{1,2}_coeffs
  //Real expectation(const RealVector2DArray& t1_coeffs,
  //		   const RealMatrix2DArray& t2_coeffs,
  //		   const UShort3DArray& pt_partition);

  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using active weights from the HierarchSparseGridDriver
  Real expectation(const RealVector& x, const RealVector2DArray& t1_coeffs,
		   const RealMatrix2DArray& t2_coeffs,
		   const UShort2DArray& set_partition = UShort2DArray());
  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using partial weights determined from sm_mi and colloc_key
  Real expectation(const RealVector& x, const RealVector2DArray& t1_coeffs,
		   const RealMatrix2DArray& t2_coeffs,
		   const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
		   const UShort2DArray& set_partition = UShort2DArray());

  /// compute the expected value of the gradient interpolant given by
  /// t1_coeff_grads using weights from the HierarchSparseGridDriver
  const RealVector& expectation_gradient(
    const RealMatrix2DArray& t1_coeff_grads);
  /// compute the expected value of the gradient interpolant given by
  /// t1_coeff_grads using t1_wts
  const RealVector& expectation_gradient(
    const RealMatrix2DArray& t1_coeff_grads, const RealVector2DArray& t1_wts);

  /// compute the expectation of t1_coeff_grads for index t1cg_index
  Real expectation_gradient(const RealVector& x,
			    const RealMatrix2DArray& t1_coeff_grads,
			    size_t t1cg_index);
  /// compute the expectation of t1_coeff_grads for index t1cg_index
  Real expectation_gradient(const RealVector& x,
			    const RealMatrix2DArray& t1_coeff_grads,
			    const UShort3DArray& sm_mi,
			    const UShort4DArray& colloc_key, size_t t1cg_index);

  /// compute the expectation of the gradient of {t1,t2}_coeffs for
  /// index deriv_index
  Real expectation_gradient(const RealVector& x,
			    const RealVector2DArray& t1_coeffs,
			    const RealMatrix2DArray& t2_coeffs,
			    size_t deriv_index);
  /// compute the expectation of the gradient of {t1,t2}_coeffs for
  /// index deriv_index
  Real expectation_gradient(const RealVector& x,
			    const RealVector2DArray& t1_coeffs,
			    const RealMatrix2DArray& t2_coeffs,
			    const UShort3DArray& sm_mi,
			    const UShort4DArray& colloc_key,
			    size_t deriv_index);

  /// increment expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  /// for a single index_set
  void increment_coefficients(const UShortArray& index_set);

  /// move the expansion coefficients for push_set from
  /// poppedExp{T1Coeffs,T2Coeffs,T1CoeffGrads} to
  /// expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void push_coefficients(const UShortArray& push_set);

  /// compute member expansion for Sobol' index integration
  void member_coefficients_weights(const BitArray&    member_bits,
    RealVector2DArray& member_t1_coeffs, RealVector2DArray& member_t1_wts,
    RealMatrix2DArray& member_t2_coeffs, RealMatrix2DArray& member_t2_wts,
    UShort4DArray& member_colloc_key,    Sizet3DArray& member_colloc_index);
  /// form hierarchical interpolant of (h-mean)^2 from member-variable
  /// expansion of h
  void central_product_member_coefficients(const BitArray& m_bits,
    const RealVector2DArray& m_t1_coeffs, const RealMatrix2DArray& m_t2_coeffs,
    const UShort4DArray&    m_colloc_key, const Sizet3DArray&   m_colloc_index,
    Real mean, RealVector2DArray& cprod_m_t1_coeffs,
    RealMatrix2DArray& cprod_m_t2_coeffs);

  //
  //- Heading: Data
  //

  /// bookkeeping to track computation of reference mean to avoid
  /// unnecessary recomputation
  short computedRefMean;
  /// bookkeeping to track computation of mean increment to avoid
  /// unnecessary recomputation
  short computedDeltaMean;
  /// bookkeeping to track computation of reference variance to avoid
  /// unnecessary recomputation
  short computedRefVariance;
  /// bookkeeping to track computation of variance increment to avoid
  /// unnecessary recomputation
  short computedDeltaVariance;

  /// track previous evaluation point for all_variables reference mean
  /// to avoid unnecessary recomputation
  RealVector xPrevRefMean;
  /// track previous evaluation point for all_variables mean increment
  /// to avoid unnecessary recomputation
  RealVector xPrevDeltaMean;
  /// track previous evaluation point for all_variables reference
  /// variance to avoid unnecessary recomputation
  RealVector xPrevRefVar;
  /// track previous evaluation point for all_variables variance
  /// increment to avoid unnecessary recomputation
  RealVector xPrevDeltaVar;

  /// storage for reference mean and variance
  RealVector referenceMoments;
  /// storage for mean and variance increments
  RealVector deltaMoments;
  /// storage for reference mean gradient
  RealVector meanRefGradient;
  /// storage for reference variance gradient
  RealVector varianceRefGradient;

  /// the type1 coefficients of the expansion for interpolating values
  std::map<UShortArray, RealVector2DArray> expansionType1Coeffs;
  /// iterator pointing to active node in expansionCoeffs
  std::map<UShortArray, RealVector2DArray>::iterator expT1CoeffsIter;

  /// the type2 coefficients of the expansion for interpolating gradients
  std::map<UShortArray, RealMatrix2DArray> expansionType2Coeffs;
  /// iterator pointing to active node in expansionCoeffGrads
  std::map<UShortArray, RealMatrix2DArray>::iterator expT2CoeffsIter;

  /// the gradients of the type1 expansion coefficients
  /** may be interpreted as either the gradients of the expansion coefficients
      or the coefficients of expansions for the response gradients.  This
      array is used when sensitivities of moments are needed with respect to
      variables that do not appear in the expansion (e.g., with respect to
      design variables for an expansion only over the random variables). */
  std::map<UShortArray, RealMatrix2DArray> expansionType1CoeffGrads;
  /// iterator pointing to active node in expansionCoeffGrads
  std::map<UShortArray, RealMatrix2DArray>::iterator expT1CoeffGradsIter;

  /// type 1 expansion coefficients popped during decrement for later
  /// restoration to expansionType1Coeffs
  std::map<UShortArray, std::map<UShortArray, RealVector> > poppedExpT1Coeffs;
  /// type 2 expansion coefficients popped during decrement for later
  /// restoration to expansionType2Coeffs
  std::map<UShortArray, std::map<UShortArray, RealMatrix> > poppedExpT2Coeffs;
  /// type 1 expansion coefficient gradients popped during decrement
  /// for later restoration to expansionType1CoeffGrads
  std::map<UShortArray, std::map<UShortArray, RealMatrix> >
    poppedExpT1CoeffGrads;

  /// roll up of expansion type 1 coefficients across all keys
  RealVector2DArray combinedExpT1Coeffs;
  /// roll up of expansion type 2 coefficient gradients across all keys
  RealMatrix2DArray combinedExpT2Coeffs;
  /// roll up of expansion type 1 coefficient gradients across all keys
  RealMatrix2DArray combinedExpT1CoeffGrads;
};


inline HierarchInterpPolyApproximation::
HierarchInterpPolyApproximation(const SharedBasisApproxData& shared_data):
  InterpPolyApproximation(shared_data)
{ }


inline HierarchInterpPolyApproximation::~HierarchInterpPolyApproximation()
{ }


inline void HierarchInterpPolyApproximation::update_active_iterators()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;
  origSurrData.active_key(key);
  if (deep_copied_surrogate_data())
    surrData.active_key(key);

  expT1CoeffsIter = expansionType1Coeffs.find(key);
  if (expT1CoeffsIter == expansionType1Coeffs.end()) {
    std::pair<UShortArray, RealVector2DArray> rv_pair(key, RealVector2DArray());
    expT1CoeffsIter = expansionType1Coeffs.insert(rv_pair).first;
  }
  expT2CoeffsIter = expansionType2Coeffs.find(key);
  if (expT2CoeffsIter == expansionType2Coeffs.end()) {
    std::pair<UShortArray, RealMatrix2DArray> rm_pair(key, RealMatrix2DArray());
    expT2CoeffsIter = expansionType2Coeffs.insert(rm_pair).first;
  }
  expT1CoeffGradsIter = expansionType1CoeffGrads.find(key);
  if (expT1CoeffGradsIter == expansionType1CoeffGrads.end()) {
    std::pair<UShortArray, RealMatrix2DArray> rm_pair(key, RealMatrix2DArray());
    expT1CoeffGradsIter = expansionType1CoeffGrads.insert(rm_pair).first;
  }
}


inline Real HierarchInterpPolyApproximation::variance()
{ return covariance(this); }


inline Real HierarchInterpPolyApproximation::variance(const RealVector& x)
{ return covariance(x, this); }


inline Real HierarchInterpPolyApproximation::delta_mean()
{
  if ( !(computedDeltaMean & 1) ) {
    SharedHierarchInterpPolyApproxData* data_rep
      = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
    UShort2DArray ref_key, incr_key;
    data_rep->hsg_driver()->partition_keys(ref_key, incr_key);
    return delta_mean(incr_key);
  }
  else
    return deltaMoments[0];
}


inline Real HierarchInterpPolyApproximation::delta_mean(const RealVector& x)
{
  if ( !(computedDeltaMean & 1) ) {
    SharedHierarchInterpPolyApproxData* data_rep
      = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
    UShort2DArray ref_key, incr_key;
    data_rep->hsg_driver()->partition_keys(ref_key, incr_key);
    return delta_mean(x, incr_key);
  }
  else
    return deltaMoments[0];
}


inline Real HierarchInterpPolyApproximation::delta_std_deviation()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_std_deviation(ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
delta_std_deviation(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_std_deviation(x, ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
delta_beta(bool cdf_flag, Real z_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_beta(cdf_flag, z_bar, ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
delta_beta(const RealVector& x, bool cdf_flag, Real z_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_beta(x, cdf_flag, z_bar, ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
delta_z(bool cdf_flag, Real beta_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_z(cdf_flag, beta_bar, ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
delta_z(const RealVector& x, bool cdf_flag, Real beta_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_z(x, cdf_flag, beta_bar, ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
expectation(const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs,
	    const UShort2DArray& set_partition)
{
  // This version defaults to active type1/2 wts from HierarchSparseGridDriver
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  return expectation(t1_coeffs, hsg_driver->type1_hierarchical_weight_sets(),
		     t2_coeffs, hsg_driver->type2_hierarchical_weight_sets(),
		     set_partition);
}


inline Real HierarchInterpPolyApproximation::
expectation(const RealVector& x, const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs,
	    const UShort2DArray& set_partition)
{
  // This version defaults to active sm_mi/key from HierarchSparseGridDriver
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  return expectation(x, t1_coeffs, t2_coeffs, hsg_driver->smolyak_multi_index(),
		     hsg_driver->collocation_key(), set_partition);
}


inline const RealVector& HierarchInterpPolyApproximation::
expectation_gradient(const RealMatrix2DArray& t1_coeff_grads)
{
  // This version defaults to active type1 wts from HierarchSparseGridDriver
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  return expectation_gradient(t1_coeff_grads,
    data_rep->hsg_driver()->type1_hierarchical_weight_sets());
}


inline Real HierarchInterpPolyApproximation::
expectation_gradient(const RealVector& x,
		     const RealMatrix2DArray& t1_coeff_grads, size_t t1cg_index)
{
  // This version defaults to active sm_mi/key from HierarchSparseGridDriver
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  return expectation_gradient(x, t1_coeff_grads,
			      hsg_driver->smolyak_multi_index(),
			      hsg_driver->collocation_key(), t1cg_index);
}


inline Real HierarchInterpPolyApproximation::
expectation_gradient(const RealVector& x, const RealVector2DArray& t1_coeffs,
		     const RealMatrix2DArray& t2_coeffs, size_t deriv_index)
{
  // This version defaults to active sm_mi/key from HierarchSparseGridDriver
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  return expectation_gradient(x, t1_coeffs, t2_coeffs,
			      hsg_driver->smolyak_multi_index(),
			      hsg_driver->collocation_key(), deriv_index);
}


inline void HierarchInterpPolyApproximation::push_coefficients()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;

  // mirror changes to origSurrData for deep copied surrData
  if (deep_copied_surrogate_data())
    surrData.push(data_rep->retrieval_index());

  push_coefficients(data_rep->hsg_driver()->trial_set());
  increment_current_from_reference();
}


inline Real HierarchInterpPolyApproximation::value(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return value(x, sm_mi, hsg_driver->collocation_key(), expT1CoeffsIter->second,
	       expT2CoeffsIter->second, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_basis_variables(x, sm_mi, hsg_driver->collocation_key(),
				  expT1CoeffsIter->second,
				  expT2CoeffsIter->second, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_basis_variables(x, sm_mi, hsg_driver->collocation_key(),
				  expT1CoeffsIter->second,
				  expT2CoeffsIter->second, dvv, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_nonbasis_variables(x, sm_mi, hsg_driver->collocation_key(),
				     expT1CoeffGradsIter->second, max_level);
}


inline const RealSymMatrix& HierarchInterpPolyApproximation::
hessian_basis_variables(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return hessian_basis_variables(x, sm_mi, hsg_driver->collocation_key(),
				 expT1CoeffsIter->second, max_level);
}


inline Real HierarchInterpPolyApproximation::
stored_value(const RealVector& x, const UShortArray& key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index(key);
  unsigned short   max_level = sm_mi.size() - 1;
  return value(x, sm_mi, hsg_driver->collocation_key(key),
	       expansionType1Coeffs[key], expansionType2Coeffs[key], max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x, const UShortArray& key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index(key);
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_basis_variables(x, sm_mi, hsg_driver->collocation_key(key),
				  expansionType1Coeffs[key],
				  expansionType2Coeffs[key], max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x, const SizetArray& dvv,
				const UShortArray& key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index(key);
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_basis_variables(x, sm_mi, hsg_driver->collocation_key(key),
				  expansionType1Coeffs[key],
				  expansionType2Coeffs[key], dvv, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x, const UShortArray& key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index(key);
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_nonbasis_variables(x, sm_mi, hsg_driver->collocation_key(key),
				     expansionType1CoeffGrads[key], max_level);
}


inline const RealSymMatrix& HierarchInterpPolyApproximation::
stored_hessian_basis_variables(const RealVector& x, const UShortArray& key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index(key);
  unsigned short   max_level = sm_mi.size() - 1;
  return hessian_basis_variables(x, sm_mi, hsg_driver->collocation_key(key),
				 expansionType1Coeffs[key], max_level);
}

} // namespace Pecos

#endif
