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
  HierarchInterpPolyApproximation(short basis_type, size_t num_vars,
				  bool use_derivs);
  /// destructor
  ~HierarchInterpPolyApproximation();

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
  /// storedColloc{Key,Indices}, and storedLevMultiIndex
  void store_coefficients();
  /// augment current interpolant using
  /// storedExpType{1Coeffs,2Coeffs,1CoeffGrads}, storedColloc{Key,Indices},
  /// and storedLevMultiIndex
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

  Real delta_mean();
  Real delta_std_deviation();
  Real delta_beta(bool cdf_flag, Real z_bar);
  Real delta_z(bool cdf_flag, Real beta_bar);
  Real delta_covariance(PolynomialApproximation* poly_approx_2);
  Real delta_covariance(const RealVector& x,
			PolynomialApproximation* poly_approx_2);

  void compute_total_sobol_indices();
  void compute_partial_variance(const BitArray& set_value);

private:

  //
  //- Heading: Convenience functions
  //

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

  /// update bookkeeping when adding a grid increment relative to the
  /// grid reference
  void increment_current_from_reference();
  /// update bookkeeping when removing a grid increment and returning
  /// to the grid reference
  void decrement_current_to_reference();

  /// compute the reference mean, excluding the current grid increment
  Real reference_mean();
  /// compute the reference mean, excluding the current grid
  /// increment, using ref_key
  Real reference_mean(const UShort2DArray& ref_key);
  /// compute the reference variance, excluding the current grid increment
  Real reference_variance();
  /// compute the reference variance, excluding the current grid
  /// increment, using ref_key
  Real reference_variance(const UShort2DArray& ref_key);

  /// compute the mean increment due to the current grid increment
  Real delta_mean(const UShort2DArray& incr_key);
  /// compute the variance increment due to the current grid increment
  Real delta_variance();
  /// compute the variance increment due to the current grid increment
  Real delta_variance(const UShort2DArray& ref_key,
		      const UShort2DArray& incr_key);
  /// compute the standard deviation increment due to the current grid increment
  Real delta_std_deviation(const UShort2DArray& ref_key,
			   const UShort2DArray& incr_key);
  /// compute the reliability index increment due to the current grid increment
  Real delta_beta(bool cdf_flag, Real z_bar, const UShort2DArray& ref_key,
		  const UShort2DArray& incr_key);
  /// compute the response level increment due to the current grid increment
  Real delta_z(bool cdf_flag, Real beta_bar, const UShort2DArray& ref_key,
	       const UShort2DArray& incr_key);

  /// compute the covariance increment due to the current grid increment
  Real delta_covariance(PolynomialApproximation* poly_approx_2,
			const UShort2DArray& ref_key,
			const UShort2DArray& incr_key);

  /// form type 1/2 coefficients for interpolation of R_1 R_2
  void product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
    RealVector2DArray& r1r2_t1_coeffs, RealMatrix2DArray& r1r2_t2_coeffs,
    const UShort2DArray& reference_key = UShort2DArray());
  /// form type 1/2 coefficients for interpolation of (R_1 - mu_1)(R_2 - mu_2)
  void central_product_interpolant(
    HierarchInterpPolyApproximation* hip_approx_2, Real mean_1, Real mean_2,
    RealVector2DArray& cov_t1_coeffs, RealMatrix2DArray& cov_t2_coeffs,
    const UShort2DArray& reference_key = UShort2DArray());

  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using weights from the HierarchSparseGridDriver
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

  /// increment expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  /// for a single index_set
  void increment_expansion_coefficients(const UShortArray& index_set);

  /// move the expansion coefficients for restore_set from
  /// savedExp{T1Coeffs,T2Coeffs,T1CoeffGrads} to
  /// expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void restore_expansion_coefficients(const UShortArray& restore_set);

  /// compute member expansion for Sobol' index integration
  void member_coefficients_weights(const BitArray&    member_bits,
				   UShort4DArray&     member_colloc_key,
				   Sizet3DArray&      member_colloc_index,
				   RealVector2DArray& member_t1_coeffs,
				   RealVector2DArray& member_t1_wts,
				   RealMatrix2DArray& member_t2_coeffs,
				   RealMatrix2DArray& member_t2_wts);
  /// form hierarchical interpolant of h^2 from member-variable
  /// expansion of h
  void product_member_coefficients(const BitArray& m_bits,
				   const UShort4DArray& m_colloc_key,
				   const Sizet3DArray&  m_colloc_index,
				   const RealVector2DArray& m_t1_coeffs,
				   const RealMatrix2DArray& m_t2_coeffs,
				   RealVector2DArray& prod_m_t1_coeffs,
				   RealMatrix2DArray& prod_m_t2_coeffs);
  /// form hierarchical interpolant of (h-\mu)^2 from member-variable
  /// expansion of h
  void central_product_member_coefficients(const BitArray& m_bits,
    const UShort4DArray& m_colloc_key, const Sizet3DArray&  m_colloc_index,
    const RealVector2DArray& m_t1_coeffs, const RealMatrix2DArray& m_t2_coeffs,
    Real mean, RealVector2DArray& cprod_m_t1_coeffs,
    RealMatrix2DArray& cprod_m_t2_coeffs);

  //
  //- Heading: Data
  //

  /// Pecos:PIECEWISE_INTERP_POLYNOMIAL or Pecos:PIECEWISE_CUBIC_INTERP
  short polyType;

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
  /// storage for reference mean and variance
  RealVector referenceMoments;
  /// storage for mean and variance increments
  RealVector deltaMoments;
  /// storage for reference mean gradient
  RealVector meanRefGradient;
  /// storage for reference variance gradient
  RealVector varianceRefGradient;

  /// the type1 coefficients of the expansion for interpolating values
  RealVector2DArray expansionType1Coeffs;
  /// the type2 coefficients of the expansion for interpolating gradients
  RealMatrix2DArray expansionType2Coeffs;
  /// the gradients of the type1 expansion coefficients
  /** may be interpreted as either the gradients of the expansion coefficients
      or the coefficients of expansions for the response gradients.  This
      array is used when sensitivities of moments are needed with respect to
      variables that do not appear in the expansion (e.g., with respect to
      design variables for an expansion only over the random variables). */
  RealMatrix2DArray expansionType1CoeffGrads;

  /// saved type 1 expansion coefficients for restoration to
  /// expansionType1Coeffs
  std::map<UShortArray, RealVector> savedExpT1Coeffs;
  /// saved type 2 expansion coefficients for restoration to
  /// expansionType2Coeffs
  std::map<UShortArray, RealMatrix> savedExpT2Coeffs;
  /// saved type 1 expansion coefficient gradients for restoration to
  /// expansionType1CoeffGrads
  std::map<UShortArray, RealMatrix> savedExpT1CoeffGrads;

  /// storage of expansionType1Coeffs state for subsequent restoration
  RealVector2DArray storedExpType1Coeffs;
  /// storage of expansionType2Coeffs state for subsequent restoration
  RealMatrix2DArray storedExpType2Coeffs;
  /// storage of expansionType1CoeffGrads state for subsequent restoration
  RealMatrix2DArray storedExpType1CoeffGrads;
  /// storage of level multi-index (levels for tensor or sparse grids)
  /// for subsequent restoration
  UShort3DArray storedLevMultiIndex;
  /// storage of IntegrationDriver collocation key state for
  /// subsequent restoration
  UShort4DArray storedCollocKey;
  // storage of IntegrationDriver collocation indices state for
  // subsequent restoration
  //Sizet3DArray storedCollocIndices;
};


inline HierarchInterpPolyApproximation::
HierarchInterpPolyApproximation(short basis_type, size_t num_vars,
				bool use_derivs):
  InterpPolyApproximation(basis_type, num_vars, use_derivs)
{}


inline HierarchInterpPolyApproximation::~HierarchInterpPolyApproximation()
{}


inline Real HierarchInterpPolyApproximation::
expectation(const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs,
	    const UShort2DArray& set_partition)
{
  // This version defaults to type1/2 weights from HierarchSparseGridDriver
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  return expectation(t1_coeffs, hsg_driver->type1_weight_set_arrays(),
		     t2_coeffs, hsg_driver->type2_weight_set_arrays(),
		     set_partition);
}


inline void HierarchInterpPolyApproximation::restore_expansion_coefficients()
{
  SparseGridDriver* sg_driver = (SparseGridDriver*)driverRep;
  restore_expansion_coefficients(sg_driver->trial_set());
  increment_current_from_reference();
}


inline Real HierarchInterpPolyApproximation::value(const RealVector& x)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return value(x, sm_mi, hsg_driver->collocation_key(), expansionType1Coeffs,
	       expansionType2Coeffs, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_basis_variables(x, sm_mi, hsg_driver->collocation_key(),
				  expansionType1Coeffs, expansionType2Coeffs,
				  max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_basis_variables(x, sm_mi, hsg_driver->collocation_key(),
				  expansionType1Coeffs, expansionType2Coeffs,
				  dvv, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_nonbasis_variables(x, sm_mi, hsg_driver->collocation_key(),
				     expansionType1CoeffGrads, max_level);
}


inline Real HierarchInterpPolyApproximation::stored_value(const RealVector& x)
{
  unsigned short max_level = storedLevMultiIndex.size() - 1;
  return value(x, storedLevMultiIndex, storedCollocKey, storedExpType1Coeffs,
	       storedExpType2Coeffs, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x)
{
  unsigned short max_level = storedLevMultiIndex.size() - 1;
  return gradient_basis_variables(x, storedLevMultiIndex, storedCollocKey,
				  storedExpType1Coeffs, storedExpType2Coeffs,
				  max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x)
{
  unsigned short max_level = storedLevMultiIndex.size() - 1;
  return gradient_nonbasis_variables(x, storedLevMultiIndex, storedCollocKey,
				     storedExpType1CoeffGrads, max_level);
}

} // namespace Pecos

#endif
