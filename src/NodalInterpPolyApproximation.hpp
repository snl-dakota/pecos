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
#include "SharedNodalInterpPolyApproxData.hpp"

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
  NodalInterpPolyApproximation(const SharedBasisApproxData& shared_data);
  /// destructor
  ~NodalInterpPolyApproximation();

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

  void combine_coefficients();

  bool update_active_iterators(const UShortArray& key);
  void combined_to_active(bool clear_combined = true);
  void clear_inactive();

  void integrate_response_moments(size_t num_moments, bool combined_stats);
  void integrate_expansion_moments(size_t num_moments, bool combined_stats);

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
  const RealVector& mean_gradient(const RealVector& x, const SizetArray& dvv);
  const RealVector& variance_gradient();
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

  Real covariance(PolynomialApproximation* poly_approx_2);
  Real covariance(const RealVector& x, PolynomialApproximation* poly_approx_2);

  Real combined_mean();
  Real combined_mean(const RealVector& x);

  Real combined_covariance(PolynomialApproximation* poly_approx_2);
  Real combined_covariance(const RealVector& x,
			   PolynomialApproximation* poly_approx_2);

  void compute_total_sobol_indices();
  void compute_partial_variance(const BitArray& set_value);

  RealVector approximation_coefficients(bool normalized) const;
  void approximation_coefficients(const RealVector& approx_coeffs,
				  bool normalized);

private:

  //
  //- Heading: Convenience functions
  //

  /// update expansionType{1Coeffs,2Coeffs,1CoeffGrads} following
  /// changes to surrogate data
  void update_expansion_coefficients();

  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using weights from the CombinedSparseGridDriver
  Real mean(const RealVector& t1_coeffs, const RealMatrix& t2_coeffs);
  /// helper function to evaluate mean(x) for passed coefficients
  Real mean(const RealVector& x, const RealVector& exp_t1_coeffs, 
	    const RealMatrix& exp_t2_coeffs);
  /// helper function to evaluate mean_gradient() for passed coefficient grads
  const RealVector& mean_gradient(const RealMatrix& exp_t1_coeff_grads);
  /// helper function to evaluate mean_gradient(x) for passed coefficients
  const RealVector& mean_gradient(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const RealMatrix& exp_t1_coeff_grads, const SizetArray& dvv);
  /// helper function to evaluate covariance() for passed coefficients
  Real covariance(Real mean_1, Real mean_2,    const RealVector& exp_t1c_1,
		  const RealMatrix& exp_t2c_1, const RealVector& exp_t1c_2,
		  const RealMatrix& exp_t2c_2);
  /// helper function to evaluate covariance(x) for passed coefficients
  Real covariance(const RealVector& x, Real mean_1, Real mean_2,
		  const RealVector& exp_t1c_1, const RealMatrix& exp_t2c_1,
		  const RealVector& exp_t1c_2, const RealMatrix& exp_t2c_2);
  /// helper function to evaluate variance_gradient() for passed coefficients
  const RealVector& variance_gradient(Real mean,
				      const RealVector& exp_t1_coeffs,
				      const RealMatrix& exp_t1_coeff_grads);
  /// helper function to evaluate variance_gradient(x) for passed coefficients
  const RealVector& variance_gradient(const RealVector& x, Real mean,
    const RealVector& mean_grad,      const RealVector& exp_t1_coeffs,
    const RealMatrix& exp_t2_coeffs,  const RealMatrix& exp_t1_coeff_grads,
    const SizetArray& dvv);
  
  /// compute the mean of a tensor interpolant on a tensor grid;
  /// contributes to mean(x)
  Real tensor_product_mean(const RealVector& x, const RealVector& exp_t1_coeffs,
    const RealMatrix& exp_t2_coeffs, const UShortArray& lev_index,
    const UShort2DArray& key, const SizetArray& colloc_index);

  /// compute the gradient of the mean of a tensor interpolant on a
  /// tensor grid; contributes to mean_gradient(x)
  const RealVector& tensor_product_mean_gradient(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const RealMatrix& exp_t1_coeff_grads, const UShortArray& lev_index,
    const UShort2DArray& key, const SizetArray& colloc_index,
    const SizetArray& dvv);

  /// compute the covariance of two tensor interpolants on the same
  /// tensor grid using an interpolation of products or product of
  /// interpolants approach; contributes to covariance(x, poly_approx_2)
  Real tensor_product_covariance(const RealVector& x, Real mean_1, Real mean_2,
    const RealVector& exp_t1c_1, const RealMatrix& exp_t2c_1,
    const RealVector& exp_t1c_2, const RealMatrix& exp_t2c_2,
    const UShortArray& lev_index, const UShort2DArray& key,
    const SizetArray& colloc_index);

  /// compute the covariance of two tensor interpolants on the same grid;
  /// contributes to covariance(x, poly_approx_2)
  Real product_of_interpolants(const RealVector& x, Real mean_1, Real mean_2,
			       const RealVector& exp_t1c_1,
			       const RealMatrix& exp_t2c_1,
			       const RealVector& exp_t1c_2,
			       const RealMatrix& exp_t2c_2,
			       const UShortArray& lev_index,
			       const UShort2DArray& colloc_key,
			       const SizetArray& colloc_index);
  /// compute the covariance of two tensor interpolants on different tensor
  /// grids using a PRODUCT_OF_INTERPOLANTS_FULL approach; contributes to
  /// covariance(x, poly_approx_2)
  Real product_of_interpolants(const RealVector& x, Real mean_1, Real mean_2,
			       const RealVector& exp_t1c_1,
			       const RealMatrix& exp_t2c_1,
			       const RealVector& exp_t1c_2,
			       const RealMatrix& exp_t2c_2,
			       const UShortArray& lev_index_1,
			       const UShort2DArray& colloc_key_1,
			       const SizetArray& colloc_index_1,
			       const UShortArray& lev_index_2,
			       const UShort2DArray& colloc_key_2,
			       const SizetArray& colloc_index_2);

  /// compute the gradient of the variance of a tensor interpolant on
  /// a tensor grid using an interpolation of products or product of
  /// interpolants approach; contributes to variance_gradient(x)
  const RealVector& tensor_product_variance_gradient(const RealVector& x,
    Real mean, const RealVector& mean_grad, const RealVector& exp_t1_coeffs,
    const RealMatrix& exp_t2_coeffs, const RealMatrix& exp_t1_coeff_grads,
    const UShortArray& lev_index,    const UShort2DArray& key,
    const SizetArray& colloc_index,  const SizetArray& dvv);

  /// compute the value of an expansion
  Real value(const RealVector& x, const RealVector& exp_t1_coeffs,
	     const RealMatrix& exp_t2_coeffs);
  /// compute the value of a tensor-product interpolant
  Real value(const RealVector& x, const RealVector& exp_t1_coeffs,
	     const RealMatrix& exp_t2_coeffs, const UShortArray& lev_index,
	     const UShort2DArray& colloc_key);
  /// compute the value of a sparse interpolant
  Real value(const RealVector& x, const RealVector& exp_t1_coeffs,
	     const RealMatrix& exp_t2_coeffs, const UShort2DArray& sm_mi,
	     const IntArray& sm_coeffs, const UShort3DArray& colloc_key,
	     const Sizet2DArray& colloc_index);
  /// compute value of reduced-dimension interpolant
  Real value(const RealVector& x, const RealVectorArray& t1_coeffs,
	     const RealMatrixArray& t2_coeffs, const UShort3DArray& colloc_key,
	     const SizetList& subset_indices);

  /// compute the gradient of an expansion with respect to basis variables
  const RealVector& gradient_basis_variables(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs);
  /// compute the gradient of an expansion with respect to basis variables
  /// for given DVV
  const RealVector& gradient_basis_variables(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const SizetArray& dvv);
  /// compute the gradient of a tensor-product interpolant with respect
  /// to basis variables
  const RealVector& gradient_basis_variables(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShortArray& lev_index, const UShort2DArray& colloc_key);
  /// compute the gradient of a tensor-product interpolant with respect
  /// to basis variables for given DVV
  const RealVector& gradient_basis_variables(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShortArray& lev_index, const UShort2DArray& colloc_key,
    const SizetArray& dvv);
  /// compute the gradient of a sparse interpolant with respect to
  /// basis variables
  const RealVector& gradient_basis_variables(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShort2DArray& sm_mi, const IntArray& sm_coeffs,
    const UShort3DArray& colloc_key, const Sizet2DArray& colloc_index);
  /// compute the gradient of a sparse interpolant with respect to
  /// basis variables for given DVV
  const RealVector& gradient_basis_variables(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShort2DArray& sm_mi, const IntArray& sm_coeffs,
    const UShort3DArray& colloc_key, const Sizet2DArray& colloc_index,
    const SizetArray& dvv);
  /// compute gradient of reduced-dimension interpolant with respect
  /// to basis variables
  const RealVector& gradient_basis_variables(const RealVector& x,
    const RealVectorArray& t1_coeffs, const RealMatrixArray& t2_coeffs,
    const UShort3DArray& colloc_key, const SizetList& subset_indices);

  /// compute the gradient of an expansion with respect to non-basis variables
  const RealVector& gradient_nonbasis_variables(const RealVector& x,
    const RealMatrix& exp_t1_coeff_grads);
  /// compute the gradient of a tensor-product interpolant with respect
  /// to non-basis variables
  const RealVector& gradient_nonbasis_variables(const RealVector& x,
    const RealMatrix& exp_t1_coeff_grads, const UShortArray& lev_index,
    const UShort2DArray& colloc_key);
  /// compute the gradient of a sparse interpolant with respect to
  /// non-basis variables
  const RealVector& gradient_nonbasis_variables(const RealVector& x,
    const RealMatrix& exp_t1_coeff_grads, const UShort2DArray& sm_mi,
    const IntArray& sm_coeffs, const UShort3DArray& colloc_key,
    const Sizet2DArray& colloc_index);

  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using t{1,2}_wts
  Real expectation(const RealVector& t1_coeffs, const RealVector& t1_wts,
		   const RealMatrix& t2_coeffs, const RealMatrix& t2_wts);

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

  /// the type1 coefficients of the expansion for interpolating values
  std::map<UShortArray, RealVector> expansionType1Coeffs;
  /// iterator pointing to active node in expansionCoeffs
  std::map<UShortArray, RealVector>::iterator expT1CoeffsIter;

  /// the type2 coefficients of the expansion for interpolating gradients
  std::map<UShortArray, RealMatrix> expansionType2Coeffs;
  /// iterator pointing to active node in expansionCoeffGrads
  std::map<UShortArray, RealMatrix>::iterator expT2CoeffsIter;

  /// the gradients of the type1 expansion coefficients
  /** may be interpreted as either the gradients of the expansion
      coefficients or the coefficients of expansions for the response
      gradients.  This array is used when sensitivities of moments are
      needed with respect to variables that do not appear in the
      expansion (e.g., with respect to design variables for an
      expansion only over the random variables). */
  std::map<UShortArray, RealMatrix> expansionType1CoeffGrads;
  /// iterator pointing to active node in expansionCoeffGrads
  std::map<UShortArray, RealMatrix>::iterator expT1CoeffGradsIter;

  /// roll up of expansion type 1 coefficients across all keys
  RealVector combinedExpT1Coeffs;
  /// roll up of expansion type 2 coefficient gradients across all keys
  RealMatrix combinedExpT2Coeffs;
  /// roll up of expansion type 1 coefficient gradients across all keys
  RealMatrix combinedExpT1CoeffGrads;

  /// the gradient of the mean of a tensor-product interpolant; a
  /// contributor to meanGradient
  RealVector tpMeanGrad;     // TO DO: move to shared data? (2nd pass tuning)
  /// the gradient of the variance of a tensor-product interpolant; a
  /// contributor to varianceGradient
  RealVector tpVarianceGrad; // TO DO: move to shared data? (2nd pass tuning)
};


inline NodalInterpPolyApproximation::
NodalInterpPolyApproximation(const SharedBasisApproxData& shared_data):
  InterpPolyApproximation(shared_data),
  expT1CoeffsIter(expansionType1Coeffs.end())
{ }


inline NodalInterpPolyApproximation::~NodalInterpPolyApproximation()
{ }


inline bool NodalInterpPolyApproximation::
update_active_iterators(const UShortArray& key)
{
  // Test for change
  if (expT1CoeffsIter != expansionType1Coeffs.end() &&
      expT1CoeffsIter->first == key)
    return false;
  
  expT1CoeffsIter = expansionType1Coeffs.find(key);
  if (expT1CoeffsIter == expansionType1Coeffs.end()) {
    std::pair<UShortArray, RealVector> rv_pair(key, RealVector());
    expT1CoeffsIter = expansionType1Coeffs.insert(rv_pair).first;
  }
  expT2CoeffsIter = expansionType2Coeffs.find(key);
  if (expT2CoeffsIter == expansionType2Coeffs.end()) {
    std::pair<UShortArray, RealMatrix> rm_pair(key, RealMatrix());
    expT2CoeffsIter = expansionType2Coeffs.insert(rm_pair).first;
  }
  expT1CoeffGradsIter = expansionType1CoeffGrads.find(key);
  if (expT1CoeffGradsIter == expansionType1CoeffGrads.end()) {
    std::pair<UShortArray, RealMatrix> rm_pair(key, RealMatrix());
    expT1CoeffGradsIter = expansionType1CoeffGrads.insert(rm_pair).first;
  }

  InterpPolyApproximation::update_active_iterators(key);
  return true;
}


inline void NodalInterpPolyApproximation::increment_coefficients()
{
  // TO DO: partial sync for new TP data set, e.g. update_surrogate_data() ?
  synchronize_surrogate_data();

  update_expansion_coefficients(); // updates iterators, clears computed bits

  allocate_component_sobol();
}


inline void NodalInterpPolyApproximation::push_coefficients()
{ update_expansion_coefficients(); } // updates iterators, clears computed bits


inline void NodalInterpPolyApproximation::finalize_coefficients()
{ update_expansion_coefficients(); } // updates iterators, clears computed bits


inline RealVector NodalInterpPolyApproximation::
approximation_coefficients(bool normalized) const
{
  if (normalized)
    PCerr << "Warning: normalized coefficients not supported in "
	  << "InterpPolyApproximation export." << std::endl;
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  if (data_rep->basisConfigOptions.useDerivs) {
    PCerr << "Error: approximation_coefficients() not supported in "
	  << "InterpPolyApproximation for type2 coefficients." << std::endl;
    return abort_handler_t<const RealVector&>(-1);
  }
  else {
    RealVector& exp_t1c = expT1CoeffsIter->second;
    return RealVector(Teuchos::View, exp_t1c.values(), exp_t1c.length());
  }
}


inline void NodalInterpPolyApproximation::
approximation_coefficients(const RealVector& approx_coeffs, bool normalized)
{
  if (normalized)
    PCerr << "Warning: normalized coefficients not supported in "
	  << "InterpPolyApproximation import." << std::endl;
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  if (data_rep->basisConfigOptions.useDerivs) {
    PCerr << "Error: approximation_coefficients() not supported in "
	  << "InterpPolyApproximation for type2 coefficients." << std::endl;
    abort_handler(-1);
  }
  else {
    expT1CoeffsIter->second = approx_coeffs;

    allocate_total_sobol();
    allocate_component_sobol();

    if (numericalMoments.empty()) {
      SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
      size_t num_moments = (data_rep->nonRandomIndices.empty()) ? 4 : 2;
      numericalMoments.sizeUninitialized(num_moments);
    }
  }
}


inline Real NodalInterpPolyApproximation::mean()
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  //IntegrationDriver* driver_rep = data_rep->driverRep;
  //if (!driver_rep->track_unique_product_weights()) {
  //  PCerr << "Error: unique product weights required in "
  //	  << "NodalInterpPolyApproximation::mean()" << std::endl;
  //  abort_handler(-1);
  //}
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedMean & 1))
    return numericalMoments[0];

  Real mu = mean(expT1CoeffsIter->second, expT2CoeffsIter->second);
  if (std_mode)
    { numericalMoments[0] = mu; computedMean |= 1; }
  return mu;
}


inline Real NodalInterpPolyApproximation::mean(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedMean & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevMean))
    return numericalMoments[0];

  Real mu = mean(x, expT1CoeffsIter->second, expT2CoeffsIter->second);
  if (all_mode)
    { numericalMoments[0] = mu; computedMean |= 1; xPrevMean = x; }
  return mu;
}


inline const RealVector& NodalInterpPolyApproximation::mean_gradient()
{
  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	  << "InterpPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedMean & 2))
    return meanGradient;

  if (std_mode) computedMean |=  2; //   activate 2-bit
  else          computedMean &= ~2; // deactivate 2-bit: protect mixed usage
  return mean_gradient(expT1CoeffGradsIter->second);
}


inline const RealVector& NodalInterpPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  // if already computed, return previous result
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if ( all_mode && (computedMean & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevMeanGrad) ) // && dvv == dvvPrev)
    switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE:
      return tpMeanGrad;   break;
    case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID:
      return meanGradient; break;
    }

  // compute the gradient of the mean
  if (all_mode) { computedMean |=  2; xPrevMeanGrad = x; }
  else            computedMean &= ~2; // deactivate 2-bit: protect mixed usage
  return mean_gradient(x, expT1CoeffsIter->second, expT2CoeffsIter->second,
		       expT1CoeffGradsIter->second, dvv);
}


inline const RealVector& NodalInterpPolyApproximation::variance_gradient()
{
  // Error check for required data
  if (!expansionCoeffFlag || !expansionCoeffGradFlag) {
    PCerr << "Error: insufficient expansion coefficient data in NodalInterp"
	  << "PolyApproximation::variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedVariance & 2))
    return varianceGradient;

  if (std_mode) computedVariance |=  2;
  else          computedVariance &= ~2; // deactivate 2-bit: protect mixed usage
  return variance_gradient(mean(), expT1CoeffsIter->second,
			   expT1CoeffGradsIter->second);
}


inline const RealVector& NodalInterpPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  // if already computed, return previous result
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if ( all_mode && (computedVariance & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevVarGrad) ) // && dvv == dvvPrev)
    switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE:
      return tpVarianceGrad;   break;
    case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID:
      return varianceGradient; break;
    }

  if (all_mode) { computedVariance |=  2; xPrevVarGrad = x; }
  else            computedVariance &= ~2;//deactivate 2-bit: protect mixed usage
  // don't compute expansion mean/mean_grad for case where tensor
  // means/mean_grads will be used
  return (data_rep->momentInterpType == PRODUCT_OF_INTERPOLANTS_FAST) ?
    variance_gradient(x, 0., meanGradient, // dummy values (not used)
		      expT1CoeffsIter->second, expT2CoeffsIter->second,
		      expT1CoeffGradsIter->second, dvv) :
    variance_gradient(x, mean(x), mean_gradient(x, dvv),
		      expT1CoeffsIter->second, expT2CoeffsIter->second,
		      expT1CoeffGradsIter->second,dvv);
}


inline Real NodalInterpPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  NodalInterpPolyApproximation* nip_approx_2
    = (NodalInterpPolyApproximation*)poly_approx_2;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  bool same  = (this == nip_approx_2),
    std_mode = data_rep->nonRandomIndices.empty();

  // Error check for required data
  if ( !expansionCoeffFlag || ( !same && !nip_approx_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  // TO DO:
  //IntegrationDriver* driver_rep = data_rep->driverRep;
  //if (!driver_rep->track_unique_product_weights()) {
  //  PCerr << "Error: unique product weights required in "
  //	    << "NodalInterpPolyApproximation::covariance()" << std::endl;
  //  abort_handler(-1);
  //}

  if (same && std_mode && (computedVariance & 1))
    return numericalMoments[1];

  Real mean_1 = mean(), mean_2 = (same) ? mean_1 : nip_approx_2->mean(),
    covar = covariance(mean_1, mean_2, expT1CoeffsIter->second,
		       expT2CoeffsIter->second,
		       nip_approx_2->expT1CoeffsIter->second,
		       nip_approx_2->expT2CoeffsIter->second);
  if (same && std_mode)
    { numericalMoments[1] = covar; computedVariance |= 1; }
  return covar;
}


inline Real NodalInterpPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  NodalInterpPolyApproximation* nip_approx_2
    = (NodalInterpPolyApproximation*)poly_approx_2;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  bool same = (this == nip_approx_2),
    all_mode = !data_rep->nonRandomIndices.empty();

  // Error check for required data
  if ( !expansionCoeffFlag || ( !same && !nip_approx_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if ( same && all_mode && (computedVariance & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar) )
    return numericalMoments[1];

  Real mean_1, mean_2;
  if (data_rep->momentInterpType == PRODUCT_OF_INTERPOLANTS_FAST)
    // don't compute means here since not used by SSG with _FAST approximation
    // (instead, tensor means are computed for each tensor grid)
    mean_1 = mean_2 = 0.;
  else {
    // TQP and SSG with PRODUCT_OF_INTERPOLANTS_FULL: compute global mean once
    mean_1 = mean(x);
    mean_2 = (same) ? mean_1 : nip_approx_2->mean(x);
  }
  Real covar = covariance(x, mean_1, mean_2, expT1CoeffsIter->second,
			  expT2CoeffsIter->second,
			  nip_approx_2->expT1CoeffsIter->second,
			  nip_approx_2->expT2CoeffsIter->second);
  if (same && all_mode)
    { numericalMoments[1] = covar; computedVariance |= 1; xPrevVar = x; }
  return covar;
}


inline Real NodalInterpPolyApproximation::combined_mean()
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedMean & 1))
    return numericalMoments[0];

  Real mu = mean(combinedExpT1Coeffs, combinedExpT2Coeffs);
  if (std_mode)
    { numericalMoments[0] = mu; computedMean |= 1; }
  return mu;
}


inline Real NodalInterpPolyApproximation::combined_mean(const RealVector& x)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedMean & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevMean))
    return numericalMoments[0];

  Real mu = mean(x, combinedExpT1Coeffs, combinedExpT2Coeffs);
  if (all_mode)
    { numericalMoments[0] = mu; computedMean |= 1; xPrevMean = x; }
  return mu;
}


inline Real NodalInterpPolyApproximation::
combined_covariance(PolynomialApproximation* poly_approx_2)
{
  NodalInterpPolyApproximation* nip_approx_2
    = (NodalInterpPolyApproximation*)poly_approx_2;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  bool same  = (this == nip_approx_2),
    std_mode = data_rep->nonRandomIndices.empty();

  if (same && std_mode && (computedVariance & 1))
    return numericalMoments[1];

  Real mean_1 = combined_mean(),
       mean_2 = (this == nip_approx_2) ? mean_1 : nip_approx_2->combined_mean();
  Real covar
    = covariance(mean_1, mean_2, combinedExpT1Coeffs, combinedExpT2Coeffs,
		 nip_approx_2->combinedExpT1Coeffs,
		 nip_approx_2->combinedExpT2Coeffs);

  if (same && std_mode)
    { numericalMoments[1] = covar; computedVariance |= 1; }
  return covar;
}


inline Real NodalInterpPolyApproximation::
combined_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  NodalInterpPolyApproximation* nip_approx_2
    = (NodalInterpPolyApproximation*)poly_approx_2;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  bool same = (this == nip_approx_2),
    all_mode = !data_rep->nonRandomIndices.empty();

  if ( same && all_mode && (computedVariance & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar) )
    return numericalMoments[1];

  Real mean_1, mean_2;
  if (data_rep->momentInterpType == PRODUCT_OF_INTERPOLANTS_FAST)
    mean_1 = mean_2 = 0.; // don't compute exp mean since tensor means used
  else {
    mean_1 = combined_mean(x);
    mean_2 = (this == nip_approx_2) ? mean_1 : nip_approx_2->combined_mean(x);
  }
  Real covar
    = covariance(x, mean_1, mean_2, combinedExpT1Coeffs, combinedExpT2Coeffs,
		 nip_approx_2->combinedExpT1Coeffs,
		 nip_approx_2->combinedExpT2Coeffs);

  if (same && all_mode)
    { numericalMoments[1] = covar; computedVariance |= 1; xPrevVar = x; }
  return covar;
}


inline Real NodalInterpPolyApproximation::value(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  return value(x, expT1CoeffsIter->second, expT2CoeffsIter->second);
}


inline const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  return gradient_basis_variables(x, expT1CoeffsIter->second,
				  expT2CoeffsIter->second);
}


inline const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  return gradient_basis_variables(x, expT1CoeffsIter->second,
				  expT2CoeffsIter->second, dvv);
}


inline const RealVector& NodalInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::gradient_nonbasis_variables()" << std::endl;
    abort_handler(-1);
  }

  return gradient_nonbasis_variables(x, expT1CoeffGradsIter->second);
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

} // namespace Pecos

#endif
