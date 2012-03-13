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
			       bool use_derivs);
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
  void compute_partial_variance(int set_value);

  const RealVector& approximation_coefficients() const;
  void approximation_coefficients(const RealVector& approx_coeffs);

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

  /// compute the covariance of a tensor interpolant on a tensor grid;
  /// contributes to covariance(x, poly_approx_2)
  Real tensor_product_covariance(const RealVector& x,
    const UShortArray& lev_index,   const UShort2DArray& key,
    const SizetArray& colloc_index, NodalInterpPolyApproximation* nip_approx_2);

  /// compute the gradient of the variance of a tensor interpolant on
  /// a tensor grid; contributes to variance_gradient(x)
  const RealVector& tensor_product_variance_gradient(const RealVector& x,
    const UShortArray& lev_index, const UShort2DArray& key,
    const SizetArray& colloc_index, const SizetArray& dvv);

  //
  //- Heading: Data
  //

  /// type of interpolation for variance and variance gradient
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
			     bool use_derivs):
  InterpPolyApproximation(basis_type, num_vars, use_derivs),
  //momentInterpType(INTERPOLATION_OF_PRODUCTS_GLOBAL_MEAN)
  //momentInterpType(INTERPOLATION_OF_PRODUCTS_TENSOR_MEAN)
  //momentInterpType(PRODUCT_OF_INTERPOLANTS_GLOBAL_MEAN)
  momentInterpType(PRODUCT_OF_INTERPOLANTS_TENSOR_MEAN)
{ }


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

} // namespace Pecos

#endif
