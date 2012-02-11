/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchInterpPolyApproximation
//- Description:  Class for Nodal Interpolation Polynomial Approximation
//-               
//- Owner:        Chris Miller

#ifndef HIERARCH_INTERP_POLY_APPROXIMATION_HPP
#define HIERARCH_INTERP_POLY_APPROXIMATION_HPP

#include "InterpPolyApproximation.hpp"
#include "HierarchSparseGridDriver.hpp"

namespace Pecos {


/// Derived approximation class for piecewise linear and cubic hierarchical 
/// interpolation polynomials (local approximation interpolating function 
/// values and potentially gradients at collocation points).

/** The HierarchInterpPolyApproximation class provides a local piecewise 
    polynomial approximation based on the hierarchical approach described in 
    X. Ma and N. Zabaras "An adaptive hierarchical sparse grid collocation 
    algorithm for the solution of stochastic differential equations", Journal 
    of Computational Physics, 228 (2009), 3084-3113.  Both piecewise linear 
    basis functions using function values at the collocation points and cubic 
    Hermite basis functions using both values and derivatives are available.  
    It is used primarily for stochastic collocation approaches to uncertainty 
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

  /// derived portion of allocate_arrays()
  void allocate_expansion_coefficients();
  /// derived portion of compute_coefficients()
  void compute_expansion_coefficients();
  /// store current state within storedExpType{1Coeffs,2Coeffs,1CoeffGrads},
  /// storedColloc{Key,Indices}, and storedLevMultiIndex
  void store_coefficients();
  /// augment current interpolant using
  /// storedExpType{1Coeffs,2Coeffs,1CoeffGrads}, storedColloc{Key,Indices},
  /// and storedLevMultiIndex
  void combine_coefficients(short combine_type);
  void restore_expansion_coefficients();

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

  void compute_total_sobol_indices();
  void compute_partial_variance(int set_value);
  void member_coefficients_weights(int set_value,
    const UShortArray& quad_order, const UShortArray& lev_index,
    const UShort2DArray& key, const SizetArray& colloc_index,
    RealVector& member_coeffs, RealVector& member_wts);

private:

  //
  //- Heading: Convenience functions
  //

  /// compute the value for a particular level and a particular index set
  Real tensor_product_value(const RealVector& x, const RealVector& t1_coeffs,
			    const RealMatrix& t2_coeffs,
			    const UShortArray& basis_index,
			    const UShort2DArray& key);

  /// compute the value at a point for a particular interpolation level
  Real value(const RealVector& x, unsigned short level);

  /// compute the approximate gradient with respect to the basis variables
  /// at a particular point for a particular interpolation level
  const RealVector& gradient_basis_variables(const RealVector& x,
					     unsigned short level);
  /// compute the approximate gradient with respect to the basis variables
  /// for a particular point, interpolation level, and DVV
  const RealVector& gradient_basis_variables(const RealVector& x,
					     const SizetArray& dvv,
					     unsigned short level);
  /// compute the approximate gradient with respect to the nonbasis
  /// variables at a particular point for a particular interpolation level
  const RealVector& gradient_nonbasis_variables(const RealVector& x,
						unsigned short level);

  //
  //- Heading: Data
  //

  /// Pecos:PIECEWISE_INTERP_POLYNOMIAL or Pecos:PIECEWISE_CUBIC_INTERP
  short polyType;

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
  /// storage of IntegrationDriver collocation indices state for
  /// subsequent restoration
  Sizet3DArray storedCollocIndices;
};


inline HierarchInterpPolyApproximation::
HierarchInterpPolyApproximation(short basis_type, size_t num_vars,
				bool use_derivs):
  InterpPolyApproximation(basis_type, num_vars, use_derivs)
{}


inline HierarchInterpPolyApproximation::~HierarchInterpPolyApproximation()
{}


inline Real HierarchInterpPolyApproximation::value(const RealVector& x)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  unsigned short max_level = hsg_driver->smolyak_multi_index().size()-1;
  return value(x, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  unsigned short max_level = hsg_driver->smolyak_multi_index().size()-1;
  return gradient_basis_variables(x, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  unsigned short max_level = hsg_driver->smolyak_multi_index().size()-1;
  return gradient_basis_variables(x, dvv, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  unsigned short max_level = hsg_driver->smolyak_multi_index().size()-1;
  return gradient_nonbasis_variables(x, max_level);
}

} // namespace Pecos

#endif
