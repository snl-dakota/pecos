/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        RegressOrthogPolyApproximation
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        John Jakeman

#ifndef REGRESS_ORTHOG_POLY_APPROXIMATION_HPP
#define REGRESS_ORTHOG_POLY_APPROXIMATION_HPP

#include "OrthogPolyApproximation.hpp"
#include "CompressedSensing.hpp"
#include "FaultTolerance.hpp"

namespace Pecos {

/// Derived approximation class for multivariate orthogonal polynomial
/// approximation with coefficient estimation via regression.

/** The RegressOrthogPolyApproximation class provides a global
    approximation based on multivariate orthogonal polynomials, where
    the coefficients are computed using regression approaches such as
    least squares (L2) or compressed sensing (L1).  It is used
    primarily for polynomial chaos expansion aproaches to UQ. */

class RegressOrthogPolyApproximation: public OrthogPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  RegressOrthogPolyApproximation(const SharedBasisApproxData& shared_data);
  /// destructor
  ~RegressOrthogPolyApproximation();

  //
  //- Heading: Member functions
  //

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  int min_coefficients() const;
  void compute_coefficients();
  void allocate_arrays();
  void allocate_component_sobol();

  size_t expansion_terms() const;

  Real value(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x,
					     const SizetArray& dvv);
  const RealVector& gradient_nonbasis_variables(const RealVector& x);

  Real stored_value(const RealVector& x);
  const RealVector& stored_gradient_basis_variables(const RealVector& x);
  const RealVector& stored_gradient_nonbasis_variables(const RealVector& x);

  Real mean(const RealVector& x);
  const RealVector& mean_gradient(const RealVector& x, const SizetArray& dvv);

  Real covariance(PolynomialApproximation* poly_approx_2);
  Real covariance(const RealVector& x, PolynomialApproximation* poly_approx_2);
  const RealVector& variance_gradient();
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

private:

  //
  //- Heading: Member functions
  //

  /// selects the solver for L1 or L2 minimization based on user input
  void select_solver();

  /// Run the regression method set in select_solver() to compute the
  /// expansion coefficients using L1 or L2 minimization
  void run_regression();

  /// set the information needed to ensure fault tolerance
  void set_fault_info();

  /// Use cross validation to find the 'best' PCE degree
  void run_cross_validation( RealMatrix &A, RealMatrix &B, RealMatrix &points,
			     size_t num_data_pts_fn );

  /// For a specific vandermonde matrix find the compressed sennsing 
  // options that produce the best PCE
  void estimate_compressed_sensing_options_via_cross_validation(
    RealMatrix &vandermonde_matrix, RealMatrix &rhs,
    std::vector<CompressedSensingOptions> &best_cs_opts,
    RealVector &best_predictor_indicators,
    RealMatrixArray &predictor_options_history,
    RealMatrixArray &predictor_indicators_history,
    RealMatrixArray &predictor_partition_indicators_history,
    size_t num_data_pts_fn );

  /// define multiIndex and expansionCoeffs from nonzero dense_coeffs
  void update_sparse(Real* dense_coeffs, size_t num_dense_terms);
  /// augment sparse_indices based on nonzero dense_coeffs
  void update_sparse_indices(Real* dense_coeffs, size_t num_dense_terms);
  /// define sparse expansionCoeffs from dense_coeffs and sparse_indices
  void update_sparse_coeffs(Real* dense_coeffs);
  /// define a row of sparse expansionCoeffGrads from dense_coeffs and
  /// sparse_indices
  void update_sparse_coeff_grads(Real* dense_coeffs, int row);
  // define sparse multiIndex from sparseIndices
  //void update_sparse_multi_index();
  /// define sparseSobolIndices from sparseIndices and multiIndex
  void update_sparse_sobol();

  /**
   * \brief Define the set of options used in the cross validation grid search
   * 
   * \param opts (output) the options to be used in the grid search
   * \param M The number of rows of the vandermonde matrix
   * \param N The number of columns of the vandermonde matrix
   */
  void gridSearchFunction( RealMatrix &opts, int M, int N, 
			   int num_function_samples );

  void least_interpolation( RealMatrix &pts, 
			    RealMatrix &vals );

  void transform_least_interpolant( RealMatrix &L,
				    RealMatrix &U,
				    RealMatrix &H,
				    IntVector &p,
				    RealVector &v,
				    RealMatrix &vals );

  void least_factorization( RealMatrix &x,
			    UShort2DArray &basis_indices,
			    RealMatrix &l,
			    RealMatrix &u, 
			    RealMatrix &H, 
			    IntVector &p,
			    RealVector &v, 
			    IntVector &k );

  void get_least_polynomial_coefficients( RealVector &v, IntVector &k,
					  UShort2DArray &basis_indices,
					  int num_dims, int num_pts,
					  RealMatrix &H );

  //
  //- Heading: Data
  //

  /// order of orthogonal best polynomial expansion found using cross validation
  IntVector bestApproxOrder;

  /// Stuct use to define the options of a compressed sensing solve
  CompressedSensingOptions CSOpts;

  /// store the fault info about the response data
  FaultInfo faultInfo;

  /// tracks sparse terms within multiIndex and expansion{Coeffs,CoeffGrads}
  /// that are retained from an original candidate set
  SizetSet sparseIndices;
  /// copy of sparseIndices stored in store_coefficients() for use in
  /// combine_coefficients()
  SizetSet storedSparseIndices;
  /// tracks sparse Sobol index terms, as a subset of sobolIndexMap
  SizetSet sparseSobolIndices;
};


inline RegressOrthogPolyApproximation::
RegressOrthogPolyApproximation(const SharedBasisApproxData& shared_data):
  OrthogPolyApproximation(shared_data)
{ }


inline RegressOrthogPolyApproximation::~RegressOrthogPolyApproximation()
{ }


inline size_t RegressOrthogPolyApproximation::expansion_terms() const
{
  return (sparseIndices.empty()) ? OrthogPolyApproximation::expansion_terms()
                                 : sparseIndices.size();
}

} // namespace Pecos

#endif
