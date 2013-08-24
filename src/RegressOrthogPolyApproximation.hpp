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
  RegressOrthogPolyApproximation(const UShortArray& approx_order,
				 size_t num_vars, bool use_derivs,
				 short output_level);
  /// destructor
  ~RegressOrthogPolyApproximation();

  //
  //- Heading: Member functions
  //

  /// set crossValidation flag
  void cross_validation(bool flag);
  /// set the noise tolerance(s) for compressed sensing approaches
  void noise_tolerance(const RealVector& noise_tol);
  /// set the L2 penalty parameter for LASSO (elastic net variant)
  void l2_penalty(Real l2_pen);

  /// store the fault info about the response data
  FaultInfo faultInfo;

protected:

  //
  //- Heading: Virtual function redefinitions and member functions
  //

  int min_coefficients() const;
  void compute_coefficients();

  void allocate_arrays();

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
  void run_cross_validation( RealMatrix &A, RealMatrix &B, 
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

  /// pack polynomial contributions to Psi matrix for regression
  void pack_polynomial_data(const RealVector& c_vars, const UShortArray& mi,
			    bool add_val,  double* pack_val,  size_t& pv_cntr,
			    bool add_grad, double* pack_grad, size_t& pg_cntr);
  /// pack response contributions to RHS for regression
  void pack_response_data(const SurrogateDataResp& sdr,
			  bool add_val,  double* pack_val,  size_t& pv_cntr,
			  bool add_grad, double* pack_grad, size_t& pg_cntr);

  /// define multiIndex and expansionCoeffs from nonzero dense_coeffs
  void update_sparse(Real* dense_coeffs, size_t num_dense_terms);
  /// augment sparse_indices based on nonzero dense_coeffs
  void update_sparse_indices(Real* dense_coeffs, size_t num_dense_terms,
			     SizetSet& sparse_indices);
  /// define sparse multiIndex from sparse_indices
  void update_sparse_multi_index(const SizetSet& sparse_indices);
  /// define sparse expansionCoeffs from dense_coeffs and sparse_indices
  void update_sparse_coeffs(Real* dense_coeffs, const SizetSet& sparse_indices);
  /// define a row of sparse expansionCoeffGrads from dense_coeffs and
  /// sparse_indices
  void update_sparse_coeff_grads(Real* dense_coeffs, int row,
				 const SizetSet& sparse_indices);

  /**
   * \brief Define the set of options used in the cross validation grid search
   * 
   * \param opts (output) the options to be used in the grid search
   * \param M The number of rows of the vandermonde matrix
   * \param N The number of columns of the vandermonde matrix
   */
  void gridSearchFunction( RealMatrix &opts, int M, int N, 
			   int num_function_samples );

  //
  //- Heading: Data
  //

  /// order of orthogonal best polynomial expansion found using cross validation
  IntVector bestApproxOrder;

  /// Wrapper class that is used to solve regression problems
  CompressedSensingTool CSTool;
  /// Stuct use to define the options of a compressed sensing solve
  CompressedSensingOptions CSOpts;

  /// flag for use of automatic cross-validation for parameter
  /// selection in regression approaches
  bool crossValidation;
  /// noise tolerance(s) for compressed sensing algorithms; vector form
  /// used in cross-validation
  RealVector noiseTols;
  /// the L2 penalty parameter for LASSO (elastic net variant)
  Real l2Penalty;
};


inline RegressOrthogPolyApproximation::
RegressOrthogPolyApproximation(const UShortArray& approx_order, size_t num_vars,
			       bool use_derivs, short output_level):
  OrthogPolyApproximation(approx_order, num_vars, use_derivs, output_level)
{ }


inline RegressOrthogPolyApproximation::~RegressOrthogPolyApproximation()
{ }


inline void RegressOrthogPolyApproximation::cross_validation(bool flag)
{ crossValidation = flag; }


inline void RegressOrthogPolyApproximation::
noise_tolerance(const RealVector& noise_tol)
{ noiseTols = noise_tol; }


inline void RegressOrthogPolyApproximation::l2_penalty(Real l2_pen)
{ l2Penalty = l2_pen; }


inline void RegressOrthogPolyApproximation::
pack_polynomial_data(const RealVector& c_vars, const UShortArray& mi,
		     bool add_val,  double* pack_val,  size_t& pv_cntr,
		     bool add_grad, double* pack_grad, size_t& pg_cntr)
{
  if (add_val)
    { pack_val[pv_cntr] = multivariate_polynomial(c_vars, mi); ++pv_cntr; }
  if (add_grad) {
    const RealVector& mvp_grad
      = multivariate_polynomial_gradient_vector(c_vars, mi);
    for (size_t j=0; j<numVars; ++j, ++pg_cntr)
      pack_grad[pg_cntr] = mvp_grad[j];
  }
}


inline void RegressOrthogPolyApproximation::
pack_response_data(const SurrogateDataResp& sdr,
		   bool add_val,  double* pack_val,  size_t& pv_cntr,
		   bool add_grad, double* pack_grad, size_t& pg_cntr)
{
  if (add_val)
    { pack_val[pv_cntr] = sdr.response_function(); ++pv_cntr; }
  if (add_grad) {
    const RealVector& resp_grad = sdr.response_gradient();
    for (size_t j=0; j<numVars; ++j, ++pg_cntr)
      pack_grad[pg_cntr] = resp_grad[j];
  }
}

} // namespace Pecos

#endif
