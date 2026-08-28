/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef NATAF_TRANSFORMATION_HPP
#define NATAF_TRANSFORMATION_HPP

#include "ProbabilityTransformation.hpp"

namespace Pecos {


/// Class for Nataf nonlinear distribution transformation.

/** The Nataf transformation occurs in two steps: (1) transformation
    from the original correlated distributions (x-space) to correlated
    standard normals (z-space) using CDF matching and from correlated
    standard normals to uncorrelated standard normals (u-space) using
    the inverse Cholesky factor of a modified correlation matrix. */

class NatafTransformation: public ProbabilityTransformation
{
public:

  //
  //- Heading: Constructors and destructor
  //

  NatafTransformation();  ///< constructor
  ~NatafTransformation(); ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// Transformation routine from u-space of uncorrelated standard normal
  /// variables to x-space of correlated random variables
  void trans_U_to_X(const RealVector& u_vars, SizetMultiArrayConstView u_cv_ids,
		    RealVector& x_vars, SizetMultiArrayConstView x_cv_ids);

  /// Transformation routine from x-space of correlated random variables 
  /// to u-space of uncorrelated standard normal variables
  void trans_X_to_U(const RealVector& x_vars, SizetMultiArrayConstView x_cv_ids,
		    RealVector& u_vars, SizetMultiArrayConstView u_cv_ids);

  /// As part of the Nataf distribution model (Der Kiureghian & Liu, 1986),
  /// this procedure modifies the user-specified correlation matrix
  /// (corrMatrixX) to account for correlation warping from the nonlinear
  /// X->Z transformation and performs a Cholesky factorization to create
  /// corrCholeskyFactorZ.
  void transform_correlations();

  /// Transformation routine for gradient vector from x-space to u-space
  void trans_grad_X_to_U(const RealVector& fn_grad_x,
			 SizetMultiArrayConstView x_cv_ids,
			 RealVector& fn_grad_u,
			 SizetMultiArrayConstView u_cv_ids,
			 const RealVector& x_vars, const SizetArray& x_dvv);
  /// Transformation routine for gradient vector from x-space to u-space
  void trans_grad_X_to_U(const RealVector& fn_grad_x,
			 SizetMultiArrayConstView x_cv_ids,
			 RealVector& fn_grad_u, const RealMatrix& jacobian_xu,
			 const SizetArray& x_dvv);

  /// Transformation routine from x-space gradient vector to design space
  void trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
			 const RealVector& x_vars, const SizetArray& x_dvv,
			 SizetMultiArrayConstView x_cv_ids,
			 SizetMultiArrayConstView u_cv_ids,
			 SizetMultiArrayConstView x_acv_ids,
			 const SizetArray& x_acv_map1_indices,
			 const ShortArray& x_acv_map2_targets);
  /// Transformation routine from x-space gradient vector to design space
  void trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
			 const RealMatrix& jacobian_xs, const SizetArray& x_dvv,
			 SizetMultiArrayConstView x_cv_ids,
			 SizetMultiArrayConstView u_cv_ids,
			 SizetMultiArrayConstView x_acv_ids,
			 const SizetArray& x_acv_map1_indices,
			 const ShortArray& x_acv_map2_targets);

  /// Transformation routine for gradient vector from u-space to x-space
  void trans_grad_U_to_X(const RealVector& fn_grad_u,
			 SizetMultiArrayConstView u_cv_ids,
			 RealVector& fn_grad_x,
			 SizetMultiArrayConstView x_cv_ids,
			 const RealVector& x_vars, const SizetArray& x_dvv);
  /// Transformation routine for gradient vector from u-space to x-space
  void trans_grad_U_to_X(const RealVector& fn_grad_u,
			 RealVector& fn_grad_x,
			 SizetMultiArrayConstView x_cv_ids,
			 const RealMatrix& jacobian_ux,
			 const SizetArray& x_dvv);

  /// Transformation routine for Hessian matrix from x-space to u-space
  void trans_hess_X_to_U(const RealSymMatrix& fn_hess_x,
			 SizetMultiArrayConstView x_cv_ids,
			 RealSymMatrix& fn_hess_u,
			 SizetMultiArrayConstView u_cv_ids,
			 const RealVector& x_vars, const RealVector& fn_grad_x,
			 const SizetArray& x_dvv);
  /// Transformation routine for Hessian matrix from x-space to u-space
  void trans_hess_X_to_U(const RealSymMatrix& fn_hess_x,
			 SizetMultiArrayConstView x_cv_ids,
			 RealSymMatrix& fn_hess_u,
			 const RealMatrix& jacobian_xu,
			 const RealSymMatrixArray& hessian_xu,
			 const RealVector& fn_grad_x, const SizetArray& x_dvv);

  /// Transformation routine for Hessian matrix from u-space to x-space
  void trans_hess_U_to_X(const RealSymMatrix& fn_hess_u,
			 SizetMultiArrayConstView u_cv_ids,
			 RealSymMatrix& fn_hess_x,
			 SizetMultiArrayConstView x_cv_ids,
			 const RealVector& x_vars, const RealVector& fn_grad_u,
			 const SizetArray& x_dvv);
  /// Transformation routine for Hessian matrix from u-space to x-space
  void trans_hess_U_to_X(const RealSymMatrix& fn_hess_u,
			 RealSymMatrix& fn_hess_x,
			 SizetMultiArrayConstView x_cv_ids,
			 const RealMatrix& jacobian_ux,
			 const RealSymMatrixArray& hessian_ux,
			 const RealVector& fn_grad_u, const SizetArray& x_dvv);

  /// Jacobian of x(u) mapping obtained from dX/dZ dZ/dU
  void jacobian_dX_dU(const RealVector& x_vars,
		      SizetMultiArrayConstView x_cv_ids,
		      SizetMultiArrayConstView u_cv_ids,
		      RealMatrix& jacobian_xu);

  /// Jacobian of u(x) mapping obtained from dU/dZ dZ/dX
  void jacobian_dU_dX(const RealVector& x_vars,
		      SizetMultiArrayConstView x_cv_ids,
		      SizetMultiArrayConstView u_cv_ids,
		      RealMatrix& jacobian_ux);

  /// Design Jacobian of x(u,s) mapping obtained from differentiation of
  /// trans_U_to_X() with respect to distribution parameters S
  void jacobian_dX_dS(const RealVector& x_vars, RealMatrix& jacobian_xs,
		      SizetMultiArrayConstView x_cv_ids,
		      SizetMultiArrayConstView u_cv_ids,
		      SizetMultiArrayConstView x_acv_ids,
		      const SizetArray& x_acv_map1_indices,
		      const ShortArray& x_acv_map2_targets);

  /// Computes numerical dx/ds and dz/ds Jacobians as requested by xs
  /// and zs booleans
  void numerical_design_jacobian(const RealVector& x_vars,
                                 bool xs, RealMatrix& num_jacobian_xs,
                                 bool zs, RealMatrix& num_jacobian_zs,
				 SizetMultiArrayConstView x_cv_ids,
				 SizetMultiArrayConstView u_cv_ids,
				 SizetMultiArrayConstView x_acv_ids,
				 const SizetArray& x_acv_map1_indices,
				 const ShortArray& x_acv_map2_targets);

  /// Hessian of x(u) mapping obtained from dZ/dU^T d^2X/dZ^2 dZ/dU
  void hessian_d2X_dU2(const RealVector& x_vars,
		       SizetMultiArrayConstView x_cv_ids,
		       SizetMultiArrayConstView u_cv_ids,
		       RealSymMatrixArray& hessian_xu);

  /// Hessian of x(u) mapping obtained from dZ/dU^T d^2X/dZ^2 dZ/dU
  void hessian_d2U_dX2(const RealVector& x_vars,
		       SizetMultiArrayConstView x_cv_ids,
		       SizetMultiArrayConstView u_cv_ids,
		       RealSymMatrixArray& hessian_ux);

private:

  //
  //- Heading: Utility routines
  //

  /// Transformation routine from u-space of uncorrelated standard normal
  /// variables to z-space of correlated standard normal variables
  void trans_U_to_Z(const RealVector& u_vars, RealVector& z_vars);

  /// Transformation routine from z-space of correlated standard normal
  /// variables to x-space of correlated random variables
  void trans_Z_to_X(const RealVector& z_vars, SizetMultiArrayConstView u_cv_ids,
		    RealVector& x_vars, SizetMultiArrayConstView x_cv_ids);
  /// Transformation routine from a single z-space variable to a
  /// corresponding x-space variable
  void trans_Z_to_X(Real z, size_t u_rv_index, Real& x, size_t x_rv_index);

  /// Transformation routine from x-space of correlated random variables
  /// to z-space of correlated standard normal variables
  void trans_X_to_Z(const RealVector& x_vars, SizetMultiArrayConstView x_cv_ids,
		    RealVector& z_vars, SizetMultiArrayConstView u_cv_ids);
  /// Transformation routine from a single x-space random variable to
  /// a corresponding z-space variable
  void trans_X_to_Z(Real x, size_t x_rv_index, Real& z, size_t u_rv_index);

  /// Transformation routine from z-space of correlated standard normal
  /// variables to u-space of uncorrelated standard normal variables
  void trans_Z_to_U(RealVector& z_vars, RealVector& u_vars);

  /// Jacobian of x(z) mapping obtained from differentiation of trans_Z_to_X()
  void jacobian_dX_dZ(const RealVector& x_vars,
		      SizetMultiArrayConstView x_cv_ids,
		      SizetMultiArrayConstView u_cv_ids,
		      RealMatrix& jacobian_xz);

  /// Jacobian of z(x) mapping obtained from differentiation of trans_X_to_Z()
  void jacobian_dZ_dX(const RealVector& x_vars,
		      SizetMultiArrayConstView x_cv_ids,
		      SizetMultiArrayConstView u_cv_ids,
		      RealMatrix& jacobian_zx);

  /// Hessian of x(z) mapping obtained from differentiation of jacobian_dX_dZ()
  void hessian_d2X_dZ2(const RealVector& x_vars,
		       SizetMultiArrayConstView x_cv_ids,
		       SizetMultiArrayConstView u_cv_ids,
		       RealSymMatrixArray& hessian_xz);

  /// Hessian of x(z) mapping obtained from differentiation of jacobian_dX_dZ()
  void hessian_d2Z_dX2(const RealVector& x_vars,
		       SizetMultiArrayConstView x_cv_ids,
		       SizetMultiArrayConstView u_cv_ids,
		       RealSymMatrixArray& hessian_zx);

  bool nonlinear_variables_mapping(SizetMultiArrayConstView x_cv_ids,
				   SizetMultiArrayConstView u_cv_ids);

  void compute_reindexing(const SizetArray& dvv,
			  SizetMultiArrayConstView cv_ids,
			  SizetArray& dvv_index_array);
  void gradient_reindex(const SizetArray& dvv_index_array,
			const RealVector& fn_grad,
			RealVector& fn_grad_reindex);
  void hessian_reindex(const SizetArray& dvv_index_array,
		       const RealSymMatrix& fn_hess,
		       RealSymMatrix& fn_hess_reindex);
  void hessian_index_restore(const SizetArray& dvv_index_array,
			     size_t num_deriv_v,
			     const RealSymMatrix& fn_hess_reindex,
			     RealSymMatrix& fn_hess);

  //
  //- Heading: Data
  //

  /// Cholesky factor of the modified correlation matrix; computed in
  /// transform_correlations().
  /** Note that this is not a component of a MarginalsCorrDistribution, but
      is rather a product of the Nataf transformation since it is defined
      by the original correlation matrix + transformation context/targets. */
  RealMatrix corrCholeskyFactorZ;
};


inline NatafTransformation::NatafTransformation():
  ProbabilityTransformation(BaseConstructor())
{ }


inline NatafTransformation::~NatafTransformation()
{ }


inline bool NatafTransformation::
nonlinear_variables_mapping(SizetMultiArrayConstView x_cv_ids,
			    SizetMultiArrayConstView u_cv_ids)
{
  size_t i, num_v = x_cv_ids.size(), x_cv_index, u_cv_index;
  short x_type, u_type;
  for (i=0; i<num_v; ++i) {
    x_cv_index = x_cv_ids[i] - 1;  u_cv_index = u_cv_ids[i] - 1;
    x_type = xDist.random_variable_type(x_cv_index);
    u_type = uDist.random_variable_type(u_cv_index);
    if ( ( ( x_type == CONTINUOUS_RANGE || x_type == UNIFORM ||
	     x_type == CONTINUOUS_INTERVAL_UNCERTAIN ) &&
	   u_type   != STD_UNIFORM ) ||
	 ( x_type   == NORMAL      && u_type != STD_NORMAL ) ||
	 ( x_type   == EXPONENTIAL && u_type != STD_EXPONENTIAL ) ||
	 ( x_type   == BETA        && u_type != STD_BETA   ) ||
	 ( x_type   == GAMMA       && u_type != STD_GAMMA  ) ||
	 ( ( x_type == BOUNDED_NORMAL    || x_type == LOGNORMAL  ||
	     x_type == BOUNDED_LOGNORMAL || x_type == LOGUNIFORM ||
	     x_type == TRIANGULAR        || x_type == GUMBEL     ||
	     x_type == FRECHET           || x_type == WEIBULL ) &&
	   x_type   != u_type ) ||
	 ( x_type   == HISTOGRAM_BIN && u_type != STD_UNIFORM &&
	   x_type   != u_type ) )
      return true;
  }
  return false;
}


inline void NatafTransformation::
compute_reindexing(const SizetArray& dvv, SizetMultiArrayConstView cv_ids,
		   SizetArray& dvv_index_array)
{
  size_t i, num_v = cv_ids.size();
  dvv_index_array.resize(num_v);
  for (i=0; i<num_v; ++i)
    dvv_index_array[i] = find_index(dvv, cv_ids[i]);
}

inline void NatafTransformation::
gradient_reindex(const SizetArray& dvv_index_array, const RealVector& fn_grad,
		 RealVector& fn_grad_reindex)
{
  size_t i, dvv_index, num_v = dvv_index_array.size();
  if (fn_grad_reindex.length() != num_v) fn_grad_reindex.size(num_v);
  else                                   fn_grad_reindex = 0.;
  for (i=0; i<num_v; ++i) {
    dvv_index = dvv_index_array[i];
    if (dvv_index != _NPOS)
      fn_grad_reindex[i] = fn_grad(dvv_index);
  }
}


inline void NatafTransformation::
hessian_reindex(const SizetArray& dvv_index_array, const RealSymMatrix& fn_hess,
		RealSymMatrix& fn_hess_reindex)
{
  size_t i, j, dvv_index_i, dvv_index_j, num_v = dvv_index_array.size();
  if (fn_hess_reindex.numRows() != num_v) fn_hess_reindex.shape(num_v);
  else                                    fn_hess_reindex = 0.;
  for (i=0; i<num_v; ++i) {
    dvv_index_i = dvv_index_array[i];
    if (dvv_index_i != _NPOS) {
      for (j=0; j<num_v; ++j) {
	dvv_index_j = dvv_index_array[j];
	if (dvv_index_j != _NPOS)
	  fn_hess_reindex(i, j) = fn_hess(dvv_index_i, dvv_index_j);
      }
    }
  }
}


inline void NatafTransformation::
hessian_index_restore(const SizetArray& dvv_index_array, size_t num_deriv_v,
		      const RealSymMatrix& fn_hess_reindex,
		      RealSymMatrix& fn_hess)
{
  if (fn_hess.numRows() != num_deriv_v) fn_hess.shape(num_deriv_v);
  else                                  fn_hess = 0.;
  size_t i, j, dvv_index_i, dvv_index_j, num_v = dvv_index_array.size();
  for (i=0; i<num_v; ++i) {
    dvv_index_i = dvv_index_array[i];
    if (dvv_index_i != _NPOS) {
      for (j=0; j<num_v; ++j) {
	dvv_index_j = dvv_index_array[j];
	if (dvv_index_j != _NPOS)
	  fn_hess(dvv_index_i, dvv_index_j) = fn_hess_reindex(i, j);
      }
    }
  }
}

} // namespace Pecos

#endif
