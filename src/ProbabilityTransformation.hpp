/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PROBABILITY_TRANSFORMATION_HPP
#define PROBABILITY_TRANSFORMATION_HPP

#include "pecos_data_types.hpp"

namespace Pecos {


/// Base class for all nonlinear distribution transformations

/** The base class for nonlinear distribution transformations,
    including Nataf, Rosenblatt, et al. */

class ProbabilityTransformation
{
public:

  /// default constructor
  ProbabilityTransformation();
  /// standard constructor for envelope
  ProbabilityTransformation(const String& prob_trans_type);
  /// copy constructor
  ProbabilityTransformation(const ProbabilityTransformation& prob_trans);

  /// destructor
  virtual ~ProbabilityTransformation();

  /// assignment operator
  ProbabilityTransformation
    operator=(const ProbabilityTransformation& prob_trans);

  //
  //- Heading: Virtual functions
  //

  /// Transformation routine from u-space of uncorrelated standard normal
  /// variables to x-space of correlated random variables
  virtual void trans_U_to_X(const RealVector& u_vars, RealVector& x_vars);

  /// Transformation routine from x-space of correlated random variables 
  /// to u-space of uncorrelated standard normal variables
  virtual void trans_X_to_U(const RealVector& x_vars, RealVector& u_vars);

  /// As part of the Nataf distribution model (Der Kiureghian & Liu, 1986),
  /// this procedure modifies the user-specified correlation matrix
  /// (corrMatrixX) to account for correlation warping from the nonlinear
  /// X->Z transformation and performs a Cholesky factorization to create
  /// corrCholeskyFactorZ.
  virtual void trans_correlations();

  /// Transformation routine for gradient vector from x-space to u-space
  virtual void trans_grad_X_to_U(const RealVector& fn_grad_x,
				 RealVector& fn_grad_u,
				 const RealVector& x_vars,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);
  /// Transformation routine for gradient vector from x-space to u-space
  virtual void trans_grad_X_to_U(const RealVector& fn_grad_x,
				 RealVector& fn_grad_u,
				 const RealMatrix& jacobian_xu,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);

  /// Transformation routine from x-space gradient vector to design space
  virtual void trans_grad_X_to_S(const RealVector& fn_grad_x,
				 RealVector& fn_grad_s,
				 const RealVector& x_vars,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids,
				 SizetMultiArrayConstView acv_ids,
				 const SizetArray& acv_map1_indices,
				 const ShortArray& acv_map2_targets);
  /// Transformation routine from x-space gradient vector to design space
  virtual void trans_grad_X_to_S(const RealVector& fn_grad_x,
				 RealVector& fn_grad_s,
				 const RealMatrix& jacobian_xs,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids,
				 SizetMultiArrayConstView acv_ids,
				 const SizetArray& acv_map1_indices,
				 const ShortArray& acv_map2_targets);

  /// Transformation routine for gradient vector from u-space to x-space
  virtual void trans_grad_U_to_X(const RealVector& fn_grad_u,
				 RealVector& fn_grad_x,
				 const RealVector& x_vars,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);
  /// Transformation routine for gradient vector from u-space to x-space
  virtual void trans_grad_U_to_X(const RealVector& fn_grad_u,
				 RealVector& fn_grad_x,
				 const RealMatrix& jacobian_ux,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);

  /// Transformation routine for Hessian matrix from x-space to u-space
  virtual void trans_hess_X_to_U(const RealSymMatrix& fn_hess_x,
				 RealSymMatrix& fn_hess_u,
				 const RealVector& x_vars,
				 const RealVector& fn_grad_x,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);
  /// Transformation routine for Hessian matrix from x-space to u-space
  virtual void trans_hess_X_to_U(const RealSymMatrix& fn_hess_x,
				 RealSymMatrix& fn_hess_u,
				 const RealMatrix& jacobian_xu,
				 const RealSymMatrixArray& hessian_xu,
				 const RealVector& fn_grad_x,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);

  /// Jacobian of x(u) mapping obtained from dX/dZ dZ/dU
  virtual void jacobian_dX_dU(const RealVector& x_vars,
			      RealMatrix& jacobian_xu);

  /// Jacobian of u(x) mapping obtained from dU/dZ dZ/dX
  virtual void jacobian_dU_dX(const RealVector& x_vars,
			      RealMatrix& jacobian_ux);

  /// Design Jacobian of x(u,s) mapping obtained from differentiation of
  /// trans_U_to_X() with respect to distribution parameters S
  virtual void jacobian_dX_dS(const RealVector& x_vars, RealMatrix& jacobian_xs,
			      SizetMultiArrayConstView cv_ids,
			      SizetMultiArrayConstView acv_ids,
			      const SizetArray& acv_map1_indices,
			      const ShortArray& acv_map2_targets);

  /// Hessian of x(u) mapping obtained from dZ/dU^T d^2X/dZ^2 dZ/dU
  virtual void hessian_d2X_dU2(const RealVector& x_vars,
			       RealSymMatrixArray& hessian_xu);

  //
  //- Heading: Member functions
  //

  /// set ranVarTypesX/U, ranVarMeansX, ranVarStdDevsX, ranVarLowerBndsX,
  /// ranVarUpperBndsX, ranVarAddtlParamsX, corrMatrixX, and correlationFlagX
  /// based on incoming data
  void initialize_random_variables(const ProbabilityTransformation& prob_trans);
  /// initializes ranVarTypesX and ranVarTypesU
  void initialize_random_variable_types(const ShortArray& x_types,
					const ShortArray& u_types);
  /// initializes ranVarMeansX, ranVarStdDevsX, ranVarLowerBndsX,
  /// ranVarUpperBndsX, and ranVarAddtlParamsX
  void initialize_random_variable_parameters(const RealVector& x_means,
					     const RealVector& x_std_devs,
					     const RealVector& x_l_bnds,
					     const RealVector& x_u_bnds,
					     const RealVectorArray& x_addtl);
  /// initializes corrMatrixX and correlationFlagX
  void initialize_random_variable_correlations(const RealSymMatrix& x_corr);

  /// reshape corrMatrixX for an all_variables specification
  void reshape_correlation_matrix(size_t num_leading_vars,
				  size_t num_probabilistic_vars,
				  size_t num_trailing_vars);

  /// return ranVarTypesX
  const ShortArray& x_types() const;
  /// return ranVarTypesU
  const ShortArray& u_types() const;
  /// return ranVarMeansX
  const RealVector& x_means() const;
  /// return ranVarStdDevsX
  const RealVector& x_std_deviations() const;
  /// return ranVarLowerBndsX
  const RealVector& x_lower_bounds() const;
  /// return ranVarUpperBndsX
  const RealVector& x_upper_bounds() const;
  /// return ranVarAddtlParamsX
  const RealVectorArray& x_additional_parameters() const;
  /// return correlationFlagX
  bool x_correlation() const;
  /// return corrMatrixX
  const RealSymMatrix& x_correlation_matrix() const;
  /// return corrCholeskyFactorZ
  const RealMatrix& z_correlation_factor() const;

  /// function to check modelRep (does this envelope contain a letter)
  bool is_null() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  ProbabilityTransformation(BaseConstructor);

  //
  //- Heading: Member functions
  //

  /// Computes numerical dx/ds and dz/ds Jacobians as requested by xs
  /// and zs booleans
  void numerical_design_jacobian(const RealVector& x_vars,
                                 bool xs, RealMatrix& num_jacobian_xs,
                                 bool zs, RealMatrix& num_jacobian_zs,
				 SizetMultiArrayConstView cv_ids,
				 SizetMultiArrayConstView acv_ids,
				 const SizetArray& acv_map1_indices,
				 const ShortArray& acv_map2_targets);

#ifdef DERIV_DEBUG
  /// routine for verification of transformation Jacobian/Hessian terms
  void verify_trans_jacobian_hessian(const RealVector& v0);

  /// routine for verification of design Jacobian terms
  void verify_design_jacobian(const RealVector& u0);
#endif // DERIV_DEBUG

  //
  //- Heading: Data members
  //

  // the following attributes are the required data for performing
  // transformations from X -> Z -> U and back.

  /// vector of indices indicating the type of each x-space uncertain variable
  ShortArray ranVarTypesX;
  /// vector of indices indicating the type of standard uncertain variable to
  /// which each x-space variable is transformed
  ShortArray ranVarTypesU;
  /// vector of means for all x-space uncertain variables
  RealVector ranVarMeansX;
  /// vector of standard deviations for all x-space uncertain variables
  RealVector ranVarStdDevsX;
  /// vector of distribution lower bounds for selected x-space uncertain vars
  RealVector ranVarLowerBndsX;
  /// vector of distribution upper bounds for selected x-space uncertain vars
  RealVector ranVarUpperBndsX;
  /// vector of additional distribution parameters (e.g., alphas, betas, modes)
  /// for selected x-space uncertain variables
  RealVectorArray ranVarAddtlParamsX;
  /// flag for indicating if correlation exists among the x-space
  /// uncertain variables
  bool correlationFlagX;
  /// matrix of random variable correlation coefficients
  RealSymMatrix corrMatrixX;
  /// cholesky factor of a modified correlation matrix (#corrMatrixX
  /// is modified in trans_correlations() for use in z-space)
  RealMatrix corrCholeskyFactorZ;

private:

  //
  //- Heading: Member functions
  //

  /// return a particular random variable distribution parameter
  const Real& distribution_parameter(size_t index, short target);
  /// set a particular random variable distribution parameter and
  /// update derived quantities
  void distribution_parameter(size_t index, short target, const Real& param);

  /// Used only by the standard envelope constructor to initialize
  /// probTransRep to the appropriate derived type.
  static ProbabilityTransformation*
    get_prob_trans(const String& prob_trans_type);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  ProbabilityTransformation* probTransRep;
  /// number of objects sharing probTransRep
  int referenceCount;
};


inline const ShortArray& ProbabilityTransformation::x_types() const
{ return (probTransRep) ? probTransRep->ranVarTypesX : ranVarTypesX; }


inline const ShortArray& ProbabilityTransformation::u_types() const
{ return (probTransRep) ? probTransRep->ranVarTypesU : ranVarTypesU; }


inline const RealVector& ProbabilityTransformation::x_means() const
{ return (probTransRep) ? probTransRep->ranVarMeansX : ranVarMeansX; }


inline const RealVector& ProbabilityTransformation::x_std_deviations() const
{ return (probTransRep) ? probTransRep->ranVarStdDevsX : ranVarStdDevsX; }


inline const RealVector& ProbabilityTransformation::x_lower_bounds() const
{ return (probTransRep) ? probTransRep->ranVarLowerBndsX : ranVarLowerBndsX; }


inline const RealVector& ProbabilityTransformation::x_upper_bounds() const
{ return (probTransRep) ? probTransRep->ranVarUpperBndsX : ranVarUpperBndsX; }


inline const RealVectorArray& ProbabilityTransformation::
x_additional_parameters() const
{
  return (probTransRep) ? probTransRep->ranVarAddtlParamsX : ranVarAddtlParamsX;
}


inline bool ProbabilityTransformation::x_correlation() const
{ return (probTransRep) ? probTransRep->correlationFlagX : correlationFlagX; }


inline const RealSymMatrix& ProbabilityTransformation::
x_correlation_matrix() const
{ return (probTransRep) ? probTransRep->corrMatrixX : corrMatrixX; }


inline const RealMatrix& ProbabilityTransformation::z_correlation_factor() const
{
  return (probTransRep) ? probTransRep->corrCholeskyFactorZ :
                          corrCholeskyFactorZ;
}


inline bool ProbabilityTransformation::is_null() const
{ return (probTransRep) ? false : true; }

} // namespace Pecos

#endif
