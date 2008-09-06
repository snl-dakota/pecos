/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 Transformation
//- Description: Base class for nonlinear distribution transformations
//- Owner:	 Mike Eldred
//- Checked by:
//- Version:

#ifndef TRANSFORMATION_HPP
#define TRANSFORMATION_HPP

#include "pecos_global_defs.hpp"
#ifdef HAVE_GSL
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_gamma.h"
#endif

#include <cmath>


namespace Pecos {


/// Base class for all nonlinear distribution transformations

/** The base class for nonlinear distribution transformations,
    including Nataf, Rosenblatt, et al. */

class Transformation
{
public:

  /// default constructor
  Transformation();
  /// standard constructor for envelope
  Transformation(const String& trans_type);
  /// copy constructor
  Transformation(const Transformation& trans);

  /// destructor
  virtual ~Transformation();

  /// assignment operator
  Transformation operator=(const Transformation& trans);

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
				 const UIntArray& x_dvv,
				 const UIntArray&  cv_ids);
  /// Transformation routine for gradient vector from x-space to u-space
  virtual void trans_grad_X_to_U(const RealVector& fn_grad_x,
				 RealVector& fn_grad_u,
				 const RealMatrix& jacobian_xu,
				 const UIntArray& x_dvv,
				 const UIntArray& cv_ids);

  /// Transformation routine from x-space gradient vector to design space
  virtual void trans_grad_X_to_S(const RealVector& fn_grad_x,
				 RealVector& fn_grad_s,
				 const RealVector& x_vars,
				 const UIntArray& x_dvv,
				 const UIntArray& cv_ids,
				 const UIntArray& acv_ids,
				 const SizetArray& acv_map1_indices,
				 const ShortArray& acv_map2_targets);
  /// Transformation routine from x-space gradient vector to design space
  virtual void trans_grad_X_to_S(const RealVector& fn_grad_x,
				 RealVector& fn_grad_s,
				 const RealMatrix& jacobian_xs,
				 const UIntArray& x_dvv,
				 const UIntArray& cv_ids,
				 const UIntArray& acv_ids,
				 const SizetArray& acv_map1_indices,
				 const ShortArray& acv_map2_targets);

  /// Transformation routine for gradient vector from u-space to x-space
  virtual void trans_grad_U_to_X(const RealVector& fn_grad_u,
				 RealVector& fn_grad_x,
				 const RealVector& x_vars,
				 const UIntArray& x_dvv,
				 const UIntArray& cv_ids);
  /// Transformation routine for gradient vector from u-space to x-space
  virtual void trans_grad_U_to_X(const RealVector& fn_grad_u,
				 RealVector& fn_grad_x,
				 const RealMatrix& jacobian_ux,
				 const UIntArray& x_dvv,
				 const UIntArray& cv_ids);

  /// Transformation routine for Hessian matrix from x-space to u-space
  virtual void trans_hess_X_to_U(const RealSymMatrix& fn_hess_x,
				 RealSymMatrix& fn_hess_u,
				 const RealVector& x_vars,
				 const RealVector& fn_grad_x,
				 const UIntArray& x_dvv,
				 const UIntArray& cv_ids);
  /// Transformation routine for Hessian matrix from x-space to u-space
  virtual void trans_hess_X_to_U(const RealSymMatrix& fn_hess_x,
				 RealSymMatrix& fn_hess_u,
				 const RealMatrix& jacobian_xu,
				 const RealSymMatrixArray& hessian_xu,
				 const RealVector& fn_grad_x,
				 const UIntArray& x_dvv,
				 const UIntArray& cv_ids);

  /// Jacobian of x(u) mapping obtained from dX/dZ dZ/dU
  virtual void jacobian_dX_dU(const RealVector& x_vars,
			      RealMatrix& jacobian_xu);

  /// Jacobian of u(x) mapping obtained from dU/dZ dZ/dX
  virtual void jacobian_dU_dX(const RealVector& x_vars,
			      RealMatrix& jacobian_ux);

  /// Design Jacobian of x(u,s) mapping obtained from differentiation of
  /// trans_U_to_X() with respect to distribution parameters S
  virtual void jacobian_dX_dS(const RealVector& x_vars, RealMatrix& jacobian_xs,
			      const UIntArray& cv_ids, const UIntArray& acv_ids,
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
  void initialize_random_variables(const Transformation& trans);
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
  void reshape_correlation_matrix(size_t num_design_vars,
				  size_t num_uncertain_vars,
				  size_t num_state_vars);

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

  /// Standard normal density function
  Real phi(const Real& beta);
  /// Standard normal cumulative distribution function
  Real Phi(const Real& beta);
  /// Inverse of standard normal cumulative distribution function
  Real Phi_inverse(const Real& p);

  /// compute std deviation from lognormal error factor specification
  void moments_from_lognormal_params(const Real& mean, const Real& err_fact,
				     Real& std_dev);
  /// compute mean and std deviation from uniform bounds specification
  void moments_from_uniform_params(const Real& lwr, const Real& upr, Real& mean,
				   Real& std_dev);
  /// compute mean and std deviation from loguniform bounds specification
  void moments_from_loguniform_params(const Real& lwr, const Real& upr,
				      Real& mean, Real& std_dev);
  /// compute mean and std deviation from triangular mode/bounds specification
  void moments_from_triangular_params(const Real& lwr, const Real& upr,
				      const Real& mode, Real& mean,
				      Real& std_dev);
  /// compute mean and std deviation from exponential beta specification
  void moments_from_exponential_params(const Real& beta, Real& mean,
				       Real& std_dev);
  /// compute mean and std deviation from beta parameter specification
  void moments_from_beta_params(const Real& lwr, const Real& upr,
				const Real& alpha, const Real& beta,
				Real& mean, Real& std_dev);
  /// compute mean and std deviation from gamma parameter specification
  void moments_from_gamma_params(const Real& alpha, const Real& beta,
				 Real& mean, Real& std_dev);
  /// compute mean and std deviation from gumbel parameter specification
  void moments_from_gumbel_params(const Real& alpha, const Real& beta,
				  Real& mean, Real& std_dev);
  /// compute mean and std deviation from frechet parameter specification
  void moments_from_frechet_params(const Real& alpha, const Real& beta,
				   Real& mean, Real& std_dev);
  /// compute mean and std deviation from weibull parameter specification
  void moments_from_weibull_params(const Real& alpha, const Real& beta,
				   Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  Transformation(BaseConstructor);

  //
  //- Heading: Member functions
  //

  /// Computes numerical dx/ds and dz/ds Jacobians as requested by xs
  /// and zs booleans
  void numerical_design_jacobian(const RealVector& x_vars,
                                 bool xs, RealMatrix& num_jacobian_xs,
                                 bool zs, RealMatrix& num_jacobian_zs,
				 const UIntArray& cv_ids,
				 const UIntArray& acv_ids,
				 const SizetArray& acv_map1_indices,
				 const ShortArray& acv_map2_targets);

#ifdef DERIV_DEBUG
  /// routine for verification of transformation Jacobian/Hessian terms
  void verify_trans_jacobian_hessian(const RealVector& v0);

  /// routine for verification of design Jacobian terms
  void verify_design_jacobian(const RealVector& u0);
#endif // DERIV_DEBUG

#ifdef HAVE_GSL
  /// Inverse of standard beta CDF (not supported by GSL)
  Real cdf_beta_Pinv(const Real& normcdf, const Real& alpha, const Real& beta);
#endif // HAVE_GSL

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

  /// the value for Pi used in several numerical routines
  static const Real Pi;

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
  /// transRep to the appropriate derived type.
  Transformation* get_trans(const String& trans_type);

#ifndef HAVE_GSL
  /// Inverse of error function used in Phi_inverse()
  Real erf_inverse(const Real& p);
#endif // HAVE_GSL

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  Transformation* transRep;
  /// number of objects sharing transRep
  int referenceCount;
};


inline const ShortArray& Transformation::x_types() const
{ return (transRep) ? transRep->ranVarTypesX : ranVarTypesX; }


inline const ShortArray& Transformation::u_types() const
{ return (transRep) ? transRep->ranVarTypesU : ranVarTypesU; }


inline const RealVector& Transformation::x_means() const
{ return (transRep) ? transRep->ranVarMeansX : ranVarMeansX; }


inline const RealVector& Transformation::x_std_deviations() const
{ return (transRep) ? transRep->ranVarStdDevsX : ranVarStdDevsX; }


inline const RealVector& Transformation::x_lower_bounds() const
{ return (transRep) ? transRep->ranVarLowerBndsX : ranVarLowerBndsX; }


inline const RealVector& Transformation::x_upper_bounds() const
{ return (transRep) ? transRep->ranVarUpperBndsX : ranVarUpperBndsX; }


inline const RealVectorArray& Transformation::x_additional_parameters() const
{ return (transRep) ? transRep->ranVarAddtlParamsX : ranVarAddtlParamsX; }


inline bool Transformation::x_correlation() const
{ return (transRep) ? transRep->correlationFlagX : correlationFlagX; }


inline const RealSymMatrix& Transformation::x_correlation_matrix() const
{ return (transRep) ? transRep->corrMatrixX : corrMatrixX; }


inline const RealMatrix& Transformation::z_correlation_factor() const
{ return (transRep) ? transRep->corrCholeskyFactorZ : corrCholeskyFactorZ; }


inline Real Transformation::phi(const Real& beta)
{
#ifdef HAVE_GSL
  return gsl_ran_ugaussian_pdf(beta);
#else
  return exp(-beta*beta/2.)/sqrt(2.*Pi);
#endif // HAVE_GSL
}


/** returns a probability < 0.5 for negative beta and a probability > 0.5
    for positive beta. */
inline Real Transformation::Phi(const Real& beta)
{
#ifdef HAVE_GSL
  return gsl_cdf_ugaussian_P(beta);
#else
  return .5 + .5*erf(beta/sqrt(2.));
#endif // HAVE_GSL
}


/** returns a negative beta for probability < 0.5 and a positive beta for
    probability > 0.5. */
inline Real Transformation::Phi_inverse(const Real& p)
{
#ifdef HAVE_GSL
  return gsl_cdf_ugaussian_Pinv(p);
#else
  return sqrt(2.)*erf_inverse(2.*p - 1.);
#endif // HAVE_GSL
}


inline void Transformation::
moments_from_lognormal_params(const Real& mean, const Real& err_fact,
			      Real& std_dev)
{
  Real zeta = log(err_fact)/1.645;
  std_dev   = mean*sqrt(exp(zeta*zeta)-1.);
}


inline void Transformation::
moments_from_uniform_params(const Real& lwr, const Real& upr, Real& mean,
			    Real& std_dev)
{ mean = (lwr + upr)/2.; std_dev = (upr - lwr)/sqrt(12.); }


inline void Transformation::
moments_from_loguniform_params(const Real& lwr, const Real& upr, Real& mean,
			       Real& std_dev)
{
  Real range = upr - lwr, log_range = log(upr) - log(lwr);
  mean       = range/log_range;
  std_dev    = sqrt(range*(log_range*(upr+lwr)-2.*range)/2.)/log_range;
}


inline void Transformation::
moments_from_triangular_params(const Real& lwr, const Real& upr,
			       const Real& mode, Real& mean, Real& std_dev)
{
  mean    = (lwr + mode + upr)/3.;
  std_dev = sqrt((lwr*(lwr - mode) + mode*(mode - upr) + upr*(upr - lwr))/18.);
}


inline void Transformation::
moments_from_exponential_params(const Real& beta, Real& mean, Real& std_dev)
{ mean = beta; std_dev = beta; }


inline void Transformation::
moments_from_beta_params(const Real& lwr, const Real& upr, const Real& alpha,
			 const Real& beta, Real& mean, Real& std_dev)
{
  Real range = upr - lwr;
  mean       = lwr + alpha/(alpha+beta)*range;
  std_dev    = sqrt(alpha*beta/(alpha+beta+1.))/(alpha+beta)*range;
}


inline void Transformation::
moments_from_gamma_params(const Real& alpha, const Real& beta, Real& mean,
			  Real& std_dev)
{ mean = alpha*beta; std_dev = sqrt(alpha)*beta; }


inline void Transformation::
moments_from_gumbel_params(const Real& alpha, const Real& beta, Real& mean,
			   Real& std_dev)
{ mean = beta + 0.5772/alpha; std_dev = Pi/sqrt(6.)/alpha; }


inline void Transformation::
moments_from_frechet_params(const Real& alpha, const Real& beta, Real& mean,
			    Real& std_dev)
{
#ifdef HAVE_GSL
  // See Haldar and Mahadevan, pp. 91-92
  Real gam = gsl_sf_gamma(1.-1./alpha);
  mean     = beta*gam;
  std_dev  = beta*sqrt(gsl_sf_gamma(1.-2./alpha)-gam*gam);
#else
  Cerr << "Error: frechet distributions only suported in executables "
       << "configured with the GSL library." << endl;
  abort_handler(-1);
#endif // HAVE_GSL
}


inline void Transformation::
moments_from_weibull_params(const Real& alpha, const Real& beta, Real& mean,
			    Real& std_dev)
{
#ifdef HAVE_GSL
  // See Haldar and Mahadevan, p. 97
  Real gam = gsl_sf_gamma(1.+1./alpha),
    cf_var = sqrt(gsl_sf_gamma(1.+2./alpha)/gam/gam - 1.);
  mean     = beta*gam;
  std_dev  = cf_var*beta*gam;
#else
  Cerr << "Error: weibull distributions only suported in executables "
       << "configured with the GSL library." << endl;
  abort_handler(-1);
#endif // HAVE_GSL
}

} // namespace Pecos

#endif
