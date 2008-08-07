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

#include "pecos_data_types.h"
#ifdef PECOS_GSL
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#endif


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
  //- Heading: Member functions
  //

  /// initialize ranVarTypesX/U, ranVarMeansX, ranVarStdDevsX, ranVarLowerBndsX,
  /// ranVarUpperBndsX, and ranVarAddtlParamsX based on distribution data from
  /// iteratedModel (corrCholeskyFactorZ and correlationFlagX set separately)
  void initialize_random_variables();
  /// set ranVarTypesX/U, ranVarMeansX, ranVarStdDevsX, ranVarLowerBndsX,
  /// ranVarUpperBndsX, ranVarAddtlParamsX, corrCholeskyFactorZ,
  /// correlationFlagX, and extendedUSpace based on incoming data
  void initialize_random_variables(const ShortArray& x_types,
				   const ShortArray& u_types,
				   const RealVector& x_means,
				   const RealVector& x_std_devs,
				   const RealVector& x_l_bnds,
				   const RealVector& x_u_bnds,
				   const RealVectorArray& x_addtl,
				   const RealSymMatrix& x_corr,
				   const RealMatrix& z_chol_fact,
				   bool x_corr_flag, bool ext_u_space);

  /// set distParamDerivs
  void distribution_parameter_derivatives(bool dist_param_derivs);

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

  /// initializes ranVarTypesX and ranVarTypesU
  void initialize_random_variable_types();
  /// initializes ranVarMeansX, ranVarStdDevsX, ranVarLowerBndsX,
  /// ranVarUpperBndsX, and ranVarAddtlParamsX
  void initialize_random_variable_parameters();

  /// reshape corrMatrixX for an all_variables specification
  void reshape_correlation_matrix();

  /// Standard normal density function
  Real phi(const Real& beta);

  /// Standard normal cumulative distribution function
  Real Phi(const Real& beta);

  /// Inverse of standard normal cumulative distribution function
  Real Phi_inverse(const Real& p);

#ifdef DERIV_DEBUG
  /// routine for verification of transformation Jacobian/Hessian terms
  void verify_trans_jacobian_hessian(const RealVector& v0);

  /// routine for verification of design Jacobian terms
  void verify_design_jacobian(const RealVector& u0);
#endif // DERIV_DEBUG

#ifdef PECOS_GSL
  /// Inverse of standard beta CDF (not supported by GSL)
  Real cdf_beta_Pinv(const Real& normcdf, const Real& alpha, const Real& beta);
#endif // PECOS_GSL

  //
  //- Heading: Data members
  //

  // the following attributes are the required data for performing
  // transformations from X -> Z -> U and back.

  /// Flag for extended u-space variable definitions.  If false
  /// (reliability, AIS), then u-space is defined exclusively with std
  /// normals.  If true (PCE), then u-space is defined by std normals,
  /// std uniforms, std exponentials, std betas, and std gammas.
  bool extendedUSpace;
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
  /// Epetra copy of Model::uncertainCorrelations
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

  /// Used only by the standard envelope constructor to initialize
  /// transRep to the appropriate derived type.
  Transformation* get_trans(const String& trans_type);

  /// As part of the Nataf distribution model (Der Kiureghian & Liu, 1986),
  /// this procedure modifies the user-specified correlation matrix
  /// (corrMatrixX) to account for correlation warping from the nonlinear
  /// X->Z transformation and performs a Cholesky factorization to create
  /// corrCholeskyFactorZ.
  void trans_correlations();

  /// Computes numerical dx/ds and dz/ds Jacobians as requested by xs
  /// and zs booleans
  void numerical_design_jacobian(const RealVector& x_vars,
                                 bool xs, RealMatrix& num_jacobian_xs,
                                 bool zs, RealMatrix& num_jacobian_zs);

#ifndef PECOS_GSL
  /// Inverse of error function used in Phi_inverse()
  Real erf_inverse(const Real& p);
#endif // PECOS_GSL

  /// return a particular random variable distribution parameter
  const Real& distribution_parameter(const size_t& outer_cv_index);
  /// set a particular random variable distribution parameter and
  /// update derived quantities
  void distribution_parameter(const size_t& outer_cv_index, const Real& param);

  //
  //- Heading: Data members
  //

  /// flags calculation of derivatives with respect to distribution
  /// parameters s within resp_x_to_u_mapping() using the chain rule
  /// df/dx dx/ds.  The default is to calculate derivatives with respect
  /// to standard random variables u using the chain rule df/dx dx/du.
  bool distParamDerivs;

  /// pointer to the letter (initialized only for the envelope)
  Transformation* transRep;
  /// number of objects sharing transRep
  int referenceCount;
};


inline void Transformation::
distribution_parameter_derivatives(bool dist_param_derivs)
{ distParamDerivs = dist_param_derivs; }


inline Real Transformation::phi(const Real& beta)
{
#ifdef PECOS_GSL
  return gsl_ran_ugaussian_pdf(beta);
#else
  return exp(-beta*beta/2.)/sqrt(2.*Pi);
#endif // PECOS_GSL
}


/** returns a probability < 0.5 for negative beta and a probability > 0.5
    for positive beta. */
inline Real Transformation::Phi(const Real& beta)
{
#ifdef PECOS_GSL
  return gsl_cdf_ugaussian_P(beta);
#else
  return .5 + .5*erf(beta/sqrt(2.));
#endif // PECOS_GSL
}


/** returns a negative beta for probability < 0.5 and a positive beta for
    probability > 0.5. */
inline Real Transformation::Phi_inverse(const Real& p)
{
#ifdef PECOS_GSL
  return gsl_cdf_ugaussian_Pinv(p);
#else
  return sqrt(2.)*erf_inverse(2.*p - 1.);
#endif // PECOS_GSL
}

} // namespace Pecos

#endif
