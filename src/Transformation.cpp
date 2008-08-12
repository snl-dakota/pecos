/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 Transformation
//- Description: Base class for nonlinear distribution transformations
//- Owner:       Mike Eldred
//- Checked by:
//- Version:

#include "pecos_global_defs.h"
#include "NatafTransformation.hpp"
#ifdef PECOS_GSL
#include "gsl/gsl_sf_gamma.h"
#endif

const Pecos::Real Pecos::Transformation::Pi = 3.1415926535897932385;

static const char rcsId[]="@(#) $Id: Transformation.C 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_trans() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_trans() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~Transformation). */
Transformation::Transformation(BaseConstructor):
  extendedUSpace(false), correlationFlagX(false),
  transRep(NULL), referenceCount(1)
{
  const RealMatrix& uncertain_corr = model.uncertain_correlations();
  if (!uncertain_corr.empty()) {
    // set correlationFlagX
    for (size_t i=1; i<numUncertainVars; i++)
      for (size_t j=0; j<i; j++)
	if (fabs(uncertain_corr[i][j]) > 1.e-25)
	  correlationFlagX = true;
    // Copy the correlation matrix to a Teuchos array
    if (correlationFlagX)
      copy_data(uncertain_corr, corrMatrixX);
  }

#ifdef REFCOUNT_DEBUG
  Cout << "Transformation::Transformation(BaseConstructor) called to build "
       << "base class for letter." << std::endl;
#endif
}


/** The default constructor: transRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
Transformation::Transformation(): transRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  Cout << "Transformation::Transformation() called to build empty "
       << "transformation object." << std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_trans, since Transformation(BaseConstructor) builds
    the actual base class data for the derived transformations. */
Transformation::Transformation(const std::string& trans_type):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  Cout << "Transformation::Transformation(ProblemDescDB&) called to "
       << "instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  transRep = get_trans(trans_type);
  if ( !transRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize transRep to the 
    appropriate derived type. */
Transformation* Transformation::get_trans(const std::string& trans_type)
{
#ifdef REFCOUNT_DEBUG
  Cout << "Envelope instantiating letter in get_trans(string&)." << std::endl;
#endif

  if (trans_type == "nataf")
    return new NatafTransformation();
  else {
    Cerr << "Error: Transformation type " << trans_type << " not available."
	 << std::endl;
    return NULL;
  }
}


/** Copy constructor manages sharing of transRep and incrementing
    of referenceCount. */
Transformation::Transformation(const Transformation& trans)
{
  // Increment new (no old to decrement)
  transRep = trans.transRep;
  if (transRep) // Check for an assignment of NULL
    transRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  Cout << "Transformation::Transformation(Transformation&)" << std::endl;
  if (transRep)
    Cout << "transRep referenceCount = " << transRep->referenceCount << std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old transRep, assigns
    new transRep, and increments referenceCount for new transRep. */
Transformation Transformation::operator=(const Transformation& trans)
{
  // Decrement old
  if (transRep) // Check for null pointer
    if (--transRep->referenceCount == 0) 
      delete transRep;
  // Increment new
  transRep = trans.transRep;
  if (transRep) // Check for an assignment of NULL
    transRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  Cout << "Transformation::operator=(Transformation&)" << std::endl;
  if (transRep)
    Cout << "transRep referenceCount = " << transRep->referenceCount << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes transRep
    when referenceCount reaches zero. */
Transformation::~Transformation()
{ 
  // Check for NULL pointer 
  if (transRep) {
    --transRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    Cout << "transRep referenceCount decremented to " 
	 << transRep->referenceCount << std::endl;
#endif
    if (transRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      Cout << "deleting transRep" << std::endl;
#endif
      delete transRep;
    }
  }
}


void Transformation::
reshape_correlation_matrix(size_t num_design_vars, size_t num_uncertain_vars,
			   size_t num_state_vars)
{
  size_t i, j, offset, num_corr_vars = corrMatrixX.M(),
    num_active_vars = num_design_vars + num_uncertain_vars + num_state_vars;
  if (num_corr_vars != num_active_vars) {
    if (num_corr_vars != num_uncertain_vars) {
      Cerr << "\nError: unknown symmetric matrix dimension (" << num_corr_vars
	   << ") in Transformation::reshape_correlation_matrix()." << std::endl;
      abort_handler(-1);
    }
    RealSymMatrix old_corr_matrix(corrMatrixX);
    corrMatrixX.Shape(num_active_vars); // initializes to zero
    for (i=0; i<num_design_vars; i++)
      corrMatrixX(i,i) = 1.;
    offset = num_design_vars;
    for (i=0; i<num_uncertain_vars; i++)
      for (j=0; j<num_uncertain_vars; j++)
	corrMatrixX(i+offset,j+offset) = old_corr_matrix(i,j);
    offset += num_uncertain_vars;
    for (i=0; i<num_state_vars; i++)
      corrMatrixX(i+offset,i+offset) = 1.;
  }
}


/** This function is commonly used to publish tranformation data when
    the Model variables are in a transformed space (e.g., u-space) and
    ranVarTypes et al. may not be generated directly.  This allows for
    the use of inverse transformations to return the transformed space
    variables to their original states. */
void Transformation::
initialize_random_variables(const ShortArray& x_types,
			    const ShortArray& u_types,
			    const RealVector& x_means,
			    const RealVector& x_std_devs,
			    const RealVector& x_l_bnds,
			    const RealVector& x_u_bnds,
			    const RealVectorArray& x_addtl,
			    const RealSymMatrix& x_corr,
			    const RealMatrix& z_chol_fact,
			    bool x_corr_flag, bool ext_u_space)
{
  initialize_random_variables_types(x_types, u_types);
  initialize_random_variables_parameters(x_means, x_std_devs, x_l_bnds,
					 x_u_bnds, x_addtl, x_corr, z_chol_fact,
					 x_corr_flag, ext_u_space);
}


void Transformation::
initialize_random_variables_types(const ShortArray& x_types,
				  const ShortArray& u_types)
{
  ranVarTypesX  = x_types;
  ranVarTypesU  = u_types;
}


void Transformation::
initialize_random_variables_parameters(const RealVector& x_means,
				       const RealVector& x_std_devs,
				       const RealVector& x_l_bnds,
				       const RealVector& x_u_bnds,
				       const RealVectorArray& x_addtl,
				       const RealSymMatrix& x_corr,
				       const RealMatrix& z_chol_fact,
				       bool x_corr_flag, bool ext_u_space)
{
  ranVarMeansX        = x_means;
  ranVarStdDevsX      = x_std_devs;
  ranVarLowerBndsX    = x_l_bnds;
  ranVarUpperBndsX    = x_u_bnds;
  ranVarAddtlParamsX  = x_addtl;
  corrMatrixX         = x_corr;
  corrCholeskyFactorZ = z_chol_fact;
  correlationFlagX    = x_corr_flag;
  extendedUSpace      = ext_u_space;
}


/** This procedure computes numerical derivatives of x and/or z with respect to
    distribution parameters s, and is used by jacobian_dX_dS() to provide data
    that is not available analytically.  Numerical dz/ds involves dL/ds
    (z(s) = L(s) u and dz/ds = dL/ds u) and is needed to evaluate dx/ds
    semi-analytically for correlated variables.  Numerical dx/ds is needed for
    distributions lacking simple closed-form CDF expressions (beta and gamma
    distributions). */
void Transformation::
numerical_design_jacobian(const RealVector& x_vars,
			  bool xs, RealMatrix& num_jacobian_xs,
			  bool zs, RealMatrix& num_jacobian_zs,
			  const UIntArray& cv_ids, const UIntArray& acv_ids)
{
  // For correlated vars, correlation matrix C = C(s) due to Nataf modifications
  //   z(s) = L(s) u  ->  dz/ds = dL/ds u  ->  need dL/ds
  //   C(s) = L(s) L(s)^T
  //   dC/ds (which could be derived analytically) = L dL/ds^T + dL/ds L^T
  // This takes the form dC/ds = A + A^T where A = L dL/ds^T
  // Unfortunately, solution of this equation for general A (which could
  // provide dL/ds) given symmetric dC/ds is not possible since it is nonunique.
  // Since we will not be differentiating the Cholesky solver, we will use
  // semi-analytic design sensitivities with numerical dL/ds.  Note that
  // numerical dz/ds is simpler and would likely be just as effective, but in
  // general, semi-analytic sensitivities should minimize the numerical portion.

  // Rectangular Jacobians = Gradient^T = num_Z x num_S where num_S is the total
  // number of active continuous vars flowed down from a higher iteration level.
  // The number of distribution parameter insertions is <= num_S.
  size_t i, j, k, num_var_map_1c = primaryACVarMapIndices.length();
  int x_len = x_vars.Length();
  if (xs && (num_jacobian_xs.M() != x_len ||
	     num_jacobian_xs.N() != num_var_map_1c) )
    num_jacobian_xs.Shape(x_len, num_var_map_1c);
  if (zs && (num_jacobian_zs.M() != x_len ||
	     num_jacobian_zs.N() != num_var_map_1c) )
    num_jacobian_zs.Shape(x_len, num_var_map_1c);

  RealMatrix L_s_plus_h, dL_dsi;
  RealVector dz_dsi;
  //RealVector z_vars_s_plus_h, z_vars_s_minus_h;
  RealVector x_vars_s_plus_h, x_vars_s_minus_h;
  if (zs) {
    L_s_plus_h.Shape(x_len, x_len);
    dL_dsi.Shape(x_len, x_len);
    dz_dsi.Size(x_len);
  }

  RealVector u_vars;
  trans_X_to_U(x_vars, u_vars);

  Real fd_grad_ss = 1.e-4;
  for (i=0; i<num_var_map_1c; i++) {

    size_t cv_index = cv_ids.index(acv_ids[primaryACVarMapIndices[i]]);
    if (cv_index != _NPOS && secondaryACVarMapTargets[i] != NO_TARGET) {

      Real s0 = distribution_parameter(i);

      // Compute the offset for the ith gradient variable.
      // Enforce a minimum delta of fdgss*.01
      Real h_mag = fd_grad_ss * max(fabs(s0), .01);
      Real h = (s0 < 0.0) ? -h_mag : h_mag; // h has same sign as s0

      // -----------------------------------
      // Evaluate (L/z_vars/x_vars)_s_plus_h
      // -----------------------------------
      Real s1 = s0 + h;
      distribution_parameter(i, s1); // updates ranVars & corrCholeskyFactorZ
      if (zs) {
	L_s_plus_h = corrCholeskyFactorZ;        // L
	//trans_U_to_Z(u_vars, z_vars_s_plus_h); // z
      }
      if (xs)
	trans_U_to_X(u_vars, x_vars_s_plus_h);   // x

      // ------------------------------------
      // Evaluate (L/z_vars/x_vars)_s_minus_h
      // ------------------------------------
      s1 = s0 - h;
      distribution_parameter(i, s1); // updates ranVars & corrCholeskyFactorZ
      //if (zs) {
        // utilize corrCholeskyFactorZ below      // L
        //trans_U_to_Z(u_vars, z_vars_s_minus_h); // z
      //}
      if (xs)
	trans_U_to_X(u_vars, x_vars_s_minus_h);   // x

      // -------------------------------
      // Compute the central differences
      // -------------------------------
      if (zs) {
	for (j=0; j<x_len; j++)                            // dL/ds
	  for (k=0; k<=j; k++)
	    dL_dsi(j, k) = (L_s_plus_h(j, k) - corrCholeskyFactorZ(j, k))/2./h;
	dz_dsi.Multiply('N', 'N', 1., dL_dsi, u_vars, 0.); // dz/ds
	for (j=0; j<x_len; j++)
	  num_jacobian_zs(j, i) = dz_dsi(j);
	//for (j=0; j<x_len; j++)                          // dz/ds (alt)
	//  num_jacobian_zs(j, i)=(z_vars_s_plus_h(j)-z_vars_s_minus_h(j))/2./h;
      }
      if (xs)
	for (j=0; j<x_len; j++)                            // dx/ds
	  num_jacobian_xs(j,i) = (x_vars_s_plus_h(j)-x_vars_s_minus_h(j))/2./h;

      distribution_parameter(i, s0); // resets s0 (& corrCholeskyFactorZ!)
    }
  }
}


#ifndef PECOS_GSL
Real Transformation::erf_inverse(const Real& p)
{
  // Adapted from ltqnorm.m see URL below for more info

  // The algorithm uses a minimax approximation by rational functions and
  // the result has a relative error less than 1.15e-9.  A last refinement
  // by Halley's rational method is applied to achieve full machine precision.

  // Author:      Peter J. Acklam
  // Time-stamp:  2000-07-19 16:44:07
  // E-mail:      pjacklam@online.no
  // URL:         http://home.online.no/~pjacklam/notes/invnorm/

  Real a[7] = { 0.,         -3.969683028665376e+01,  2.209460984245205e+02,
    -2.759285104469687e+02,  1.383577518672690e+02, -3.066479806614716e+01,
     2.506628277459239e+00 };
  Real b[6] = { 0.,         -5.447609879822406e+01,  1.615858368580409e+02,
    -1.556989798598866e+02,  6.680131188771972e+01, -1.328068155288572e+01 };
  Real c[7] = { 0.,         -7.784894002430293e-03, -3.223964580411365e-01,
    -2.400758277161838e+00, -2.549732539343734e+00,  4.374664141464968e+00,
     2.938163982698783e+00 };
  Real d[5] = { 0.,          7.784695709041462e-03,  3.224671290700398e-01,
     2.445134137142996e+00,  3.754408661907416e+00 };

  // Modify p to since this is the error function
  // 1/sqrt(2.)*ltqnorm((p+1)/2)=erf_inverse(p)
  // redefine p here
  Real p_new = 0.5*(p + 1.);
  // scale at the end

  // Define break-points.
  Real plow  = 0.02425;
  Real phigh = 1. - plow;

  // Rational approximation for central region
  Real r, q, z;
  if (p_new > plow && p_new < phigh) { // -plow <= p_new <= phigh 
    q = p_new - 0.5;
    r = q*q;
    z = (((((a[1]*r+a[2])*r+a[3])*r+a[4])*r+a[5])*r+a[6])*q/ 
        (((((b[1]*r+b[2])*r+b[3])*r+b[4])*r+b[5])*r+1);
  }
  // Rational approximation for lower region
  else if (p_new < plow && p_new > 0.) { // plow > p_new > 0.
    q = sqrt(-2*log(p_new));
    z = (((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6])/
        ((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1);  
  }
  else if (p_new > phigh && p_new < 1.) { // phigh < p_new < 1.
    q = sqrt(-2*(log(1-p_new)));
    z = -(((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6])/
         ((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1); 
  }
  else if (p_new == 0.)
    z = -1.*pow(10., 150.);
  else if (p_new == 1.)
    z = pow(10., 150.); 
  else if (p_new < 0. || p_new > 1.) {
    Cerr << "Error: probability greater than 1 or less than 0 in erf_inverse()."
         << std::endl;
    abort_handler(-1);
  }
  // user erf instead of erfc

  Real e = 0.5*(1. - erf(-z/sqrt(2.))) - p_new;
  Real u = e*sqrt(2.*Pi)*exp(z*z/2.);
  z = z - u/(1. + z*u/2.);
  // scale since this is the erf inverse not gaussian inverse
  // see above
  z = 1/sqrt(2.)*z;
  return z;
}
#endif // PECOS_GSL


#ifdef PECOS_GSL
/** Solve is performed in scaled space (for the standard beta distribution). */
Real Transformation::
cdf_beta_Pinv(const Real& normcdf, const Real& alpha, const Real& beta)
{
  // F(x) = Phi(z) = normcdf
  // F(x) - normcdf = 0
  // x^{i+1} = x - (F(x) - normcdf)/f(x)

  // Initial guess: model as uniform
  // (linear CDF accumulation between [0,1] bounds)
  Real scaled_x = normcdf;

  // evaluate residual F(x) - normcdf
  Real res = gsl_cdf_beta_P(scaled_x, alpha, beta) - normcdf;

  size_t newton_iters = 0, max_iters = 20;
  bool converged = false, terminate = false;
  while (!terminate && !converged) {

    // evaluate residual derivative f(x)
    Real dres_dx = gsl_ran_beta_pdf(scaled_x, alpha, beta);

    // compute Newton step
    Real delta_scaled_x;
    if (fabs(dres_dx) > DBL_MIN) {
      delta_scaled_x = -res/dres_dx; // full Newton step
      if (fabs(delta_scaled_x) < convergenceTol)
	converged = true; // but go ahead and take the step, if beneficial
    }
    else
      terminate = true;

    // Simple backtracking line search globalization
    bool reduction = false;
    size_t backtrack_iters = 0;
    while (!reduction && !terminate) {
      Real scaled_x_step = scaled_x + delta_scaled_x;
      // scaled_x must lie in [0,1]
      if (scaled_x_step < 0.) scaled_x_step = 0.;
      if (scaled_x_step > 1.) scaled_x_step = 1.;

      // evaluate residual at scaled_x_step
      Real res_step = gsl_cdf_beta_P(scaled_x_step, alpha, beta) - normcdf;

      // perform backtracking line search to enforce decrease in res
      if ( fabs(res_step) < fabs(res) ) { // accept step
	reduction = true;
	scaled_x = scaled_x_step;
	res      = res_step;
	//Cout << "residual = " << res << " delta = " << delta_scaled_x
	//     << " scaled_x = " << scaled_x << '\n';
      }
      else if (converged)
	terminate = true; // kick out of inner while
      else {
	//Cout << "Backtracking\n";
	delta_scaled_x /= 2.; // backtrack
	if (backtrack_iters++ >= max_iters) // backtrack iter must complete
	  terminate = true;
      }
    }
    if (++newton_iters >= max_iters) // Newton iteration has completed
      terminate = true;
  }
  return scaled_x;
}
#endif // PECOS_GSL


#ifdef DERIV_DEBUG
void Transformation::verify_trans_jacobian_hessian(const RealVector& v0)
{
  size_t i, j, k;
  bool fd_grad_flag = true, fd_hess_flag = true, fd_hess_by_fn_flag = false,
       fd_hess_by_grad_flag = true;

  Real fd_grad_ss = 1.e-8, fd_hess_by_fn_ss = 2.e-8, fd_hess_by_grad_ss = 1.e-8;

  RealVector trans_vars_v0;
  //trans_X_to_U(v0, trans_vars_v0); // v = x, trans_vars_v = u
  trans_U_to_X(v0, trans_vars_v0); // v = u, trans_vars_v = x
  int num_v = v0.Length(), num_tv = trans_vars_v0.Length();

  RealMatrix num_jac_dtv_dv(num_tv, num_v);
  RealSymMatrixArray num_hess_d2tv_dv2(num_tv);
  for (i=0; i<num_tv; i++)
    num_hess_d2tv_dv2[i].Shape(num_v);

  // ------------------------------
  // Estimate numerical derivatives
  // ------------------------------
  if (fd_grad_flag || fd_hess_flag) {
    RealVector v1 = v0; // for perturbed values
    // ---------------
    // Loop over num_v
    // ---------------
    for (j=0; j<num_v; j++) { // difference the 1st num_v vars
      if (fd_grad_flag) {

	// Compute the offset for the ith gradient variable.
	// Enforce a minimum delta of fdgss*.01
	Real h_mag = fd_grad_ss * max(fabs(v0(j)), .01);
	Real h = (v1(j) < 0.0) ? -h_mag : h_mag; // h has same sign as v1(j)

	// ----------------------------
	// Evaluate trans_vars_v_plus_h
	// ----------------------------
	RealVector trans_vars_v_plus_h;
	v1(j) = v0(j) + h;
	Cout << ">>>>> Pecos finite difference gradient evaluation for v["
	     << j+1 << "] + h:\n";
	//trans_X_to_U(v1, trans_vars_v_plus_h);
	trans_U_to_X(v1, trans_vars_v_plus_h);

	// -----------------------------
	// Evaluate trans_vars_v_minus_h
	// -----------------------------
	RealVector trans_vars_v_minus_h;
	v1(j) = v0(j) - h;
	Cout << ">>>>> Pecos finite difference gradient evaluation for v["
	     << j+1 << "] - h:\n";
	//trans_X_to_U(v1, trans_vars_v_minus_h);
	trans_U_to_X(v1, trans_vars_v_minus_h);

	// always use central diffs for verification purposes
	for (i=0; i<num_tv; i++)
	  num_jac_dtv_dv(i,j)
	    = (trans_vars_v_plus_h(i) - trans_vars_v_minus_h(i))/2./h;
      }

      if (fd_hess_flag) {

	if (fd_hess_by_fn_flag) {
	  RealVector trans_vars_v_plus_2h, trans_vars_v_minus_2h;

	  // Compute the 2nd-order Hessian offset for the ith variable.
	  // Enforce a minimum delta of fdhss*.01
	  Real h_mag = fd_hess_by_fn_ss * max(fabs(v0(j)), .01);
	  Real h = (v1(j) < 0.0) ? -h_mag : h_mag; // h has same sign as v1(j)

	  // evaluate diagonal term

	  // -----------------------------
	  // Evaluate trans_vars_v_plus_2h
	  // -----------------------------
	  v1(j) = v0(j) + 2.*h;
	  Cout << ">>>>> Pecos finite difference Hessian evaluation for v["
	       << j+1 << "] + 2h:\n";
	  //trans_X_to_U(v1, trans_vars_v_plus_2h);
	  trans_U_to_X(v1, trans_vars_v_plus_2h);

	  // ------------------------------
	  // Evaluate trans_vars_v_minus_2h
	  // ------------------------------
	  v1(j) = v0(j) - 2.*h;
	  Cout << ">>>>> Pecos finite difference Hessian evaluation for v["
	       << j+1 << "] - 2h:\n";
	  //trans_X_to_U(v1, trans_vars_v_minus_2h);
	  trans_U_to_X(v1, trans_vars_v_minus_2h);

	  for (i=0; i<num_tv; i++)
	    num_hess_d2tv_dv2[i](j,j)
	      = (trans_vars_v_plus_2h(i) - 2.*trans_vars_v0(i) +
		 trans_vars_v_minus_2h(i))/(4.*h*h);

	  // evaluate off-diagonal terms

	  for (k=j+1; k<num_v; k++) {
	    RealVector trans_vars_v_plus_h_plus_h,
	      trans_vars_v_plus_h_minus_h, trans_vars_v_minus_h_plus_h,
	      trans_vars_v_minus_h_minus_h;

	    // -----------------------------------
	    // Evaluate trans_vars_v_plus_h_plus_h
	    // -----------------------------------
	    v1(j) = v0(j) + h;
	    v1(k) = v0(k) + h;
	    Cout << ">>>>> Pecos finite difference Hessian evaluation for v["
		 << j+1 << "] + h, v[" << k+1 << "] + h:\n";
	    //trans_X_to_U(v1, trans_vars_v_plus_h_plus_h);
	    trans_U_to_X(v1, trans_vars_v_plus_h_plus_h);
	    // ------------------------------------
	    // Evaluate trans_vars_v_plus_h_minus_h
	    // ------------------------------------
	    //v1(j) = v0(j) + h;
	    v1(k) = v0(k) - h;
	    Cout << ">>>>> Pecos finite difference Hessian evaluation for v["
		 << j+1 << "] + h, v[" << k+1 << "] - h:\n";
	    //trans_X_to_U(v1, trans_vars_v_plus_h_minus_h);
	    trans_U_to_X(v1, trans_vars_v_plus_h_minus_h);
	    // ------------------------------------
	    // Evaluate trans_vars_v_minus_h_plus_h
	    // ------------------------------------
	    v1(j) = v0(j) - h;
	    v1(k) = v0(k) + h;
	    Cout << ">>>>> Pecos finite difference Hessian evaluation for v["
		 << j+1 << "] - h, v[" << k+1 << "] + h:\n";
	    //trans_X_to_U(v1, trans_vars_v_minus_h_plus_h);
	    trans_U_to_X(v1, trans_vars_v_minus_h_plus_h);
	    // -------------------------------------
	    // Evaluate trans_vars_v_minus_h_minus_h
	    // -------------------------------------
	    //v1(j) = v0(j) - h;
	    v1(k) = v0(k) - h;
	    Cout << ">>>>> Pecos finite difference Hessian evaluation for v["
		 << j+1 << "] - h, v[" << k+1 << "] - h:\n";
	    //trans_X_to_U(v1, trans_vars_v_minus_h_minus_h);
	    trans_U_to_X(v1, trans_vars_v_minus_h_minus_h);

	    for (i=0; i<num_tv; i++)
	      num_hess_d2tv_dv2[i](j,k) = num_hess_d2tv_dv2[i](k,j)
		= (trans_vars_v_plus_h_plus_h(i)
		- trans_vars_v_plus_h_minus_h(i)
		- trans_vars_v_minus_h_plus_h(i)
		+ trans_vars_v_minus_h_minus_h(i)) / (4.*h*h);

	    v1(k) = v0(k);
	  }
	}

	if (fd_hess_by_grad_flag) {

	  // Compute the 1st-order Hessian offset for the ith variable.
	  // Enforce a minimum delta of fdhss*.01
	  Real h_mag = fd_hess_by_grad_ss * max(fabs(v0(j)), .01);
	  Real h = (v1(j) < 0.0) ? -h_mag : h_mag; // h has same sign as v1(j)

	  // --------------------------
	  // Evaluate fn_grads_v_plus_h
	  // --------------------------
	  v1(j) = v0(j) + h;
	  Cout << ">>>>> Pecos finite difference Hessian evaluation for v["
	       << j+1 << "] + h:\n";
	  RealVector trans_vars_v_plus_h;
	  trans_U_to_X(v1, trans_vars_v_plus_h);
	  RealMatrix jac_v0, jac_v_plus_h;
	  // jacobian routines use x_vars:
	  jacobian_dX_dU(trans_vars_v0,       jac_v0);
	  jacobian_dX_dU(trans_vars_v_plus_h, jac_v_plus_h);
	  for (i=0; i<num_tv; i++)
	    for (k=0; k<num_v; k++)
	      num_hess_d2tv_dv2[i](j,k)	= (jac_v_plus_h(i,k) - jac_v0(i,k))/h;
	}
      }
      v1(j) = v0(j);
    }
  }

  // Enforce symmetry in the case of FD Hessians from 1st-order gradient
  // differences by averaging off-diagonal terms: H' = 1/2 (H + H^T)
  if (fd_hess_by_grad_flag)
    for (i=0; i<num_tv; i++)
      for (j=0; j<num_v; j++)
	for (k=j+1; k<num_v; k++)
	  num_hess_d2tv_dv2[i](j,k) = num_hess_d2tv_dv2[i](k,j)
	    = (num_hess_d2tv_dv2[i](j,k) + num_hess_d2tv_dv2[i](k,j))/2.;

  // Print out numerical and analytic:
  RealVector x0(num_tv);
  //x0 = v0;
  //RealMatrix jacobian_ux;
  //jacobian_dU_dX(x0, jacobian_ux);
  trans_U_to_X(v0, x0);
  RealMatrix jacobian_xu;
  jacobian_dX_dU(x0, jacobian_xu);
  Cout << "\nNumerical jacobian:" << num_jac_dtv_dv
       << "\nAnalytic jacobian:"  << jacobian_xu; //jacobian_ux;
  RealSymMatrixArray hessian_xu(num_tv);
  hessian_d2X_dU2(x0, hessian_xu);
  for (i=0; i<num_tv; i++)
    Cout << "\nNumerical Hessian:" << num_hess_d2tv_dv2[i]
	 << "\nAnalytic Hessian:"  << hessian_xu[i];
}


void Transformation::verify_design_jacobian(const RealVector& u0)
{
  RealVector x0;
  trans_U_to_X(u0, x0);

  RealMatrix num_jac_dx_ds, num_jac_dz_ds;
  numerical_design_jacobian(x0, true, num_jac_dx_ds, false, num_jac_dz_ds);

  RealMatrix jacobian_xs;
  jacobian_dX_dS(x0, jacobian_xs);

  // Print out numerical and analytic:
  Cout << "\nNumerical jacobian:" << num_jac_dx_ds
       << "\nAnalytic jacobian:"  << jacobian_xs;
}
#endif // DERIV_DEBUG

} // namespace Pecos
