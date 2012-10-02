/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "pecos_stat_util.hpp"
#include "NatafTransformation.hpp"

static const char rcsId[]="@(#) $Id: ProbabilityTransformation.cpp 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_prob_trans() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_prob_trans() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~ProbabilityTransformation). */
ProbabilityTransformation::ProbabilityTransformation(BaseConstructor):
  correlationFlagX(false), probTransRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "ProbabilityTransformation::ProbabilityTransformation(Base"
        << "Constructor) called to build base class for letter." << std::endl;
#endif
}


/** The default constructor: probTransRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
ProbabilityTransformation::ProbabilityTransformation():
  probTransRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "ProbabilityTransformation::ProbabilityTransformation() called to "
        << "build empty envelope." << std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_prob_trans, since ProbabilityTransformation(BaseConstructor)
    builds the actual base class data for the derived transformations. */
ProbabilityTransformation::
ProbabilityTransformation(const String& prob_trans_type):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "ProbabilityTransformation::ProbabilityTransformation(string&) "
        << "called to instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  probTransRep = get_prob_trans(prob_trans_type);
  if ( !probTransRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize probTransRep to the 
    appropriate derived type. */
ProbabilityTransformation* ProbabilityTransformation::
get_prob_trans(const String& prob_trans_type)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_prob_trans(string&)."
        << std::endl;
#endif

  if (prob_trans_type == "nataf")
    return new NatafTransformation();
  else {
    PCerr << "Error: ProbabilityTransformation type " << prob_trans_type
	  << " not available." << std::endl;
    return NULL;
  }
}


/** Copy constructor manages sharing of probTransRep and incrementing
    of referenceCount. */
ProbabilityTransformation::
ProbabilityTransformation(const ProbabilityTransformation& prob_trans)
{
  // Increment new (no old to decrement)
  probTransRep = prob_trans.probTransRep;
  if (probTransRep) // Check for an assignment of NULL
    probTransRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  PCout << "ProbabilityTransformation::ProbabilityTransformation("
        << "ProbabilityTransformation&)" << std::endl;
  if (probTransRep)
    PCout << "probTransRep referenceCount = " << probTransRep->referenceCount
	  << std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old probTransRep, assigns
    new probTransRep, and increments referenceCount for new probTransRep. */
ProbabilityTransformation ProbabilityTransformation::
operator=(const ProbabilityTransformation& prob_trans)
{
  if (probTransRep != prob_trans.probTransRep) { // normal case: old != new
    // Decrement old
    if (probTransRep) // Check for null pointer
      if (--probTransRep->referenceCount == 0) 
	delete probTransRep;
    // Assign and increment new
    probTransRep = prob_trans.probTransRep;
    if (probTransRep) // Check for an assignment of NULL
      probTransRep->referenceCount++;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "ProbabilityTransformation::operator=(ProbabilityTransformation&)"
        << std::endl;
  if (probTransRep)
    PCout << "probTransRep referenceCount = " << probTransRep->referenceCount
	  << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes probTransRep
    when referenceCount reaches zero. */
ProbabilityTransformation::~ProbabilityTransformation()
{ 
  // Check for NULL pointer 
  if (probTransRep) {
    --probTransRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "probTransRep referenceCount decremented to " 
	  << probTransRep->referenceCount << std::endl;
#endif
    if (probTransRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting probTransRep" << std::endl;
#endif
      delete probTransRep;
    }
  }
}


/** This function is commonly used to publish tranformation data when
    the Model variables are in a transformed space (e.g., u-space) and
    ranVarTypes et al. may not be generated directly.  This allows for
    the use of inverse transformations to return the transformed space
    variables to their original states. */
void ProbabilityTransformation::
initialize_random_variables(const ProbabilityTransformation& prob_trans)
{
  if (probTransRep) // target is envelope
    probTransRep->initialize_random_variables(prob_trans);
  else {
    if (prob_trans.probTransRep) { // source is envelope
      ranVarTypesX        = prob_trans.probTransRep->ranVarTypesX;
      ranVarTypesU        = prob_trans.probTransRep->ranVarTypesU;
      ranVarMeansX        = prob_trans.probTransRep->ranVarMeansX;
      ranVarStdDevsX      = prob_trans.probTransRep->ranVarStdDevsX;
      ranVarLowerBndsX    = prob_trans.probTransRep->ranVarLowerBndsX;
      ranVarUpperBndsX    = prob_trans.probTransRep->ranVarUpperBndsX;
      ranVarAddtlParamsX  = prob_trans.probTransRep->ranVarAddtlParamsX;
      correlationFlagX    = prob_trans.probTransRep->correlationFlagX;
      corrMatrixX         = prob_trans.probTransRep->corrMatrixX;
      corrCholeskyFactorZ = prob_trans.probTransRep->corrCholeskyFactorZ;
    }
    else { // source is letter
      ranVarTypesX        = prob_trans.ranVarTypesX;
      ranVarTypesU        = prob_trans.ranVarTypesU;
      ranVarMeansX        = prob_trans.ranVarMeansX;
      ranVarStdDevsX      = prob_trans.ranVarStdDevsX;
      ranVarLowerBndsX    = prob_trans.ranVarLowerBndsX;
      ranVarUpperBndsX    = prob_trans.ranVarUpperBndsX;
      ranVarAddtlParamsX  = prob_trans.ranVarAddtlParamsX;
      correlationFlagX    = prob_trans.correlationFlagX;
      corrMatrixX         = prob_trans.corrMatrixX;
      corrCholeskyFactorZ = prob_trans.corrCholeskyFactorZ;
    }
  }
}


void ProbabilityTransformation::
initialize_random_variable_types(const ShortArray& x_types,
				 const ShortArray& u_types)
{
  if (probTransRep) {
    probTransRep->ranVarTypesX  = x_types;
    probTransRep->ranVarTypesU  = u_types;
  }
  else {
    ranVarTypesX  = x_types;
    ranVarTypesU  = u_types;
  }
}


void ProbabilityTransformation::
initialize_random_variable_parameters(const RealVector& x_means,
				      const RealVector& x_std_devs,
				      const RealVector& x_l_bnds,
				      const RealVector& x_u_bnds,
				      const RealVectorArray& x_addtl)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->initialize_random_variable_parameters(x_means,  x_std_devs,
							x_l_bnds, x_u_bnds,
							x_addtl);
  else {
    ranVarMeansX       = x_means;
    ranVarStdDevsX     = x_std_devs;
    ranVarLowerBndsX   = x_l_bnds;
    ranVarUpperBndsX   = x_u_bnds;
    ranVarAddtlParamsX = x_addtl;
  }
}


void ProbabilityTransformation::
initialize_random_variable_correlations(const RealSymMatrix& x_corr)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->initialize_random_variable_correlations(x_corr);
  else {
    corrMatrixX = x_corr;
    size_t num_ran_vars = x_corr.numRows();
    correlationFlagX = false;
    for (size_t i=1; i<num_ran_vars; i++)
      for (size_t j=0; j<i; j++)
	if (std::fabs(x_corr(i,j)) > ZERO_TOL)
	  correlationFlagX = true;
  }
}


void ProbabilityTransformation::
reshape_correlation_matrix(size_t num_leading_vars,
			   size_t num_probabilistic_vars,
			   size_t num_trailing_vars)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->reshape_correlation_matrix(num_leading_vars,
					     num_probabilistic_vars,
					     num_trailing_vars);
  else {
    if (!correlationFlagX)
      return;

    size_t i, j, offset, num_corr_vars = corrMatrixX.numRows(),
      num_active_vars = num_leading_vars + num_probabilistic_vars +
      num_trailing_vars;
    if (num_corr_vars != num_active_vars) {
      if (num_corr_vars != num_probabilistic_vars) {
	PCerr << "\nError: unknown symmetric matrix dim (" << num_corr_vars
	      << ") in ProbabilityTransformation::reshape_correlation_matrix()."
	      << std::endl;
	abort_handler(-1);
      }
      RealSymMatrix old_corr_matrix(corrMatrixX);
      corrMatrixX.shape(num_active_vars); // initializes to zero
      for (i=0; i<num_leading_vars; i++)
	corrMatrixX(i,i) = 1.;
      offset = num_leading_vars;
      for (i=0; i<num_probabilistic_vars; i++)
	for (j=0; j<num_probabilistic_vars; j++)
	  corrMatrixX(i+offset,j+offset) = old_corr_matrix(i,j);
      offset += num_probabilistic_vars;
      for (i=0; i<num_trailing_vars; i++)
	corrMatrixX(i+offset,i+offset) = 1.;
    }
  }
}


void ProbabilityTransformation::
trans_U_to_X(const RealVector& u_vars, RealVector& x_vars)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_U_to_X(u_vars, x_vars);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_U_to_X() virtual fn."
	  << "\nNo default defined at ProbabilityTransformation base class.\n"
	  << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_X_to_U(const RealVector& x_vars, RealVector& u_vars)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_X_to_U(x_vars, u_vars);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_X_to_U() virtual fn."
	  << "\nNo default defined at ProbabilityTransformation base class.\n"
	  << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::transform_correlations()
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->transform_correlations();
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine transform_correlations() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_X_to_U(const RealVector& fn_grad_x, RealVector& fn_grad_u,
		  const RealVector& x_vars,    const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_X_to_U(fn_grad_x, fn_grad_u, x_vars, x_dvv,
				    cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_X_to_U() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_X_to_U(const RealVector& fn_grad_x,   RealVector& fn_grad_u,
		  const RealMatrix& jacobian_xu, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_X_to_U(fn_grad_x, fn_grad_u, jacobian_xu, x_dvv,
				    cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_X_to_U() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
		  const RealVector& x_vars, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids,
		  SizetMultiArrayConstView acv_ids,
		  const SizetArray& acv_map1_indices,
		  const ShortArray& acv_map2_targets)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_X_to_S(fn_grad_x, fn_grad_s, x_vars, x_dvv, cv_ids,
				    acv_ids, acv_map1_indices,
				    acv_map2_targets);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_X_to_S() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << "class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
		  const RealMatrix& jacobian_xs, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids,
		  SizetMultiArrayConstView acv_ids,
		  const SizetArray& acv_map1_indices,
		  const ShortArray& acv_map2_targets)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_X_to_S(fn_grad_x, fn_grad_s, jacobian_xs, x_dvv,
				    cv_ids, acv_ids, acv_map1_indices,
				    acv_map2_targets);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_X_to_S() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_U_to_X(const RealVector& fn_grad_u, RealVector& fn_grad_x,
		  const RealVector& x_vars,    const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_U_to_X(fn_grad_u, fn_grad_x, x_vars, x_dvv,
				    cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_U_to_X() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_U_to_X(const RealVector& fn_grad_u,   RealVector& fn_grad_x,
		  const RealMatrix& jacobian_ux, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_U_to_X(fn_grad_u, fn_grad_x, jacobian_ux, x_dvv,
				    cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_U_to_X() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_hess_X_to_U(const RealSymMatrix& fn_hess_x, RealSymMatrix& fn_hess_u,
		  const RealVector& x_vars, const RealVector& fn_grad_x,
		  const SizetArray& x_dvv, SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_hess_X_to_U(fn_hess_x, fn_hess_u, x_vars, fn_grad_x,
				    x_dvv, cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_hess_X_to_U() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_hess_X_to_U(const RealSymMatrix& fn_hess_x, RealSymMatrix& fn_hess_u,
		  const RealMatrix& jacobian_xu,
		  const RealSymMatrixArray& hessian_xu,
		  const RealVector& fn_grad_x, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_hess_X_to_U(fn_hess_x, fn_hess_u, jacobian_xu,
				    hessian_xu, fn_grad_x, x_dvv, cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_hess_X_to_U() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
jacobian_dX_dU(const RealVector& x_vars, RealMatrix& jacobian_xu)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->jacobian_dX_dU(x_vars, jacobian_xu);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine jacobian_dX_dU() virtual "
          << "fn.\nNo default defined at ProbabilityTransformation base class."
	  << "\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
jacobian_dU_dX(const RealVector& x_vars, RealMatrix& jacobian_ux)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->jacobian_dU_dX(x_vars, jacobian_ux);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine jacobian_dU_dX() virtual "
          << "fn.\nNo default defined at ProbabilityTransformation base class."
	  << "\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
jacobian_dX_dS(const RealVector& x_vars, RealMatrix& jacobian_xs,
	       SizetMultiArrayConstView cv_ids,
	       SizetMultiArrayConstView acv_ids,
	       const SizetArray& acv_map1_indices,
	       const ShortArray& acv_map2_targets)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->jacobian_dX_dS(x_vars, jacobian_xs, cv_ids, acv_ids,
				 acv_map1_indices, acv_map2_targets);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine jacobian_dX_dS() virtual "
          << "fn.\nNo default defined at ProbabilityTransformation base class."
	  << "\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
hessian_d2X_dU2(const RealVector& x_vars, RealSymMatrixArray& hessian_xu)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->hessian_d2X_dU2(x_vars, hessian_xu);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine hessian_d2X_dU2() virtual "
          << "fn.\nNo default defined at ProbabilityTransformation base class."
	  << "\n" << std::endl;
    abort_handler(-1);
  }
}


/** This procedure computes numerical derivatives of x and/or z with respect to
    distribution parameters s, and is used by jacobian_dX_dS() to provide data
    that is not available analytically.  Numerical dz/ds involves dL/ds
    (z(s) = L(s) u and dz/ds = dL/ds u) and is needed to evaluate dx/ds
    semi-analytically for correlated variables.  Numerical dx/ds is needed for
    distributions lacking simple closed-form CDF expressions (beta and gamma
    distributions). */
void ProbabilityTransformation::
numerical_design_jacobian(const RealVector& x_vars,
			  bool xs, RealMatrix& num_jacobian_xs,
			  bool zs, RealMatrix& num_jacobian_zs,
			  SizetMultiArrayConstView cv_ids,
			  SizetMultiArrayConstView acv_ids,
			  const SizetArray& acv_map1_indices,
			  const ShortArray& acv_map2_targets)
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
  size_t i, j, k, num_var_map_1c = acv_map1_indices.size();
  int x_len = x_vars.length();
  if (xs && (num_jacobian_xs.numRows() != x_len ||
	     num_jacobian_xs.numCols() != num_var_map_1c) )
    num_jacobian_xs.shape(x_len, num_var_map_1c);
  if (zs && (num_jacobian_zs.numRows() != x_len ||
	     num_jacobian_zs.numCols() != num_var_map_1c) )
    num_jacobian_zs.shape(x_len, num_var_map_1c);

  RealMatrix L_s_plus_h, dL_dsi;
  RealVector dz_dsi;
  //RealVector z_vars_s_plus_h, z_vars_s_minus_h;
  RealVector x_vars_s_plus_h, x_vars_s_minus_h;
  if (zs) {
    L_s_plus_h.shape(x_len, x_len);
    dL_dsi.shape(x_len, x_len);
    dz_dsi.size(x_len);
  }

  RealVector u_vars;
  trans_X_to_U(x_vars, u_vars);

  Real fd_grad_ss = 1.e-4;
  for (i=0; i<num_var_map_1c; i++) {

    size_t cv_index        = find_index(cv_ids, acv_ids[acv_map1_indices[i]]);
    short  acv_map2_target = acv_map2_targets[i];
    if (cv_index != _NPOS && acv_map2_target != NO_TARGET) {

      Real s0 = distribution_parameter(cv_index, acv_map2_target);

      // Compute the offset for the ith gradient variable.
      // Enforce a minimum delta of fdgss*.01
      Real h_mag = fd_grad_ss * std::max(std::fabs(s0), .01);
      Real h = (s0 < 0.0) ? -h_mag : h_mag; // h has same sign as s0

      // -----------------------------------
      // Evaluate (L/z_vars/x_vars)_s_plus_h
      // -----------------------------------
      Real s1 = s0 + h;
      // updates ranVars & corrCholeskyFactorZ:
      distribution_parameter(cv_index, acv_map2_target, s1);
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
      // updates ranVars & corrCholeskyFactorZ:
      distribution_parameter(cv_index, acv_map2_target, s1);
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
	dz_dsi.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1., dL_dsi,
			u_vars, 0.); // dz/ds
	for (j=0; j<x_len; j++)
	  num_jacobian_zs(j, i) = dz_dsi(j);
	//for (j=0; j<x_len; j++)                          // dz/ds (alt)
	//  num_jacobian_zs(j, i)=(z_vars_s_plus_h(j)-z_vars_s_minus_h(j))/2./h;
      }
      if (xs)
	for (j=0; j<x_len; j++)                            // dx/ds
	  num_jacobian_xs(j,i) = (x_vars_s_plus_h(j)-x_vars_s_minus_h(j))/2./h;

      // resets s0 & corrCholeskyFactorZ:
      distribution_parameter(cv_index, acv_map2_target, s0);
    }
  }
}


/** This function accommodates the native Model space (X, Z, or U)
    through the use of num[Distribution]Vars counts and
    iteratedModel.some_distribution_parameter(), but is tied to
    x-space through the VarMapIndices user specifications. */
const Real& ProbabilityTransformation::
distribution_parameter(size_t index, short target)
{
  switch (target) {
  case CDV_LWR_BND: case N_LWR_BND: case LN_LWR_BND: case   U_LWR_BND:
  case  LU_LWR_BND: case T_LWR_BND: case BE_LWR_BND: case CSV_LWR_BND:
    return ranVarLowerBndsX[index];      break;
  case CDV_UPR_BND: case N_UPR_BND: case LN_UPR_BND: case   U_UPR_BND:
  case  LU_UPR_BND: case T_UPR_BND: case BE_UPR_BND: case CSV_UPR_BND:
    return ranVarUpperBndsX[index];      break;
  case N_MEAN:      case LN_MEAN:
    return ranVarMeansX[index];          break;
  case N_STD_DEV:   case LN_STD_DEV:
    return ranVarStdDevsX[index];        break;
  case LN_LAMBDA:   case T_MODE:    case E_BETA:  case BE_ALPHA: case GA_ALPHA:
  case GU_ALPHA:    case F_ALPHA:   case W_ALPHA:
    return ranVarAddtlParamsX[index][0]; break;
  case LN_ZETA:     case BE_BETA:   case GA_BETA: case GU_BETA:  case F_BETA:
  case W_BETA:
    return ranVarAddtlParamsX[index][1]; break;
  case LN_ERR_FACT:
    if (ranVarAddtlParamsX[index].length() < 3) {
      PCerr << "Error: LN_ERR_FACT cannot be returned in Probability"
	    << "Transformation::distribution_parameter()." << std::endl;
      abort_handler(-1);
    }
    return ranVarAddtlParamsX[index][2];
    break;
  }
}


/** This function accommodates the native Model space (X, Z, or U)
    through the use of num[Distribution]Vars counts and
    iteratedModel.some_distribution_parameter(), but is tied to
    x-space through the VarMapIndices user specifications. */
void ProbabilityTransformation::
distribution_parameter(size_t index, short target, const Real& param)
{
  switch (target) {
  // -------------------------
  // Distribution lower bounds
  // -------------------------
  case CDV_LWR_BND: case U_LWR_BND: case CSV_LWR_BND:
    ranVarLowerBndsX[index] = param;
    moments_from_uniform_params(param, ranVarUpperBndsX[index],
				ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case N_LWR_BND:
    if (ranVarTypesX[index] != BOUNDED_NORMAL) { // protect this for now
      PCerr << "Error: setting normal bounds only allowed for BOUNDED_NORMAL."
	    << std::endl;
      abort_handler(-1);
    }
    //if (ranVarLowerBndsX[index] == -DBL_MAX && param != -DBL_MAX) {
    //  PCerr << "Error: for BOUNDED_NORMAL, activating an inactive lower bound"
    //        << " is not allowed." << std::endl;
    //  abort_handler(-1);
    //}
    //else
    ranVarLowerBndsX[index] = param; break;
  case LN_LWR_BND:
    if (ranVarTypesX[index] != BOUNDED_LOGNORMAL) { // protect this for now
      PCerr << "Error: setting lognormal bounds only allowed for "
	    << "BOUNDED_LOGNORMAL." << std::endl;
      abort_handler(-1);
    }
    //if (ranVarLowerBndsX[index] == 0. && param != 0.) {
    //  PCerr << "Error: for BOUNDED_LOGNORMAL, activating an inactive lower "
    //        << "bound is not allowed." << std::endl;
    //  abort_handler(-1);
    //}
    //else
    ranVarLowerBndsX[index] = param; break;
  case LU_LWR_BND:
    ranVarLowerBndsX[index] = param;
    moments_from_loguniform_params(param, ranVarUpperBndsX[index],
				   ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case T_LWR_BND:
    ranVarLowerBndsX[index] = param;
    moments_from_triangular_params(param, ranVarUpperBndsX[index],
				   ranVarAddtlParamsX[index][0],
				   ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case BE_LWR_BND:
    ranVarLowerBndsX[index] = param;
    moments_from_beta_params(param, ranVarUpperBndsX[index],
			     ranVarAddtlParamsX[index][0],
			     ranVarAddtlParamsX[index][1],
			     ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  // -------------------------
  // Distribution upper bounds
  // -------------------------
  case CDV_UPR_BND: case U_UPR_BND: case CSV_UPR_BND:
    ranVarUpperBndsX[index] = param;
    moments_from_uniform_params(ranVarLowerBndsX[index], param,
				ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case N_UPR_BND:
    if (ranVarTypesX[index] != BOUNDED_NORMAL) { // protect this for now
      PCerr << "Error: setting normal bounds only allowed for BOUNDED_NORMAL."
	    << std::endl;
      abort_handler(-1);
    }
    //if (ranVarUpperBndsX[index] == DBL_MAX && param != DBL_MAX) {
    //  PCerr << "Error: for BOUNDED_NORMAL, activating an inactive upper "
    //        << "bound is not allowed." << std::endl;
    //  abort_handler(-1);
    //}
    //else
    ranVarUpperBndsX[index] = param; break;
  case LN_UPR_BND:
    if (ranVarTypesX[index] != BOUNDED_LOGNORMAL) { // protect this for now
      PCerr << "Error: setting lognormal bounds only allowed for "
	    << "BOUNDED_LOGNORMAL." << std::endl;
      abort_handler(-1);
    }
    //if (ranVarUpperBndsX[index] == DBL_MAX && param != DBL_MAX) {
    //  PCerr << "Error: for BOUNDED_LOGNORMAL, activating an inactive upper "
    //        << "bound is not allowed." << std::endl;
    //  abort_handler(-1);
    //}
    //else
    ranVarUpperBndsX[index] = param; break;
  case LU_UPR_BND:
    ranVarUpperBndsX[index] = param;
    moments_from_loguniform_params(ranVarLowerBndsX[index], param,
				   ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case T_UPR_BND:
    ranVarUpperBndsX[index] = param;
    moments_from_triangular_params(ranVarLowerBndsX[index], param,
				   ranVarAddtlParamsX[index][0],
				   ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case BE_UPR_BND:
    ranVarUpperBndsX[index] = param;
    moments_from_beta_params(ranVarLowerBndsX[index], param,
			     ranVarAddtlParamsX[index][0],
			     ranVarAddtlParamsX[index][1],
			     ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  // -------------------------
  // Distribution shape params
  // -------------------------
  case N_MEAN:
    ranVarMeansX[index]   = param; break;
  case N_STD_DEV:
    ranVarStdDevsX[index] = param; break;
  case LN_MEAN:
    ranVarMeansX[index] = param;
    lognormal_params_from_moments(param, ranVarStdDevsX[index],
				  ranVarAddtlParamsX[index][0],
				  ranVarAddtlParamsX[index][1]);
    if (ranVarAddtlParamsX[index].length() > 2)
      lognormal_err_factor_from_std_deviation(param, ranVarStdDevsX[index],
					      ranVarAddtlParamsX[index][2]);
    break;
  case LN_STD_DEV:
    ranVarStdDevsX[index] = param;
    lognormal_params_from_moments(ranVarMeansX[index], param,
				  ranVarAddtlParamsX[index][0],
				  ranVarAddtlParamsX[index][1]);
    if (ranVarAddtlParamsX[index].length() > 2)
      lognormal_err_factor_from_std_deviation(ranVarMeansX[index], param,
					      ranVarAddtlParamsX[index][2]);
    break;
  case LN_LAMBDA:
    ranVarAddtlParamsX[index][0] = param;
    moments_from_lognormal_params(param, ranVarAddtlParamsX[index][1],
				  ranVarMeansX[index], ranVarStdDevsX[index]);
    if (ranVarAddtlParamsX[index].length() > 2)
      lognormal_err_factor_from_std_deviation(ranVarMeansX[index],
					      ranVarStdDevsX[index],
					      ranVarAddtlParamsX[index][2]);
    break;
  case LN_ZETA:
    ranVarAddtlParamsX[index][1] = param;
    moments_from_lognormal_params(ranVarAddtlParamsX[index][0], param,
				  ranVarMeansX[index], ranVarStdDevsX[index]);
    if (ranVarAddtlParamsX[index].length() > 2)
      lognormal_err_factor_from_std_deviation(ranVarMeansX[index],
					      ranVarStdDevsX[index],
					      ranVarAddtlParamsX[index][2]);
    break;
  case LN_ERR_FACT:
    if (ranVarAddtlParamsX[index].length() < 3) {
      PCerr << "Error: LN_ERR_FACT cannot be set in Probability"
	    << "Transformation::distribution_parameter()." << std::endl;
      abort_handler(-1);
    }
    ranVarAddtlParamsX[index][2] = param;
    lognormal_std_deviation_from_err_factor(ranVarMeansX[index], param,
					    ranVarStdDevsX[index]);
    lognormal_params_from_moments(ranVarMeansX[index], ranVarStdDevsX[index],
				  ranVarAddtlParamsX[index][0],
				  ranVarAddtlParamsX[index][1]);
    break;
  case T_MODE:
    ranVarAddtlParamsX[index][0] = param;
    moments_from_triangular_params(ranVarLowerBndsX[index],
				   ranVarUpperBndsX[index], param,
				   ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case E_BETA:
    ranVarAddtlParamsX[index][0] = param;
    moments_from_exponential_params(param, ranVarMeansX[index],
				    ranVarStdDevsX[index]);
    break;
  case BE_ALPHA:
    ranVarAddtlParamsX[index][0] = param;
    moments_from_beta_params(ranVarLowerBndsX[index], ranVarUpperBndsX[index],
			     param, ranVarAddtlParamsX[index][1],
			     ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case GA_ALPHA:
    ranVarAddtlParamsX[index][0] = param;
    moments_from_gamma_params(param, ranVarAddtlParamsX[index][1],
			      ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case GU_ALPHA:
    ranVarAddtlParamsX[index][0] = param;
    moments_from_gumbel_params(param, ranVarAddtlParamsX[index][1],
			       ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case F_ALPHA:
    ranVarAddtlParamsX[index][0] = param;
    moments_from_frechet_params(param, ranVarAddtlParamsX[index][1],
				ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case W_ALPHA:
    ranVarAddtlParamsX[index][0] = param;
    moments_from_weibull_params(param, ranVarAddtlParamsX[index][1],
				ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case BE_BETA:
    ranVarAddtlParamsX[index][1] = param;
    moments_from_beta_params(ranVarLowerBndsX[index], ranVarUpperBndsX[index],
			     ranVarAddtlParamsX[index][0], param,
			     ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case GA_BETA:
    ranVarAddtlParamsX[index][1] = param;
    moments_from_gamma_params(ranVarAddtlParamsX[index][0], param,
			      ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case GU_BETA:
    ranVarAddtlParamsX[index][1] = param;
    moments_from_gumbel_params(ranVarAddtlParamsX[index][0], param,
			       ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  case F_BETA:
    ranVarAddtlParamsX[index][1] = param;
    moments_from_frechet_params(ranVarAddtlParamsX[index][0], param,
				ranVarMeansX[index], ranVarStdDevsX[index]);
    break;  
  case W_BETA:
    ranVarAddtlParamsX[index][1] = param;
    moments_from_weibull_params(ranVarAddtlParamsX[index][0], param,
				ranVarMeansX[index], ranVarStdDevsX[index]);
    break;
  }

  // update corrCholeskyFactorZ for new ranVarMeans/ranVarStdDevs
  transform_correlations();
}


#ifdef DERIV_DEBUG
void ProbabilityTransformation::
verify_trans_jacobian_hessian(const RealVector& v0)
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
	Real h_mag = fd_grad_ss * std::max(std::fabs(v0(j)), .01);
	Real h = (v1(j) < 0.0) ? -h_mag : h_mag; // h has same sign as v1(j)

	// ----------------------------
	// Evaluate trans_vars_v_plus_h
	// ----------------------------
	RealVector trans_vars_v_plus_h;
	v1(j) = v0(j) + h;
	PCout << ">>>>> Pecos finite difference gradient evaluation for v["
	      << j+1 << "] + h:\n";
	//trans_X_to_U(v1, trans_vars_v_plus_h);
	trans_U_to_X(v1, trans_vars_v_plus_h);

	// -----------------------------
	// Evaluate trans_vars_v_minus_h
	// -----------------------------
	RealVector trans_vars_v_minus_h;
	v1(j) = v0(j) - h;
	PCout << ">>>>> Pecos finite difference gradient evaluation for v["
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
	  Real h_mag = fd_hess_by_fn_ss * std::max(std::fabs(v0(j)), .01);
	  Real h = (v1(j) < 0.0) ? -h_mag : h_mag; // h has same sign as v1(j)

	  // evaluate diagonal term

	  // -----------------------------
	  // Evaluate trans_vars_v_plus_2h
	  // -----------------------------
	  v1(j) = v0(j) + 2.*h;
	  PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
	        << j+1 << "] + 2h:\n";
	  //trans_X_to_U(v1, trans_vars_v_plus_2h);
	  trans_U_to_X(v1, trans_vars_v_plus_2h);

	  // ------------------------------
	  // Evaluate trans_vars_v_minus_2h
	  // ------------------------------
	  v1(j) = v0(j) - 2.*h;
	  PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
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
	    PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
		  << j+1 << "] + h, v[" << k+1 << "] + h:\n";
	    //trans_X_to_U(v1, trans_vars_v_plus_h_plus_h);
	    trans_U_to_X(v1, trans_vars_v_plus_h_plus_h);
	    // ------------------------------------
	    // Evaluate trans_vars_v_plus_h_minus_h
	    // ------------------------------------
	    //v1(j) = v0(j) + h;
	    v1(k) = v0(k) - h;
	    PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
		  << j+1 << "] + h, v[" << k+1 << "] - h:\n";
	    //trans_X_to_U(v1, trans_vars_v_plus_h_minus_h);
	    trans_U_to_X(v1, trans_vars_v_plus_h_minus_h);
	    // ------------------------------------
	    // Evaluate trans_vars_v_minus_h_plus_h
	    // ------------------------------------
	    v1(j) = v0(j) - h;
	    v1(k) = v0(k) + h;
	    PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
		  << j+1 << "] - h, v[" << k+1 << "] + h:\n";
	    //trans_X_to_U(v1, trans_vars_v_minus_h_plus_h);
	    trans_U_to_X(v1, trans_vars_v_minus_h_plus_h);
	    // -------------------------------------
	    // Evaluate trans_vars_v_minus_h_minus_h
	    // -------------------------------------
	    //v1(j) = v0(j) - h;
	    v1(k) = v0(k) - h;
	    PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
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
	  Real h_mag = fd_hess_by_grad_ss * std::max(std::fabs(v0(j)), .01);
	  Real h = (v1(j) < 0.0) ? -h_mag : h_mag; // h has same sign as v1(j)

	  // --------------------------
	  // Evaluate fn_grads_v_plus_h
	  // --------------------------
	  v1(j) = v0(j) + h;
	  PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
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
  PCout << "\nNumerical jacobian:" << num_jac_dtv_dv
        << "\nAnalytic jacobian:"  << jacobian_xu; //jacobian_ux;
  RealSymMatrixArray hessian_xu(num_tv);
  hessian_d2X_dU2(x0, hessian_xu);
  for (i=0; i<num_tv; i++)
    PCout << "\nNumerical Hessian:" << num_hess_d2tv_dv2[i]
	  << "\nAnalytic Hessian:"  << hessian_xu[i];
}


void ProbabilityTransformation::verify_design_jacobian(const RealVector& u0)
{
  RealVector x0;
  trans_U_to_X(u0, x0);

  RealMatrix num_jac_dx_ds, num_jac_dz_ds;
  numerical_design_jacobian(x0, true, num_jac_dx_ds, false, num_jac_dz_ds);

  RealMatrix jacobian_xs;
  jacobian_dX_dS(x0, jacobian_xs);

  // Print out numerical and analytic:
  PCout << "\nNumerical jacobian:" << num_jac_dx_ds
        << "\nAnalytic jacobian:"  << jacobian_xs;
}
#endif // DERIV_DEBUG

} // namespace Pecos
