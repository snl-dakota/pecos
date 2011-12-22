/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        NodalInterpPolyApproximation
//- Description:  Implementation code for NodalInterpPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "NodalInterpPolyApproximation.hpp"
#include "TensorProductDriver.hpp"
#include "SparseGridDriver.hpp"

//#define DEBUG


namespace Pecos {


Real NodalInterpPolyApproximation::
tensor_product_value(const RealVector& x, const RealVector& exp_t1_coeffs,
		     const RealMatrix& exp_t2_coeffs,
		     const UShortArray& basis_index, const UShort2DArray& key,
		     const SizetArray& colloc_index)
{
  Real tp_val = 0.; size_t i, num_colloc_pts = key.size();
  if (exp_t2_coeffs.empty()) {
    if (colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[i] * 
	          type1_interpolant_value(x, key[i], basis_index);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[colloc_index[i]] *
	          type1_interpolant_value(x, key[i], basis_index);
  }
  else {
    size_t j, c_index;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      tp_val += exp_t1_coeffs[c_index] *
	        type1_interpolant_value(x, key_i, basis_index);
      const Real* exp_t2_coeff_i = exp_t2_coeffs[c_index];
      for (j=0; j<numVars; ++j)
	tp_val += exp_t2_coeff_i[j] *
	          type2_interpolant_value(x, j, key_i, basis_index);
    }
  }
  return tp_val;
}


const RealVector& NodalInterpPolyApproximation::
tensor_product_gradient_basis_variables(const RealVector& x,
					const RealVector& exp_t1_coeffs,
					const RealMatrix& exp_t2_coeffs,
					const UShortArray& basis_index,
					const UShort2DArray& key,
					const SizetArray& colloc_index)
{
  if (tpGradient.length() != numVars)
    tpGradient.size(numVars); // init to 0
  else
    tpGradient = 0.;
  size_t i, j, num_colloc_pts = key.size();
  if (exp_t2_coeffs.empty()) {
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += exp_t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index);
    }
  }
  else {
    size_t k;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      const Real* exp_t2_coeff_i = (colloc_index.empty()) ?
	exp_t2_coeffs[i] : exp_t2_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j) { // ith contribution to jth grad component
	tpGradient[j] += exp_t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index);
	for (k=0; k<numVars; ++k) // type2 interpolant for kth grad comp
	  tpGradient[j] += exp_t2_coeff_i[k] *
	    type2_interpolant_gradient(x, j, k, key_i, basis_index);
      }
    }
  }
  return tpGradient;
}


const RealVector& NodalInterpPolyApproximation::
tensor_product_gradient_basis_variables(const RealVector& x,
					const RealVector& exp_t1_coeffs,
					const RealMatrix& exp_t2_coeffs,
					const UShortArray& basis_index,
					const UShort2DArray& key,
					const SizetArray& colloc_index,
					const SizetArray& dvv)
{
  size_t i, j, deriv_index, num_colloc_pts = key.size(),
    num_deriv_vars = dvv.size();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.size(num_deriv_vars); // init to 0
  else
    tpGradient = 0.;
  if (exp_t2_coeffs.empty()) {
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i   = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += exp_t1_coeff_i *
	  type1_interpolant_gradient(x, deriv_index, key_i, basis_index);
      }
    }
  }
  else {
    size_t k;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real& exp_t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      const Real* exp_t2_coeff_i = (colloc_index.empty()) ?
	exp_t2_coeffs[i] : exp_t2_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += exp_t1_coeff_i *
	  type1_interpolant_gradient(x, deriv_index, key_i, basis_index);
	for (k=0; k<numVars; ++k)
	  tpGradient[j] += exp_t2_coeff_i[k] *
	    type2_interpolant_gradient(x, deriv_index, k, key_i, basis_index);
      }
    }
  }
  return tpGradient;
}


const RealVector& NodalInterpPolyApproximation::
tensor_product_gradient_nonbasis_variables(const RealVector& x,
					   const RealMatrix& exp_t1_coeff_grads,
					   const UShortArray& basis_index,
					   const UShort2DArray& key,
					   const SizetArray& colloc_index)
{
  size_t i, j, num_colloc_pts = key.size(),
    num_deriv_vars = exp_t1_coeff_grads.numRows();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.size(num_deriv_vars); // init to 0
  else
    tpGradient = 0.;
  for (i=0; i<num_colloc_pts; ++i) {
    const Real* exp_t1_coeff_grad_i = (colloc_index.empty()) ?
      exp_t1_coeff_grads[i] : exp_t1_coeff_grads[colloc_index[i]];
    Real t1_val = type1_interpolant_value(x, key[i], basis_index);
    for (j=0; j<numVars; ++j)
      tpGradient[j] += exp_t1_coeff_grad_i[j] * t1_val;
  }
  return tpGradient;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
Real NodalInterpPolyApproximation::
tensor_product_mean(const RealVector&    x,   const UShortArray& lev_index,
		    const UShort2DArray& key, const SizetArray&  colloc_index)
{
  const Real3DArray& colloc_wts_1d
    = driverRep->type1_collocation_weights_array();
  Real tp_mean = 0.;
  size_t i, j, c_index, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    Real prod_i = 1.;
    for (j=0; j<numVars; ++j)
      prod_i *= (randomVarsKey[j]) ? colloc_wts_1d[lev_index[j]][j][key_i[j]] :
	polynomialBasis[lev_index[j]][j].type1_value(x[j], key_i[j]);
    c_index = (colloc_index.empty()) ? i : colloc_index[i];
    tp_mean += expansionType1Coeffs[c_index] * prod_i;
  }
  return tp_mean;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_mean_gradient(const RealVector& x, const UShortArray& lev_index,
			     const UShort2DArray& key,
			     const SizetArray& colloc_index,
			     const SizetArray& dvv)
{
  // -------------------------------------------------------------------
  // Mixed variable key:
  //   xi = ran vars, sa = augmented des vars, si = inserted design vars
  // Active variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
  //   mu(sa, si)    = Sum_i r_i(sa, si) wt_prod_i
  //   dmu/ds        = Sum_i dr_i/ds wt_prod_i
  // All variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
  //   mu(sa, si)    = Sum_i r_i(si) Lsa_i wt_prod_i
  //   dmu/dsa       = Sum_i r_i(si) dLsa_i/dsa wt_prod_i
  //   dmu/dsi       = Sum_i dr_i/dsi Lsa_i wt_prod_i
  // -------------------------------------------------------------------
  size_t i, j, k, c_index, deriv_index, cntr = 0, 
    num_deriv_vars = dvv.size(), num_colloc_pts = key.size();
  if (tpMeanGrad.length() != num_deriv_vars)
    tpMeanGrad.sizeUninitialized(num_deriv_vars);
  tpMeanGrad = 0.;
  SizetList::iterator it;
  const Real3DArray& colloc_wts_1d
    = driverRep->type1_collocation_weights_array();
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    // Error check for required data
    if (randomVarsKey[deriv_index] &&
	!expConfigOptions.expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	    << "InterpPolyApproximation::tensor_product_mean_gradient()."
	    << std::endl;
      abort_handler(-1);
    }
    else if (!randomVarsKey[deriv_index] &&
	     !expConfigOptions.expansionCoeffFlag) {
      PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	    << "Approximation::tensor_product_mean_gradient()" << std::endl;
      abort_handler(-1);
    }
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      c_index = (colloc_index.empty()) ? j : colloc_index[j];
      Real wt_prod_j = 1.;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	{ k = *it; wt_prod_j *= colloc_wts_1d[lev_index[k]][k][key_j[k]]; }
      if (randomVarsKey[deriv_index]) {
	// --------------------------------------------------------------------
	// derivative of All var expansion w.r.t. random var (design insertion)
	// --------------------------------------------------------------------
	Real Lsa_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  Lsa_j *= polynomialBasis[lev_index[k]][k].type1_value(x[k], key_j[k]);
	}
	tpMeanGrad[i] += wt_prod_j * Lsa_j *
	                 expansionType1CoeffGrads(cntr, c_index);
      }
      else {
	// ---------------------------------------------------------------------
	// deriv of All var expansion w.r.t. nonrandom var (design augmentation)
	// ---------------------------------------------------------------------
	Real dLsa_j_dsa_i = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    polynomialBasis[lev_index[k]][k].type1_gradient(x[k], key_j[k]) :
	    polynomialBasis[lev_index[k]][k].type1_value(x[k],    key_j[k]);
	}
	tpMeanGrad[i] += expansionType1Coeffs[c_index]
	              *  wt_prod_j * dLsa_j_dsa_i;
      }
    }
    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      cntr++;
  }

  return tpMeanGrad;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
Real NodalInterpPolyApproximation::
tensor_product_covariance(const RealVector& x, const UShortArray& lev_index,
			  const UShort2DArray& key,
			  const SizetArray& colloc_index,
			  const RealVector& exp_coeffs_2)
{
  size_t i, j, k, c_index_i, c_index_j, num_colloc_pts = key.size();
  Real mean1 = 0., mean2 = 0., tp_covar = 0.; SizetList::iterator it;
  const Real3DArray& colloc_wts_1d
    = driverRep->type1_collocation_weights_array();
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    Real wt_prod_i = 1., Ls_prod_i = 1.;
    for (k=0; k<numVars; ++k)
      if (randomVarsKey[k])
	wt_prod_i *= colloc_wts_1d[lev_index[k]][k][key_i[k]];
      else
	Ls_prod_i *=
	  polynomialBasis[lev_index[k]][k].type1_value(x[k], key_i[k]);
    c_index_i = (colloc_index.empty()) ? i : colloc_index[i];
    const Real& exp_coeff_i = expansionType1Coeffs[c_index_i];
    Real wt_Ls_prod_i = wt_prod_i    * Ls_prod_i;
    mean1 += exp_coeff_i             * wt_Ls_prod_i;
    mean2 += exp_coeffs_2[c_index_i] * wt_Ls_prod_i;
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      // to include the ij-th term,  basis i must be the same as basis j for the
      // random variable subset.  In this case, wt_prod_i may be reused.  Note
      // that it is not necessary to collapse terms with the same random basis
      // subset, since cross term in (a+b)(a+b) = a^2+2ab+b^2 gets included.
      // If terms were collapsed (following eval of non-random portions), the
      // nested loop could be replaced with a single loop to evaluate (a+b)^2.
      bool include = true;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	if (key_i[*it] != key_j[*it])
	  { include = false; break; }
      if (include) {
	Real Ls_prod_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  Ls_prod_j *=
	    polynomialBasis[lev_index[k]][k].type1_value(x[k], key_j[k]);
	}
	c_index_j = (colloc_index.empty()) ? j : colloc_index[j];
	tp_covar += wt_prod_i * Ls_prod_i * Ls_prod_j * exp_coeff_i
	         *  exp_coeffs_2[c_index_j];
      }
    }
  }
  tp_covar -= mean1*mean2;
  return tp_covar;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_variance_gradient(const RealVector& x,
				 const UShortArray& lev_index,
				 const UShort2DArray& key,
				 const SizetArray& colloc_index,
				 const SizetArray& dvv)
{
  // -------------------------------------------------------------------
  // Mixed variable key:
  //   xi = ran vars, sa = augmented des vars, si = inserted design vars
  // Active variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
  //   var(sa, si)   = Sum_i Sum_j r_i(sa, si) r_j(sa, si) wt_prod_ij - mu^2
  //   dvar/ds       = Sum_i Sum_j (r_i dr_j/ds + dr_i/ds r_j) wt_prod_ij -
  //                   2 mu dmu/ds
  // All variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
  //   var(sa, si)   = Sum_i Sum_j r_i(si) r_j(si) Lsa_i Lsa_j wt_prod_ij -
  //                   mu^2
  //   dvar/dsa      = Sum_i Sum_j r_i(si) r_j(si) (Lsa_i dLsa_j/dsa +
  //                   dLsa_i/dsa Lsa_j) wt_prod_ij - 2 mu dmu/dsa
  //   dvar/dsi      = Sum_i Sum_j (r_i dr_j/dsi + dr_i/dsi r_j)
  //                   Lsa_i Lsa_j wt_prod_ij - 2 mu dmu/dsi
  // -------------------------------------------------------------------
  size_t i, j, k, l, deriv_index, c_index_j, c_index_k, cntr = 0,
    num_deriv_vars = dvv.size(), num_colloc_pts = key.size();
  if (tpVarianceGrad.length() != num_deriv_vars)
    tpVarianceGrad.sizeUninitialized(num_deriv_vars);
  tpVarianceGrad = 0.;
  Real mean = 0.; SizetList::iterator it;
  const Real3DArray& colloc_wts_1d
    = driverRep->type1_collocation_weights_array();
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (randomVarsKey[deriv_index] &&
	!expConfigOptions.expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	    << "InterpPolyApproximation::tensor_product_variance_gradient()."
	    << std::endl;
      abort_handler(-1);
    }
    Real mean_grad_i = 0.;
    // first loop of double sum
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      c_index_j = (colloc_index.empty()) ? j : colloc_index[j];
      // compute wt_prod_j and Lsa_j
      Real wt_prod_j = 1., Lsa_j = 1., dLsa_j_dsa_i = 1.;
      for (k=0; k<numVars; ++k)
	if (randomVarsKey[k])
	  wt_prod_j *= colloc_wts_1d[lev_index[k]][k][key_j[k]];
	else
	  Lsa_j *= polynomialBasis[lev_index[k]][k].type1_value(x[k], key_j[k]);
      // update mean (once) and mean_grad_i
      if (i == 0)
	mean += expansionType1Coeffs[c_index_j] * wt_prod_j * Lsa_j;
      if (randomVarsKey[deriv_index])
	mean_grad_i += wt_prod_j * Lsa_j
	            *  expansionType1CoeffGrads(cntr, c_index_j);
      else {
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    polynomialBasis[lev_index[k]][k].type1_gradient(x[k], key_j[k]) :
	    polynomialBasis[lev_index[k]][k].type1_value(x[k],    key_j[k]);
	}
	mean_grad_i += wt_prod_j * dLsa_j_dsa_i
	            *  expansionType1Coeffs[c_index_j];
      }
      // second loop of double sum
      for (k=0; k<num_colloc_pts; ++k) {
	const UShortArray& key_k = key[k];
	c_index_k = (colloc_index.empty()) ? k : colloc_index[k];
	// to include the jk-th term, colloc pts xi_j must be the same as xi_k
	// for the random var subset.  In this case, wt_prod_j may be reused.
	bool include = true;
	for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	  if (key_j[*it] != key_k[*it])
	    { include = false; break; }
	if (include) {
	  Real Lsa_k = 1.;
	  for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	    l = *it;
	    Lsa_k *=
	      polynomialBasis[lev_index[l]][l].type1_value(x[l], key_k[l]);
	  }
	  if (randomVarsKey[deriv_index])
	    // ---------------------------------------------------------------
	    // deriv of All var expansion w.r.t. random var (design insertion)
	    // ---------------------------------------------------------------
	    tpVarianceGrad[i] += wt_prod_j * Lsa_j * Lsa_k *
	      (expansionType1Coeffs[c_index_j] *
	       expansionType1CoeffGrads(cntr, c_index_k) +
	       expansionType1CoeffGrads(cntr, c_index_j) *
	       expansionType1Coeffs[c_index_k]);
	  else {
	    // ---------------------------------------------------------------
	    // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	    // ---------------------------------------------------------------
	    Real dLsa_k_dsa_i = 1.;
	    for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it){
	      l = *it;
	      dLsa_k_dsa_i *= (l == deriv_index) ?
		polynomialBasis[lev_index[l]][l].type1_gradient(x[l], key_k[l]):
		polynomialBasis[lev_index[l]][l].type1_value(x[l],    key_k[l]);
	    }
	    tpVarianceGrad[i] += wt_prod_j * expansionType1Coeffs[c_index_j] *
	      expansionType1Coeffs[c_index_k] * (Lsa_j * dLsa_k_dsa_i +
						 dLsa_j_dsa_i * Lsa_k);
	  }
	}
      }
    }

    // subtract 2 mu dmu/ds
    tpVarianceGrad[i] -= 2. * mean * mean_grad_i;

    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      cntr++;
  }
  return tpVarianceGrad;
}


Real NodalInterpPolyApproximation::value(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_value(x, expansionType1Coeffs, expansionType2Coeffs,
				tpq_driver->level_index(),
				tpq_driver->collocation_key(), colloc_index);
    break;
  }
  case SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    Real approx_val = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	approx_val += sm_coeffs[i] *
	  tensor_product_value(x, expansionType1Coeffs, expansionType2Coeffs,
			       sm_mi[i], colloc_key[i], colloc_indices[i]);
    return approx_val;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_gradient_basis_variables(x, expansionType1Coeffs,
      expansionType2Coeffs, tpq_driver->level_index(),
      tpq_driver->collocation_key(), colloc_index);
    break;
  }
  case SPARSE_GRID: {
    if (approxGradient.length() != numVars)
      approxGradient.sizeUninitialized(numVars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient_basis_variables(x,
	  expansionType1Coeffs, expansionType2Coeffs, sm_mi[i], colloc_key[i],
	  colloc_indices[i]);
	for (j=0; j<numVars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_gradient_basis_variables(x, expansionType1Coeffs,
      expansionType2Coeffs, tpq_driver->level_index(),
      tpq_driver->collocation_key(), colloc_index, dvv);
    break;
  }
  case SPARSE_GRID: {
    size_t num_deriv_vars = dvv.size();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.sizeUninitialized(num_deriv_vars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient_basis_variables(x,
	  expansionType1Coeffs, expansionType2Coeffs, sm_mi[i], colloc_key[i],
	  colloc_indices[i], dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in NodalInterp"
	  << "PolyApproximation::gradient_nonbasis_variables()" << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_gradient_nonbasis_variables(x,
      expansionType1CoeffGrads, tpq_driver->level_index(),
      tpq_driver->collocation_key(), colloc_index);
    break;
  }
  case SPARSE_GRID: {
    size_t i, j, num_deriv_vars = expansionType1CoeffGrads.numRows();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.size(num_deriv_vars); // init to 0
    else
      approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient_nonbasis_variables(
	  x, expansionType1CoeffGrads, sm_mi[i], colloc_key[i],
	  colloc_indices[i]);
	for (j=0; j<num_deriv_vars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


Real NodalInterpPolyApproximation::stored_value(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not available in "
	  << "NodalInterpPolyApproximation::stored_value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_value(x, storedExpType1Coeffs, storedExpType2Coeffs,
				storedLevMultiIndex[0], storedCollocKey[0],
				colloc_index);
    break;
  }
  case SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    Real approx_val = 0.;
    size_t i, num_smolyak_indices = storedLevCoeffs.size();
    for (i=0; i<num_smolyak_indices; ++i)
      if (storedLevCoeffs[i])
	approx_val += storedLevCoeffs[i] *
	  tensor_product_value(x, storedExpType1Coeffs, storedExpType2Coeffs,
			       storedLevMultiIndex[i], storedCollocKey[i],
			       storedCollocIndices[i]);
    return approx_val;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not available in NodalInterpPoly"
	  << "Approximation::stored_gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_gradient_basis_variables(x, storedExpType1Coeffs,
      storedExpType2Coeffs, storedLevMultiIndex[0], storedCollocKey[0],
      colloc_index);
    break;
  }
  case SPARSE_GRID: {
    if (approxGradient.length() != numVars)
      approxGradient.size(numVars); // init to 0
    else
      approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    size_t i, j, num_smolyak_indices = storedLevCoeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = storedLevCoeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient_basis_variables(x,
	  storedExpType1Coeffs, storedExpType2Coeffs, storedLevMultiIndex[i],
	  storedCollocKey[i], storedCollocIndices[i]);
	for (j=0; j<numVars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not available in Nodal"
	  << "InterpPolyApproximation::stored_gradient_nonbasis_variables()"
	  << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_gradient_nonbasis_variables(x,
      storedExpType1CoeffGrads, storedLevMultiIndex[0], storedCollocKey[0],
      colloc_index);
    break;
  }
  case SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    size_t i, j, num_smolyak_indices = storedLevCoeffs.size(),
      num_deriv_vars = storedExpType1CoeffGrads.numRows();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.size(num_deriv_vars); // init to 0
    else
      approxGradient = 0.;
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = storedLevCoeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient_nonbasis_variables(
	  x, storedExpType1CoeffGrads, storedLevMultiIndex[i],
	  storedCollocKey[i], storedCollocIndices[i]);
	for (j=0; j<num_deriv_vars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


/** In this case, all expansion variables are random variables and the
    mean of the expansion is simply the sum over i of r_i w_i. */
Real NodalInterpPolyApproximation::mean()
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  // TO DO:
  //if (!driverRep->track_unique_product_weights()) {
  //  PCerr << "Error: unique product weights required in "
  //	  << "NodalInterpPolyApproximation::mean()" << std::endl;
  //  abort_handler(-1);
  //}

  if (numericalMoments.empty())
    numericalMoments.sizeUninitialized(4); // standard mode
  Real& mean = numericalMoments[0]; mean = 0.;
  const RealVector& t1_wts = driverRep->type1_weight_sets();
  switch (basisConfigOptions.useDerivs) {
  case false:
    for (size_t i=0; i<numCollocPts; ++i)
      mean += expansionType1Coeffs[i] * t1_wts[i];
    break;
  case true: {
    size_t i, j;
    const RealMatrix& t2_wts = driverRep->type2_weight_sets();
    for (i=0; i<numCollocPts; ++i) {
      mean += expansionType1Coeffs[i] * t1_wts[i];
      const Real* coeff2_i = expansionType2Coeffs[i];
      const Real*  t2_wt_i = t2_wts[i];
      for (j=0; j<numVars; ++j)
	mean += coeff2_i[j] * t2_wt_i[j];
    }
    break;
  }
  }

  return mean;
}


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves evaluating the
    expectation over this subset. */
Real NodalInterpPolyApproximation::mean(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_mean(x, tpq_driver->level_index(),
			       tpq_driver->collocation_key(), colloc_index);
    break;
  }
  case SPARSE_GRID: {
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    if (numericalMoments.empty())
      numericalMoments.sizeUninitialized(2); // all_variables mode
    Real& mean = numericalMoments[0]; mean = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	mean += sm_coeffs[i] *
	  tensor_product_mean(x, sm_mi[i], colloc_key[i], colloc_indices[i]);
    return mean;
    break;
  }
  }
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  In
    this case, the derivative of the expectation is the expectation of
    the derivative.  The mixed derivative case (some design variables
    are inserted and some are augmented) requires no special treatment. */
const RealVector& NodalInterpPolyApproximation::mean_gradient()
{
  // d/ds <R> = <dR/ds>

  // Error check for required data
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	  << "InterpPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  const RealVector& t1_wts = driverRep->type1_weight_sets();
  size_t i, j, num_deriv_vars = expansionType1CoeffGrads.numRows();
  if (meanGradient.length() != num_deriv_vars)
    meanGradient.sizeUninitialized(num_deriv_vars);
  meanGradient = 0.;
  for (i=0; i<numCollocPts; ++i) {
    const Real& t1_wt_i = t1_wts[i];
    for (j=0; j<num_deriv_vars; ++j)
      meanGradient[j] += expansionType1CoeffGrads(j,i) * t1_wt_i;
  }
  return meanGradient;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  In this case, the mean of the expansion is the
    expectation over the random subset and the derivative of the mean
    is the derivative of the remaining expansion over the non-random
    subset.  This function must handle the mixed case, where some
    design/state variables are augmented (and are part of the
    expansion: derivatives are evaluated as described above) and some
    are inserted (derivatives are obtained from expansionType1CoeffGrads). */
const RealVector& NodalInterpPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_mean_gradient(x, tpq_driver->level_index(),
					tpq_driver->collocation_key(),
					colloc_index, dvv);
    break;
  }
  case SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (meanGradient.length() != num_deriv_vars)
      meanGradient.sizeUninitialized(num_deriv_vars);
    meanGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      if (coeff) {
	const RealVector& tpm_grad
	  = tensor_product_mean_gradient(x, sm_mi[i], colloc_key[i],
					 colloc_indices[i], dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  meanGradient[j] += coeff * tpm_grad[j];
      }
    }
    return meanGradient;
    break;
  }
}


/** In this case, all expansion variables are random variables and the
    variance of the expansion is the sum over all but the first term
    of the coefficients squared times the polynomial norms squared. */
Real NodalInterpPolyApproximation::variance()
{
  if (numericalMoments.empty())
    numericalMoments.sizeUninitialized(4); // standard mode
  numericalMoments[1] = covariance(this);
  return numericalMoments[1];
}


/** In this case, a subset of the expansion variables are random
    variables and the variance of the expansion involves summations
    over this subset. */
Real NodalInterpPolyApproximation::variance(const RealVector& x)
{
  if (numericalMoments.empty())
    numericalMoments.sizeUninitialized(2); // all_variables mode
  numericalMoments[1] = covariance(x, this);
  return numericalMoments[1];
}


Real NodalInterpPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  // TO DO:
  //if (!driverRep->track_unique_product_weights()) {
  //  PCerr << "Error: unique product weights required in "
  //	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
  //  abort_handler(-1);
  //}

  // compute mean_1,mean_2 first, then compute covariance as
  // wt_prod*(coeff1-mean_1)*(coeff2-mean_2) in order to avoid precision
  // loss from computing covariance as <R_i R_j> - \mu_i \mu_j
  NodalInterpPolyApproximation* nip_approx_2
    = (NodalInterpPolyApproximation*)poly_approx_2;
  // Note: compute_statistics() in dakota/src/NonDExpansion.C orders calls
  //       to reduce repetition in moment calculations.
  Real  mean_1 = mean(), mean_2 = nip_approx_2->mean();
  const RealVector& t1_coeffs_2 = nip_approx_2->expansionType1Coeffs;
  const RealVector& t1_wts      = driverRep->type1_weight_sets();
  Real covar = 0.; size_t i, j;
  switch (basisConfigOptions.useDerivs) {
  case false: // type1 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
    for (i=0; i<numCollocPts; ++i)
      covar += (expansionType1Coeffs[i] - mean_1) * (t1_coeffs_2[i] - mean_2)
	    *  t1_wts[i];
    break;
  case true: {
    const RealMatrix& t2_coeffs_2 = nip_approx_2->expansionType2Coeffs;
    const RealMatrix& t2_wts      = driverRep->type2_weight_sets();
    for (i=0; i<numCollocPts; ++i) {
      // type1 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      Real coeff1_i_mm1 = expansionType1Coeffs[i] - mean_1,
	  coeff1_2i_mm2 = t1_coeffs_2[i]          - mean_2;
      covar += coeff1_i_mm1 * coeff1_2i_mm2 * t1_wts[i];
      // type2 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      // --> interpolated gradients are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
      const Real *coeff2_i  = expansionType2Coeffs[i],
	         *coeff2_2i = t2_coeffs_2[i], *t2_wt_i = t2_wts[i];
      for (j=0; j<numVars; ++j)
	covar  += (coeff1_i_mm1 * coeff2_2i[j] + coeff1_2i_mm2 * coeff2_i[j])
	       *  t2_wt_i[j];
    }
    break;
  }
  }
  return covar;
}


/** In this case, a subset of the expansion variables are random
    variables and the variance of the expansion involves summations
    over this subset. */
Real NodalInterpPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  const RealVector& t1_coeffs_2
    = ((NodalInterpPolyApproximation*)poly_approx_2)->expansionType1Coeffs;
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_covariance(x, tpq_driver->level_index(),
				     tpq_driver->collocation_key(),
				     colloc_index, t1_coeffs_2);
    break;
  }
  case SPARSE_GRID:
    // *** TO DO: verify correctness of TP summation approach
    // *** More rigorous: collapse into unique multiIndex of random terms
    //     prior to sum squared [or use a nested loop over Smolyak indices
    //     with extended tensor_product_covariance(key1,key2,etc.) definition
    //     --> good way to explore why cross-terms appear to cancel out].
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, num_smolyak_indices = sm_coeffs.size(); Real covar = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	covar += sm_coeffs[i] *
	  tensor_product_covariance(x, sm_mi[i], colloc_key[i],
				    colloc_indices[i], t1_coeffs_2);
    return covar;
    break;
  }
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  The
    mixed derivative case (some design variables are inserted and some
    are augmented) requires no special treatment. */
const RealVector& NodalInterpPolyApproximation::variance_gradient()
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::variance_gradient()" << std::endl;
    abort_handler(-1);
  }
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	  << "InterpPolyApproximation::variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  const RealVector& t1_wts = driverRep->type1_weight_sets();
  size_t i, j, num_deriv_vars = expansionType1CoeffGrads.numRows();
  if (varianceGradient.length() != num_deriv_vars)
    varianceGradient.sizeUninitialized(num_deriv_vars);
  varianceGradient = 0.;

  Real mean = 0.;
  for (i=0; i<numCollocPts; ++i)
    mean += expansionType1Coeffs[i] * t1_wts[i];
  for (i=0; i<numCollocPts; ++i) {
    Real term_i = 2. * (expansionType1Coeffs[i] - mean) * t1_wts[i];
    for (j=0; j<num_deriv_vars; ++j)
      varianceGradient[j] += term_i * expansionType1CoeffGrads(j,i);
  }

  return varianceGradient;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  This function must handle the mixed case, where
    some design/state variables are augmented (and are part of the
    expansion) and some are inserted (derivatives are obtained from
    expansionType1CoeffGrads). */
const RealVector& NodalInterpPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::variance_gradient()" << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_variance_gradient(x, tpq_driver->level_index(),
					    tpq_driver->collocation_key(),
					    colloc_index, dvv);
    break;
  }
  case SPARSE_GRID:
    // *** TO DO: verify correctness of TP summation approach
    // *** More rigorous: collapse into unique multiIndex of random terms
    //     prior to sum squared.
    size_t num_deriv_vars = dvv.size();
    if (varianceGradient.length() != num_deriv_vars)
      varianceGradient.sizeUninitialized(num_deriv_vars);
    varianceGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      if (coeff) {
	const RealVector& tpv_grad
	  = tensor_product_variance_gradient(x, sm_mi[i], colloc_key[i],
					     colloc_indices[i], dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  varianceGradient[j] += coeff * tpv_grad[j];
      }
    }
    return varianceGradient;
    break;
  }
}

} // namespace Pecos
