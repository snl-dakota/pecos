/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchInterpPolyApproximation
//- Description:  Implementation code for InterpPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "HierarchInterpPolyApproximation.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG

namespace Pecos {


void HierarchInterpPolyApproximation::allocate_expansion_coefficients()
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort4DArray& key = hsg_driver->collocation_key();
  size_t i, j, k, num_levels = key.size(), num_sets, num_tp_pts,
    num_deriv_vars = surrData.num_derivative_variables();

  if (expansionType1Coeffs.size() != num_levels)
    expansionType1Coeffs.resize(num_levels);
  if ( expansionType2Coeffs.size() != num_levels)
    expansionType2Coeffs.resize(num_levels);
  if (expansionType1CoeffGrads.size() != num_levels)
    expansionType1CoeffGrads.resize(num_levels);
  for (i=0; i<num_levels; ++i) {
    const UShort3DArray& key_i = key[i];
    num_sets = key_i.size();
    if (expansionType1Coeffs[i].size() != num_sets)
      expansionType1Coeffs[i].resize(num_sets);
    if (expansionType2Coeffs[i].size() != num_sets)
      expansionType2Coeffs[i].resize(num_sets);
    if (expansionType1CoeffGrads[i].size() != num_sets)
      expansionType1CoeffGrads[i].resize(num_sets);
    for (j=0; j<num_sets; ++j) {
      num_tp_pts = key_i[j].size();
      for (k=0; k<num_tp_pts; ++k) {
	if (expConfigOptions.expansionCoeffFlag) {
	  expansionType1Coeffs[i][j].sizeUninitialized(num_tp_pts);
	  if (basisConfigOptions.useDerivs)
	    expansionType2Coeffs[i][j].shapeUninitialized(num_deriv_vars,
							  num_tp_pts);
	}
	if (expConfigOptions.expansionCoeffGradFlag)
	  expansionType1CoeffGrads[i][j].shapeUninitialized(num_deriv_vars,
							    num_tp_pts);
      }
    }
  }

  // checking numCollocPts is insufficient due to anisotropy --> changes in
  // anisotropic weights could move points around without changing the total.
  //bool update_exp_form =
  //  ( (expConfigOptions.expansionCoeffFlag &&
  //     expansionType1Coeffs.length()      != numCollocPts) ||
  //    (expConfigOptions.expansionCoeffGradFlag &&
  //     expansionType1CoeffGrads.numCols() != numCollocPts ) );
}


void HierarchInterpPolyApproximation::compute_expansion_coefficients()
{
  if (surrData.anchor()) {
    PCerr << "Error: anchor point not supported in HierarchInterpPoly"
	  << "Approximation::compute_expansion_coefficients" << std::endl;
    abort_handler(-1);
    /*
    if (expConfigOptions.expansionCoeffFlag) {
      expansionType1Coeffs[0][0][0] = surrData.anchor_function();
      if (basisConfigOptions.useDerivs)
	Teuchos::setCol(surrData.anchor_gradient(), 0,
			expansionType2Coeffs[0][0]);
    }
    if (expConfigOptions.expansionCoeffGradFlag)
      Teuchos::setCol(surrData.anchor_gradient(), 0,
		      expansionType1CoeffGrads[0][0]);
    */
  }

  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&      sm_mi      = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      key        = hsg_driver->collocation_key();
  size_t lev, set, pt, v, num_levels = key.size(), num_sets, num_tp_pts,
    cntr = 0, num_deriv_vars = surrData.num_derivative_variables();

  // level 0
  if (expConfigOptions.expansionCoeffFlag) {
    expansionType1Coeffs[0][0][0] = surrData.response_function(cntr);
    if (basisConfigOptions.useDerivs)
      Teuchos::setCol(surrData.response_gradient(cntr), 0,
		      expansionType2Coeffs[0][0]);
  }
  if (expConfigOptions.expansionCoeffGradFlag)
    Teuchos::setCol(surrData.response_gradient(cntr), 0,
		    expansionType1CoeffGrads[0][0]);
  ++cntr;
  // levels 1 to num_levels
  for (lev=1; lev<num_levels; ++lev) {
    const UShort3DArray& key_l = key[lev];
    num_sets = key_l.size();
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = key_l[set].size();
      for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	const RealVector& c_vars = surrData.continuous_variables(cntr);
	// coefficients are hierarchical surpluses
	if (expConfigOptions.expansionCoeffFlag) {
	  expansionType1Coeffs[lev][set][pt] = surrData.response_function(cntr)
	    - value(c_vars, sm_mi, key, expansionType1Coeffs,
		    expansionType2Coeffs, lev-1);
	  if (basisConfigOptions.useDerivs) {
	    const RealVector& data_grad = surrData.response_gradient(cntr);
	    const RealVector& prev_grad = gradient_basis_variables(c_vars,
	      sm_mi, key, expansionType1Coeffs, expansionType2Coeffs, lev-1);
	    Real* hier_grad = expansionType2Coeffs[lev][set][pt];
	    for (v=0; v<num_deriv_vars; ++v)
	      hier_grad[v] = data_grad[v] - prev_grad[v];
	  }
	}
	if (expConfigOptions.expansionCoeffGradFlag) {
	  const RealVector& data_grad = surrData.response_gradient(cntr);
	  const RealVector& prev_grad = gradient_nonbasis_variables(c_vars,
	    sm_mi, key, expansionType1CoeffGrads, lev-1);
	  Real* hier_grad = expansionType1CoeffGrads[lev][set][pt];
	  for (v=0; v<num_deriv_vars; ++v)
	    hier_grad[v] = data_grad[v] - prev_grad[v];
	}
      }
    }
  }
}


void HierarchInterpPolyApproximation::store_coefficients()
{
  if (expConfigOptions.expansionCoeffFlag) {
    storedExpType1Coeffs   = expansionType1Coeffs;
    if (basisConfigOptions.useDerivs)
      storedExpType2Coeffs = expansionType2Coeffs;
  }
  if (expConfigOptions.expansionCoeffGradFlag)
    storedExpType1CoeffGrads = expansionType1CoeffGrads;

  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  storedLevMultiIndex = hsg_driver->smolyak_multi_index();
  storedCollocKey     = hsg_driver->collocation_key();
  //storedCollocIndices = hsg_driver->collocation_indices();
}


void HierarchInterpPolyApproximation::combine_coefficients(short combine_type)
{
  // update expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads} by adding or
  // multiplying stored expansion evaluated at current collocation points
  size_t i, j, num_pts = surrData.size();
  Real lf_val, discrep_val;
  /*
  for (i=0; i<num_pts; ++i) {
    const RealVector& c_vars = (anchor_pt && i == 0) ?
      surrData.anchor_continuous_variables() :
      surrData.continuous_variables(i);
    if (combine_type == MULT_COMBINE) { // eval once for both Coeffs/CoeffGrads
      discrep_val = stored_value(c_vars);
      lf_val = expansionType1Coeffs[i]; // copy prior to update
    }
    if (expConfigOptions.expansionCoeffFlag) {
      // split up type1/type2 contribs so increments are performed properly
      if (combine_type == ADD_COMBINE)
	expansionType1Coeffs[i] += stored_value(c_vars);
      else if (combine_type == MULT_COMBINE)
	expansionType1Coeffs[i] *= discrep_val;
      if (basisConfigOptions.useDerivs) {
	const RealVector& discrep_grad
	  = stored_gradient_basis_variables(c_vars);
	Real* exp_t2_coeffs_i = expansionType2Coeffs[i];
	size_t num_deriv_vars = discrep_grad.length();
	if (combine_type == ADD_COMBINE)
	  for (j=0; j<num_deriv_vars; ++j)
	    exp_t2_coeffs_i[j] += discrep_grad[j];
	else if (combine_type == MULT_COMBINE)
	  // hf = lf*discrep --> dhf/dx = dlf/dx*discrep + lf*ddiscrep/dx
	  for (j=0; j<num_deriv_vars; ++j)
	    exp_t2_coeffs_i[j] = exp_t2_coeffs_i[j] * discrep_val
	                       + discrep_grad[j]    * lf_val;
      }
    }
    if (expConfigOptions.expansionCoeffGradFlag) {
      Real* exp_t1_grad_i = expansionType1CoeffGrads[i];
      const RealVector& discrep_grad
	= stored_gradient_nonbasis_variables(c_vars);
      size_t num_deriv_vars = discrep_grad.length();
      if (combine_type == ADD_COMBINE)
	for (j=0; j<num_deriv_vars; ++j)
	  exp_t1_grad_i[j] += discrep_grad[j];
      else if (combine_type == MULT_COMBINE)
	for (j=0; j<num_deriv_vars; ++j)
	  exp_t1_grad_i[j] = exp_t1_grad_i[j] * discrep_val
	                   + discrep_grad[j]  * lf_val;
    }
  }
  */

  // clear stored data now that it has been combined
  storedExpType1Coeffs.clear();
  storedExpType2Coeffs.clear();
  storedExpType1CoeffGrads.clear();
  storedLevMultiIndex.clear();
  storedCollocKey.clear();
  //storedCollocIndices.clear();
}


void HierarchInterpPolyApproximation::restore_expansion_coefficients()
{
  size_t new_colloc_pts = surrData.size();

  /*
  if (expConfigOptions.expansionCoeffFlag) {
    expansionType1Coeffs.resize(new_colloc_pts);
    if (basisConfigOptions.useDerivs) {
      size_t num_deriv_vars = expansionType2Coeffs.numRows();
      expansionType2Coeffs.reshape(num_deriv_vars, new_colloc_pts);
    }
  }
  if (expConfigOptions.expansionCoeffGradFlag) {
    size_t num_deriv_vars = expansionType1CoeffGrads.numRows();
    expansionType1CoeffGrads.reshape(num_deriv_vars, new_colloc_pts);
  }

  for (int i=numCollocPts; i<new_colloc_pts; ++i) {
    if (expConfigOptions.expansionCoeffFlag) {
      expansionType1Coeffs[i] = surrData.response_function(i);
      if (basisConfigOptions.useDerivs)
	Teuchos::setCol(surrData.response_gradient(i), i, expansionType2Coeffs);
    }
    if (expConfigOptions.expansionCoeffGradFlag)
      Teuchos::setCol(surrData.response_gradient(i),i,expansionType1CoeffGrads);
  }
  */

  numCollocPts = new_colloc_pts;
}


Real HierarchInterpPolyApproximation::
value(const RealVector& x, const UShort3DArray& sm_mi, const UShort4DArray& key,
      const RealVector2DArray& t1_coeffs, const RealMatrix2DArray& t2_coeffs,
      unsigned short level)
{
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  Real approx_val = 0.;
  SizetArray colloc_index; // empty -> default indexing
  size_t lev, set, num_sets;
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    num_sets = sm_mi_l.size();
    for (set=0; set<num_sets; ++set)
      approx_val +=
	tensor_product_value(x, t1_coeffs_l[set], t2_coeffs_l[set],
			     sm_mi_l[set], key_l[set], colloc_index);
  }
  return approx_val;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			 const UShort4DArray& key,
			 const RealVector2DArray& t1_coeffs,
			 const RealMatrix2DArray& t2_coeffs,
			 unsigned short level)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != numVars)
    approxGradient.sizeUninitialized(numVars);
  approxGradient = 0.;

  SizetArray colloc_index; // empty -> default indexing
  size_t lev, set, num_sets;
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    num_sets = sm_mi_l.size();
    for (set=0; set<num_sets; ++set)
      approxGradient +=
	tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			 const UShort4DArray& key,
			 const RealVector2DArray& t1_coeffs,
			 const RealMatrix2DArray& t2_coeffs,
			 const SizetArray& dvv, unsigned short level)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t lev, set, num_sets, num_deriv_vars = dvv.size();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.;

  SizetArray colloc_index; // empty -> default indexing
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    num_sets = sm_mi_l.size();
    for (set=0; set<num_sets; ++set)
      approxGradient +=
	tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index, dvv);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			    const UShort4DArray& key,
			    const RealMatrix2DArray& t1_coeff_grads,
			    unsigned short level)
{
  // Error check for required data
  size_t lev, set, num_sets, num_deriv_vars;
  if (expConfigOptions.expansionCoeffGradFlag) {
    if (t1_coeff_grads.size() > level && t1_coeff_grads[level].size())
      num_deriv_vars = t1_coeff_grads[level][0].numRows();
    else {
      PCerr << "Error: insufficient size in type1 expansion coefficient "
	    << "gradients in\n       HierarchInterpPolyApproximation::"
	    << "gradient_nonbasis_variables()" << std::endl;
      abort_handler(-1);
    }
  }
  else {
    PCerr << "Error: expansion coefficient gradients not defined in Hierarch"
	  << "InterpPolyApproximation::gradient_nonbasis_variables()"
	  << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.;

  SizetArray colloc_index; // empty -> default indexing
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&            sm_mi_l = sm_mi[lev];
    const UShort3DArray&              key_l = key[lev];
    const RealMatrixArray& t1_coeff_grads_l = t1_coeff_grads[lev];
    num_sets = sm_mi_l.size();
    for (set=0; set<num_sets; ++set)
      approxGradient +=
	tensor_product_gradient_nonbasis_variables(x, t1_coeff_grads_l[set],
	  sm_mi_l[set], key_l[set], colloc_index);
  }

  return approxGradient;
}


Real HierarchInterpPolyApproximation::mean()
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  return expectation(expansionType1Coeffs, expansionType2Coeffs);
}



Real HierarchInterpPolyApproximation::mean(const RealVector& x)
{
  PCerr << "TODO: mean in all variables mode";
  return numericalMoments[0];
}


const RealVector& HierarchInterpPolyApproximation::mean_gradient()
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Hierarch"
	  << "InterpPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const RealVector2DArray& t1_wts = hsg_driver->type1_weight_set_arrays();
  /*
  size_t i, j, num_deriv_vars = expansionType1CoeffGrads.numRows();
  if (meanGradient.length() != num_deriv_vars)
    meanGradient.sizeUninitialized(num_deriv_vars);
  meanGradient = 0.;
  for (i=0; i<numCollocPts; ++i) {
    const Real& t1_wt_i = t1_wts[i];
    for (j=0; j<num_deriv_vars; ++j)
      meanGradient[j] += expansionType1CoeffGrads(j,i) * t1_wt_i;
  }
  */
  return meanGradient;
}


const RealVector& HierarchInterpPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  PCerr << "TODO: mean_gradient in all variables mode" << std::endl;
  return meanGradient;
}


Real HierarchInterpPolyApproximation::variance()
{ return covariance(this); }


Real HierarchInterpPolyApproximation::variance(const RealVector& x)
{
  //TODO
  PCerr << "TODO: variance in all variables mode" << std::endl;
  return numericalMoments[1];
}


const RealVector& HierarchInterpPolyApproximation::variance_gradient()
{
  //TODO
  PCerr << "TODO: variance_gradient()" << std::endl;
  return varianceGradient;
}


const RealVector& HierarchInterpPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  //TODO
  PCerr << "TODO: variance_gradient in all variables mode" << std::endl;
  return varianceGradient;
}


Real HierarchInterpPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  // form hierarchical t1/t2 coeffs for (R_1 - \mu_1) (R_2 - \mu_2)
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&      sm_mi      = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      key        = hsg_driver->collocation_key();
  size_t lev, set, pt, num_levels = expansionType1Coeffs.size(),
    num_sets, num_tp_pts, cntr = 0;

  HierarchInterpPolyApproximation* hip_approx_2 = 
    static_cast<HierarchInterpPolyApproximation*>(poly_approx_2);
  const SurrogateData& s_data_2 = hip_approx_2->surrData;
  Real mean_1 = mean(),  mean_2 = hip_approx_2->mean();

  RealVector2DArray cov_t1_coeffs(num_levels); cov_t1_coeffs[0].resize(1);
  RealMatrix2DArray cov_t2_coeffs(num_levels); cov_t2_coeffs[0].resize(1);
  cov_t1_coeffs[0][0].sizeUninitialized(1);
  switch (basisConfigOptions.useDerivs) {
  case false:
    // level 0
    cov_t1_coeffs[0][0][0] = (surrData.response_function(cntr) - mean_1) *
                             (s_data_2.response_function(cntr) - mean_2);
    ++cntr;
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = key[lev].size();
      cov_t1_coeffs[lev].resize(num_sets); cov_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = key[lev][set].size();
	RealVector& cov_t1_coeffs_ls = cov_t1_coeffs[lev][set];
	cov_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr)
	  cov_t1_coeffs_ls[pt]
	    = (surrData.response_function(cntr) - mean_1)
	    * (s_data_2.response_function(cntr) - mean_2)
	    - value(surrData.continuous_variables(cntr), sm_mi, key,
		    cov_t1_coeffs, cov_t2_coeffs, lev-1);
      }
    }
    break;
  case true:
    size_t v;
    // level 0
    Real data_fn1_mm1 = surrData.response_function(cntr) - mean_1;
    Real data_fn2_mm2 = s_data_2.response_function(cntr) - mean_2;
    cov_t1_coeffs[0][0][0] = data_fn1_mm1 * data_fn2_mm2;
    cov_t2_coeffs[0][0].shapeUninitialized(numVars, 1);
    Real *cov_t2_coeffs_000 = cov_t2_coeffs[0][0][0];
    const RealVector& data_grad1 = surrData.response_gradient(cntr);
    const RealVector& data_grad2 = s_data_2.response_gradient(cntr);
    for (v=0; v<numVars; ++v)
      cov_t2_coeffs_000[v] = (data_fn1_mm1 * data_grad2[v] +
			      data_fn2_mm2 * data_grad1[v]);
    ++cntr;
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = key[lev].size();
      cov_t1_coeffs[lev].resize(num_sets); cov_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& cov_t1_coeffs_ls = cov_t1_coeffs[lev][set];
	RealMatrix& cov_t2_coeffs_ls = cov_t2_coeffs[lev][set];
	num_tp_pts = key[lev][set].size();
	cov_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	cov_t2_coeffs_ls.shapeUninitialized(numVars, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  const RealVector& c_vars = surrData.continuous_variables(cntr);
	  // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  data_fn1_mm1 = surrData.response_function(cntr) - mean_1;
	  data_fn2_mm2 = s_data_2.response_function(cntr) - mean_2;
	  cov_t1_coeffs_ls[pt] = data_fn1_mm1 * data_fn2_mm2 -
	    value(c_vars, sm_mi, key, cov_t1_coeffs, cov_t2_coeffs, lev-1);
	  // type2 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  // --> interpolated grads are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
	  Real* cov_t2_coeffs_lsp = cov_t2_coeffs_ls[pt];
	  const RealVector& data_grad1 = surrData.response_gradient(cntr);
	  const RealVector& data_grad2 = s_data_2.response_gradient(cntr);
	  const RealVector& prev_grad  = gradient_basis_variables(c_vars,
	    sm_mi, key, cov_t1_coeffs, cov_t2_coeffs, lev-1);
	  for (v=0; v<numVars; ++v)
	    cov_t2_coeffs_lsp[v]
	      = (data_fn1_mm1 * data_grad2[v] + data_fn2_mm2 * data_grad1[v])
	      - prev_grad[v];
	}
      }
    }
    break;
  }

  // pass these exp coefficients into a general expectation fn
  return expectation(cov_t1_coeffs, cov_t2_coeffs);
}


Real HierarchInterpPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  //TODO
  PCerr << "TODO: covariance in all variables" << std::endl;
  Real covar = 0.;
  return covar;
}


Real HierarchInterpPolyApproximation::
expectation(const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs)
{
  Real integral = 0.;
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const RealVector2DArray& t1_wts = hsg_driver->type1_weight_set_arrays();
  size_t lev, set, pt, num_levels = t1_coeffs.size(),
    num_sets, num_tp_pts;
  switch (basisConfigOptions.useDerivs) {
  case false:
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      num_sets = t1_coeffs_l.size();
      for (set=0; set<num_sets; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt)
	  integral += t1_coeffs_ls[pt] * t1_wts_ls[pt];
      }
    }
    break;
  case true: {
    const RealMatrix2DArray& t2_wts = hsg_driver->type2_weight_set_arrays();
    size_t v;
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      num_sets = t1_coeffs_l.size();
      for (set=0; set<num_sets; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealMatrix& t2_coeffs_ls = t2_coeffs[lev][set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	const RealMatrix&    t2_wts_ls = t2_wts[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) {
	  integral += t1_coeffs_ls[pt] * t1_wts_ls[pt];
	  const Real* t2_coeffs_lsp = t2_coeffs_ls[pt];
	  const Real* t2_wts_lsp    = t2_wts_ls[pt];
	  for (v=0; v<numVars; ++v)
	    integral += t2_coeffs_lsp[v] * t2_wts_lsp[v];
	}
      }
    }
    break;
  }
  }

  return integral;
}


void HierarchInterpPolyApproximation::
compute_numerical_response_moments(size_t num_moments)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in InterpPoly"
	  << "Approximation::compute_numerical_response_moments()" << std::endl;
    abort_handler(-1);
  }

  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&      sm_mi      = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      key        = hsg_driver->collocation_key();
  size_t lev, set, pt, num_levels = key.size(), num_sets, num_tp_pts, cntr;
  int m_index, moment;

  if (numericalMoments.length() != num_moments)
    numericalMoments.sizeUninitialized(num_moments);
  Real& mean = numericalMoments[0];
  mean = expectation(expansionType1Coeffs, expansionType2Coeffs);

  // size moment coefficient arrays
  RealVector2DArray mom_t1_coeffs(num_levels);
  RealMatrix2DArray mom_t2_coeffs(num_levels);
  for (lev=0; lev<num_levels; ++lev) {
    num_sets = key[lev].size();
    mom_t1_coeffs[lev].resize(num_sets);
    mom_t2_coeffs[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = key[lev][set].size();
      mom_t1_coeffs[lev][set].sizeUninitialized(num_tp_pts);
      if (basisConfigOptions.useDerivs)
	mom_t2_coeffs[lev][set].shapeUninitialized(numVars, num_tp_pts);
    }
  }

  for (m_index=1; m_index<num_moments; ++m_index) {
    moment = m_index+1; cntr = 0;
    switch (basisConfigOptions.useDerivs) {
    case false:
      // level 0
      mom_t1_coeffs[0][0][0]
	= std::pow(surrData.response_function(cntr) - mean, moment);
      ++cntr;
      // levels 1:w
      for (lev=1; lev<num_levels; ++lev) {
	num_sets = key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  RealVector& mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  num_tp_pts = key[lev][set].size();
	  // type1 hierarchical interpolation of (R - \mu)^moment
	  for (pt=0; pt<num_tp_pts; ++pt, ++cntr)
	    mom_t1_coeffs_ls[pt]
	      = std::pow(surrData.response_function(cntr) - mean, moment)
	      - value(surrData.continuous_variables(cntr), sm_mi, key,
		      mom_t1_coeffs, mom_t2_coeffs, lev-1);
	}
      }
      break;
    case true:
      size_t v;
      // level 0
      Real data_fn_mm         = surrData.response_function(cntr) - mean;
      mom_t1_coeffs[0][0][0]  = std::pow(data_fn_mm, moment);
      Real* mom_t2_coeffs_000 = mom_t2_coeffs[0][0][0];
      Real deriv = moment * std::pow(data_fn_mm, m_index);
      const RealVector& data_grad = surrData.response_gradient(cntr);
      for (v=0; v<numVars; ++v)
	mom_t2_coeffs_000[v] = deriv * data_grad[v];
      // levels 1:w
      for (lev=1, cntr=1; lev<num_levels; ++lev) {
	num_sets = key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  num_tp_pts = key[lev][set].size();
	  RealVector& mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  RealMatrix& mom_t2_coeffs_ls = mom_t2_coeffs[lev][set];
	  for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	    const RealVector& c_vars = surrData.continuous_variables(cntr);
	    // type1 hierarchical interpolation of (R - \mu)^moment
	    data_fn_mm = surrData.response_function(cntr) - mean;
	    mom_t1_coeffs_ls[pt] = std::pow(data_fn_mm, moment) -
	      value(c_vars, sm_mi, key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	    // type2 hierarchical interpolation of (R - \mu)^moment
	    // --> interpolated grads are moment(R-\mu)^{moment-1} R'
	    Real* mom_t2_coeffs_lsp = mom_t2_coeffs_ls[pt];
	    deriv = moment * std::pow(data_fn_mm, m_index);
	    const RealVector& data_grad = surrData.response_gradient(cntr);
	    const RealVector& prev_grad = gradient_basis_variables(c_vars,
	      sm_mi, key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	    for (v=0; v<numVars; ++v)
	      mom_t2_coeffs_lsp[v] = deriv * data_grad[v] - prev_grad[v];
	  }
	}
      }
      break;
    }

    // pass these exp coefficients into a general expectation fn
    numericalMoments[m_index] = expectation(mom_t1_coeffs, mom_t2_coeffs);
  }

  // standardize third and higher central moments, if present
  standardize_moments(numericalMoments);

  /*
  if (numericalMoments.size() != num_moments)
    numericalMoments.size(num_moments);
  if (basisConfigOptions.useDerivs)
    compute_numerical_moments(expansionType1Coeffs, expansionType2Coeffs,
			      hsg_driver->type1_weight_set_arrays(),
			      hsg_driver->type2_weight_set_arrays(),
			      numericalMoments);
  else
    compute_numerical_moments(expansionType1Coeffs,
			      hsg_driver->type1_weight_set_arrays(),
			      numericalMoments);
  */
}


void HierarchInterpPolyApproximation::
compute_numerical_expansion_moments(size_t num_moments)
{
  // for now: nested interpolation is exact
  expansionMoments = numericalMoments;

  // a couple different ways to go with this in the future:
  // (1) evaluate hierarchical value(lev) - value(lev-1) with HSGDriver wts
  // (2) evaluate value() with CSGDriver wts
  //  > promote Nodal implementation of this function to base class
  //  > redefine HierarchSparseGridDriver::type1_weight_sets() to generate
  //    from 1D weights array in CSG-style approach (not simple concatenation)
}


/** Computes the variance of component functions. Assumes that all
    subsets of set_value have been computed in advance which will be
    true so long as the partial_variance is called following
    appropriate enumeration of set value  */
void HierarchInterpPolyApproximation::compute_partial_variance(int set_value)
{
  Real& variance = partialVariance[sobolIndexMap[set_value]];
  // Computes the integral first

  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&    sm_index = hsg_driver->smolyak_multi_index();
  const UShort4DArray&  colloc_key = hsg_driver->collocation_key();
  const Sizet3DArray& colloc_index = hsg_driver->collocation_indices();
  // Smolyak recursion of anisotropic tensor products
  size_t i, num_levels = sm_index.size();
  /*
  UShortArray quad_order;
  for (i=0; i<num_levels; ++i) {
    for (j; j<num_sets; ++j) {
      hsg_driver->level_to_order(sm_index[i][j], quad_order);
      variance += partial_variance_integral(set_value, quad_order, sm_index[i],
					    colloc_key[i], colloc_index[i]);
    }
  }
  */

  // manage constituentSets
  InterpPolyApproximation::compute_partial_variance(set_value);
}


void HierarchInterpPolyApproximation::compute_total_sobol_indices()
{
  const Real& m1 = numericalMoments[1]; // standardized, if not num exception
  Real total_variance = (m1 > 0.) ? m1*m1 : m1;
  int j, set_value;

  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&    sm_index = hsg_driver->smolyak_multi_index();
  const UShort4DArray&  colloc_key = hsg_driver->collocation_key();
  const Sizet3DArray& colloc_index = hsg_driver->collocation_indices();

    // Smolyak recursion of anisotropic tensor products
    size_t i, num_levels = sm_index.size();
    /*
    UShortArray quad_order;
    // iterate each variable 
    for (j=0; j<numVars; ++j) {
      set_value = (int)std::pow(2.,(int)numVars) - (int)std::pow(2.,j) - 1; 
      for (i=0; i<num_smolyak_indices; ++i) {
	hsg_driver->level_to_order(sm_index[i], quad_order);
	totalSobolIndices[j] += 
	  total_effects_integral(set_value, quad_order, sm_index[i],
				 colloc_key[i], colloc_index[i]);
      }
      totalSobolIndices[j] = std::abs(1. - totalSobolIndices[j]/total_variance);
    }
    */
}


void HierarchInterpPolyApproximation::
member_coefficients_weights(int set_value, const UShortArray& quad_order,
			    const UShortArray& lev_index,
			    const UShort2DArray& key,
			    const SizetArray& colloc_index,
			    RealVector& member_coeffs, RealVector& member_wts)
{
  // create member variable key and get number of expansion coeffs in
  // member-variable-only expansion
  BoolDeque nonmember_vars(numVars); // distinguish set members from non-members
  int num_member_coeffs = 1; // # exp coeffs in member-variable-only expansion
  IntVector indexing_factor(numVars, false); // factors indexing member vars 
  for (int k=0; k<numVars; ++k) {
    // if subset contains variable k, set key for variable k to true
    if (set_value & (1 << k)) {
      nonmember_vars[k]  = false;	
      indexing_factor[k] = num_member_coeffs; // for indexing of member_coeffs
      num_member_coeffs *= quad_order[k];
    }
    else {
      nonmember_vars[k]  = true;	
      indexing_factor[k] = 1;
    }
  }

  // Size vectors to store new coefficients
  member_coeffs.sizeUninitialized(num_member_coeffs);
  member_wts.sizeUninitialized(num_member_coeffs);

  // Perform integration over non-member variables and store indices
  // of new expansion
  size_t i, j, num_colloc_pts = key.size();
  const Real3DArray& colloc_wts_1d
    = driverRep->type1_collocation_weights_array();
  /*
  for (i=0; i <num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    size_t member_coeffs_index = 0;	
    Real prod_i_nonmembers = 1., prod_i_members = 1.;
    for (j=0; j<numVars; ++j)
      // Save the product of the weights of the member and non-member variables 
      if (nonmember_vars[j])
	prod_i_nonmembers *= colloc_wts_1d[lev_index[j]][j][key_i[j]];
      else {
	// Convert key to corresponding index on member_coeffs
	member_coeffs_index += key_i[j]*indexing_factor[j];
	prod_i_members      *= colloc_wts_1d[lev_index[j]][j][key_i[j]];
      }

    // member_wts is performed more time than necessary here, but it
    // seems to be the simplest place to put it
    member_wts[member_coeffs_index] = prod_i_members;
    // sort coefficients by the "signature" of the member variables
    // (i.e. member_coeffs_index)
    unsigned short c_index = (colloc_index.empty()) ? i : colloc_index[i];
    member_coeffs[member_coeffs_index]
      += expansionType1Coeffs[c_index]*prod_i_nonmembers;
  }
  */
}

}
