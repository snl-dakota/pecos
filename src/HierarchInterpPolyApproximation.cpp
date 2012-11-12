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
#include "pecos_stat_util.hpp"

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

  if (expConfigOptions.refinementControl) {
    size_t num_moments = (nonRandomIndices.empty()) ? 4 : 2;
    if (referenceMoments.empty())
      referenceMoments.sizeUninitialized(num_moments);
    if (deltaMoments.empty())
      deltaMoments.sizeUninitialized(num_moments);
  }
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

  HierarchSparseGridDriver* hsg_driver   = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      key          = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, v, num_levels = key.size(), num_sets, num_tp_pts,
    cntr = 0, index, num_deriv_vars = surrData.num_derivative_variables();

  // level 0
  index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
  if (expConfigOptions.expansionCoeffFlag) {
    expansionType1Coeffs[0][0][0] = surrData.response_function(index);
    if (basisConfigOptions.useDerivs)
      Teuchos::setCol(surrData.response_gradient(index), 0,
		      expansionType2Coeffs[0][0]);
  }
  if (expConfigOptions.expansionCoeffGradFlag)
    Teuchos::setCol(surrData.response_gradient(index), 0,
		    expansionType1CoeffGrads[0][0]);
  ++cntr;
  // levels 1 to num_levels
  for (lev=1; lev<num_levels; ++lev) {
    const UShort3DArray& key_l = key[lev];
    num_sets = key_l.size();
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = key_l[set].size();
      for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	const RealVector& c_vars = surrData.continuous_variables(index);
	// coefficients are hierarchical surpluses
	if (expConfigOptions.expansionCoeffFlag) {
	  expansionType1Coeffs[lev][set][pt] = surrData.response_function(index)
	    - value(c_vars, sm_mi, key, expansionType1Coeffs,
		    expansionType2Coeffs, lev-1);
	  if (basisConfigOptions.useDerivs) {
	    const RealVector& data_grad = surrData.response_gradient(index);
	    const RealVector& prev_grad = gradient_basis_variables(c_vars,
	      sm_mi, key, expansionType1Coeffs, expansionType2Coeffs, lev-1);
	    Real* hier_grad = expansionType2Coeffs[lev][set][pt];
	    for (v=0; v<num_deriv_vars; ++v)
	      hier_grad[v] = data_grad[v] - prev_grad[v];
	  }
	}
	if (expConfigOptions.expansionCoeffGradFlag) {
	  const RealVector& data_grad = surrData.response_gradient(index);
	  const RealVector& prev_grad = gradient_nonbasis_variables(c_vars,
	    sm_mi, key, expansionType1CoeffGrads, lev-1);
	  Real* hier_grad = expansionType1CoeffGrads[lev][set][pt];
	  for (v=0; v<num_deriv_vars; ++v)
	    hier_grad[v] = data_grad[v] - prev_grad[v];
	}
      }
    }
  }

  computedMean = computedVariance
    = computedRefMean = computedDeltaMean
    = computedRefVariance = computedDeltaVariance = 0;
}


void HierarchInterpPolyApproximation::increment_expansion_coefficients()
{
  increment_current_from_reference();

  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  switch (expConfigOptions.refinementControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: // generalized sparse grids
    increment_expansion_coefficients(hsg_driver->trial_set());
    break;
  default: {
    const UShort3DArray&   sm_mi = hsg_driver->smolyak_multi_index();
    const UShortArray& incr_sets = hsg_driver->increment_sets();
    size_t lev, num_lev = sm_mi.size(), set, start_set, num_sets;
    for (lev=0; lev<num_lev; ++lev) {
      start_set = incr_sets[lev]; num_sets = sm_mi[lev].size();
      for (set=start_set; set<num_sets; ++set)
	increment_expansion_coefficients(sm_mi[lev][set]);
    }
    break;
  }
  }
}


// ******************************************************************
// TO DO: verify that decrement/restore is always valid for surpluses
// ******************************************************************


void HierarchInterpPolyApproximation::decrement_expansion_coefficients()
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShortArray&        trial_set  = hsg_driver->trial_set();
  size_t lev = hsg_driver->index_norm(trial_set);

  if (expConfigOptions.expansionCoeffFlag) {
    savedExpT1Coeffs[trial_set] = expansionType1Coeffs[lev].back();
    expansionType1Coeffs[lev].pop_back();
    if (basisConfigOptions.useDerivs) {
      savedExpT2Coeffs[trial_set] = expansionType2Coeffs[lev].back();
      expansionType2Coeffs[lev].pop_back();
    }
  }
  if (expConfigOptions.expansionCoeffGradFlag) {
    savedExpT1CoeffGrads[trial_set] = expansionType1CoeffGrads[lev].back();
    expansionType1CoeffGrads[lev].pop_back();
  }

  decrement_current_to_reference();
}


void HierarchInterpPolyApproximation::finalize_expansion_coefficients()
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&      sm_mi      = hsg_driver->smolyak_multi_index();

  size_t lev, set, num_sets, num_levels = sm_mi.size(), num_smolyak_sets,
    num_coeff_sets;
  for (lev=0; lev<num_levels; ++lev) {
    const UShort2DArray& sm_mi_l = sm_mi[lev];
    num_smolyak_sets = sm_mi_l.size();
    num_coeff_sets = (expConfigOptions.expansionCoeffFlag) ?
      expansionType1Coeffs[lev].size() : expansionType1CoeffGrads[lev].size();
    for (set=num_coeff_sets; set<num_smolyak_sets; ++set)
      restore_expansion_coefficients(sm_mi_l[set]);
  }
  savedExpT1Coeffs.clear(); savedExpT2Coeffs.clear();
  savedExpT1CoeffGrads.clear();

  computedMean = computedVariance = 0;
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

  computedMean = computedVariance = 0;
}


/** Lower level helper function to process a single index set. */
void HierarchInterpPolyApproximation::
increment_expansion_coefficients(const UShortArray& index_set)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  size_t lev = hsg_driver->index_norm(index_set);

  if (lev >= expansionType1Coeffs.size()) {
    expansionType1Coeffs.resize(lev+1);
    expansionType2Coeffs.resize(lev+1);
    expansionType1CoeffGrads.resize(lev+1);
  }
  size_t set = expansionType1Coeffs[lev].size();
  // append empty and update in place
  RealVector fns; RealMatrix grads;
  expansionType1Coeffs[lev].push_back(fns);
  expansionType2Coeffs[lev].push_back(grads);
  expansionType1CoeffGrads[lev].push_back(grads);
  RealVector& t1_coeffs      = expansionType1Coeffs[lev][set];
  RealMatrix& t2_coeffs      = expansionType2Coeffs[lev][set];
  RealMatrix& t1_coeff_grads = expansionType1CoeffGrads[lev][set];

  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  const UShort4DArray& key   = hsg_driver->collocation_key();
  size_t index, pt, num_trial_pts = key[lev][set].size(), v, num_deriv_vars = 0;
  if (expConfigOptions.expansionCoeffFlag) {
    t1_coeffs.sizeUninitialized(num_trial_pts);
    if (basisConfigOptions.useDerivs) {
      num_deriv_vars = expansionType2Coeffs[0][0].numRows();
      t2_coeffs.shapeUninitialized(num_deriv_vars, num_trial_pts);
    }
  }
  if (expConfigOptions.expansionCoeffGradFlag) {
    num_deriv_vars = expansionType1CoeffGrads[0][0].numRows();
    t1_coeff_grads.shapeUninitialized(num_deriv_vars, num_trial_pts);
  }
 
  for (pt=0, index=numCollocPts; pt<num_trial_pts; ++pt, ++index) {
    const RealVector& c_vars = surrData.continuous_variables(index);
    if (expConfigOptions.expansionCoeffFlag) {
      t1_coeffs[pt] = surrData.response_function(index) - value(c_vars, sm_mi,
	key, expansionType1Coeffs, expansionType2Coeffs, lev-1);
      if (basisConfigOptions.useDerivs) {
	const RealVector& data_grad = surrData.response_gradient(index);
	const RealVector& prev_grad = gradient_basis_variables(c_vars, sm_mi,
	  key, expansionType1Coeffs, expansionType2Coeffs, lev-1);
	Real* hier_grad = t2_coeffs[pt];
	for (v=0; v<num_deriv_vars; ++v)
	  hier_grad[v] = data_grad[v] - prev_grad[v];
      }
    }
    if (expConfigOptions.expansionCoeffGradFlag) {
      const RealVector& data_grad = surrData.response_gradient(index);
      const RealVector& prev_grad = gradient_nonbasis_variables(c_vars, sm_mi,
	key, expansionType1CoeffGrads, lev-1);
      Real* hier_grad = t1_coeff_grads[pt];
      for (v=0; v<num_deriv_vars; ++v)
	hier_grad[v] = data_grad[v] - prev_grad[v];
    }
  }
  numCollocPts += num_trial_pts;
}


/** Lower level helper function to process a single index set. */
void HierarchInterpPolyApproximation::
restore_expansion_coefficients(const UShortArray& restore_set)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  size_t lev = hsg_driver->index_norm(restore_set);
  if (expConfigOptions.expansionCoeffFlag) {
    expansionType1Coeffs[lev].push_back(savedExpT1Coeffs[restore_set]);
    savedExpT1Coeffs.erase(restore_set);
    if (basisConfigOptions.useDerivs) {
      expansionType2Coeffs[lev].push_back(savedExpT2Coeffs[restore_set]);
      savedExpT2Coeffs.erase(restore_set);
    }
  }
  if (expConfigOptions.expansionCoeffGradFlag) {
    expansionType1CoeffGrads[lev].push_back(savedExpT1CoeffGrads[restore_set]);
    savedExpT1CoeffGrads.erase(restore_set);
  }
}


void HierarchInterpPolyApproximation::increment_current_from_reference()
{
  computedRefMean     = computedMean;
  computedRefVariance = computedVariance;

  if ( (computedMean & 1) || (computedVariance & 1) )
    referenceMoments = numericalMoments;
  if (computedMean & 2)
    meanRefGradient = meanGradient;
  if (computedVariance & 2)
    varianceRefGradient = varianceGradient;

  // clear current and delta
  computedMean = computedVariance =
    computedDeltaMean = computedDeltaVariance = 0;
}


void HierarchInterpPolyApproximation::decrement_current_to_reference()
{
  computedMean     = computedRefMean;
  computedVariance = computedRefVariance;

  if ( (computedRefMean & 1) || (computedRefVariance & 1) )
    numericalMoments = referenceMoments;
  if (computedRefMean & 2)
    meanGradient = meanRefGradient;
  if (computedRefVariance & 2)
    varianceGradient = varianceRefGradient;

  // leave reference settings, but clear delta settings
  computedDeltaMean = computedDeltaVariance = 0;
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
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  size_t lev, set, num_sets;
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    // sm_mi and key include all current index sets, whereas the t1/t2
    // coeffs may refect a partial state derived from a reference_key
    num_sets = t1_coeffs_l.size();
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

  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  size_t lev, set, num_sets;
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    // sm_mi and key include all current index sets, whereas the t1/t2
    // coeffs may refect a partial state derived from a reference_key
    num_sets = t1_coeffs_l.size();
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

  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    // sm_mi and key include all current index sets, whereas the t1/t2
    // coeffs may refect a partial state derived from a reference_key
    num_sets = t1_coeffs_l.size();
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

  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&            sm_mi_l = sm_mi[lev];
    const UShort3DArray&              key_l = key[lev];
    const RealMatrixArray& t1_coeff_grads_l = t1_coeff_grads[lev];
    // sm_mi and key include all current index sets, whereas the t1/t2
    // coeffs may refect a partial state derived from a reference_key
    num_sets = t1_coeff_grads_l.size();
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

  Real& mean = numericalMoments[0];
  if ( !(computedMean & 1) ) {
    mean = expectation(expansionType1Coeffs, expansionType2Coeffs);
    computedMean |= 1;
  }
  return mean;
}



Real HierarchInterpPolyApproximation::mean(const RealVector& x)
{
  Real& mean = numericalMoments[0];
  if ( !(computedMean & 1) || !match_nonrandom_vars(x, xPrevMean) ) {

    // TO DO
    PCerr << "\nError: mean not yet implemented for all variables mode."
	  << std::endl;
    abort_handler(-1);

    computedMean |= 1; xPrevMean = x;
  }
  return mean;
}


const RealVector& HierarchInterpPolyApproximation::mean_gradient()
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Hierarch"
	  << "InterpPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  if ( !(computedMean & 2) ) {

    HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
    const RealVector2DArray& t1_wts = hsg_driver->type1_weight_set_arrays();

    // TO DO
    PCerr << "\nError: mean_gradient not yet implemented." << std::endl;
    abort_handler(-1);
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
    computedMean |= 2;
  }
  return meanGradient;
}


const RealVector& HierarchInterpPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  if ( !(computedMean & 2) || !match_nonrandom_vars(x, xPrevMeanGrad) ) {
    // TO DO
    PCerr << "\nError: mean_gradient not yet implemented for all variables "
	  << "mode." << std::endl;
    abort_handler(-1);

    computedMean |= 2; xPrevMeanGrad = x;
  }
  return meanGradient;
}


Real HierarchInterpPolyApproximation::variance()
{
  Real& var = numericalMoments[1];
  if ( !(computedVariance & 1) ) {
    var = covariance(this);
    computedVariance |= 1;
  }
  return var;
}


Real HierarchInterpPolyApproximation::variance(const RealVector& x)
{
  Real& var = numericalMoments[1];
  if ( !(computedVariance & 1) || !match_nonrandom_vars(x, xPrevVar) ) {
    var = covariance(x, this);
    computedVariance |= 1; xPrevVar = x;
  }
  return var;
}


const RealVector& HierarchInterpPolyApproximation::variance_gradient()
{
  if ( !(computedVariance & 2) ) {
    // TO DO
    PCerr << "\nError: variance_gradient not yet implemented." << std::endl;
    abort_handler(-1);

    computedVariance |= 2;
  }
  return varianceGradient;
}


const RealVector& HierarchInterpPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  if ( !(computedVariance & 2) || !match_nonrandom_vars(x, xPrevVarGrad) ) {
    // TO DO
    PCerr << "\nError: variance_gradient not yet implemented for all variables "
	  << "mode." << std::endl;
    abort_handler(-1);

    computedVariance |= 2; xPrevVarGrad = x;
  }
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

  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  bool same = (this == hip_approx_2);
  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  Real mean_1 = mean(), mean_2 = (same) ? mean_1 : hip_approx_2->mean();
  central_product_interpolant(hip_approx_2, mean_1, mean_2,
			      cov_t1_coeffs, cov_t2_coeffs);

  // evaluate expectation of these t1/t2 coefficients
  Real covar = expectation(cov_t1_coeffs, cov_t2_coeffs);
  // Note: separation of reference and increment using cov_t{1,2}_coeffs
  // with {ref,incr}_key would provide an increment of a central moment
  // around an invariant center.  For hierarchical covariance, one must
  // also account for the change in mean as in delta_covariance().

  if (same)
    { numericalMoments[1] = covar; computedVariance |= 1; }
  return covar;
}


Real HierarchInterpPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  // TO DO
  PCerr << "\nError: covariance not yet implemented for all variables mode."
	<< std::endl;
  abort_handler(-1);

  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  bool same = (this == hip_approx_2);
  Real covar = 0.;

  if (same)
    { numericalMoments[1] = covar; computedVariance |= 1; }
  return covar;
}


Real HierarchInterpPolyApproximation::reference_mean()
{
  if ( !(computedRefMean & 1) ) {
    HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
    UShort2DArray ref_key, incr_key;
    hsg_driver->partition_keys(ref_key, incr_key);
    return reference_mean(ref_key);
  }
  else
    return referenceMoments[0];
}


Real HierarchInterpPolyApproximation::
reference_mean(const UShort2DArray& ref_key)
{
  Real& ref_mean = referenceMoments[0];
  if ( !(computedRefMean & 1) ) {
    ref_mean = expectation(expansionType1Coeffs, expansionType2Coeffs, ref_key);
    computedRefMean |= 1;
  }
  return ref_mean;
}


Real HierarchInterpPolyApproximation::reference_variance()
{
  if ( !(computedRefVariance & 1) ) {
    HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
    UShort2DArray ref_key, incr_key;
    hsg_driver->partition_keys(ref_key, incr_key);
    return reference_variance(ref_key);
  }
  else
    return referenceMoments[1];
}


Real HierarchInterpPolyApproximation::
reference_variance(const UShort2DArray& ref_key)
{
  Real& ref_var = referenceMoments[1];
  if ( !(computedRefVariance & 1) ) {
    Real ref_mean = reference_mean(ref_key);
    RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
    central_product_interpolant(this, ref_mean, ref_mean, cov_t1_coeffs,
				cov_t2_coeffs, ref_key);
    ref_var = expectation(cov_t1_coeffs, cov_t2_coeffs, ref_key),
    computedRefVariance |= 1;
  }
  return ref_var;
}


Real HierarchInterpPolyApproximation::delta_mean()
{
  if ( !(computedDeltaMean & 1) ) {
    HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
    UShort2DArray ref_key, incr_key;
    hsg_driver->partition_keys(ref_key, incr_key);
    return delta_mean(incr_key);
  }
  else
    return deltaMoments[0];
}


Real HierarchInterpPolyApproximation::delta_mean(const UShort2DArray& incr_key)
{
  Real& delta_var = deltaMoments[0];
  if ( !(computedDeltaMean & 1) ) {
    delta_var
      = expectation(expansionType1Coeffs, expansionType2Coeffs, incr_key);
    computedDeltaMean |= 1;
  }
  return delta_var;
}


Real HierarchInterpPolyApproximation::delta_variance()
{
  if ( !(computedDeltaVariance & 1) ) {
    HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
    UShort2DArray ref_key, incr_key;
    hsg_driver->partition_keys(ref_key, incr_key);
    return delta_variance(ref_key, incr_key);
  }
  else
    return deltaMoments[1];
}


Real HierarchInterpPolyApproximation::
delta_variance(const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  Real& delta_var = deltaMoments[1];
  if ( !(computedDeltaVariance & 1) ) {
    delta_var = delta_covariance(this, ref_key, incr_key),
    computedDeltaVariance |= 1;
  }
  return delta_var;
}


Real HierarchInterpPolyApproximation::delta_std_deviation()
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  UShort2DArray ref_key, incr_key;
  hsg_driver->partition_keys(ref_key, incr_key);

  return delta_std_deviation(ref_key, incr_key);
}


Real HierarchInterpPolyApproximation::
delta_std_deviation(const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  // delta-sigma = sqrt( var0 + delta-var ) - sigma0
  //             = [ sqrt(1 + delta_var / var0) - 1 ] * sigma0
  //             = sqrt1pm1(delta_var / var0) * sigma0
  // where sqrt1pm1(x) = expm1[ log1p(x) / 2 ]

  Real delta_var = delta_variance(ref_key, incr_key),
       var0      = reference_variance(ref_key),
       sigma0    = std::sqrt(var0);

  return (delta_var < var0) ?
    bmth::sqrt1pm1(delta_var / var0) * sigma0 :            // preserve precision
    std::sqrt(var0 + delta_var) - sigma0; // precision OK; prevent div by var0=0
}


Real HierarchInterpPolyApproximation::delta_beta(bool cdf_flag, Real z_bar)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  UShort2DArray ref_key, incr_key;
  hsg_driver->partition_keys(ref_key, incr_key);

  return delta_beta(cdf_flag, z_bar, ref_key, incr_key);
}


Real HierarchInterpPolyApproximation::
delta_beta(bool cdf_flag, Real z_bar, const UShort2DArray& ref_key,
	   const UShort2DArray& incr_key)
{
  //  CDF delta-beta = (mu1 - z-bar)/sigma1 - (mu0 - z-bar)/sigma0
  //    = (mu1 sigma0 - z-bar sigma0 - mu0 sigma1 + z-bar sigma1)/sigma1/sigma0
  //    = (delta-mu sigma0 - mu0 delta-sigma + z-bar delta-sigma)/sigma1/sigma0
  //    = (delta-mu - delta-sigma beta0)/sigma1
  // CCDF delta-beta = (z-bar - mu1)/sigma1 - (z-bar - mu0)/sigma0
  //    = (z-bar sigma0 - mu1 sigma0 - z-bar sigma1 + mu0 sigma1)/sigma1/sigma0
  //    = (mu0 delta-sigma - delta-mu sigma0 - z-bar delta-sigma)/sigma1/sigma0
  //    = -delta-mu/sigma1 - delta_sigma (z-bar - mu0) / sigma0 / sigma1
  //    = (-delta-mu - delta-sigma beta0)/sigma1

  Real beta0, mu0 = reference_mean(ref_key), delta_mu = delta_mean(incr_key),
    var0 = reference_variance(ref_key), sigma0 = std::sqrt(var0),
    delta_sigma = delta_std_deviation(ref_key, incr_key),
    sigma1 = sigma0 + delta_sigma;

  // Error traps are needed for zero variance: a single point ref grid
  // (level=0 sparse or m=1 tensor) has zero variance.  Unchanged response
  // values along an index set could then cause sigma1 also = 0.
  if (cdf_flag) {
    if (sigma0 > SMALL_NUMBER && sigma1 > SMALL_NUMBER) {
      beta0 = (mu0 - z_bar) / sigma0;
      return ( delta_mu - delta_sigma * beta0) / sigma1;
    }
    else if (sigma1 > SMALL_NUMBER)// neglect beta0 term (zero init reliability)
      return delta_mu / sigma1; // or delta = beta1 = (mu1 - z_bar) / sigma1 ?
    else if (sigma0 > SMALL_NUMBER) // assume beta1 = 0 -> delta = -beta0
      return (z_bar - mu0) / sigma0;
    else                      // assume beta0 = beta1 = 0
      return 0;
  }
  else {
    if (sigma0 > SMALL_NUMBER && sigma1 > SMALL_NUMBER) {
      beta0 = (z_bar - mu0) / sigma0;
      return (-delta_mu - delta_sigma * beta0) / sigma1;
    }
    else if (sigma1 > SMALL_NUMBER)// neglect beta0 term (zero init reliability)
      return -delta_mu / sigma1;
    else if (sigma0 > SMALL_NUMBER) // assume beta1 = 0 -> delta = -beta0
      return (mu0 - z_bar) / sigma0;
    else                      // assume beta0 = beta1 = 0
      return 0;
  }
}


Real HierarchInterpPolyApproximation::delta_z(bool cdf_flag, Real beta_bar)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  UShort2DArray ref_key, incr_key;
  hsg_driver->partition_keys(ref_key, incr_key);

  return delta_z(cdf_flag, beta_bar, ref_key, incr_key);
}


Real HierarchInterpPolyApproximation::
delta_z(bool cdf_flag, Real beta_bar, const UShort2DArray& ref_key,
	const UShort2DArray& incr_key)
{
  //  CDF delta-z = (mu1 - sigma1 beta-bar) - (mu0 - sigma0 beta-bar)
  //              = delta-mu - delta-sigma * beta-bar
  // CCDF delta-z = (mu1 + sigma1 beta-bar) - (mu0 + sigma0 beta-bar)
  //              = delta-mu + delta-sigma * beta-bar

  Real delta_mu = delta_mean(incr_key),
    delta_sigma = delta_std_deviation(ref_key, incr_key);
  return (cdf_flag) ? delta_mu - delta_sigma * beta_bar :
                      delta_mu + delta_sigma * beta_bar;
}


Real HierarchInterpPolyApproximation::
delta_covariance(PolynomialApproximation* poly_approx_2)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  UShort2DArray ref_key, incr_key;
  hsg_driver->partition_keys(ref_key, incr_key);

  return delta_covariance(poly_approx_2, ref_key, incr_key);
}


Real HierarchInterpPolyApproximation::
delta_covariance(PolynomialApproximation* poly_approx_2,
		 const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::delta_covariance()" << std::endl;
    abort_handler(-1);
  }

  // Supports multiple grid increments in discerning nominal from delta based
  // on isotropic/anisotropic/generalized index set increments.  In current
  // use, 2D keys with set ranges are sufficient: level -> {start,end} set.
  // In the future, may need 3D keys for level/set/point.
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  bool same = (this == hip_approx_2);
  RealVector2DArray r1r2_t1_coeffs; RealMatrix2DArray r1r2_t2_coeffs;
  product_interpolant(hip_approx_2, r1r2_t1_coeffs, r1r2_t2_coeffs);

  // Compute surplus for r1, r2, and r1r2 and retrieve reference mean values
  Real ref_mean_r1 = reference_mean(ref_key),
    delta_mean_r1 = delta_mean(incr_key),
    ref_mean_r2 = (same) ? ref_mean_r1 :
      expectation(hip_approx_2->expansionType1Coeffs,
		  hip_approx_2->expansionType2Coeffs, ref_key),
    delta_mean_r2 = (same) ? delta_mean_r1 :
      expectation(hip_approx_2->expansionType1Coeffs,
		  hip_approx_2->expansionType2Coeffs, incr_key),
    delta_mean_r1r2 = expectation(r1r2_t1_coeffs, r1r2_t2_coeffs, incr_key);

  // Hierarchical increment to covariance:
  // \Delta\Sigma_ij = \Sigma^1_ij - \Sigma^0_ij
  //   = ( E[Ri Rj]^1 - E[Ri]^1 E[Rj]^1 ) - ( E[Ri Rj]^0 - E[Ri]^0 E[Rj]^0 )
  //   = E[Ri Rj]^0 + \DeltaE[Ri Rj]
  //     - (E[Ri]^0 + \DeltaE[Ri]) (E[Rj]^0 + \DeltaE[Rj])
  //     - E[Ri Rj]^0 + E[Ri]^0 E[Rj]^0
  //   = \DeltaE[Ri Rj] - \DeltaE[Ri] E[Rj]^0 - E[Ri]^0 \DeltaE[Rj]
  //     - \DeltaE[Ri] \DeltaE[Rj]
  Real delta_covar = delta_mean_r1r2 - ref_mean_r1 * delta_mean_r2
     - ref_mean_r2 * delta_mean_r1 - delta_mean_r1 * delta_mean_r2;
  if (same)
    { deltaMoments[1] = delta_covar; computedDeltaVariance |= 1; }
  return delta_covar;
}


Real HierarchInterpPolyApproximation::
delta_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  // TO DO
  PCerr << "\nError: delta_covariance not yet implemented for all variables "
	<< "mode." << std::endl;
  abort_handler(-1);

  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  bool same = (this == hip_approx_2);

  Real delta_covar = 0.;

  if (same)
    { deltaMoments[1] = delta_covar; computedDeltaVariance |= 1; }
  return delta_covar;
}


Real HierarchInterpPolyApproximation::
expectation(const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs,
	    const UShort2DArray& set_partition)
{
  Real integral = 0.;
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const RealVector2DArray& t1_wts = hsg_driver->type1_weight_set_arrays();
  size_t lev, set, pt, num_levels = t1_coeffs.size(), set_start = 0, set_end,
    num_tp_pts;
  bool partial = !set_partition.empty();
  switch (basisConfigOptions.useDerivs) {
  case false:
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1_coeffs_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) // omitted if empty surplus vector
	  integral += t1_coeffs_ls[pt] * t1_wts_ls[pt];
      }
    }
    break;
  case true: {
    const RealMatrix2DArray& t2_wts = hsg_driver->type2_weight_set_arrays();
    size_t v;
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1_coeffs_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealMatrix& t2_coeffs_ls = t2_coeffs[lev][set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	const RealMatrix&    t2_wts_ls = t2_wts[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
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


Real HierarchInterpPolyApproximation::
expectation(const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs,
	    const UShort3DArray& pt_partition)
{
  Real integral = 0.;
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const RealVector2DArray& t1_wts = hsg_driver->type1_weight_set_arrays();
  size_t lev, set, pt, num_levels = t1_coeffs.size(), num_sets,
    tp_pt_start = 0, tp_pt_end;
  bool partial = !pt_partition.empty();
  switch (basisConfigOptions.useDerivs) {
  case false:
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      num_sets = t1_coeffs_l.size();
      for (set=0; set<num_sets; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	if (partial) {
	  tp_pt_start = pt_partition[lev][set][0];
	  tp_pt_end   = pt_partition[lev][set][1];
	}
	else
	  tp_pt_end   = t1_coeffs_ls.length();
	for (pt=tp_pt_start; pt<tp_pt_end; ++pt)
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
	if (partial) {
	  tp_pt_start = pt_partition[lev][set][0];
	  tp_pt_end   = pt_partition[lev][set][1];
	}
	else
	  tp_pt_end   = t1_coeffs_ls.length();
	for (pt=tp_pt_start; pt<tp_pt_end; ++pt) {
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


/** Whereas expectation() supports either a reference or increment key
    (passed as generic set_partition), functions forming hierarchical
    interpolant coefficients support only a reference key (starting
    point must be set 0; end point can be controlled). */
void HierarchInterpPolyApproximation::
product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
		    RealVector2DArray& r1r2_t1_coeffs,
		    RealMatrix2DArray& r1r2_t2_coeffs,
		    const UShort2DArray& reference_key)
{
  HierarchSparseGridDriver* hsg_driver   = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      key          = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, num_levels = expansionType1Coeffs.size(),
    num_sets, num_tp_pts, cntr = 0, index;
  bool partial = !reference_key.empty();
  const SurrogateData& s_data_2 = hip_approx_2->surrData;

  // form hierarchical t1/t2 coeffs for raw moment R1 R2
  r1r2_t1_coeffs.resize(num_levels); r1r2_t1_coeffs[0].resize(1);
  r1r2_t2_coeffs.resize(num_levels); r1r2_t2_coeffs[0].resize(1);
  r1r2_t1_coeffs[0][0].sizeUninitialized(1);
  switch (basisConfigOptions.useDerivs) {
  case false:
    // level 0 (assume this is always contained in a partial reference_key)
    index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
    r1r2_t1_coeffs[0][0][0]
      = surrData.response_function(index) * s_data_2.response_function(index);
    ++cntr;
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? reference_key[lev][1] : key[lev].size();
      r1r2_t1_coeffs[lev].resize(num_sets);
      r1r2_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = key[lev][set].size();
	RealVector& r1r2_t1_coeffs_ls = r1r2_t1_coeffs[lev][set];
	r1r2_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of R1 R2
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  r1r2_t1_coeffs_ls[pt] = surrData.response_function(index)
	    * s_data_2.response_function(index)
	    - value(surrData.continuous_variables(index), sm_mi, key,
		    r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1);
	}
      }
    }
    break;
  case true:
    size_t v;
    // level 0 (assume this is always contained in a partial reference_key)
    index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
    Real data_fn1 = surrData.response_function(index);
    Real data_fn2 = s_data_2.response_function(index);
    r1r2_t1_coeffs[0][0][0] = data_fn1 * data_fn2;
    r1r2_t2_coeffs[0][0].shapeUninitialized(numVars, 1);
    Real *r1r2_t2_coeffs_000 = r1r2_t2_coeffs[0][0][0];
    const RealVector& data_grad1 = surrData.response_gradient(index);
    const RealVector& data_grad2 = s_data_2.response_gradient(index);
    for (v=0; v<numVars; ++v)
      r1r2_t2_coeffs_000[v]
	= data_fn1 * data_grad2[v] + data_fn2 * data_grad1[v];
    ++cntr;
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? reference_key[lev][1] : key[lev].size();
      r1r2_t1_coeffs[lev].resize(num_sets);
      r1r2_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& r1r2_t1_coeffs_ls = r1r2_t1_coeffs[lev][set];
	RealMatrix& r1r2_t2_coeffs_ls = r1r2_t2_coeffs[lev][set];
	num_tp_pts = key[lev][set].size();
	r1r2_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	r1r2_t2_coeffs_ls.shapeUninitialized(numVars, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  const RealVector& c_vars = surrData.continuous_variables(index);
	  // type1 hierarchical interpolation of R1 R2
	  data_fn1 = surrData.response_function(index);
	  data_fn2 = s_data_2.response_function(index);
	  r1r2_t1_coeffs_ls[pt] = data_fn1 * data_fn2 -
	    value(c_vars, sm_mi, key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1);
	  // type2 hierarchical interpolation of R1 R2
	  // --> interpolated grads are R1 * R2' + R2 * R1'
	  Real* r1r2_t2_coeffs_lsp = r1r2_t2_coeffs_ls[pt];
	  const RealVector& data_grad1 = surrData.response_gradient(index);
	  const RealVector& data_grad2 = s_data_2.response_gradient(index);
	  const RealVector& prev_grad  = gradient_basis_variables(c_vars,
	    sm_mi, key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1);
	  for (v=0; v<numVars; ++v)
	    r1r2_t2_coeffs_lsp[v] = data_fn1 * data_grad2[v]
	      + data_fn2 * data_grad1[v] - prev_grad[v];
	}
      }
    }
    break;
  }
}


/** Whereas expectation() supports either a reference or increment key
    (passed as generic set_partition), functions forming hierarchical
    interpolant coefficients support only a reference key (starting
    point must be set 0; end point can be controlled). */
void HierarchInterpPolyApproximation::
central_product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
			    Real mean_1, Real mean_2,
			    RealVector2DArray& cov_t1_coeffs,
			    RealMatrix2DArray& cov_t2_coeffs,
			    const UShort2DArray& reference_key)
{
  // form hierarchical t1/t2 coeffs for (R_1 - \mu_1) (R_2 - \mu_2)
  HierarchSparseGridDriver* hsg_driver   = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      key          = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, num_levels = expansionType1Coeffs.size(),
    num_sets, num_tp_pts, cntr = 0, index;
  bool partial = !reference_key.empty();
  const SurrogateData& s_data_2 = hip_approx_2->surrData;

  cov_t1_coeffs.resize(num_levels); cov_t1_coeffs[0].resize(1);
  cov_t2_coeffs.resize(num_levels); cov_t2_coeffs[0].resize(1);
  cov_t1_coeffs[0][0].sizeUninitialized(1);
  switch (basisConfigOptions.useDerivs) {
  case false:
    // level 0 (assume this is always contained in a partial reference_key)
    index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
    cov_t1_coeffs[0][0][0] = (surrData.response_function(index) - mean_1) *
                             (s_data_2.response_function(index) - mean_2);
    ++cntr;
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? reference_key[lev][1] : key[lev].size();
      cov_t1_coeffs[lev].resize(num_sets); cov_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = key[lev][set].size();
	RealVector& cov_t1_coeffs_ls = cov_t1_coeffs[lev][set];
	cov_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  cov_t1_coeffs_ls[pt]
	    = (surrData.response_function(index) - mean_1)
	    * (s_data_2.response_function(index) - mean_2)
	    - value(surrData.continuous_variables(index), sm_mi, key,
		    cov_t1_coeffs, cov_t2_coeffs, lev-1);
	}
      }
    }
    break;
  case true:
    size_t v;
    // level 0 (assume this is always contained in a partial reference_key)
    index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
    Real data_fn1_mm1 = surrData.response_function(index) - mean_1;
    Real data_fn2_mm2 = s_data_2.response_function(index) - mean_2;
    cov_t1_coeffs[0][0][0] = data_fn1_mm1 * data_fn2_mm2;
    cov_t2_coeffs[0][0].shapeUninitialized(numVars, 1);
    Real *cov_t2_coeffs_000 = cov_t2_coeffs[0][0][0];
    const RealVector& data_grad1 = surrData.response_gradient(index);
    const RealVector& data_grad2 = s_data_2.response_gradient(index);
    for (v=0; v<numVars; ++v)
      cov_t2_coeffs_000[v]
	= data_fn1_mm1 * data_grad2[v] + data_fn2_mm2 * data_grad1[v];
    ++cntr;
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? reference_key[lev][1] : key[lev].size();
      cov_t1_coeffs[lev].resize(num_sets); cov_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& cov_t1_coeffs_ls = cov_t1_coeffs[lev][set];
	RealMatrix& cov_t2_coeffs_ls = cov_t2_coeffs[lev][set];
	num_tp_pts = key[lev][set].size();
	cov_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	cov_t2_coeffs_ls.shapeUninitialized(numVars, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  const RealVector& c_vars = surrData.continuous_variables(index);
	  // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  data_fn1_mm1 = surrData.response_function(index) - mean_1;
	  data_fn2_mm2 = s_data_2.response_function(index) - mean_2;
	  cov_t1_coeffs_ls[pt] = data_fn1_mm1 * data_fn2_mm2 -
	    value(c_vars, sm_mi, key, cov_t1_coeffs, cov_t2_coeffs, lev-1);
	  // type2 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  // --> interpolated grads are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
	  Real* cov_t2_coeffs_lsp = cov_t2_coeffs_ls[pt];
	  const RealVector& data_grad1 = surrData.response_gradient(index);
	  const RealVector& data_grad2 = s_data_2.response_gradient(index);
	  const RealVector& prev_grad  = gradient_basis_variables(c_vars,
	    sm_mi, key, cov_t1_coeffs, cov_t2_coeffs, lev-1);
	  for (v=0; v<numVars; ++v)
	    cov_t2_coeffs_lsp[v] = data_fn1_mm1 * data_grad2[v]
	      + data_fn2_mm2 * data_grad1[v] - prev_grad[v];
	}
      }
    }
    break;
  }
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

  HierarchSparseGridDriver* hsg_driver   = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      key          = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, num_levels = key.size(), num_sets, num_tp_pts,
    cntr, index;
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
      index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
      mom_t1_coeffs[0][0][0]
	= std::pow(surrData.response_function(index) - mean, moment);
      ++cntr;
      // levels 1:w
      for (lev=1; lev<num_levels; ++lev) {
	num_sets = key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  RealVector& mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  num_tp_pts = key[lev][set].size();
	  // type1 hierarchical interpolation of (R - \mu)^moment
	  for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	    index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	    mom_t1_coeffs_ls[pt]
	      = std::pow(surrData.response_function(index) - mean, moment)
	      - value(surrData.continuous_variables(index), sm_mi, key,
		      mom_t1_coeffs, mom_t2_coeffs, lev-1);
	  }
	}
      }
      break;
    case true:
      size_t v;
      // level 0
      index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
      Real data_fn_mm         = surrData.response_function(index) - mean;
      mom_t1_coeffs[0][0][0]  = std::pow(data_fn_mm, moment);
      Real* mom_t2_coeffs_000 = mom_t2_coeffs[0][0][0];
      Real deriv = moment * std::pow(data_fn_mm, m_index);
      const RealVector& data_grad = surrData.response_gradient(index);
      for (v=0; v<numVars; ++v)
	mom_t2_coeffs_000[v] = deriv * data_grad[v];
      ++cntr;
      // levels 1:w
      for (lev=1; lev<num_levels; ++lev) {
	num_sets = key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  num_tp_pts = key[lev][set].size();
	  RealVector& mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  RealMatrix& mom_t2_coeffs_ls = mom_t2_coeffs[lev][set];
	  for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	    index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	    const RealVector& c_vars = surrData.continuous_variables(index);
	    // type1 hierarchical interpolation of (R - \mu)^moment
	    data_fn_mm = surrData.response_function(index) - mean;
	    mom_t1_coeffs_ls[pt] = std::pow(data_fn_mm, moment) -
	      value(c_vars, sm_mi, key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	    // type2 hierarchical interpolation of (R - \mu)^moment
	    // --> interpolated grads are moment(R-\mu)^{moment-1} R'
	    Real* mom_t2_coeffs_lsp = mom_t2_coeffs_ls[pt];
	    deriv = moment * std::pow(data_fn_mm, m_index);
	    const RealVector& data_grad = surrData.response_gradient(index);
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
  //standardize_moments(numericalMoments);

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


/** Computes the variance of component functions.  Assumes that all
    subsets of set_value have been computed in advance which will be
    true so long as the partial_variance is called following
    appropriate enumeration of set value. */
void HierarchInterpPolyApproximation::
compute_partial_variance(const BitArray& set_value)
{
  Real& variance = partialVariance[sobolIndexMap[set_value]];

  // Compute the partial integral corresponding to set_value
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&        sm_index = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key = hsg_driver->collocation_key();
  const Sizet3DArray&     colloc_index = hsg_driver->collocation_indices();

  // Smolyak recursion of anisotropic tensor products
  size_t lev, set, num_levels = colloc_key.size(), num_sets;
  UShortArray quad_order; SizetArray empty_array;
  variance = 0.;
  for (lev=0; lev<num_levels; ++lev) {
    num_sets = colloc_key[lev].size();
    for (set=0; set<num_sets; ++set) {
      hsg_driver->level_to_order(sm_index[lev][set], quad_order);
      const SizetArray& colloc_index_ls = (colloc_index.empty()) ?
	empty_array : colloc_index[lev][set];
      variance +=
	partial_variance_integral(set_value, quad_order, sm_index[lev][set],
				  colloc_key[lev][set], colloc_index_ls,
				  expansionType1Coeffs[lev][set],
				  expansionType2Coeffs[lev][set]);
    }
  }

  // compute proper subsets and subtract their contributions
  InterpPolyApproximation::compute_partial_variance(set_value);
}


void HierarchInterpPolyApproximation::compute_total_sobol_indices()
{
  // Compute the total expansion variance.  For standard mode, the full variance
  // is likely already available, as managed by computedVariance in variance().
  // For all variables mode, we use covariance(this) without passing x for the
  // nonRandomIndices (bypass computedVariance checks by not using variance()).
  Real total_variance = (nonRandomIndices.empty()) ? variance() : // std mode
                        covariance(this);                    // all vars mode

  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray&        sm_index = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key = hsg_driver->collocation_key();
  const Sizet3DArray&     colloc_index = hsg_driver->collocation_indices();

  // Smolyak recursion of anisotropic tensor products
  size_t v, lev, set, num_levels = colloc_key.size(), num_sets;
  UShortArray quad_order; SizetArray empty_array; Real complement_variance;
  // iterate each variable 
  BitArray set_value(numVars);
  for (v=0; v<numVars; ++v) {
    // define set_value that includes all but index of interest
    set_value.set(); set_value[v].flip();
    complement_variance = 0.;
    for (lev=0; lev<num_levels; ++lev) {
      num_sets = colloc_key[lev].size();
      for (set=0; set<num_sets; ++set) {
	hsg_driver->level_to_order(sm_index[lev][set], quad_order);
	const SizetArray& colloc_index_ls = (colloc_index.empty()) ?
	  empty_array : colloc_index[lev][set];
	complement_variance +=
	  total_effects_integral(set_value, quad_order, sm_index[lev][set],
				 colloc_key[lev][set], colloc_index_ls,
				 expansionType1Coeffs[lev][set],
				 expansionType2Coeffs[lev][set]);
      }
    }
    totalSobolIndices[v] = std::abs(1. - complement_variance / total_variance);
  }
}

}
