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
//#define INTERPOLATION_TEST

namespace Pecos {


void HierarchInterpPolyApproximation::allocate_expansion_coefficients()
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort4DArray& key = hsg_driver->collocation_key();
  size_t i, j, k, num_levels = key.size(), num_sets, num_tp_pts,
    num_deriv_vars = surrData.num_derivative_variables();
  bool t1c = false, t2c = false, t1g = false;
  if (expConfigOptions.expansionCoeffFlag) {
    if (expansionType1Coeffs.size() != num_levels)
      { expansionType1Coeffs.resize(num_levels); t1c = true; }
    if (basisConfigOptions.useDerivs &&
	expansionType2Coeffs.size() != num_levels)
      { expansionType2Coeffs.resize(num_levels); t2c = true; }
  }
  if (expConfigOptions.expansionCoeffGradFlag &&
      expansionType1CoeffGrads.size() != num_levels)
    { expansionType1CoeffGrads.resize(num_levels); t1g = true; }

  if (t1c || t2c || t1g)
    for (i=0; i<num_levels; ++i) {
      const UShort3DArray& key_i = key[i];
      num_sets = key_i.size();
      if (t1c) expansionType1Coeffs[i].resize(num_sets);
      if (t2c) expansionType2Coeffs[i].resize(num_sets);
      if (t1g) expansionType1CoeffGrads[i].resize(num_sets);
      for (j=0; j<num_sets; ++j) {
	num_tp_pts = key_i[j].size();
	for (k=0; k<num_tp_pts; ++k) {
	  if (t1c) expansionType1Coeffs[i][j].sizeUninitialized(num_tp_pts);
	  if (t2c) expansionType2Coeffs[i][j].shapeUninitialized(
	    num_deriv_vars, num_tp_pts);
	  if (t1g) expansionType1CoeffGrads[i][j].shapeUninitialized(
	    num_deriv_vars, num_tp_pts);
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
  const UShort4DArray& key = hsg_driver->collocation_key();
  size_t i, j, k, l, num_levels = key.size(), num_sets, num_tp_pts, cntr = 0,
    num_deriv_vars = surrData.num_derivative_variables();

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
  for (i=1; i<num_levels; ++i) {
    const UShort3DArray& key_i = key[i];
    num_sets = key_i.size();
    for (j=0; j<num_sets; ++j) {
      num_tp_pts = key_i[j].size();
      for (k=0; k<num_tp_pts; ++k, ++cntr) {
	const RealVector& c_vars = surrData.continuous_variables(cntr);
	// coefficients are hierarchical surpluses
	if (expConfigOptions.expansionCoeffFlag) {
	  expansionType1Coeffs[i][j][k]
	    = surrData.response_function(cntr) - value(c_vars, i-1);
	  if (basisConfigOptions.useDerivs) {
	    const RealVector& data_grad = surrData.response_gradient(cntr);
	    const RealVector& prev_grad = gradient_basis_variables(c_vars, i-1);
	    Real* hier_grad = expansionType2Coeffs[i][j][k];
	    for (l=0; l<num_deriv_vars; ++l)
	      hier_grad[l] = data_grad[l] - prev_grad[l];
	    //RealVector hier_grad = surrData.response_gradient(cntr);
	    //hier_grad -= gradient_basis_variables(c_vars, i-1);
	    //Teuchos::setCol(hier_grad, (int)k, expansionType2Coeffs[i][j]);
	  }
	}
	if (expConfigOptions.expansionCoeffGradFlag) {
	  const RealVector& data_grad = surrData.response_gradient(cntr);
	  const RealVector& prev_grad = gradient_nonbasis_variables(c_vars,i-1);
	  Real* hier_grad = expansionType1CoeffGrads[i][j][k];
	  for (l=0; l<num_deriv_vars; ++l)
	    hier_grad[l] = data_grad[l] - prev_grad[l];
	  //RealVector hier_grad = surrData.response_gradient(cntr);
	  //hier_grad -= gradient_nonbasis_variables(c_vars, i-1);
	  //Teuchos::setCol(hier_grad, (int)k, expansionType1CoeffGrads[i][j]);
	}
      }
    }
  }

#ifdef INTERPOLATION_TEST

  // TO DO

  // SC should accurately interpolate the collocation data for TPQ and
  // SSG with fully nested rules, but will exhibit interpolation error
  // for SSG with other rules.
  index = 0;
  Real val, err, val_max_err = 0., grad_max_err = 0.,
       val_rmse = 0., grad_rmse = 0.;
  PCout << std::scientific << std::setprecision(WRITE_PRECISION);
  for (size_t i=0; i<numCollocPts; ++i, ++index) {
    const Real&       coeff1 = expansionType1Coeffs[i];
    const RealVector& c_vars = surrData.continuous_variables(index);
    val = value(c_vars);
    err = (std::abs(coeff1) > DBL_MIN) ? std::abs(1. - val/coeff1) :
                                         std::abs(coeff1 - val);
    PCout << "Colloc pt " << std::setw(3) << i+1
	  << ": truth value  = "  << std::setw(WRITE_PRECISION+7) << coeff1
	  << " interpolant = "    << std::setw(WRITE_PRECISION+7) << val
	  << " relative error = " << std::setw(WRITE_PRECISION+7) << err <<'\n';
    if (err > val_max_err) val_max_err = err; val_rmse += err * err;
    if (basisConfigOptions.useDerivs) {
      const Real*     coeff2 = expansionType2Coeffs[i];
      const RealVector& grad = gradient_basis_variables(c_vars);
      for (size_t j=0; j<numVars; ++j) {
	err = (std::abs(coeff2[j]) > DBL_MIN) ?
	  std::abs(1. - grad[j]/coeff2[j]) : std::abs(coeff2[j] - grad[j]);
	PCout << "               " << "truth grad_" << j+1 << " = "
	      << std::setw(WRITE_PRECISION+7) << coeff2[j] << " interpolant = "
	      << std::setw(WRITE_PRECISION+7) << grad[j] << " relative error = "
	      << std::setw(WRITE_PRECISION+7) << err << '\n';
	if (err > grad_max_err) grad_max_err = err; grad_rmse += err * err;
      }
    }
  }
  val_rmse = std::sqrt(val_rmse/numCollocPts);
  PCout << "\nValue interpolation errors:    " << std::setw(WRITE_PRECISION+7)
	<< val_max_err << " (max) "            << std::setw(WRITE_PRECISION+7)
	<< val_rmse    << " (RMS)\n";
  if (basisConfigOptions.useDerivs) {
    grad_rmse = std::sqrt(grad_rmse/numCollocPts/numVars);
    PCout << "Gradient interpolation errors: " << std::setw(WRITE_PRECISION+7)
	  << grad_max_err << " (max) "         << std::setw(WRITE_PRECISION+7)
	  << grad_rmse    << " (RMS)\n";
  }
#endif // INTERPOLATION_TEST
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
value(const RealVector& x, unsigned short level)
{
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::get_value()" << std::endl;
    abort_handler(-1);
  }

  Real approx_val = 0.;
  // loop over outer level indices
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  const UShort4DArray& key   = hsg_driver->collocation_key();
  SizetArray colloc_index; // empty -> default indexing
  size_t lev, set, num_sets;
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = expansionType1Coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = expansionType2Coeffs[lev];
    num_sets = sm_mi_l.size();
    for (set=0; set<num_sets; ++set)
      approx_val +=
	tensor_product_value(x, t1_coeffs_l[set], t2_coeffs_l[set],
			     sm_mi_l[set], key_l[set], colloc_index);
  }
  return approx_val;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, unsigned short level)
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
  // loop over outer level indices
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  const UShort4DArray& key   = hsg_driver->collocation_key();
  SizetArray colloc_index; // empty -> default indexing
  size_t lev, set, num_sets;
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = expansionType1Coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = expansionType2Coeffs[lev];
    num_sets = sm_mi_l.size();
    for (set=0; set<num_sets; ++set)
      approxGradient +=
	tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv,
			 unsigned short level)
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
  // loop over outer level indices
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  const UShort4DArray& key   = hsg_driver->collocation_key();
  SizetArray colloc_index; // empty -> default indexing
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = expansionType1Coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = expansionType2Coeffs[lev];
    num_sets = sm_mi_l.size();
    for (set=0; set<num_sets; ++set)
      approxGradient +=
	tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index, dvv);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x, unsigned short level)
{
  // Error check for required data
  size_t lev, set, num_sets, num_deriv_vars;
  if (expConfigOptions.expansionCoeffGradFlag) {
    if (expansionType1CoeffGrads.size() > level &&
	expansionType1CoeffGrads[level].size())
      num_deriv_vars = expansionType1CoeffGrads[level][0].numRows();
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
  // loop over outer level indices
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  const UShort4DArray& key   = hsg_driver->collocation_key();
  SizetArray colloc_index; // empty -> default indexing
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&            sm_mi_l = sm_mi[lev];
    const UShort3DArray&              key_l = key[lev];
    const RealMatrixArray& t1_coeff_grads_l = expansionType1CoeffGrads[lev];
    num_sets = sm_mi_l.size();
    for (set=0; set<num_sets; ++set)
      approxGradient +=
	tensor_product_gradient_nonbasis_variables(x, t1_coeff_grads_l[set],
	  sm_mi_l[set], key_l[set], colloc_index);
  }

  return approxGradient;
}


Real HierarchInterpPolyApproximation::stored_value(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not available in "
	  << "HierarchInterpPolyApproximation::stored_value()" << std::endl;
    abort_handler(-1);
  }

  Real approx_val = 0.;
  SizetArray colloc_index; // empty -> default indexing
  size_t lev, set, num_sets, num_levels = storedLevMultiIndex.size();
  for (lev=0; lev<num_levels; ++lev) {
    const UShort2DArray&       sm_mi_l = storedLevMultiIndex[lev];
    const UShort3DArray&         key_l = storedCollocKey[lev];
    const RealVectorArray& t1_coeffs_l = storedExpType1Coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = storedExpType2Coeffs[lev];
    num_sets = sm_mi_l.size();
    for (set=0; set<num_sets; ++set)
      approx_val +=
	tensor_product_value(x, t1_coeffs_l[set], t2_coeffs_l[set],
			     sm_mi_l[set], key_l[set], colloc_index);
  }
  return approx_val;
}


const RealVector& HierarchInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not available in HierarchInterpPoly"
	  << "Approximation::stored_gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != numVars)
    approxGradient.sizeUninitialized(numVars);
  approxGradient = 0.;
  SizetArray colloc_index; // empty -> default indexing
  size_t lev, set, num_sets, num_levels = storedLevMultiIndex.size();
  for (lev=0; lev<num_levels; ++lev) {
    const UShort2DArray&       sm_mi_l = storedLevMultiIndex[lev];
    const UShort3DArray&         key_l = storedCollocKey[lev];
    const RealVectorArray& t1_coeffs_l = storedExpType1Coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = storedExpType2Coeffs[lev];
    num_sets = sm_mi_l.size();
    for (set=0; set<num_sets; ++set)
      approxGradient +=
	tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x)
{
  // Error check for required data
  size_t lev, set, num_sets, num_deriv_vars,
    num_levels = storedLevMultiIndex.size();
  if (expConfigOptions.expansionCoeffGradFlag) {
    lev = num_levels - 1;
    if (storedExpType1CoeffGrads.size() == num_levels &&
	storedExpType1CoeffGrads[lev].size())
      num_deriv_vars = storedExpType1CoeffGrads[lev][0].numRows();
    else {
      PCerr << "Error: insufficient size in stored type1 expansion coefficient "
	    << "gradients in\n       HierarchInterpPolyApproximation::stored_"
	    << "gradient_nonbasis_variables()" << std::endl;
      abort_handler(-1);
    }
  }
  else {
    PCerr << "Error: expansion coefficient gradients not available in Hierarch"
	  << "InterpPolyApproximation::stored_gradient_nonbasis_variables()"
	  << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.;
  SizetArray colloc_index; // empty -> default indexing
  for (lev=0; lev<num_levels; ++lev) {
    const UShort2DArray&            sm_mi_l = storedLevMultiIndex[lev];
    const UShort3DArray&              key_l = storedCollocKey[lev];
    const RealMatrixArray& t1_coeff_grads_l = storedExpType1CoeffGrads[lev];
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
  Real mean = 0.;
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  const RealVector2DArray& t1_wts = hsg_driver->type1_weight_set_arrays();
  switch (basisConfigOptions.useDerivs) {
  case false:
    //for ( unsigned int i = 0; i<numCollocPts; ++i)
    //  mean += expansionType1Coeffs[i] * t1_wts[i];
    break;
  case true:
    const RealMatrix2DArray& t2_wts = hsg_driver->type2_weight_set_arrays();
    /*
    for ( unsigned int i = 0; i< numCollocPts; ++i) {
      mean += expansionType1Coeffs[i] * t1_wts[i];
      const Real* coeff2_i = expansionType2Coeffs[i];
      const Real* t2_wt_i = t2_wts[i];
      for (unsigned int j = 0; j< numVars; ++j)
	mean += coeff2_i[j] * t2_wt_i[j];
    }
    */
    break;
  }
  return mean;
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
  Real covar = 0.;
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  HierarchInterpPolyApproximation* hip_approx_2 = 
    static_cast<HierarchInterpPolyApproximation*>(poly_approx_2);
  Real mean_1 = mean(), mean_2 = hip_approx_2->mean();
  const RealVector2DArray& t1_coeffs_2 = hip_approx_2->expansionType1Coeffs;
  const RealVector2DArray& t1_wts = hsg_driver->type1_weight_set_arrays();
  switch (basisConfigOptions.useDerivs) {
  case false:
    /*
    for ( unsigned int i=0; i<numCollocPts; ++i)
      covar += (expansionType1Coeffs[i] - mean_1) * (t1_coeffs_2[i] - mean_2)
	*  t1_wts[i];
    */
    break;
  case true:
    const RealMatrix2DArray& t2_coeffs_2 = hip_approx_2->expansionType2Coeffs;
    const RealMatrix2DArray& t2_wts = hsg_driver->type2_weight_set_arrays();
    /*
    for ( unsigned int i = 0; i < numCollocPts; ++i) {
      // type1 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      Real coeff1_i_mm1 = expansionType1Coeffs[i] - mean_1;
      Real coeff1_2i_mm2 = t1_coeffs_2[i]          - mean_2;
      covar += coeff1_i_mm1 * coeff1_2i_mm2 * t1_wts[i];
      // type2 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      // --> interpolated gradients are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
      const Real *coeff2_i  = expansionType2Coeffs[i];
      const Real *coeff2_2i = t2_coeffs_2[i], *t2_wt_i = t2_wts[i];
      for (unsigned int j=0; j<numVars; ++j)
	covar  += (coeff1_i_mm1 * coeff2_2i[j] + coeff1_2i_mm2 * coeff2_i[j])
	  *  t2_wt_i[j];
    }
    */
    break;
  }
  return covar;
}


Real HierarchInterpPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  //TODO
  PCerr << "TODO: covariance in all variables" << std::endl;
  Real covar = 0.;
  return covar;
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
  /* TO DO: add another PolynomialApproximation formulation
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
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in InterpPoly"
	  << "Approximation::compute_numerical_expansion_moments()"<< std::endl;
    abort_handler(-1);
  }
  if (expansionMoments.length() != num_moments)
    expansionMoments.sizeUninitialized(num_moments);

  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  /* TO DO: add another PolynomialApproximation formulation
  size_t i, num_pts = surrData.size();
  RealVector t1_exp(num_pts);
  if (basisConfigOptions.useDerivs) {
    RealMatrix t2_exp(numVars, num_pts);
    for (i=0; i<num_pts; ++i) {
      const RealVector& c_vars = surrData.continuous_variables(i);
      t1_exp[i] = value(c_vars);
      Teuchos::setCol(gradient_basis_variables(c_vars), (int)i, t2_exp);
    }
    compute_numerical_moments(t1_exp, t2_exp,
                              hsg_driver->type1_weight_set_arrays(),
			      hsg_driver->type2_weight_set_arrays(),
			      expansionMoments);
  }
  else {
    for (i=0; i<num_pts; ++i)
      t1_exp[i] = value(surrData.continuous_variables(i));
    compute_numerical_moments(t1_exp, hsg_driver->type1_weight_set_arrays(),
			      expansionMoments);
  }
  */
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
      set_value = (int)std::pow(2.,int(numVars)) - (int)std::pow(2.,j) - 1; 
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
