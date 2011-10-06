/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        InterpPolyApproximation
//- Description:  Implementation code for InterpPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "InterpPolyApproximation.hpp"
#include "TensorProductDriver.hpp"
#include "SparseGridDriver.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG
//#define INTERPOLATION_TEST


namespace Pecos {

int InterpPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface of multiple dimensions
  return (configOptions.expansionCoeffFlag ||
	  configOptions.expansionCoeffGradFlag) ? 1 : 0;
}


void InterpPolyApproximation::
distribution_types(short& poly_type_1d, short& rule)
{
  switch (basisType) {
  case PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    poly_type_1d = (configOptions.useDerivs) ?
      PIECEWISE_CUBIC_INTERP : PIECEWISE_LINEAR_INTERP;
    rule = NEWTON_COTES;                    break;
  case GLOBAL_INTERPOLATION_POLYNOMIAL:
    poly_type_1d = (configOptions.useDerivs) ? HERMITE_INTERP : LAGRANGE_INTERP;
    rule = NO_RULE;                         break;
  default:
    poly_type_1d = NO_POLY; rule = NO_RULE; break;
  }
}


void InterpPolyApproximation::allocate_arrays()
{
  allocate_component_effects();
  allocate_total_effects();

  size_t num_deriv_vars = surrData.num_derivative_variables();
  if (configOptions.expansionCoeffFlag) {
    if (expansionType1Coeffs.length() != numCollocPts)
      expansionType1Coeffs.sizeUninitialized(numCollocPts);
    if ( configOptions.useDerivs &&
	 ( expansionType2Coeffs.numRows() != num_deriv_vars ||
	   expansionType2Coeffs.numCols() != numCollocPts ) )
      expansionType2Coeffs.shapeUninitialized(num_deriv_vars, numCollocPts);
  }
  if ( configOptions.expansionCoeffGradFlag &&
       ( expansionType1CoeffGrads.numRows() != num_deriv_vars ||
	 expansionType1CoeffGrads.numCols() != numCollocPts ) )
    expansionType1CoeffGrads.shapeUninitialized(num_deriv_vars, numCollocPts);

  // checking numCollocPts is insufficient due to anisotropy --> changes in
  // anisotropic weights could move points around without changing the total.
  //bool update_exp_form =
  //  ( (configOptions.expansionCoeffFlag &&
  //     expansionType1Coeffs.length()      != numCollocPts) ||
  //    (configOptions.expansionCoeffGradFlag &&
  //     expansionType1CoeffGrads.numCols() != numCollocPts ) );

  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray&   quad_order = tpq_driver->quadrature_order();

    // verify total number of collocation pts (should not include anchor pt)
    size_t i, j, num_colloc_pts = 1;
    for (i=0; i<numVars; ++i)
      num_colloc_pts *= quad_order[i];
    if (num_colloc_pts != numCollocPts) {
      PCerr << "Error: inconsistent total collocation point count in "
	    << "InterpPolyApproximation::allocate_arrays()" << std::endl;
      abort_handler(-1);
    }

    //bool update_exp_form = (quad_order != quadOrderPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials
    //if (update_exp_form) {
    //}

    // can't use quad_order > quadOrderPrev logic since only 1 pt set is stored
    bool update_basis_form = (quad_order != quadOrderPrev);
    if (update_basis_form)
      update_tensor_interpolation_basis();

    quadOrderPrev = quad_order;
    break;
  }
  case SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    unsigned short    ssg_level  = ssg_driver->level();
    const RealVector& aniso_wts  = ssg_driver->anisotropic_weights();

    //bool update_exp_form
    //  = (ssg_level != ssgLevelPrev || aniso_wts != ssgAnisoWtsPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials
    //if (update_exp_form) {
    //}

    // Ignore weights since they only reduce the interpolation depth from the
    // level and the basis update uses a coarse increment based on level.  This
    // matches isotropic sparse grids, but forces fewer and larger updates in
    // the case of anisotropic or generalized grids.
    bool update_basis_form
      = (ssgLevelPrev == USHRT_MAX || ssg_level > ssgLevelPrev);
    if (update_basis_form)
      update_sparse_interpolation_basis(ssg_level);

    ssgLevelPrev = ssg_level; ssgAnisoWtsPrev = aniso_wts;
    break;
  }
  }
}


void InterpPolyApproximation::compute_coefficients()
{
  if (!configOptions.expansionCoeffFlag &&
      !configOptions.expansionCoeffGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in "
	  << "InterpPolyApproximation::compute_coefficients().\n         "
	  << "Bypassing approximation construction." << std::endl;
    return;
  }

  // For testing of anchor point logic:
  //size_t index = surrData.size() - 1;
  //surrData.anchor_point(surrData.variables_data()[index],
  //                      surrData.response_data()[index]);
  //surrData.pop(1);

  size_t offset = 0;
  numCollocPts = surrData.size();
  // anchor point, if present, is the first expansionSample.
  if (surrData.anchor())
    { offset = 1; ++numCollocPts; }
  if (!numCollocPts) {
    PCerr << "Error: nonzero number of sample points required in "
	  << "InterpPolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }

  allocate_arrays();

  if (surrData.anchor()) {
    if (configOptions.expansionCoeffFlag) {
      expansionType1Coeffs[0] = surrData.anchor_function();
      if (configOptions.useDerivs)
	Teuchos::setCol(surrData.anchor_gradient(), 0, expansionType2Coeffs);
    }
    if (configOptions.expansionCoeffGradFlag)
      Teuchos::setCol(surrData.anchor_gradient(), 0, expansionType1CoeffGrads);
  }

  size_t index = 0;
  for (int i=offset; i<numCollocPts; ++i, ++index) {
    if (configOptions.expansionCoeffFlag) {
      expansionType1Coeffs[i] = surrData.response_function(index);
      if (configOptions.useDerivs)
	Teuchos::setCol(surrData.response_gradient(index), i,
			expansionType2Coeffs);
    }
    if (configOptions.expansionCoeffGradFlag)
      Teuchos::setCol(surrData.response_gradient(index), i,
		      expansionType1CoeffGrads);
  }

#ifdef INTERPOLATION_TEST
  // SC should accurately interpolate the collocation data for TPQ and
  // SSG with fully nested rules, but will exhibit interpolation error
  // for SSG with other rules.
  index = 0;
  for (size_t i=offset; i<numCollocPts; ++i, ++index) {
    const Real& coeff1 = expansionType1Coeffs[i];
    const Real&    val = value(surrData.continuous_variables(index));
    PCout << "Colloc pt " << std::setw(3) << i+1
	  << ": truth value  = " << std::setw(WRITE_PRECISION+7) << coeff1
	  << " interpolant = "   << std::setw(WRITE_PRECISION+7) << val
	  << " error = " << std::setw(WRITE_PRECISION+7)
	  << std::abs(coeff1 - val) << '\n';
    if (configOptions.useDerivs) {
      const Real*     coeff2 = expansionType2Coeffs[i];
      const RealVector& grad = gradient(surrData.continuous_variables(index));
      for (size_t j=0; j<numVars; ++j)
	PCout << "               " << "truth grad_" << j+1 << " = "
	      << std::setw(WRITE_PRECISION+7) << coeff2[j] << " interpolant = "
	      << std::setw(WRITE_PRECISION+7) << grad[j] << " error = "
	      << std::setw(WRITE_PRECISION+7) << std::abs(coeff2[j] - grad[j])
	      << '\n';
    }
  }
#endif // INTERPOLATION_TEST
}


void InterpPolyApproximation::increment_coefficients()
{
  bool err_flag = false;
  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID: {
    // As for allocate_arrays(), increments are performed in coarser steps
    // than may be strictly necessary: all increments are filled in for all
    // vars for a step in level (ignoring anisotropy or generalized indices).
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    const UShortArray& trial_set = ssg_driver->trial_set();
    unsigned short max_trial_index = 0;
    for (size_t i=0; i<numVars; ++i)
      if (trial_set[i] > max_trial_index)
	max_trial_index = trial_set[i];
    update_sparse_interpolation_basis(max_trial_index);
    break;
  }
  default:
    err_flag = true; break;
  }
  if (err_flag) {
    PCerr << "Error: unsupported grid definition in InterpPolyApproximation::"
	  << "increment_coefficients()" << std::endl;
    abort_handler(-1);
  }

  restore_expansion_coefficients();
}


void InterpPolyApproximation::decrement_coefficients()
{
  // leave polynomialBasis as is

  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID: {
    // move previous expansion data to current expansion
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    savedSmolyakMultiIndex.push_back(ssg_driver->trial_set());
    break;
  }
  }

  // not necessary to prune; next increment/restore/finalize takes care of this
  //if (configOptions.expansionCoeffFlag) {
  //  expansionType1Coeffs.resize(numCollocPts);
  //  if (configOptions.useDerivs) {
  //    size_t num_deriv_vars = expansionType2Coeffs.numRows();
  //    expansionType2Coeffs.reshape(num_deriv_vars, numCollocPts);
  //  }
  //}
  //if (configOptions.expansionCoeffGradFlag) {
  //  size_t num_deriv_vars = expansionType1CoeffGrads.numRows();
  //  expansionType1CoeffGrads.reshape(num_deriv_vars, numCollocPts);
  //}

  numCollocPts = surrData.size(); // data already decremented
  if (surrData.anchor())
    ++numCollocPts;
}


void InterpPolyApproximation::restore_coefficients()
{
  // leave polynomialBasis as is (a previous increment is being restored)

  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID: {
    // move previous expansion data to current expansion
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    std::deque<UShortArray>::iterator sit
      = std::find(savedSmolyakMultiIndex.begin(), savedSmolyakMultiIndex.end(),
		  ssg_driver->trial_set());
    if (sit != savedSmolyakMultiIndex.end())
      savedSmolyakMultiIndex.erase(sit);
    break;
  }
  }

  restore_expansion_coefficients();
}


void InterpPolyApproximation::finalize_coefficients()
{
  // leave polynomialBasis as is (all previous increments are being restored)

  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID:
    // move previous expansion data to current expansion
    savedSmolyakMultiIndex.clear();
    break;
  }

  restore_expansion_coefficients();
}


void InterpPolyApproximation::store_coefficients()
{
  if (configOptions.expansionCoeffFlag) {
    storedExpType1Coeffs   = expansionType1Coeffs;
    if (configOptions.useDerivs)
      storedExpType2Coeffs = expansionType2Coeffs;
  }
  if (configOptions.expansionCoeffGradFlag)
    storedExpType1CoeffGrads = expansionType1CoeffGrads;
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    storedCollocKey.resize(1); storedLevMultiIndex.resize(1);
    storedCollocKey[0]     = tpq_driver->collocation_key();
    storedLevMultiIndex[0] = tpq_driver->level_index();
    break;
  }
  case SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    storedLevMultiIndex = ssg_driver->smolyak_multi_index();
    storedLevCoeffs     = ssg_driver->smolyak_coefficients();
    storedCollocKey     = ssg_driver->collocation_key();
    storedCollocIndices = ssg_driver->collocation_indices();
    break;
  }
  }
}


void InterpPolyApproximation::combine_coefficients(short corr_type)
{
#ifdef DEBUG
  PCout << "Original type1 expansion coefficients prior to combination:\n";
  write_data(PCout, expansionType1Coeffs);
#endif // DEBUG

  // update expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads} by adding or
  // multiplying stored expansion evaluated at current collocation points
  size_t i, j, offset = 0, num_pts = surrData.size();
  bool anchor_pt = surrData.anchor();
  if (anchor_pt) { offset = 1; ++num_pts; }
  for (i=0; i<num_pts; ++i) {
    const RealVector& c_vars = (anchor_pt && i == 0) ?
      surrData.anchor_continuous_variables() :
      surrData.continuous_variables(i-offset);
    if (configOptions.expansionCoeffFlag) {
      if (corr_type == ADD_COMBINE) {
	// split up type1/type2 contribs so increments are performed properly
	expansionType1Coeffs[i] += stored_value(c_vars);
	if (configOptions.useDerivs) {
	  const RealVector& discrep_grad = stored_gradient(c_vars);
	  Real*          exp_t2_coeffs_i = expansionType2Coeffs[i];
	  size_t num_discrep_vars = discrep_grad.length();
	  for (j=0; j<num_discrep_vars; ++j)
	    exp_t2_coeffs_i[j] += discrep_grad[j];
	}
      }
      else if (corr_type == MULT_COMBINE) {
	Real discrep_val = stored_value(c_vars),
	          lf_val = expansionType1Coeffs[i]; // copy prior to update
	expansionType1Coeffs[i] *= discrep_val;
	if (configOptions.useDerivs) {
	  // hf = lf*discrep --> dhf/dx = dlf/dx*discrep + lf*ddiscrep/dx
	  const RealVector& discrep_grad = stored_gradient(c_vars);
	  Real*          exp_t2_coeffs_i = expansionType2Coeffs[i];
	  size_t num_discrep_vars = discrep_grad.length();
	  for (j=0; j<num_discrep_vars; ++j)
	    exp_t2_coeffs_i[j] = exp_t2_coeffs_i[j] * discrep_val
	                       + discrep_grad[j]    * lf_val;
	}
      }
    }
    if (configOptions.expansionCoeffGradFlag) {
      Real*   exp_grad_i = expansionType1CoeffGrads[i];
      /* TO DO
      const Real* grad_i = gradient(vars_set[i]);
      //if (corr_type == ADD_COMBINE)
      for (j=0; j<num_deriv_vars; ++j)
	exp_grad_i[j] += grad_i[j];
      //else if (corr_type == MULT_COMBINE)
      //for (j=0; j<num_deriv_vars; ++j)
      //  exp_grad_i[j] = ...;
      */
    }
  }
#ifdef DEBUG
  PCout << "Updated type1 expansion coefficients following combination:\n";
  write_data(PCout, expansionType1Coeffs);
#endif // DEBUG

  // clear stored data now that it has been combined
  if (configOptions.expansionCoeffFlag) {
    storedExpType1Coeffs.resize(0);
    if (configOptions.useDerivs) storedExpType2Coeffs.reshape(0,0);
  }
  if (configOptions.expansionCoeffGradFlag)
    storedExpType1CoeffGrads.reshape(0,0);
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    storedCollocKey.clear(); break;
  case SPARSE_GRID:
    storedLevMultiIndex.clear(); storedLevCoeffs.clear();
    storedCollocKey.clear();     storedCollocIndices.clear(); break;
  }
}


void InterpPolyApproximation::restore_expansion_coefficients()
{
  size_t offset = 0, new_colloc_pts = surrData.size();
  if (surrData.anchor())
    { offset = 1; ++new_colloc_pts; }

  if (configOptions.expansionCoeffFlag) {
    expansionType1Coeffs.resize(new_colloc_pts);
    if (configOptions.useDerivs) {
      size_t num_deriv_vars = expansionType2Coeffs.numRows();
      expansionType2Coeffs.reshape(num_deriv_vars, new_colloc_pts);
    }
  }
  if (configOptions.expansionCoeffGradFlag) {
    size_t num_deriv_vars = expansionType1CoeffGrads.numRows();
    expansionType1CoeffGrads.reshape(num_deriv_vars, new_colloc_pts);
  }

  size_t index = numCollocPts - offset;
  for (int i=numCollocPts; i<new_colloc_pts; ++i, ++index) {
    if (configOptions.expansionCoeffFlag) {
      expansionType1Coeffs[i] = surrData.response_function(index);
      if (configOptions.useDerivs)
	Teuchos::setCol(surrData.response_gradient(index), i,
			expansionType2Coeffs);
    }
    if (configOptions.expansionCoeffGradFlag)
      Teuchos::setCol(surrData.response_gradient(index), i,
		      expansionType1CoeffGrads);
  }

  numCollocPts = new_colloc_pts;
}


void InterpPolyApproximation::update_tensor_interpolation_basis()
{
  TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
  const UShortArray&   quad_order = tpq_driver->quadrature_order();
  const UShortArray&    lev_index = tpq_driver->level_index();

  // resize if needed (leaving previous levels unmodified)
  size_t i, j, k, basis_size = polynomialBasis.size();
  unsigned short max_order = quad_order[0];
  for (i=1; i<numVars; ++i)
    if (quad_order[i] > max_order)
      max_order = quad_order[i];
  // quad_order range is 1:m; quad_index range is 0:m-1
  if (max_order > basis_size) {
    polynomialBasis.resize(max_order);
    for (i=basis_size; i<max_order; ++i)
      polynomialBasis[i].resize(numVars);
  }

  // fill any required gaps in polynomialBasis.
  const Real3DArray& colloc_pts_1d = driverRep->collocation_points_array();
  short poly_type_1d; short rule; bool found; unsigned short l_index;
  distribution_types(poly_type_1d, rule);
  for (j=0; j<numVars; ++j) {
    l_index = lev_index[j];
    const RealArray& colloc_pts_1d_ij          =   colloc_pts_1d[l_index][j];
    std::vector<BasisPolynomial>& poly_basis_i = polynomialBasis[l_index];
    BasisPolynomial& poly_basis_ij = poly_basis_i[j];
    if (poly_basis_ij.is_null()) { // does not account for parametric changes
                                   // resulting in new pts for existing orders
      found = false;
      for (k=0; k<numVars; ++k)
	if (k != j && !poly_basis_i[k].is_null() &&
	    colloc_pts_1d_ij == poly_basis_i[k].interpolation_points())
	  { found = true; break; }
      if (found) // reuse previous instance via shared representation
	poly_basis_ij = poly_basis_i[k]; // shared rep
      else { // instantiate and initialize a new unique instance
	poly_basis_ij = BasisPolynomial(poly_type_1d, rule);
	poly_basis_ij.interpolation_points(colloc_pts_1d_ij);
      }
    }
  }
}


void InterpPolyApproximation::
update_sparse_interpolation_basis(unsigned short max_level)
{
  // resize if needed (leaving previous levels unmodified)
  size_t i, j, k, basis_size = polynomialBasis.size();
  // j range is 0:w inclusive; i range is 1:w+1 inclusive
  unsigned short num_levels = max_level + 1;
  if (num_levels > basis_size) {
    polynomialBasis.resize(num_levels);
    for (i=basis_size; i<num_levels; ++i)
      polynomialBasis[i].resize(numVars);
  }

  // fill gaps that may exist within any level (SparseGridDriver::
  // update_1d_collocation_points_weights() updates in an unstructured manner)
  const Real3DArray& colloc_pts_1d = driverRep->collocation_points_array();
  short poly_type_1d; short rule; bool found;
  distribution_types(poly_type_1d, rule);
  for (i=0; i<num_levels; ++i) { // i -> 0:num_levels-1 -> 0:ssg_level
    const Real2DArray& colloc_pts_1d_i = colloc_pts_1d[i];
    std::vector<BasisPolynomial>& poly_basis_i = polynomialBasis[i];
    for (j=0; j<numVars; ++j) {
      const RealArray& colloc_pts_1d_ij = colloc_pts_1d_i[j];
      BasisPolynomial&    poly_basis_ij =    poly_basis_i[j];
      if ( poly_basis_ij.is_null() &&  // doesn't account for parametric changes
	  !colloc_pts_1d_ij.empty()) { // resulting in new pts for existing levs
	found = false;
	for (k=0; k<numVars; ++k)
	  if (k != j && !poly_basis_i[k].is_null() &&
	      colloc_pts_1d_ij == colloc_pts_1d_i[k])  // vector equality
	    { found = true; break; }
	if (found) // reuse previous instance via shared representation
	  poly_basis_ij = poly_basis_i[k]; // shared rep
	else { // instantiate and initialize a new unique instance
	  poly_basis_ij = BasisPolynomial(poly_type_1d, rule);
	  poly_basis_ij.interpolation_points(colloc_pts_1d_ij);
	}
      }
    }
  }
}


void InterpPolyApproximation::
compute_numerical_response_moments(size_t num_moments)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in InterpPoly"
	  << "Approximation::compute_numerical_response_moments()" << std::endl;
    abort_handler(-1);
  }

  if (configOptions.useDerivs)
    compute_numerical_moments(num_moments, expansionType1Coeffs,
			      expansionType2Coeffs, numericalMoments);
  else
    compute_numerical_moments(num_moments, expansionType1Coeffs,
			      numericalMoments);
}


void InterpPolyApproximation::
compute_numerical_expansion_moments(size_t num_moments)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in InterpPoly"
	  << "Approximation::compute_numerical_expansion_moments()"<< std::endl;
    abort_handler(-1);
  }

  size_t i, offset = 0, num_pts = surrData.size();
  bool anchor_pt = surrData.anchor();
  if (anchor_pt) { offset = 1; ++num_pts; }
  RealVector t1_exp(num_pts);
  if (configOptions.useDerivs) {
    RealMatrix t2_exp(numVars, num_pts);
    for (i=0; i<num_pts; ++i) {
      const RealVector& c_vars = (anchor_pt && i == 0) ?
	surrData.anchor_continuous_variables() :
	surrData.continuous_variables(i-offset);
      t1_exp[i] = value(c_vars);
      Teuchos::setCol(gradient(c_vars), (int)i, t2_exp);
    }
    compute_numerical_moments(num_moments, t1_exp, t2_exp, expansionMoments);
  }
  else {
    for (i=0; i<num_pts; ++i) {
      const RealVector& c_vars = (anchor_pt && i == 0) ?
	surrData.anchor_continuous_variables() :
	surrData.continuous_variables(i-offset);
      t1_exp[i] = value(c_vars);
    }
    compute_numerical_moments(num_moments, t1_exp, expansionMoments);
  }
}


void InterpPolyApproximation::compute_component_effects()
{
  // perform subset sort
  constituentSets.resize(sobolIndices.length());
  get_subsets();

  // initialize partialVariance
  if (partialVariance.empty())
    partialVariance.sizeUninitialized(sobolIndices.length());
  partialVariance = 0.;

  const Real& mean           = numericalMoments[0];
  const Real& total_variance = numericalMoments[1];
  partialVariance[0]         = mean*mean; // init with mean sq

  // Solve for partial variance
  for (IntIntMIter map_iter=sobolIndexMap.begin();
       map_iter!=sobolIndexMap.end(); ++map_iter) {
    // partialVariance[0] stores the mean; it is not a component function
    // and does not follow the procedures for obtaining variance 
    if (map_iter->first) {
      partial_variance(map_iter->first);
      sobolIndices[map_iter->second]
	= partialVariance[map_iter->second]/total_variance;
      // total indices simply identify the membership of the sobolIndices 
      // and adds it to the appropriate bin
    }
  }
#ifdef DEBUG
  PCout << "In InterpPolyApproximation::compute_component_effects(), "
	<< "sobolIndices =\n"; write_data(PCout, sobolIndices);
#endif // DEBUG
}


void InterpPolyApproximation::compute_total_effects()
{
  // iterate through existing indices if all component indices are available
  totalSobolIndices = 0.; // init total indices
  if (configOptions.vbdControl == ALL_VBD)
    for (IntIntMIter itr=sobolIndexMap.begin(); itr!=sobolIndexMap.end(); ++itr)
      for (int k=0; k<numVars; ++k) {
        if (itr->first & (1 << k))
          totalSobolIndices[k] += sobolIndices[itr->second];
        totalSobolIndices[k] = std::abs(totalSobolIndices[k]);
      }

  // If not available, compute total indices independently.
  // This approach parallels partial_variance_integral where the algorithm is 
  // separated by integration approach.
  else {
    const Real& total_variance = numericalMoments[1];
    int j, set_value;
    switch (configOptions.expCoeffsSolnApproach) {
    case QUADRATURE: {
      TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
      const UShortArray&   quad_order = tpq_driver->quadrature_order();
      const UShortArray&    lev_index = tpq_driver->level_index();
      const UShort2DArray& colloc_key = tpq_driver->collocation_key();
      SizetArray colloc_index; // empty -> default indexing
      for (j=0; j<numVars; ++j) {
	// define set_value that includes all but index of interest
	set_value = (int)std::pow(2.,int(numVars)) - (int)std::pow(2.,j) - 1;
	totalSobolIndices[j] = std::abs(1. -
	  total_effects_integral(set_value, quad_order, lev_index,
				 colloc_key, colloc_index) / total_variance);
      }
      break;
    }
    case SPARSE_GRID: {
      SparseGridDriver*     ssg_driver = (SparseGridDriver*)driverRep;
      const IntArray&        sm_coeffs = ssg_driver->smolyak_coefficients();
      const UShort2DArray&    sm_index = ssg_driver->smolyak_multi_index();
      const UShort3DArray&  colloc_key = ssg_driver->collocation_key();
      const Sizet2DArray& colloc_index = ssg_driver->collocation_indices();
      // Smolyak recursion of anisotropic tensor products
      size_t i, num_smolyak_indices = sm_coeffs.size();
      UShortArray quad_order;
      // iterate each variable 
      for (j=0; j<numVars; ++j) {
	set_value = (int)std::pow(2.,int(numVars)) - (int)std::pow(2.,j) - 1; 
	for (i=0; i<num_smolyak_indices; ++i)
	  if (sm_coeffs[i]) {
	    ssg_driver->level_to_order(sm_index[i], quad_order);
	    totalSobolIndices[j] += sm_coeffs[i] *
	      total_effects_integral(set_value, quad_order, sm_index[i],
				     colloc_key[i], colloc_index[i]);
	  }
	totalSobolIndices[j]
	  = std::abs(1. - totalSobolIndices[j]/total_variance);
      }
      break;
    }
    }
  }
#ifdef DEBUG
  PCout << "In InterpPolyApproximation::compute_total_effects(), "
	<< "totalSobolIndices =\n"; write_data(PCout, totalSobolIndices);
#endif // DEBUG
}


Real InterpPolyApproximation::
total_effects_integral(int set_value, const UShortArray& quad_order,
		       const UShortArray& lev_index, const UShort2DArray& key,
		       const SizetArray& colloc_index)
{
  // Distinguish between non-members and members of the given set, set_value
  BoolDeque nonmember_vars(numVars,true); 
  int num_mem_exp_coeffs = 1; // number of expansion coeffs in
                              // member-variable-only expansion 
  IntVector indexing_factor; // factors indexing member variables 
  indexing_factor.sizeUninitialized(numVars);

  // create member variable key and get number of expansion coeffs in
  // member-variable-only expansion
  indexing_factor = 1;
  for (int k=0; k<numVars; ++k) {
    // if subset contains variable k, set key for variable k to true
    if (set_value & (1 << k)) {
      nonmember_vars[k] = false;	
      // information to properly index mem_exp_coeffs
      indexing_factor[k] = num_mem_exp_coeffs;
      num_mem_exp_coeffs *= quad_order[k];
    }	
  }
        
  // Create vector to store new coefficients
  RealVector mem_exp_coeffs(num_mem_exp_coeffs),
             mem_weights(num_mem_exp_coeffs);
 
  // Perform integration over non-member variables and store indices
  // of new expansion
  size_t i, j, num_colloc_pts = key.size();
  const Real3DArray& colloc_wts_1d
    = driverRep->type1_collocation_weights_array();
  for (i=0; i <num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    size_t mem_exp_coeffs_index = 0;	
    Real prod_i_nonmembers = 1, prod_i_members = 1;
    for (j=0; j<numVars; ++j)
      // Save the product of the weights of the member and non-member variables 
      if (nonmember_vars[j])
        prod_i_nonmembers    *= colloc_wts_1d[lev_index[j]][j][key_i[j]];
      else {
	// Convert key to corresponding index on mem_exp_coeffs
	mem_exp_coeffs_index += key_i[j]*indexing_factor[j];
        prod_i_members       *= colloc_wts_1d[lev_index[j]][j][key_i[j]];
      }
 
    // mem_weights is performed more time than necessary here, but it
    // seems to be the simplest place to put it
    mem_weights[mem_exp_coeffs_index] = prod_i_members;
    // sort coefficients by the "signature" of the member variables
    // (i.e. mem_exp_coeffs_index)
    unsigned short c_index = (colloc_index.empty()) ? i : colloc_index[i];
    mem_exp_coeffs[mem_exp_coeffs_index]
      += expansionType1Coeffs[c_index]*prod_i_nonmembers;
  }
 
  // Now integrate over the remaining variables	
  Real  integral   = 0;
  const Real& total_mean = numericalMoments[0];
  for (i=0; i<num_mem_exp_coeffs; ++i)
    integral += std::pow(mem_exp_coeffs[i] - total_mean, 2.)*mem_weights[i];
  return integral;
}


/** Find constituent subsets. */
void InterpPolyApproximation::get_subsets()
{
  // includes the "zero" set
  //int num_subsets = sobolIndices.length(); 

  // Here we want to utilize the integer representation of the subset
  // but we want to store it in a size appropriate container
  // so finding lower sets is given the argument of the integer rep (->first)
  // and stored in constituentSets in size-appropriate-index-map (->second)
  for (IntIntMIter map_iter=sobolIndexMap.begin();
       map_iter!=sobolIndexMap.end(); ++map_iter) {
    lower_sets(map_iter->first, constituentSets[map_iter->second]);
    constituentSets[map_iter->second].erase(map_iter->first);
  }
}


/** For input set, recursively finds constituent subsets with one
    fewer element */
void InterpPolyApproximation::
lower_sets(int plus_one_set, IntSet& top_level_set)
{
  // if this set has been stored before, stop
  if (top_level_set.count(plus_one_set))
    return;
  // otherwise store current set
  else
    top_level_set.insert(plus_one_set);
  // and find lower level sets
  for (int k=0; k<numVars; ++k)
    // this performs a bitwise comparison by shifting 1 by k spaces 
    // and comparing that to a binary form of plus_one_set; this allows 
    // the variable membership using integers instead of a d-array of bools
    if (plus_one_set & (1 << k)) 
      // if subset i contains variable k, remove that variable from the set 
      // by converting the bit-form of (1<<k) to an integer and subtract from
      // the plus_one_set
      lower_sets(plus_one_set-(int)std::pow(2.0,k),top_level_set);
}


/** Computes the variance of component functions. Assumes that all
    subsets of set_value have been computed in advance which will be
    true so long as the partial_variance is called following
    appropriate enumeration of set value  */
void InterpPolyApproximation::partial_variance(int set_value)
{
  // Computes the integral first
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray&   quad_order = tpq_driver->quadrature_order();
    const UShortArray&    lev_index = tpq_driver->level_index();
    const UShort2DArray& colloc_key = tpq_driver->collocation_key();
    SizetArray colloc_index; // empty -> default indexing
    partialVariance[sobolIndexMap[set_value]]
      = partial_variance_integral(set_value, quad_order, lev_index,
				  colloc_key, colloc_index);
    break;
  }
  case SPARSE_GRID: {
    SparseGridDriver*     ssg_driver = (SparseGridDriver*)driverRep;
    const IntArray&        sm_coeffs = ssg_driver->smolyak_coefficients();
    const UShort2DArray&    sm_index = ssg_driver->smolyak_multi_index();
    const UShort3DArray&  colloc_key = ssg_driver->collocation_key();
    const Sizet2DArray& colloc_index = ssg_driver->collocation_indices();
    // Smolyak recursion of anisotropic tensor products
    size_t i, num_smolyak_indices = sm_coeffs.size();
    UShortArray quad_order;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i]) {
	ssg_driver->level_to_order(sm_index[i], quad_order);
	partialVariance[sobolIndexMap[set_value]] += sm_coeffs[i]
	  * partial_variance_integral(set_value, quad_order, sm_index[i],
				      colloc_key[i], colloc_index[i]);
      }
    break;
  }
  }
  
  // Now subtract the contributions from constituent subsets
  IntSet::iterator itr;
  for (itr  = constituentSets[sobolIndexMap[set_value]].begin();
       itr != constituentSets[sobolIndexMap[set_value]].end(); ++itr) 
    partialVariance[sobolIndexMap[set_value]]
      -= partialVariance[sobolIndexMap[*itr]];
}


/** Forms an interpolant over variables that are members of the given set.
    Finds the variance of the interpolant w.r.t. the variables in the set.
    Overloaded version supporting Smolyak sparse grids. */
Real InterpPolyApproximation::
partial_variance_integral(int set_value, const UShortArray& quad_order,
			  const UShortArray& lev_index,
			  const UShort2DArray& key,
			  const SizetArray& colloc_index)
{
  // Distinguish between non-members and members of the given set, set_value
  BoolDeque nonmember_vars(numVars,true); 
  int num_mem_exp_coeffs = 1; // number of expansion coeffs in
                              // member-variable-only expansion 
  IntVector indexing_factor; // factors indexing member variables 
  indexing_factor.sizeUninitialized(numVars);

  // create member variable key and get number of expansion coeffs in
  // member-variable-only expansion
  indexing_factor = 1;
  for (int k=0; k<numVars; ++k) {
    // if subset contains variable k, set key for variable k to true
    if (set_value & (1 << k)) {
      nonmember_vars[k] = false;	
      // information to properly index mem_exp_coeffs
      indexing_factor[k] = num_mem_exp_coeffs;
      num_mem_exp_coeffs *= quad_order[k];
    }	
  }
	
  // Create vector to store new coefficients
  RealVector mem_exp_coeffs(num_mem_exp_coeffs), 
             mem_weights(num_mem_exp_coeffs);

  // Perform integration over non-member variables and store indices
  // of new expansion
  size_t i, j, num_colloc_pts = key.size();
  const Real3DArray& colloc_wts_1d
    = driverRep->type1_collocation_weights_array();
  for (i=0; i <num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    size_t mem_exp_coeffs_index = 0;	
    Real prod_i_nonmembers = 1, prod_i_members = 1;
    for (j=0; j<numVars; ++j)
      // Save the product of the weights of the member and non-member variables 
      if (nonmember_vars[j])
	prod_i_nonmembers    *= colloc_wts_1d[lev_index[j]][j][key_i[j]];
      else {
	// Convert key to corresponding index on mem_exp_coeffs
	mem_exp_coeffs_index += key_i[j]*indexing_factor[j];
	prod_i_members       *= colloc_wts_1d[lev_index[j]][j][key_i[j]];
      }

    // mem_weights is performed more time than necessary here, but it
    // seems to be the simplest place to put it
    mem_weights[mem_exp_coeffs_index] = prod_i_members;
    // sort coefficients by the "signature" of the member variables
    // (i.e. mem_exp_coeffs_index)
    unsigned short c_index = (colloc_index.empty()) ? i : colloc_index[i];
    mem_exp_coeffs[mem_exp_coeffs_index]
      += expansionType1Coeffs[c_index]*prod_i_nonmembers;
  }

  // Now integrate over the remaining variables	
  Real integral = 0;
  for (i=0; i<num_mem_exp_coeffs; ++i)
    integral += std::pow(mem_exp_coeffs[i], 2.)*mem_weights[i];
  return integral;	
}

} // namespace Pecos
