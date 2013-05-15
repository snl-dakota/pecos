/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        ProjectOrthogPolyApproximation
//- Description:  Implementation code for ProjectOrthogPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "ProjectOrthogPolyApproximation.hpp"
#include "TensorProductDriver.hpp"
#include "CombinedSparseGridDriver.hpp"
#include "CubatureDriver.hpp"
#include "pecos_global_defs.hpp"

//#define DEBUG


namespace Pecos {

int ProjectOrthogPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface in multiple dimensions
  if (expConfigOptions.expansionCoeffFlag ||
      expConfigOptions.expansionCoeffGradFlag)
    return 1; // quadrature: (int)pow((Real)MIN_GAUSS_PTS, numVars);
  else
    return 0;
}


void ProjectOrthogPolyApproximation::allocate_arrays()
{
  // if base class version not invoked, invoke here
  if (expConfigOptions.expCoeffsSolnApproach != SAMPLING) {
    allocate_component_effects();
    allocate_total_effects();

    if (expansionMoments.empty())
      expansionMoments.sizeUninitialized(2);
  }
  // no combination by default, even if stored{MultiIndex,ExpCoeffs,
  // ExpCoeffGrads} are defined.  Redefined by attribute passed in
  // combine_coefficients(short).
  storedExpCombineType = NO_COMBINE;

  // update_exp_form controls when to update (refinement) and when not to
  // update (subIterator execution) an expansion's multiIndex definition.
  // Simple logic of updating if previous number of points != current number
  // is not robust enough for anisotropic updates --> track using Prev arrays.
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray& quad_order = tpq_driver->quadrature_order();
    bool update_exp_form = (quad_order != quadOrderPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials?

    if (update_exp_form) {
      UShortArray integrand_order(numVars);
      quadrature_order_to_integrand_order(quad_order, integrand_order);
      integrand_order_to_expansion_order(integrand_order, approxOrder);
      tensor_product_multi_index(approxOrder, multiIndex);
      size_expansion();
      quadOrderPrev = quad_order; // update reference point
    }
    PCout << "Orthogonal polynomial approximation order = { ";
    for (size_t i=0; i<numVars; ++i) PCout << approxOrder[i] << ' ';
    PCout << "} using tensor-product expansion of " << numExpansionTerms
	  << " terms\n";
    break;
  }
  case CUBATURE: {
    CubatureDriver* cub_driver = (CubatureDriver*)driverRep;
    //unsigned short cub_int_order = cub_driver->integrand_order();
    //bool update_exp_form = (cub_int_order != cubIntOrderPrev);

    //if (update_exp_form) {
      UShortArray integrand_order(numVars, cub_driver->integrand_order());
      integrand_order_to_expansion_order(integrand_order, approxOrder);
      total_order_multi_index(approxOrder, multiIndex);
      size_expansion();
      //cubIntOrderPrev = cub_int_order; // update reference point
    //}
    PCout << "Orthogonal polynomial approximation order = { ";
    for (size_t i=0; i<numVars; ++i)
      PCout << approxOrder[i] << ' ';
    PCout << "} using total-order expansion of " << numExpansionTerms
	  << " terms\n";
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    unsigned short    ssg_level = csg_driver->level();
    const RealVector& aniso_wts = csg_driver->anisotropic_weights();
    bool update_exp_form = (ssg_level != ssgLevelPrev ||
      aniso_wts != ssgAnisoWtsPrev || expConfigOptions.refinementControl ==
      DIMENSION_ADAPTIVE_CONTROL_GENERALIZED);
    // *** TO DO: capture updates to parameterized/numerical polynomials?

    if (update_exp_form) { // compute and output number of terms
      sparse_grid_multi_index(multiIndex);
      size_expansion();
      ssgLevelPrev = ssg_level; ssgAnisoWtsPrev = aniso_wts; // update ref pts
    }
    PCout << "Orthogonal polynomial approximation level = " << ssg_level
	  << " using tensor integration and tensor sum expansion of "
	  << numExpansionTerms << " terms\n"; break;
    break;
  }
  default: // SAMPLING
    OrthogPolyApproximation::allocate_arrays(); // default implementation
    break;
  }
}


void ProjectOrthogPolyApproximation::compute_coefficients()
{
  if (!expConfigOptions.expansionCoeffFlag &&
      !expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in "
	  << "ProjectOrthogPolyApproximation::compute_coefficients().\n"
	  << "         Bypassing approximation construction." << std::endl;
    return;
  }

  // For testing of anchor point logic:
  //size_t index = surrData.size() - 1;
  //surrData.anchor_point(surrData.variables_data()[index],
  //                      surrData.response_data()[index]);
  //surrData.pop(1);

  // anchor point, if present, is handled differently for different
  // expCoeffsSolnApproach settings:
  //   SAMPLING:   treat it as another data point
  //   QUADRATURE/CUBATURE/COMBINED_SPARSE_GRID: error
  //   LEAST_SQ_REGRESSION: use equality-constrained least squares
  size_t i, j, num_total_pts = surrData.size();
  if (surrData.anchor())
    ++num_total_pts;
  if (!num_total_pts) {
    PCerr << "Error: nonzero number of sample points required in ProjectOrthog"
	  << "PolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }

  // Array sizing can be divided into two parts:
  // > data used in all cases (size in allocate_arrays())
  // > data not used in expansion import case (size here)
  allocate_arrays();
#ifdef DEBUG
  gradient_check();
#endif // DEBUG

  // calculate polynomial chaos coefficients
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    // verify quad_order stencil matches num_total_pts
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray& quad_order = tpq_driver->quadrature_order();
    if (quad_order.size() != numVars) {
      PCerr << "Error: quadrature order array is not consistent with number of "
	    << "variables (" << numVars << ")\n       in ProjectOrthogPoly"
	    << "Approximation::compute_coefficients()." << std::endl;
      abort_handler(-1);
    }
    size_t num_colloc_pts = 1;
    for (i=0; i<numVars; ++i)
      num_colloc_pts *= quad_order[i];
    if (num_total_pts != num_colloc_pts) {
      PCerr << "Error: number of current points (" << num_total_pts
	    << ") is not consistent with\n       quadrature data in Project"
	    << "OrthogPolyApproximation::compute_coefficients()." << std::endl;
      abort_handler(-1);
    }

#ifdef DEBUG
    for (i=0; i<numVars; ++i) {
      OrthogonalPolynomial* poly_rep
	= (OrthogonalPolynomial*)polynomialBasis[i].polynomial_rep();
      for (j=1; j<=quad_order[i]; ++j)
	poly_rep->gauss_check(j);
    }
#endif // DEBUG

    // single expansion integration
    integration_checks();
    integrate_expansion(multiIndex, surrData.variables_data(),
			surrData.response_data(),
			driverRep->type1_weight_sets(),
			expansionCoeffs, expansionCoeffGrads);
    break;
  }
  case CUBATURE:
    // single expansion integration
    integration_checks();
    integrate_expansion(multiIndex, surrData.variables_data(),
			surrData.response_data(),
			driverRep->type1_weight_sets(),
			expansionCoeffs, expansionCoeffGrads);
    break;
  case COMBINED_SPARSE_GRID: {
    // multiple tensor expansion integrations
    // Note: allocate_arrays() calls sparse_grid_multi_index() which uses
    // append_multi_index() to build multiIndex.
    if (expConfigOptions.expansionCoeffFlag)     expansionCoeffs = 0.;
    if (expConfigOptions.expansionCoeffGradFlag) expansionCoeffGrads = 0.;
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)driverRep;
    const IntArray& sm_coeffs = csg_driver->smolyak_coefficients();
    size_t i, num_tensor_grids = tpMultiIndex.size(); int coeff;
    SDVArray tp_data_vars; SDRArray tp_data_resp;
    RealVector tp_wts, tp_coeffs; RealMatrix tp_coeff_grads;
    bool store_tp = (expConfigOptions.refinementControl ==
		     DIMENSION_ADAPTIVE_CONTROL_GENERALIZED);
    // loop over tensor-products, forming sub-expansions, and sum them up
    for (i=0; i<num_tensor_grids; ++i) {
      // form tp_data_vars, tp_data_resp, tp_wts using collocKey et al.
      integration_data(i, tp_data_vars, tp_data_resp, tp_wts);

      // form tp_multi_index from tpMultiIndexMap
      //map_tensor_product_multi_index(tp_multi_index, i);

      // form tp expansion coeffs
      RealVector& tp_coeffs_i = (store_tp) ? tpExpansionCoeffs[i] : tp_coeffs;
      RealMatrix& tp_grads_i
	= (store_tp) ? tpExpansionCoeffGrads[i] : tp_coeff_grads;
      integrate_expansion(tpMultiIndex[i], tp_data_vars, tp_data_resp, tp_wts,
			  tp_coeffs_i, tp_grads_i);

      // sum tensor product coeffs/grads into expansion coeffs/grads
      coeff = sm_coeffs[i];
      if (coeff)
	overlay_expansion(tpMultiIndexMap[i], tp_coeffs_i, tp_grads_i, coeff);
    }
    //if (!reEntrantFlag) {
    //  csg_driver->clear_smolyak_arrays();
    //  csg_driver->clear_collocation_arrays();
    //  tpMultiIndex.clear(); tpMultiIndexMap.clear();
    //}
    break;
  }
  case SAMPLING:
    surrData.data_checks();
    expectation();
    break;
  default:
    PCerr << "Error: unsupported expCoeffsSolnApproach in ProjectOrthogPoly"
	  << "Approximation::compute_coefficients()" << std::endl;
    abort_handler(-1);
    break;
  }

  computedMean = computedVariance = 0;
}


void ProjectOrthogPolyApproximation::increment_coefficients()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in ProjectOrthogPoly"
	  << "Approximation::increment_coefficients()" << std::endl;
    abort_handler(-1);
  }

  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  size_t last_index = tpMultiIndex.size();
  // size tpMultiIndex and tpExpansion{Coeffs,CoeffGrads}
  size_t new_size = last_index+1;
  tpMultiIndex.resize(new_size);
  tpExpansionCoeffs.resize(new_size);
  tpExpansionCoeffGrads.resize(new_size);
  // update tpMultiIndex
  UShortArray exp_order(numVars);
  sparse_grid_level_to_expansion_order(csg_driver->trial_set(), exp_order);
  tensor_product_multi_index(exp_order, tpMultiIndex[last_index]);
  // update multiIndex and numExpansionTerms
  append_multi_index(tpMultiIndex[last_index], multiIndex, tpMultiIndexMap,
		     tpMultiIndexMapRef);
  resize_expansion();
  // form tp_data_pts, tp_wts using collocKey et al.
  SDVArray tp_data_vars; SDRArray tp_data_resp; RealVector tp_wts;
  integration_data(last_index, tp_data_vars, tp_data_resp, tp_wts);
  // form trial expansion coeffs/grads
  integrate_expansion(tpMultiIndex[last_index], tp_data_vars, tp_data_resp,
		      tp_wts, tpExpansionCoeffs[last_index],
		      tpExpansionCoeffGrads[last_index]);
  // sum trial expansion into expansionCoeffs/expansionCoeffGrads
  append_tensor_expansions(last_index);
  // cleanup
  //if (!reEntrantFlag) {
  //  csg_driver->clear_smolyak_arrays();
  //  csg_driver->clear_collocation_arrays();
  //  tpMultiIndex.clear(); tpMultiIndexMap.clear();
  //}

  computedMean = computedVariance = 0;
}


void ProjectOrthogPolyApproximation::decrement_coefficients()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in ProjectOrthogPoly"
	  << "Approximation::decrement_coefficients()" << std::endl;
    abort_handler(-1);
  }

  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  // reset expansion{Coeffs,CoeffGrads}: (set in append_tensor_expansions())
  expansionCoeffs     = prevExpCoeffs;
  expansionCoeffGrads = prevExpCoeffGrads;

  // reset multiIndex and numExpansionTerms:
  numExpansionTerms   = tpMultiIndexMapRef.back();
  multiIndex.resize(numExpansionTerms); // truncate previous increment
  // resize not necessary since (1) already updated from prevExp and 
  // (2) not updating expansion on decrement (next increment updates).
  //resize_expansion();

  // reset tensor-product bookkeeping and save restorable data
  savedLevMultiIndex.push_back(csg_driver->trial_set());
  savedTPMultiIndex.push_back(tpMultiIndex.back());
  savedTPMultiIndexMap.push_back(tpMultiIndexMap.back());
  savedTPMultiIndexMapRef.push_back(numExpansionTerms);
  savedTPExpCoeffs.push_back(tpExpansionCoeffs.back());
  savedTPExpCoeffGrads.push_back(tpExpansionCoeffGrads.back());
  tpMultiIndex.pop_back();       tpMultiIndexMap.pop_back();
  tpMultiIndexMapRef.pop_back(); tpExpansionCoeffs.pop_back();
  tpExpansionCoeffGrads.pop_back();

  computedMean = computedVariance = 0;
}


void ProjectOrthogPolyApproximation::restore_coefficients()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in ProjectOrthogPoly"
	  << "Approximation::restore_coefficients()" << std::endl;
    abort_handler(-1);
  }

  // move previous expansion data to current expansion
  size_t last_index = tpMultiIndex.size();
  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  std::deque<UShortArray>::iterator sit
    = std::find(savedLevMultiIndex.begin(), savedLevMultiIndex.end(),
		csg_driver->trial_set());
  size_t index_star = std::distance(savedLevMultiIndex.begin(), sit);
  savedLevMultiIndex.erase(sit);
  std::deque<UShort2DArray>::iterator iit = savedTPMultiIndex.begin();
  std::deque<SizetArray>::iterator    mit = savedTPMultiIndexMap.begin();
  std::deque<size_t>::iterator        rit = savedTPMultiIndexMapRef.begin();
  std::deque<RealVector>::iterator    cit = savedTPExpCoeffs.begin();
  std::deque<RealMatrix>::iterator    git = savedTPExpCoeffGrads.begin();
  std::advance(iit, index_star);      std::advance(mit, index_star);
  std::advance(rit, index_star);      std::advance(cit, index_star);
  std::advance(git, index_star);
  tpMultiIndex.push_back(*iit);          savedTPMultiIndex.erase(iit);
  tpMultiIndexMap.push_back(*mit);       savedTPMultiIndexMap.erase(mit);
  tpMultiIndexMapRef.push_back(*rit);    savedTPMultiIndexMapRef.erase(rit);
  tpExpansionCoeffs.push_back(*cit);     savedTPExpCoeffs.erase(cit);
  tpExpansionCoeffGrads.push_back(*git); savedTPExpCoeffGrads.erase(git);
  // update multiIndex and numExpansionTerms
  append_multi_index(tpMultiIndex[last_index], tpMultiIndexMap[last_index],
		     tpMultiIndexMapRef[last_index], multiIndex);
  resize_expansion();
  // sum trial expansion into expansionCoeffs/expansionCoeffGrads
  append_tensor_expansions(last_index);

  computedMean = computedVariance = 0;
}


void ProjectOrthogPolyApproximation::finalize_coefficients()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in ProjectOrthogPoly"
	  << "Approximation::finalize_coefficients()" << std::endl;
    abort_handler(-1);
  }

  size_t start_index = tpMultiIndex.size();
  // update multiIndex and numExpansionTerms
  std::deque<UShort2DArray>::iterator iit = savedTPMultiIndex.begin();
  std::deque<SizetArray>::iterator    mit = savedTPMultiIndexMap.begin();
  std::deque<size_t>::iterator        rit = savedTPMultiIndexMapRef.begin();
  for (; iit!=savedTPMultiIndex.end(); ++iit, ++mit, ++rit)
    append_multi_index(*iit, *mit, *rit, multiIndex);
  resize_expansion();
  // move previous expansion data to current expansion
  tpMultiIndex.insert(tpMultiIndex.end(), savedTPMultiIndex.begin(),
    savedTPMultiIndex.end());
  tpMultiIndexMap.insert(tpMultiIndexMap.end(), savedTPMultiIndexMap.begin(),
    savedTPMultiIndexMap.end());
  tpMultiIndexMapRef.insert(tpMultiIndexMapRef.end(),
    savedTPMultiIndexMapRef.begin(), savedTPMultiIndexMapRef.end());
  tpExpansionCoeffs.insert(tpExpansionCoeffs.end(), savedTPExpCoeffs.begin(),
    savedTPExpCoeffs.end());
  tpExpansionCoeffGrads.insert(tpExpansionCoeffGrads.end(),
    savedTPExpCoeffGrads.begin(), savedTPExpCoeffGrads.end());
  savedLevMultiIndex.clear();     savedTPMultiIndex.clear();
  savedTPMultiIndexMap.clear();   savedTPMultiIndexMapRef.clear();
  savedTPExpCoeffs.clear();       savedTPExpCoeffGrads.clear();
  // sum remaining trial expansions into expansionCoeffs/expansionCoeffGrads
  append_tensor_expansions(start_index);

  computedMean = computedVariance = 0;
}


void ProjectOrthogPolyApproximation::store_coefficients()
{
  // Store the aggregated expansion data.  This approach is preferred to
  // appending to savedTP{MultiIndex,Coeffs,CoeffGrads} since the savedTP
  // approach is less general (TP and sum of TP only), less memory
  // efficient (tensor redundancies in sparse grids), and causes ambiguity
  // in finalize_coefficients() for generalized sparse grids.
  storedMultiIndex = multiIndex;
  if (expConfigOptions.expansionCoeffFlag)
    storedExpCoeffs = expansionCoeffs;
  if (expConfigOptions.expansionCoeffGradFlag)
    storedExpCoeffGrads = expansionCoeffGrads;

  switch (expConfigOptions.expCoeffsSolnApproach) { // approach-specific storage
  case COMBINED_SPARSE_GRID: { // sum of tensor-product expansions
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    storedLevMultiIndex = csg_driver->smolyak_multi_index(); break;
  }
  default: // tensor and total-order expansions
    storedApproxOrder = approxOrder;                         break;
  }
}


void ProjectOrthogPolyApproximation::combine_coefficients(short combine_type)
{
  // based on incoming combine_type, combine the data stored previously
  // by store_coefficients()

  // storedExpCombineType used later in compute_numerical_response_moments()
  storedExpCombineType = combine_type;

  switch (combine_type) {
  case ADD_COMBINE: {
    // Note: would like to preserve tensor indexing (at least for QUADRATURE
    // case) so that Horner's rule performance opt could be used within
    // tensor_product_value()).  However, a tensor result in the overlay
    // will not occur unless one expansion order dominates the other (partial
    // domination results in sum of tensor expansions as for sparse grids).
    // Therefore, stick with the general-purpose expansion overlay and exclude
    // tensor_product_value() usage for combined coefficient sets.

    // base class version is sufficient; no specialization based on exp form
    OrthogPolyApproximation::combine_coefficients(combine_type);
    break;
  }
  case MULT_COMBINE: {
    // compute form of product expansion
    switch (expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE: { // product of two tensor-product expansions
      for (size_t i=0; i<numVars; ++i)
	approxOrder[i] += storedApproxOrder[i];
      UShort2DArray multi_index_prod;
      tensor_product_multi_index(approxOrder, multi_index_prod);

      // perform the multiplication of current and stored expansions
      multiply_expansion(storedMultiIndex, storedExpCoeffs, storedExpCoeffGrads,
			 multi_index_prod);
      computedMean = computedVariance = 0;
      break;
    }
    case COMBINED_SPARSE_GRID: { // product of two sums of tensor-product exp.
      CombinedSparseGridDriver* csg_driver
	= (CombinedSparseGridDriver*)driverRep;
      // filter out dominated Smolyak multi-indices that don't contribute
      // to the definition of the product expansion
      UShort2DArray curr_pareto, stored_pareto;
      update_pareto(csg_driver->smolyak_multi_index(), curr_pareto);
      update_pareto(storedLevMultiIndex,             stored_pareto);
      size_t i, j, k, num_stored_mi = stored_pareto.size(),
	num_curr_mi = curr_pareto.size();
      // overlay each product expansion from the tensor-product combinations
      UShortArray exp_order_i, exp_order_j, exp_order_prod(numVars);
      UShort2DArray multi_index_prod, tp_multi_index_prod;
      for (i=0; i<num_stored_mi; ++i) {
	sparse_grid_level_to_expansion_order(stored_pareto[i], exp_order_i);
	for (j=0; j<num_curr_mi; ++j) {
	  sparse_grid_level_to_expansion_order(curr_pareto[j], exp_order_j);
	  for (k=0; k<numVars; ++k)
	    exp_order_prod[k] = exp_order_i[k] + exp_order_j[k];
	  tensor_product_multi_index(exp_order_prod, tp_multi_index_prod);
	  append_multi_index(tp_multi_index_prod, multi_index_prod);
	}
      }

      // perform the multiplication of current and stored expansions
      multiply_expansion(storedMultiIndex, storedExpCoeffs, storedExpCoeffGrads,
			 multi_index_prod);
      computedMean = computedVariance = 0;
      break;
    }
    default:
      // base class version supports product of two total-order expansions
      OrthogPolyApproximation::combine_coefficients(combine_type);
      break;
    }
    break;
  }
  case ADD_MULT_COMBINE:
    // base class manages this placeholder
    OrthogPolyApproximation::combine_coefficients(combine_type);
    break;
  }
}


void ProjectOrthogPolyApproximation::
sparse_grid_multi_index(UShort2DArray& multi_index)
{
  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  const UShort2DArray&  sm_multi_index = csg_driver->smolyak_multi_index();
  size_t i, num_smolyak_indices = sm_multi_index.size();

  // assemble a complete list of individual polynomial coverage
  // defined from the linear combination of mixed tensor products
  multi_index.clear(); tpMultiIndexMap.clear(); tpMultiIndexMapRef.clear();
  tpMultiIndex.resize(num_smolyak_indices);
  if (expConfigOptions.refinementControl ==
      DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
    tpExpansionCoeffs.resize(num_smolyak_indices);
    tpExpansionCoeffGrads.resize(num_smolyak_indices);
  }
  UShortArray exp_order(numVars);
  for (i=0; i<num_smolyak_indices; ++i) {
    // regenerate i-th exp_order as collocKey[i] cannot be used in general case
    // (i.e., for nested rules GP, CC, F2, or GK).  Rather, collocKey[i] is to
    // be used only as the key to the collocation pts.
    sparse_grid_level_to_expansion_order(sm_multi_index[i], exp_order);
    tensor_product_multi_index(exp_order, tpMultiIndex[i]);
    append_multi_index(tpMultiIndex[i], multi_index, tpMultiIndexMap,
		       tpMultiIndexMapRef);
#ifdef DEBUG
    PCout << "level =\n" << sm_multi_index[i] << "expansion_order =\n"
	  << exp_order << "tp_multi_index =\n" << tpMultiIndex[i]
	  << "multi_index =\n" << multi_index << '\n';
#endif // DEBUG
  }

  /*
  case SPARSE_INT_TENSOR_SUM_EXP: case SPARSE_INT_RESTR_TENSOR_SUM_EXP: {
    // assemble a complete list of individual polynomial coverage
    // defined from the linear combination of mixed tensor products
    multi_index.clear();
    UShort2DArray tp_multi_index;
    UShortArray int_order(numVars), exp_order(numVars);
    // Note: restricted rule growth within the sparse grid point set is separate
    // from restricted definition of the expansion terms.  The former makes
    // integrand precision more uniform by delaying exponential sequences, but
    // may still contain some nonuniformity.  The latter may enforce expansion
    // uniformity that would not otherwise be present based on integrand
    // precision alone, in order to reduce the possibility of a response-basis
    // product landing in the concave interior of the integrand resolution.
    short exp_growth = (sparseGridExpansion == SPARSE_INT_RESTR_TENSOR_SUM_EXP)
      ? MODERATE_RESTRICTED_GROWTH : UNRESTRICTED_GROWTH;
    for (i=0; i<num_smolyak_indices; ++i) {
      sparse_grid_level_to_expansion_order(sm_multi_index[i], exp_order,
					   exp_growth);
      tensor_product_multi_index(exp_order, tp_multi_index);
      append_multi_index(tp_multi_index, multi_index);
#ifdef DEBUG
      PCout << "level =\n" << sm_multi_index[i] << "integrand order =\n"
	    << int_order << "expansion order =\n" << exp_order << '\n';
	  //<< "tp_multi_index =\n" << tp_multi_index
	  //<< "multi_index =\n" << multi_index << '\n';
#endif // DEBUG
    }
    break;
  }
  case SPARSE_INT_TOTAL_ORD_EXP: {
    // back out approxOrder & use total_order_multi_index()
    UShortArray quad_order(numVars), integrand_order(numVars);
    UShort2DArray pareto(1), total_pareto;
    for (i=0; i<num_smolyak_indices; ++i) {
      csg_driver->level_to_order(sm_multi_index[i], quad_order);
      quadrature_order_to_integrand_order(quad_order, integrand_order);
      // maintain an n-dimensional Pareto front of nondominated multi-indices
      pareto[0] = integrand_order;
      update_pareto(pareto, total_pareto);
#ifdef DEBUG
      PCout << "level =\n" << sm_multi_index[i] << "\nquad_order =\n"
	    << quad_order << "\nintegrand_order =\n" << integrand_order << '\n';
#endif // DEBUG
    }
#ifdef DEBUG
    PCout << "total_pareto =\n" << total_pareto << '\n';
#endif // DEBUG

    // first pass: compute max isotropic integrand that fits within Pareto front
    unsigned short order = 0;
    integrand_order.assign(numVars, order);
    bool total_order_dominated = true;
    while (total_order_dominated) {
      // calculate all nondominated polynomials for a total-order expansion
      pareto.clear();
      total_order_multi_index(integrand_order, pareto, 0);
      total_order_dominated = assess_dominance(pareto, total_pareto);
#ifdef DEBUG
      PCout << "integrand_order =\n" << integrand_order << "pareto =\n"
	    << pareto << "total_order_dominated = " << total_order_dominated
	    << '\n';
#endif // DEBUG
      // could increment/decrement by 2's due to expansion_order conversion,
      // but the actual resolvable integrand order is typically odd.
      if (total_order_dominated)
	++order; // advance to next test level
      else
	--order; // exiting loop: rewind to last successful
      integrand_order.assign(numVars, order);
    }
#ifdef DEBUG
    PCout << "Isotropic integrand_order =\n" << integrand_order << '\n';
#endif // DEBUG
    integrand_order_to_expansion_order(integrand_order, approxOrder);
    total_order_multi_index(approxOrder, multi_index);
    break;
  }
  case SPARSE_INT_HEUR_TOTAL_ORD_EXP: // early heuristic
    heuristic_sparse_grid_level_to_expansion_order(csg_driver->level(),
						   approxOrder);
    total_order_multi_index(approxOrder, multi_index);
    break;
  }
  */
}


/* This approach reduces memory requirements but must perform additional
   calculation to regenerate the tp_multi_index instances (previously
   generated in sparse_grid_multi_index()).  Currently, these tp_multi_index
   instances are stored in tpMultiIndex for later use in compute_coefficients().
void ProjectOrthogPolyApproximation::
map_tensor_product_multi_index(UShort2DArray& tp_multi_index, size_t tp_index)
{
  const SizetArray& tp_mi_map = tpMultiIndexMap[tp_index];
  size_t i, num_tp_terms = tp_mi_map.size();
  tp_multi_index.resize(num_tp_terms);
  for (i=0; i<num_tp_terms; ++i)
    tp_multi_index[i] = multiIndex[tp_mi_map[i]];
}
*/


/** The optional growth_rate supports the option of forcing the computed
    integrand order to be conservative in the presence of exponential growth
    due to nested quadrature rules.  This avoids aggressive formulation of PCE
    expansion orders when an exponential rule takes a large jump that is not
    balanced by the other index set component mappings.  Note that restricting
    the expansion growth directly from the level (*_RESTRICTED_GROWTH cases
    below used by SPARSE_INT_RESTR_TENSOR_SUM_EXP) is similar but not identical
    to restricting the quadrature order growth from the level and then computing
    the integrand and expansion orders from the restricted quadrature order
    (default UNRESTRICTED_GROWTH case below used by SPARSE_INT_TENSOR_SUM_EXP
    and TENSOR_INT_TENSOR_SUM_EXP, where quadrature rule restriction happens
    elsewhere).  In particular, these approaches differ in granularity of
    control, since the former approach grows linearly and the latter approach
    selects the minimal quadrature order (from nonlinear growth or lookup) that
    meets a linear target. */
void ProjectOrthogPolyApproximation::
sparse_grid_level_to_expansion_order(const UShortArray& level,
				     UShortArray& exp_order)
                                     //,short growth_rate)
{
  size_t n = level.size();
  UShortArray int_order(n);
  //switch (growth_rate) {
  //case UNRESTRICTED_GROWTH: { // used for {SPARSE,TENSOR}_INT_TENSOR_SUM_EXP
    // Best option for TENSOR_INT_TENSOR_SUM_EXP, but SPARSE_INT_TENSOR_SUM_EXP
    // is generally too aggressive for nested rules and exponential growth
    // (SPARSE_INT_RESTR_TENSOR_SUM_EXP is preferred).
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    UShortArray quad_order(n);
    csg_driver->level_to_order(level, quad_order);
    quadrature_order_to_integrand_order(quad_order, int_order);
    //break;
  //}
  /*
  case SLOW_RESTRICTED_GROWTH: // not currently used
    for (size_t i=0; i<n; ++i) // synch with slow linear growth: i = 2l + 1
      int_order[i] =  2*level[i] + 1;
    break;
  case MODERATE_RESTRICTED_GROWTH: // used for SPARSE_INT_RESTR_TENSOR_SUM_EXP
    // mitigate uneven integrand coverage due to exponential rule growth by
    // enforcing moderate linear expansion growth.
    for (size_t i=0; i<n; ++i) // synch with moderate linear growth: i = 4l + 1
      int_order[i] =  4*level[i] + 1;
    break;
  }
  */
  integrand_order_to_expansion_order(int_order, exp_order);
}


void ProjectOrthogPolyApproximation::
quadrature_order_to_integrand_order(const UShortArray& quad_order,
				    UShortArray& int_order)
{
  // Need to know exact polynomial resolution for each mixed tensor grid:
  //   Gaussian integrand resolution:        2m-1
  //   Gauss-Patterson integrand resolution: 2m-1 - previous constraints + 1
  //   Clenshaw-Curtis integrand resolution: m (odd m), m-1 (even m)

  // Burkardt monomial test logic:
  //   sparse_grid_monomial_test: resolve monomials of total degree 2*level + 1
  //   for all rules --> doesn't make sense for exponential growth rules where
  //   order grows faster for Gauss than CC (level_to_order exponential is
  //   2^{w+1}-1 for Gauss and 2^w+1 for CC) --> estimate appears valid for CC
  //   (although it does not define a crisp boundary, since some monomials above
  //   the boundary are resolved) but overly conservative for Gauss (whole
  //   orders above the boundary estimate are resolved).

  size_t i, n = quad_order.size();
  if (int_order.size() != n)
    int_order.resize(n);
  const ShortArray& colloc_rules = driverRep->collocation_rules();
  if (colloc_rules.empty()) // use basisTypes with default modes
    for (i=0; i<n; ++i)
      switch (basisTypes[i]) {
      case CHEBYSHEV_ORTHOG: // default mode is Clenshaw-Curtis
	int_order[i] = (quad_order[i] % 2) ? quad_order[i] : quad_order[i] - 1;
	break;
      default: // default mode is standard non-nested Gauss rules
	int_order[i] =  2*quad_order[i] - 1; // i = 2m - 1
	break;
      }
  else {
    const UShortArray& gk_order = driverRep->genz_keister_order();
    const UShortArray& gk_prec  = driverRep->genz_keister_precision();
    for (i=0; i<n; ++i)
      switch (colloc_rules[i]) {
      case CLENSHAW_CURTIS: case FEJER2:
	// i = m (odd m), m-1 (even m).  Note that growth rule enforces odd.
	// TO DO: verify FEJER2 same as CC
	int_order[i] = (quad_order[i] % 2) ? quad_order[i] : quad_order[i] - 1;
	break;
      case GAUSS_PATTERSON: {
	// for o(l)=2^{l+1}-1, o(l-1) = (o(l)-1)/2
	unsigned short prev_o = std::max(1,(quad_order[i] - 1)/2);
	int_order[i] = 2*quad_order[i] - prev_o;
	break;
      }
      case GENZ_KEISTER: {
	// same relationship as Gauss-Patterson, except prev_o does not follow
	// simple pattern and requires lookup
	unsigned short lev = 0, max_lev = 5;
	for (lev=0; lev<=max_lev; ++lev)
	  if (gk_order[lev] == quad_order[i])
	    { int_order[i] = gk_prec[lev]; break; }
	/*
	int lev, o, prev_o = 1, max_lev = 4, i_rule = GENZ_KEISTER,
	  g_rule = FULL_EXPONENTIAL; // map l->o directly without restriction
	for (lev=0; lev<=max_lev; ++lev) {
	  webbur::level_growth_to_order(1, &lev, &i_rule, &g_rule, &o);
	  if (o == quad_order[i])
	    { int_order[i] = 2*quad_order[i] - prev_o; break; }
	  else
	    prev_o = o;
	}
	*/
	if (lev > max_lev) {
	  PCerr << "Error: maximum GENZ_KEISTER level exceeded in ProjectOrthog"
		<< "PolyApproximation::quadrature_order_to_integrand_order()."
		<< std::endl;
	  abort_handler(-1);
	}
	break;
      }
      default: // standard non-nested Gauss rules
	int_order[i] =  2*quad_order[i] - 1; break; // i = 2m - 1
      }
  }
}


void ProjectOrthogPolyApproximation::
integrand_order_to_expansion_order(const UShortArray& int_order,
				   UShortArray& exp_order)
{
  // reserve half of the integrand order for the expansion and half for the
  // response function (integrand = 2p)
  size_t i, n = int_order.size();
  if (exp_order.size() != n)
    exp_order.resize(n);
  for (i=0; i<n; ++i)
    exp_order[i] = int_order[i] / 2; // remainder truncated
}


void ProjectOrthogPolyApproximation::
append_tensor_expansions(size_t start_index)
{
  // for use in decrement_coefficients()
  prevExpCoeffs = expansionCoeffs; prevExpCoeffGrads = expansionCoeffGrads;

  // update expansion{Coeffs,CoeffGrads} using a hierarchical update
  // rather than building from scratch
  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  const IntArray&     sm_coeffs = csg_driver->smolyak_coefficients();
  const IntArray& sm_coeffs_ref = csg_driver->smolyak_coefficients_reference();
#ifdef DEBUG
  PCout << "In ProjectOrthogPolyApproximation::append_tensor_expansions() with "
	<< "start index " << start_index << "\nsm_coeffs:\n" << sm_coeffs
	<< "sm_coeffs_ref:\n" << sm_coeffs_ref << std::endl;
#endif // DEBUG

  // add trial expansions
  size_t index, num_tensor_grids = sm_coeffs.size();
  int coeff, delta_coeff;
  for (index=start_index; index<num_tensor_grids; ++index) {
    coeff = sm_coeffs[index];
    if (coeff)
      overlay_expansion(tpMultiIndexMap[index], tpExpansionCoeffs[index],
			tpExpansionCoeffGrads[index], coeff);
#ifdef DEBUG
    PCout << "Trial set sm_coeff = " << coeff << "\ntpExpansionCoeffs:\n";
    write_data(PCout, tpExpansionCoeffs[index]);
    PCout << "\ntpMultiIndexMap:\n" << tpMultiIndexMap[index] << '\n';
#endif // DEBUG
  }
  // update other expansion contributions with a changed smolyak coefficient
  for (index=0; index<start_index; ++index) {
    // add new, subtract previous
    delta_coeff = sm_coeffs[index] - sm_coeffs_ref[index];
#ifdef DEBUG
    PCout << "Old set delta_coeff = " << delta_coeff
	  << "\ntpExpansionCoeffs:\n";
    write_data(PCout, tpExpansionCoeffs[index]);
    PCout << "\ntpMultiIndexMap:\n" << tpMultiIndexMap[index] << '\n';
#endif // DEBUG
    if (delta_coeff)
      overlay_expansion(tpMultiIndexMap[index], tpExpansionCoeffs[index],
			tpExpansionCoeffGrads[index], delta_coeff);
  }
}


void ProjectOrthogPolyApproximation::
integration_data(size_t tp_index, SDVArray& tp_data_vars,
		 SDRArray& tp_data_resp, RealVector& tp_weights)
{
  // extract tensor vars/resp from surrData and tensor wts from type1CollocWts1D
  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  const UShortArray&    sm_index = csg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&       key = csg_driver->collocation_key()[tp_index];
  const SizetArray&  colloc_index = csg_driver->collocation_indices()[tp_index];
  const Real3DArray& colloc_wts_1d = csg_driver->type1_collocation_weights_1d();
  const SDVArray& data_vars = surrData.variables_data();
  const SDRArray& data_resp = surrData.response_data();
  size_t i, j, index, num_tp_pts = colloc_index.size();
  tp_data_vars.resize(num_tp_pts); tp_data_resp.resize(num_tp_pts);
  tp_weights.resize(num_tp_pts);
  for (i=0; i<num_tp_pts; ++i) {
    // tensor-product vars/resp
    index = colloc_index[i];
    tp_data_vars[i] = data_vars[index];
    tp_data_resp[i] = data_resp[index];
    // tensor-product weight
    Real& tp_wts_i = tp_weights[i]; tp_wts_i = 1.;
    const UShortArray& key_i = key[i];
    for (j=0; j<numVars; ++j)
      tp_wts_i *= colloc_wts_1d[sm_index[j]][j][key_i[j]];
  }
}


/** The coefficients of the PCE for the response are calculated using a
    spectral projection of the response against each multivariate orthogonal
    polynomial basis fn using the inner product ratio <f,Psi>/<Psi^2>, where
    inner product <a,b> is the n-dimensional integral of a*b*weighting over
    the support range of the n-dimensional (composite) weighting function.
    1-D quadrature rules are defined for specific 1-D weighting functions
    and support ranges and approximate the integral of f*weighting as the
    Sum_i of w_i f_i.  To extend this to n-dimensions, a tensor product
    quadrature rule, cubature, or Smolyak sparse grid rule is applied.  
    It is not necessary to approximate the integral for the denominator
    numerically, since this is available analytically. */
void ProjectOrthogPolyApproximation::
integrate_expansion(const UShort2DArray& multi_index,
		    const SDVArray& data_vars, const SDRArray& data_resp,
		    const RealVector& wt_sets, RealVector& exp_coeffs,
		    RealMatrix& exp_coeff_grads)
{
  // Perform numerical integration via tensor-product quadrature/cubature/
  // Smolyak sparse grids.  Quadrature/cubature use a single application of
  // point and weight sets computed by TensorProductDriver/CubatureDriver, and
  // sparse grids could do this as well, but it is better to integrate the
  // sparse grid on a per-tensor-product basis folowed by summing the
  // corresponding PC expansions.
  if (data_resp[0].is_null()) {
    PCerr << "Error: null SDR in ProjectOrthogPolyApproximation::"
	  << "integrate_expansion()" << std::endl;
    abort_handler(-1);
  }
  size_t i, j, k, num_exp_terms = multi_index.size(),
    num_pts = std::min(data_vars.size(), data_resp.size()),
    num_deriv_vars = data_resp[0].response_gradient().length();
  Real wt_resp_fn_i, Psi_ij; Real* exp_grad;
  RealVector wt_resp_grad_i;
  if (expConfigOptions.expansionCoeffFlag) { // shape if needed and zero out
    if (exp_coeffs.length() != num_exp_terms)
      exp_coeffs.size(num_exp_terms); // init to 0
    else
      exp_coeffs = 0.;
  }
  if (expConfigOptions.expansionCoeffGradFlag) {
    if (exp_coeff_grads.numRows() != num_deriv_vars ||
	exp_coeff_grads.numCols() != num_exp_terms)
      exp_coeff_grads.shape(num_deriv_vars, num_exp_terms); // init to 0
    else
      exp_coeff_grads = 0.;
    wt_resp_grad_i.sizeUninitialized(num_deriv_vars);
  }
  for (i=0; i<num_pts; ++i) {
    if (expConfigOptions.expansionCoeffFlag)
      wt_resp_fn_i = wt_sets[i] * data_resp[i].response_function();
    if (expConfigOptions.expansionCoeffGradFlag) {
      wt_resp_grad_i = data_resp[i].response_gradient(); // copy
      wt_resp_grad_i.scale(wt_sets[i]);
    }
#ifdef DEBUG
    PCout << "wt = " << wt_sets[i] << " resp = "
	  << data_resp[i].response_function() << std::endl;
#endif //DEBUG
    const RealVector& c_vars_i = data_vars[i].continuous_variables();
    for (j=0; j<num_exp_terms; ++j) {
      Psi_ij = multivariate_polynomial(c_vars_i, multi_index[j]);
      if (expConfigOptions.expansionCoeffFlag) {
	exp_coeffs[j] += Psi_ij * wt_resp_fn_i;
#ifdef DEBUG
	PCout << "Psi[" << i << "][" << j << "] = " << Psi_ij
	      << " exp_coeffs[" << j << "] = " << exp_coeffs[j] << std::endl;
#endif //DEBUG
      }
      if (expConfigOptions.expansionCoeffGradFlag) {
	exp_grad = exp_coeff_grads[j];
	for (k=0; k<num_deriv_vars; ++k)
	  exp_grad[k] += Psi_ij * wt_resp_grad_i[k];
      }
    }
  }

  for (i=0; i<num_exp_terms; ++i) {
    Real norm_sq = norm_squared(multi_index[i]);
    if (expConfigOptions.expansionCoeffFlag)
      exp_coeffs[i] /= norm_sq;
    if (expConfigOptions.expansionCoeffGradFlag) {
      exp_grad = exp_coeff_grads[i];
      for (k=0; k<num_deriv_vars; ++k)
	exp_grad[k] /= norm_sq;
    }
  }
#ifdef DEBUG
  PCout << "expansion_coeffs:\n"; write_data(PCout, exp_coeffs);
  if (exp_coeff_grads.numRows()) {
    PCout << "expansion_coeff_grads:\n";
    write_data(PCout, exp_coeff_grads, true, true, true);
  }
  PCout << "\n\n";
#endif // DEBUG
}


/** The coefficients of the PCE for the response are calculated using a
    spectral projection of the response against each multivariate orthogonal
    polynomial basis fn using the inner product ratio <f,Psi>/<Psi^2>,
    where inner product <a,b> is the n-dimensional integral of a*b*weighting
    over the support range of the n-dimensional (composite) weighting
    function.  When interpreting the weighting function as a probability
    density function, <a,b> = expected value of a*b, which can be evaluated
    by sampling from the probability density function and computing the mean
    statistic.  It is not necessary to compute the mean statistic for the
    denominator, since this is available analytically. */
void ProjectOrthogPolyApproximation::expectation()
{
  // "lhs" or "random", no weights needed
  size_t i, j, k, num_surr_data_pts = surrData.size(), num_failed_surr_fn = 0,
    num_failed_surr_grad = 0, num_deriv_vars = expansionCoeffGrads.numRows();
  SizetShortMap::const_iterator fit;
  short                failed_anchor_data = surrData.failed_anchor_data();
  const SizetShortMap& failed_resp_data   = surrData.failed_response_data();
  for (fit=failed_resp_data.begin(); fit!=failed_resp_data.end(); ++fit) {
    if (fit->second & 1) ++num_failed_surr_fn;
    if (fit->second & 2) ++num_failed_surr_grad;
  }
  size_t num_data_pts_fn = num_surr_data_pts - num_failed_surr_fn,
    num_data_pts_grad    = num_surr_data_pts - num_failed_surr_grad,
    num_total_pts_fn = num_data_pts_fn, num_total_pts_grad = num_data_pts_grad;
  bool anchor_fn = false, anchor_grad = false;
  if (surrData.anchor()) {
    if (expConfigOptions.expansionCoeffFlag     && !(failed_anchor_data & 1))
      { anchor_fn   = true; ++num_total_pts_fn; }
    if (expConfigOptions.expansionCoeffGradFlag && !(failed_anchor_data & 2))
      { anchor_grad = true; ++num_total_pts_grad; }
  }
  if (expConfigOptions.expansionCoeffFlag)
    PCout << "Expectations of " << numExpansionTerms << " chaos coefficients "
	  << "using " << num_total_pts_fn << " observations.\n";
  if (expConfigOptions.expansionCoeffGradFlag)
    PCout << "Expectations of gradients of " << numExpansionTerms << " chaos "
	  << "coefficients using " << num_total_pts_grad << " observations.\n";

  /*
  // The following implementation evaluates all PCE coefficients
  // using a consistent expectation formulation
  for (i=0; i<numExpansionTerms; ++i) {
    Real& exp_coeff_i = expansionCoeffs[i];
    exp_coeff_i = (anchor_fn) ?
      surrData.anchor_function() * multivariate_polynomial(
        surrData.anchor_continuous_variables(), multiIndex[i]) : 0.0;
    for (j=0; j<num_data_pts; ++j)
      exp_coeff_i += surrData.response_function(j) * 
        multivariate_polynomial(surrData.continuous_variables(j),multiIndex[i]);
    exp_coeff_i /= num_total_pts * norm_squared(multiIndex[i]);
#ifdef DEBUG
    PCout << "coeff[" << i << "] = " << exp_coeff_i
	  << " norm squared[" << i <<"] = " << norm_squared(multiIndex[i])
	  << '\n';
#endif // DEBUG
  }
  */

  // This alternate implementation evaluates the first PCE coefficient (the
  // response mean) as an expectation and then removes the mean from the
  // expectation evaluation of all subsequent coefficients.  This approach
  // has been observed to result in better results for small sample sizes.
  Real empty_r;
  Real& mean      = (expConfigOptions.expansionCoeffFlag) ?
    expansionCoeffs[0] : empty_r;
  Real* mean_grad = (expConfigOptions.expansionCoeffGradFlag) ?
    expansionCoeffGrads[0] : NULL;
  if (expConfigOptions.expansionCoeffFlag) {
    if (anchor_fn)   mean = surrData.anchor_function();
    else             expansionCoeffs = 0.;
  }
  if (expConfigOptions.expansionCoeffGradFlag) {
    if (anchor_grad) copy_data(surrData.anchor_gradient().values(),
			       num_deriv_vars, mean_grad);
    else             expansionCoeffGrads = 0.;
  }
  for (k=0, fit=failed_resp_data.begin(); k<num_surr_data_pts; ++k) {
    bool add_val  = expConfigOptions.expansionCoeffFlag,
         add_grad = expConfigOptions.expansionCoeffGradFlag;
    fail_booleans(fit, k, add_val, add_grad);
    if (add_val)
      mean += surrData.response_function(k);
    if (add_grad) {
      const RealVector& curr_pt_grad = surrData.response_gradient(k);
      for (j=0; j<num_deriv_vars; ++j)
	mean_grad[j] += curr_pt_grad[j];
    }
  }
  if (expConfigOptions.expansionCoeffFlag)
    mean /= num_total_pts_fn;
  if (expConfigOptions.expansionCoeffGradFlag)
    for (j=0; j<num_deriv_vars; ++j)
      mean_grad[j] /= num_total_pts_grad;

  Real chaos_sample, resp_fn_minus_mean, norm_sq; Real* exp_grad_i;
  RealVector resp_grad_minus_mean;
  if (expConfigOptions.expansionCoeffGradFlag)
    resp_grad_minus_mean.sizeUninitialized(num_deriv_vars);
  if (anchor_fn || anchor_grad) {
    if (anchor_fn)
      resp_fn_minus_mean = surrData.anchor_function() - mean;
    if (anchor_grad) {
      const RealVector& anch_grad = surrData.anchor_gradient();
      for (j=0; j<num_deriv_vars; ++j)
	resp_grad_minus_mean[j] = anch_grad[j] - mean_grad[j];
    }
    const RealVector& c_vars = surrData.anchor_continuous_variables();
    for (i=1; i<numExpansionTerms; ++i) {
      chaos_sample = multivariate_polynomial(c_vars, multiIndex[i]);
      if (anchor_fn)
	expansionCoeffs[i] = resp_fn_minus_mean * chaos_sample;
      if (anchor_grad) {
	exp_grad_i = expansionCoeffGrads[i];
	for (j=0; j<num_deriv_vars; ++j)
	  exp_grad_i[j] = resp_grad_minus_mean[j] * chaos_sample;
      }
    }
  }
  for (k=0, fit=failed_resp_data.begin(); k<num_surr_data_pts; ++k) {
    bool add_val  = expConfigOptions.expansionCoeffFlag,
         add_grad = expConfigOptions.expansionCoeffGradFlag;
    fail_booleans(fit, k, add_val, add_grad);
    if (add_val)
      resp_fn_minus_mean = surrData.response_function(k) - mean;
    if (add_grad) {
      const RealVector& resp_grad = surrData.response_gradient(k);
      for (j=0; j<num_deriv_vars; ++j)
	resp_grad_minus_mean[j] = resp_grad[j] - mean_grad[j];
    }
    const RealVector& c_vars = surrData.continuous_variables(k);
    for (i=1; i<numExpansionTerms; ++i) {
      chaos_sample = multivariate_polynomial(c_vars, multiIndex[i]);
      if (add_val)
	expansionCoeffs[i] += resp_fn_minus_mean * chaos_sample;
      if (add_grad) {
	exp_grad_i = expansionCoeffGrads[i];
	for (j=0; j<num_deriv_vars; ++j)
	  exp_grad_i[j] += resp_grad_minus_mean[j] * chaos_sample;
      }
    }
  }
  for (i=1; i<numExpansionTerms; ++i) {
    norm_sq = norm_squared(multiIndex[i]);
    if (expConfigOptions.expansionCoeffFlag)
      expansionCoeffs[i] /= norm_sq * num_total_pts_fn;
    if (expConfigOptions.expansionCoeffGradFlag) {
      exp_grad_i = expansionCoeffGrads[i];
      for (j=0; j<num_deriv_vars; ++j)
	exp_grad_i[j] /= norm_sq * num_total_pts_grad;
    }
#ifdef DEBUG
    PCout << "coeff[" << i << "] = " << expansionCoeffs[i]
        //<< "coeff_grad[" << i <<"] = " << exp_grad_i
	  << " norm squared[" << i <<"] = " << norm_sq << '\n';
#endif // DEBUG
  }
}


void ProjectOrthogPolyApproximation::
update_pareto(const UShort2DArray& new_pareto, UShort2DArray& total_pareto)
{
  //if (total_pareto.empty())
  //  total_pareto = new_pareto;
  //else {
    size_t i, num_new_p = new_pareto.size(),
      num_total_p = total_pareto.size();
    std::list<UShort2DArray::iterator> rm_iters; // sorted, unique
    UShort2DArray::iterator jit;
    for (i=0; i<num_new_p; ++i) {
      const UShortArray& new_i = new_pareto[i];
      bool new_i_dominated = false;
      for (jit=total_pareto.begin(); jit!=total_pareto.end(); ++jit) {
	bool new_i_dominated_by_j, total_j_dominated;
	assess_dominance(new_i, *jit, new_i_dominated_by_j, total_j_dominated);
	if (new_i_dominated_by_j)
	  { new_i_dominated = true; break; }
	if (total_j_dominated)
	  rm_iters.push_back(jit);
      }
      // 
      // prune newly dominated in reverse order (vector iterators
      // following a point of insertion or deletion are invalidated)
      while (!rm_iters.empty())
	{ total_pareto.erase(rm_iters.back()); rm_iters.pop_back(); }
      // add nondominated
      if (!new_i_dominated)
	total_pareto.push_back(new_i);
    }
  //}
}


bool ProjectOrthogPolyApproximation::
assess_dominance(const UShort2DArray& new_pareto,
		 const UShort2DArray& total_pareto)
{
  bool new_dominated = true;
  size_t i, j, num_new_p = new_pareto.size(),
    num_total_p = total_pareto.size();
  for (i=0; i<num_new_p; ++i) {
    const UShortArray& new_i = new_pareto[i];
    bool new_i_dominated = false;
    for (j=0; j<num_total_p; ++j) {
      bool i_dominated_by_j, j_dominated;
      assess_dominance(new_i, total_pareto[j], i_dominated_by_j, j_dominated);
      if (i_dominated_by_j)
	{ new_i_dominated = true; break; }
    }
    if (!new_i_dominated) {
      new_dominated = false;
#ifdef DEBUG
      PCout << "Nondominated new pareto member =\n" << new_i;
#else
      break;
#endif // DEBUG
    }
  }
  return new_dominated;
}


void ProjectOrthogPolyApproximation::
assess_dominance(const UShortArray& new_order,
		 const UShortArray& existing_order,
		 bool& new_dominated, bool& existing_dominated)
{
  // can't use std::vector::operator< (used for component-wise sorting)
  size_t i, n = new_order.size();
  bool equal = true, existing_dominated_temp = true;
  new_dominated = true;
  for (i=0; i<n; ++i)
    if (new_order[i] > existing_order[i])
      { equal = false; new_dominated = false; }
    else if (existing_order[i] > new_order[i])
      { equal = false; existing_dominated_temp = false; }
  // asymmetric logic since incumbent wins a tie
  existing_dominated = (!equal && existing_dominated_temp);
}


void ProjectOrthogPolyApproximation::
compute_numerical_response_moments(size_t num_moments)
{
  size_t i, num_pts = surrData.size();
  bool anchor_pt = surrData.anchor();
  if (anchor_pt) ++num_pts;

  // define data_coeffs
  RealVector data_coeffs(num_pts);
  if (anchor_pt) {
    data_coeffs[0] = surrData.anchor_function();
    for (i=1; i<num_pts; ++i)
      data_coeffs[i] = surrData.response_function(i-1);
  }
  else
    for (i=0; i<num_pts; ++i)
      data_coeffs[i] = surrData.response_function(i);

  if (storedExpCombineType && !storedExpCoeffs.empty()) {
    // update data_coeffs using evaluations from stored expansions
    switch (storedExpCombineType) {
    case ADD_COMBINE:
      if (anchor_pt) {
	data_coeffs[0] += stored_value(surrData.anchor_continuous_variables());
	for (i=1; i<num_pts; ++i)
	  data_coeffs[i] += stored_value(surrData.continuous_variables(i-1));
      }
      else
	for (i=0; i<num_pts; ++i)
	  data_coeffs[i] += stored_value(surrData.continuous_variables(i));
      break;
    case MULT_COMBINE:
      if (anchor_pt) {
	data_coeffs[0] *= stored_value(surrData.anchor_continuous_variables());
	for (i=1; i<num_pts; ++i)
	  data_coeffs[i] *= stored_value(surrData.continuous_variables(i-1));
      }
      else
	for (i=0; i<num_pts; ++i)
	  data_coeffs[i] *= stored_value(surrData.continuous_variables(i));
      break;
    }
    // stored data may now be cleared
    storedMultiIndex.clear();
    if (expConfigOptions.expansionCoeffFlag)
      storedExpCoeffs.resize(0);
    if (expConfigOptions.expansionCoeffGradFlag)
      storedExpCoeffGrads.reshape(0,0);
  }

  // update numericalMoments based on data_coeffs
  if (numericalMoments.length() != num_moments)
    numericalMoments.sizeUninitialized(num_moments);
  compute_numerical_moments(data_coeffs, driverRep->type1_weight_sets(),
			    numericalMoments);
}


Real ProjectOrthogPolyApproximation::value(const RealVector& x)
{
  // sum expansion to get response value prediction

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    if (storedExpCombineType) // not guaranteed to use tensor indexing
      return OrthogPolyApproximation::value(x);
    else { // Horner's rule approach applicable for tensor indexing
      // Error check for required data
      if (!expConfigOptions.expansionCoeffFlag) {
	PCerr << "Error: expansion coefficients not defined in "
	      << "ProjectOrthogPolyApproximation::value()" << std::endl;
	abort_handler(-1);
      }
      TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
      RealVector accumulator(numVars); // init to 0.
      return tensor_product_value(x, expansionCoeffs, approxOrder, multiIndex,
				  accumulator);
    }
    break;
  /*
  case COMBINED_SPARSE_GRID: {
    // Horner's rule approach requires storage of tpExpansionCoeffs in
    // compute_coefficients().  For now, leave store_tp as is and use
    // default approach if tpExpansionCoeffs is empty.  In addition,
    // tp arrays are not currently updated for expansion combinations.
    if (tpExpansionCoeffs.empty() || storedExpCombineType) // most cases
      return OrthogPolyApproximation::value(x);
    else { // generalized sparse grid case
      // Error check for required data
      if (!expConfigOptions.expansionCoeffFlag) {
	PCerr << "Error: expansion coefficients not defined in "
	      << "ProjectOrthogPolyApproximation::value()" << std::endl;
	abort_handler(-1);
      }
      CombinedSparseGridDriver* csg_driver
	= (CombinedSparseGridDriver*)driverRep;
      const UShort2DArray& sm_mi     = csg_driver->smolyak_multi_index();
      const IntArray&      sm_coeffs = csg_driver->smolyak_coefficients();
      RealVector accumulator(numVars); // init to 0.
      Real approx_val = 0.;
      size_t i, num_sm_mi = sm_mi.size(); int sm_coeff;
      for (i=0; i<num_sm_mi; ++i) {
	sm_coeff = sm_coeffs[i];
	if (sm_coeff)
	  approx_val += sm_coeff *
	    tensor_product_value(x, tpExpansionCoeffs[i],
				 tpApproxOrders[i], // TO DO
				 tpMultiIndex[i], accumulator);
      }
      return approx_val;
    }
    break;
  }
  */
  default: // other cases are total-order expansions
    return OrthogPolyApproximation::value(x);
    break;
  }
}


Real ProjectOrthogPolyApproximation::stored_value(const RealVector& x)
{
  // sum expansion to get response value prediction

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: { // Horner's rule approach
    // Error check for required data
    size_t i, num_stored_terms = storedMultiIndex.size();
    if (!num_stored_terms || storedExpCoeffs.length() != num_stored_terms) {
      PCerr << "Error: stored expansion coefficients not available in "
	    << "ProjectOrthogPolyApproximation::stored_value()" << std::endl;
      abort_handler(-1);
    }
    // Note: requires tensor indexing in storedMultiIndex (see OPA::value(x)),
    // which is safe to assume prior to support of >2 levels of fidelity.
    RealVector accumulator(numVars); // init to 0.
    return tensor_product_value(x, storedExpCoeffs, storedApproxOrder,
				storedMultiIndex, accumulator);
    break;
  }
  // Horner's rule approach would require storage of tensor product components
  //case COMBINED_SPARSE_GRID:
    //break;
  default: // other cases are total-order expansions
    return OrthogPolyApproximation::stored_value(x);
    break;
  }
}


Real ProjectOrthogPolyApproximation::
tensor_product_value(const RealVector& x, const RealVector& tp_coeffs,
		     const UShortArray& approx_order,
		     const UShort2DArray& tp_mi, RealVector& accumulator)
{
  unsigned short ao_0 = approx_order[0], ao_j, mi_i0, mi_ij;
  size_t i, j, num_tp_coeffs = tp_coeffs.length();
  BasisPolynomial& poly_0 = polynomialBasis[0]; Real x0 = x[0];
  for (i=0; i<num_tp_coeffs; ++i) {
    const UShortArray& tp_mi_i = tp_mi[i]; mi_i0 = tp_mi_i[0];
    if (ao_0)
      accumulator[0] += (mi_i0) ? tp_coeffs[i] * poly_0.type1_value(x0, mi_i0)
	                        : tp_coeffs[i];
    else
      accumulator[0]  = tp_coeffs[i];
    if (mi_i0 == ao_0) {
      // accumulate sums over variables with max key value
      for (j=1; j<numVars; ++j) {
	mi_ij = tp_mi_i[j]; ao_j = approx_order[j];
	if (ao_j)
	  accumulator[j] += (mi_ij) ? accumulator[j-1] *
	    polynomialBasis[j].type1_value(x[j], mi_ij) : accumulator[j-1];
	else
	  accumulator[j]  = accumulator[j-1];
	accumulator[j-1] = 0.;
	if (mi_ij != ao_j)
	  break;
      }
    }
  }
  Real tp_val = accumulator[numVars-1];
  accumulator[numVars-1] = 0.;
  return tp_val;
}

} // namespace Pecos
