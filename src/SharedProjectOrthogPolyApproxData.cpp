/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:       SharedProjectOrthogPolyApproxData
//- Description: Implementation code for SharedProjectOrthogPolyApproxData class
//-               
//- Owner:       Mike Eldred

#include "SharedProjectOrthogPolyApproxData.hpp"
#include "TensorProductDriver.hpp"
#include "CombinedSparseGridDriver.hpp"
#include "CubatureDriver.hpp"
#include "pecos_global_defs.hpp"

//#define DEBUG

namespace Pecos {


void SharedProjectOrthogPolyApproxData::allocate_data()
{
  // update_exp_form controls when to update (refinement) and when not to
  // update (subIterator execution) an expansion's multiIndex definition.
  // Simple logic of updating if previous number of points != current number
  // is not robust enough for anisotropic updates --> track using Prev arrays.
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray& quad_order = tpq_driver->quadrature_order();
    bool update_exp_form
      = (quad_order != quadOrderPrev || activeKey != activeKeyPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials?

    if (update_exp_form) {
      UShortArray int_order(numVars);
      quadrature_order_to_integrand_order(driverRep, quad_order, int_order);
      UShortArray& ao = approxOrdIter->second;
      integrand_order_to_expansion_order(int_order, ao);
      UShort2DArray& mi = multiIndexIter->second;
      tensor_product_multi_index(ao, mi); // include upper bound

      // precomputation performed by tpqDriver prior to allocate_data()
      //precompute_maximal_rules(ao);

      allocate_component_sobol(mi);
      quadOrderPrev = quad_order;
      activeKeyPrev = activeKey;
    }

#ifdef DEBUG
    for (i=0; i<numVars; ++i) {
      OrthogonalPolynomial* poly_rep
	= (OrthogonalPolynomial*)polynomialBasis[i].polynomial_rep();
      for (j=1; j<=quad_order[i]; ++j)
	poly_rep->gauss_check(j);
    }
#endif // DEBUG

    PCout << "Orthogonal polynomial approximation order = { ";
    for (size_t i=0; i<numVars; ++i) PCout << approxOrdIter->second[i] << ' ';
    PCout << "} using tensor-product expansion of " << expansion_terms()
	  << " terms\n";
    break;
  }
  case CUBATURE: {
    CubatureDriver* cub_driver = (CubatureDriver*)driverRep;
    //unsigned short cub_int_order = cub_driver->integrand_order();
    //bool update_exp_form
    //  = (cub_int_order != cubIntOrderPrev || activeKey != activeKeyPrev);

    //if (update_exp_form) {
      UShortArray integrand_order(numVars, cub_driver->integrand_order());
      UShortArray& ao = approxOrdIter->second;
      integrand_order_to_expansion_order(integrand_order, ao);
      UShort2DArray& mi = multiIndexIter->second;
      total_order_multi_index(ao, mi);

      // See special logic in CubatureDriver::compute_grid() for GOLUB_WELSCH
      //precompute_maximal_rules(ao);

      allocate_component_sobol(mi);
      //cubIntOrderPrev = cub_int_order; // update reference point
      //activeKeyPrev = activeKey;
    //}

    PCout << "Orthogonal polynomial approximation order = { ";
    for (size_t i=0; i<numVars; ++i) PCout << ao[i] << ' ';
    PCout << "} using total-order expansion of " << expansion_terms()
	  << " terms\n";
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    unsigned short    ssg_level = csg_driver->level();
    const RealVector& aniso_wts = csg_driver->anisotropic_weights();
    bool update_exp_form
      = (ssg_level != ssgLevelPrev || aniso_wts != ssgAnisoWtsPrev ||
	 activeKey != activeKeyPrev || expConfigOptions.refinementControl ==
	 DIMENSION_ADAPTIVE_CONTROL_GENERALIZED);
    // *** TO DO: capture updates to parameterized/numerical polynomials?

    UShort2DArray& mi = multiIndexIter->second;
    if (update_exp_form) {
      sparse_grid_multi_index(csg_driver, mi);

      // precomputation performed by ssgDriver prior to allocate_data()
      //precompute_maximal_rules(multiIndex);

      allocate_component_sobol(mi);
      ssgLevelPrev = ssg_level; ssgAnisoWtsPrev = aniso_wts;
      activeKeyPrev = activeKey;
    }
    PCout << "Orthogonal polynomial approximation level = " << ssg_level
	  << " using tensor integration and tensor sum expansion of "
	  << mi.size() << " terms\n"; break;
    break;
  }
  default: // SAMPLING
    SharedOrthogPolyApproxData::allocate_data();
    break;
  }
}


void SharedProjectOrthogPolyApproxData::increment_data()
{
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: case CUBATURE: {
    // Use same approach as SharedRegressOrthogPolyApproxData --> elevate?
    // *** TO DO ***
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    switch (expConfigOptions.refinementControl) {
    case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED:
      // increment tpMultiIndex{,Map,MapRef} arrays, update tpMultiIndex,
      // update multiIndex and append bookkeeping
      increment_trial_set(csg_driver, multiIndexIter->second);
      break;
    case UNIFORM_CONTROL:
      increment_sparse_grid_multi_index(csg_driver, multiIndexIter->second);
      break;
    }
    break;
  }
  }
  // update Sobol' array sizes to pick up new interaction terms
  increment_component_sobol();
}


void SharedProjectOrthogPolyApproxData::increment_component_sobol()
{
  if (expConfigOptions.vbdFlag && expConfigOptions.vbdOrderLimit != 1) {
    reset_sobol_index_map_values();
    multi_index_to_sobol_index_map(tpMultiIndex[activeKey].back());
    assign_sobol_index_map_values();
  }
}


void SharedProjectOrthogPolyApproxData::decrement_data()
{
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: case CUBATURE: {
    // Use same approach as SharedRegressOrthogPolyApproxData --> elevate?
    // *** TO DO ***
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    switch (expConfigOptions.refinementControl) {
    case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED:
      decrement_trial_set(csg_driver->trial_set(), multiIndexIter->second);
      break;
    case UNIFORM_CONTROL:
      decrement_sparse_grid_multi_index(csg_driver, multiIndexIter->second);
      break;
    }
  }
  }
}


bool SharedProjectOrthogPolyApproxData::push_available()
{
  switch (expConfigOptions.refinementControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
    SparseGridDriver* sg_driver = (SparseGridDriver*)driverRep;
    return push_trial_available(sg_driver->trial_set());
    break;
  }
  //case UNIFORM_CONTROL:  case DIMENSION_ADAPTIVE_CONTROL_SOBOL:
  //case DIMENSION_ADAPTIVE_CONTROL_DECAY:
  default:
    return !poppedTPMultiIndex[activeKey].empty(); break;
  }
}


void SharedProjectOrthogPolyApproxData::pre_push_data()
{
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: case CUBATURE: {
    // Use same approach as SharedRegressOrthogPolyApproxData --> elevate?
    // *** TO DO ***
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    switch (expConfigOptions.refinementControl) {
    case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED:
      pre_push_trial_set(csg_driver->trial_set(), multiIndexIter->second);
      break;
    case UNIFORM_CONTROL:
      push_sparse_grid_multi_index(csg_driver, multiIndexIter->second);
      break;
    }
  }
  }
}


void SharedProjectOrthogPolyApproxData::post_push_data()
{
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: case CUBATURE: {
    // Use same approach as SharedRegressOrthogPolyApproxData --> elevate?
    // *** TO DO ***
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    switch (expConfigOptions.refinementControl) {
    case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED:
      post_push_trial_set(csg_driver->trial_set(), multiIndexIter->second);
      break;
    }
    //case UNIFORM_CONTROL: // no-op
  }
  }
}


void SharedProjectOrthogPolyApproxData::pre_finalize_data()
{
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: case CUBATURE: {
    // Use same approach as SharedRegressOrthogPolyApproxData --> elevate?
    // *** TO DO ***
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    switch (expConfigOptions.refinementControl) {
    case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
      // update multiIndex
      std::deque<UShort2DArray>& popped_tp_mi  = poppedTPMultiIndex[activeKey];
      std::deque<SizetArray>& popped_tp_mi_map
	= poppedTPMultiIndexMap[activeKey];
      std::deque<size_t>& popped_tp_mi_map_ref
	= poppedTPMultiIndexMapRef[activeKey];

      std::deque<UShort2DArray>::iterator iit = popped_tp_mi.begin();
      std::deque<SizetArray>::iterator    mit = popped_tp_mi_map.begin();
      std::deque<size_t>::iterator        rit = popped_tp_mi_map_ref.begin();
      UShort2DArray& mi = multiIndexIter->second;
      for (; iit!=popped_tp_mi.end(); ++iit, ++mit, ++rit)
	append_multi_index(*iit, *mit, *rit, mi);
      // move previous expansion data to current expansion
      UShort3DArray& tp_mi         = tpMultiIndex[activeKey];
      Sizet2DArray&  tp_mi_map     = tpMultiIndexMap[activeKey];
      SizetArray&    tp_mi_map_ref = tpMultiIndexMapRef[activeKey];
      tp_mi.insert(tp_mi.end(), popped_tp_mi.begin(), popped_tp_mi.end());
      tp_mi_map.insert(tp_mi_map.end(), popped_tp_mi_map.begin(),
		       popped_tp_mi_map.end());
      tp_mi_map_ref.insert(tp_mi_map_ref.end(), popped_tp_mi_map_ref.begin(),
			   popped_tp_mi_map_ref.end());
      break;
    }
    case UNIFORM_CONTROL:
      // *** TO DO ***
      break;
    }
    }
  }
}


void SharedProjectOrthogPolyApproxData::post_finalize_data()
{
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: case CUBATURE: {
    // Use same approach as SharedRegressOrthogPolyApproxData --> elevate?
    // *** TO DO ***
    break;
  }
  case COMBINED_SPARSE_GRID:
    switch (expConfigOptions.refinementControl) {
    case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED:
      poppedLevMultiIndex[activeKey].clear();//.erase(activeKey);
      poppedTPMultiIndex[activeKey].clear();//.erase(activeKey);
      poppedTPMultiIndexMap[activeKey].clear();//.erase(activeKey);
      poppedTPMultiIndexMapRef[activeKey].clear();//.erase(activeKey);
      break;
    case UNIFORM_CONTROL:
      // *** TO DO ***
      break;
    }
  }
}


void SharedProjectOrthogPolyApproxData::pre_combine_data()
{
  switch (expConfigOptions.combineType) {
  case ADD_COMBINE:
    // Note: would like to preserve tensor indexing (at least for QUADRATURE
    // case) so that Horner's rule performance opt could be used within
    // tensor_product_value()).  However, a tensor result in the overlay
    // will not occur unless one expansion order dominates the other (partial
    // domination results in sum of tensor expansions as for sparse grids).
    // Therefore, stick with the general-purpose expansion overlay and exclude
    // tensor_product_value() usage for combined coefficient sets.

    // base class version is sufficient; no specialization based on exp form
    SharedOrthogPolyApproxData::pre_combine_data();  break;
  case MULT_COMBINE:
    // compute form of product expansion
    switch (expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE: { // product of two tensor-product expansions
      //active_key(driverRep->maximal_grid());
      // roll up approxOrders to define combinedMultiIndex
      size_t cntr, j, num_seq = approxOrder.size() - 2; // bridge first to last
      if (num_seq) combinedMultiIndexSeq.resize(num_seq);
      std::map<UShortArray, UShortArray>::iterator ao_it = approxOrder.begin();
      UShortArray combined_ao = ao_it->second;   ++ao_it; // copy
      for (cntr=0; ao_it!=approxOrder.end(); ++ao_it) {
	const UShortArray& ao = ao_it->second;
	for (size_t i=0; i<numVars; ++i)
	  combined_ao[i] += ao[i];
	UShort2DArray& combined_mi = (cntr < num_seq) ?
	  combinedMultiIndexSeq[cntr] : combinedMultiIndex;
	tensor_product_multi_index(combined_ao, combined_mi);
      }
      //allocate_component_sobol(combinedMultiIndex); // defer
      break;
    }
    case COMBINED_SPARSE_GRID: { // product of two sums of tensor-product exp.
      //active_key(driverRep->maximal_grid());
      // filter out dominated Smolyak multi-indices that don't contribute
      // to the definition of the product expansion
      CombinedSparseGridDriver* csg_driver
	= (CombinedSparseGridDriver*)driverRep;
      const std::map<UShortArray, UShort2DArray>& sm_mi_map
	= csg_driver->smolyak_multi_index_map();
      std::map<UShortArray, UShort2DArray>::const_iterator sm_it;

      // Define Pareto expansion orders for first level
      sm_it = sm_mi_map.begin();
      const UShort2DArray& sm_mi_1 = sm_it->second;
      size_t i, j, v, num_keys = sm_mi_map.size(), num_sm_mi = sm_mi_1.size();
      UShortArray exp_order, exp_order_prod(numVars);
      UShort2DArray pareto_eo_1, pareto_eo_2, pareto_eo_prod, tp_mi;
      for (i=0; i<num_sm_mi; ++i) {
	sparse_grid_level_to_expansion_order(csg_driver, sm_mi_1[i], exp_order);
	update_pareto_set(exp_order, pareto_eo_1);
      }

      // Define product multi-index arrays across sequence
      size_t num_pareto_1, num_pareto_2, num_pareto_prod, s,
	num_seq = num_keys - 2; // bridge from first to last
      if (num_seq) combinedMultiIndexSeq.resize(num_seq);

      for (s=0; s<=num_seq; ++s) {

	// define Pareto expansion orders for new level
	++sm_it;  const UShort2DArray& sm_mi_2 = sm_it->second;
	num_sm_mi = sm_mi_2.size();  pareto_eo_2.clear();
	for (i=0; i<num_sm_mi; ++i) {
	  sparse_grid_level_to_expansion_order(csg_driver,sm_mi_2[i],exp_order);
	  update_pareto_set(exp_order, pareto_eo_2);
	}

	// enumerate combinations of Pareto expansion orders
        num_pareto_1 = pareto_eo_1.size();  num_pareto_2 = pareto_eo_2.size();
        pareto_eo_prod.clear();
	for (i=0; i<num_pareto_1; ++i) {
	  const UShortArray& eo_1 = pareto_eo_1[i];
	  for (j=0; j<num_pareto_2; ++j) {
	    const UShortArray& eo_2 = pareto_eo_2[j];
	    for (v=0; v<numVars; ++v)
	      exp_order_prod[v] = eo_1[v] + eo_2[v];
	    // refilter products for Pareto exp_order_prod
	    update_pareto_set(exp_order_prod, pareto_eo_prod);
	  }
	}

	// overlay each product expansion from the tensor-product combinations
	UShort2DArray& combined_mi = (s < num_seq) ?
	  combinedMultiIndexSeq[s] : combinedMultiIndex;
	num_pareto_prod = pareto_eo_prod.size();  combined_mi.clear();
	for (i=0; i<num_pareto_prod; ++i) {
	  tensor_product_multi_index(pareto_eo_prod[i], tp_mi);
	  append_multi_index(tp_mi, combined_mi);
	}

	// This is a running product -> subsequent multi-index terms must build
	// on previous results (carry forward Pareto sets of exp_order from
	// prev level for combination with exp_orders for new level).
	if (s < num_seq)
	  pareto_eo_1 = pareto_eo_prod;
      }

      //allocate_component_sobol(combinedMultiIndex); // defer
      break;
    }
    default:
      // base class version supports product of two total-order expansions
      SharedOrthogPolyApproxData::pre_combine_data();  break;
    }
    break;
  case ADD_MULT_COMBINE:
    // base class manages this placeholder
    SharedOrthogPolyApproxData::pre_combine_data();  break;
  }
}


void SharedProjectOrthogPolyApproxData::
sparse_grid_multi_index(CombinedSparseGridDriver* csg_driver,
			UShort2DArray& multi_index)
{
  const UShort2DArray& sm_mi = csg_driver->smolyak_multi_index();
  size_t i, num_smolyak_indices = sm_mi.size();

  // assemble a complete list of individual polynomial coverage
  // defined from the linear combination of mixed tensor products
  multi_index.clear();
  UShort3DArray& tp_mi         = tpMultiIndex[activeKey];
  Sizet2DArray&  tp_mi_map     = tpMultiIndexMap[activeKey];
  SizetArray&    tp_mi_map_ref = tpMultiIndexMapRef[activeKey];
  tp_mi.resize(num_smolyak_indices);
  tp_mi_map.resize(num_smolyak_indices);
  tp_mi_map_ref.resize(num_smolyak_indices);
  UShortArray exp_order(numVars);
  for (i=0; i<num_smolyak_indices; ++i) {
    // regenerate i-th exp_order as collocKey[i] cannot be used in general case
    // (i.e., for nested rules GP, CC, F2, or GK).  Rather, collocKey[i] is to
    // be used only as the key to the collocation pts.
    sparse_grid_level_to_expansion_order(csg_driver, sm_mi[i], exp_order);
    tensor_product_multi_index(exp_order, tp_mi[i]);
    append_multi_index(tp_mi[i], multi_index, tp_mi_map[i], tp_mi_map_ref[i]);
#ifdef DEBUG
    PCout << "level =\n" << sm_mi[i] << "expansion_order =\n" << exp_order
	  << "tp_multi_index =\n" << tp_mi[i] << "multi_index =\n"
	  << multi_index << '\n';
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
      sparse_grid_level_to_expansion_order(csg_driver, sm_multi_index[i],
                                           exp_order,  exp_growth);
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
      quadrature_order_to_integrand_order(driverRep, quad_order,
                                          integrand_order);
      // maintain an n-dimensional Pareto front of nondominated multi-indices
      pareto[0] = integrand_order;
      update_pareto_set(pareto, total_pareto);
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


void SharedProjectOrthogPolyApproxData::
increment_sparse_grid_multi_index(CombinedSparseGridDriver* csg_driver,
				  UShort2DArray& multi_index)
{
  // assemble a complete list of individual polynomial coverage
  // defined from the linear combination of mixed tensor products
  UShort3DArray& tp_mi         = tpMultiIndex[activeKey];
  Sizet2DArray&  tp_mi_map     = tpMultiIndexMap[activeKey];
  SizetArray&    tp_mi_map_ref = tpMultiIndexMapRef[activeKey];
  size_t num_tp_mi = tp_mi.size();

  const UShort2DArray& sm_mi = csg_driver->smolyak_multi_index();
  size_t i, num_smolyak_indices = sm_mi.size();

  tp_mi.resize(num_smolyak_indices);
  tp_mi_map.resize(num_smolyak_indices);
  tp_mi_map_ref.resize(num_smolyak_indices);
  UShortArray exp_order(numVars);
  for (i=num_tp_mi; i<num_smolyak_indices; ++i) {
    // regenerate i-th exp_order as collocKey[i] cannot be used in general case
    // (i.e., for nested rules GP, CC, F2, or GK).  Rather, collocKey[i] is to
    // be used only as the key to the collocation pts.
    sparse_grid_level_to_expansion_order(csg_driver, sm_mi[i], exp_order);
    tensor_product_multi_index(exp_order, tp_mi[i]);
    append_multi_index(tp_mi[i], multi_index, tp_mi_map[i], tp_mi_map_ref[i]);
  }
}


void SharedProjectOrthogPolyApproxData::
decrement_sparse_grid_multi_index(CombinedSparseGridDriver* csg_driver,
				  UShort2DArray& multi_index)
{
  // assemble a complete list of individual polynomial coverage
  // defined from the linear combination of mixed tensor products
  UShort3DArray& tp_mi         = tpMultiIndex[activeKey];
  Sizet2DArray&  tp_mi_map     = tpMultiIndexMap[activeKey];
  SizetArray&    tp_mi_map_ref = tpMultiIndexMapRef[activeKey];
  size_t num_tp_mi = tp_mi.size();

  const UShort2DArray& sm_mi = csg_driver->smolyak_multi_index();
  size_t i, num_smolyak_indices = sm_mi.size();
  std::deque<UShort2DArray>&  pop_tp_mi = poppedTPMultiIndex[activeKey];
  std::deque<SizetArray>& pop_tp_mi_map = poppedTPMultiIndexMap[activeKey];
  std::deque<size_t>& pop_tp_mi_map_ref = poppedTPMultiIndexMapRef[activeKey];
  for (i=num_smolyak_indices; i<num_tp_mi; ++i) {
    pop_tp_mi.push_back(tp_mi[i]);
    pop_tp_mi_map.push_back(tp_mi_map[i]);
    pop_tp_mi_map_ref.push_back(tp_mi_map_ref[i]);
  }

  tp_mi.resize(num_smolyak_indices);         // prune
  tp_mi_map.resize(num_smolyak_indices);     // prune
  tp_mi_map_ref.resize(num_smolyak_indices); // prune
}


void SharedProjectOrthogPolyApproxData::
push_sparse_grid_multi_index(CombinedSparseGridDriver* csg_driver,
			     UShort2DArray& multi_index)
{
  // assemble a complete list of individual polynomial coverage
  // defined from the linear combination of mixed tensor products
  UShort3DArray& tp_mi         = tpMultiIndex[activeKey];
  Sizet2DArray&  tp_mi_map     = tpMultiIndexMap[activeKey];
  SizetArray&    tp_mi_map_ref = tpMultiIndexMapRef[activeKey];

  std::deque<UShort2DArray>&  pop_tp_mi = poppedTPMultiIndex[activeKey];
  std::deque<SizetArray>& pop_tp_mi_map = poppedTPMultiIndexMap[activeKey];
  std::deque<size_t>& pop_tp_mi_map_ref = poppedTPMultiIndexMapRef[activeKey];
  size_t i, num_pop = pop_tp_mi.size();
  for (i=0; i<num_pop; ++i) {
    tp_mi.push_back(pop_tp_mi[i]);
    tp_mi_map.push_back(pop_tp_mi_map[i]);
    tp_mi_map_ref.push_back(pop_tp_mi_map_ref[i]);
  }

  pop_tp_mi.clear();  pop_tp_mi_map.clear();  pop_tp_mi_map_ref.clear();
}


/* This approach reduces memory requirements but must perform additional
   calculation to regenerate the tp_multi_index instances (previously
   generated in sparse_grid_multi_index()).  Currently, these tp_multi_index
   instances are stored in tpMultiIndex for later use in compute_coefficients().
void SharedProjectOrthogPolyApproxData::
map_tensor_product_multi_index(UShort2DArray& tp_multi_index, size_t tp_index)
{
  const SizetArray& tp_mi_map = tpMultiIndexMap[activeKey][tp_index];
  size_t i, num_tp_terms = tp_mi_map.size();
  tp_multi_index.resize(num_tp_terms);
  for (i=0; i<num_tp_terms; ++i)
    tp_multi_index[i] = multiIndex[tp_mi_map[i]];
}
*/


Real SharedProjectOrthogPolyApproxData::
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


/*
Real SharedProjectOrthogPolyApproxData::
tensor_product_value(const RealVector& x, const RealVector& tp_coeffs,
		     const UShortArray& approx_order,
		     const UShort2DArray& tp_mi, RealVector& accumulator)
{
  //PCout << "test\n";
  unsigned short ao_0 = approx_order[0], ao_j, mi_i0, mi_ij;
  size_t i, j, num_tp_coeffs = tp_coeffs.length();
  BasisPolynomial& poly_0 = polynomialBasis[0]; Real x0 = x[0];
  Teuchos::SerialDenseVector<unsigned short,unsigned short> 
    max_order_1d( numVars );
  std::vector< std::set<unsigned short> > orders_1d( numVars );
  for (i=0; i<num_tp_coeffs; ++i) {
    const UShortArray& tp_mi_i = tp_mi[i];
    for (j=0; j<numVars; ++j) {
      max_order_1d[j] = std::max( max_order_1d[j], tp_mi_i[j] );
      orders_1d[j].insert( tp_mi_i[j] );
    }
  }
  std::vector< RealVector > bases_1d( numVars );
  std::set<unsigned short>::iterator it;
  for (j=0; j<numVars; ++j) {
    bases_1d[j].size( max_order_1d[j] );
    for ( it = orders_1d[j].begin(); it != orders_1d[j].end(); ++it )
      bases_1d[j][*it] = polynomialBasis[j].type1_value( x[j], *it );
  }

  for (i=0; i<num_tp_coeffs; ++i) {
    const UShortArray& tp_mi_i = tp_mi[i]; mi_i0 = tp_mi_i[0];
    if (ao_0)
    //accumulator[0] += (mi_i0) ? tp_coeffs[i] * poly_0.type1_value(x0, mi_i0)
      //: tp_coeffs[i];
      accumulator[0] += (mi_i0) ? tp_coeffs[i] * bases_1d[0][mi_i0] :
	tp_coeffs[i];
    else
      accumulator[0]  = tp_coeffs[i];
    if (mi_i0 == ao_0) {
      // accumulate sums over variables with max key value
      for (j=1; j<numVars; ++j) {
	mi_ij = tp_mi_i[j]; ao_j = approx_order[j];
	if (ao_j)
	  accumulator[j] += (mi_ij) ? accumulator[j-1] *
	    //polynomialBasis[j].type1_value(x[j], mi_ij) : accumulator[j-1];
	    bases_1d[j][mi_ij] : accumulator[j-1];
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
*/

} // namespace Pecos
