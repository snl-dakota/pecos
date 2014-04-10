/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedOrthogPolyApproxData
//- Description:  Implementation code for SharedOrthogPolyApproxData class
//-               
//- Owner:        Mike Eldred

#include "SharedOrthogPolyApproxData.hpp"
#include "pecos_global_defs.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG
//#define DECAY_DEBUG

namespace Pecos {


/** This version supports only orthogonal polynomial types.  In this case,
    the polynomial types needed for an orthogonal basis and for computing
    collocation points and weights in an integration driver are the same. */
bool SharedOrthogPolyApproxData::
initialize_orthog_poly_basis_types(const ShortArray& u_types,
				   ShortArray& opb_types)
{
  bool extra_dist_params = false;

  // Initialize opb_types and extra_dist_params from u_types.
  size_t i, num_vars = u_types.size();
  if (opb_types.size() != num_vars)
    opb_types.resize(num_vars);
  for (i=0; i<num_vars; ++i) {
    switch (u_types[i]) {
    case STD_NORMAL:      opb_types[i] = HERMITE_ORTHOG;            break;
    case STD_UNIFORM:     opb_types[i] = LEGENDRE_ORTHOG;           break;
    // weight fn = 1/sqrt(1-x^2); same as BETA/Jacobi for alpha=beta=-1/2
    //case xxx:           opb_types[i] = CHEBYSHEV_ORTHOG;          break;
    case STD_EXPONENTIAL: opb_types[i] = LAGUERRE_ORTHOG;           break;
    case STD_BETA:
      opb_types[i] = JACOBI_ORTHOG;       extra_dist_params = true; break;
    case STD_GAMMA:
      opb_types[i] = GEN_LAGUERRE_ORTHOG; extra_dist_params = true; break;
    default:
      opb_types[i] = NUM_GEN_ORTHOG;      extra_dist_params = true; break;
    }
  }

  return extra_dist_params;
}


void SharedOrthogPolyApproxData::allocate_data()
{
  // detect changes since previous construction
  bool update_exp_form = (approxOrder != approxOrderPrev);
  //bool restore_exp_form = (multiIndex.size() != t*_*_terms(approxOrder));

  if (update_exp_form) { //|| restore_exp_form) {
    switch (expConfigOptions.expBasisType) {
    // ADAPTED starts from total-order basis, like generalized sparse grid
    case ADAPTED_BASIS: {
      // We could assume that:
      // (1) exp_order defines the upper bound of the basis (not the starting
      //     point) --> e.g., we are in 100D and would like to recover terms up
      //     to order 5, but can't form a candidate multiIndex that large. So we
      //     start from SGL 0 and select components up to this upper bnd.  In 
      //     this case, we specify either colloc_points or colloc_ratio << 1.
      // (2) exp_order defines the starting point and there is no explicit
      //     upper bound --> colloc_ratio is more (initially) meaningful.
      // (3) neither: exp_order only used to define colloc_pts from colloc_ratio
      //     and we hard-wire the starting point (e.g., level 0) for adaptation
      // For now, we use case 3 as it is the simplest

      // initialize the sparse grid driver (lightweight mode) for generating
      // candidate index sets
      csgDriver.initialize_grid(numVars, referenceSGLevel);

      // define reference multiIndex and tpMultiIndex{,Map,MapRef} from 
      // initial sparse grid level
      //sparse_grid_multi_index(&csgDriver, multiIndex); // heavyweight mapping
      const UShort2DArray& sm_multi_index = csgDriver.smolyak_multi_index();
      size_t i, num_sm_mi = sm_multi_index.size();
      multiIndex.clear();
      tpMultiIndex.clear(); tpMultiIndexMap.clear(); tpMultiIndexMapRef.clear();
      for (i=0; i<num_sm_mi; ++i)
	increment_trial_set(sm_multi_index[i], multiIndex); // lightwt mapping
      break;
    }
    case DEFAULT_BASIS: // should not occur (reassigned in NonDPCE ctor)
    case TOTAL_ORDER_BASIS:
      allocate_total_order(); // defines approxOrder and (candidate) multiIndex
      break;
    case TENSOR_PRODUCT_BASIS:
      inflate_scalar(approxOrder, numVars); // promote scalar->vector, if needed
      tensor_product_multi_index(approxOrder, multiIndex);
      break;
    }
    allocate_component_sobol(multiIndex);
    // Note: defer this if update_exp_form is needed downstream
    approxOrderPrev = approxOrder;
  }

  // output (candidate) expansion form
  PCout << "Orthogonal polynomial approximation order = { ";
  for (size_t i=0; i<numVars; ++i)
    PCout << approxOrder[i] << ' ';
  switch (expConfigOptions.expBasisType) {
  case ADAPTED_BASIS:
    PCout << "} using adapted expansion initiated from "; break;
  case DEFAULT_BASIS: // should not occur (reassigned in NonDPCE ctor)
  case TOTAL_ORDER_BASIS:
    PCout << "} using total-order expansion of ";         break;
  case TENSOR_PRODUCT_BASIS:
    PCout << "} using tensor-product expansion of ";      break;
  }
  PCout << multiIndex.size() << " terms\n";
}


void SharedOrthogPolyApproxData::allocate_data(const UShort2DArray& mi)
{
  multiIndex = mi;
  allocate_component_sobol(mi);

  // output form of imported expansion
  PCout << "Orthogonal polynomial approximation using imported expansion of "
	<< multiIndex.size() << " terms\n";
}


void SharedOrthogPolyApproxData::allocate_total_order()
{
  // For uniform refinement, all refinements are based off of approxOrder.
  // For PCBDO, approxOrder is invariant and a multiIndex update is prevented
  // by update_exp_form.

  // promote a scalar input into an isotropic vector
  inflate_scalar(approxOrder, numVars);

  // define (candidate) multiIndex
  total_order_multi_index(approxOrder, multiIndex);
}


void SharedOrthogPolyApproxData::
allocate_component_sobol(const UShort2DArray& multi_index)
{
  if (expConfigOptions.vbdFlag) {
    if (expConfigOptions.vbdOrderLimit == 1) // main effects only
      allocate_main_sobol();
    else { // main + interaction effects
      sobolIndexMap.clear();
      multi_index_to_sobol_index_map(multi_index);
      assign_sobol_index_map_values();
      /*
      unsigned short max_order = approxOrder[0];
      size_t v, num_v = sharedDataRep->numVars;
      for (v=1; v<num_v; ++v)
	if (approxOrder[v] > max_order)
	  max_order = approxOrder[v];
      if (max_order >= num_v)	{
	if (sobolIndices.empty())
	  allocate_main_interaction_sobol(num_v); // all n-way interactions
      }
      else {
	bool update_exp_form = (approxOrder != approxOrderPrev);
	if (update_exp_form)
	  allocate_main_interaction_sobol(max_order);
      }
      */
    }
  }
}


void SharedOrthogPolyApproxData::
increment_trial_set(CombinedSparseGridDriver* csg_driver,
		    UShort2DArray& aggregated_mi)
{
  size_t last_index = tpMultiIndex.size();
  // increment tpMultiIndex{,Map,MapRef} arrays
  UShort2DArray new_us2a; SizetArray new_sa;
  tpMultiIndex.push_back(new_us2a);
  tpMultiIndexMap.push_back(new_sa); tpMultiIndexMapRef.push_back(0);
  // update tpMultiIndex
  UShortArray exp_order(numVars);
  sparse_grid_level_to_expansion_order(csg_driver, csg_driver->trial_set(),
				       exp_order);
  tensor_product_multi_index(exp_order, tpMultiIndex[last_index]);
  // update multiIndex and append bookkeeping
  append_multi_index(tpMultiIndex[last_index], aggregated_mi,
		     tpMultiIndexMap[last_index],
		     tpMultiIndexMapRef[last_index]);
}


void SharedOrthogPolyApproxData::
increment_trial_set(const UShortArray& trial_set, UShort2DArray& aggregated_mi)
{
  size_t i, last_index = tpMultiIndex.size();
  // increment tpMultiIndex{,Map,MapRef} arrays
  UShort2DArray new_us2a; SizetArray new_sa;
  tpMultiIndex.push_back(new_us2a);
  tpMultiIndexMap.push_back(new_sa); tpMultiIndexMapRef.push_back(0);
  // update tpMultiIndex
  UShortArray exp_order(numVars);
  // linear growth in Gaussian rules would normally result in a factor of 2:
  //   m = 2l+1 -> i = 2m-1 = 4l+1 -> o = i/2 = 2l
  // This is the default, but finer and coarser grain growth can be used.
  for (i=0; i<numVars; ++i)
    exp_order[i] = multiIndexGrowthFactor * trial_set[i];
  tensor_product_multi_index(exp_order, tpMultiIndex[last_index]);
  // update multiIndex and append bookkeeping
  append_multi_index(tpMultiIndex[last_index], aggregated_mi,
		     tpMultiIndexMap[last_index],
		     tpMultiIndexMapRef[last_index]);
}


void SharedOrthogPolyApproxData::
decrement_trial_set(const UShortArray& trial_set, UShort2DArray& aggregated_mi)
{
  // reset the aggregated multi-index
  size_t num_exp_terms = tpMultiIndexMapRef.back();
  aggregated_mi.resize(num_exp_terms); // truncate previous increment

  // reset tensor-product bookkeeping and save restorable data
  savedLevMultiIndex.push_back(trial_set);
  savedTPMultiIndex.push_back(tpMultiIndex.back());
  savedTPMultiIndexMap.push_back(tpMultiIndexMap.back());
  savedTPMultiIndexMapRef.push_back(num_exp_terms);

  tpMultiIndex.pop_back();
  tpMultiIndexMap.pop_back();
  tpMultiIndexMapRef.pop_back();
}


void SharedOrthogPolyApproxData::
pre_restore_trial_set(const UShortArray& trial_set,
		      UShort2DArray& aggregated_mi)
{
  restoreIndex = find_index(savedLevMultiIndex, trial_set);

  std::deque<UShort2DArray>::iterator iit = savedTPMultiIndex.begin();
  std::deque<SizetArray>::iterator    mit = savedTPMultiIndexMap.begin();
  std::deque<size_t>::iterator        rit = savedTPMultiIndexMapRef.begin();

  size_t last_index = tpMultiIndex.size();
  std::advance(iit, restoreIndex); tpMultiIndex.push_back(*iit);
  std::advance(mit, restoreIndex); tpMultiIndexMap.push_back(*mit);
  std::advance(rit, restoreIndex); tpMultiIndexMapRef.push_back(*rit);

  // update multiIndex
  append_multi_index(tpMultiIndex[last_index], tpMultiIndexMap[last_index],
		     tpMultiIndexMapRef[last_index], aggregated_mi);
}


void SharedOrthogPolyApproxData::
post_restore_trial_set(const UShortArray& trial_set,
		       UShort2DArray& aggregated_mi)
{
  std::deque<UShortArray>::iterator   sit = savedLevMultiIndex.begin();
  std::deque<UShort2DArray>::iterator iit = savedTPMultiIndex.begin();
  std::deque<SizetArray>::iterator    mit = savedTPMultiIndexMap.begin();
  std::deque<size_t>::iterator        rit = savedTPMultiIndexMapRef.begin();
  std::advance(sit, restoreIndex); savedLevMultiIndex.erase(sit);
  std::advance(iit, restoreIndex); savedTPMultiIndex.erase(iit);
  std::advance(mit, restoreIndex); savedTPMultiIndexMap.erase(mit);
  std::advance(rit, restoreIndex); savedTPMultiIndexMapRef.erase(rit);
}


void SharedOrthogPolyApproxData::
restore_best_solution(UShort2DArray& aggregated_mi)
{
  // reset the aggregated multi-index
  aggregated_mi.resize(bestExpTerms); // truncate previous increments

  // reset tensor-product bookkeeping and save restorable data
  savedLevMultiIndex.clear();   savedTPMultiIndex.clear();
  savedTPMultiIndexMap.clear(); savedTPMultiIndexMapRef.clear();

  tpMultiIndex.clear();       // will be rebuilt each time in allocate_data()
  tpMultiIndexMap.clear();    // will be rebuilt each time in allocate_data()
  tpMultiIndexMapRef.clear(); // will be rebuilt each time in allocate_data()
}


void SharedOrthogPolyApproxData::
update_component_sobol(const UShort2DArray& multi_index)
{
  if (expConfigOptions.vbdFlag && expConfigOptions.vbdOrderLimit != 1) {
    reset_sobol_index_map_values();
    multi_index_to_sobol_index_map(multi_index);
    assign_sobol_index_map_values();
  }
}


/** Default storage, specialized in derived classes. */
void SharedOrthogPolyApproxData::store_data()
{ storedApproxOrder = approxOrder; storedMultiIndex = multiIndex; }


void SharedOrthogPolyApproxData::pre_combine_data(short combine_type)
{
  // based on incoming combine_type, combine the data stored previously
  // by store_coefficients()

  switch (combine_type) {
  case ADD_COMBINE: {
    // update multiIndex with any storedMultiIndex terms not yet included.
    // An update in place is sufficient.
    size_t stored_mi_map_ref;
    //append_multi_index(multiIndex, storedMultiIndex, combinedMultiIndex,
    //                   storedMultiIndexMap, stored_mi_map_ref);
    append_multi_index(storedMultiIndex, multiIndex, storedMultiIndexMap,
		       stored_mi_map_ref);
    // update sobolIndexMap with any storedMultiIndex terms not yet included
    update_component_sobol(storedMultiIndex);
    break;
  }
  case MULT_COMBINE: {
    // default implementation: product of two total-order expansions
    // (specialized in SharedProjectOrthogPolyApproxData::pre_combine_data())

    // update approxOrder and define combinedMultiIndex
    for (size_t i=0; i<numVars; ++i)
      approxOrder[i] += storedApproxOrder[i];
    UShort2DArray multi_index_prod;
    total_order_multi_index(approxOrder, combinedMultiIndex);
    // define sobolIndexMap from combinedMultiIndex
    allocate_component_sobol(combinedMultiIndex);
    break;
  }
  case ADD_MULT_COMBINE:
    PCerr << "Error : additive+multiplicative combination not yet implemented "
	  << "in SharedOrthogPolyApproxData::combine_data()" << std::endl;
    abort_handler(-1);
    break;
  }
}


void SharedOrthogPolyApproxData::post_combine_data(short combine_type)
{
  storedApproxOrder.clear();
  storedMultiIndex.clear(); storedMultiIndexMap.clear();

  switch (combine_type) {
  case MULT_COMBINE:
    std::swap(multiIndex, combinedMultiIndex); // pointer swap for efficiency
    combinedMultiIndex.clear();
    break;
  }
}


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
void SharedOrthogPolyApproxData::
sparse_grid_level_to_expansion_order(CombinedSparseGridDriver* csg_driver,
				     const UShortArray& level,
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
    UShortArray quad_order(n);
    csg_driver->level_to_order(level, quad_order);
    quadrature_order_to_integrand_order(csg_driver, quad_order, int_order);
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


void SharedOrthogPolyApproxData::
quadrature_order_to_integrand_order(IntegrationDriver* int_driver,
				    const UShortArray& quad_order,
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
  const ShortArray& colloc_rules = int_driver->collocation_rules();
  if (colloc_rules.empty()) // use orthogPolyTypes with default modes
    for (i=0; i<n; ++i)
      switch (orthogPolyTypes[i]) {
      case CHEBYSHEV_ORTHOG: // default mode is Clenshaw-Curtis
	int_order[i] = (quad_order[i] % 2) ? quad_order[i] : quad_order[i] - 1;
	break;
      default: // default mode is standard non-nested Gauss rules
	int_order[i] =  2*quad_order[i] - 1; // i = 2m - 1
	break;
      }
  else {
    const UShortArray& gk_order = int_driver->genz_keister_order();
    const UShortArray& gk_prec  = int_driver->genz_keister_precision();
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


void SharedOrthogPolyApproxData::
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


/** Append to combined_mi based on append_mi. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& append_mi, UShort2DArray& combined_mi)
{
  if (combined_mi.empty())
    combined_mi = append_mi;
  else {
    size_t i, num_app_mi = append_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      if (std::find(combined_mi.begin(), combined_mi.end(), search_mi) ==
	  combined_mi.end())
	combined_mi.push_back(search_mi);
    }
  }
}


/** Append append_mi to combined_mi, and update append_mi_map
    (SizetArray) and append_mi_map_ref to facilitate related
    aggregations without repeated searching. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& append_mi, UShort2DArray& combined_mi,
		   SizetArray& append_mi_map, size_t& append_mi_map_ref)
{
  size_t i, num_app_mi = append_mi.size();
  append_mi_map.resize(num_app_mi);
  if (combined_mi.empty()) {
    combined_mi = append_mi;
    append_mi_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      append_mi_map[i] = i;
  }
  else {
    append_mi_map_ref = combined_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      size_t index = find_index(combined_mi, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	append_mi_map[i] = combined_mi.size();
	combined_mi.push_back(search_mi);
      }
      else
	append_mi_map[i] = index;
    }
  }
}


/** Append append_mi to combined_mi, and update append_mi_map
    (SizetSet) and append_mi_map_ref to facilitate related
    aggregations without repeated searching. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& append_mi, UShort2DArray& combined_mi,
		   SizetSet& append_mi_map, size_t& append_mi_map_ref)
{
  size_t i, num_app_mi = append_mi.size();
  append_mi_map.clear();
  if (combined_mi.empty()) {
    combined_mi = append_mi;
    append_mi_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      append_mi_map.insert(i);
  }
  else {
    append_mi_map_ref = combined_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      size_t index = find_index(combined_mi, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	append_mi_map.insert(combined_mi.size());
	combined_mi.push_back(search_mi);
      }
      else
	append_mi_map.insert(index);
    }
  }
}


/** Append append_mi to combined_mi, and update append_mi_map (SizetSet)
    and append_mi_map_ref to facilitate related aggregations without
    repeated searching.  This case is used when append_mi and
    combined_mi follow a consistent order without gaps. */
void SharedOrthogPolyApproxData::
append_leading_multi_index(const UShort2DArray& append_mi,
			   UShort2DArray& combined_mi,
			   SizetSet& append_mi_map, size_t& append_mi_map_ref)
{
  size_t i, num_app_mi = append_mi.size();
  append_mi_map.clear();
  if (combined_mi.empty()) {
    combined_mi = append_mi;
    append_mi_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      append_mi_map.insert(i);
  }
  else {
    append_mi_map_ref = combined_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      append_mi_map.insert(i);
      if (i < append_mi_map_ref) {
	// verify that append_mi is a leading subset with consistent ordering
	if (append_mi[i] != combined_mi[i]) {
	  PCerr << "Error: leading subset assumption violated in SharedOrthog"
		<< "PolyApproxData::append_leading_multi_index()." << std::endl;
	  abort_handler(-1);
	}
      }
      else
	combined_mi.push_back(append_mi[i]);
    }
  }
}


/** Append to combined_mi based on append_mi and previously defined
    append_mi_map and append_mi_map_ref.  If necessary, update
    append_mi_map and append_mi_map_ref. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& append_mi, SizetArray& append_mi_map,
		   size_t& append_mi_map_ref,   UShort2DArray& combined_mi)
{
  if (combined_mi.empty())
    combined_mi = append_mi; // assume append_mi_map{,_ref} are up to date
  else {
    size_t i, num_app_mi = append_mi.size(), num_mi = combined_mi.size();
    if (num_mi == append_mi_map_ref) { // current mi corresponds to ref
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref)
	  combined_mi.push_back(append_mi[i]);
    }
    else if (num_mi > append_mi_map_ref) { // mi has grown since ref taken
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref) { // previously appended
	  const UShortArray& search_mi = append_mi[i];
	  // search from reference pt forward
	  UShort2DArray::iterator it, it_start = combined_mi.begin();
	  std::advance(it_start, append_mi_map_ref);
	  it = std::find(it_start, combined_mi.end(), search_mi);
	  if (it == combined_mi.end()) { // still an append: update map, append
	    append_mi_map[i] = combined_mi.size();
	    combined_mi.push_back(append_mi[i]);
	  }
	  else // no longer an append: only update map
	    append_mi_map[i] = append_mi_map_ref + std::distance(it_start, it);
	}
      append_mi_map_ref = num_mi; // reference point now updated
    }
    else { // combined_mi is not allowed to shrink since ref taken
      PCerr << "Error: combined_mi inconsistent with reference size in "
	    << "OrthogPolyApproximation::append_multi_index()." << std::endl;
      abort_handler(-1);
    }
  }
}


/** Append to combined_mi based on append_mi and previously defined
    append_mi_map and append_mi_map_ref.  If necessary, update
    append_mi_map and append_mi_map_ref. */
void SharedOrthogPolyApproxData::
append_multi_index(SizetSet& sparse_indices, const UShort2DArray& append_mi,
		   UShort2DArray& combined_mi, RealVector& exp_coeffs,
		   RealMatrix& exp_coeff_grads)
{
  if (combined_mi.empty())
    combined_mi = append_mi; // sparse indices & exp coeffs are up to date
  else { // merge multi-indices; update sparse_indices and exp coeffs

    bool sparse_append = !sparse_indices.empty(); // empty if over-determined LS
    bool coeff_flag = !exp_coeffs.empty(), grad_flag = !exp_coeff_grads.empty();
    RealVector old_exp_coeffs; RealMatrix old_exp_coeff_grads;
    if (coeff_flag) old_exp_coeffs      = exp_coeffs;
    if (grad_flag)  old_exp_coeff_grads = exp_coeff_grads;
 
    size_t i, combined_index, coeff_index, num_app_mi = append_mi.size(),
      num_coeff = (sparse_append) ? sparse_indices.size() : num_app_mi;
    SizetArray append_mi_map(num_app_mi);
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      combined_index = find_index(combined_mi, search_mi);
      if (combined_index == _NPOS) { // search_mi does not exist in combined_mi
	combined_index = combined_mi.size();
	combined_mi.push_back(search_mi);
      }
      append_mi_map[i] = combined_index;
      if (!sparse_append)
	sparse_indices.insert(combined_index); // becomes resorted
    }

    SizetSet old_sparse_indices; SizetSet::iterator it;
    if (sparse_append) {
      old_sparse_indices = sparse_indices;
      sparse_indices.clear();
      for (it=old_sparse_indices.begin(); it!=old_sparse_indices.end(); ++it)
	sparse_indices.insert(append_mi_map[*it]); // becomes resorted
      it = old_sparse_indices.begin(); // reset for loop to follow
    }

    // now that resorting is completed, reorder exp_coeff{s,_grads} to match
    for (i=0; i<num_coeff; ++i) {
      if (sparse_append) { combined_index = append_mi_map[*it]; ++it; }
      else                 combined_index = append_mi_map[i];
      coeff_index = std::distance(sparse_indices.begin(),
				  sparse_indices.find(combined_index));
      if (coeff_flag) exp_coeffs[coeff_index] = old_exp_coeffs[i];
      if (grad_flag) {
	Real *exp_coeff_grad     = exp_coeff_grads[coeff_index],
	     *old_exp_coeff_grad = old_exp_coeff_grads[i];
	for (size_t j=0; j<numVars; ++j)
	  exp_coeff_grad[j] = old_exp_coeff_grad[j];
      }
    }
  }
}


// The following variants maintain a separation between ref_mi and combined_mi,
// rather than updating in place.  In current use cases, append_mi_map has
// provided sufficient bookkeeping to allow in-place updating.

/*  Create combined_mi by appending append_mi to ref_mi.
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& ref_mi, const UShort2DArray& append_mi,
		   UShort2DArray& combined_mi)
{
  if (ref_mi.empty())
    combined_mi = append_mi;
  else {
    combined_mi = ref_mi;
    size_t i, num_app_mi = append_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      if (std::find(combined_mi.begin(), combined_mi.end(),
		    search_mi) == combined_mi.end())
	combined_mi.push_back(search_mi);
    }
  }
}
*/

/*  Append append_mi to ref_mi to create combined_mi, and update
    append_mi_map and append_mi_map_ref to facilitate related
    aggregations without repeated searching.
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& ref_mi, const UShort2DArray& append_mi,
		   UShort2DArray& combined_mi,  SizetArray& append_mi_map,
		   size_t& append_mi_map_ref)
{
  size_t i, num_app_mi = append_mi.size();
  append_mi_map.resize(num_app_mi);
  if (ref_mi.empty()) {
    combined_mi = append_mi;
    append_mi_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      append_mi_map[i] = i;
  }
  else {
    combined_mi = ref_mi;
    append_mi_map_ref = combined_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = app_mi[i];
      size_t index = find_index(combined_mi, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	append_mi_map[i] = combined_mi.size();
	combined_mi.push_back(search_mi);
      }
      else
	append_mi_map[i] = index;
    }
  }
}
*/

/*  Append append_mi to ref_mi to create combined_mi using previously
    defined append_mi_map and append_mi_map_ref.  If necessary, update
    append_mi_map and append_mi_map_ref.
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& ref_mi, const UShort2DArray& append_mi,
		   SizetArray& append_mi_map, size_t& append_mi_map_ref,
		   UShort2DArray& combined_mi)
{
  if (ref_mi.empty())
    combined_mi = append_mi; // assume append_mi_map{,_ref} are up to date
  else {
    combined_mi = ref_mi;
    size_t i, num_app_mi = append_mi.size(), num_mi = combined_mi.size();
    if (num_mi == append_mi_map_ref) { // current mi corresponds to ref
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref)
	  combined_mi.push_back(append_mi[i]);
    }
    else if (num_mi > append_mi_map_ref) { // mi has grown since ref taken
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref) { // previously appended
	  const UShortArray& search_mi = append_mi[i];
	  // search from reference pt forward
	  UShort2DArray::iterator it, it_start = combined_mi.begin();
	  std::advance(it_start, append_mi_map_ref);
	  it = std::find(it_start, combined_mi.end(), search_mi);
	  if (it == combined_mi.end()) { // still an append: update map, append
	    append_mi_map[i] = combined_mi.size();
	    combined_mi.push_back(append_mi[i]);
	  }
	  else // no longer an append: only update map
	    append_mi_map[i] = append_mi_map_ref + std::distance(it_start, it);
	}
      append_mi_map_ref = num_mi; // reference point now updated
    }
    else { // combined_mi is not allowed to shrink since ref taken
      PCerr << "Error: ref_mi inconsistent with reference size in "
	    << "OrthogPolyApproximation::append_multi_index()." << std::endl;
      abort_handler(-1);
    }
  }
}
*/


/** This test works in combination with DEBUG settings in
    (Legendre,Laguerre,Jacobi,GenLaguerre)OrthogPolynomial::type1_gradient(). */
void SharedOrthogPolyApproxData::gradient_check()
{
  BasisPolynomial hermite_poly(HERMITE_ORTHOG), legendre_poly(LEGENDRE_ORTHOG),
    laguerre_poly(LAGUERRE_ORTHOG), jacobi_poly(JACOBI_ORTHOG),
    gen_laguerre_poly(GEN_LAGUERRE_ORTHOG), chebyshev_poly(CHEBYSHEV_ORTHOG);
  // alpha/beta selections mirror dakota_uq_rosenbrock_pce.in
  jacobi_poly.alpha_stat(1.5);
  jacobi_poly.beta_stat(2.);
  gen_laguerre_poly.alpha_stat(2.5);

  Real x = 0.5; // valid for all support ranges: [-1,1], [0,Inf], [-Inf, Inf]
  PCout << "-------------------------------------------------\n";
  for (size_t n=0; n<=10; n++) {
    PCout << "Gradients at " << x << " for order " << n << '\n';
    hermite_poly.type1_gradient(x, n);
    legendre_poly.type1_gradient(x, n);
    laguerre_poly.type1_gradient(x, n);
    jacobi_poly.type1_gradient(x, n);
    gen_laguerre_poly.type1_gradient(x, n);
    chebyshev_poly.type1_gradient(x, n);
    PCout << "-------------------------------------------------\n";
  }
}


void SharedOrthogPolyApproxData::
get_tag(char* tag, size_t j, unsigned short order) const
{
  switch (orthogPolyTypes[j]) {
  case HERMITE_ORTHOG:
    std::sprintf(tag,  "He%i", order); break;
  case LEGENDRE_ORTHOG:
    std::sprintf(tag,   "P%i", order); break;
  case LAGUERRE_ORTHOG:
    std::sprintf(tag,   "L%i", order); break;
  case JACOBI_ORTHOG:
    std::sprintf(tag, "Pab%i", order); break;
  case GEN_LAGUERRE_ORTHOG:
    std::sprintf(tag,  "La%i", order); break;
  case CHEBYSHEV_ORTHOG:
    std::sprintf(tag,   "T%i", order); break;
  case NUM_GEN_ORTHOG:
    std::sprintf(tag, "Num%i", order); break;
  default:
    PCerr << "Error: bad polynomial type = " << orthogPolyTypes[j]
	  << " in SharedOrthogPolyApproxData::get_tag()." << std::endl;
    abort_handler(-1);
    break; 
  }
}

} // namespace Pecos
