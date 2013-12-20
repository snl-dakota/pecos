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
  // no combination by default, even if stored{MultiIndex,ExpCoeffs,
  // ExpCoeffGrads} are defined.  Redefined by attribute passed in
  // combine_coefficients(short).
  storedExpCombineType = NO_COMBINE; // reset to initial state (if needed)

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
      allocate_component_sobol();
      quadOrderPrev = quad_order;
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
    for (size_t i=0; i<numVars; ++i) PCout << approxOrder[i] << ' ';
    PCout << "} using tensor-product expansion of " << multiIndex.size()
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
      allocate_component_sobol();
      //cubIntOrderPrev = cub_int_order; // update reference point
    //}

    PCout << "Orthogonal polynomial approximation order = { ";
    for (size_t i=0; i<numVars; ++i)
      PCout << approxOrder[i] << ' ';
    PCout << "} using total-order expansion of " << multiIndex.size()
	  << " terms\n";
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    unsigned short    ssg_level = csg_driver->level();
    const RealVector& aniso_wts = csg_driver->anisotropic_weights();
    bool update_exp_form
      = (ssg_level != ssgLevelPrev || aniso_wts != ssgAnisoWtsPrev ||
	 expConfigOptions.refinementControl ==
	 DIMENSION_ADAPTIVE_CONTROL_GENERALIZED);
    // *** TO DO: capture updates to parameterized/numerical polynomials?

    if (update_exp_form) {
      sparse_grid_multi_index(multiIndex);
      allocate_component_sobol();
      ssgLevelPrev = ssg_level; ssgAnisoWtsPrev = aniso_wts;
    }
    PCout << "Orthogonal polynomial approximation level = " << ssg_level
	  << " using tensor integration and tensor sum expansion of "
	  << multiIndex.size() << " terms\n"; break;
    break;
  }
  default: // SAMPLING
    SharedOrthogPolyApproxData::allocate_data();
    break;
  }
}


void SharedProjectOrthogPolyApproxData::increment_data()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in SharedProjectOrthogPoly"
	  << "ApproxData::increment_data()" << std::endl;
    abort_handler(-1);
  }

  size_t last_index = tpMultiIndex.size();
  // increment tpMultiIndex{,Map,MapRef} arrays
  UShort2DArray new_us2a; SizetArray new_sa;
  tpMultiIndex.push_back(new_us2a);
  tpMultiIndexMap.push_back(new_sa); tpMultiIndexMapRef.push_back(0);
  // update tpMultiIndex
  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  UShortArray exp_order(numVars);
  sparse_grid_level_to_expansion_order(csg_driver->trial_set(), exp_order);
  tensor_product_multi_index(exp_order, tpMultiIndex[last_index]);
  // update multiIndex and append bookkeeping
  append_multi_index(tpMultiIndex[last_index], multiIndex,
		     tpMultiIndexMap[last_index],
		     tpMultiIndexMapRef[last_index]);
  // update Sobol' array sizes to pick up new interaction terms
  increment_component_sobol();

  // cleanup
  //if (!reEntrantFlag) {
  //  csg_driver->clear_smolyak_arrays();
  //  csg_driver->clear_collocation_arrays();
  //  tpMultiIndex.clear(); tpMultiIndexMap.clear();
  //}
}


void SharedProjectOrthogPolyApproxData::increment_component_sobol()
{
  if (expConfigOptions.vbdFlag && expConfigOptions.vbdOrderLimit != 1) {
    reset_sobol_index_map_values();
    multi_index_to_sobol_index_map(tpMultiIndex.back());
    assign_sobol_index_map_values();
  }
}


void SharedProjectOrthogPolyApproxData::decrement_data()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in SharedProjectOrthogPoly"
	  << "ApproxData::decrement_data()" << std::endl;
    abort_handler(-1);
  }

  // reset multiIndex
  size_t num_exp_terms = tpMultiIndexMapRef.back();
  multiIndex.resize(num_exp_terms); // truncate previous increment

  // reset tensor-product bookkeeping and save restorable data
  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  savedLevMultiIndex.push_back(csg_driver->trial_set());
  savedTPMultiIndex.push_back(tpMultiIndex.back());
  savedTPMultiIndexMap.push_back(tpMultiIndexMap.back());
  savedTPMultiIndexMapRef.push_back(num_exp_terms);

  tpMultiIndex.pop_back();
  tpMultiIndexMap.pop_back();
  tpMultiIndexMapRef.pop_back();
}


void SharedProjectOrthogPolyApproxData::restore_data()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in SharedProjectOrthogPoly"
	  << "ApproxDataP::restore_coefficients()" << std::endl;
    abort_handler(-1);
  }

  // move previous expansion data to current expansion
  size_t last_index = tpMultiIndex.size();
  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  std::deque<UShortArray>::iterator sit
    = std::find(savedLevMultiIndex.begin(), savedLevMultiIndex.end(),
		csg_driver->trial_set());
  restoreIndex = std::distance(savedLevMultiIndex.begin(), sit);
  savedLevMultiIndex.erase(sit);
  std::deque<UShort2DArray>::iterator iit = savedTPMultiIndex.begin();
  std::deque<SizetArray>::iterator    mit = savedTPMultiIndexMap.begin();
  std::deque<size_t>::iterator        rit = savedTPMultiIndexMapRef.begin();
  std::advance(iit, restoreIndex);    std::advance(mit, restoreIndex);
  std::advance(rit, restoreIndex);

  tpMultiIndex.push_back(*iit);          savedTPMultiIndex.erase(iit);
  tpMultiIndexMap.push_back(*mit);       savedTPMultiIndexMap.erase(mit);
  tpMultiIndexMapRef.push_back(*rit);    savedTPMultiIndexMapRef.erase(rit);

  // update multiIndex
  append_multi_index(tpMultiIndex[last_index], tpMultiIndexMap[last_index],
		     tpMultiIndexMapRef[last_index], multiIndex);
}


void SharedProjectOrthogPolyApproxData::finalize_data()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in SharedProjectOrthogPoly"
	  << "ApproxData" << "::finalize_data()" << std::endl;
    abort_handler(-1);
  }

  size_t start_index = tpMultiIndex.size();
  // update multiIndex
  std::deque<UShort2DArray>::iterator iit = savedTPMultiIndex.begin();
  std::deque<SizetArray>::iterator    mit = savedTPMultiIndexMap.begin();
  std::deque<size_t>::iterator        rit = savedTPMultiIndexMapRef.begin();
  for (; iit!=savedTPMultiIndex.end(); ++iit, ++mit, ++rit)
    append_multi_index(*iit, *mit, *rit, multiIndex);
  // move previous expansion data to current expansion
  tpMultiIndex.insert(tpMultiIndex.end(), savedTPMultiIndex.begin(),
    savedTPMultiIndex.end());
  tpMultiIndexMap.insert(tpMultiIndexMap.end(), savedTPMultiIndexMap.begin(),
    savedTPMultiIndexMap.end());
  tpMultiIndexMapRef.insert(tpMultiIndexMapRef.end(),
    savedTPMultiIndexMapRef.begin(), savedTPMultiIndexMapRef.end());
  savedLevMultiIndex.clear();     savedTPMultiIndex.clear();
  savedTPMultiIndexMap.clear();   savedTPMultiIndexMapRef.clear();
}


void SharedProjectOrthogPolyApproxData::store_data()
{
  // Store the aggregated expansion data.  This approach is preferred to
  // appending to savedTP{MultiIndex,Coeffs,CoeffGrads} since the savedTP
  // approach is less general (TP and sum of TP only), less memory
  // efficient (tensor redundancies in sparse grids), and causes ambiguity
  // in finalize_coefficients() for generalized sparse grids.
  storedMultiIndex = multiIndex;

  // approach-specific storage
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: { // sum of tensor-product expansions
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    storedLevMultiIndex = csg_driver->smolyak_multi_index(); break;
  }
  default: // tensor and total-order expansions
    storedApproxOrder = approxOrder;                         break;
  }
}


void SharedProjectOrthogPolyApproxData::combine_data(short combine_type)
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
    SharedOrthogPolyApproxData::combine_data(combine_type);
    break;
  }
  case MULT_COMBINE: {
    // compute form of product expansion
    switch (expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE: { // product of two tensor-product expansions
      for (size_t i=0; i<numVars; ++i)
	approxOrder[i] += storedApproxOrder[i];
      UShort2DArray multi_index_prod;
      tensor_product_multi_index(approxOrder, combinedMultiIndex);
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
	  append_multi_index(tp_multi_index_prod, combinedMultiIndex);
	}
      }
      break;
    }
    default:
      // base class version supports product of two total-order expansions
      SharedOrthogPolyApproxData::combine_data(combine_type);
      break;
    }
    break;
  }
  case ADD_MULT_COMBINE:
    // base class manages this placeholder
    SharedOrthogPolyApproxData::combine_data(combine_type);
    break;
  }
}


void SharedProjectOrthogPolyApproxData::
sparse_grid_multi_index(UShort2DArray& multi_index)
{
  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  const UShort2DArray&  sm_multi_index = csg_driver->smolyak_multi_index();
  size_t i, num_smolyak_indices = sm_multi_index.size();

  // assemble a complete list of individual polynomial coverage
  // defined from the linear combination of mixed tensor products
  multi_index.clear();
  tpMultiIndex.resize(num_smolyak_indices);
  tpMultiIndexMap.resize(num_smolyak_indices);
  tpMultiIndexMapRef.resize(num_smolyak_indices);
  UShortArray exp_order(numVars);
  for (i=0; i<num_smolyak_indices; ++i) {
    // regenerate i-th exp_order as collocKey[i] cannot be used in general case
    // (i.e., for nested rules GP, CC, F2, or GK).  Rather, collocKey[i] is to
    // be used only as the key to the collocation pts.
    sparse_grid_level_to_expansion_order(sm_multi_index[i], exp_order);
    tensor_product_multi_index(exp_order, tpMultiIndex[i]);
    append_multi_index(tpMultiIndex[i], multi_index, tpMultiIndexMap[i],
		       tpMultiIndexMapRef[i]);
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
void SharedProjectOrthogPolyApproxData::
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
void SharedProjectOrthogPolyApproxData::
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


void SharedProjectOrthogPolyApproxData::
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


void SharedProjectOrthogPolyApproxData::
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


void SharedProjectOrthogPolyApproxData::
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


bool SharedProjectOrthogPolyApproxData::
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


void SharedProjectOrthogPolyApproxData::
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
