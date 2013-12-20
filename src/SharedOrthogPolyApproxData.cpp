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
  // multiIndex is now a shared reference --> sparsity pruning is no longer
  // an issue at this level
  //bool restore_exp_form = (multiIndex.size()!=total_order_terms(approxOrder));

  if (update_exp_form) { //|| restore_exp_form) {
    allocate_total_order(); // defines approxOrder and (candidate) multiIndex
    allocate_component_sobol();
    // Note: defer this if update_exp_form is needed downstream
    approxOrderPrev = approxOrder;
  }

  // output (candidate) expansion form
  PCout << "Orthogonal polynomial approximation order = { ";
  for (size_t i=0; i<numVars; ++i)
    PCout << approxOrder[i] << ' ';
  PCout << "} using total-order expansion of " << multiIndex.size()
	<< " terms\n";
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


void SharedOrthogPolyApproxData::allocate_component_sobol()
{
  if (expConfigOptions.vbdFlag) {
    if (expConfigOptions.vbdOrderLimit == 1) // main effects only
      allocate_main_sobol();
    else { // main + interaction effects
      sobolIndexMap.clear();
      multi_index_to_sobol_index_map(multiIndex);
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


void SharedOrthogPolyApproxData::store_data()
{ storedApproxOrder = approxOrder; storedMultiIndex = multiIndex; }


void SharedOrthogPolyApproxData::combine_data(short combine_type)
{
  // based on incoming combine_type, combine the data stored previously
  // by store_coefficients()

  switch (combine_type) {
  case ADD_COMBINE: {
    // update multiIndex with any storedMultiIndex terms not yet included
    size_t stored_mi_map_ref;
    append_multi_index(storedMultiIndex, multiIndex, storedMultiIndexMap,
		       stored_mi_map_ref);
    break;
  }
  case MULT_COMBINE: {
    // default implementation: product of two total-order expansions
    // (specialized in ProjectOrthogPolyApproximation::combine_coefficients())
    for (size_t i=0; i<numVars; ++i)
      approxOrder[i] += storedApproxOrder[i];
    UShort2DArray multi_index_prod;
    total_order_multi_index(approxOrder, combinedMultiIndex);
    break;
  }
  case ADD_MULT_COMBINE:
    PCerr << "Error : additive+multiplicative combination not yet implemented "
	  << "in SharedOrthogPolyApproxData::combine_data()" << std::endl;
    abort_handler(-1);
    break;
  }
}


/** Append to multi_index based on app_multi_index. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& app_multi_index,
		   UShort2DArray& multi_index)
{
  if (multi_index.empty())
    multi_index = app_multi_index;
  else {
    size_t i, num_app_mi = app_multi_index.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = app_multi_index[i];
      if (std::find(multi_index.begin(), multi_index.end(), search_mi) ==
	  multi_index.end()) // search_mi does not yet exist in multi_index
	multi_index.push_back(search_mi);
    }
  }
}


/** Append to multi_index, multi_index_map, and multi_index_map_ref
    based on app_multi_index. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& app_multi_index,
		   UShort2DArray& multi_index, SizetArray& multi_index_map,
		   size_t& multi_index_map_ref)
{
  size_t i, num_app_mi = app_multi_index.size();
  multi_index_map.resize(num_app_mi);
  if (multi_index.empty()) {
    multi_index = app_multi_index;
    multi_index_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      multi_index_map[i] = i;
  }
  else {
    multi_index_map_ref = multi_index.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = app_multi_index[i];
      size_t index = find_index(multi_index, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	multi_index_map[i] = multi_index.size();
	multi_index.push_back(search_mi);
      }
      else
	multi_index_map[i] = index;
    }
  }
}


/** Append to multi_index based on app_multi_index and previously
    defined multi_index_map and multi_index_map_ref.  If necessary,
    update multi_index_map and multi_index_map_ref. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& app_multi_index,
		   SizetArray& multi_index_map, size_t& multi_index_map_ref,
		   UShort2DArray& multi_index)
{
  if (multi_index.empty())
    multi_index = app_multi_index;
  else {
    size_t i, num_app_mi = app_multi_index.size(), num_mi = multi_index.size();
    if (num_mi == multi_index_map_ref) { // current mi corresponds to ref
      for (i=0; i<num_app_mi; ++i)
	if (multi_index_map[i] >= multi_index_map_ref)
	  multi_index.push_back(app_multi_index[i]);
    }
    else if (num_mi > multi_index_map_ref) { // mi has grown since ref taken
      for (i=0; i<num_app_mi; ++i)
	if (multi_index_map[i] >= multi_index_map_ref) { // previously appended
	  const UShortArray& search_mi = app_multi_index[i];
	  // search from reference pt forward
	  UShort2DArray::iterator it, it_start = multi_index.begin();
	  std::advance(it_start, multi_index_map_ref);
	  it = std::find(it_start, multi_index.end(), search_mi);
	  if (it == multi_index.end()) { // still an append: update map, append
	    multi_index_map[i] = multi_index.size();
	    multi_index.push_back(app_multi_index[i]);
	  }
	  else // no longer an append: only update map
	    multi_index_map[i] = multi_index_map_ref
	                       + std::distance(it_start, it);
	}
      multi_index_map_ref = num_mi; // reference point now updated
    }
    else { // multi_index is not allowed to shrink since ref taken
      PCerr << "Error: multi_index inconsistent with reference size in "
	    << "OrthogPolyApproximation::append_multi_index()." << std::endl;
      abort_handler(-1);
    }
  }
}


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
