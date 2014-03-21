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
    case TOTAL_ORDER_BASIS: case ADAPTED_BASIS:
      allocate_total_order(); // defines approxOrder and (candidate) multiIndex
      break;
    case TENSOR_PRODUCT_BASIS:
      inflate_scalar(approxOrder, numVars); // promote scalar->vector, if needed
      tensor_product_multi_index(approxOrder, multiIndex);
      break;
    //case DEFAULT_BASIS: // reassigned in NonDPolynomialChaos ctor
    //  break;
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
    PCout << "} using adapted total-order expansion of "; break;
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


/** Append append_mi to combined_mi, and update append_mi_map and
    append_mi_map_ref to facilitate related aggregations without
    repeated searching. */
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
