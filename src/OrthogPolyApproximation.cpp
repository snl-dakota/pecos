/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        OrthogPolyApproximation
//- Description:  Implementation code for OrthogPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "OrthogPolyApproximation.hpp"
#include "pecos_global_defs.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG
//#define DECAY_DEBUG


namespace Pecos {

/** This version supports only orthogonal polynomial types.  In this
    case, the polynomial types needed for an orthogonal basis and for
    computing collocation points and weights in an integration driver
    are the same. */
bool OrthogPolyApproximation::
initialize_basis_types(const ShortArray& u_types, ShortArray& basis_types)
{
  bool extra_dist_params = false;

  // Initialize basis_types and extra_dist_params from u_types.
  size_t i, num_vars = u_types.size();
  if (basis_types.size() != num_vars)
    basis_types.resize(num_vars);
  for (i=0; i<num_vars; ++i) {
    switch (u_types[i]) {
    case STD_NORMAL:      basis_types[i] = HERMITE_ORTHOG;            break;
    case STD_UNIFORM:     basis_types[i] = LEGENDRE_ORTHOG;           break;
    // weight fn = 1/sqrt(1-x^2); same as BETA/Jacobi for alpha=beta=-1/2
    //case xxx:           basis_types[i] = CHEBYSHEV_ORTHOG;          break;
    case STD_EXPONENTIAL: basis_types[i] = LAGUERRE_ORTHOG;           break;
    case STD_BETA:
      basis_types[i] = JACOBI_ORTHOG;       extra_dist_params = true; break;
    case STD_GAMMA:
      basis_types[i] = GEN_LAGUERRE_ORTHOG; extra_dist_params = true; break;
    default:
      basis_types[i] = NUM_GEN_ORTHOG;      extra_dist_params = true; break;
    }
  }

  return extra_dist_params;
}


int OrthogPolyApproximation::min_coefficients() const
{ return 0; } // coefficient import case


void OrthogPolyApproximation::allocate_arrays()
{
  // default implementation employs a total-order expansion (needed for
  // PCE import case which instantiates an OrthogPolyApproximation).

  allocate_total_order();
  allocate_total_sobol();
  allocate_component_sobol(); // needs multiIndex from allocate_total_order()

  if (expansionMoments.empty())
    expansionMoments.sizeUninitialized(2);

  // size expansion even if !update_exp_form due to possibility of change
  // to expansion{Coeff,GradFlag} settings
  size_expansion();

  // output expansion form
  PCout << "Orthogonal polynomial approximation order = { ";
  for (size_t i=0; i<numVars; ++i)
    PCout << approxOrder[i] << ' ';
  PCout << "} using total-order expansion of " << numExpansionTerms
	<< " terms\n";
}


void OrthogPolyApproximation::allocate_total_order()
{
  // For uniform refinement, all refinements are based off of approxOrder.
  // For PCBDO, numExpansionTerms and approxOrder are invariant and a
  // multiIndex update is prevented by update_exp_form.

  // promote a scalar input into an isotropic vector
  inflate_scalar(approxOrder, numVars);

  // capture changes due to order increments or sparsity pruning
  bool update_exp_form = (approxOrder       != approxOrderPrev),
      restore_exp_form = (numExpansionTerms != total_order_terms(approxOrder));
  if (update_exp_form || restore_exp_form) {
    total_order_multi_index(approxOrder, multiIndex);
    numExpansionTerms = multiIndex.size();

    // Note: defer this if update_exp_form is needed downstream
    approxOrderPrev = approxOrder;
  }
}


void OrthogPolyApproximation::allocate_component_sobol()
{
  if (expConfigOptions.vbdFlag && expConfigOptions.expansionCoeffFlag) {
    if (expConfigOptions.vbdOrderLimit == 1) // main effects only
      { if (sobolIndices.empty()) allocate_main_sobol(); }
    else { // main + interaction effects
      sobolIndexMap.clear();
      multi_index_to_sobol_index_map(multiIndex);
      sobol_index_map_to_sobol_indices();
      /*
      unsigned short max_order = approxOrder[0]; size_t v;
      for (v=1; v<numVars; ++v)
	if (approxOrder[v] > max_order)
	  max_order = approxOrder[v];
      if (max_order >= numVars)	{
	if (sobolIndices.empty())
	  allocate_main_interaction_sobol(numVars); // all n-way interactions
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


void OrthogPolyApproximation::store_coefficients()
{
  // store the aggregated expansion data
  storedApproxOrder = approxOrder; storedMultiIndex = multiIndex;
  if (expConfigOptions.expansionCoeffFlag)
    storedExpCoeffs = expansionCoeffs;
  if (expConfigOptions.expansionCoeffGradFlag)
    storedExpCoeffGrads = expansionCoeffGrads;
}


void OrthogPolyApproximation::combine_coefficients(short combine_type)
{
  // based on incoming combine_type, combine the data stored previously
  // by store_coefficients()

  switch (combine_type) {
  case ADD_COMBINE: {
    // update multiIndex with any storedMultiIndex terms not yet included
    Sizet2DArray stored_mi_map; SizetArray stored_mi_map_ref;
    append_multi_index(storedMultiIndex, multiIndex, stored_mi_map,
		       stored_mi_map_ref);
    // resize expansion{Coeffs,CoeffGrads} based on updated multiIndex
    resize_expansion();
    // update expansion{Coeffs,CoeffGrads}
    overlay_expansion(stored_mi_map.back(), storedExpCoeffs,
		      storedExpCoeffGrads, 1);
    break;
  }
  case MULT_COMBINE: {
    // default implementation: product of two total-order expansions
    // (specialized in ProjectOrthogPolyApproximation::combine_coefficients())
    for (size_t i=0; i<numVars; ++i)
      approxOrder[i] += storedApproxOrder[i];
    UShort2DArray multi_index_prod;
    total_order_multi_index(approxOrder, multi_index_prod);

    // perform the multiplication of current and stored expansions
    multiply_expansion(storedMultiIndex, storedExpCoeffs, storedExpCoeffGrads,
		       multi_index_prod);
    break;
  }
  case ADD_MULT_COMBINE:
    //overlay_expansion(storedMultiIndex, storedExpCoeffs,
    //                  storedExpCoeffGrads, addCoeffs, addCoeffGrads);
    //multiply_expansion(storedMultiIndex, storedExpCoeffs,
    //                   storedExpCoeffGrads, multCoeffs, multCoeffGrads);
    //compute_combine_factors(addCoeffs, multCoeffs);
    //apply_combine_factors();
    PCerr << "Error : additive+multiplicative combination not yet implemented "
	  << "in OrthogPolyApproximation::combine_coefficients()" << std::endl;
    abort_handler(-1);
    break;
  }

  /* Code below moved to compute_numerical_response_moments()
  storedMultiIndex.clear();
  if (expConfigOptions.expansionCoeffFlag)     storedExpCoeffs.resize(0);
  if (expConfigOptions.expansionCoeffGradFlag) storedExpCoeffGrads.reshape(0,0);
  */

  computedMean = computedVariance = 0;
}


/** Append to multi_index based on app_multi_index. */
void OrthogPolyApproximation::
append_multi_index(const UShort2DArray& app_multi_index,
		   UShort2DArray& multi_index)
{
  if (multi_index.empty())
    multi_index = app_multi_index;
  else {
    size_t i, num_app_mi = app_multi_index.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = app_multi_index[i];
      // TO DO: make this process more efficient
      if (std::find(multi_index.begin(), multi_index.end(), search_mi) ==
	  multi_index.end()) // search_mi does not yet exist in multi_index
	multi_index.push_back(search_mi);
    }
  }
}


/** Append to multi_index, multi_index_map, and multi_index_map_ref
    based on app_multi_index. */
void OrthogPolyApproximation::
append_multi_index(const UShort2DArray& app_multi_index,
		   UShort2DArray& multi_index, Sizet2DArray& multi_index_map,
		   SizetArray& multi_index_map_ref)
{
  size_t i, num_app_mi = app_multi_index.size();
  if (multi_index.empty()) {
    multi_index = app_multi_index;
    multi_index_map_ref.push_back(0);
    multi_index_map.resize(1); multi_index_map[0].resize(num_app_mi);
    for (i=0; i<num_app_mi; ++i)
      multi_index_map[0][i] = i;
  }
  else {
    multi_index_map_ref.push_back(multi_index.size());
    SizetArray app_mi_map(num_app_mi);
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = app_multi_index[i];
      // TO DO: make this process more efficient
      size_t index = find_index(multi_index, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	app_mi_map[i] = multi_index.size();
	multi_index.push_back(search_mi);
      }
      else
	app_mi_map[i] = index;
    }
    multi_index_map.push_back(app_mi_map);
  }
}


/** Append to multi_index based on app_multi_index and previously
    defined multi_index_map and multi_index_map_ref.  If necessary,
    update multi_index_map and multi_index_map_ref. */
void OrthogPolyApproximation::
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


void OrthogPolyApproximation::
overlay_expansion(const SizetArray& multi_index_map,
		  const RealVector& exp_coeffs, const RealMatrix& exp_grads,
		  int coeff)
{
  size_t i, j, index, num_terms = multi_index_map.size(), 
    num_deriv_vars = expansionCoeffGrads.numRows();
  for (i=0; i<num_terms; ++i) {
    index = multi_index_map[i];
    if (expConfigOptions.expansionCoeffFlag)
      expansionCoeffs[index] += coeff * exp_coeffs[i];
    if (expConfigOptions.expansionCoeffGradFlag) {
      Real*       exp_grad_ndx = expansionCoeffGrads[index];
      const Real* grad_i       = exp_grads[i];
      for (j=0; j<num_deriv_vars; ++j)
	exp_grad_ndx[j] += coeff * grad_i[j];
    }
  }
}


void OrthogPolyApproximation::
multiply_expansion(const UShort2DArray& multi_index_b,
		   const RealVector&    exp_coeffs_b,
		   const RealMatrix&    exp_grads_b,
		   const UShort2DArray& multi_index_c)
{
  UShort2DArray multi_index_a = multiIndex;  // copy (both expConfigOptions)
  RealVector exp_coeffs_a = expansionCoeffs; // copy (both expConfigOptions)
  RealMatrix exp_grads_a;
  if (expConfigOptions.expansionCoeffGradFlag)
    exp_grads_a = expansionCoeffGrads;       // copy (CoeffGrads only)
  size_t i, j, k, v, num_a = multi_index_a.size(), num_b = multi_index_b.size(),
    num_c = multi_index_c.size(), num_deriv_vars = exp_grads_a.numRows();

  // precompute 1D basis triple products required
  unsigned short max_a, max_b, max_c; UShortMultiSet max_abc;
  OrthogonalPolynomial* poly_rep_v;
  for (v=0; v<numVars; ++v) {
    max_a = max_b = max_c = 0; max_abc.clear();
    // could track max_abc within combine_coefficients() and pass in, but would
    // need max orders for each dimension for both factors plus their product.
    // Since this would be awkward and only marginally more efficient, just
    // compute them here from the available multi-index arrays.
    for (i=0; i<num_a; ++i)
      if (multi_index_a[i][v] > max_a)
	max_a = multi_index_a[i][v];
    for (i=0; i<num_b; ++i)
      if (multi_index_b[i][v] > max_b)
	max_b = multi_index_b[i][v];
    for (i=0; i<num_c; ++i)
      if (multi_index_c[i][v] > max_c)
	max_c = multi_index_c[i][v];
    max_abc.insert(max_a); max_abc.insert(max_b); max_abc.insert(max_c); 
    poly_rep_v = (OrthogonalPolynomial*)polynomialBasis[v].polynomial_rep();
    poly_rep_v->precompute_triple_products(max_abc);
  }

  // For c = a * b, compute coefficient of product expansion as:
  // \Sum_k c_k \Psi_k = \Sum_i \Sum_j a_i b_j \Psi_i \Psi_j
  //    c_k <\Psi_k^2> = \Sum_i \Sum_j a_i b_j <\Psi_i \Psi_j \Psi_k>
  if (expConfigOptions.expansionCoeffFlag)
    expansionCoeffs.size(num_c);                      // init to 0
  if (expConfigOptions.expansionCoeffGradFlag)
    expansionCoeffGrads.shape(num_deriv_vars, num_c); // init to 0
  Real trip_prod, trip_prod_v, norm_sq_k; bool non_zero;
  for (k=0; k<num_c; ++k) {
    for (i=0; i<num_a; ++i) {
      for (j=0; j<num_b; ++j) {
	trip_prod = 1.;
	for (v=0; v<numVars; ++v) {
	  poly_rep_v=(OrthogonalPolynomial*)polynomialBasis[v].polynomial_rep();
	  non_zero = poly_rep_v->triple_product(multi_index_a[i][v],
	    multi_index_b[j][v], multi_index_c[k][v], trip_prod_v);
	  if (non_zero) trip_prod *= trip_prod_v;
	  else          break;
	}
	if (non_zero) {
	  if (expConfigOptions.expansionCoeffFlag)
	    expansionCoeffs[k] += exp_coeffs_a[i] * exp_coeffs_b[j] * trip_prod;
	  if (expConfigOptions.expansionCoeffGradFlag)
	    for (v=0; v<num_deriv_vars; ++v)
	      expansionCoeffGrads(v,k) += (exp_coeffs_a[i] * exp_grads_b(v,j)
		+ exp_coeffs_b[j] * exp_grads_a(v,i)) * trip_prod;
	}
      }
    }
    norm_sq_k = norm_squared(multi_index_c[k]);
    if (expConfigOptions.expansionCoeffFlag)
      expansionCoeffs[k] /= norm_sq_k;
    if (expConfigOptions.expansionCoeffGradFlag)
      for (v=0; v<num_deriv_vars; ++v)
	expansionCoeffGrads(v,k) /= norm_sq_k;
  }
  multiIndex = multi_index_c; numExpansionTerms = multiIndex.size();
}


void OrthogPolyApproximation::integration_checks()
{
  if (surrData.anchor()) {
    PCerr << "Error: anchor point not supported for numerical integration in "
	  << "OrthogPolyApproximation::integration()." << std::endl;
    abort_handler(-1);
  }
  if (!driverRep) {
    PCerr << "Error: pointer to integration driver required in "
	  << "OrthogPolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }
  size_t num_data_pts = surrData.size(), num_grid_pts = driverRep->grid_size();
  if (num_data_pts != num_grid_pts) {
    PCerr << "Error: number of current points (" << num_data_pts << ") is "
	  << "not consistent with\n       number of points/weights ("
	  << num_grid_pts << ") from integration driver in\n       "
	  << "OrthogPolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }
}


Real OrthogPolyApproximation::value(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  Real approx_val = 0.;
  for (size_t i=0; i<numExpansionTerms; ++i)
    approx_val += expansionCoeffs[i] * multivariate_polynomial(x,multiIndex[i]);
  return approx_val;
}


const RealVector& OrthogPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  // could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != numVars)
    approxGradient.size(numVars); // init to 0
  else
    approxGradient = 0.;

  // sum expansion to get response gradient prediction
  size_t i, j;
  for (i=0; i<numExpansionTerms; ++i) {
    const RealVector& term_i_grad
      = multivariate_polynomial_gradient_vector(x, multiIndex[i]);
    Real& coeff_i = expansionCoeffs[i];
    for (j=0; j<numVars; ++j)
      approxGradient[j] += coeff_i * term_i_grad[j];
  }
  return approxGradient;
}


const RealVector& OrthogPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t i, j, num_deriv_vars = dvv.size();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.size(num_deriv_vars); // init to 0
  else
    approxGradient = 0.;

  // sum expansion to get response gradient prediction
  for (i=0; i<numExpansionTerms; ++i) {
    const RealVector& term_i_grad
      = multivariate_polynomial_gradient_vector(x, multiIndex[i], dvv);
    Real& coeff_i = expansionCoeffs[i];
    for (j=0; j<num_deriv_vars; ++j)
      approxGradient[j] += coeff_i * term_i_grad[j];
  }
  return approxGradient;
}


const RealVector& OrthogPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in OrthogPoly"
	  << "Approximation::gradient_coefficient_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t i, j, num_deriv_vars = expansionCoeffGrads.numRows();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.size(num_deriv_vars); // init to 0
  else
    approxGradient = 0.;

  // sum expansion to get response gradient prediction
  for (i=0; i<numExpansionTerms; ++i) {
    Real term_i = multivariate_polynomial(x, multiIndex[i]);
    const Real* exp_coeff_grad_i = expansionCoeffGrads[i];
    for (j=0; j<num_deriv_vars; ++j)
      approxGradient[j] += exp_coeff_grad_i[j] * term_i;
  }
  return approxGradient;
}


Real OrthogPolyApproximation::stored_value(const RealVector& x)
{
  // Error check for required data
  size_t i, num_stored_terms = storedMultiIndex.size();
  if (!num_stored_terms || storedExpCoeffs.length() != num_stored_terms) {
    PCerr << "Error: stored expansion coefficients not available in "
	  << "OrthogPolyApproximation::stored_value()" << std::endl;
    abort_handler(-1);
  }

  Real approx_val = 0.;
  for (size_t i=0; i<num_stored_terms; ++i)
    approx_val += storedExpCoeffs[i]
               *  multivariate_polynomial(x, storedMultiIndex[i]);
  return approx_val;
}


const RealVector& OrthogPolyApproximation::
stored_gradient_basis_variables(const RealVector& x)
{
  // Error check for required data
  size_t i, j, num_stored_terms = storedMultiIndex.size();
  if (!num_stored_terms || storedExpCoeffs.length() != num_stored_terms) {
    PCerr << "Error: stored expansion coefficients not available in OrthogPoly"
	  << "Approximation::stored_gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != numVars)
    approxGradient.size(numVars); // init to 0
  else
    approxGradient = 0.;

  // sum expansion to get response gradient prediction
  for (i=0; i<num_stored_terms; ++i) {
    const RealVector& term_i_grad
      = multivariate_polynomial_gradient_vector(x, storedMultiIndex[i]);
    Real& coeff_i = storedExpCoeffs[i];
    for (j=0; j<numVars; ++j)
      approxGradient[j] += coeff_i * term_i_grad[j];
  }
  return approxGradient;
}


const RealVector& OrthogPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x)
{
  // Error check for required data
  size_t i, j, num_stored_terms = storedMultiIndex.size(),
    num_deriv_vars = storedExpCoeffGrads.numRows();
  if (!num_stored_terms || storedExpCoeffGrads.numCols() != num_stored_terms) {
    PCerr << "Error: stored expansion coeff grads not available in OrthogPoly"
	  << "Approximation::stored_gradient_nonbasis_variables()" << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != num_deriv_vars)
    approxGradient.size(num_deriv_vars); // init to 0
  else
    approxGradient = 0.;

  // sum expansion to get response gradient prediction
  for (i=0; i<num_stored_terms; ++i) {
    Real term_i = multivariate_polynomial(x, storedMultiIndex[i]);
    const Real* coeff_grad_i = storedExpCoeffGrads[i];
    for (j=0; j<num_deriv_vars; ++j)
      approxGradient[j] += coeff_grad_i[j] * term_i;
  }
  return approxGradient;
}


/** In this case, all expansion variables are random variables and the
    mean of the expansion is simply the first chaos coefficient. */
Real OrthogPolyApproximation::mean()
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  bool std_mode = nonRandomIndices.empty();
  if (std_mode && (computedMean & 1))
    return expansionMoments[0];

  Real mean = expansionCoeffs[0];
  if (std_mode)
    { expansionMoments[0] = mean; computedMean |= 1; }
  return mean;
}


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves evaluating the
    expectation over this subset. */
Real OrthogPolyApproximation::mean(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  bool all_mode = !nonRandomIndices.empty();
  if (all_mode && (computedMean & 1) && match_nonrandom_vars(x, xPrevMean))
    return expansionMoments[0];

  Real mean = expansionCoeffs[0];
  for (size_t i=1; i<numExpansionTerms; ++i)
    // expectations are zero for expansion terms with nonzero random indices
    if (zero_random(multiIndex[i])) {
      mean += expansionCoeffs[i]
	   *  multivariate_polynomial(x, multiIndex[i], nonRandomIndices);
#ifdef DEBUG
      PCout << "Mean estimate inclusion: term index = " << i << " Psi = "
	    << multivariate_polynomial(x, multiIndex[i], nonRandomIndices)
	    << " PCE coeff = " << expansionCoeffs[i] << " total = " << mean
	    << std::endl;
#endif // DEBUG
    }

  if (all_mode)
    { expansionMoments[0] = mean; computedMean |= 1; xPrevMean = x; }
  return mean;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  In
    this case, the derivative of the expectation is the expectation of
    the derivative.  The mixed derivative case (some design variables
    are inserted and some are augmented) requires no special treatment. */
const RealVector& OrthogPolyApproximation::mean_gradient()
{
  // d/ds \mu_R = d/ds \alpha_0 = <dR/ds>

  // Error check for required data
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in "
	  << "OrthogPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  bool std_mode = nonRandomIndices.empty();
  if (std_mode && (computedMean & 2))
    return meanGradient;

  meanGradient = Teuchos::getCol(Teuchos::Copy, expansionCoeffGrads, 0);
  if (std_mode) computedMean |=  2; //   activate 2-bit
  else          computedMean &= ~2; // deactivate 2-bit: protect mixed usage
  return meanGradient;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  In this case, the mean of the expansion is the
    expectation over the random subset and the derivative of the mean
    is the derivative of the remaining expansion over the non-random
    subset.  This function must handle the mixed case, where some
    design/state variables are augmented (and are part of the
    expansion: derivatives are evaluated as described above) and some
    are inserted (derivatives are obtained from expansionCoeffGrads). */
const RealVector& OrthogPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  // if already computed, return previous result
  bool all_mode = !nonRandomIndices.empty();
  if ( all_mode && (computedMean & 2) &&
       match_nonrandom_vars(x, xPrevMeanGrad) ) // && dvv == dvvPrev)
    return meanGradient;

  size_t i, j, deriv_index, num_deriv_vars = dvv.size(),
    cntr = 0; // insertions carried in order within expansionCoeffGrads
  if (meanGradient.length() != num_deriv_vars)
    meanGradient.sizeUninitialized(num_deriv_vars);
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    Real& grad_i = meanGradient[i];
    if (randomVarsKey[deriv_index]) { // deriv w.r.t. design var insertion
      // Error check for required data
      if (!expConfigOptions.expansionCoeffGradFlag) {
	PCerr << "Error: expansion coefficient gradients not defined in "
	      << "OrthogPolyApproximation::mean_gradient()." << std::endl;
	abort_handler(-1);
      }
      grad_i = expansionCoeffGrads[0][cntr];
    }
    else {
      grad_i = 0.;
      if (!expConfigOptions.expansionCoeffFlag) { // check for reqd data
	PCerr << "Error: expansion coefficients not defined in "
	      << "OrthogPolyApproximation::mean_gradient()" << std::endl;
	abort_handler(-1);
      }
    }
    for (j=1; j<numExpansionTerms; ++j) {
      // expectations are zero for expansion terms with nonzero random indices
      if (zero_random(multiIndex[j])) {
	// In both cases below, term to differentiate is alpha_j(s) Psi_j(s)
	// since <Psi_j>_xi = 1 for included terms.  The difference occurs
	// based on whether a particular s_i dependence appears in alpha
	// (for inserted) or Psi (for augmented).
	if (randomVarsKey[deriv_index])
	  // -------------------------------------------
	  // derivative w.r.t. design variable insertion
	  // -------------------------------------------
	  grad_i += expansionCoeffGrads[j][cntr]
	    * multivariate_polynomial(x, multiIndex[j], nonRandomIndices);
	else
	  // ----------------------------------------------
	  // derivative w.r.t. design variable augmentation
	  // ----------------------------------------------
	  grad_i += expansionCoeffs[j] *
	    multivariate_polynomial_gradient(x, deriv_index, multiIndex[j],
					     nonRandomIndices);
      }
    }
    if (randomVarsKey[deriv_index]) // derivative w.r.t. design var insertion
      ++cntr;
  }
  if (all_mode) { computedMean |=  2; xPrevMeanGrad = x; }
  else            computedMean &= ~2; // deactivate 2-bit: protect mixed usage
  return meanGradient;
}


Real OrthogPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  OrthogPolyApproximation* opa_2 = (OrthogPolyApproximation*)poly_approx_2;
  bool same = (opa_2 == this), std_mode = nonRandomIndices.empty();

  // Error check for required data
  if ( !expConfigOptions.expansionCoeffFlag ||
       ( !same && !opa_2->expConfigOptions.expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if (same) {
    if (std_mode && (computedVariance & 1))
      return expansionMoments[1];
    else {
      Real covar = 0.;
      for (size_t i=1; i<numExpansionTerms; ++i)
	covar += expansionCoeffs[i] * expansionCoeffs[i]
	      *  norm_squared(multiIndex[i]);
      if (std_mode)
	{ expansionMoments[1] = covar; computedVariance |= 1; }
      return covar;
    }
  }
  else {
    const RealVector& exp_coeffs_2 = opa_2->expansionCoeffs;
    Real covar = 0.;
    if (sparseIndices.empty()) // dense indexing is consistent
      for (size_t i=1; i<numExpansionTerms; ++i)
	covar += expansionCoeffs[i] * exp_coeffs_2[i]
	      *  norm_squared(multiIndex[i]);
    else { // for sparse PCE, multiIndex sequences may differ
      const SizetSet& sparse_ind_2 = opa_2->sparseIndices;
      SizetSet::const_iterator cit1 = ++sparseIndices.begin(),
	cit2 = ++sparse_ind_2.begin();
      size_t si1, si2, i1 = 1, i2 = 1;
      while (cit1 != sparseIndices.end() && cit2 != sparse_ind_2.end()) {
	si1 = *cit1; si2 = *cit2;
	if (si1 == si2) {// equality in sparse index implies multiIndex equality
	  covar += expansionCoeffs[i1] * exp_coeffs_2[i2]
	        *  norm_squared(multiIndex[i1]); // or multi_index_2[i2]
	  ++cit1; ++cit2; ++i1; ++i2;
	}
	else if (si1 < si2) { ++cit1; ++i1; }
	else                { ++cit2; ++i2; }
      }
    }
    return covar;
  }
}


Real OrthogPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  OrthogPolyApproximation* opa_2 = (OrthogPolyApproximation*)poly_approx_2;
  bool same = (this == opa_2), all_mode = !nonRandomIndices.empty();

  // Error check for required data
  if ( !expConfigOptions.expansionCoeffFlag ||
       ( !same && !opa_2->expConfigOptions.expansionCoeffFlag )) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if ( same && all_mode && (computedVariance & 1) &&
       match_nonrandom_vars(x, xPrevVar) )
    return expansionMoments[1];

  const RealVector&     exp_coeffs_2 = opa_2->expansionCoeffs;
  const UShort2DArray& multi_index_2 = opa_2->multiIndex;
  size_t i1, i2, num_i1 = multiIndex.size(), num_i2 = multi_index_2.size();
  Real covar = 0.;
  for (i1=1; i1<num_i1; ++i1) {
    // For r = random_vars and nr = non_random_vars,
    // sigma^2_R(nr) = < (R(r,nr) - \mu_R(nr))^2 >_r
    // -> only include terms from R(r,nr) which don't appear in \mu_R(nr)
    if (!zero_random(multiIndex[i1])) {
      Real norm_sq_i = norm_squared(multiIndex[i1], randomIndices);
      for (i2=1; i2<num_i2; ++i2)
	// random polynomial part must be identical to contribute to variance
	// (else orthogonality drops term).  Note that it is not necessary to
	// collapse terms with the same random basis subset, since cross term
	// in (a+b)(a+b) = a^2+2ab+b^2 gets included.  If terms were collapsed
	// (following eval of non-random portions), the nested loop could be
	// replaced with a single loop to evaluate (a+b)^2.
	if (match_random_key(multiIndex[i1], multi_index_2[i2]))
	  covar += expansionCoeffs[i1]  * exp_coeffs_2[i2] * norm_sq_i *
	    multivariate_polynomial(x,    multiIndex[i1], nonRandomIndices) *
	    multivariate_polynomial(x, multi_index_2[i2], nonRandomIndices);
    }
  }
  if (same && all_mode)
    { expansionMoments[1] = covar; computedVariance |= 1; xPrevVar = x; }
  return covar;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  The
    mixed derivative case (some design variables are inserted and some
    are augmented) requires no special treatment. */
const RealVector& OrthogPolyApproximation::variance_gradient()
{
  // d/ds \sigma^2_R = Sum_{j=1}^P <Psi^2_j> d/ds \alpha^2_j
  //                 = 2 Sum_{j=1}^P \alpha_j <dR/ds, Psi_j>

  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag ||
      !expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: insufficient expansion coefficient data in "
	  << "OrthogPolyApproximation::variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  bool std_mode = nonRandomIndices.empty();
  if (std_mode && (computedVariance & 2))
    return varianceGradient;

  size_t i, j, num_deriv_vars = expansionCoeffGrads.numRows();
  if (varianceGradient.length() != num_deriv_vars)
    varianceGradient.sizeUninitialized(num_deriv_vars);
  varianceGradient = 0.;
  for (i=1; i<numExpansionTerms; ++i) {
    Real term_i = 2. * expansionCoeffs[i] * norm_squared(multiIndex[i]);
    for (j=0; j<num_deriv_vars; ++j)
      varianceGradient[j] += term_i * expansionCoeffGrads[i][j];
  }
  if (std_mode) computedVariance |=  2;
  else          computedVariance &= ~2; // deactivate 2-bit: protect mixed usage
  return varianceGradient;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  This function must handle the mixed case, where
    some design/state variables are augmented (and are part of the
    expansion) and some are inserted (derivatives are obtained from
    expansionCoeffGrads). */
const RealVector& OrthogPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::variance_gradient()" << std::endl;
    abort_handler(-1);
  }

  // if already computed, return previous result
  bool all_mode = !nonRandomIndices.empty();
  if ( all_mode && (computedVariance & 2) &&
       match_nonrandom_vars(x, xPrevVarGrad) ) // && dvv == dvvPrev)
    return varianceGradient;

  size_t i, j, k, deriv_index, num_deriv_vars = dvv.size(),
    cntr = 0; // insertions carried in order within expansionCoeffGrads
  if (varianceGradient.length() != num_deriv_vars)
    varianceGradient.sizeUninitialized(num_deriv_vars);
  varianceGradient = 0.;

  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (randomVarsKey[deriv_index] && !expConfigOptions.expansionCoeffGradFlag){
      PCerr << "Error: expansion coefficient gradients not defined in "
	    << "OrthogPolyApproximation::variance_gradient()." << std::endl;
      abort_handler(-1);
    }
    for (j=1; j<numExpansionTerms; ++j) {
      if (!zero_random(multiIndex[j])) {
	Real norm_sq_j = norm_squared(multiIndex[j], randomIndices);
	for (k=1; k<numExpansionTerms; ++k) {
	  // random part of polynomial must be identical to contribute to
	  // variance (else orthogonality drops term)
	  if (match_random_key(multiIndex[j], multiIndex[k])) {
	    // In both cases below, the term to differentiate is
	    // alpha_j(s) alpha_k(s) <Psi_j^2>_xi Psi_j(s) Psi_k(s) and the
	    // difference occurs based on whether a particular s_i dependence
	    // appears in alpha (for inserted) or Psi (for augmented).
	    if (randomVarsKey[deriv_index])
	      // -------------------------------------------
	      // derivative w.r.t. design variable insertion
	      // -------------------------------------------
	      varianceGradient[i] += norm_sq_j *
		(expansionCoeffs[j] * expansionCoeffGrads[k][cntr] +
		 expansionCoeffs[k] * expansionCoeffGrads[j][cntr]) *
		multivariate_polynomial(x, multiIndex[j], nonRandomIndices) *
		multivariate_polynomial(x, multiIndex[k], nonRandomIndices);
	    else
	      // ----------------------------------------------
	      // derivative w.r.t. design variable augmentation
	      // ----------------------------------------------
	      varianceGradient[i] +=
		expansionCoeffs[j] * expansionCoeffs[k] * norm_sq_j *
		// Psi_j * dPsi_k_ds_i + dPsi_j_ds_i * Psi_k
		(multivariate_polynomial(x, multiIndex[j], nonRandomIndices) *
		 multivariate_polynomial_gradient(x, deriv_index, multiIndex[k],
						  nonRandomIndices) +
		 multivariate_polynomial_gradient(x, deriv_index, multiIndex[j],
						  nonRandomIndices) *
		 multivariate_polynomial(x, multiIndex[k], nonRandomIndices));
	  }
	}
      }
    }
    if (randomVarsKey[deriv_index]) // derivative w.r.t. design var insertion
      ++cntr;
  }
  if (all_mode) { computedVariance |=  2; xPrevVarGrad = x; }
  else            computedVariance &= ~2;//deactivate 2-bit: protect mixed usage
  return varianceGradient;
}


/** This test works in combination with DEBUG settings in
    (Legendre,Laguerre,Jacobi,GenLaguerre)OrthogPolynomial::type1_gradient(). */
void OrthogPolyApproximation::gradient_check()
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


void OrthogPolyApproximation::compute_component_sobol()
{
  // sobolIndices are indexed via a bit array, one bit per variable.
  // A bit is turned on for an expansion term if there is a variable
  // dependence (i.e., its multi-index component is non-zero).  Since
  // the Sobol' indices involve a consolidation of variance contributions
  // from the expansion terms, we define a bit array from the multIndex
  // and then use a lookup within sobolIndexMap to assign the expansion
  // term contribution to the correct Sobol' index.

  // compute and sum the variance contributions for each expansion term.  For
  // all_vars mode, this approach picks up the total expansion variance, which
  // is the desired reference pt for type-agnostic global sensitivity analysis.
  size_t i, j, k, sobol_len = sobolIndices.length();
  Real sum_p_var = 0.; RealVector p_var(numExpansionTerms-1, false);
  for (i=1, k=0; i<numExpansionTerms; ++i, ++k) {
    p_var[k] = expansionCoeffs(i) * expansionCoeffs(i)
             * norm_squared(multiIndex[i]);
    sum_p_var += p_var[k];
  }

  // iterate through multiIndex and store sensitivities.  Note: sobolIndices[0]
  // (corresponding to constant exp term with no variable dependence) is unused.
  sobolIndices = 0.; // initialize
  if (sum_p_var > SMALL_NUMBER) { // don't attribute variance if zero/negligible
    BitArray set(numVars);
    for (i=1, k=0; i<numExpansionTerms; ++i, ++k) {

      // determine the bit set corresponding to this expansion term
      for (j=0; j<numVars; ++j)
	if (multiIndex[i][j]) set.set(j);   //   activate bit j
	else                  set.reset(j); // deactivate bit j

      // lookup the bit set within sobolIndexMap --> increment the correct
      // Sobol' index with the variance contribution from this expansion term.
      BAULMIter it = sobolIndexMap.find(set);
      if (it != sobolIndexMap.end()) // may not be found if vbdOrderLimit
	sobolIndices[it->second] += p_var[k];// divide by sum_p_var outside loop
    }
    for (i=0; i<sobol_len; ++i)
      sobolIndices[i] /= sum_p_var;
  }
#ifdef DEBUG
  PCout << "In OrthogPolyApproximation::compute_component_sobol(), "
	<< "sobolIndices =\n"; write_data(PCout, sobolIndices);
#endif // DEBUG
}


void OrthogPolyApproximation::compute_total_sobol() 
{
  totalSobolIndices = 0.;

  if (expConfigOptions.vbdOrderLimit) {
    // all component indices may not be available, so compute total indices
    // from scratch by computing and summing variance contributions for each
    // expansion term
    size_t i, j, k;
    Real sum_p_var = 0., ratio_i;
    RealVector p_var(numExpansionTerms-1, false);
    for (i=1, k=0; i<numExpansionTerms; ++i, ++k) {
      p_var[k] = expansionCoeffs(i) * expansionCoeffs(i)
               * norm_squared(multiIndex[i]);
      sum_p_var += p_var[k];
    }
    // if negligible variance (e.g., a deterministic test fn), then attribution
    // of this variance is suspect.  Defaulting totalSobolIndices to zero is a
    // good choice since it drops out from anisotropic refinement based on the
    // response-average of these indices.
    if (sum_p_var > SMALL_NUMBER) { // avoid division by zero
      Real p_var_k;
      for (i=1, k=0; i<numExpansionTerms; ++i, ++k) {
	p_var_k = p_var[k];
	const UShortArray& mi_i = multiIndex[i];
	// for any constituent variable j in exansion term i, the expansion
	// term contributes to the total sensitivity of variable j
	for (j=0; j<numVars; ++j)
	  if (mi_i[j])
	    totalSobolIndices[j] += p_var_k; // divide by sum_p_var outside loop
      }
      for (j=0; j<numVars; ++j)
	totalSobolIndices[j] /= sum_p_var;
    }
  }
  else {
    // all component effects are present, so simply add them up:
    // totalSobolIndices parses the bit sets of each of the sobolIndices
    // and adds them to each matching variable bin
    for (BAULMIter it=sobolIndexMap.begin(); it!=sobolIndexMap.end(); ++it)
      for (size_t k=0; k<numVars; ++k) 
        if (it->first[k]) // var k is present in this Sobol' index
          totalSobolIndices[k] += sobolIndices[it->second];
  }

#ifdef DEBUG
  PCout << "In OrthogPolyApproximation::compute_total_sobol(), "
	<< "totalSobolIndices =\n"; write_data(PCout, totalSobolIndices);
#endif // DEBUG
}


const RealVector& OrthogPolyApproximation::dimension_decay_rates()
{
  if (decayRates.empty())
    decayRates.sizeUninitialized(numVars);

  // define max_orders for each var for sizing LLS matrices/vectors
  size_t i, j;
  UShortArray max_orders(numVars, 0);
  for (i=0; i<numExpansionTerms; ++i)
    for (j=0; j<numVars; ++j)
      if (multiIndex[i][j] > max_orders[j])
	max_orders[j] = multiIndex[i][j];
  // size A_vectors and b_vectors
  RealVectorArray A_vectors(numVars), b_vectors(numVars);
  for (i=0; i<numVars; ++i) {
    A_vectors[i].sizeUninitialized(max_orders[i]);
    b_vectors[i].sizeUninitialized(max_orders[i]);
  }

  // populate A_vectors and b_vectors
  unsigned short order, non_zero, var_index, order_index;
  bool univariate;
  for (i=1; i<numExpansionTerms; ++i) {
    univariate = true; non_zero = 0;
    for (j=0; j<numVars; ++j) {
      if (multiIndex[i][j]) {
	++non_zero;
	if (non_zero > 1) { univariate = false; break; }
	else { order = multiIndex[i][j]; var_index = j; order_index = order-1; }
      }
    }
    if (univariate) {
      // find a for y = ax + b with x = term order, y = log(coeff), and
      // b = known intercept for order x = 0
      Real norm   = std::sqrt(polynomialBasis[var_index].norm_squared(order)),
	abs_coeff = std::abs(expansionCoeffs[i]);
#ifdef DECAY_DEBUG
      PCout << "Univariate contribution: order = " << order << " coeff = "
	    << abs_coeff << " norm = " << norm << '\n';
#endif // DECAY_DEBUG
      A_vectors[var_index][order_index] = (Real)order;
      b_vectors[var_index][order_index] = (abs_coeff > 1.e-25) ?
	std::log10(abs_coeff * norm) : std::log10(norm) - 25.;
    }
  }
#ifdef DECAY_DEBUG
  PCout << "raw b_vectors:\n";
  for (i=0; i<numVars; ++i)
    { PCout << "Variable " << i+1 << '\n'; write_data(PCout, b_vectors[i]); }
#endif // DECAY_DEBUG

  // first coefficient is used in each of the LLS solves
  Real log_coeff0 = std::log10(std::abs(expansionCoeffs[0])), tol = -10.;
  short last_index_above = -1, new_size;

  for (i=0; i<numVars; ++i) {
    RealVector& A_i = A_vectors[i]; RealVector& b_i = b_vectors[i];

    // Handle case of flatline at numerical precision by ignoring subsequent
    // values below a tolerance (allow first, prune second)
    // > high decay rate will de-emphasize refinement, but consider zeroing
    //   out refinement for a direction that is converged below tol (?)
    // > for now, truncate max_orders and scale back {A,b}_vectors
    order = max_orders[i];
    for (j=0; j<order; ++j)
      if (b_i[j] > tol)
	last_index_above = j;
    new_size = last_index_above+2; // include one value below tolerance
    if (new_size < order) {
      max_orders[i] = order = new_size;
      A_i.resize(order); // needed for psuedo-inv but not LAPACK
      b_i.resize(order); // needed for psuedo-inv but not LAPACK
    }

    // subtract intercept b for y = Ax+b  ==>  Ax = y-b
    for (j=0; j<order; ++j)
      b_i[j] -= log_coeff0;

    // employ simple 1-D pseudo inverse for LLS:
    //   A^T A x = A^T(y-b)  ==>  x = A^T(y-b) / A^T A
    // negate negative slope in log space such that large>0 is fast
    // convergence, small>0 is slow convergence, and <0 is divergence
    decayRates[i] = -A_i.dot(b_i) / A_i.dot(A_i);
  }

#ifdef DECAY_DEBUG
  PCout << "Intercept log(abs(coeff0)) = " << log_coeff0 << '\n';
  PCout << "b_vectors after truncation & intercept subtraction:\n";
  for (i=0; i<numVars; ++i)
    { PCout << "Variable " << i+1 << '\n'; write_data(PCout, b_vectors[i]); }
  PCout << "Individual approximation decay:\n"; write_data(PCout, decayRates);
#endif // DECAY_DEBUG

  return decayRates;
}


void OrthogPolyApproximation::
print_coefficients(std::ostream& s, bool normalized)
{
  size_t i, j;
  char tag[10];

  // terms and term identifiers
  for (i=0; i<numExpansionTerms; ++i) {
    s << "\n  " << std::setw(WRITE_PRECISION+7);
    if (normalized) // basis is divided by norm, so coeff is multiplied by norm
      s << expansionCoeffs[i] * std::sqrt(norm_squared(multiIndex[i]));
    else
      s << expansionCoeffs[i];
    for (j=0; j<numVars; ++j) {
      get_tag(tag, i, j);
      s << std::setw(5) << tag;
    }
  }
  s << '\n';
}


void OrthogPolyApproximation::
coefficient_labels(std::vector<std::string>& coeff_labels) const
{
  size_t i, j;
  char tag[10];

  coeff_labels.reserve(numExpansionTerms);

  // terms and term identifiers
  for (i=0; i<numExpansionTerms; ++i) {
    std::string tags;
    for (j=0; j<numVars; ++j) {
      if (j!=0)
	tags += ' ';
      get_tag(tag, i, j);
      tags += tag;
    }
    coeff_labels.push_back(tags);
  }
}


void OrthogPolyApproximation::get_tag(char* tag, size_t i, size_t j) const
{
  switch (basisTypes[j]) {
  case HERMITE_ORTHOG:
    std::sprintf(tag,  "He%i", multiIndex[i][j]);
    break;
  case LEGENDRE_ORTHOG:
    std::sprintf(tag,   "P%i", multiIndex[i][j]);
    break;
  case LAGUERRE_ORTHOG:
    std::sprintf(tag,   "L%i", multiIndex[i][j]);
    break;
  case JACOBI_ORTHOG:
    std::sprintf(tag, "Pab%i", multiIndex[i][j]);
    break;
  case GEN_LAGUERRE_ORTHOG:
    std::sprintf(tag,  "La%i", multiIndex[i][j]);
    break;
  case CHEBYSHEV_ORTHOG:
    std::sprintf(tag,   "T%i", multiIndex[i][j]);
    break;
  case NUM_GEN_ORTHOG:
    std::sprintf(tag, "Num%i", multiIndex[i][j]);
    break;
  default:
    PCerr << "Error: bad polynomial type = " << basisTypes[j]
	  << " in OrthogPolyApproximation::get_tag()." << std::endl;
    abort_handler(-1);
    break; 
  }
}

} // namespace Pecos
