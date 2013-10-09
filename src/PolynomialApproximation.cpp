/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        PolynomialApproximation
//- Description:  Implementation code for PolynomialApproximation class
//-               
//- Owner:        Mike Eldred

#include "PolynomialApproximation.hpp"
#include "BasisPolynomial.hpp"
#include "SparseGridDriver.hpp"
#include "DistributionParams.hpp"
#include "NumericGenOrthogPolynomial.hpp"

//#define DEBUG


namespace Pecos {

void PolynomialApproximation::
initialize_collocation_rules(const ShortArray& u_types,
			     const BasisConfigOptions& bc_options,
			     ShortArray& colloc_rules)
{
  size_t num_vars = u_types.size();
  colloc_rules.resize(num_vars);

  // set colloc_rules based on u_types: open Gauss rules are used for all
  // global cases; nested closed rules are used for piecewise approximations.
  for (size_t i=0; i<num_vars; ++i) {
    switch (u_types[i]) {
    case STD_NORMAL:
      colloc_rules[i] = (bc_options.nestedRules) ? GENZ_KEISTER : GAUSS_HERMITE;
      break;
    case STD_UNIFORM:
      if (bc_options.piecewiseBasis) // closed nested rules required
	colloc_rules[i] = (bc_options.equidistantRules) ? NEWTON_COTES :
	  CLENSHAW_CURTIS;
      else
	colloc_rules[i] = (bc_options.nestedRules) ? GAUSS_PATTERSON :
	  GAUSS_LEGENDRE;
      // For tensor-product quadrature without refinement, Gauss-Legendre
      // is preferred due to greater polynomial exactness since nesting is
      // not a concern.  For sparse grids and quadrature with refinement,
      // Gauss-Patterson or Clenshaw-Curtis can be better options.
      break;
    case STD_EXPONENTIAL: colloc_rules[i] = GAUSS_LAGUERRE;     break;
    case STD_BETA:        colloc_rules[i] = GAUSS_JACOBI;       break;
    case STD_GAMMA:       colloc_rules[i] = GEN_GAUSS_LAGUERRE; break;
    default:              colloc_rules[i] = GOLUB_WELSCH;       break;
    }
  }
}


void PolynomialApproximation::
initialize_polynomial_basis(const ShortArray& basis_types,
			    const ShortArray& colloc_rules,
			    std::vector<BasisPolynomial>& poly_basis)
{
  size_t i, num_vars = basis_types.size(), num_rules = colloc_rules.size();

  // instantiate poly_basis using basis_types and colloc_rules
  if (poly_basis.size() != num_vars) {
    poly_basis.resize(num_vars);
    if (num_rules == num_vars)
      for (i=0; i<num_vars; ++i)
	poly_basis[i] = BasisPolynomial(basis_types[i], colloc_rules[i]);
    else if (num_rules == 0)   // use default rules for each basis type
      for (i=0; i<num_vars; ++i)
	poly_basis[i] = BasisPolynomial(basis_types[i]);
    else if (num_rules == 1) { // cubature utilizes a single rule
      short colloc_rule = colloc_rules[0];
      for (i=0; i<num_vars; ++i)
	poly_basis[i] = BasisPolynomial(basis_types[i], colloc_rule);
    }

    /*
    // Could reuse objects as in InterpPolyApproximation, but this would require
    // caching dist params and moving code from distribution_parameters() to
    // here, which isn't yet well motivated.
    size_t i, j;
    for (i=0; i<numVars; ++i) {
      // reuse prev instance via shared rep or instantiate new unique instance
      short basis_type_i = basis_types[i];
      bool found = false;
      for (j=0; j<i; ++j)
	if ( basis_type_i == basis_types[j] && ( basis_type_i <= LAGUERRE ||
	     ( basis_type_i == JACOBI && jacobiAlphas[] == jacobiAlphas[] &&
	       jacobiBetas[] == jacobiBetas[] ) ||
	     ( basis_type_i == GENERALIZED_LAGUERRE &&
	       genLagAlphas[] == genLagAlphas[] ) ) )
	  { found = true; break; }
      polynomialBasis[i]
	= (found) ? polynomialBasis[j] : BasisPolynomial(basis_type_i);
    }
    */
  }
}


void PolynomialApproximation::
update_basis_distribution_parameters(const ShortArray& u_types,
				     const AleatoryDistParams& adp,
				     std::vector<BasisPolynomial>& poly_basis)
{
  size_t i, num_vars = u_types.size(), nuv_cntr = 0, lnuv_cntr = 0,
    luuv_cntr = 0, tuv_cntr = 0, beuv_cntr = 0, gauv_cntr = 0, guuv_cntr = 0,
    fuv_cntr = 0, wuv_cntr = 0, hbuv_cntr = 0;

  // update poly_basis using distribution data from dp
  for (i=0; i<num_vars; ++i)
    switch (u_types[i]) {
    case STD_NORMAL:
      ++nuv_cntr; break;
    case STD_UNIFORM: case STD_EXPONENTIAL:
      break;
    case STD_BETA:
      poly_basis[i].alpha_stat(adp.beta_alpha(beuv_cntr));
      poly_basis[i].beta_stat(adp.beta_beta(beuv_cntr));
      ++beuv_cntr; break;
    case STD_GAMMA:
      poly_basis[i].alpha_stat(adp.gamma_alpha(gauv_cntr));
      ++gauv_cntr; break;
    case BOUNDED_NORMAL:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	bounded_normal_distribution(adp.normal_mean(nuv_cntr),
				    adp.normal_std_deviation(nuv_cntr),
				    adp.normal_lower_bound(nuv_cntr),
				    adp.normal_upper_bound(nuv_cntr));
      ++nuv_cntr; break;
    case LOGNORMAL: {
      Real mean, stdev;
      moments_from_lognormal_spec(adp.lognormal_means(),
				  adp.lognormal_std_deviations(),
				  adp.lognormal_lambdas(),
				  adp.lognormal_zetas(),
				  adp.lognormal_error_factors(), lnuv_cntr,
				  mean, stdev);
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	lognormal_distribution(mean, stdev);
      ++lnuv_cntr; break;
    }
    case BOUNDED_LOGNORMAL: {
      Real mean, stdev;
      moments_from_lognormal_spec(adp.lognormal_means(),
				  adp.lognormal_std_deviations(),
				  adp.lognormal_lambdas(),
				  adp.lognormal_zetas(),
				  adp.lognormal_error_factors(), lnuv_cntr,
				  mean, stdev);
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	bounded_lognormal_distribution(mean, stdev,
				       adp.lognormal_lower_bound(lnuv_cntr),
				       adp.lognormal_upper_bound(lnuv_cntr));
      ++lnuv_cntr; break;
    }
    case LOGUNIFORM:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	loguniform_distribution(adp.loguniform_lower_bound(luuv_cntr),
				adp.loguniform_upper_bound(luuv_cntr));
      ++luuv_cntr; break;
    case TRIANGULAR:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	triangular_distribution(adp.triangular_mode(tuv_cntr),
				adp.triangular_lower_bound(tuv_cntr),
				adp.triangular_upper_bound(tuv_cntr));
      ++tuv_cntr; break;
    case GUMBEL:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	gumbel_distribution(adp.gumbel_alpha(guuv_cntr),
			    adp.gumbel_beta(guuv_cntr));
      ++guuv_cntr; break;
    case FRECHET:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	frechet_distribution(adp.frechet_alpha(fuv_cntr),
			     adp.frechet_beta(fuv_cntr));
      ++fuv_cntr; break;
    case WEIBULL:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	weibull_distribution(adp.weibull_alpha(wuv_cntr),
			     adp.weibull_beta(wuv_cntr));
      ++wuv_cntr; break;
    case HISTOGRAM_BIN:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	histogram_bin_distribution(adp.histogram_bin_pairs(hbuv_cntr));
      ++hbuv_cntr; break;
    default:
      PCerr << "Error: unsupported u-space type in PolynomialApproximation::"
	    << "distribution_parameters()" << std::endl;
      abort_handler(-1);
      break;
    }
}


/*
void PolynomialApproximation::allocate_component_sobol()
{
  // default implementation is reasonable for tensor expansions, but is
  // wasteful (and should be overridden) for total-order and sparse 
  // sum-of-tensor (which emulate total-order) expansions.
  if (expConfigOptions.vbdControl && expConfigOptions.expansionCoeffFlag &&
      sobolIndices.empty()) {
    switch (expConfigOptions.vbdControl) {
    case UNIVARIATE_VBD: // main effects only
      allocate_main_sobol();                    break;
    case ALL_VBD: // main + all n-way interactions
      allocate_main_interaction_sobol(numVars); break;
    }
  }
}
*/


void PolynomialApproximation::increment_component_sobol()
{
  PCerr << "Error: derived class does not redefine increment_component_sobol()."
	<< std::endl;
  abort_handler(-1);
}


void PolynomialApproximation::allocate_main_sobol()
{
  // define binary sets corresponding to main effects
  BitArray set(numVars, 0);
  // prepend the 0-way interaction for indexing consistency
  sobolIndexMap[set] = 0;
  for (size_t v=0; v<numVars; ++v) // activate bit for variable v
    { set.set(v); sobolIndexMap[set] = v+1; set.reset(v); }
  sobolIndices.sizeUninitialized(numVars+1);
}


/*
void PolynomialApproximation::
allocate_main_interaction_sobol(unsigned short max_order)
{
  // include m-way interactions for m <= max_order.  Note: sobol index 0 is
  // unused (corresponds to constant exp term with no variable dependence).
  unsigned long sobol_len = 0; size_t v;
  BitArray set(numVars, 0);
  if (max_order >= 1) {
    sobol_len += numVars;
    for (v=0; v<numVars; ++v)
      { set.set(v); sobolIndexMap[set] = v; set.reset(v); }
  }
  unsigned long cntr = sobol_len;
  for (unsigned short ord=2; ord<=max_order; ++ord) {
    // compute number of terms from n_choose_k() for each exp order
    sobol_len += BasisPolynomial::n_choose_k(numVars, ord);
    // use increment_terms() and exclude term equalities to get interactions
    UShortArray terms(ord, 1); // # of terms = level
    bool order_complete = false;
    while (!order_complete) {
      size_t last_index = ord - 1, prev_index = ord - 2;
      for (terms[last_index]=1; terms[last_index]<terms[prev_index];
	   ++terms[last_index]) {
	// convert orders (within terms) to variable indices (within set)
	set.reset();
	for (size_t i=0; i<ord; ++i)
	  set.set(terms[i]-1);
	sobolIndexMap[set] = cntr; ++cntr;
      }
      increment_terms(terms, last_index, prev_index, numVars,
		      order_complete, true);
    }
  }

  sobolIndices.sizeUninitialized(sobol_len);
}
*/


void PolynomialApproximation::allocate_total_sobol()
{
  // number of total indices independent of number of component indices
  if (expConfigOptions.vbdControl && expConfigOptions.expansionCoeffFlag &&
      totalSobolIndices.empty())
    totalSobolIndices.sizeUninitialized(numVars);
}


void PolynomialApproximation::reset_sobol_index_map_values()
{
  for (BAULMIter it=sobolIndexMap.begin(); it!=sobolIndexMap.end(); ++it)
    it->second = it->first.count(); // order of interaction
}


void PolynomialApproximation::
multi_index_to_sobol_index_map(const UShort2DArray& mi)
{
  BitArray set(numVars); size_t i, j, num_mi = mi.size();
  for (i=0; i<num_mi; ++i) {
    // determine the bit set corresponding to this expansion term
    for (j=0; j<numVars; ++j)
      if (mi[i][j]) set.set(j);   //   activate bit j
      else          set.reset(j); // deactivate bit j

    // If set does not already exist, insert it.
    // For value in key-value pair, initially use the interaction order;
    // will be updated below in sobol_index_map_to_sobol_indices().  The
    // 0-way interaction is included to support child lookup requirements
    // in InterpPolyApproximation::compute_partial_variance().
    if (sobolIndexMap.find(set) == sobolIndexMap.end())
      sobolIndexMap[set] = set.count(); // order of interaction
  }
}


void PolynomialApproximation::sobol_index_map_to_sobol_indices()
{
  // Compute total counts for interactions within each interaction group.
  // The 0-way interaction is included to support child lookup requirements
  // in InterpPolyApproximation::compute_partial_variance().
  ULongArray counters(numVars+1, 0);
  for (BAULMIter it=sobolIndexMap.begin(); it!=sobolIndexMap.end(); ++it)
    ++counters[it->second];

  // aggregate the totals into starting indices for each interaction group
  ULongArray indices(numVars+1); indices[0] = 0;
  for (size_t i=1; i<=numVars; ++i)
    indices[i] = indices[i-1] + counters[i-1];
  unsigned long sobol_len = indices[numVars] + counters[numVars];

  // Now that the sobolIndexMap sets are fully defined and sorted, define the
  // mapping of these sets to sobolIndices based on counters.  We choose to
  // group the sobolIndices by the number of interactions to simplify output.
  unsigned long interaction_order;
  for (BAULMIter it=sobolIndexMap.begin(); it!=sobolIndexMap.end(); ++it) {
    interaction_order = it->second;
    it->second = indices[interaction_order]++;
  }

  sobolIndices.sizeUninitialized(sobol_len);
}


/** Return the number of terms in a tensor-product expansion.  For
    isotropic and anisotropic expansion orders, calculation of the
    number of expansion terms is straightforward: Prod(p_i + 1). */
size_t PolynomialApproximation::
tensor_product_terms(const UShortArray& order, bool include_upper_bound)
{
  size_t i, n = order.size(), num_terms = 1;
  if (include_upper_bound)
    for (i=0; i<n; ++i)
      num_terms *= order[i] + 1; // multi-index from expansion order p (default)
  else
    for (i=0; i<n; ++i)
      num_terms *= order[i];     // collocation index from quadrature order m
  return num_terms;
}


void PolynomialApproximation::
tensor_product_multi_index(const UShortArray& order, UShort2DArray& multi_index,
			   bool include_upper_bound)
{
  // This function may be used for defining collocation indices for quadrature
  // points (indices = 0:m-1) or polynomial order multi-indices for tensor
  // expansions (indices = 0:p).  The inclusion or exclusion of the upper bound
  // (m or p) is managed with include_upper_bound (false for m, true for p).

  // rather than inserting into multi_index, go ahead and estimate its length
  // (since its inexpensive) and use term-by-term assignment.
  size_t i, n = order.size(),
    mi_len = tensor_product_terms(order, include_upper_bound);
  if (mi_len != multi_index.size())
    multi_index.resize(mi_len);
  UShortArray indices(n, 0); multi_index[0] = indices;
  for (i=1; i<mi_len; ++i) {
    // increment the n-dimensional index set
    increment_indices(indices, order, include_upper_bound);
    multi_index[i] = indices;
  }
}


void PolynomialApproximation::
hierarchical_tensor_product_multi_index(const UShort2DArray& delta_quad,
					UShort2DArray& multi_index)
{
  // delta_quad are non-sequential indices of hierarchical collocation point
  // increments
  size_t i, j, n = delta_quad.size(), mi_len = 1;
  UShortArray index_bound(n);
  for (i=0; i<n; ++i)
    mi_len *= index_bound[i] = delta_quad[i].size();
  if (mi_len != multi_index.size())
    multi_index.resize(mi_len);
  UShortArray indices(n, 0);
  for (i=0; i<mi_len; ++i) {
    multi_index[i].resize(n);
    for (j=0; j<n; ++j)
      multi_index[i][j] = delta_quad[j][indices[j]];
    if (i < mi_len-1)
      increment_indices(indices, index_bound, false); // 0 <= index < bound
  }
}


/** Return the number of terms in a total-order expansion.  For anisotropic
    expansion order, no simple expression is currently available and the
    number of expansion terms is computed using the multiIndex recursion. */
size_t PolynomialApproximation::
total_order_terms(const UShortArray& upper_bound, short lower_bound_offset)
{
  size_t i, n = upper_bound.size();
  if (!n) {
    PCerr << "Error: empty upper_bound in PolynomialApproximation::"
	  << "total_order_terms()." << std::endl;
    abort_handler(-1);
  }

  bool isotropic = true;
  unsigned short order = upper_bound[0];
  for (i=1; i<n; ++i)
    if (upper_bound[i] != order)
      { isotropic = false; break; }

  size_t num_terms;
  if (isotropic) {
    num_terms = (size_t)std::floor(
      BasisPolynomial::n_choose_k(order+n, order)+.5); // round non-integral
    if (lower_bound_offset >= 0) { // default is -1
      int omit_order = order - lower_bound_offset - 1;
      if (omit_order >= 0)
	num_terms -= (size_t)std::floor(
	  BasisPolynomial::n_choose_k(omit_order+n, omit_order)+.5); // round
    }
  }
  else { // anisotropic: use multiIndex recursion to compute
    bool mi_lower_bound = (lower_bound_offset >= 0); // default is -1
    if (mi_lower_bound) {
      // Smolyak combinatorial form is invalid for anisotropic levels
      PCerr << "Error: anisotropic orders not currently supported with "
	    << "multi-index lower bound\n       in PolynomialApproximation::"
	    << "total_order_terms()." << std::endl;
      abort_handler(-1);
    }
    unsigned short max_order = order, order_nd;
    for (i=1; i<n; ++i)
      max_order = std::max(max_order, upper_bound[i]);
    num_terms = 1; // order 0
    if (max_order >= 1)   // order 1
      for (i=0; i<n; ++i)
	if (upper_bound[i] >= 1) // only upper bound may be anisotropic
	  ++num_terms;
    for (order_nd=2; order_nd<=max_order; ++order_nd) { // order 2 through max
      UShortArray terms(order_nd, 1); // # of terms = current order
      bool order_complete = false;
      while (!order_complete) {
	size_t last_index = order_nd - 1, prev_index = order_nd - 2;
	for (terms[last_index]=1; terms[last_index]<=terms[prev_index]; 
	     ++terms[last_index]) {
	  bool include = true;
	  for (i=0; i<n; ++i) {
	    if (std::count(terms.begin(), terms.end(), i+1) > upper_bound[i]) {
	      // only upper bound may be anisotropic
	      include = false;
	      break;
	    }
	  }
	  if (include)
	    ++num_terms;
	}
	// increment term hierarchy
	increment_terms(terms, last_index, prev_index, n, order_complete);
      }
    }
  }
  return num_terms;
}


/** Overloaded version for defining the multi-indices for a single
    scalar level.  Anisotropy is not supported, so this version is not
    usable as a kernel within other overloaded versions. */
void PolynomialApproximation::
total_order_multi_index(unsigned short level, size_t num_vars, 
			UShort2DArray& multi_index)
{
  UShortArray mi(num_vars, 0);
  multi_index.clear();
  // special logic required for level < 2 due to prev_index defn below
  switch (level) {
  case 0:
    multi_index.push_back(mi); break;
  case 1:
    for (size_t i=0; i<num_vars; ++i)
      { mi[i] = 1; multi_index.push_back(mi); mi[i] = 0; }
    break;
  default: {
    UShortArray terms(level, 1); // # of terms = level
    bool order_complete = false;
    while (!order_complete) {
      // this is the inner-most loop w/i the nested looping managed by terms
      size_t last_index = level - 1, prev_index = level - 2;
      for (terms[last_index]=1;
	   terms[last_index]<=terms[prev_index]; ++terms[last_index]) {
	for (size_t i=0; i<num_vars; ++i)
	  mi[i] = std::count(terms.begin(), terms.end(), i+1);
	multi_index.push_back(mi);
      }
      increment_terms(terms, last_index, prev_index, num_vars, order_complete);
    }
    break;
  }
  }
}


void PolynomialApproximation::
total_order_multi_index(const UShortArray& upper_bound,
			UShort2DArray& multi_index, short lower_bound_offset, 
			size_t max_terms)
{
  // populate multi_index: implementation follows ordering of Eq. 4.1 in
  // [Xiu and Karniadakis, 2002].
  // To handle anisotropy, we currently perform a complete total-order
  // recursion based on max_order and reject multi-indices with components
  // greater than those in upper_bound.
  size_t i, cntr = 0, n = upper_bound.size();
  if (!n) {
    PCerr << "Error: empty upper_bound in PolynomialApproximation::"
	  << "total_order_multi_index()." << std::endl;
    abort_handler(-1);
  }

  // Since total_order_terms() is expensive, use insertions into multi_index
  // rather than sizing with entry-by-entry assignment
  bool isotropic = true;
  unsigned short order = upper_bound[0];
  for (i=1; i<n; ++i)
    if (upper_bound[i] != order)
      { isotropic = false; break; }

  bool mi_lower_bound = (lower_bound_offset >= 0); // default is -1
  unsigned short max_order = order, min_order = 0, order_nd;
  if (isotropic) {
    if (mi_lower_bound)
      min_order = (lower_bound_offset >= order)
	? 0 : order - lower_bound_offset; // e.g., Smolyak l.b. = w-N+1
  }
  else {
    if (mi_lower_bound) {
      // Smolyak combinatorial form is invalid for anisotropic levels
      PCerr << "Error: anisotropic orders not currently supported with "
	    << "multi-index lower bound\n       in PolynomialApproximation::"
	    << "total_order_multi_index()." << std::endl;
      abort_handler(-1);
    }
    for (i=1; i<n; ++i)
      max_order = std::max(max_order, upper_bound[i]);
  }

  // special logic required for order_nd < 2 due to prev_index defn below
  UShortArray mi(n,0);
  multi_index.clear();
  if (min_order == 0) // && max_order >= 0
    multi_index.push_back(mi); // order 0
  if (min_order <= 1 && max_order >= 1) { // order 1
    for (i=0; i<n && cntr<max_terms; ++i) {
      if (upper_bound[i] >= 1) { // only upper bound may be anisotropic
	mi[i] = 1; // ith entry is nonzero
	multi_index.push_back(mi);
	mi[i] = 0; // reset
      }
    }
  }
  for (order_nd=std::max(min_order,(unsigned short)2);
       order_nd<=max_order; ++order_nd) {
    UShortArray terms(order_nd, 1);//# of terms = current order
    bool order_complete = false;
    while (!order_complete) {
      // this is the inner-most loop w/i the nested looping managed by terms
      size_t last_index = order_nd - 1, prev_index = order_nd - 2;
      for (terms[last_index]=1;
	   terms[last_index]<=terms[prev_index] && cntr<max_terms;
	   ++terms[last_index]) {
	// store the orders of the univariate polynomials to be used for
	// constructing the current multivariate basis function
	bool include = true;
	for (i=0; i<n; ++i) {
	  mi[i] = std::count(terms.begin(), terms.end(), i+1);
	  if (mi[i] > upper_bound[i]) { // only upper bound may be anisotropic
	    include = false;
	    break;
	  }
	}
	if (include)
	  multi_index.push_back(mi);
#ifdef DEBUG
	PCout << "terms:\n" << terms << std::endl;
#endif // DEBUG
      }
      if (cntr == max_terms)
	order_complete = true;
      else // increment term hierarchy
	increment_terms(terms, last_index, prev_index, n, order_complete);
    }
  }

#ifdef DEBUG
  size_t mi_len = multi_index.size();
  for (i=0; i<mi_len; ++i)
    PCout << "multiIndex[" << i << "]:\n" << multi_index[i] << std::endl;
#endif // DEBUG
}


void PolynomialApproximation::
compute_numerical_moments(const RealVector& coeffs, const RealVector& t1_wts,
			  RealVector& moments)
{
  // computes and stores the following moments:
  // > mean     (1st raw moment)
  // > variance (2nd central moment)
  // > skewness (3rd standardized moment)
  // > kurtosis (4th standardized moment with offset to eliminate "excess")

  // current support for this implementation: can't be open-ended since we
  // employ a specific combination of raw, central, and standardized moments
  size_t num_moments = moments.length();
  if (num_moments < 1 || num_moments > 4) {
    PCerr << "Error: unsupported number of moments requested in Polynomial"
	  << "Approximation::compute_numerical_moments()" << std::endl;
    abort_handler(-1);
  }
  size_t i, j, num_pts = coeffs.length();
  if (t1_wts.length() != num_pts) {
    PCerr << "Error: mismatch in array lengths between integration driver "
	  << "weights ("  << t1_wts.length() << ") and coefficients ("
	  << num_pts << ") in PolynomialApproximation::compute_numerical_"
	  << "moments()." << std::endl;
    abort_handler(-1);
  }

  // estimate 1st raw moment (mean)
  moments = 0.;
  Real& mean = moments[0];
  for (i=0; i<num_pts; ++i)
    mean += t1_wts[i] * coeffs[i];

  // estimate central moments 2 through num_moments
  if (num_moments > 1) {
    Real centered_fn, pow_fn;
    for (i=0; i<num_pts; ++i) {
      pow_fn = centered_fn = coeffs[i] - mean;
      for (j=1; j<num_moments; ++j) {
	pow_fn     *= centered_fn;
	moments[j] += t1_wts[i] * pow_fn;
      }
    }
  }

  // standardize third and higher central moments, if present
  //standardize_moments(moments);
}


void PolynomialApproximation::
compute_numerical_moments(const RealVector& t1_coeffs,
			  const RealMatrix& t2_coeffs, const RealVector& t1_wts,
			  const RealMatrix& t2_wts, RealVector& moments)
{
  // computes and stores the following moments:
  // > mean     (1st raw moment)
  // > variance (2nd central moment)
  // > skewness (3rd standardized moment)
  // > kurtosis (4th standardized moment with offset to eliminate "excess")

  // current support for this implementation: can't be open-ended since we
  // employ a specific combination of raw, central, and standardized moments
  size_t num_moments = moments.length();
  if (num_moments < 1 || num_moments > 4) {
    PCerr << "Error: unsupported number of moments requested in Polynomial"
	  << "Approximation::compute_numerical_moments()" << std::endl;
    abort_handler(-1);
  }
  size_t i, j, k, num_pts = t1_coeffs.length();
  if (t1_wts.length() != num_pts || t2_wts.numCols() != num_pts ||
      t2_coeffs.numCols() != num_pts) {
    PCerr << "Error: mismatch in array lengths among integration driver "
	  << "weights ("  << t1_wts.length() << ", " << t2_wts.numCols()
	  << ") and coefficients (" << num_pts << ", " << t2_coeffs.numCols()
	  << ") in PolynomialApproximation::compute_numerical_moments()."
	  << std::endl;
    abort_handler(-1);
  }

  // estimate 1st raw moment (mean)
  moments = 0.;
  Real& mean = moments[0];
  for (i=0; i<num_pts; ++i) {
    mean += t1_wts[i] * t1_coeffs[i];
    const Real* coeff2_i = t2_coeffs[i];
    const Real*  t2_wt_i = t2_wts[i];
    for (j=0; j<numVars; ++j)
      mean += coeff2_i[j] * t2_wt_i[j];
  }

  // estimate central moments 2 through num_moments
  if (num_moments > 1) {
    Real centered_fn, pow_fn;
    for (i=0; i<num_pts; ++i) {
      pow_fn = centered_fn = t1_coeffs[i] - mean;
      const Real* coeff2_i = t2_coeffs[i];
      const Real*  t2_wt_i = t2_wts[i];
      for (j=1; j<num_moments; ++j) {
	Real& moment_j = moments[j];
	// type2 interpolation of (R - \mu)^n
	// --> interpolated gradients are n(R - \mu)^{n-1} dR/dx
	for (k=0; k<numVars; ++k)
	  moment_j += (j+1) * pow_fn * coeff2_i[k] * t2_wt_i[k];
	// type 1 interpolation of (R - \mu)^n
	pow_fn   *= centered_fn;
	moment_j += t1_wts[i] * pow_fn;
      }
    }
  }

  // convert central moments to std deviation/skewness/kurtosis
  //standardize_moments(moments);
}


void PolynomialApproximation::
standardize_moments(const RealVector& central_moments, RealVector& std_moments)
{
  size_t num_moments = central_moments.length();
  std_moments.sizeUninitialized(num_moments);
  if (num_moments >= 1) std_moments[0] = central_moments[0]; // mean
  if (num_moments <  2) return;

  const Real& var = central_moments[1];
  Real&   std_dev = std_moments[1];
  if (var > 0.) {
    // standardized moment k is E[((X-mu)/sigma)^k] = E[(X-mu)^k]/sigma^k
    std_dev = std::sqrt(var); // not standardized (2nd standardized moment is 1)
    Real pow_fn = var;
    for (size_t i=2; i<num_moments; ++i)
      { pow_fn *= std_dev; std_moments[i] = central_moments[i] / pow_fn; }
    // offset the fourth standardized moment to eliminate excess kurtosis
    if (num_moments > 3)
      std_moments[3] -= 3.;
  }
  else {
    // don't leave uninitialized, even if undefined
    for (size_t i=1; i<num_moments; ++i)
      std_moments[i] = 0.;
    // special case of zero variance is OK for num_moments == 2, but not higher
    if ( !(num_moments == 2 && var == 0.) ) // std_dev OK for var == 0.
      PCerr << "Warning: moments cannot be standardized due to non-positive "
	    << "variance.\n         Skipping standardization." << std::endl;
  }
}


Real PolynomialApproximation::
delta_covariance(PolynomialApproximation* poly_approx_2)
{
  PCerr << "Error: delta_covariance() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
delta_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  PCerr << "Error: delta_covariance() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_mean()
{
  PCerr << "Error: delta_mean() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_mean(const RealVector& x)
{
  PCerr << "Error: delta_mean(x) not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_std_deviation()
{
  PCerr << "Error: delta_std_deviation() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_std_deviation(const RealVector& x)
{
  PCerr << "Error: delta_std_deviation(x) not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_beta(bool cdf_flag, Real z_bar)
{
  PCerr << "Error: delta_beta() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
delta_beta(const RealVector& x, bool cdf_flag, Real z_bar)
{
  PCerr << "Error: delta_beta(x) not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_z(bool cdf_flag, Real beta_bar)
{
  PCerr << "Error: delta_z() not available for this polynomial approximation "
	<< "type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
delta_z(const RealVector& x, bool cdf_flag, Real beta_bar)
{
  PCerr << "Error: delta_z(x) not available for this polynomial approximation "
	<< "type." << std::endl;
  abort_handler(-1);
  return 0.;
}


const RealVector& PolynomialApproximation::dimension_decay_rates()
{
  PCerr << "Error: dimension_decay_rates() not available for this polynomial "
	<< "approximation type." << std::endl;
  return abort_handler_t<const RealVector&>(-1);
}


void PolynomialApproximation::increment_order()
{
  PCerr << "Error: increment_order() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
}


void PolynomialApproximation::update_order(const UShortArray& order)
{
  PCerr << "Error: update_order() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
}

} // namespace Pecos
