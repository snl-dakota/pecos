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

void PolynomialApproximation::allocate_component_effects()
{
  // sobolIndices[0] is reserved for mean 
 
  // Allocate memory specific to output control
  if (configOptions.vbdControl && configOptions.expansionCoeffFlag &&
      sobolIndices.empty()) {
    int i, index_length;
    switch (configOptions.vbdControl) {
    case UNIVARIATE_VBD: // main effects only
      index_length = (int)numVars+1;
      // map binary to integer main effects form
      sobolIndexMap[0] = 0;
      for (i=1; i<index_length; ++i) 
	sobolIndexMap[int(std::pow(2.,i-1))] = i;
      break;
    case ALL_VBD: // main + interaction effects
      // don't recompute total separately; rather sum from component effects
      index_length = (int)std::pow(2.,(int)numVars);
      for (i=0; i<index_length; ++i) 
	sobolIndexMap[i] = i; // already in binary notation
      break;
    }
    sobolIndices.sizeUninitialized(index_length);
  }
}


void PolynomialApproximation::allocate_total_effects()
{
  // number of total indices independent of number of component indices
  if (configOptions.vbdControl && configOptions.expansionCoeffFlag &&
      totalSobolIndices.empty())
    totalSobolIndices.sizeUninitialized(numVars);
}


/** This static function is used by the integration drivers to define a
    basis for computing points and weights.  Thus, it only employs bases
    that support a "pull" mode (orthogonal and piecewise interpolation
    polynomials) and excludes bases that utilize a "push" mode (Lagrange
    and Hermite global interpolation polynomials). [Note: piecewise 
    interpolation polynomials could support push of points and then pull
    of weights, but for now, it is pull for both.] */
bool PolynomialApproximation::
distribution_types(const ShortArray& u_types, bool piecewise_basis,
		   bool use_derivs, ShortArray& basis_types)
{
  bool extra_dist_params = false;
  size_t i, num_vars = u_types.size();

  // Initialize basis_types from u_types.
  if (basis_types.size() != num_vars) {
    basis_types.resize(num_vars);
    for (i=0; i<num_vars; ++i) {
      switch (u_types[i]) {
      case STD_NORMAL:
	basis_types[i] = HERMITE_ORTHOG;                                  break;
      case STD_UNIFORM:
	if (piecewise_basis)
	  basis_types[i] = (use_derivs) ? PIECEWISE_CUBIC_INTERP:
	                                  PIECEWISE_LINEAR_INTERP;
	else // Legendre or Chebyshev OrthogPolynomial (Lagrange and Hermite
	     // interpolation not supported -- see doxygen comment above).
	  basis_types[i] = LEGENDRE_ORTHOG;
	  // To employ Chebyshev for uniform, have to multiply inner product
	  // integrands by the inverse of the weight function (weight fn =
	  // 1/sqrt(1-x^2); same as beta PDF/Jacobi poly for alpha=beta=-1/2).
	  //basis_types[i] = CHEBYSHEV_ORTHOG;
	break;
      case STD_EXPONENTIAL:
	basis_types[i] = LAGUERRE_ORTHOG;                                 break;
      case STD_BETA:
	basis_types[i] = JACOBI_ORTHOG;
	extra_dist_params = true;                                         break;
      case STD_GAMMA:
	basis_types[i] = GEN_LAGUERRE_ORTHOG;
	extra_dist_params = true;                                         break;
      case BOUNDED_NORMAL: case BOUNDED_LOGNORMAL: case LOGNORMAL:
      case LOGUNIFORM: case TRIANGULAR:
      case GUMBEL: case FRECHET: case WEIBULL: case HISTOGRAM_BIN:
	basis_types[i] = NUM_GEN_ORTHOG;
	extra_dist_params = true;                                         break;
      default:
	PCerr << "Error: unsupported u-space type in PolynomialApproximation::"
	      << "distribution_types()." << std::endl;
	abort_handler(-1);                                                break;
      }
    }
  }
  else
    for (i=0; i<num_vars; ++i)
      if (u_types[i] != STD_NORMAL && u_types[i] != STD_UNIFORM &&
	  u_types[i] != STD_EXPONENTIAL)
	{ extra_dist_params = true; break; }

  return extra_dist_params;
}


void PolynomialApproximation::
distribution_rules(const ShortArray& u_types, bool nested_rules,
		   bool  piecewise_basis,     bool equidistant_rules,
		   ShortArray& colloc_rules)
{
  size_t i, num_vars = u_types.size();
  colloc_rules.resize(num_vars);

  // set colloc_rules based on u_types: open Gauss rules are used for all
  // global cases; nested closed rules are used for piecewise approximations.
  for (size_t i=0; i<num_vars; ++i) {
    switch (u_types[i]) {
    case STD_NORMAL:
      colloc_rules[i] = (nested_rules) ? GENZ_KEISTER : GAUSS_HERMITE; break;
    case STD_UNIFORM:
      if (piecewise_basis) // closed nested rules required
	colloc_rules[i] = (equidistant_rules) ? NEWTON_COTES : CLENSHAW_CURTIS;
      else
	colloc_rules[i] = (nested_rules) ? GAUSS_PATTERSON : GAUSS_LEGENDRE;
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
distribution_basis(const ShortArray& basis_types,
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
distribution_parameters(const ShortArray& u_types, const DistributionParams& dp,
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
      poly_basis[i].alpha_stat(dp.beta_alpha(beuv_cntr));
      poly_basis[i].beta_stat(dp.beta_beta(beuv_cntr));
      ++beuv_cntr; break;
    case STD_GAMMA:
      poly_basis[i].alpha_stat(dp.gamma_alpha(gauv_cntr));
      ++gauv_cntr; break;
    case BOUNDED_NORMAL:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	bounded_normal_distribution(dp.normal_mean(nuv_cntr),
				    dp.normal_std_deviation(nuv_cntr),
				    dp.normal_lower_bound(nuv_cntr),
				    dp.normal_upper_bound(nuv_cntr));
      ++nuv_cntr; break;
    case LOGNORMAL: {
      Real mean, stdev;
      moments_from_lognormal_spec(dp.lognormal_means(),
				  dp.lognormal_std_deviations(),
				  dp.lognormal_lambdas(), dp.lognormal_zetas(),
				  dp.lognormal_error_factors(), lnuv_cntr,
				  mean, stdev);
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	lognormal_distribution(mean, stdev);
      ++lnuv_cntr; break;
    }
    case BOUNDED_LOGNORMAL: {
      Real mean, stdev;
      moments_from_lognormal_spec(dp.lognormal_means(),
				  dp.lognormal_std_deviations(),
				  dp.lognormal_lambdas(), dp.lognormal_zetas(),
				  dp.lognormal_error_factors(), lnuv_cntr,
				  mean, stdev);
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	bounded_lognormal_distribution(mean, stdev,
				       dp.lognormal_lower_bound(lnuv_cntr),
				       dp.lognormal_upper_bound(lnuv_cntr));
      ++lnuv_cntr; break;
    }
    case LOGUNIFORM:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	loguniform_distribution(dp.loguniform_lower_bound(luuv_cntr),
				dp.loguniform_upper_bound(luuv_cntr));
      ++luuv_cntr; break;
    case TRIANGULAR:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	triangular_distribution(dp.triangular_mode(tuv_cntr),
				dp.triangular_lower_bound(tuv_cntr),
				dp.triangular_upper_bound(tuv_cntr));
      ++tuv_cntr; break;
    case GUMBEL:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	gumbel_distribution(dp.gumbel_alpha(guuv_cntr),
			    dp.gumbel_beta(guuv_cntr));
      ++guuv_cntr; break;
    case FRECHET:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	frechet_distribution(dp.frechet_alpha(fuv_cntr),
			     dp.frechet_beta(fuv_cntr));
      ++fuv_cntr; break;
    case WEIBULL:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	weibull_distribution(dp.weibull_alpha(wuv_cntr),
			     dp.weibull_beta(wuv_cntr));
      ++wuv_cntr; break;
    case HISTOGRAM_BIN:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	histogram_bin_distribution(dp.histogram_bin_pairs(hbuv_cntr));
      ++hbuv_cntr; break;
    default:
      PCerr << "Error: unsupported u-space type in PolynomialApproximation::"
	    << "distribution_parameters()" << std::endl;
      abort_handler(-1);
      break;
    }
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
    num_terms = (size_t)BasisPolynomial::factorial_ratio(order+n, order);
    if (lower_bound_offset >= 0) { // default is -1
      int omit_order = order - lower_bound_offset - 1;
      if (omit_order >= 0)
	num_terms -= (size_t)BasisPolynomial::factorial_ratio(
	  omit_order+n, omit_order);
    }
    num_terms /= (size_t)BasisPolynomial::factorial(n);
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
compute_numerical_moments(size_t num_moments, const RealVector& coeffs,
			  RealVector& moments)
{
  // computes and stores the following moments:
  // > mean     (1st raw moment)
  // > variance (2nd central moment)
  // > skewness (3rd standardized moment)
  // > kurtosis (4th standardized moment with offset to eliminate "excess")

  // current support for this implementation: can't be open-ended since we
  // employ a specific combination of raw, central, and standardized moments
  if (num_moments < 1 || num_moments > 4) {
    PCerr << "Error: unsupported number of moments requested in Polynomial"
	  << "Approximation::compute_numerical_moments()" << std::endl;
    abort_handler(-1);
  }

  if (moments.length() == num_moments)
    moments = 0.;
  else
    moments.size(num_moments); // init to 0

  size_t i, j, num_pts = coeffs.length();
  const RealVector& t1_wts = driverRep->type1_weight_sets();
  if (t1_wts.length() != num_pts) {
    PCerr << "Error: mismatch in array lengths between integration driver "
	  << "weights ("  << t1_wts.length() << ") and coefficients ("
	  << num_pts << ") in PolynomialApproximation::compute_numerical_"
	  << "moments()." << std::endl;
    abort_handler(-1);
  }

  // estimate 1st raw moment (mean)
  Real& mean = moments[0];
  for (i=0; i<num_pts; ++i)
    mean += t1_wts[i] * coeffs[i];

  // estimate central moments 2 through num_moments
  Real centered_fn, pow_fn;
  for (i=0; i<num_pts; ++i) {
    pow_fn = centered_fn = coeffs[i] - mean;
    for (j=1; j<num_moments; ++j) {
      pow_fn     *= centered_fn;
      moments[j] += t1_wts[i] * pow_fn;
    }
  }

  // standardize third and higher central moments, if present
  standardize_moments(moments);
}


void PolynomialApproximation::
compute_numerical_moments(size_t num_moments, const RealVector& t1_coeffs,
			  const RealMatrix& t2_coeffs, RealVector& moments)
{
  // computes and stores the following moments:
  // > mean     (1st raw moment)
  // > variance (2nd central moment)
  // > skewness (3rd standardized moment)
  // > kurtosis (4th standardized moment with offset to eliminate "excess")

  // current support for this implementation: can't be open-ended since we
  // employ a specific combination of raw, central, and standardized moments
  if (num_moments < 1 || num_moments > 4) {
    PCerr << "Error: unsupported number of moments requested in Polynomial"
	  << "Approximation::compute_numerical_moments()" << std::endl;
    abort_handler(-1);
  }

  if (moments.length() == num_moments)
    moments = 0.;
  else
    moments.size(num_moments); // init to 0

  size_t i, j, k, num_pts = t1_coeffs.length();
  const RealVector& t1_wts = driverRep->type1_weight_sets();
  const RealMatrix& t2_wts = driverRep->type2_weight_sets();
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
  Real& mean = moments[0];
  for (i=0; i<num_pts; ++i) {
    mean += t1_wts[i] * t1_coeffs[i];
    const Real* coeff2_i = t2_coeffs[i];
    const Real*  t2_wt_i = t2_wts[i];
    for (j=0; j<numVars; ++j)
      mean += coeff2_i[j] * t2_wt_i[j];
  }

  // estimate central moments 2 through num_moments
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

  // standardize third and higher central moments, if present
  standardize_moments(moments);
}


void PolynomialApproximation::standardize_moments(RealVector& moments)
{
  size_t num_moments = moments.length();
  if (num_moments <= 2) return;

  const Real& var = moments[1];
  if (var > 0.) {
    // standardized moment k is E[((X-mu)/sigma)^k] = E[(X-mu)^k]/sigma^k
    Real pow_fn = var, std_dev = std::sqrt(var);
    for (size_t i=2; i<num_moments; ++i)
      { pow_fn *= std_dev; moments[i] /= pow_fn; }
    // offset the fourth standardized moment to eliminate excess kurtosis
    if (num_moments > 3)
      moments[3] -= 3.;
  }
  else
    PCerr << "Warning: skewness and kurtosis cannot be standardized due to "
	  << "non-positive variance.\n         Skipping standardization."
	  << std::endl;
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

} // namespace Pecos
