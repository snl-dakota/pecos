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


/** Return the number of terms in a tensor-product expansion.  For
    isotropic and anisotropic expansion orders, calculation of the
    number of expansion terms is straightforward: Prod(p_i + 1). */
size_t PolynomialApproximation::
tensor_product_terms(const UShortArray& order, bool exp_order_offset)
{
  size_t i, n = order.size(), num_terms = 1;
  if (exp_order_offset)
    for (i=0; i<n; ++i)
      num_terms *= order[i] + 1; // PCE: expansion order p
  else
    for (i=0; i<n; ++i)
      num_terms *= order[i];     // SC:  quadrature order m (default)
  return num_terms;
}


void PolynomialApproximation::
tensor_product_multi_index(const UShortArray& order,
			   UShort2DArray& multi_index, bool exp_order_offset)
{
  // rather than inserting into multi_index, go ahead and estimate its length
  // (since its inexpensive) and use term-by-term assignment.
  size_t i, n = order.size(),
    mi_len = tensor_product_terms(order, exp_order_offset);
  if (mi_len != multi_index.size())
    multi_index.resize(mi_len);
  UShortArray indices(n, 0); multi_index[0] = indices;
  for (i=1; i<mi_len; ++i) {
    // increment the n-dimensional index set
    increment_indices(indices, order, !exp_order_offset);
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


/// TO DO: Add overloaded function to support integration over only ran vars
void PolynomialApproximation::compute_numerical_moments(size_t num_moments)
{
  // computes and stores the following moments:
  // > mean     (1st raw moment)
  // > variance (2nd central moment)
  // > skewness (3rd standardized moment)
  // > kurtosis (4th standardized moment with offset to eliminate "excess")

  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in Polynomial"
	  << "Approximation::compute_numerical_moments()" << std::endl;
    abort_handler(-1);
  }
  // current support for this implementation: can't be open-ended since we
  // employ a specific combination of raw, central, and standardized moments
  if (num_moments < 1 || num_moments > 4) {
    PCerr << "Error: unsupported number of moments requested in Polynomial"
	  << "Approximation::compute_numerical_moments()" << std::endl;
    abort_handler(-1);
  }

  numericalMoments.size(num_moments); // init to 0

  size_t i, j, offset = 0, num_pts = dataPoints.size();
  bool anchor_pt = !anchorPoint.is_null();
  const RealVector& wt_sets = driverRep->weight_sets();

  // estimate 1st raw moment (mean)
  Real& mean = numericalMoments[0];
  if (anchor_pt) {
    offset = 1; num_pts += offset;
    mean   = wt_sets[0] * anchorPoint.response_function();
  }
  for (size_t i=offset; i<num_pts; ++i)
    mean += wt_sets[i] * dataPoints[i].response_function();

  // estimate central moments 2 through num_moments
  Real centered_fn, pow_fn;
  if (anchor_pt) {
    pow_fn = centered_fn = anchorPoint.response_function() - mean;
    for (j=1; j<num_moments; ++j)
      { pow_fn *= centered_fn; numericalMoments[j] = wt_sets[0] * pow_fn; }
  }
  for (i=offset; i<num_pts; ++i) {
    pow_fn = centered_fn = dataPoints[i].response_function() - mean;
    for (j=1; j<num_moments; ++j)
      { pow_fn *= centered_fn; numericalMoments[j] += wt_sets[i] * pow_fn; }
  }

  // standardize third and higher central moments, if present
  if (num_moments > 2) {
    // standardized moment k is E[((X-mu)/sigma)^k] = E[(X-mu)^k]/sigma^k
    Real std_dev = std::sqrt(numericalMoments[1]); pow_fn = std_dev*std_dev;
    for (j=2; j<num_moments; ++j)
      { pow_fn *= std_dev; numericalMoments[j] /= pow_fn; }

    // offset the fourth standardized moment to eliminate excess kurtosis
    if (num_moments > 3)
      numericalMoments[3] -= 3.;
  }

  //return numericalMoments;
}


const RealVector& PolynomialApproximation::dimension_decay_rates()
{
  PCerr << "Error: dimension_decay_rates() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return expansionCoeffs; // meaningless return to keep compilers happy
}


void PolynomialApproximation::increment_order()
{
  PCerr << "Error: increment_order() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
}

} // namespace Pecos
