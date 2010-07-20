/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
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

void PolynomialApproximation::allocate_arrays()
{
  // size the sensitivity member variables
  if (expansionCoeffFlag) {
    if (sobolIndices.empty())
      sobolIndices.sizeUninitialized((int)std::pow(2.,(int)numVars));
    if (totalSobolIndices.empty())
      totalSobolIndices.sizeUninitialized(numVars);
  }
}


void PolynomialApproximation::allocate_collocation_arrays()
{
  // define mapping from 1:numCollocPts to set of 1d interpolation indices
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  size_t num_smolyak_indices = smolyakMultiIndex.size();
  collocKey.resize(num_smolyak_indices);
  expansionCoeffIndices.resize(num_smolyak_indices);
  UShortArray quad_order(numVars); //, gauss_indices(numVars);
  //IntArray key(2*numVars);
  //unsigned short closed_order_max;
  //ssg_driver->level_to_order_closed_exponential(ssg_level,
  //  closed_order_max);
  //const IntArray& rules = driverRep->integration_rules();
  //const IntArraySizetMap& index_map = ssg_driver->index_map();
  const IntArray& unique_mapping = ssg_driver->unique_index_mapping();
  size_t i, j, cntr = 0;
  for (i=0; i<num_smolyak_indices; ++i) {
    //const UShortArray& sm_index_i = smolyakMultiIndex[i];
    UShort2DArray& colloc_key_i = collocKey[i];
    SizetArray&    coeff_map_i  = expansionCoeffIndices[i];

    ssg_driver->level_to_order(smolyakMultiIndex[i], quad_order);
    tensor_product_multi_index(quad_order, colloc_key_i);
    size_t num_tp_pts = colloc_key_i.size();
    coeff_map_i.resize(num_tp_pts);
    //gauss_indices = 0;
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      coeff_map_i[j] = unique_mapping[cntr];

      /*
      // update expansionCoeffIndices: sparse grid nesting of points may be
      // tracked by demoting each point to its lowest order representation
      // or by promoting each point to the highest order grid.
      for (k=0; k<numVars; k++) {
	switch (rules[k]) {
	// For open weakly-nested rules (Gauss-Hermite, Gauss-Legendre),
	// demote to order=1,index=0 if a center pt for a particular order>1
	case GAUSS_HERMITE: case GAUSS_LEGENDRE:
	  if (numVars > 1 && quad_order[k] > 1 &&
	      gauss_indices[k] == (quad_order[k]-1)/2) // demoted base/index
	    { key[k] = 1;        key[k+numVars] = 0; }
	  else                                    // original base/index
	    { key[k] = quad_order[k]; key[k+numVars] = gauss_indices[k]; }
	  break;
	// For closed nested rules (Clenshaw-Curtis), base is a dummy and
	// index is promoted to the highest order grid
	case CLENSHAW_CURTIS:
	  key[k] = closed_order_max; // promoted base
	  if (sm_index_i[k] == 0)
	    key[k+numVars] = (closed_order_max-1)/2;      // promoted index
	  else {
	    key[k+numVars] = gauss_indices[k];
	    unsigned short delta_w = ssg_level - sm_index_i[k];
	    if (delta_w)
	      key[k+numVars] *= (size_t)pow(2., delta_w); // promoted index
	  }
	  break;
	// For open non-nested rules (Gauss-Laguerre), no modification reqd
        // Note: might need to check for symmetric cases of Gauss-Jacobi
	default: // original base/index
	  key[k] = quad_order[k]; key[k+numVars] = gauss_indices[k];
	  break;
	}
      }
#ifdef DEBUG
      PCout << "lookup key:\n" << key << std::endl;
#endif // DEBUG
      IntArraySizetMap::const_iterator cit = index_map.find(key);
      if (cit == index_map.end()) {
	PCerr << "Error: lookup on sparse grid index map failed in "
	      << "InterpPolyApproximation::find_coefficients()"
	      << std::endl;
	abort_handler(-1);
      }
      else
	coeff_map_i[j] = cit->second;

      // increment the n-dimensional gauss point index set
      increment_indices(gauss_indices, quad_order, true);
      */
#ifdef DEBUG
      PCout << "collocKey[" << i << "][" << j << "]:\n" << colloc_key_i[j]
	    << "expansionCoeffIndices[" << i << "][" << j << "] = "
	    << coeff_map_i[j] << '\n';
#endif // DEBUG
    }
  }
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
			unsigned int max_terms)
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
smolyak_multi_index(UShort2DArray& multi_index, RealArray& coeffs)
{
  // Populate smolyakMultiIndex and smolyakCoeffs.  Identifies
  // use of polynomialBasis[variable][index] based on index 0:num_levels-1.
  // w = q - N = dimension-independent level.  For isotropic,
  //   w + 1 <= |i| <= w + N for i starts at 1 (used for index set defn.)
  //   w - N + 1 <= |j| <= w for j = i - 1 starts at 0 (used for generation)
  // For anisotropic, a weighted linear index set constraint is used.

  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  if (ssg_driver->isotropic()) {
    unsigned short ssg_level = ssg_driver->level();
    // initialize multi_index
    UShortArray levels(numVars, ssg_level);
    total_order_multi_index(levels, multi_index, numVars-1);
    size_t i, j, num_terms = multi_index.size();
    // initialize coeffs
    coeffs.resize(num_terms);
    for (i=0; i<num_terms; i++) {
      const UShortArray& index_set_i = multi_index[i];
      int wpNmi = ssg_level; // w + N - |i| = w - |j|
      for (j=0; j<numVars; j++) // subtract 1-norm of index set
	wpNmi -= index_set_i[j];
      coeffs[i] = std::pow(-1., wpNmi)
	* BasisPolynomial::n_choose_k(numVars - 1, wpNmi);
    }
  }
  else {
    // utilize Pecos wrapper to sgmga_vcn_{ordered,coef}
    Int2DArray pmi;
    ssg_driver->anisotropic_multi_index(pmi, coeffs);
    // copy Int2DArray -> UShort2DArray
    size_t i, j, num_pmi = pmi.size();
    multi_index.resize(num_pmi);
    for (i=0; i<num_pmi; ++i) {
      multi_index[i].resize(numVars);
      for (j=0; j<numVars; ++j)
	multi_index[i][j] = (unsigned short)pmi[i][j];
    }
  }

#ifdef DEBUG
  size_t i, num_terms = coeffs.size();
  PCout << "\nnum Smolyak terms = " << num_terms << '\n';
  for (i=0; i<num_terms; i++)
    PCout << "multi_index[" << i << "]:\n" << multi_index[i]
	  << "coeffs[" << i << "] = " << coeffs[i] << "\n\n";
#endif // DEBUG
}

} // namespace Pecos
