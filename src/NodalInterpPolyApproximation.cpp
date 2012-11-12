/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        NodalInterpPolyApproximation
//- Description:  Implementation code for NodalInterpPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "NodalInterpPolyApproximation.hpp"
#include "TensorProductDriver.hpp"
#include "CombinedSparseGridDriver.hpp"
#include "InterpolationPolynomial.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG

namespace Pecos {


void NodalInterpPolyApproximation::allocate_expansion_coefficients()
{
  size_t num_deriv_vars = surrData.num_derivative_variables();
  if (expConfigOptions.expansionCoeffFlag) {
    if (expansionType1Coeffs.length() != numCollocPts)
      expansionType1Coeffs.sizeUninitialized(numCollocPts);
    if ( basisConfigOptions.useDerivs &&
	 ( expansionType2Coeffs.numRows() != num_deriv_vars ||
	   expansionType2Coeffs.numCols() != numCollocPts ) )
      expansionType2Coeffs.shapeUninitialized(num_deriv_vars, numCollocPts);
  }
  if ( expConfigOptions.expansionCoeffGradFlag &&
       ( expansionType1CoeffGrads.numRows() != num_deriv_vars ||
	 expansionType1CoeffGrads.numCols() != numCollocPts ) )
    expansionType1CoeffGrads.shapeUninitialized(num_deriv_vars, numCollocPts);

  // checking numCollocPts is insufficient due to anisotropy --> changes in
  // anisotropic weights could move points around without changing the total.
  //bool update_exp_form =
  //  ( (expConfigOptions.expansionCoeffFlag &&
  //     expansionType1Coeffs.length()      != numCollocPts) ||
  //    (expConfigOptions.expansionCoeffGradFlag &&
  //     expansionType1CoeffGrads.numCols() != numCollocPts ) );
}


void NodalInterpPolyApproximation::compute_expansion_coefficients()
{
  if (surrData.anchor()) {
    if (expConfigOptions.expansionCoeffFlag) {
      expansionType1Coeffs[0] = surrData.anchor_function();
      if (basisConfigOptions.useDerivs)
	Teuchos::setCol(surrData.anchor_gradient(), 0, expansionType2Coeffs);
    }
    if (expConfigOptions.expansionCoeffGradFlag)
      Teuchos::setCol(surrData.anchor_gradient(), 0, expansionType1CoeffGrads);
  }

  size_t index = 0, offset = (surrData.anchor()) ? 1 : 0;
  for (int i=offset; i<numCollocPts; ++i, ++index) {
    if (expConfigOptions.expansionCoeffFlag) {
      expansionType1Coeffs[i] = surrData.response_function(index);
      // Note: gradients from DAKOTA already scaled in u-space Recast
      if (basisConfigOptions.useDerivs)
	Teuchos::setCol(surrData.response_gradient(index), i,
			expansionType2Coeffs);
    }
    if (expConfigOptions.expansionCoeffGradFlag)
      Teuchos::setCol(surrData.response_gradient(index), i,
		      expansionType1CoeffGrads);
  }

  computedMean = computedVariance = 0;
}


void NodalInterpPolyApproximation::store_coefficients()
{
  if (expConfigOptions.expansionCoeffFlag) {
    storedExpType1Coeffs   = expansionType1Coeffs;
    if (basisConfigOptions.useDerivs)
      storedExpType2Coeffs = expansionType2Coeffs;
  }
  if (expConfigOptions.expansionCoeffGradFlag)
    storedExpType1CoeffGrads = expansionType1CoeffGrads;

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    storedCollocKey.resize(1); storedLevMultiIndex.resize(1);
    storedCollocKey[0]     = tpq_driver->collocation_key();
    storedLevMultiIndex[0] = tpq_driver->level_index();
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    storedLevMultiIndex = csg_driver->smolyak_multi_index();
    storedLevCoeffs     = csg_driver->smolyak_coefficients();
    storedCollocKey     = csg_driver->collocation_key();
    storedCollocIndices = csg_driver->collocation_indices();
    break;
  }
  }
}


void NodalInterpPolyApproximation::combine_coefficients(short combine_type)
{
#ifdef DEBUG
  PCout << "Original type1 expansion coefficients prior to combination:\n";
  write_data(PCout, expansionType1Coeffs);
#endif // DEBUG

  // update expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads} by adding or
  // multiplying stored expansion evaluated at current collocation points
  size_t i, j, offset = 0, num_pts = surrData.size();
  bool anchor_pt = surrData.anchor();
  if (anchor_pt) { offset = 1; ++num_pts; }
  Real lf_val, discrep_val;
  for (i=0; i<num_pts; ++i) {
    const RealVector& c_vars = (anchor_pt && i == 0) ?
      surrData.anchor_continuous_variables() :
      surrData.continuous_variables(i-offset);
    if (combine_type == MULT_COMBINE) { // eval once for both Coeffs/CoeffGrads
      discrep_val = stored_value(c_vars);
      lf_val = expansionType1Coeffs[i]; // copy prior to update
    }
    if (expConfigOptions.expansionCoeffFlag) {
      // split up type1/type2 contribs so increments are performed properly
      if (combine_type == ADD_COMBINE)
	expansionType1Coeffs[i] += stored_value(c_vars);
      else if (combine_type == MULT_COMBINE)
	expansionType1Coeffs[i] *= discrep_val;
      if (basisConfigOptions.useDerivs) {
	const RealVector& discrep_grad
	  = stored_gradient_basis_variables(c_vars);
	Real* exp_t2_coeffs_i = expansionType2Coeffs[i];
	size_t num_deriv_vars = discrep_grad.length();
	if (combine_type == ADD_COMBINE)
	  for (j=0; j<num_deriv_vars; ++j)
	    exp_t2_coeffs_i[j] += discrep_grad[j];
	else if (combine_type == MULT_COMBINE)
	  // hf = lf*discrep --> dhf/dx = dlf/dx*discrep + lf*ddiscrep/dx
	  for (j=0; j<num_deriv_vars; ++j)
	    exp_t2_coeffs_i[j] = exp_t2_coeffs_i[j] * discrep_val
	                       + discrep_grad[j]    * lf_val;
      }
    }
    if (expConfigOptions.expansionCoeffGradFlag) {
      Real* exp_t1_grad_i = expansionType1CoeffGrads[i];
      const RealVector& discrep_grad
	= stored_gradient_nonbasis_variables(c_vars);
      size_t num_deriv_vars = discrep_grad.length();
      if (combine_type == ADD_COMBINE)
	for (j=0; j<num_deriv_vars; ++j)
	  exp_t1_grad_i[j] += discrep_grad[j];
      else if (combine_type == MULT_COMBINE)
	for (j=0; j<num_deriv_vars; ++j)
	  exp_t1_grad_i[j] = exp_t1_grad_i[j] * discrep_val
	                   + discrep_grad[j]  * lf_val;
    }
  }
#ifdef DEBUG
  PCout << "Updated type1 expansion coefficients following combination:\n";
  write_data(PCout, expansionType1Coeffs);
#endif // DEBUG

  // clear stored data now that it has been combined
  if (expConfigOptions.expansionCoeffFlag) {
    storedExpType1Coeffs.resize(0);
    if (basisConfigOptions.useDerivs) storedExpType2Coeffs.reshape(0,0);
  }
  if (expConfigOptions.expansionCoeffGradFlag)
    storedExpType1CoeffGrads.reshape(0,0);
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    storedCollocKey.clear(); break;
  case COMBINED_SPARSE_GRID:
    storedLevMultiIndex.clear(); storedLevCoeffs.clear();
    storedCollocKey.clear();     storedCollocIndices.clear(); break;
  }

  computedMean = computedVariance = 0;
}


void NodalInterpPolyApproximation::restore_expansion_coefficients()
{
  size_t offset = 0, new_colloc_pts = surrData.size();
  if (surrData.anchor())
    { offset = 1; ++new_colloc_pts; }

  if (expConfigOptions.expansionCoeffFlag) {
    expansionType1Coeffs.resize(new_colloc_pts);
    if (basisConfigOptions.useDerivs) {
      size_t num_deriv_vars = expansionType2Coeffs.numRows();
      expansionType2Coeffs.reshape(num_deriv_vars, new_colloc_pts);
    }
  }
  if (expConfigOptions.expansionCoeffGradFlag) {
    size_t num_deriv_vars = expansionType1CoeffGrads.numRows();
    expansionType1CoeffGrads.reshape(num_deriv_vars, new_colloc_pts);
  }

  size_t index = numCollocPts - offset;
  for (int i=numCollocPts; i<new_colloc_pts; ++i, ++index) {
    if (expConfigOptions.expansionCoeffFlag) {
      expansionType1Coeffs[i] = surrData.response_function(index);
      if (basisConfigOptions.useDerivs)
	Teuchos::setCol(surrData.response_gradient(index), i,
			expansionType2Coeffs);
    }
    if (expConfigOptions.expansionCoeffGradFlag)
      Teuchos::setCol(surrData.response_gradient(index), i,
		      expansionType1CoeffGrads);
  }

  computedMean = computedVariance = 0;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
Real NodalInterpPolyApproximation::
tensor_product_mean(const RealVector&    x,   const UShortArray& lev_index,
		    const UShort2DArray& key, const SizetArray&  colloc_index)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::tensor_product_mean()" << std::endl;
    abort_handler(-1);
  }

  Real tp_mean = 0.;
  size_t i, j, c_index, num_colloc_pts = key.size();
  if (basisConfigOptions.useDerivs) {
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      tp_mean += expansionType1Coeffs[c_index]
	      *  type1_interpolant_value(x, key_i, lev_index, nonRandomIndices)
	      *  type1_weight(key_i, lev_index, randomIndices);
      const Real *t2_coeff_i = expansionType2Coeffs[c_index];
      for (j=0; j<numVars; ++j)
	tp_mean += t2_coeff_i[j]
	  * type2_interpolant_value(x, j, key_i, lev_index, nonRandomIndices)
	  * type2_weight(j, key_i, lev_index, randomIndices);
    }
  }
  else
    for (i=0; i<num_colloc_pts; ++i) {
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      tp_mean += expansionType1Coeffs[c_index]
	      *  type1_interpolant_value(x, key[i], lev_index, nonRandomIndices)
	      *  type1_weight(key[i], lev_index, randomIndices);
    }

  return tp_mean;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_mean_gradient(const RealVector& x, const UShortArray& lev_index,
			     const UShort2DArray& key,
			     const SizetArray& colloc_index,
			     const SizetArray& dvv)
{
  // -------------------------------------------------------------------
  // Mixed variable key:
  //   xi = ran vars, sa = augmented des vars, si = inserted design vars
  // Active variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
  //   mu(sa, si)    = Sum_i r_i(sa, si) wt_prod_i
  //   dmu/ds        = Sum_i dr_i/ds wt_prod_i
  // All variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
  //   mu(sa, si)    = Sum_i r_i(si) Lsa_i wt_prod_i
  //   dmu/dsa       = Sum_i r_i(si) dLsa_i/dsa wt_prod_i
  //   dmu/dsi       = Sum_i dr_i/dsi Lsa_i wt_prod_i
  // -------------------------------------------------------------------
  size_t i, j, k, c_index, deriv_index, cntr = 0, 
    num_deriv_vars = dvv.size(), num_colloc_pts = key.size();
  if (tpMeanGrad.length() != num_deriv_vars)
    tpMeanGrad.sizeUninitialized(num_deriv_vars);
  tpMeanGrad = 0.;
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    // Error check for required data
    if (randomVarsKey[deriv_index] &&
	!expConfigOptions.expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	    << "InterpPolyApproximation::tensor_product_mean_gradient()."
	    << std::endl;
      abort_handler(-1);
    }
    else if (!randomVarsKey[deriv_index] &&
	     !expConfigOptions.expansionCoeffFlag) {
      PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	    << "Approximation::tensor_product_mean_gradient()" << std::endl;
      abort_handler(-1);
    }
    Real& grad_i = tpMeanGrad[i];
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      c_index = (colloc_index.empty()) ? j : colloc_index[j];
      if (randomVarsKey[deriv_index]) {
	// --------------------------------------------------------------------
	// derivative of All var expansion w.r.t. random var (design insertion)
	// --------------------------------------------------------------------
	grad_i += expansionType1CoeffGrads(cntr, c_index)
	       *  type1_interpolant_value(x, key_j, lev_index, nonRandomIndices)
	       *  type1_weight(key_j, lev_index, randomIndices);
	if (basisConfigOptions.useDerivs) {
	  PCerr << "Error: combination of coefficient gradients and "
		<< "use_derivatives in NodalInterpPolyApproximation::"
		<< "tensor_product_mean_gradient()" << std::endl;
	  abort_handler(-1);
	}
      }
      else {
	// ---------------------------------------------------------------------
	// deriv of All var expansion w.r.t. nonrandom var (design augmentation)
	// ---------------------------------------------------------------------
	grad_i += expansionType1Coeffs[c_index]
	  * type1_interpolant_gradient(x, deriv_index, key_j, lev_index,
				       nonRandomIndices)
	  * type1_weight(key_j, lev_index, randomIndices);
	if (basisConfigOptions.useDerivs) {
	  const Real *t2_coeff_j = expansionType2Coeffs[c_index];
	  for (k=0; k<numVars; ++k)
	    grad_i += t2_coeff_j[k]
	      * type2_interpolant_gradient(x, deriv_index, k, key_j, lev_index,
					   nonRandomIndices)
	      * type2_weight(k, key_j, lev_index, randomIndices);
	}
      }
    }
    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      ++cntr;
  }

  return tpMeanGrad;
}


/** Covariance of response functions for same tensor product grid.
    Supports all_variables mode and either interpolation of products
    or product of interpolants.  For interpolation of products, can be
    used for either tensor or sparse interpolants, but for product of
    interpolants, only valid for tensor interpolants. */
Real NodalInterpPolyApproximation::
tensor_product_covariance(const RealVector& x, const UShortArray& lev_index,
			  const UShort2DArray& key,
			  const SizetArray& colloc_index,
			  NodalInterpPolyApproximation* nip_approx_2)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag ||
      !nip_approx_2->expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::tensor_product_covariance()" << std::endl;
    abort_handler(-1);
  }

  /* TO DO: weights are not a function of x.  Consider pre-computing these
     for each TP for reuse among multiple integrations and multiple x.
  if ( !(computedPartialWeights & 1) )
    tensor_product_type1_partial_weights(randomIndices);
  if ( basisConfigOptions.useDerivs && !(computedPartialWeights & 2) )
    tensor_product_type2_partial_weights(randomIndices);
  */

  bool same = (this == nip_approx_2);
  size_t i, j, k, c_index_i, num_colloc_pts = key.size();
  const RealVector& t1_coeffs_2 = nip_approx_2->expansionType1Coeffs;
  const RealMatrix& t2_coeffs_2 = nip_approx_2->expansionType2Coeffs;
  Real tp_covar = 0., t1_coeff_1_mm1, t1_coeff_2_mm2, mean_1, mean_2;
  if (momentInterpType == PRODUCT_OF_INTERPOLANTS_FAST) {
    mean_1 = tensor_product_mean(x, lev_index, key, colloc_index);
    mean_2 = (same) ? mean_1 :
      nip_approx_2->tensor_product_mean(x, lev_index, key, colloc_index);
  }
  else
    { mean_1 = mean(x); mean_2 = (same) ? mean_1 : nip_approx_2->mean(x); }

  switch (momentInterpType) {
  case INTERPOLATION_OF_PRODUCTS:
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      c_index_i = (colloc_index.empty()) ? i : colloc_index[i];
      t1_coeff_1_mm1 = expansionType1Coeffs[c_index_i] - mean_1;
      t1_coeff_2_mm2 = t1_coeffs_2[c_index_i]          - mean_2;
      tp_covar += t1_coeff_1_mm1 * t1_coeff_2_mm2
	       *  type1_interpolant_value(x, key_i, lev_index, nonRandomIndices)
	       *  type1_weight(key_i, lev_index, randomIndices);
#ifdef DEBUG
      PCout << "t1_coeff_1_mm1 = " << t1_coeff_1_mm1 << " t1_coeff_2_mm2 = "
	    << t1_coeff_2_mm2 << " type1_weight() = "
	    << type1_weight(key_i, lev_index, randomIndices)
	    << " type1_interpolant_value() = "
	    << type1_interpolant_value(x, key_i, lev_index, nonRandomIndices)
	    << " sum = " << tp_covar << std::endl;
#endif // DEBUG
      if (basisConfigOptions.useDerivs) {
	const Real *t2_coeff_1 = expansionType2Coeffs[c_index_i],
	           *t2_coeff_2 = t2_coeffs_2[c_index_i];
	for (j=0; j<numVars; ++j)
	  tp_covar +=
	    (t1_coeff_1_mm1 * t2_coeff_2[j] + t1_coeff_2_mm2 * t2_coeff_1[j])
	    * type2_interpolant_value(x, j, key_i, lev_index, nonRandomIndices)
	    * type2_weight(j, key_i, lev_index, randomIndices);
      }
    }
    break;
  // "diagonal" covariance terms for products of same tensor grids
  case PRODUCT_OF_INTERPOLANTS_FAST: case PRODUCT_OF_INTERPOLANTS_FULL: {
    size_t c_index_j; Real t1_wt_Ls_prod_i;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      c_index_i = (colloc_index.empty()) ? i : colloc_index[i];
      t1_coeff_1_mm1  = expansionType1Coeffs[c_index_i] - mean_1;
      t1_wt_Ls_prod_i = type1_weight(key_i, lev_index, randomIndices) *
	type1_interpolant_value(x, key_i, lev_index, nonRandomIndices);
      for (j=0; j<num_colloc_pts; ++j) {
	const UShortArray& key_j = key[j];
	// to include the ij-th term,  basis i must be the same as basis j for
	// the random var subset.  In this case, wt_prod_i may be reused.  Note
	// that it is not necessary to collapse terms with the same random basis
	// subset, since cross term in (a+b)(a+b) = a^2+2ab+b^2 gets included.
	// If terms were collapsed (following eval of non-random portions), the
	// nested loop could be replaced with a single loop to evaluate (a+b)^2.
	if (match_random_key(key_i, key_j)) {
	  c_index_j = (colloc_index.empty()) ? j : colloc_index[j];
	  t1_coeff_2_mm2 = t1_coeffs_2[c_index_j] - mean_2;
	  tp_covar += t1_coeff_1_mm1 * t1_coeff_2_mm2 * t1_wt_Ls_prod_i *
	    type1_interpolant_value(x, key_j, lev_index, nonRandomIndices);
	  /* TO DO
	  if (basisConfigOptions.useDerivs) {
	    const Real *t2_coeff_1i = expansionType2Coeffs[c_index_i],
	               *t2_coeff_2i = t2_coeffs_2[c_index_i];
	    for (j=0; j<numVars; ++j)
	      tp_covar += (t1_coeff_1i_mm1 * t2_coeff_2i[j] +
	                   t1_coeff_2i_mm2 * t2_coeff_1i[j])
		* type2_interpolant_value(x,j,key_i,lev_index,nonRandomIndices)
		* type2_weight(j, key_i, lev_index, randomIndices);
	  }
	  */
	}
      }
    }
    break;
  }
  }

  return tp_covar;
}


/** Covariance of response functions for differing tensor products using a 
    product of interpolants approach.  Needed for sparse interpolants in
    all_variables mode. */
Real NodalInterpPolyApproximation::
tensor_product_covariance(const RealVector& x, const UShortArray& lev_index_1,
			  const UShort2DArray& key_1,
			  const SizetArray& colloc_index_1,
			  const UShortArray& lev_index_2,
			  const UShort2DArray& key_2,
			  const SizetArray& colloc_index_2,
			  NodalInterpPolyApproximation* nip_approx_2)
{
#ifdef DEBUG
  PCout << "tensor_product_covariance for lev_index_1 =\n" << lev_index_1
	<< "lev_index_2 =\n" << lev_index_2 << "key_1 =\n" << key_1
	<< "key_2 =\n" << key_2 << std::endl;
#endif // DEBUG

  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag ||
      !nip_approx_2->expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::tensor_product_covariance()" << std::endl;
    abort_handler(-1);
  }
  if (momentInterpType == INTERPOLATION_OF_PRODUCTS ||
      momentInterpType == PRODUCT_OF_INTERPOLANTS_FAST) {
    PCerr << "Error: mixed tensor product covariance only required for full "
	  << "products of interpolants. " << std::endl;
    abort_handler(-1);
  }

  bool same = (this == nip_approx_2);
  size_t i, j, k, c_index_1i, c_index_2j, num_pts_1 = key_1.size(),
    num_pts_2 = key_2.size();
  const RealVector& t1_coeffs_2 = nip_approx_2->expansionType1Coeffs;
  const RealMatrix& t2_coeffs_2 = nip_approx_2->expansionType2Coeffs;
  Real tp_covar = 0., t1_coeff_1i_mm1, t1_coeff_2j_mm2, t1_Ls_1i,
    basis_prod, basis_prod_k, mean_1 = mean(x),
    mean_2 = (same) ? mean_1 : nip_approx_2->mean(x);

  bool non_zero; unsigned short l1k, l2k; SizetList::iterator it;
  for (i=0; i<num_pts_1; ++i) {
    const UShortArray& key_1i = key_1[i];
    c_index_1i = (colloc_index_1.empty()) ? i : colloc_index_1[i];
    t1_coeff_1i_mm1 = expansionType1Coeffs[c_index_1i] - mean_1;
    t1_Ls_1i = type1_interpolant_value(x, key_1i, lev_index_1,nonRandomIndices);
    for (j=0; j<num_pts_2; ++j) {
      const UShortArray& key_2j = key_2[j];
      // to include ij-th term, expectation of the product of the interpolation
      // polynomials must be nonzero
      if (basis_product(lev_index_1, key_1i, lev_index_2, key_2j, basis_prod)) {
	c_index_2j = (colloc_index_2.empty()) ? j : colloc_index_2[j];
	t1_coeff_2j_mm2 = t1_coeffs_2[c_index_2j] - mean_2;
	tp_covar += basis_prod * t1_coeff_1i_mm1 * t1_coeff_2j_mm2 * t1_Ls_1i *
	  type1_interpolant_value(x, key_2j, lev_index_2, nonRandomIndices);

	  /* TO DO
	  if (basisConfigOptions.useDerivs) {
	    const Real *t2_coeff_1i = expansionType2Coeffs[c_index_i],
	               *t2_coeff_2i = t2_coeffs_2[c_index_i];
	    for (j=0; j<numVars; ++j)
	      tp_covar += (t1_coeff_1i_mm1 * t2_coeff_2i[j] +
	                   t1_coeff_2i_mm2 * t2_coeff_1i[j])
		* type2_interpolant_value(x,j,key_i,lev_index,nonRandomIndices)
		* type2_weight(j, key_i, lev_index, randomIndices);
	  }
	  */
      }
    }
  }

  return tp_covar;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_variance_gradient(const RealVector& x,
				 const UShortArray& lev_index,
				 const UShort2DArray& key,
				 const SizetArray& colloc_index,
				 const SizetArray& dvv)
{
  // -------------------------------------------------------------------
  // Mixed variable key:
  //   xi = ran vars, sa = augmented des vars, si = inserted design vars
  // Active variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
  //   var(sa, si)   = Sum_i Sum_j r_i(sa, si) r_j(sa, si) wt_prod_ij - mu^2
  //   dvar/ds       = Sum_i Sum_j (r_i dr_j/ds + dr_i/ds r_j) wt_prod_ij -
  //                   2 mu dmu/ds
  // All variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
  //   var(sa, si)   = Sum_i Sum_j r_i(si) r_j(si) Lsa_i Lsa_j wt_prod_ij -
  //                   mu^2
  //   dvar/dsa      = Sum_i Sum_j r_i(si) r_j(si) (Lsa_i dLsa_j/dsa +
  //                   dLsa_i/dsa Lsa_j) wt_prod_ij - 2 mu dmu/dsa
  //   dvar/dsi      = Sum_i Sum_j (r_i dr_j/dsi + dr_i/dsi r_j)
  //                   Lsa_i Lsa_j wt_prod_ij - 2 mu dmu/dsi
  // -------------------------------------------------------------------
  size_t i, j, k, l, deriv_index, c_index_j, cntr = 0,
    num_deriv_vars = dvv.size(), num_colloc_pts = key.size();
  if (tpVarianceGrad.length() != num_deriv_vars)
    tpVarianceGrad.sizeUninitialized(num_deriv_vars);
  tpVarianceGrad = 0.;
  Real t1_coeff_j_mm, mean_1 =
    (momentInterpType == PRODUCT_OF_INTERPOLANTS_FAST) ?
    tensor_product_mean(x, lev_index, key, colloc_index) : mean(x);
  const RealVector& mean_grad =
    (momentInterpType == PRODUCT_OF_INTERPOLANTS_FAST) ?
    tensor_product_mean_gradient(x, lev_index, key, colloc_index, dvv) :
    mean_gradient(x, dvv);

  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    // Error check for required data
    if (randomVarsKey[deriv_index] &&
	!expConfigOptions.expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	    << "InterpPolyApproximation::tensor_product_variance_gradient()."
	    << std::endl;
      abort_handler(-1);
    }
    else if (!randomVarsKey[deriv_index] &&
	     !expConfigOptions.expansionCoeffFlag) {
      PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	    << "Approximation::tensor_product_variance_gradient()" << std::endl;
      abort_handler(-1);
    }
    Real& grad_i = tpVarianceGrad[i];
    switch (momentInterpType) {
    case INTERPOLATION_OF_PRODUCTS:
      for (j=0; j<num_colloc_pts; ++j) {
	const UShortArray& key_j = key[j];
	c_index_j     = (colloc_index.empty()) ? j : colloc_index[j];
	t1_coeff_j_mm = expansionType1Coeffs[c_index_j] - mean_1;
	if (randomVarsKey[deriv_index]) {
	  // ---------------------------------------------------------------
	  // deriv of All var expansion w.r.t. random var (design insertion)
	  // ---------------------------------------------------------------
	  // d/dx[(R-mu)^2] = 2(R-mu)(dR/dx - dmu/dx)
	  grad_i += 2. * t1_coeff_j_mm
	    * (expansionType1CoeffGrads(cntr, c_index_j) - mean_grad[i])
	    * type1_interpolant_value(x, key_j, lev_index, nonRandomIndices)
	    * type1_weight(key_j, lev_index, randomIndices);
	  if (basisConfigOptions.useDerivs) {
	    PCerr << "Error: combination of coefficient gradients and "
		  << "use_derivatives in NodalInterpPolyApproximation::"
		  << "tensor_product_variance_gradient()" << std::endl;
	    abort_handler(-1);
	  }
	}
	else {
	  // ---------------------------------------------------------------
	  // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	  // ---------------------------------------------------------------
	  grad_i += t1_coeff_j_mm * t1_coeff_j_mm
	    * type1_interpolant_gradient(x, deriv_index, key_j, lev_index,
					 nonRandomIndices)
	    * type1_weight(key_j, lev_index, randomIndices);
	  if (basisConfigOptions.useDerivs) {
	    const Real *t2_coeff_j = expansionType2Coeffs[c_index_j];
	    for (k=0; k<numVars; ++k)
	      grad_i += 2. * (t1_coeff_j_mm * t2_coeff_j[k])
		* type2_interpolant_gradient(x, deriv_index, k, key_j,
					      lev_index, nonRandomIndices)
		* type2_weight(k, key_j, lev_index, randomIndices);
	  }
	}
      }
      break;
    case PRODUCT_OF_INTERPOLANTS_FAST: {
      size_t c_index_k;
      Real wt_prod_j, Lsa_j, dLsa_j_dsa_i, t1_coeff_k_mm, Lsa_k;
      // first loop of double sum
      for (j=0; j<num_colloc_pts; ++j) {
	const UShortArray& key_j = key[j];
	c_index_j = (colloc_index.empty()) ? j : colloc_index[j];
	// compute wt_prod_j and Lsa_j
	wt_prod_j = type1_weight(key_j, lev_index, randomIndices);
	Lsa_j = type1_interpolant_value(x, key_j, lev_index, nonRandomIndices);
	dLsa_j_dsa_i = type1_interpolant_gradient(x, deriv_index, key_j,
						  lev_index, nonRandomIndices);
	// second loop of double sum
	for (k=0; k<num_colloc_pts; ++k) {
	  const UShortArray& key_k = key[k];
	  c_index_k = (colloc_index.empty()) ? k : colloc_index[k];
	  // to include jk-th term, colloc pts xi_j must be the same as xi_k
	  // for random var subset.  In this case, wt_prod_j may be reused.
	  if (match_random_key(key_j, key_k)) {
	    t1_coeff_j_mm = expansionType1Coeffs[c_index_j] - mean_1;
	    t1_coeff_k_mm = expansionType1Coeffs[c_index_k] - mean_1;
	    Lsa_k = type1_interpolant_value(x,key_k,lev_index,nonRandomIndices);
	    if (randomVarsKey[deriv_index])
	      // ---------------------------------------------------------
	      // deriv of All var exp w.r.t. random var (design insertion)
	      // ---------------------------------------------------------
	      grad_i += wt_prod_j * Lsa_j * Lsa_k *
		( t1_coeff_j_mm * ( expansionType1CoeffGrads(cntr, c_index_k) -
				    mean_grad[i] )
		+ t1_coeff_k_mm * ( expansionType1CoeffGrads(cntr, c_index_j) -
				    mean_grad[i] ) );
	    else {
	      // ---------------------------------------------------------------
	      // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	      // ---------------------------------------------------------------
	      Real dLsa_k_dsa_i
		= type1_interpolant_gradient(x, deriv_index, key_k,
					     lev_index, nonRandomIndices);
	      grad_i += wt_prod_j * t1_coeff_j_mm * t1_coeff_k_mm *
		(Lsa_j * dLsa_k_dsa_i + dLsa_j_dsa_i * Lsa_k);
	    }
	  }
	}
      }
      break;
    }
    case PRODUCT_OF_INTERPOLANTS_FULL: {
      // TO DO
      break;
    }

    }

    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      ++cntr;
  }

  return tpVarianceGrad;
}


Real NodalInterpPolyApproximation::value(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_value(x, expansionType1Coeffs, expansionType2Coeffs,
				tpq_driver->level_index(),
				tpq_driver->collocation_key(), colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    Real approx_val = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	approx_val += sm_coeffs[i] *
	  tensor_product_value(x, expansionType1Coeffs, expansionType2Coeffs,
			       sm_mi[i], colloc_key[i], colloc_index[i]);
    return approx_val;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_gradient_basis_variables(x, expansionType1Coeffs,
      expansionType2Coeffs, tpq_driver->level_index(),
      tpq_driver->collocation_key(), colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    if (approxGradient.length() != numVars)
      approxGradient.sizeUninitialized(numVars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient_basis_variables(x,
	  expansionType1Coeffs, expansionType2Coeffs, sm_mi[i], colloc_key[i],
	  colloc_index[i]);
	for (j=0; j<numVars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_gradient_basis_variables(x, expansionType1Coeffs,
      expansionType2Coeffs, tpq_driver->level_index(),
      tpq_driver->collocation_key(), colloc_index, dvv);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    size_t num_deriv_vars = dvv.size();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.sizeUninitialized(num_deriv_vars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient_basis_variables(x,
	  expansionType1Coeffs, expansionType2Coeffs, sm_mi[i], colloc_key[i],
	  colloc_index[i], dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in NodalInterp"
	  << "PolyApproximation::gradient_nonbasis_variables()" << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_gradient_nonbasis_variables(x,
      expansionType1CoeffGrads, tpq_driver->level_index(),
      tpq_driver->collocation_key(), colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    size_t num_deriv_vars = expansionType1CoeffGrads.numRows();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.sizeUninitialized(num_deriv_vars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient_nonbasis_variables(
	  x, expansionType1CoeffGrads, sm_mi[i], colloc_key[i],colloc_index[i]);
	for (j=0; j<num_deriv_vars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


Real NodalInterpPolyApproximation::stored_value(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not available in "
	  << "NodalInterpPolyApproximation::stored_value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_value(x, storedExpType1Coeffs, storedExpType2Coeffs,
				storedLevMultiIndex[0], storedCollocKey[0],
				colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    Real approx_val = 0.;
    size_t i, num_smolyak_indices = storedLevCoeffs.size();
    for (i=0; i<num_smolyak_indices; ++i)
      if (storedLevCoeffs[i])
	approx_val += storedLevCoeffs[i] *
	  tensor_product_value(x, storedExpType1Coeffs, storedExpType2Coeffs,
			       storedLevMultiIndex[i], storedCollocKey[i],
			       storedCollocIndices[i]);
    return approx_val;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not available in NodalInterpPoly"
	  << "Approximation::stored_gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_gradient_basis_variables(x, storedExpType1Coeffs,
      storedExpType2Coeffs, storedLevMultiIndex[0], storedCollocKey[0],
      colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    if (approxGradient.length() != numVars)
      approxGradient.sizeUninitialized(numVars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    size_t i, j, num_smolyak_indices = storedLevCoeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = storedLevCoeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient_basis_variables(x,
	  storedExpType1Coeffs, storedExpType2Coeffs, storedLevMultiIndex[i],
	  storedCollocKey[i], storedCollocIndices[i]);
	for (j=0; j<numVars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not available in Nodal"
	  << "InterpPolyApproximation::stored_gradient_nonbasis_variables()"
	  << std::endl;
    abort_handler(-1);
  }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_gradient_nonbasis_variables(x,
      storedExpType1CoeffGrads, storedLevMultiIndex[0], storedCollocKey[0],
      colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    size_t i, j, num_smolyak_indices = storedLevCoeffs.size(),
      num_deriv_vars = storedExpType1CoeffGrads.numRows();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.sizeUninitialized(num_deriv_vars);
    approxGradient = 0.;
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = storedLevCoeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient_nonbasis_variables(
	  x, storedExpType1CoeffGrads, storedLevMultiIndex[i],
	  storedCollocKey[i], storedCollocIndices[i]);
	for (j=0; j<num_deriv_vars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


/** In this case, all expansion variables are random variables and the
    mean of the expansion is simply the sum over i of r_i w_i. */
Real NodalInterpPolyApproximation::mean()
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  // TO DO:
  //if (!driverRep->track_unique_product_weights()) {
  //  PCerr << "Error: unique product weights required in "
  //	  << "NodalInterpPolyApproximation::mean()" << std::endl;
  //  abort_handler(-1);
  //}

  Real& mean = numericalMoments[0];
  if ( !(computedMean & 1) ) {
    mean = 0.;
    const RealVector& t1_wts = driverRep->type1_weight_sets();
    if (basisConfigOptions.useDerivs) {
      size_t i, j;
      const RealMatrix& t2_wts = driverRep->type2_weight_sets();
      for (i=0; i<numCollocPts; ++i) {
	mean += expansionType1Coeffs[i] * t1_wts[i];
	const Real* coeff2_i = expansionType2Coeffs[i];
	const Real*  t2_wt_i = t2_wts[i];
	for (j=0; j<numVars; ++j)
	  mean += coeff2_i[j] * t2_wt_i[j];
      }
    }
    else
      for (size_t i=0; i<numCollocPts; ++i)
	mean += expansionType1Coeffs[i] * t1_wts[i];

    computedMean |= 1;
  }

  return mean;
}


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves evaluating the
    expectation over this subset. */
Real NodalInterpPolyApproximation::mean(const RealVector& x)
{
  Real& mean = numericalMoments[0];
  if ( !(computedMean & 1) || !match_nonrandom_vars(x, xPrevMean) ) {
    switch (expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE: {
      TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
      SizetArray colloc_index; // empty -> default indexing
      mean = tensor_product_mean(x, tpq_driver->level_index(),
				 tpq_driver->collocation_key(), colloc_index);
      break;
    }
    case COMBINED_SPARSE_GRID: {
      CombinedSparseGridDriver* csg_driver
	= (CombinedSparseGridDriver*)driverRep;
      const UShort2DArray& sm_mi        = csg_driver->smolyak_multi_index();
      const IntArray&      sm_coeffs    = csg_driver->smolyak_coefficients();
      const UShort3DArray& colloc_key   = csg_driver->collocation_key();
      const Sizet2DArray&  colloc_index = csg_driver->collocation_indices();
      size_t i, num_smolyak_indices = sm_coeffs.size();
      mean = 0.;
      for (i=0; i<num_smolyak_indices; ++i)
	if (sm_coeffs[i])
	  mean += sm_coeffs[i] *
	    tensor_product_mean(x, sm_mi[i], colloc_key[i], colloc_index[i]);
      break;
    }
    }
    computedMean |= 1; xPrevMean = x;
  }

  return mean;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  In
    this case, the derivative of the expectation is the expectation of
    the derivative.  The mixed derivative case (some design variables
    are inserted and some are augmented) requires no special treatment. */
const RealVector& NodalInterpPolyApproximation::mean_gradient()
{
  // d/ds <R> = <dR/ds>

  // Error check for required data
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	  << "InterpPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  if ( !(computedMean & 2) ) {
    const RealVector& t1_wts = driverRep->type1_weight_sets();
    size_t i, j, num_deriv_vars = expansionType1CoeffGrads.numRows();
    if (meanGradient.length() != num_deriv_vars)
      meanGradient.sizeUninitialized(num_deriv_vars);
    meanGradient = 0.;
    for (i=0; i<numCollocPts; ++i) {
      const Real& t1_wt_i = t1_wts[i];
      for (j=0; j<num_deriv_vars; ++j)
	meanGradient[j] += expansionType1CoeffGrads(j,i) * t1_wt_i;
    }
    computedMean |= 2;
  }

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
    are inserted (derivatives are obtained from expansionType1CoeffGrads). */
const RealVector& NodalInterpPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  // if already computed, return previous result
  if ( (computedMean & 2) &&
       match_nonrandom_vars(x, xPrevMeanGrad) ) // && dvv == dvvPrev)
    switch (expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE:           return tpMeanGrad;   break;
    case COMBINED_SPARSE_GRID: return meanGradient; break;
    }

  // compute the gradient of the mean
  computedMean |= 2; xPrevMeanGrad = x;
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_mean_gradient(x, tpq_driver->level_index(),
					tpq_driver->collocation_key(),
					colloc_index, dvv);
    break;
  }
  case COMBINED_SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (meanGradient.length() != num_deriv_vars)
      meanGradient.sizeUninitialized(num_deriv_vars);
    meanGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      if (coeff) {
	const RealVector& tpm_grad
	  = tensor_product_mean_gradient(x, sm_mi[i], colloc_key[i],
					 colloc_index[i], dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  meanGradient[j] += coeff * tpm_grad[j];
      }
    }
    return meanGradient;
    break;
  }
}


/** In this case, all expansion variables are random variables and the
    variance of the expansion is the sum over all but the first term
    of the coefficients squared times the polynomial norms squared. */
Real NodalInterpPolyApproximation::variance()
{
  Real& var = numericalMoments[1];
  if ( !(computedVariance & 1) ) {
    var = covariance(this);
    computedVariance |= 1;
  }
  return var;
}


/** In this case, a subset of the expansion variables are random
    variables and the variance of the expansion involves summations
    over this subset. */
Real NodalInterpPolyApproximation::variance(const RealVector& x)
{
  Real& var = numericalMoments[1];
  if ( !(computedVariance & 1) || !match_nonrandom_vars(x, xPrevVar) ) {
    var = covariance(x, this);
    computedVariance |= 1; xPrevVar = x;
  }
  return var;
}


Real NodalInterpPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  NodalInterpPolyApproximation* nip_approx_2
    = (NodalInterpPolyApproximation*)poly_approx_2;
  bool same = (this == nip_approx_2);

  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag ||
      !nip_approx_2->expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  // TO DO:
  //if (!driverRep->track_unique_product_weights()) {
  //  PCerr << "Error: unique product weights required in "
  //	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
  //  abort_handler(-1);
  //}

  // compute mean_1,mean_2 first, then compute covariance as
  // wt_prod*(coeff1-mean_1)*(coeff2-mean_2) in order to avoid precision
  // loss from computing covariance as <R_i R_j> - \mu_i \mu_j
  // Note: compute_statistics() in dakota/src/NonDExpansion.C orders calls
  //       to reduce repetition in moment calculations.
  Real mean_1 = mean(), mean_2  = (same) ? mean_1 : nip_approx_2->mean();
  const RealVector& t1_coeffs_2 = nip_approx_2->expansionType1Coeffs;
  const RealVector& t1_wts      = driverRep->type1_weight_sets();
  Real covar = 0.; size_t i, j;
  switch (basisConfigOptions.useDerivs) {
  case false: // type1 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
    for (i=0; i<numCollocPts; ++i)
      covar += (expansionType1Coeffs[i] - mean_1) * (t1_coeffs_2[i] - mean_2)
	    *  t1_wts[i];
    break;
  case true: {
    const RealMatrix& t2_coeffs_2 = nip_approx_2->expansionType2Coeffs;
    const RealMatrix& t2_wts      = driverRep->type2_weight_sets();
    for (i=0; i<numCollocPts; ++i) {
      // type1 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      Real coeff1_i_mm1 = expansionType1Coeffs[i] - mean_1,
	  coeff1_2i_mm2 = t1_coeffs_2[i]          - mean_2;
      covar += coeff1_i_mm1 * coeff1_2i_mm2 * t1_wts[i];
      // type2 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      // --> interpolated gradients are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
      const Real *coeff2_i  = expansionType2Coeffs[i],
	         *coeff2_2i = t2_coeffs_2[i], *t2_wt_i = t2_wts[i];
      for (j=0; j<numVars; ++j)
	covar  += (coeff1_i_mm1 * coeff2_2i[j] + coeff1_2i_mm2 * coeff2_i[j])
	       *  t2_wt_i[j];
    }
    break;
  }
  }
  if (same)
    { numericalMoments[1] = covar; computedVariance |= 1; }
  return covar;
}


/** In this case, a subset of the expansion variables are random
    variables and the variance of the expansion involves integration
    over this subset. */
Real NodalInterpPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  NodalInterpPolyApproximation* nip_approx_2
    = (NodalInterpPolyApproximation*)poly_approx_2;
  bool same = (this == nip_approx_2);
  Real covar;

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    covar = tensor_product_covariance(x, tpq_driver->level_index(),
				      tpq_driver->collocation_key(),
				      colloc_index, nip_approx_2);
    break;
  }

  // While we can collapse the Smolyak recursion and combine the weights in the
  // distinct variables case, we cannot do this here for the all_variables case
  // since the non-integrated interpolation polynomial portions are not constant
  // and are coupled with the weight combination.
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    covar = 0.;
    switch (momentInterpType) {
    // For an interpolation of products (which captures product cross-terms),
    // a sum of tensor-product covariances is correct and straightforward.
    case INTERPOLATION_OF_PRODUCTS: case PRODUCT_OF_INTERPOLANTS_FAST:
      for (i=0; i<num_smolyak_indices; ++i)
	if (sm_coeffs[i])
	  covar += sm_coeffs[i] *
	    tensor_product_covariance(x, sm_mi[i], colloc_key[i],
				      colloc_index[i], nip_approx_2);
      break;
    // For a product of interpolants, the rigorous approach must manage the
    // expectation of products of random interpolation polynomials.  The
    // PRODUCT_OF_INTERPOLANTS_TENSOR_MEAN shortcut approach has been shown to
    // be numerically equivalent and less computationally involved, but the
    // equivalence has not been proven.
    case PRODUCT_OF_INTERPOLANTS_FULL:
      update_nonzero_basis_products(sm_mi);
      for (i=0; i<num_smolyak_indices; ++i) {
	int sm_coeff_i = sm_coeffs[i];
	if (sm_coeff_i) {   // one of each diagonal term
#ifdef DEBUG
	  Real tp_covar
	    = tensor_product_covariance(x, sm_mi[i], colloc_key[i],
					colloc_index[i], nip_approx_2);
	  covar += sm_coeff_i * sm_coeff_i * tp_covar;
	  PCout << "Diagonal covar: sm_coeffs[" << i << "] = " << sm_coeff_i
		<< " tp_covar = " << tp_covar << '\n';
#else
	  covar += sm_coeff_i * sm_coeff_i *
	    tensor_product_covariance(x, sm_mi[i], colloc_key[i],
				      colloc_index[i], nip_approx_2);
#endif // DEBUG
	  for (j=0; j<i; ++j)
	    if (sm_coeffs[j]) { // two of each off-diagonal term
#ifdef DEBUG
	      Real tp_covar =
		tensor_product_covariance(x,
		  sm_mi[i], colloc_key[i], colloc_index[i],
		  sm_mi[j], colloc_key[j], colloc_index[j], nip_approx_2);
	      covar += 2. * sm_coeff_i * sm_coeffs[j] * tp_covar;
	      PCout << "Off-diagonal covar: sm_coeffs[" << i << "] = "
		    << sm_coeff_i << " sm_coeffs[" << j << "] = "
		    << sm_coeffs[j] << " tp_covar = " << tp_covar << '\n';
#else
	      covar += 2. * sm_coeff_i * sm_coeffs[j] *
		tensor_product_covariance(x,
		  sm_mi[i], colloc_key[i], colloc_index[i],
		  sm_mi[j], colloc_key[j], colloc_index[j], nip_approx_2);
#endif // DEBUG
	    }
	}
      }
      break;
    }
    break;
  }
  }
  if (same)
    { numericalMoments[1] = covar; computedVariance |= 1; xPrevVar = x; }
  return covar;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  The
    mixed derivative case (some design variables are inserted and some
    are augmented) requires no special treatment. */
const RealVector& NodalInterpPolyApproximation::variance_gradient()
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::variance_gradient()" << std::endl;
    abort_handler(-1);
  }
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	  << "InterpPolyApproximation::variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  if ( !(computedVariance & 2) ) {
    const RealVector& t1_wts = driverRep->type1_weight_sets();
    size_t i, j, num_deriv_vars = expansionType1CoeffGrads.numRows();
    if (varianceGradient.length() != num_deriv_vars)
      varianceGradient.sizeUninitialized(num_deriv_vars);
    varianceGradient = 0.;

    Real mean_1 = mean();
    for (i=0; i<numCollocPts; ++i) {
      Real term_i = 2. * (expansionType1Coeffs[i] - mean_1) * t1_wts[i];
      for (j=0; j<num_deriv_vars; ++j)
	varianceGradient[j] += term_i * expansionType1CoeffGrads(j,i);
    }
    computedVariance |= 2;
  }

  return varianceGradient;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  This function must handle the mixed case, where
    some design/state variables are augmented (and are part of the
    expansion) and some are inserted (derivatives are obtained from
    expansionType1CoeffGrads). */
const RealVector& NodalInterpPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  // if already computed, return previous result
  if ( (computedVariance & 2) &&
       match_nonrandom_vars(x, xPrevVarGrad) ) // && dvv == dvvPrev)
    switch (expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE:           return tpVarianceGrad;   break;
    case COMBINED_SPARSE_GRID: return varianceGradient; break;
    }

  computedVariance |= 2; xPrevVarGrad = x;
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_variance_gradient(x, tpq_driver->level_index(),
					    tpq_driver->collocation_key(),
					    colloc_index, dvv);
    break;
  }
  case COMBINED_SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (varianceGradient.length() != num_deriv_vars)
      varianceGradient.sizeUninitialized(num_deriv_vars);
    varianceGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    switch (momentInterpType) {
    // For an interpolation of products (which captures product cross-terms),
    // a sum of tensor-product covariances is correct and straightforward.
    case INTERPOLATION_OF_PRODUCTS: case PRODUCT_OF_INTERPOLANTS_FAST:
      for (i=0; i<num_smolyak_indices; ++i) {
	int coeff = sm_coeffs[i];
	if (coeff) {
	  const RealVector& tpv_grad
	    = tensor_product_variance_gradient(x, sm_mi[i], colloc_key[i],
					       colloc_index[i], dvv);
	  for (j=0; j<num_deriv_vars; ++j)
	    varianceGradient[j] += coeff * tpv_grad[j];
	}
      }
      break;
    case PRODUCT_OF_INTERPOLANTS_FULL:
      // TO DO
      break;
    }
    return varianceGradient;
    break;
  }
}


void NodalInterpPolyApproximation::
update_nonzero_basis_products(const UShort2DArray& sm_multi_index)
{
  UShortArray max_lev_index(numVars, 0);
  size_t i, j, num_sm_mi = sm_multi_index.size();
  for (i=0; i<num_sm_mi; ++i) {
    const UShortArray& sm_mi_i = sm_multi_index[i];
    for (j=0; j<numVars; ++j)
      if (sm_mi_i[j] > max_lev_index[j])
	max_lev_index[j] = sm_mi_i[j];
  }

  // TO DO: ensure point mode is Gaussian quadrature, not nested rules
  const Real3DArray&    pts_1d = driverRep->collocation_points_1d();
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();

  // set up bookkeeping for random variables with unique interpolants
  bool empty_nz = nonZerosMapIndices.empty();
  SizetList::iterator it;
  if (empty_nz) {
    size_t num_v = randomIndices.size(), num_nz = 0, v1, v2, v1_cntr, v2_cntr;
    unsigned short min_max_lev;
    nonZerosMapIndices.resize(num_v); nonZerosMapIndices[0] = num_nz++;
    it=randomIndices.begin(); ++it; SizetList::iterator it2; bool found;
    for (v1_cntr=1; it!=randomIndices.end(); ++it, ++v1_cntr) {
      v1 = *it;
      for (it2=randomIndices.begin(), v2_cntr=0; it2!=it; ++it2, ++v2_cntr) {
	v2 = *it2; min_max_lev = std::min(max_lev_index[v1], max_lev_index[v2]);
	if (same_basis(min_max_lev, v1, v2)) {
	  nonZerosMapIndices[v1_cntr] = nonZerosMapIndices[v2_cntr];
	  found = true; break;
	}
      }
      if (!found)
	nonZerosMapIndices[v1_cntr] = num_nz++;
    }
    nonZerosMapMaxLevels.resize(num_nz);
    nonZerosMapMaxLevels.assign(num_nz, 0);
    nonZerosMapArray.resize(num_nz);
  }
  InterpolationPolynomial *poly_rep_i, *poly_rep_j;
  UShort2DMultiSet map_key; UShortMultiSet mk1, mk2;
  size_t k, l, v, v_cntr = 0, num_pts_i, num_pts_j, nz_index;
  unsigned short j_start, max_lev_v, start; Real basis_prod_v;
  UShortMultiSet::iterator kit, lit; UShort2DMultiSet::iterator mit;
  for (it=randomIndices.begin(); it!=randomIndices.end(); ++it, ++v_cntr) {
    v = *it;
    max_lev_v = max_lev_index[v];
    nz_index  = nonZerosMapIndices[v_cntr];
    start     = nonZerosMapMaxLevels[nz_index] + 1;
    if (start <= max_lev_v) {
      const RealArray&    pts_1d_v =    pts_1d[max_lev_v][v];
      const RealArray& t1_wts_1d_v = t1_wts_1d[max_lev_v][v];
      UShort2DMultiSetRealMap& non_zeros_map = nonZerosMapArray[nz_index];
      for (i=1; i<=max_lev_v; ++i) {
	// InterpolationPolynomial is polynomialBasis[level][var]
	poly_rep_i = (InterpolationPolynomial*)
	  polynomialBasis[i][v].polynomial_rep();
	num_pts_i = poly_rep_i->interpolation_size();
	j_start = (i<start) ? start : 1;
	for (j=j_start; j<i; ++j) {
	  poly_rep_j = (InterpolationPolynomial*)
	    polynomialBasis[j][v].polynomial_rep();
	  num_pts_j = poly_rep_j->interpolation_size();
	  mk1.clear(); mk1.insert(num_pts_i);
	  mk2.clear(); mk2.insert(num_pts_j);
	  for (k=0; k<num_pts_i; ++k) {
	    kit = mk1.insert(k);
	    map_key.clear(); map_key.insert(mk1);
	    for (l=0; l<num_pts_j; ++l) {
	      if (basis_product_1d(poly_rep_i, poly_rep_j, k, l, pts_1d_v,
				   t1_wts_1d_v, basis_prod_v)) {
		lit = mk2.insert(l); mit = map_key.insert(mk2);
		non_zeros_map[map_key] = basis_prod_v;
		map_key.erase(mit); mk2.erase(lit);
#ifdef DEBUG
		PCout << "basis product var " << v << " lev_1 = " << i
		      << " lev_2 = " << j << " key_1 = " << k
		      << " key_2 = " << l << " basis_prod_v = "
		      << basis_prod_v << std::endl;
#endif // DEBUG
	      }
	    }
	    mk1.erase(kit);
	  }
	}
      }
      // Update max level tracking
      nonZerosMapMaxLevels[nz_index] = max_lev_v;
    }
  }
}


bool NodalInterpPolyApproximation::
basis_product(const UShortArray& lev_index_1, const UShortArray& key_1,
	      const UShortArray& lev_index_2, const UShortArray& key_2,
	      Real& prod)
{
  SizetList::iterator it; size_t v, v_cntr = 0;
  unsigned short l1v, l2v, k1v, k2v;
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  prod = 1.;
  for (it=randomIndices.begin(); it!=randomIndices.end(); ++it, ++v_cntr) {
    v = *it;
    l1v = lev_index_1[v]; l2v = lev_index_2[v]; k1v = key_1[v]; k2v = key_2[v];
    // Embed knowledge of special cases:
    // > if one level index is 0 (single colloc pt for which L=1=constant), then
    //   expected val is normal weight from integration of single Lagrange poly
    // > if basis order is same, then E[L_i L_j] = delta_ij * weight_i by 1/0
    //   properties at Gauss points (this fact also used in
    //   tensor_product_covariance() function over the same tensor grid)
    // > if level indices are nonzero and different (for either nested or
    //   non-nested points), then numerical evaluation is needed.
    if (l1v == 0)
      prod *= t1_wts_1d[l2v][v][k2v];
    else if (l2v == 0)
      prod *= t1_wts_1d[l1v][v][k1v];
    else if (l1v == l2v) { // same level: bases eval as delta_ij using interpPts
      if (k1v == k2v) // 1D integral is 1D collocation weight
	prod *= t1_wts_1d[l1v][v][k1v];
      else            // 1D integral is zero
	return false;
    }
    else { // nonzero differing levels: lookup precomputed product
      UShort2DMultiSetRealMap& non_zeros_map
	= nonZerosMapArray[nonZerosMapIndices[v_cntr]];
      // lookup 1D terms
      UShortMultiSet mk1, mk2; mk1.insert(k1v); mk2.insert(k2v);
      mk1.insert(polynomialBasis[l1v][v].interpolation_size());
      mk2.insert(polynomialBasis[l2v][v].interpolation_size());
      UShort2DMultiSet map_key; map_key.insert(mk1); map_key.insert(mk2);
      UShort2DMultiSetRealMap::iterator it = non_zeros_map.find(map_key);
      if (it == non_zeros_map.end())
	return false; // zeros are not stored
      else
	prod *= it->second;
    }
  }
  return true;
}


void NodalInterpPolyApproximation::
compute_numerical_response_moments(size_t num_moments)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in InterpPoly"
	  << "Approximation::compute_numerical_response_moments()" << std::endl;
    abort_handler(-1);
  }

  if (numericalMoments.length() != num_moments)
    numericalMoments.sizeUninitialized(num_moments);
  if (basisConfigOptions.useDerivs)
    compute_numerical_moments(expansionType1Coeffs, expansionType2Coeffs,
			      driverRep->type1_weight_sets(),
			      driverRep->type2_weight_sets(), numericalMoments);
  else
    compute_numerical_moments(expansionType1Coeffs,
			      driverRep->type1_weight_sets(), numericalMoments);
}


void NodalInterpPolyApproximation::
compute_numerical_expansion_moments(size_t num_moments)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in InterpPoly"
	  << "Approximation::compute_numerical_expansion_moments()"<< std::endl;
    abort_handler(-1);
  }
  if (expansionMoments.length() != num_moments)
    expansionMoments.sizeUninitialized(num_moments);

  size_t i, offset = 0, num_pts = surrData.size();
  bool anchor_pt = surrData.anchor();
  if (anchor_pt) { offset = 1; ++num_pts; }
  RealVector t1_exp(num_pts);
  if (basisConfigOptions.useDerivs) {
    RealMatrix t2_exp(numVars, num_pts);
    for (i=0; i<num_pts; ++i) {
      const RealVector& c_vars = (anchor_pt && i == 0) ?
	surrData.anchor_continuous_variables() :
	surrData.continuous_variables(i-offset);
      t1_exp[i] = value(c_vars);
      Teuchos::setCol(gradient_basis_variables(c_vars), (int)i, t2_exp);
    }
    compute_numerical_moments(t1_exp, t2_exp, driverRep->type1_weight_sets(),
			      driverRep->type2_weight_sets(), expansionMoments);
  }
  else {
    for (i=0; i<num_pts; ++i) {
      const RealVector& c_vars = (anchor_pt && i == 0) ?
	surrData.anchor_continuous_variables() :
	surrData.continuous_variables(i-offset);
      t1_exp[i] = value(c_vars);
    }
    compute_numerical_moments(t1_exp, driverRep->type1_weight_sets(),
			      expansionMoments);
#ifdef DEBUG
    PCout << "Expansion moments type 1 coefficients:\n";
    write_data(PCout, t1_exp);
    PCout << "Expansion moments type 1 weight sets:\n";
    write_data(PCout, driverRep->type1_weight_sets());
#endif // DEBUG
  }
}


/** Computes the variance of component functions. Assumes that all
    subsets of set_value have been computed in advance which will be
    true so long as the partial_variance is called following
    appropriate enumeration of set value  */
void NodalInterpPolyApproximation::
compute_partial_variance(const BitArray& set_value)
{
  Real& variance = partialVariance[sobolIndexMap[set_value]];

  // Compute the partial integral corresponding to set_value
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray&   quad_order = tpq_driver->quadrature_order();
    const UShortArray&    lev_index = tpq_driver->level_index();
    const UShort2DArray& colloc_key = tpq_driver->collocation_key();
    SizetArray colloc_index; // empty -> default indexing
    variance = partial_variance_integral(set_value, quad_order, lev_index,
      colloc_key, colloc_index, expansionType1Coeffs, expansionType2Coeffs);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    const IntArray&            sm_coeffs = csg_driver->smolyak_coefficients();
    const UShort2DArray&        sm_index = csg_driver->smolyak_multi_index();
    const UShort3DArray&      colloc_key = csg_driver->collocation_key();
    const Sizet2DArray&     colloc_index = csg_driver->collocation_indices();
    // Smolyak recursion of anisotropic tensor products
    size_t i, num_smolyak_indices = sm_coeffs.size();
    UShortArray quad_order;
    variance = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i]) {
	csg_driver->level_to_order(sm_index[i], quad_order);
	variance += sm_coeffs[i] * partial_variance_integral(set_value,
	  quad_order, sm_index[i], colloc_key[i], colloc_index[i],
	  expansionType1Coeffs, expansionType2Coeffs);
      }
    break;
  }
  }

  // compute proper subsets and subtract their contributions
  InterpPolyApproximation::compute_partial_variance(set_value);
}


void NodalInterpPolyApproximation::compute_total_sobol_indices()
{
  // Compute the total expansion variance.  For standard mode, the full variance
  // is likely already available, as managed by computedVariance in variance().
  // For all variables mode, we use covariance(this) without passing x for the
  // nonRandomIndices (bypass computedVariance checks by not using variance()).
  Real total_variance = (nonRandomIndices.empty()) ? variance() : // std mode
                        covariance(this);                    // all vars mode

  size_t j; BitArray set_value(numVars);
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray&   quad_order = tpq_driver->quadrature_order();
    const UShortArray&    lev_index = tpq_driver->level_index();
    const UShort2DArray& colloc_key = tpq_driver->collocation_key();
    SizetArray colloc_index; // empty -> default indexing
    for (j=0; j<numVars; ++j) {
      // define set_value that includes all but index of interest
      set_value.set(); set_value[j].flip();
      totalSobolIndices[j] = std::abs(1. -
	total_effects_integral(set_value, quad_order, lev_index, colloc_key,
			       colloc_index, expansionType1Coeffs,
			       expansionType2Coeffs) / total_variance);
    }
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    const IntArray&            sm_coeffs = csg_driver->smolyak_coefficients();
    const UShort2DArray&        sm_index = csg_driver->smolyak_multi_index();
    const UShort3DArray&      colloc_key = csg_driver->collocation_key();
    const Sizet2DArray&     colloc_index = csg_driver->collocation_indices();
    // Smolyak recursion of anisotropic tensor products
    size_t i, num_smolyak_indices = sm_coeffs.size();
    UShortArray quad_order; Real complement_variance;
    // iterate each variable 
    for (j=0; j<numVars; ++j) {
      // define set_value that includes all but index of interest
      set_value.set(); set_value[j].flip();
      complement_variance = 0.;
      for (i=0; i<num_smolyak_indices; ++i)
	if (sm_coeffs[i]) {
	  csg_driver->level_to_order(sm_index[i], quad_order);
	  complement_variance += sm_coeffs[i] *
	    total_effects_integral(set_value, quad_order, sm_index[i],
				   colloc_key[i], colloc_index[i],
				   expansionType1Coeffs, expansionType2Coeffs);
	}
      totalSobolIndices[j] = std::abs(1. - complement_variance/total_variance);
    }
    break;
  }
  }
}

} // namespace Pecos
