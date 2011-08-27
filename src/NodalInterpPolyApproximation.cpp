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
#include "SparseGridDriver.hpp"

//#define DEBUG


namespace Pecos {


/** Overloaded version supporting tensor-product quadrature. */
Real NodalInterpPolyApproximation::
tensor_product_value(const RealVector& x, const RealVector& exp_t1_coeffs,
		     const RealMatrix& exp_t2_coeffs, const UShort2DArray& key)
{
  Real tp_val = 0.; size_t i, num_colloc_pts = key.size();
  switch (configOptions.useDerivs) {
  case false:
    for (size_t i=0; i<num_colloc_pts; ++i)
      tp_val += exp_t1_coeffs[i] * type1_interpolant_value(x, key[i]);
    break;
  case true:
    for (size_t i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      tp_val += exp_t1_coeffs[i] * type1_interpolant_value(x, key_i);
      const Real* coeff2_i = exp_t2_coeffs[i];
      for (size_t j=0; j<numVars; ++j)
	tp_val += coeff2_i[j] * type2_interpolant_value(x, j, key_i);
    }
    break;
  }
  return tp_val;
}


/** Overloaded version supporting Smolyak sparse grids. */
Real NodalInterpPolyApproximation::
tensor_product_value(const RealVector& x, const RealVector& exp_t1_coeffs,
		     const RealMatrix& exp_t2_coeffs,
		     const UShortArray& sm_index, const UShort2DArray& key,
		     const SizetArray& colloc_index)
{
  Real tp_val = 0.; size_t i, num_colloc_pts = key.size();
  switch (configOptions.useDerivs) {
  case false:
    for (i=0; i<num_colloc_pts; ++i)
      tp_val += exp_t1_coeffs[colloc_index[i]] *
	        type1_interpolant_value(x, key[i], sm_index);
    break;
  case true:
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      tp_val += exp_t1_coeffs[colloc_index[i]] *
	        type1_interpolant_value(x, key_i, sm_index);
      const Real* coeff2_i = exp_t2_coeffs[colloc_index[i]];
      for (size_t j=0; j<numVars; ++j)
	tp_val += coeff2_i[j] * type2_interpolant_value(x, j, key_i, sm_index);
    }
    break;
  }
  return tp_val;
}


/** Overloaded version supporting tensor-product quadrature. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_gradient(const RealVector& x, const RealVector& exp_t1_coeffs,
			const RealMatrix& exp_t2_coeffs,
			const UShort2DArray& key)
{
  if (tpGradient.length() != numVars)
    tpGradient.sizeUninitialized(numVars);
  tpGradient = 0.;

  size_t i, j, num_colloc_pts = key.size();
  switch (configOptions.useDerivs) {
  case false:
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real&     coeff1_i = exp_t1_coeffs[i];
      for (j=0; j<numVars; ++j) // compute ith contribution to tpGradient[j]
	tpGradient[j] += coeff1_i * type1_interpolant_gradient(x, j, key_i);
    }
    break;
  case true:
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real&     coeff1_i = exp_t1_coeffs[i];
      const Real*     coeff2_i = exp_t2_coeffs[i]; // column vector
      for (j=0; j<numVars; ++j) { // ith contribution to jth grad component
	tpGradient[j] += coeff1_i * type1_interpolant_gradient(x, j, key_i);
	for (size_t k=0; k<numVars; ++k) // type2 interpolant for kth grad comp
	  tpGradient[j] +=
	    coeff2_i[k] * type2_interpolant_gradient(x, j, k, key_i);
      }
    }
    break;
  }
  return tpGradient;
}


/** Overloaded version supporting Smolyak sparse grids. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_gradient(const RealVector& x, const RealVector& exp_t1_coeffs,
			const RealMatrix& exp_t2_coeffs,
			const UShortArray& sm_index, const UShort2DArray& key,
			const SizetArray& colloc_index)
{
  if (tpGradient.length() != numVars)
    tpGradient.sizeUninitialized(numVars);
  tpGradient = 0.;

  size_t i, j, num_colloc_pts = key.size();
  switch (configOptions.useDerivs) {
  case false:
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real&     coeff1_i = exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] +=
	  coeff1_i * type1_interpolant_gradient(x, j, key_i, sm_index);
    }
    break;
  case true:
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real&     coeff1_i = exp_t1_coeffs[colloc_index[i]];
      const Real*     coeff2_i = exp_t2_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j) { // ith contribution to jth grad component
	tpGradient[j] +=
	  coeff1_i * type1_interpolant_gradient(x, j, key_i, sm_index);
	for (size_t k=0; k<numVars; ++k) // type2 interpolant for kth grad comp
	  tpGradient[j] +=
	    coeff2_i[k] * type2_interpolant_gradient(x, j, k, key_i, sm_index);
      }
    }
    break;
  }
  return tpGradient;
}


/** Overloaded version supporting tensor-product quadrature. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_gradient(const RealVector& x, const RealVector& exp_t1_coeffs,
			const RealMatrix& exp_t2_coeffs,
			const UShort2DArray& key, const SizetArray& dvv)
{
  size_t i, j, deriv_index, num_colloc_pts = key.size(),
    num_deriv_vars = dvv.size();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.sizeUninitialized(num_deriv_vars);
  tpGradient = 0.;

  switch (configOptions.useDerivs) {
  case false:
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i   = key[i];
      const Real&        coeff1_i = exp_t1_coeffs[i];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += coeff1_i *
	  type1_interpolant_gradient(x, deriv_index, key_i);
      }
    }
    break;
  case true:
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real&     coeff1_i = exp_t1_coeffs[i];
      const Real*     coeff2_i = exp_t2_coeffs[i]; // column vector
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += coeff1_i *
	  type1_interpolant_gradient(x, deriv_index, key_i);
	for (size_t k=0; k<numVars; ++k)
	  tpGradient[j] +=
	    coeff2_i[k] * type2_interpolant_gradient(x, deriv_index, k, key_i);
      }
    }
    break;
  }
  return tpGradient;
}


/** Overloaded version supporting Smolyak sparse grids. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_gradient(const RealVector& x, const RealVector& exp_t1_coeffs,
			const RealMatrix& exp_t2_coeffs,
			const UShortArray& sm_index, const UShort2DArray& key,
			const SizetArray& colloc_index, const SizetArray& dvv)
{
  size_t i, j, k, deriv_index, num_colloc_pts = key.size(),
    num_deriv_vars = dvv.size();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.sizeUninitialized(num_deriv_vars);
  tpGradient = 0.;

  switch (configOptions.useDerivs) {
  case false:
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i   = key[i];
      const Real&        coeff_i = exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += coeff_i *
	  type1_interpolant_gradient(x, deriv_index, key_i, sm_index);
      }
    }
    break;
  case true:
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      const Real&     coeff1_i = exp_t1_coeffs[colloc_index[i]];
      const Real*     coeff2_i = exp_t2_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += coeff1_i *
	  type1_interpolant_gradient(x, deriv_index, key_i, sm_index);
	for (size_t k=0; k<numVars; ++k)
	  tpGradient[j] += coeff2_i[k] *
	    type2_interpolant_gradient(x, deriv_index, k, key_i, sm_index);
      }
    }
    break;
  }
  return tpGradient;
}


/** Overloaded all_variables version supporting tensor-product quadrature. */
Real NodalInterpPolyApproximation::
tensor_product_mean(const RealVector& x, const UShort2DArray& key)
{
  TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
  const Real2DArray&   colloc_wts_1d
    = tpq_driver->type1_collocation_weights_array();

  std::vector<BasisPolynomial>& basis_0 = polynomialBasis[0];
  size_t i, j, num_colloc_pts = key.size(); Real tp_mean = 0.;
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    Real prod_i = 1.;
    for (j=0; j<numVars; ++j)
      prod_i *= (randomVarsKey[j]) ? colloc_wts_1d[j][key_i[j]] :
	basis_0[j].type1_value(x[j], key_i[j]);
    tp_mean += expansionType1Coeffs[i] * prod_i;
  }
  return tp_mean;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
Real NodalInterpPolyApproximation::
tensor_product_mean(const RealVector& x, const UShortArray& sm_index,
		    const UShort2DArray& key, const SizetArray& colloc_index)
{
  SparseGridDriver*  ssg_driver = (SparseGridDriver*)driverRep;
  const Real3DArray& colloc_wts_1d
    = ssg_driver->type1_collocation_weights_array();

  Real tp_mean = 0.;
  size_t i, j, num_colloc_pts = key.size();
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    Real prod_i = 1.;
    for (j=0; j<numVars; ++j)
      prod_i *= (randomVarsKey[j]) ? colloc_wts_1d[sm_index[j]][j][key_i[j]] :
	polynomialBasis[sm_index[j]][j].type1_value(x[j], key_i[j]);
    tp_mean += expansionType1Coeffs[colloc_index[i]] * prod_i;
  }
  return tp_mean;
}


/** Overloaded all_variables version supporting tensor-product quadrature. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_mean_gradient(const RealVector& x, const UShort2DArray& key,
			     const SizetArray& dvv)
{
  TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
  const Real2DArray&   colloc_wts_1d
    = tpq_driver->type1_collocation_weights_array();

  size_t i, j, k, deriv_index, num_colloc_pts = key.size(),
    num_deriv_vars = dvv.size(),
    cntr = 0; // insertions in cntr order w/i tpCoeffGrads
  if (tpMeanGrad.length() != num_deriv_vars)
    tpMeanGrad.sizeUninitialized(num_deriv_vars);
  tpMeanGrad = 0.;
  std::vector<BasisPolynomial>& basis_0 = polynomialBasis[0];
  SizetList::iterator it;
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    // Error check for required data
    if (randomVarsKey[deriv_index] && !configOptions.expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	    << "InterpPolyApproximation::tensor_product_mean_gradient()."
	    << std::endl;
      abort_handler(-1);
    }
    else if (!randomVarsKey[deriv_index] && !configOptions.expansionCoeffFlag) {
      PCerr << "Error: expansion coefficients not defined in NodalInterp"
	    << "PolyApproximation::tensor_product_mean_gradient()" << std::endl;
      abort_handler(-1);
    }
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      Real wt_prod_j = 1.;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	{ k = *it; wt_prod_j *= colloc_wts_1d[k][key_j[k]]; }
      if (randomVarsKey[deriv_index]) {
	// --------------------------------------------------------------------
	// derivative of All var expansion w.r.t. random var (design insertion)
	// --------------------------------------------------------------------
	Real Lsa_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it)
	  { k = *it; Lsa_j *= basis_0[k].type1_value(x[k], key_j[k]); }
	tpMeanGrad[i] += wt_prod_j * Lsa_j * expansionType1CoeffGrads(cntr, j);
      }
      else {
	// ---------------------------------------------------------------------
	// deriv of All var expansion w.r.t. nonrandom var (design augmentation)
	// ---------------------------------------------------------------------
	Real dLsa_j_dsa_i = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    basis_0[k].type1_gradient(x[k], key_j[k]) :
	    basis_0[k].type1_value(x[k],    key_j[k]);
	}
	tpMeanGrad[i] += expansionType1Coeffs[j] * wt_prod_j * dLsa_j_dsa_i;
      }
    }
    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      cntr++;
  }

  return tpMeanGrad;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_mean_gradient(const RealVector& x, const UShortArray& sm_index,
			     const UShort2DArray& key,
			     const SizetArray& colloc_index,
			     const SizetArray& dvv)
{
  SparseGridDriver*  ssg_driver = (SparseGridDriver*)driverRep;
  const Real3DArray& colloc_wts_1d
    = ssg_driver->type1_collocation_weights_array();

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
  size_t i, j, k, deriv_index,
    cntr = 0, // insertions in cntr order w/i tpCoeffGrads
    num_deriv_vars = dvv.size(), num_colloc_pts = key.size();
  if (tpMeanGrad.length() != num_deriv_vars)
    tpMeanGrad.sizeUninitialized(num_deriv_vars);
  tpMeanGrad = 0.;
  SizetList::iterator it;
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    // Error check for required data
    if (randomVarsKey[deriv_index] && !configOptions.expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	    << "InterpPolyApproximation::tensor_product_mean_gradient()."
	    << std::endl;
      abort_handler(-1);
    }
    else if (!randomVarsKey[deriv_index] && !configOptions.expansionCoeffFlag) {
      PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	    << "Approximation::tensor_product_mean_gradient()" << std::endl;
      abort_handler(-1);
    }
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      Real wt_prod_j = 1.;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it) {
	k = *it;
	wt_prod_j *= colloc_wts_1d[sm_index[k]][k][key_j[k]];
      }
      if (randomVarsKey[deriv_index]) {
	// --------------------------------------------------------------------
	// derivative of All var expansion w.r.t. random var (design insertion)
	// --------------------------------------------------------------------
	Real Lsa_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  Lsa_j *= polynomialBasis[sm_index[k]][k].type1_value(x[k], key_j[k]);
	}
	tpMeanGrad[i] += wt_prod_j * Lsa_j *
	  expansionType1CoeffGrads(cntr, colloc_index[j]);
      }
      else {
	// ---------------------------------------------------------------------
	// deriv of All var expansion w.r.t. nonrandom var (design augmentation)
	// ---------------------------------------------------------------------
	Real dLsa_j_dsa_i = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    polynomialBasis[sm_index[k]][k].type1_gradient(x[k], key_j[k]) :
	    polynomialBasis[sm_index[k]][k].type1_value(x[k],    key_j[k]);
	}
	tpMeanGrad[i]
	  += expansionType1Coeffs[colloc_index[j]] * wt_prod_j * dLsa_j_dsa_i;
      }
    }
    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      cntr++;
  }

  return tpMeanGrad;
}


/** Overloaded all_variables version supporting tensor-product quadrature. */
Real NodalInterpPolyApproximation::
tensor_product_covariance(const RealVector& x, const UShort2DArray& key,
			  const RealVector& exp_coeffs_2)
{
  TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
  const Real2DArray&   colloc_wts_1d
    = tpq_driver->type1_collocation_weights_array();

  Real tp_covar = 0.;
  size_t i, j, k, num_colloc_pts = key.size();
  Real mean1 = 0., mean2 = 0.;
  SizetList::iterator it;
  std::vector<BasisPolynomial>& basis_0 = polynomialBasis[0];
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    Real wt_prod_i = 1., Ls_prod_i = 1.;
    for (k=0; k<numVars; ++k)
      if (randomVarsKey[k])
	wt_prod_i *= colloc_wts_1d[k][key_i[k]];
      else
	Ls_prod_i *= basis_0[k].type1_value(x[k], key_i[k]);
    const Real& exp_coeff_i = expansionType1Coeffs[i];
    Real wt_Ls_prod_i = wt_prod_i * Ls_prod_i;
    mean1 += exp_coeff_i     * wt_Ls_prod_i;
    mean2 += exp_coeffs_2[i] * wt_Ls_prod_i;
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      // to include the ij-th term, colloc pts xi_i must be the same as xi_j
      // for the ran var subset.  In this case, wt_prod_i may be reused.
      bool include = true;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	if (key_i[*it] != key_j[*it])
	  { include = false; break; }
      if (include) {
	Real Ls_prod_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  Ls_prod_j *= basis_0[k].type1_value(x[k], key_j[k]);
	}
	tp_covar += wt_prod_i * Ls_prod_i * Ls_prod_j * exp_coeff_i
	  * exp_coeffs_2[j];
      }
    }
  }
  tp_covar -= mean1*mean2;
  return tp_covar;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
Real NodalInterpPolyApproximation::
tensor_product_covariance(const RealVector& x, const UShortArray& sm_index,
			  const UShort2DArray& key,
			  const SizetArray& colloc_index,
			  const RealVector& exp_coeffs_2)
{
  SparseGridDriver*  ssg_driver = (SparseGridDriver*)driverRep;
  const Real3DArray& colloc_wts_1d
    = ssg_driver->type1_collocation_weights_array();

  size_t i, j, k, index, num_colloc_pts = key.size();
  Real mean1 = 0., mean2 = 0.;
  SizetList::iterator it;
  Real tp_covar = 0.;
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = key[i];
    Real wt_prod_i = 1., Ls_prod_i = 1.;
    for (k=0; k<numVars; ++k)
      if (randomVarsKey[k])
	wt_prod_i *= colloc_wts_1d[sm_index[k]][k][key_i[k]];
      else
	Ls_prod_i *= polynomialBasis[sm_index[k]][k].type1_value(x[k],key_i[k]);
    index = colloc_index[i];
    const Real& exp_coeff_i = expansionType1Coeffs[index];
    Real wt_Ls_prod_i = wt_prod_i * Ls_prod_i;
    mean1 += exp_coeff_i         * wt_Ls_prod_i;
    mean2 += exp_coeffs_2[index] * wt_Ls_prod_i;
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      // to include the ij-th term, colloc pts xi_i must be the same as xi_j
      // for the ran var subset.  In this case, wt_prod_i may be reused.
      bool include = true;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	if (key_i[*it] != key_j[*it])
	  { include = false; break; }
      if (include) {
	Real Ls_prod_j = 1.;
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  Ls_prod_j *=
	    polynomialBasis[sm_index[k]][k].type1_value(x[k], key_j[k]);
	}
	tp_covar += wt_prod_i * Ls_prod_i * Ls_prod_j * exp_coeff_i
	  * exp_coeffs_2[colloc_index[j]];
      }
    }
  }
  tp_covar -= mean1*mean2;
  return tp_covar;
}


/** Overloaded all_variables version supporting tensor-product quadrature. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_variance_gradient(const RealVector& x, const UShort2DArray& key,
				 const SizetArray& dvv)
{
  TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
  const Real2DArray&   colloc_wts_1d
    = tpq_driver->type1_collocation_weights_array();

  Real mean = 0.;
  size_t i, j, k, l, deriv_index, num_colloc_pts = key.size(),
    num_deriv_vars = dvv.size(),
    cntr = 0; // insertions in cntr order w/i expansionType1CoeffGrads
  if (tpVarianceGrad.length() != num_deriv_vars)
    tpVarianceGrad.sizeUninitialized(num_deriv_vars);
  tpVarianceGrad = 0.;
  SizetList::iterator it;
  std::vector<BasisPolynomial>& basis_0 = polynomialBasis[0];
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (randomVarsKey[deriv_index] && !configOptions.expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	    << "InterpPolyApproximation::tensor_product_variance_gradient()."
	    << std::endl;
      abort_handler(-1);
    }
    Real mean_grad_i = 0.;
    // first loop of double sum
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      // compute wt_prod_j and Lsa_j
      Real wt_prod_j = 1., Lsa_j = 1., dLsa_j_dsa_i = 1.;
      for (k=0; k<numVars; ++k)
	if (randomVarsKey[k])
	  wt_prod_j *= colloc_wts_1d[k][key_j[k]];
	else
	  Lsa_j *= basis_0[k].type1_value(x[k], key_j[k]);
      // update mean (once) and mean_grad_i
      if (i == 0)
	mean += expansionType1Coeffs[j] * wt_prod_j * Lsa_j;
      if (randomVarsKey[deriv_index])
	mean_grad_i += wt_prod_j * Lsa_j * expansionType1CoeffGrads(cntr, j);
      else {
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    basis_0[k].type1_gradient(x[k], key_j[k]) :
	    basis_0[k].type1_value(x[k],    key_j[k]);
	}
	mean_grad_i += wt_prod_j * dLsa_j_dsa_i * expansionType1Coeffs[j];
      }
      // second loop of double sum
      for (k=0; k<num_colloc_pts; ++k) {
	const UShortArray& key_k = key[k];
	// to include the jk-th term, colloc pts xi_j must be the same as xi_k
	// for the random var subset.  In this case, wt_prod_j may be reused.
	bool include = true;
	for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	  if (key_j[*it] != key_k[*it])
	    { include = false; break; }
	if (include) {
	  Real Lsa_k = 1.;
	  for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	    l = *it;
	    Lsa_k *= basis_0[l].type1_value(x[l], key_k[l]);
	  }
	  if (randomVarsKey[deriv_index])
	    // ---------------------------------------------------------------
	    // deriv of All var expansion w.r.t. random var (design insertion)
	    // ---------------------------------------------------------------
	    tpVarianceGrad[i] += wt_prod_j * Lsa_j * Lsa_k *
	      (expansionType1Coeffs[j] * expansionType1CoeffGrads(cntr, k) +
	       expansionType1CoeffGrads(cntr, j) * expansionType1Coeffs[k]);
	  else {
	    // ---------------------------------------------------------------
	    // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	    // ---------------------------------------------------------------
	    Real dLsa_k_dsa_i = 1.;
	    for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it){
	      l = *it;
	      dLsa_k_dsa_i *= (l == deriv_index) ?
		basis_0[l].type1_gradient(x[l], key_k[l]):
		basis_0[l].type1_value(x[l],    key_k[l]);
	    }
	    tpVarianceGrad[i] +=
	      wt_prod_j * expansionType1Coeffs[j] * expansionType1Coeffs[k] *
	      (Lsa_j * dLsa_k_dsa_i + dLsa_j_dsa_i * Lsa_k);
	  }
	}
      }
    }

    // subtract 2 mu dmu/ds
    tpVarianceGrad[i] -= 2. * mean * mean_grad_i;

    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      cntr++;
  }
  return tpVarianceGrad;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_variance_gradient(const RealVector& x,
				 const UShortArray& sm_index,
				 const UShort2DArray& key,
				 const SizetArray& colloc_index,
				 const SizetArray& dvv)
{
  SparseGridDriver*  ssg_driver = (SparseGridDriver*)driverRep;
  const Real3DArray& colloc_wts_1d
    = ssg_driver->type1_collocation_weights_array();

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
  size_t i, j, k, l, deriv_index,
    cntr = 0, // insertions in cntr order w/i expansionType1CoeffGrads
    num_deriv_vars = dvv.size(), num_colloc_pts = key.size();
  if (tpVarianceGrad.length() != num_deriv_vars)
    tpVarianceGrad.sizeUninitialized(num_deriv_vars);
  tpVarianceGrad = 0.;
  Real mean = 0.;
  SizetList::iterator it;
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (randomVarsKey[deriv_index] && !configOptions.expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	    << "InterpPolyApproximation::tensor_product_variance_gradient()."
	    << std::endl;
      abort_handler(-1);
    }
    Real mean_grad_i = 0.;
    // first loop of double sum
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = key[j];
      // compute wt_prod_j and Lsa_j
      Real wt_prod_j = 1., Lsa_j = 1., dLsa_j_dsa_i = 1.;
      for (k=0; k<numVars; ++k)
	if (randomVarsKey[k])
	  wt_prod_j *= colloc_wts_1d[sm_index[k]][k][key_j[k]];
	else
	  Lsa_j *= polynomialBasis[sm_index[k]][k].type1_value(x[k], key_j[k]);
      // update mean (once) and mean_grad_i
      if (i == 0)
	mean += expansionType1Coeffs[colloc_index[j]] * wt_prod_j * Lsa_j;
      if (randomVarsKey[deriv_index])
	mean_grad_i
	  += wt_prod_j * Lsa_j * expansionType1CoeffGrads(cntr,colloc_index[j]);
      else {
	for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	  k = *it;
	  dLsa_j_dsa_i *= (k == deriv_index) ?
	    polynomialBasis[sm_index[k]][k].type1_gradient(x[k], key_j[k]) :
	    polynomialBasis[sm_index[k]][k].type1_value(x[k],    key_j[k]);
	}
	mean_grad_i
	  += wt_prod_j * dLsa_j_dsa_i * expansionType1Coeffs[colloc_index[j]];
      }
      // second loop of double sum
      for (k=0; k<num_colloc_pts; ++k) {
	const UShortArray& key_k = key[k];
	// to include the jk-th term, colloc pts xi_j must be the same as xi_k
	// for the random var subset.  In this case, wt_prod_j may be reused.
	bool include = true;
	for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	  if (key_j[*it] != key_k[*it])
	    { include = false; break; }
	if (include) {
	  Real Lsa_k = 1.;
	  for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	    l = *it;
	    Lsa_k *= polynomialBasis[sm_index[l]][l].type1_value(x[l],key_k[l]);
	  }
	  if (randomVarsKey[deriv_index])
	    // ---------------------------------------------------------------
	    // deriv of All var expansion w.r.t. random var (design insertion)
	    // ---------------------------------------------------------------
	    tpVarianceGrad[i] += wt_prod_j * Lsa_j * Lsa_k *
	      (expansionType1Coeffs[colloc_index[j]] *
	       expansionType1CoeffGrads(cntr,colloc_index[k]) +
	       expansionType1CoeffGrads(cntr,colloc_index[j]) *
	       expansionType1Coeffs[colloc_index[k]]);
	  else {
	    // ---------------------------------------------------------------
	    // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	    // ---------------------------------------------------------------
	    Real dLsa_k_dsa_i = 1.;
	    for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it){
	      l = *it;
	      dLsa_k_dsa_i *= (l == deriv_index) ?
		polynomialBasis[sm_index[l]][l].type1_gradient(x[l], key_k[l]):
		polynomialBasis[sm_index[l]][l].type1_value(x[l],    key_k[l]);
	    }
	    tpVarianceGrad[i] +=
	      wt_prod_j * expansionType1Coeffs[colloc_index[j]] *
	      expansionType1Coeffs[colloc_index[k]] *
	      (Lsa_j * dLsa_k_dsa_i + dLsa_j_dsa_i * Lsa_k);
	  }
	}
      }
    }

    // subtract 2 mu dmu/ds
    tpVarianceGrad[i] -= 2. * mean * mean_grad_i;

    if (randomVarsKey[deriv_index]) // deriv w.r.t. design var insertion
      cntr++;
  }
  return tpVarianceGrad;
}


Real NodalInterpPolyApproximation::value(const RealVector& x)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    return tensor_product_value(x, expansionType1Coeffs, expansionType2Coeffs,
				tpq_driver->collocation_key());
    break;
  }
  case SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    Real approx_val = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	approx_val += sm_coeffs[i] *
	  tensor_product_value(x, expansionType1Coeffs, expansionType2Coeffs,
			       sm_mi[i], colloc_key[i], colloc_indices[i]);
    return approx_val;
    break;
  }
  }
}


Real NodalInterpPolyApproximation::stored_value(const RealVector& x)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_value(x, storedExpType1Coeffs, storedExpType2Coeffs,
				storedCollocKey[0]);
    break;
  case SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    Real approx_val = 0.;
    size_t i, num_smolyak_indices = storedSmolyakCoeffs.size();
    for (i=0; i<num_smolyak_indices; ++i)
      if (storedSmolyakCoeffs[i])
	approx_val += storedSmolyakCoeffs[i] *
	  tensor_product_value(x, storedExpType1Coeffs, storedExpType2Coeffs,
			       storedSmolyakMultiIndex[i], storedCollocKey[i],
			       storedCollocIndices[i]);
    return approx_val;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::gradient(const RealVector& x)
{
  // this could define a default_dvv and call gradient(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::gradient()" << std::endl;
    abort_handler(-1);
  }

  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    return tensor_product_gradient(x, expansionType1Coeffs,
      expansionType2Coeffs, tpq_driver->collocation_key());
    break;
  }
  case SPARSE_GRID: {
    if (approxGradient.length() != numVars)
      approxGradient.sizeUninitialized(numVars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient(x,
	  expansionType1Coeffs, expansionType2Coeffs, sm_mi[i],
	  colloc_key[i], colloc_indices[i]);
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
stored_gradient(const RealVector& x)
{
  // this could define a default_dvv and call gradient(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::gradient()" << std::endl;
    abort_handler(-1);
  }

  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    return tensor_product_gradient(x, storedExpType1Coeffs,
      storedExpType2Coeffs, storedCollocKey[0]);
    break;
  case SPARSE_GRID: {
    if (approxGradient.length() != numVars)
      approxGradient.sizeUninitialized(numVars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    size_t i, j, num_smolyak_indices = storedSmolyakCoeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = storedSmolyakCoeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient(x,
	  storedExpType1Coeffs, storedExpType2Coeffs,
	  storedSmolyakMultiIndex[i], storedCollocKey[i],
	  storedCollocIndices[i]);
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
gradient(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::gradient()" << std::endl;
    abort_handler(-1);
  }

  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    return tensor_product_gradient(x, expansionType1Coeffs,
      expansionType2Coeffs, tpq_driver->collocation_key(), dvv);
    break;
  }
  case SPARSE_GRID: {
    size_t num_deriv_vars = dvv.size();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.sizeUninitialized(num_deriv_vars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = tensor_product_gradient(x,
	  expansionType1Coeffs, expansionType2Coeffs, sm_mi[i],
	  colloc_key[i], colloc_indices[i], dvv);
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
  if (!configOptions.expansionCoeffFlag) {
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

  if (numericalMoments.empty())
    numericalMoments.sizeUninitialized(4); // standard mode
  Real& mean = numericalMoments[0]; mean = 0.;
  const RealVector& t1_wts = driverRep->type1_weight_sets();
  switch (configOptions.useDerivs) {
  case false:
    for (size_t i=0; i<numCollocPts; ++i)
      mean += expansionType1Coeffs[i] * t1_wts[i];
    break;
  case true: {
    size_t i, j;
    const RealMatrix& t2_wts = driverRep->type2_weight_sets();
    for (i=0; i<numCollocPts; ++i) {
      mean += expansionType1Coeffs[i] * t1_wts[i];
      const Real* coeff2_i = expansionType2Coeffs[i];
      const Real*  t2_wt_i = t2_wts[i];
      for (j=0; j<numVars; ++j)
	mean += coeff2_i[j] * t2_wt_i[j];
    }
    break;
  }
  }

  return mean;
}


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves evaluating the
    expectation over this subset. */
Real NodalInterpPolyApproximation::mean(const RealVector& x)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    return tensor_product_mean(x, tpq_driver->collocation_key());
    break;
  }
  case SPARSE_GRID: {
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    if (numericalMoments.empty())
      numericalMoments.sizeUninitialized(2); // all_variables mode
    Real& mean = numericalMoments[0]; mean = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	mean += sm_coeffs[i] *
	  tensor_product_mean(x, sm_mi[i], colloc_key[i], colloc_indices[i]);
    return mean;
    break;
  }
  }
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
  if (!configOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	  << "InterpPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }

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
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    return tensor_product_mean_gradient(x, tpq_driver->collocation_key(), dvv);
    break;
  }
  case SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (meanGradient.length() != num_deriv_vars)
      meanGradient.sizeUninitialized(num_deriv_vars);
    meanGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      if (coeff) {
	const RealVector& tpm_grad
	  = tensor_product_mean_gradient(x, sm_mi[i], colloc_key[i],
					 colloc_indices[i], dvv);
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
  if (numericalMoments.empty())
    numericalMoments.sizeUninitialized(4); // standard mode
  numericalMoments[1] = covariance(this);
  return numericalMoments[1];
}


/** In this case, a subset of the expansion variables are random
    variables and the variance of the expansion involves summations
    over this subset. */
Real NodalInterpPolyApproximation::variance(const RealVector& x)
{
  if (numericalMoments.empty())
    numericalMoments.sizeUninitialized(2); // all_variables mode
  numericalMoments[1] = covariance(x, this);
  return numericalMoments[1];
}


Real NodalInterpPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
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
  NodalInterpPolyApproximation* nip_approx_2
    = (NodalInterpPolyApproximation*)poly_approx_2;
  // Note: compute_statistics() in dakota/src/NonDExpansion.C orders calls
  //       to reduce repetition in moment calculations.
  Real  mean_1 = mean(), mean_2 = nip_approx_2->mean();
  const RealVector& t1_coeffs_2 = nip_approx_2->expansionType1Coeffs;
  const RealVector& t1_wts      = driverRep->type1_weight_sets();
  Real covar = 0.; size_t i, j;
  switch (configOptions.useDerivs) {
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
  return covar;
}


/** In this case, a subset of the expansion variables are random
    variables and the variance of the expansion involves summations
    over this subset. */
Real NodalInterpPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  const RealVector& t1_coeffs_2
    = ((NodalInterpPolyApproximation*)poly_approx_2)->expansionType1Coeffs;
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    return tensor_product_covariance(x, tpq_driver->collocation_key(),
				     t1_coeffs_2);
    break;
  }
  case SPARSE_GRID:
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, num_smolyak_indices = sm_coeffs.size(); Real covar = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	covar += sm_coeffs[i] *
	  tensor_product_covariance(x, sm_mi[i], colloc_key[i],
				    colloc_indices[i], t1_coeffs_2);
    return covar;
    break;
  }
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  The
    mixed derivative case (some design variables are inserted and some
    are augmented) requires no special treatment. */
const RealVector& NodalInterpPolyApproximation::variance_gradient()
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::variance_gradient()" << std::endl;
    abort_handler(-1);
  }
  if (!configOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	  << "InterpPolyApproximation::variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  const RealVector& t1_wts = driverRep->type1_weight_sets();
  size_t i, j, num_deriv_vars = expansionType1CoeffGrads.numRows();
  if (varianceGradient.length() != num_deriv_vars)
    varianceGradient.sizeUninitialized(num_deriv_vars);
  varianceGradient = 0.;

  Real mean = 0.;
  for (i=0; i<numCollocPts; ++i)
    mean += expansionType1Coeffs[i] * t1_wts[i];
  for (i=0; i<numCollocPts; ++i) {
    Real term_i = 2. * (expansionType1Coeffs[i] - mean) * t1_wts[i];
    for (j=0; j<num_deriv_vars; ++j)
      varianceGradient[j] += term_i * expansionType1CoeffGrads(j,i);
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
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::variance_gradient()" << std::endl;
    abort_handler(-1);
  }

  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    return
      tensor_product_variance_gradient(x, tpq_driver->collocation_key(), dvv);
    break;
  }
  case SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (varianceGradient.length() != num_deriv_vars)
      varianceGradient.sizeUninitialized(num_deriv_vars);
    varianceGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
    const UShort2DArray& sm_mi          = ssg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs      = ssg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key     = ssg_driver->collocation_key();
    const Sizet2DArray&  colloc_indices = ssg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      if (coeff) {
	const RealVector& tpv_grad
	  = tensor_product_variance_gradient(x, sm_mi[i], colloc_key[i],
					     colloc_indices[i], dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  varianceGradient[j] += coeff * tpv_grad[j];
      }
    }
    return varianceGradient;
    break;
  }
}

} // namespace Pecos
