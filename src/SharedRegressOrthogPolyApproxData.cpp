/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedRegressOrthogPolyApproxData
//- Description:  Implementation code for SharedRegressOrthogPolyApproxData class
//-               
//- Owner:        John Jakeman

#include "SharedRegressOrthogPolyApproxData.hpp"
#include "pecos_global_defs.hpp"
#include "SurrogateData.hpp"


namespace Pecos {

void SharedRegressOrthogPolyApproxData::allocate_data()
{
  if (expConfigOptions.expCoeffsSolnApproach == ORTHOG_LEAST_INTERPOLATION) {
    // clear history from previous expansion; new pts -> new least interpolant
    approxOrder.clear();   // for update_approx_order() -> exp combination logic
    multiIndex.clear();    // for reuse check in ROPA::least_interpolation()
    sobolIndexMap.clear(); // for update_component_sobol()

    if (expConfigOptions.vbdFlag && expConfigOptions.vbdOrderLimit == 1)
      allocate_main_sobol(); // main effects only

    PCout << "Orthogonal polynomial approximation of least order\n";
  }
  else if (expConfigOptions.expBasisType == ADAPTED_BASIS) {
    // *** TO DO ***: build up from tensor index sets (generalized sparse grid)
  }
  else if (expConfigOptions.expBasisType == TENSOR_PRODUCT_BASIS /* || 
	   expConfigOptions.expBasisType == DEFAULT_BASIS && tensor_grid */ ) {
    // for tensor samples, seems like a good idea to use tensor candidate exp
    // *** only for CS candidate? or true basis for Least sq as well? ***

    inflate_scalar(approxOrder, numVars); // promote scalar to vector, if needed
    tensor_product_multi_index(approxOrder, multiIndex);
    allocate_component_sobol(multiIndex);
    approxOrderPrev = approxOrder;// Note: defer if update_exp_form needed later
  }
  else // use default total-order for unstructured data sets w/o adaptation
    SharedOrthogPolyApproxData::allocate_data();
}


void SharedRegressOrthogPolyApproxData::
update_approx_order(unsigned short new_order)
{
  if (approxOrder.empty() || new_order > approxOrder[0])
    approxOrder.assign(numVars, new_order);
}


void SharedRegressOrthogPolyApproxData::
pack_polynomial_data(const RealVector& c_vars, const UShortArray& mi,
		     bool add_val,  double* pack_val,  size_t& pv_cntr,
		     bool add_grad, double* pack_grad, size_t& pg_cntr)
{
  if (add_val)
    { pack_val[pv_cntr] = multivariate_polynomial(c_vars, mi); ++pv_cntr; }
  if (add_grad) {
    const RealVector& mvp_grad
      = multivariate_polynomial_gradient_vector(c_vars, mi);
    for (size_t j=0; j<numVars; ++j, ++pg_cntr)
      pack_grad[pg_cntr] = mvp_grad[j];
  }
}


void SharedRegressOrthogPolyApproxData::
pack_response_data(const SurrogateDataResp& sdr,
		   bool add_val,  double* pack_val,  size_t& pv_cntr,
		   bool add_grad, double* pack_grad, size_t& pg_cntr)
{
  if (add_val)
    { pack_val[pv_cntr] = sdr.response_function(); ++pv_cntr; }
  if (add_grad) {
    const RealVector& resp_grad = sdr.response_gradient();
    for (size_t j=0; j<numVars; ++j, ++pg_cntr)
      pack_grad[pg_cntr] = resp_grad[j];
  }
}

} // namespace Pecos
