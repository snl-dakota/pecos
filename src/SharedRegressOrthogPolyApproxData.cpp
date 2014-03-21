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
  else 
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
