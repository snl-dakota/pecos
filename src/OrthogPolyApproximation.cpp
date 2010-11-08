/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        OrthogPolyApproximation
//- Description:  Implementation code for OrthogPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "OrthogPolyApproximation.hpp"
#include "TensorProductDriver.hpp"
#include "SparseGridDriver.hpp"
#include "CubatureDriver.hpp"
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"
#include "pecos_global_defs.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG


namespace Pecos {

bool OrthogPolyApproximation::
distribution_types(const ShortArray& u_types,
		   const IntArray& int_rules, ShortArray& basis_types,
		   ShortArray& gauss_modes)
{
  bool extra_dist_params = false;
  size_t i, num_vars = u_types.size(), num_rules = int_rules.size();

  // Inactive code below is overkill: just copy values since
  // BasisPolynomial::get_polynomial() does not pass modes to types that
  // don't support them.  These modes allow control of Gauss-Legendre vs.
  // Gauss-Patterson points/weights for Legendre polynomials and
  // Clenshaw-Curtis vs. Fejer points/weights for Chebyshev polynomials.
  if (gauss_modes.size() != num_rules) {
    gauss_modes.resize(num_rules);
    for (i=0; i<num_rules; i++)
      //switch (u_types[i]) {
      //case STD_UNIFORM:
        gauss_modes[i] = (short)int_rules[i];// break;
      //default:
      //  gauss_modes[i] = 0;                   break;
      //}
  }

  if (basis_types.size() != num_vars) {
    basis_types.resize(num_vars);
    for (i=0; i<num_vars; i++) {
      switch (u_types[i]) {
      case STD_NORMAL:
	basis_types[i] = HERMITE;                                  break;
      case STD_UNIFORM: // Legendre or Chebyshev OrthogPolynomial
	// To employ Chebyshev for uniform, have to multiply inner product
	// integrands by the inverse of the weight function (weight fn =
	// 1/sqrt(1-x^2); same as beta PDF/Jacobi poly for alpha=beta=-1/2).
	//basis_types[i] = CHEBYSHEV;                              break;
	basis_types[i] = LEGENDRE;                                 break;
      case STD_EXPONENTIAL:
	basis_types[i] = LAGUERRE;                                 break;
      case STD_BETA:
	basis_types[i] = JACOBI;
	extra_dist_params = true;                                         break;
      case STD_GAMMA:
	basis_types[i] = GENERALIZED_LAGUERRE;
	extra_dist_params = true;                                         break;
      case BOUNDED_NORMAL: case BOUNDED_LOGNORMAL: case LOGNORMAL:
      case LOGUNIFORM: case TRIANGULAR:
      case GUMBEL: case FRECHET: case WEIBULL: case HISTOGRAM_BIN:
	basis_types[i] = NUMERICALLY_GENERATED;
	extra_dist_params = true;                                         break;
      default:
	PCerr << "Error: unsupported u-space type in OrthogPolyApproximation::"
	      << "distribution_types()." << std::endl;
	abort_handler(-1);                                                break;
      }
    }
  }
  else
    for (i=0; i<num_vars; i++)
      if (u_types[i] != STD_NORMAL && u_types[i] != STD_UNIFORM &&
	  u_types[i] != STD_EXPONENTIAL)
	{ extra_dist_params = true; break; }

  return extra_dist_params;
}


void OrthogPolyApproximation::
distribution_basis(const ShortArray& basis_types, const ShortArray& gauss_modes,
		   std::vector<BasisPolynomial>& poly_basis)
{
  size_t i, num_vars = basis_types.size();
  bool modes = (!gauss_modes.empty());
  if (poly_basis.size() != num_vars) {
    poly_basis.resize(num_vars);
    for (i=0; i<num_vars; ++i)
      poly_basis[i] = (modes) ?	BasisPolynomial(basis_types[i], gauss_modes[i])
	                      :	BasisPolynomial(basis_types[i]);

    /*
    // Could reuse objects as in InterpPolyApproximation, but this would require
    // caching dist params and moving code from distribution_parameters() to
    // here, which isn't yet well motivated.
    size_t i, j;
    for (i=0; i<numVars; ++i) {
      // reuse prev instance via shared rep or instantiate new unique instance
      short basis_type_i = basisTypes[i];
      bool found = false;
      for (j=0; j<i; ++j)
	if ( basis_type_i == basisTypes[j] && ( basis_type_i <= LAGUERRE ||
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


void OrthogPolyApproximation::
distribution_parameters(const ShortArray& u_types, const DistributionParams& dp,
			std::vector<BasisPolynomial>& poly_basis)
{
  size_t i, num_vars = u_types.size(), nuv_cntr = 0, lnuv_cntr = 0,
    luuv_cntr = 0, tuv_cntr = 0, beuv_cntr = 0, gauv_cntr = 0, guuv_cntr = 0,
    fuv_cntr = 0, wuv_cntr = 0, hbuv_cntr = 0;
  for (i=0; i<num_vars; i++)
    switch (u_types[i]) {
    case STD_NORMAL:
      ++nuv_cntr; break;
    case STD_UNIFORM:
    case STD_EXPONENTIAL:
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
      PCerr << "Error: unsupported u-space type in OrthogPolyApproximation::"
	    << "distribution_parameters()" << std::endl;
      abort_handler(-1);
      break;
    }
}


int OrthogPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface of multiple dimensions
  if (configOptions.expansionCoeffFlag || configOptions.expansionCoeffGradFlag)
    switch (configOptions.expCoeffsSolnApproach) {
    case QUADRATURE: case CUBATURE: case SPARSE_GRID: case SAMPLING:
      return 1; // quadrature: (int)pow((Real)MIN_GAUSS_PTS, numVars);
      break;
    case REGRESSION:
      // numExpansionTerms is either set from the NonDPolynomialChaos ctor if
      // expansion_terms is specified or computed by the allocate_arrays() call
      // in compute_coefficients() if expansion_order is specified.  The latter
      // case is too late for use of this fn by ApproximationInterface::
      // minimum_samples() in DataFitSurrModel::build_global(), so
      // numExpansionTerms must be calculated for this case.
      return (numExpansionTerms) ?
	numExpansionTerms : total_order_terms(approxOrder);
      break;
    case -1: default: // coefficient import 
      return 0;
      break;
    }
  else
    return 0;
}


void OrthogPolyApproximation::allocate_arrays()
{
  allocate_component_effects();
  allocate_total_effects();

  // Infer expansion formulation from quadrature_order or sparse_grid_level
  // spec, as in SC.  Preserve previous capability (quadrature_order and
  // sparse_grid_level with total-order expansions) for paper results via
  // quadratureExpansion/sparseGridExpansion (compile-time) switches.

  // update_exp_form controls when to update (refinement) and when not to
  // update (subIterator execution) an expansion's multiIndex definition.
  // Simple logic of updating if previous number of points != current number
  // is not robust enough for anisotropic updates --> track using Prev arrays.
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray& quad_order = tpq_driver->quadrature_order();
    bool update_exp_form = (quad_order != quadOrderPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials?

    if (update_exp_form) {
      UShortArray integrand_order(numVars);
      quadrature_order_to_integrand_order(quad_order, integrand_order);
      integrand_order_to_expansion_order(integrand_order, approxOrder);
    }
    PCout << "Orthogonal polynomial approximation order = { ";
    for (size_t i=0; i<numVars; ++i)
      PCout << approxOrder[i] << ' ';
    if (quadratureExpansion == TENSOR_INT_TENSOR_EXP) {
      if (update_exp_form)
	tensor_product_multi_index(approxOrder, multiIndex, true);
      PCout << "} using tensor-product expansion of ";
    }
    else if (quadratureExpansion == TENSOR_INT_TOTAL_ORD_EXP) {
      if (update_exp_form)
	total_order_multi_index(approxOrder, multiIndex);
      PCout << "} using total-order expansion of ";
    }
    else {
      PCerr << "\nError: unsupported setting for quadratureExpansion in "
	    << "OrthogPolyApproximation::allocate_arrays()." << std::endl;
      abort_handler(-1);
    }
    if (update_exp_form)
      numExpansionTerms = multiIndex.size();
    PCout << numExpansionTerms << " terms\n";

    quadOrderPrev = quad_order;
    break;
  }
  case CUBATURE: {
    CubatureDriver* cub_driver = (CubatureDriver*)driverRep;
    //unsigned short cub_int_order = cub_driver->integrand_order();
    //bool update_exp_form = (cub_int_order != cubIntOrderPrev);

    UShortArray integrand_order(numVars, cub_driver->integrand_order());
    integrand_order_to_expansion_order(integrand_order, approxOrder);
    total_order_multi_index(approxOrder, multiIndex);
    numExpansionTerms = multiIndex.size();
    PCout << "Orthogonal polynomial approximation order = { ";
    for (size_t i=0; i<numVars; ++i)
      PCout << approxOrder[i] << ' ';
    PCout << "} using total-order expansion of " << numExpansionTerms
	  << " terms\n";

    //cubIntOrderPrev = cub_int_order;
    break;
  }
  case SPARSE_GRID: {
    const RealVector& aniso_wts  = driverRep->anisotropic_weights();
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    unsigned short    ssg_level  = ssg_driver->level();
    bool update_exp_form
      = (ssg_level != ssgLevelPrev || aniso_wts != ssgAnisoWtsPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials?

    if (update_exp_form) { // compute and output number of terms
      sparse_grid_multi_index(multiIndex);
      numExpansionTerms = multiIndex.size();
    }
    if (sparseGridExpansion == TENSOR_INT_TENSOR_SUM_EXP)
      PCout << "Orthogonal polynomial approximation level = " << ssg_level
	    << " using tensor integration and tensor sum expansion of "
	    << numExpansionTerms << " terms\n";
    else if (sparseGridExpansion == SPARSE_INT_TENSOR_SUM_EXP)
      PCout << "Orthogonal polynomial approximation level = " << ssg_level
	    << " using sparse integration and tensor sum expansion of "
	    << numExpansionTerms << " terms\n";
    else if (sparseGridExpansion == SPARSE_INT_TOTAL_ORD_EXP ||
	     sparseGridExpansion == SPARSE_INT_HEUR_TOTAL_ORD_EXP) {
      PCout << "Orthogonal polynomial approximation order = { ";
      for (size_t i=0; i<numVars; ++i)
	PCout << approxOrder[i] << ' ';
      PCout << "} using sparse integration and total-order expansion of "
	    << numExpansionTerms << " terms\n";
    }
    else {
      PCerr << "Error: unsupported setting for sparseGridExpansion in "
	    << "OrthogPolyApproximation::allocate_arrays()." << std::endl;
      abort_handler(-1);
    }

    ssgLevelPrev = ssg_level; ssgAnisoWtsPrev = aniso_wts;
    break;
  }
  default: { // SAMPLING and REGRESSION
    // Support efficiency in PCBDO for now.
    bool update_exp_form
      = (expansionCoeffs.empty() && expansionCoeffGrads.empty());
    // TO DO: support uniform/adaptive refinement.
    //bool update_exp_form = (approxOrder       != approxOrderPrev ||
    //                        numExpansionTerms != numExpTermsPrev);

    // TO DO: this logic is not reentrant since approxOrder updates
    //        numExpansionTerms and vice versa
    if (!approxOrder.empty()) {
      if (update_exp_form) {
	size_t order_len = approxOrder.size();
	if (order_len != numVars) {
	  if (order_len == 1) {
	    unsigned short order = approxOrder[0];
	    approxOrder.assign(numVars, order);
	  }
	  else {
	    PCerr << "Error: expansion_order specification length does not "
		  << "match number of active variables." << std::endl;
	    abort_handler(-1);
	  }
	}
	total_order_multi_index(approxOrder, multiIndex);
	numExpansionTerms = multiIndex.size();
      }
      PCout << "Orthogonal polynomial approximation order = { ";
      for (size_t i=0; i<numVars; ++i)
	PCout << approxOrder[i] << ' ';
      PCout << "} using total-order expansion of " << numExpansionTerms
	    << " terms\n";
    }
    else if (numExpansionTerms) { // default is 0
      unsigned short order = 0;
      size_t full_order_expansion = 1;
      while (numExpansionTerms > full_order_expansion) {
	++order;
	// do in 2 steps rather than 1 to avoid truncation from int division
	full_order_expansion *= numVars+order;
	full_order_expansion /= order;
      }
      PCout << "Orthogonal polynomial approximation order = " << order;
      if (numExpansionTerms == full_order_expansion)
	PCout << " using full total-order expansion of ";
      else
	PCout << " using partial total-order expansion of ";
      PCout << numExpansionTerms << " terms\n";
      if (update_exp_form) {
	approxOrder.assign(numVars, order);
	total_order_multi_index(approxOrder, multiIndex, -1, numExpansionTerms);
      }
    }
    else {
      PCerr << "Error: bad expansion specification in "
	    << "OrthogPolyApproximation::allocate_arrays()." << std::endl;
      abort_handler(-1);
    }

    //approxOrderPrev = approxOrder; numExpTermsPrev = numExpansionTerms;
    break;
  }
  }

  // now that terms & order are known, shape some arrays.  This is done here,
  // rather than in compute_coefficients(), in order to support array sizing
  // for the data import case.
  if (configOptions.expansionCoeffFlag &&
      expansionCoeffs.length() != numExpansionTerms)
    expansionCoeffs.sizeUninitialized(numExpansionTerms);
  if (configOptions.expansionCoeffGradFlag) {
    const SurrogateDataPoint& sdp
      = (anchorPoint.is_null()) ? dataPoints[0] : anchorPoint;
    size_t num_deriv_vars = sdp.response_gradient().length();
    if (expansionCoeffGrads.numRows() != num_deriv_vars ||
	expansionCoeffGrads.numCols() != numExpansionTerms)
      expansionCoeffGrads.shapeUninitialized(num_deriv_vars, numExpansionTerms);
  }
}


void OrthogPolyApproximation::compute_coefficients()
{
  if (!configOptions.expansionCoeffFlag &&
      !configOptions.expansionCoeffGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in "
	  << "OrthogPolyApproximation::compute_coefficients().\n         "
	  << "Bypassing approximation construction." << std::endl;
    return;
  }

  // For testing of anchorPoint logic:
  //anchorPoint = dataPoints.front();
  //dataPoints.pop_front();

  // anchorPoint, if present, is handled differently for different
  // expCoeffsSolnApproach settings:
  //   SAMPLING:   treat it as another currentPoint
  //   QUADRATURE/CUBATURE/SPARSE_GRID: error
  //   REGRESSION: use equality-constrained least squares
  size_t i, j, num_pts = dataPoints.size();
  if (!anchorPoint.is_null())
    ++num_pts;
  if (!num_pts) {
    PCerr << "Error: nonzero number of sample points required in "
	  << "OrthogPolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }

  // Array sizing can be divided into two parts:
  // > data used in all cases (size in allocate_arrays())
  // > data not used in expansion import case (size here)
  allocate_arrays();
#ifdef DEBUG
  gradient_check();
#endif // DEBUG

  // calculate polynomial chaos coefficients
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    // verify quad_order stencil matches num_pts
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray& quad_order = tpq_driver->quadrature_order();
    if (quad_order.size() != numVars) {
      PCerr << "Error: quadrature order array is not consistent with number of "
	    << "variables (" << numVars << ")\n       in "
	    << "OrthogPolyApproximation::compute_coefficients()." << std::endl;
      abort_handler(-1);
    }
    size_t num_gauss_pts = 1;
    for (i=0; i<numVars; ++i)
      num_gauss_pts *= quad_order[i];
    if (num_pts != num_gauss_pts) {
      PCerr << "Error: number of current points (" << num_pts << ") is not "
	    << "consistent with\n       quadrature data in "
	    << "OrthogPolyApproximation::compute_coefficients()." << std::endl;
      abort_handler(-1);
    }

#ifdef DEBUG
    for (i=0; i<numVars; ++i) {
      OrthogonalPolynomial* poly_rep
	= (OrthogonalPolynomial*)polynomialBasis[i].polynomial_rep();
      for (j=1; j<=quad_order[i]; ++j)
	poly_rep->gauss_check(j);
    }
#endif // DEBUG

    // single expansion integration
    integration_checks();
    integrate_expansion(multiIndex, dataPoints, driverRep->weight_sets(),
			expansionCoeffs, expansionCoeffGrads);
    break;
  }
  case CUBATURE:
    // single expansion integration
    integration_checks();
    integrate_expansion(multiIndex, dataPoints, driverRep->weight_sets(),
			expansionCoeffs, expansionCoeffGrads);
    break;
  case SPARSE_GRID:
    integration_checks();
    switch (sparseGridExpansion) {
    case TENSOR_INT_TENSOR_SUM_EXP: {
      // multiple expansion integration
      // Note: allocate_arrays() calls sparse_grid_multi_index() which uses
      // append_multi_index() to build multiIndex.
      if (configOptions.expansionCoeffFlag)     expansionCoeffs = 0.;
      if (configOptions.expansionCoeffGradFlag) expansionCoeffGrads = 0.;
      SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
      const IntArray& sm_coeffs = ssg_driver->smolyak_coefficients();
      size_t i, num_tensor_grids = tpMultiIndex.size();
      std::vector<SurrogateDataPoint> tp_data_pts;
      RealVector tp_wts, tp_coeffs; RealMatrix tp_coeff_grads;
      bool store_tp = (configOptions.refinementControl == GENERALIZED_SPARSE);
      // loop over tensor-products, forming sub-expansions, and sum them up
      for (i=0; i<num_tensor_grids; ++i) {
	// form tp_data_pts, tp_wts using collocKey et al.
	integration_data(i, tp_data_pts, tp_wts);

	// form tp_multi_index from tpMultiIndexMap
	//map_tensor_product_multi_index(tp_multi_index, i);

	// form tp expansion coeffs
	RealVector& tp_coeffs_i = (store_tp) ? tpExpansionCoeffs[i] : tp_coeffs;
	RealMatrix& tp_grads_i
	  = (store_tp) ? tpExpansionCoeffGrads[i] : tp_coeff_grads;
	integrate_expansion(tpMultiIndex[i], tp_data_pts, tp_wts, tp_coeffs_i,
			    tp_grads_i);

	// sum tensor product coeffs/grads into expansion coeffs/grads
	combine_expansion(tpMultiIndexMap[i], tp_coeffs_i, tp_grads_i,
			  sm_coeffs[i]);
      }
      //if (!reEntrantFlag) {
      //  ssg_driver->clear_smolyak_arrays();
      //  ssg_driver->clear_collocation_arrays();
      //  tpMultiIndex.clear(); tpMultiIndexMap.clear();
      //}
      break;
    }
    default: // SPARSE_INT_*
      // single expansion integration
      integrate_expansion(multiIndex, dataPoints, driverRep->weight_sets(),
			  expansionCoeffs, expansionCoeffGrads);
      break;
    }
    break;
  case REGRESSION:
    regression();
    break;
  case SAMPLING:
    expectation();
    break;
  }
}


void OrthogPolyApproximation::increment_coefficients()
{
  bool err_flag = false;
  size_t last_index;
  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    switch (sparseGridExpansion) {
    case TENSOR_INT_TENSOR_SUM_EXP: {
      last_index = tpMultiIndex.size();
      // size tpMultiIndex and tpExpansion{Coeffs,CoeffGrads}
      size_t new_size = last_index+1;
      tpMultiIndex.resize(new_size);
      tpExpansionCoeffs.resize(new_size);
      tpExpansionCoeffGrads.resize(new_size);
      // update tpMultiIndex
      UShortArray quad_order(numVars), int_order(numVars), exp_order(numVars);
      ssg_driver->level_to_order(ssg_driver->trial_index_set(), quad_order);
      quadrature_order_to_integrand_order(quad_order, int_order);
      integrand_order_to_expansion_order(int_order, exp_order);
      tensor_product_multi_index(exp_order, tpMultiIndex[last_index], true);
      // update multiIndex and numExpansionTerms
      append_multi_index(tpMultiIndex[last_index], multiIndex, true);
      numExpansionTerms = multiIndex.size();
      resize_expansion();
      // form tp_data_pts, tp_wts using collocKey et al.
      std::vector<SurrogateDataPoint> tp_data_pts; RealVector tp_wts;
      integration_data(last_index, tp_data_pts, tp_wts);
      // form trial expansion coeffs/grads
      integrate_expansion(tpMultiIndex[last_index], tp_data_pts, tp_wts,
			  tpExpansionCoeffs[last_index],
			  tpExpansionCoeffGrads[last_index]);
      // sum trial expansion into expansionCoeffs/expansionCoeffGrads
      append_expansions(last_index);
      // cleanup
      //if (!reEntrantFlag) {
      //  ssg_driver->clear_smolyak_arrays();
      //  ssg_driver->clear_collocation_arrays();
      //  tpMultiIndex.clear(); tpMultiIndexMap.clear();
      //}
      break;
    }
    default:
      err_flag = true; break;
    }
    break;
  }
  default:
    err_flag = true; break;
  }

  if (err_flag) {
    PCerr << "Error: unsupported grid definition in OrthogPolyApproximation::"
	  << "increment_coefficients()" << std::endl;
    abort_handler(-1);
  }
}


void OrthogPolyApproximation::decrement_coefficients()
{
  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID:
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    switch (sparseGridExpansion) {
    case TENSOR_INT_TENSOR_SUM_EXP: {
      // reset expansion{Coeffs,CoeffGrads}: (set in append_expansions())
      expansionCoeffs     = prevExpCoeffs;
      expansionCoeffGrads = prevExpCoeffGrads;
      // reset multiIndex and numExpansionTerms:
      numExpansionTerms = tpMultiIndexMapRef.back();
      multiIndex.resize(numExpansionTerms); // truncate previous increment
      // reset tensor-product bookkeeping and save restorable data
      savedSmolyakMultiIndex.push_back(ssg_driver->trial_index_set());
      savedTPMultiIndex.push_back(tpMultiIndex.back());
      savedTPMultiIndexMap.push_back(tpMultiIndexMap.back());
      savedTPMultiIndexMapRef.push_back(numExpansionTerms);
      savedTPExpCoeffs.push_back(tpExpansionCoeffs.back());
      savedTPExpCoeffGrads.push_back(tpExpansionCoeffGrads.back());
      tpMultiIndex.pop_back();
      tpMultiIndexMap.pop_back();
      tpMultiIndexMapRef.pop_back();
      tpExpansionCoeffs.pop_back();
      tpExpansionCoeffGrads.pop_back();
      // resize not necessary since not updating expansion on decrement;
      // rather, next increment takes care of this.
      //resize_expansion();
      break;
    }
    }
    break;
  }
}


void OrthogPolyApproximation::restore_coefficients()
{
  size_t last_index;
  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID:
    switch (sparseGridExpansion) {
    case TENSOR_INT_TENSOR_SUM_EXP: {
      // move previous expansion data to current expansion
      last_index = tpMultiIndex.size();
      SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
      UShort2DArray::iterator  sit = std::find(savedSmolyakMultiIndex.begin(),
					       savedSmolyakMultiIndex.end(),
					       ssg_driver->trial_index_set());
      size_t index_star = std::distance(savedSmolyakMultiIndex.begin(), sit);
      savedSmolyakMultiIndex.erase(sit);
      UShort3DArray::iterator   iit = savedTPMultiIndex.begin();
      Sizet2DArray::iterator    mit = savedTPMultiIndexMap.begin();
      SizetArray::iterator      rit = savedTPMultiIndexMapRef.begin();
      RealVectorArray::iterator cit = savedTPExpCoeffs.begin();
      RealMatrixArray::iterator git = savedTPExpCoeffGrads.begin();
      std::advance(iit, index_star); std::advance(mit, index_star);
      std::advance(rit, index_star); std::advance(cit, index_star);
      std::advance(git, index_star);
      tpMultiIndex.push_back(*iit);          savedTPMultiIndex.erase(iit);
      tpMultiIndexMap.push_back(*mit);       savedTPMultiIndexMap.erase(mit);
      tpMultiIndexMapRef.push_back(*rit);    savedTPMultiIndexMapRef.erase(rit);
      tpExpansionCoeffs.push_back(*cit);     savedTPExpCoeffs.erase(cit);
      tpExpansionCoeffGrads.push_back(*git); savedTPExpCoeffGrads.erase(git);
      // update multiIndex and numExpansionTerms
      append_multi_index(tpMultiIndex[last_index], tpMultiIndexMap[last_index],
			 tpMultiIndexMapRef[last_index], multiIndex);
      numExpansionTerms = multiIndex.size();
      resize_expansion();
      // sum trial expansion into expansionCoeffs/expansionCoeffGrads
      append_expansions(last_index);
      break;
    }
    break;
  }
  }
}


void OrthogPolyApproximation::finalize_coefficients()
{
  size_t start_index;
  switch (configOptions.expCoeffsSolnApproach) {
  case SPARSE_GRID:
    switch (sparseGridExpansion) {
    case TENSOR_INT_TENSOR_SUM_EXP: {
      start_index = tpMultiIndex.size();
      // update multiIndex and numExpansionTerms
      UShort3DArray::iterator iit = savedTPMultiIndex.begin();
      Sizet2DArray::iterator  mit = savedTPMultiIndexMap.begin();
      SizetArray::iterator    rit = savedTPMultiIndexMapRef.begin();
      for (; iit!=savedTPMultiIndex.end(); ++iit, ++mit, ++rit)
	append_multi_index(*iit, *mit, *rit, multiIndex);
      numExpansionTerms = multiIndex.size();
      resize_expansion();
      // move previous expansion data to current expansion
      tpMultiIndex.insert(tpMultiIndex.end(),
	savedTPMultiIndex.begin(), savedTPMultiIndex.end());
      tpMultiIndexMap.insert(tpMultiIndexMap.end(),
	savedTPMultiIndexMap.begin(), savedTPMultiIndexMap.end());
      tpMultiIndexMapRef.insert(tpMultiIndexMapRef.end(),
	savedTPMultiIndexMapRef.begin(), savedTPMultiIndexMapRef.end());
      tpExpansionCoeffs.insert(tpExpansionCoeffs.end(),
	savedTPExpCoeffs.begin(), savedTPExpCoeffs.end());
      tpExpansionCoeffGrads.insert(tpExpansionCoeffGrads.end(),
	savedTPExpCoeffGrads.begin(), savedTPExpCoeffGrads.end());
      savedSmolyakMultiIndex.clear();
      savedTPMultiIndex.clear();
      savedTPMultiIndexMap.clear();
      savedTPMultiIndexMapRef.clear();
      savedTPExpCoeffs.clear();
      savedTPExpCoeffGrads.clear();
      // sum remaining trial expansions into expansionCoeffs/expansionCoeffGrads
      append_expansions(start_index);
      break;
    }
    break;
  }
  }
}


void OrthogPolyApproximation::
sparse_grid_multi_index(UShort2DArray& multi_index)
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  const UShort2DArray& sm_multi_index = ssg_driver->smolyak_multi_index();
  size_t i, num_smolyak_indices = sm_multi_index.size();

  switch (sparseGridExpansion) {
  case TENSOR_INT_TENSOR_SUM_EXP: {
    // assemble a complete list of individual polynomial coverage
    // defined from the linear combination of mixed tensor products
    multi_index.clear();
    tpMultiIndex.resize(num_smolyak_indices);
    if (configOptions.refinementControl == GENERALIZED_SPARSE) {
      tpExpansionCoeffs.resize(num_smolyak_indices);
      tpExpansionCoeffGrads.resize(num_smolyak_indices);
    }
    UShortArray quad_order(numVars), integrand_order(numVars),
      expansion_order(numVars);
    for (i=0; i<num_smolyak_indices; ++i) {
      // regenerate i-th expansion_order as collocKey[i] cannot be used in
      // general case (i.e., for nested rules GP, CC, F2, or GK).  Rather,
      // collocKey[i] is to be used only as the key to the collocation pts.
      ssg_driver->level_to_order(sm_multi_index[i], quad_order);
      quadrature_order_to_integrand_order(quad_order, integrand_order);
      integrand_order_to_expansion_order(integrand_order, expansion_order);
      tensor_product_multi_index(expansion_order, tpMultiIndex[i], true);
      append_multi_index(tpMultiIndex[i], multi_index, true);
#ifdef DEBUG
      PCout << "level =\n" << sm_multi_index[i] << "quad_order =\n"
	    << quad_order << "integrand_order =\n" << integrand_order
	    << "expansion_order =\n" << expansion_order << "tp_multi_index =\n"
	    << tpMultiIndex[i] << "multi_index =\n" << multi_index << '\n';
#endif // DEBUG
    }
    break;
  }
  case SPARSE_INT_TENSOR_SUM_EXP: {
    // assemble a complete list of individual polynomial coverage
    // defined from the linear combination of mixed tensor products
    multi_index.clear();
    UShort2DArray tp_multi_index;
    UShortArray integrand_order(numVars), expansion_order(numVars);
    for (i=0; i<num_smolyak_indices; ++i) {
      // Mitigate uneven integrand coverage due to exponential growth by
      // enforcing moderate linear growth (TO DO: SLOW_RESTRICTED_GROWTH):
      level_growth_to_integrand_order(sm_multi_index[i],
        MODERATE_RESTRICTED_GROWTH, integrand_order);
      integrand_order_to_expansion_order(integrand_order, expansion_order);
      tensor_product_multi_index(expansion_order, tp_multi_index, true);
      append_multi_index(tp_multi_index, multi_index, false);
#ifdef DEBUG
      PCout << "level =\n" << sm_multi_index[i] << "integrand_order =\n"
	    << integrand_order << "expansion_order =\n" << expansion_order
	    << "tp_multi_index =\n" << tp_multi_index << "multi_index =\n"
	    << multi_index << '\n';
#endif // DEBUG
    }
    break;
  }
  case SPARSE_INT_TOTAL_ORD_EXP: {
    // back out approxOrder & use total_order_multi_index()
    UShortArray quad_order(numVars), integrand_order(numVars);
    UShort2DArray pareto(1), total_pareto;
    for (i=0; i<num_smolyak_indices; ++i) {
      ssg_driver->level_to_order(sm_multi_index[i], quad_order);
      quadrature_order_to_integrand_order(quad_order, integrand_order);
      // maintain an n-dimensional Pareto front of nondominated multi-indices
      pareto[0] = integrand_order;
      update_pareto(pareto, total_pareto);
#ifdef DEBUG
      PCout << "level =\n" << sm_multi_index[i] << "\nquad_order =\n"
	    << quad_order << "\nintegrand_order =\n" << integrand_order << '\n';
#endif // DEBUG
    }
#ifdef DEBUG
    PCout << "total_pareto =\n" << total_pareto << '\n';
#endif // DEBUG

    // first pass: compute max isotropic integrand that fits within Pareto front
    unsigned short order = 0;
    integrand_order.assign(numVars, order);
    bool total_order_dominated = true;
    while (total_order_dominated) {
      // calculate all nondominated polynomials for a total-order expansion
      pareto.clear();
      total_order_multi_index(integrand_order, pareto, 0);
      total_order_dominated = assess_dominance(pareto, total_pareto);
#ifdef DEBUG
      PCout << "integrand_order =\n" << integrand_order << "pareto =\n"
	    << pareto << "total_order_dominated = " << total_order_dominated
	    << '\n';
#endif // DEBUG
      // could increment/decrement by 2's due to expansion_order conversion,
      // but the actual resolvable integrand order is typically odd.
      if (total_order_dominated)
	++order; // advance to next test level
      else
	--order; // exiting loop: rewind to last successful
      integrand_order.assign(numVars, order);
    }
#ifdef DEBUG
    PCout << "Isotropic integrand_order =\n" << integrand_order << '\n';
#endif // DEBUG
    /* Since anisotropic total-order expansions only prune corners of the Pascal
       hyper-pyramid and polynomial gaps from Smolyak (with nonlinear growth
       rules) are on the interior, the following approach is not effective.

    // second pass: if integrand order is anisotropic, assess dimensional
    // increments that may fit within Pareto front.  An important current
    // use case is PCBDO in all_variables mode.
    UShortArray sgl_to_int_order(n); // sg level -> quad order -> int order
    ssg_driver->level_to_order(ssg_driver->level(), quad_order);
    quadrature_order_to_integrand_order(quad_order, sgl_to_int_order);
    bool isotropic_integrand = true;
    order = sgl_to_int_order[0];
    for (i=1; i<n; ++i)
      if (sgl_to_int_order[i] != order)
	{ isotropic_integrand = false; break; }
    if (!isotropic_integrand) {
      bool remaining_dimensions = true;
      BoolDeque dominated(n, true);
      UShort2DArray to_multi_index;
      while (remaining_dimensions) {
	remaining_dimensions = false;
	for (i=0; i<n; ++i) {
	  if (dominated[i]) {
	    remaining_dimensions = true;
	    ++integrand_order[i];   // advance to next test level
	    // anisotropic total_order doesn't support lower bound offset,
	    // so must resort to defining the Pareto set in two steps
	    to_multi_index.clear();
	    total_order_multi_index(integrand_order, to_multi_index);//, 0);
	    pareto.clear();
	    update_pareto(to_multi_index, pareto);
	    dominated[i] = assess_dominance(pareto, total_pareto);
#ifdef DEBUG
	    PCout << "Anisotropic integrand_order =\n" << integrand_order
	        //<< "pareto =\n" << pareto
		  << "dominated[" << i << "] = " << dominated[i] << '\n';
#endif // DEBUG
	    if (!dominated[i])
	      --integrand_order[i]; // rewind to last successful
	  }
	}
      }
    }
    */

    integrand_order_to_expansion_order(integrand_order, approxOrder);
    total_order_multi_index(approxOrder, multi_index);
    break;
  }
  case SPARSE_INT_HEUR_TOTAL_ORD_EXP: // early heuristic
    sparse_grid_level_to_expansion_order(ssg_driver->level(), approxOrder);
    total_order_multi_index(approxOrder, multi_index);
    break;
  default: // possible fallback: revert to [optional] expansion_order spec.
    PCerr << "Error: unsupported sparseGridExpansion option in "
	  << "OrthogPolyApproximation::sparse_grid_multi_index()" << std::endl;
    abort_handler(-1);
    break;
  }
}


/* This approach reduces memory requirements but must perform additional
   calculation to regenerate the tp_multi_index instances (previously
   generated in sparse_grid_multi_index()).  Currently, these tp_multi_index
   instances are stored in tpMultiIndex for later use in compute_coefficients().
void OrthogPolyApproximation::
map_tensor_product_multi_index(UShort2DArray& tp_multi_index, size_t tp_index)
{
  const SizetArray& tp_mi_map = tpMultiIndexMap[tp_index];
  size_t i, num_tp_terms = tp_mi_map.size();
  tp_multi_index.resize(num_tp_terms);
  for (i=0; i<num_tp_terms; ++i)
    tp_multi_index[i] = multiIndex[tp_mi_map[i]];
}
*/


/* These mappings reflect simplified heuristics (expected to be exact
   for CC, Gauss-Patterson, and linear growth Gaussian, but not for
   exponential Gaussian). */
void OrthogPolyApproximation::
sparse_grid_level_to_expansion_order(unsigned short ssg_level,
				     UShortArray& exp_order)
{
  if (exp_order.size() != numVars)
    exp_order.resize(numVars);
  const IntArray& int_rules
    = ((SparseGridDriver*)driverRep)->integration_rules();
  const IntArray& growth_rules
    = ((SparseGridDriver*)driverRep)->growth_rules();
  for (size_t i=0; i<numVars; ++i) {
    switch (growth_rules[i]) {
    case FULL_EXPONENTIAL:
      switch (int_rules[i]) {
      case CLENSHAW_CURTIS: case FEJER2: { // TO DO: verify F2
	// integrand order = 2*level+1 ("sharp" result from Novak & Ritter,1996)
	unsigned short integrand = 2*ssg_level + 1;
	exp_order[i] = integrand / 2; // remainder truncated
	// results in exp_order = level
	break;
      }
      default:
	PCerr << "Error: unsupported integration rule for "
	      << "FULL_EXPONENTIAL growth in OrthogPolyApproximation::"
	      << "sparse_grid_level_to_expansion_order()." << std::endl;
	abort_handler(-1);
      }
      break;
    case MODERATE_LINEAR: {
      // linear growth Gauss rules & matching moderate growth exponential rules

      // For std exponential growth: use total integrand order = m (based on
      // experimental observations in PCE_logratio.m for n=2 and levels<5).
      // This heuristic has been shown to be nonconservative for higher
      // dimensions (e.g., short column w/ n=3 and cantilever beam w/ n=4).
      //SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
      //ssg_driver->level_to_order_open_exponential(ssg_level, integrand);

      // For linear growth (quad_order = 2*level+1):
      // This relationship derived from Pareto PCE eo for n={2,4} and w={0,8}
      // and verified for n=5 and n=9.
      // > Log ratio  n = 2: For w={0,8}, eo = {0,1,3,5,7,9,11,13,15}
      // > Short col  n = 3: For w={0,8}, eo = {0,1,2,4,6,8,10,12,14}
      // > Cantilever n = 4: For w={0,8}, eo = {0,1,2,3,5,7, 9,11,13}
      // > Gen Rosen  n = 5: For w={0,6}, eo = {0,1,2,3,4,6, 8}
      // > Steel col  n = 9: For w={0,5}, eo = {0,1,2,3,4,5}
      // -> Increases by 1 up to wmNp1=0 and increases by 2 thereafter.
      short wmNp1 = ssg_level - numVars + 1;
      exp_order[i] = (wmNp1 > 0) ? wmNp1 + ssg_level : ssg_level;
      break;
    }
    default:
      PCerr << "Error: unsupported growth rule in OrthogPolyApproximation::"
	    << "sparse_grid_level_to_expansion_order()," << std::endl;
      abort_handler(-1);
    }
  }
}


/** The computed integrand_order may be conservative (by design) in the
    presence of exponential growth.  This avoids aggressive optimization
    of PCE expansion orders when an exponential rule takes a large jump
    that is not balanced by the other index set component mappings. */
void OrthogPolyApproximation::
level_growth_to_integrand_order(const UShortArray& levels, short growth_rate,
				UShortArray& int_order)
{
  size_t i, n = levels.size();
  if (int_order.size() != n)
    int_order.resize(n);
  switch (growth_rate) {
  case SLOW_RESTRICTED_GROWTH:
    for (i=0; i<n; ++i) // synch with slow linear growth: i = 2l + 1
      int_order[i] =  2*levels[i] + 1;
    break;
  case MODERATE_RESTRICTED_GROWTH:
    for (i=0; i<n; ++i) // synch with moderate linear growth: i = 4l + 1
      int_order[i] =  4*levels[i] + 1;
    break;
  case UNRESTRICTED_GROWTH: {
    // This was previously the only option for SPARSE_INT_TENSOR_SUM_EXP
    // expansions, but was too aggressive for use with exponential growth.
    // It is retained for completeness, but is not recommended (this member
    // fn generalizes, but still provides access to the previous behavior).
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    UShortArray quad_order(numVars);
    ssg_driver->level_to_order(levels, quad_order);
    quadrature_order_to_integrand_order(quad_order, int_order);
    break;
  }
  }
}


void OrthogPolyApproximation::
quadrature_order_to_integrand_order(const UShortArray& quad_order,
				    UShortArray& int_order)
{
  // Need to know exact polynomial resolution for each mixed tensor grid:
  //   Gaussian integrand resolution:        2m-1
  //   Gauss-Patterson integrand resolution: 2m-1 - previous constraints + 1
  //   Clenshaw-Curtis integrand resolution: m (odd m), m-1 (even m)

  // Burkardt monomial test logic:
  //   sparse_grid_monomial_test: resolve monomials of total degree 2*level + 1
  //   for all rules --> doesn't make sense for exponential growth rules where
  //   order grows faster for Gauss than CC (level_to_order exponential is
  //   2^{w+1}-1 for Gauss and 2^w+1 for CC) --> estimate appears valid for CC
  //   (although it does not define a crisp boundary, since some monomials above
  //   the boundary are resolved) but overly conservative for Gauss (whole
  //   orders above the boundary estimate are resolved).

  size_t i, n = quad_order.size();
  if (int_order.size() != n)
    int_order.resize(n);
  if (gaussModes.empty()) // use basisTypes with default modes
    for (i=0; i<n; ++i)
      switch (basisTypes[i]) {
      case CHEBYSHEV: // default mode is Clenshaw-Curtis
	int_order[i] = (quad_order[i] % 2) ? quad_order[i] : quad_order[i] - 1;
	break;
      default: // default mode is standard non-nested Gauss rules
	int_order[i] =  2*quad_order[i] - 1; // i = 2m - 1
	break;
      }
  else {
    const UShortArray& gk_order = driverRep->genz_keister_order();
    const UShortArray& gk_prec  = driverRep->genz_keister_precision();
    for (i=0; i<n; ++i)
      switch (gaussModes[i]) {
      case CLENSHAW_CURTIS: case FEJER2:
	// i = m (odd m), m-1 (even m).  Note that growth rule enforces odd.
	// TO DO: verify FEJER2 same as CC
	int_order[i] = (quad_order[i] % 2) ? quad_order[i] : quad_order[i] - 1;
	break;
      case GAUSS_PATTERSON: {
	// for o(l)=2^{l+1}-1, o(l-1) = (o(l)-1)/2
	unsigned short prev_o = std::max(1,(quad_order[i] - 1)/2);
	int_order[i] = 2*quad_order[i] - prev_o;
	break;
      }
      case GENZ_KEISTER: {
	// same relationship as Gauss-Patterson, except prev_o does not follow
	// simple pattern and requires lookup
	unsigned short lev = 0, max_lev = 4;
	for (lev=0; lev<=max_lev; ++lev)
	  if (gk_order[lev] == quad_order[i])
	    { int_order[i] = gk_prec[lev]; break; }
	/*
	int lev, o, prev_o = 1, max_lev = 4, i_rule = GENZ_KEISTER,
	  g_rule = FULL_EXPONENTIAL; // map l->o directly without restriction
	for (lev=0; lev<=max_lev; ++lev) {
	  webbur::level_growth_to_order(1, &lev, &i_rule, &g_rule, &o);
	  if (o == quad_order[i])
	    { int_order[i] = 2*quad_order[i] - prev_o; break; }
	  else
	    prev_o = o;
	}
	*/
	if (lev>max_lev) {
	  PCerr << "Error: maximum GENZ_KEISTER level exceeded in OrthogPoly"
		<< "Approximation::quadrature_order_to_integrand_order()"
		<< std::endl;
	  abort_handler(-1);
	}
	break;
      }
      default: // standard non-nested Gauss rules
	int_order[i] =  2*quad_order[i] - 1; break; // i = 2m - 1
      }
  }
}


void OrthogPolyApproximation::
integrand_order_to_expansion_order(const UShortArray& int_order,
				   UShortArray& exp_order)
{
  // reserve half of the integrand order for the expansion and half for the
  // response function (integrand = 2p)
  size_t i, n = int_order.size();
  if (exp_order.size() != n)
    exp_order.resize(n);
  for (i=0; i<n; ++i)
    exp_order[i] = int_order[i] / 2; // remainder truncated
}


void OrthogPolyApproximation::
append_multi_index(const UShort2DArray& tp_multi_index,
		   UShort2DArray& multi_index, bool define_tp_mi_map)
{
  size_t i, num_tp_mi = tp_multi_index.size();
  if (multi_index.empty()) {
    multi_index = tp_multi_index;
    if (define_tp_mi_map) {
      tpMultiIndexMapRef.push_back(0);
      tpMultiIndexMap.resize(1); tpMultiIndexMap[0].resize(num_tp_mi);
      for (i=0; i<num_tp_mi; ++i)
	tpMultiIndexMap[0][i] = i;
    }
  }
  else if (define_tp_mi_map) {
    tpMultiIndexMapRef.push_back(multi_index.size());
    SizetArray tp_mi_map(num_tp_mi);
    for (i=0; i<num_tp_mi; ++i) {
      const UShortArray& search_mi = tp_multi_index[i];
      // TO DO: make this process more efficient
      size_t index = find_index(multi_index, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	tp_mi_map[i] = multi_index.size();
	multi_index.push_back(search_mi);
      }
      else
	tp_mi_map[i] = index;
    }
    tpMultiIndexMap.push_back(tp_mi_map);
  }
  else
    for (i=0; i<num_tp_mi; ++i) {
      const UShortArray& search_mi = tp_multi_index[i];
      // TO DO: make this process more efficient
      if (std::find(multi_index.begin(), multi_index.end(), search_mi) ==
	  multi_index.end()) // search_mi does not yet exist in multi_index
	multi_index.push_back(search_mi);
    }
}


void OrthogPolyApproximation::
append_multi_index(const UShort2DArray& tp_multi_index,
		   SizetArray& tp_mi_map, size_t& tp_mi_map_ref,
		   UShort2DArray& multi_index)
{
  if (multi_index.empty())
    multi_index = tp_multi_index;
  else {
    size_t i, num_tp_mi = tp_multi_index.size(), num_mi = multi_index.size();
    if (num_mi == tp_mi_map_ref) { // current multi_index corresponds to ref
      for (i=0; i<num_tp_mi; ++i)
	if (tp_mi_map[i] >= tp_mi_map_ref)
	  multi_index.push_back(tp_multi_index[i]);
    }
    else if (num_mi > tp_mi_map_ref) { // multi_index has grown since ref taken
      for (i=0; i<num_tp_mi; ++i)
	if (tp_mi_map[i] >= tp_mi_map_ref) { // was previously appended
	  const UShortArray& search_mi = tp_multi_index[i];
	  // search from reference pt forward
	  UShort2DArray::iterator it, it_start = multi_index.begin();
	  std::advance(it_start, tp_mi_map_ref);
	  it = std::find(it_start, multi_index.end(), search_mi);
	  if (it == multi_index.end()) { // still an append: update map, append
	    tp_mi_map[i] = multi_index.size();
	    multi_index.push_back(tp_multi_index[i]);
	  }
	  else // no longer an append: only update map
	    tp_mi_map[i] = tp_mi_map_ref + std::distance(it_start, it);
	}
      tp_mi_map_ref = num_mi; // reference point now updated
    }
    else { // multi_index is not allowed to shrink since ref taken
      PCerr << "Error: multi_index inconsistent with reference size in "
	    << "OrthogPolyApproximation::append_multi_index()." << std::endl;
      abort_handler(-1);
    }
  }
}


void OrthogPolyApproximation::
combine_expansion(const SizetArray& tp_mi_map,
		  const RealVector& tp_expansion_coeffs,
		  const RealMatrix& tp_expansion_grads, int coeff)
{
  size_t i, j, index, num_tp_terms = tp_mi_map.size(), 
    num_deriv_vars = expansionCoeffGrads.numRows();
  for (i=0; i<num_tp_terms; ++i) {
    index = tp_mi_map[i];
    if (configOptions.expansionCoeffFlag)
      expansionCoeffs[index] += coeff * tp_expansion_coeffs[i];
    if (configOptions.expansionCoeffGradFlag) {
      Real*       exp_grad_ndx = expansionCoeffGrads[index];
      const Real* tp_grad_i    = tp_expansion_grads[i];
      for (j=0; j<num_deriv_vars; ++j)
	exp_grad_ndx[j] += coeff * tp_grad_i[j];
    }
  }
}


void OrthogPolyApproximation::append_expansions(size_t start_index)
{
  // for use in decrement_coefficients()
  prevExpCoeffs = expansionCoeffs; prevExpCoeffGrads = expansionCoeffGrads;

  // update expansion{Coeffs,CoeffGrads} using a hierarchical update
  // rather than building from scratch
  SparseGridDriver*  ssg_driver = (SparseGridDriver*)driverRep;
  const IntArray&     sm_coeffs = ssg_driver->smolyak_coefficients();
  const IntArray& ref_sm_coeffs = ssg_driver->reference_smolyak_coefficients();
#ifdef DEBUG
  PCout << "In OrthogPolyApproximation::append_expansions() with start index "
	<< start_index << "\nsm_coeffs:\n" << sm_coeffs << "ref_sm_coeffs:\n"
	<< ref_sm_coeffs << std::endl;
#endif // DEBUG

  // add trial expansions
  size_t index, num_tensor_grids = sm_coeffs.size();
  for (index=start_index; index<num_tensor_grids; ++index) {
    combine_expansion(tpMultiIndexMap[index], tpExpansionCoeffs[index],
		      tpExpansionCoeffGrads[index], sm_coeffs[index]);
#ifdef DEBUG
    PCout << "Trial set sm_coeff = " << sm_coeffs[index]
	  << "\ntpExpansionCoeffs:\n";
    write_data(PCout, tpExpansionCoeffs[index]);
    PCout << "\ntpMultiIndexMap:\n" << tpMultiIndexMap[index] << '\n';
#endif // DEBUG
  }
  // update other expansion contributions with a changed smolyak coefficient
  int delta_coeff;
  for (index=0; index<start_index; ++index) {
    // add new, subtract previous
    delta_coeff = sm_coeffs[index] - ref_sm_coeffs[index];
#ifdef DEBUG
    PCout << "Old set delta_coeff = " << delta_coeff
	  << "\ntpExpansionCoeffs:\n";
    write_data(PCout, tpExpansionCoeffs[index]);
    PCout << "\ntpMultiIndexMap:\n" << tpMultiIndexMap[index] << '\n';
#endif // DEBUG
    if (delta_coeff)
      combine_expansion(tpMultiIndexMap[index], tpExpansionCoeffs[index],
			tpExpansionCoeffGrads[index], delta_coeff);
  }
}


void OrthogPolyApproximation::
update_pareto(const UShort2DArray& new_pareto, UShort2DArray& total_pareto)
{
  //if (total_pareto.empty())
  //  total_pareto = new_pareto;
  //else {
    size_t i, num_new_p = new_pareto.size(),
      num_total_p = total_pareto.size();
    std::list<UShort2DArray::iterator> rm_iters; // sorted, unique
    UShort2DArray::iterator jit;
    for (i=0; i<num_new_p; ++i) {
      const UShortArray& new_i = new_pareto[i];
      bool new_i_dominated = false;
      for (jit=total_pareto.begin(); jit!=total_pareto.end(); ++jit) {
	bool new_i_dominated_by_j, total_j_dominated;
	assess_dominance(new_i, *jit, new_i_dominated_by_j, total_j_dominated);
	if (new_i_dominated_by_j)
	  new_i_dominated = true;
	if (total_j_dominated)
	  rm_iters.push_back(jit);
      }
      // 
      // prune newly dominated in reverse order (vector iterators
      // following a point of insertion or deletion are invalidated)
      while (!rm_iters.empty())
	{ total_pareto.erase(rm_iters.back()); rm_iters.pop_back(); }
      // add nondominated
      if (!new_i_dominated)
	total_pareto.push_back(new_i);
    }
  //}
}


bool OrthogPolyApproximation::
assess_dominance(const UShort2DArray& new_pareto,
		 const UShort2DArray& total_pareto)
{
  bool new_dominated = true;
  size_t i, j, num_new_p = new_pareto.size(),
    num_total_p = total_pareto.size();
  for (i=0; i<num_new_p; ++i) {
    const UShortArray& new_i = new_pareto[i];
    bool new_i_dominated = false;
    for (j=0; j<num_total_p; ++j) {
      bool i_dominated_by_j, j_dominated;
      assess_dominance(new_i, total_pareto[j], i_dominated_by_j, j_dominated);
      if (i_dominated_by_j)
	{ new_i_dominated = true; break; }
    }
    if (!new_i_dominated) {
      new_dominated = false;
#ifdef DEBUG
      PCout << "Nondominated new pareto member =\n" << new_i;
#else
      break;
#endif // DEBUG
    }
  }
  return new_dominated;
}


void OrthogPolyApproximation::
assess_dominance(const UShortArray& new_order,
		 const UShortArray& existing_order,
		 bool& new_dominated, bool& existing_dominated)
{
  // can't use std::vector::operator< (used for component-wise sorting)
  size_t i, n = new_order.size();
  bool equal = true, existing_dominated_temp = true;
  new_dominated = true;
  for (i=0; i<n; ++i)
    if (new_order[i] > existing_order[i])
      { equal = false; new_dominated = false; }
    else if (existing_order[i] > new_order[i])
      { equal = false; existing_dominated_temp = false; }
  existing_dominated = (!equal && existing_dominated_temp);
}


void OrthogPolyApproximation::integration_checks()
{
  if (!anchorPoint.is_null()) { // TO DO: verify this creates a problem
    PCerr << "Error: anchor point not supported for numerical integration in "
	  << "OrthogPolyApproximation::integration()." << std::endl;
    abort_handler(-1);
  }
  if (!driverRep) {
    PCerr << "Error: pointer to integration driver required in "
	  << "OrthogPolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }
  size_t num_pts = dataPoints.size();
  if (num_pts != driverRep->weight_sets().length()) {
    PCerr << "Error: number of current points (" << num_pts << ") is not "
	  << "consistent with\n       number of points/weights from "
	  << "integration driver in\n       OrthogPolyApproximation::"
	  << "compute_coefficients()." << std::endl;
    abort_handler(-1);
  }
}


void OrthogPolyApproximation::
integration_data(size_t tp_index,
		 std::vector<SurrogateDataPoint>& tp_data_points,
		 RealVector& tp_weights)
{
  // extract tensor points from dataPoints and tensor weights from gaussWts1D
  SparseGridDriver*    ssg_driver = (SparseGridDriver*)driverRep;
  const UShortArray&     sm_index = ssg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&        key = ssg_driver->collocation_key()[tp_index];
  const SizetArray&  colloc_index = ssg_driver->collocation_indices()[tp_index];
  const Real3DArray& gauss_wts_1d = ssg_driver->gauss_weights_array();
  size_t i, j, num_tp_pts = colloc_index.size();
  tp_data_points.resize(num_tp_pts); tp_weights.resize(num_tp_pts);
  for (i=0; i<num_tp_pts; ++i) {
    // tensor-product point
    tp_data_points[i] = dataPoints[colloc_index[i]];
    // tensor-product weight
    Real& tp_wts_i = tp_weights[i]; tp_wts_i = 1.;
    const UShortArray& key_i = key[i];
    for (j=0; j<numVars; ++j)
      tp_wts_i *= gauss_wts_1d[sm_index[j]][j][key_i[j]];
  }
}


/** The coefficients of the PCE for the response are calculated using a
    spectral projection of the response against each multivariate orthogonal
    polynomial basis fn using the inner product ratio <f,Psi>/<Psi^2>, where
    inner product <a,b> is the n-dimensional integral of a*b*weighting over
    the support range of the n-dimensional (composite) weighting function.
    1-D quadrature rules are defined for specific 1-D weighting functions
    and support ranges and approximate the integral of f*weighting as the
    Sum_i of w_i f_i.  To extend this to n-dimensions, a tensor product
    quadrature rule, cubature, or Smolyak sparse grid rule is applied.  
    It is not necessary to approximate the integral for the denominator
    numerically, since this is available analytically. */
void OrthogPolyApproximation::
integrate_expansion(const UShort2DArray& multi_index,
		    const std::vector<SurrogateDataPoint>& data_pts,
		    const RealVector& wt_sets, RealVector& exp_coeffs,
		    RealMatrix& exp_coeff_grads)
{
  // Perform numerical integration via tensor-product quadrature/cubature/
  // Smolyak sparse grids.  Quadrature/cubature use a single application of
  // point and weight sets computed by TensorProductDriver/CubatureDriver, and
  // sparse grids could do this as well, but it is better to integrate the
  // sparse grid on a per-tensor-product basis folowed by summing the
  // corresponding PC expansions.
  std::vector<SurrogateDataPoint>::const_iterator cit;
  size_t i, j, k, num_exp_terms = multi_index.size(),
    num_deriv_vars = data_pts[0].response_gradient().length();
  Real wt_resp_fn_i, Psi_ij; Real* exp_grad;
  RealVector wt_resp_grad_i;
  if (configOptions.expansionCoeffFlag) { // shape if needed and zero out
    if (exp_coeffs.length() != num_exp_terms)
      exp_coeffs.size(num_exp_terms); // init to 0
    else
      exp_coeffs = 0.;
  }
  if (configOptions.expansionCoeffGradFlag) {
    if (exp_coeff_grads.numRows() != num_deriv_vars ||
	exp_coeff_grads.numCols() != num_exp_terms)
      exp_coeff_grads.shape(num_deriv_vars, num_exp_terms); // init to 0
    else
      exp_coeff_grads = 0.;
    wt_resp_grad_i.sizeUninitialized(num_deriv_vars);
  }
  for (i=0, cit=data_pts.begin(); cit!=data_pts.end(); ++i, ++cit) {
    if (configOptions.expansionCoeffFlag)
      wt_resp_fn_i = wt_sets[i] * cit->response_function();
    if (configOptions.expansionCoeffGradFlag) {
      wt_resp_grad_i = cit->response_gradient(); // copy
      wt_resp_grad_i.scale(wt_sets[i]);
    }
    const RealVector& c_vars_i = cit->continuous_variables();
    for (j=0; j<num_exp_terms; ++j) {
      Psi_ij = multivariate_polynomial(c_vars_i, multi_index[j]);
      if (configOptions.expansionCoeffFlag)
	exp_coeffs[j] += Psi_ij * wt_resp_fn_i;
      if (configOptions.expansionCoeffGradFlag) {
	exp_grad = exp_coeff_grads[j];
	for (k=0; k<num_deriv_vars; ++k)
	  exp_grad[k] += Psi_ij * wt_resp_grad_i[k];
      }
    }
  }

  for (i=0; i<num_exp_terms; ++i) {
    const Real& norm_sq = norm_squared(multi_index[i]);
    if (configOptions.expansionCoeffFlag)
      exp_coeffs[i] /= norm_sq;
    if (configOptions.expansionCoeffGradFlag) {
      exp_grad = exp_coeff_grads[i];
      for (k=0; k<num_deriv_vars; ++k)
	exp_grad[k] /= norm_sq;
    }
  }
#ifdef DEBUG
  PCout << "expansion_coeffs:\n" << exp_coeffs
	<< "expansion_coeff_grads:\n" << exp_coeff_grads<< "\n\n";
#endif // DEBUG
}


// Compute moments by projecting data points along the normalized identity function
void OrthogPolyApproximation::
integrate_moments(const UShort2DArray& multi_index,
		    const std::vector<SurrogateDataPoint>& data_pts,
                    RealVector& function_vals,
		    const RealVector& wt_sets, RealVector& moments)
		    
{
  std::vector<SurrogateDataPoint>::const_iterator cit;
  size_t i, j, k;
  size_t m = moments.length(); 
  RealVector wt_resp_fn_i(m); Real Psi_ij; 
  moments = 0;

  // construct vector of wts and modified data points
  for (i=0, cit=data_pts.begin(); cit!=data_pts.end(); ++i, ++cit) {
    for (j=0; j<m; ++j) {
      // do not center the first moment
      if (j == 0)
        wt_resp_fn_i[j] = wt_sets[i] * cit->response_function(); 
      else
        wt_resp_fn_i[j] = wt_sets[i] * std::pow(function_vals[i],Real(j+1));
    }
    const RealVector& c_vars_i = cit->continuous_variables();
    // only need to access the first basis function (i.e. identity function)
    Psi_ij = multivariate_polynomial(c_vars_i, multi_index[0]);
    for (j=0; j<m; ++j)
      moments[j] += Psi_ij * wt_resp_fn_i[j];
  }

}

/** In this case, regression is used in place of spectral projection.  That
    is, instead of calculating the PCE coefficients using inner product
    ratios, linear least squares is used to estimate the PCE coefficients
    which best match a set of response samples.  This approach is also known
    as stochastic response surfaces.  The least squares estimation is
    performed using DGELSS (SVD) or DGGLSE (equality-constrained) from
    LAPACK, based on anchorPoint and derivative data availability. */
void OrthogPolyApproximation::regression()
{
  size_t i, j, num_pts = dataPoints.size();
  bool err_flag = false, anchor_pt = !anchorPoint.is_null();

  // compute order of data contained within anchorPoint/dataPoints
  SurrogateDataPoint& pt0 = (anchor_pt) ? anchorPoint : dataPoints[0];
  unsigned short data_order = 1;
  if (!pt0.response_gradient().empty()) data_order |= 2;
  if (!pt0.response_hessian().empty())  data_order |= 4;

  // use_grads_flag indicates usage of gradient data with respect to the
  // expansion variables (aleatory or combined) within the expansion
  // coefficient calculation process, which must be distinguished from the
  // expansionCoeffGradFlag case.
  bool use_grads_flag
    = ((data_order & 2) && !configOptions.expansionCoeffGradFlag),
    fn_constrained_lls = (use_grads_flag && num_pts < numExpansionTerms);
  size_t eqns_per_pt = (use_grads_flag) ? 1 + numVars : 1;
  std::vector<SurrogateDataPoint>::iterator it;
  Teuchos::LAPACK<int, Real> la;
  double *A_matrix, *work;
  int info       = 0, // output flag from DGELSS/DGGLSE subroutines
      num_cols_A = numExpansionTerms; // Number of columns in matrix A

  if (anchor_pt) {
    // Use DGGLSE for equality-constrained LLS soln using GRQ factorization.
    // Solves min ||b - Ax||_2 s.t. Cx = d (Note: b,C switched from LAPACK docs)
    // where {b,d} are single vectors (multiple RHS not supported).

    int num_cons     = eqns_per_pt; // constraints from one anchor point
    int num_rows_A   = num_pts*eqns_per_pt; // Number of rows in matrix A
    // Matrix of polynomial terms in a contiguous block of memory packed in
    // column-major ordering as required by F77 LAPACK subroutines.
    A_matrix         = new double [num_rows_A*num_cols_A]; // "A" in A*x = b
    // Vector of response values that correspond to the samples in matrix A.
    double* b_vector = new double [num_rows_A]; // "b" in A*x = b
    // Matrix of constraints unrolled into a vector
    double* C_matrix = new double [num_cols_A]; // "C" in C*x = d
    // RHS of constraints
    double* d_vector = new double [num_cons];   // "d" in C*x = d
    // Solution vector
    double* x_vector = new double [num_cols_A]; // "x" in b - A*x, C*x = d

    // Get the optimal work array size
    int lwork        = -1; // special code for workspace query
    work             = new double [1]; // temporary work array
    la.GGLSE(num_rows_A, num_cols_A, num_cons, A_matrix, num_rows_A, C_matrix,
	     num_cons, b_vector, d_vector, x_vector, work, lwork, &info);
    lwork = (int)work[0]; // optimal work array size returned by query
    delete [] work;
    work             = new double [lwork]; // Optimal work array

    if (configOptions.expansionCoeffFlag) {
      // pack dataPoints values/gradients in A,b and anchor point
      // values/gradients in C,d.
      size_t a_cntr = 0, b_cntr = 0, c_cntr = 0, d_cntr = 0;
      for (i=0; i<numExpansionTerms; ++i) {
	const UShortArray& mi = multiIndex[i];
	for (it=dataPoints.begin(); it!=dataPoints.end(); ++it) {
	  const RealVector& c_vars = it->continuous_variables();
	  A_matrix[a_cntr] = multivariate_polynomial(c_vars, mi); ++a_cntr;
	  if (use_grads_flag) {
	    const RealVector& mvp_grad
	      = multivariate_polynomial_gradient(c_vars, mi);
	    for (j=0; j<numVars; ++j, ++a_cntr)
	      A_matrix[a_cntr] = mvp_grad[j];
	  }
	}
	const RealVector& ap_c_vars = anchorPoint.continuous_variables();
	C_matrix[c_cntr] = multivariate_polynomial(ap_c_vars, mi); ++c_cntr;
	if (use_grads_flag) {
	  const RealVector& mvp_grad
	    = multivariate_polynomial_gradient(ap_c_vars, mi);
	  for (j=0; j<numVars; ++j, ++c_cntr)
	    C_matrix[c_cntr] = mvp_grad[j];
	}
      }
      for (it=dataPoints.begin(); it!=dataPoints.end(); ++it) {
	b_vector[b_cntr] = it->response_function(); ++b_cntr;
	if (use_grads_flag) {
	  const RealVector& resp_grad = it->response_gradient();
	  for (j=0; j<numVars; ++j, ++b_cntr)
	    b_vector[b_cntr] = resp_grad[j];
	}
      }
      d_vector[d_cntr] = anchorPoint.response_function(); ++d_cntr;
      if (use_grads_flag) {
	const RealVector& resp_grad = anchorPoint.response_gradient();
	for (j=0; j<numVars; ++j, ++d_cntr)
	  d_vector[d_cntr] = resp_grad[j];
      }
      // Least squares computation using LAPACK's DGGLSE subroutine which uses a
      // GRQ factorization method for solving the least squares problem
      info = 0;
      la.GGLSE(num_rows_A, num_cols_A, num_cons, A_matrix, num_rows_A,
	       C_matrix, num_cons, b_vector, d_vector, x_vector, work,
	       lwork, &info);
      if (info)
	err_flag = true;
      copy_data(x_vector, numExpansionTerms, expansionCoeffs);
    }

    if (configOptions.expansionCoeffGradFlag) {
      size_t num_deriv_vars = expansionCoeffGrads.numRows();
      for (i=0; i<num_deriv_vars; ++i) {
	// must be recomputed each time since DGGLSE solves in place
	size_t a_cntr = 0;
	for (j=0; j<numExpansionTerms; ++j) {
	  const UShortArray& mi = multiIndex[j];
	  for (it=dataPoints.begin(); it!=dataPoints.end(); ++it, ++a_cntr)
	    A_matrix[a_cntr]
	      = multivariate_polynomial(it->continuous_variables(), mi);
	  C_matrix[j]
	    = multivariate_polynomial(anchorPoint.continuous_variables(), mi);
	}
	// the Ax=b RHS is the dataPoints values for the i-th grad component
	for (j=0, it=dataPoints.begin(); j<num_rows_A; ++j, ++it)
	  b_vector[j] = it->response_gradient()[i];
	// the Cx=d RHS is the anchorPoint value for the i-th gradient component
	d_vector[0] = anchorPoint.response_gradient()[i];
	// solve for the PCE coefficients for the i-th gradient component
	info = 0;
	la.GGLSE(num_rows_A, num_cols_A, num_cons, A_matrix, num_rows_A,
		 C_matrix, num_cons, b_vector, d_vector, x_vector, work,
		 lwork, &info);
	if (info)
	  err_flag = true;
	for (j=0; j<numExpansionTerms; ++j)
	  expansionCoeffGrads(i,j) = x_vector[j];
      }
    }

    delete [] b_vector;
    delete [] C_matrix;
    delete [] d_vector;
    delete [] x_vector;
  }
  else if (fn_constrained_lls) {

    int num_rows_A = num_pts*numVars; // Number of rows in matrix A
    int num_cons   = num_pts; // constraints from response function values
    // Matrix of polynomial terms in a contiguous block of memory packed in
    // column-major ordering as required by F77 LAPACK subroutines.
    A_matrix         = new double [num_rows_A*num_cols_A]; // "A" in A*x = b
    // Vector of response values that correspond to the samples in matrix A.
    double* b_vector = new double [num_rows_A]; // "b" in A*x = b
    // Matrix of constraints unrolled into a vector
    double* C_matrix = new double [num_cons*num_cols_A]; // "C" in C*x = d
    // RHS of constraints
    double* d_vector = new double [num_cons];   // "d" in C*x = d
    // Solution vector
    double* x_vector = new double [num_cols_A]; // "x" in b - A*x, C*x = d

    // Get the optimal work array size
    int lwork        = -1; // special code for workspace query
    work             = new double [1]; // temporary work array
    la.GGLSE(num_rows_A, num_cols_A, num_cons, A_matrix, num_rows_A, C_matrix,
	     num_cons, b_vector, d_vector, x_vector, work, lwork, &info);
    lwork = (int)work[0]; // optimal work array size returned by query
    delete [] work;
    work             = new double [lwork]; // Optimal work array

    // pack dataPoints gradients in A,b and dataPoints values in C,d.
    size_t a_cntr = 0, b_cntr = 0, c_cntr = 0, d_cntr = 0;
    for (i=0; i<numExpansionTerms; ++i) {
      const UShortArray& mi = multiIndex[i];
      for (it=dataPoints.begin(); it!=dataPoints.end(); ++it, ++c_cntr) {
	const RealVector& c_vars = it->continuous_variables();
	// hard constraint on response values
	C_matrix[c_cntr] = multivariate_polynomial(c_vars, mi);
	// LLS on response gradients
	const RealVector& mvp_grad
	  = multivariate_polynomial_gradient(c_vars, mi);
	for (j=0; j<numVars; ++j, ++a_cntr)
	  A_matrix[a_cntr] = mvp_grad[j];
      }
    }
    for (it=dataPoints.begin(); it!=dataPoints.end(); ++it, ++d_cntr) {
      d_vector[d_cntr] = it->response_function();
      const RealVector& resp_grad = it->response_gradient();
      for (j=0; j<numVars; ++j, ++b_cntr)
	b_vector[b_cntr] = resp_grad[j];
    }
    // Least squares computation using LAPACK's DGGLSE subroutine which uses a
    // GRQ factorization method for solving the least squares problem
    info = 0;
    la.GGLSE(num_rows_A, num_cols_A, num_cons, A_matrix, num_rows_A,
	     C_matrix, num_cons, b_vector, d_vector, x_vector, work,
	     lwork, &info);
    if (info)
      err_flag = true;
    copy_data(x_vector, numExpansionTerms, expansionCoeffs);

    delete [] b_vector;
    delete [] C_matrix;
    delete [] d_vector;
    delete [] x_vector;
  }
  else {
    // Use DGELSS for LLS soln using SVD.  Solves min ||b - Ax||_2
    // where {b} may have multiple RHS -> multiple {x} solutions.

    int    num_rows_A    = num_pts*eqns_per_pt; // Number of rows in matrix A
    size_t num_grad_rhs  = (configOptions.expansionCoeffGradFlag) ?
      expansionCoeffGrads.numRows() : 0;
    size_t num_coeff_rhs = (configOptions.expansionCoeffFlag) ? 1 : 0;
    int    num_rhs = num_coeff_rhs + num_grad_rhs; // input: # of RHS vectors
    double rcond = -1.; // input:  use macheps to rank singular vals of A
    int    rank  =  0;  // output: effective rank of matrix A
    // Matrix of polynomial terms unrolled into a vector.
    A_matrix          = new double [num_rows_A*num_cols_A]; // "A" in A*x = b
    // input: vectors of response data that correspond to samples in matrix A.
    double* b_vectors = new double [num_rows_A*num_rhs]; // "b" in A*x = b
    // output: vector of singular values, dimensioned for overdetermined system
    double* s_vector  = new double [num_cols_A];

    // Get the optimal work array size
    int lwork         = -1; // special code for workspace query
    work              = new double [1]; // temporary work array
    la.GELSS(num_rows_A, num_cols_A, num_rhs, A_matrix, num_rows_A, b_vectors,
	     num_rows_A, s_vector, rcond, &rank, work, lwork, &info);
    lwork = (int)work[0]; // optimal work array size returned by query
    delete [] work;
    work              = new double [lwork]; // Optimal work array

    // The "A" matrix is a contiguous block of memory packed in column-major
    // ordering as required by F77 for the GELSS subroutine from LAPACK.  For
    // example, the 6 elements of A(2,3) are stored in the order A(1,1), A(2,1),
    // A(1,2), A(2,2), A(1,3), A(2,3).
    size_t a_cntr = 0, b_cntr = 0;
    for (i=0; i<numExpansionTerms; ++i) {
      const UShortArray& mi = multiIndex[i];
      for (it=dataPoints.begin(); it!=dataPoints.end(); ++it) {
	const RealVector& c_vars = it->continuous_variables();
	A_matrix[a_cntr] = multivariate_polynomial(c_vars, mi); ++a_cntr;
	if (use_grads_flag) {
	  const RealVector& mvp_grad
	    = multivariate_polynomial_gradient(c_vars, mi);
	  for (j=0; j<numVars; ++j, ++a_cntr)
	    A_matrix[a_cntr] = mvp_grad[j];
	}
      }
    }

    // response data (values/gradients) define the multiple RHS which are
    // matched in the LS soln.  b_vectors is num_pts (rows) x num_rhs (cols),
    // arranged in column-major order.
    if (use_grads_flag)
      for (it=dataPoints.begin(); it!=dataPoints.end(); ++it) {
	b_vectors[b_cntr] = it->response_function(); ++b_cntr;
	const RealVector& resp_grad = it->response_gradient();
	for (j=0; j<numVars; ++j, ++b_cntr)
	  b_vectors[b_cntr] = resp_grad[j];
      }
    else
      for (i=0, it=dataPoints.begin(); i<num_pts; ++it, ++i) {
	if (configOptions.expansionCoeffFlag)
	  b_vectors[i] = it->response_function();
	if (configOptions.expansionCoeffGradFlag) {
	  const RealVector& resp_grad = it->response_gradient();
	  for (j=0; j<num_grad_rhs; ++j) // i-th point, j-th gradient component
	    b_vectors[(j+num_coeff_rhs)*num_pts+i] = resp_grad[j];
	}
      }

    // Least squares computation using LAPACK's DGELSS subroutine which uses a
    // SVD method for solving the least squares problem
    info = 0;
    la.GELSS(num_rows_A, num_cols_A, num_rhs, A_matrix, num_rows_A, b_vectors,
	     num_rows_A, s_vector, rcond, &rank, work, lwork, &info);
    if (info)
      err_flag = true;

    // DGELSS returns the x solution in the numExpansionTerms x num_rhs
    // submatrix of b_vectors
    if (configOptions.expansionCoeffFlag)
      copy_data(b_vectors, numExpansionTerms, expansionCoeffs);
    if (configOptions.expansionCoeffGradFlag)
      for (i=0; i<numExpansionTerms; ++i)
	for (j=0; j<num_grad_rhs; ++j)
	  expansionCoeffGrads(j,i) = b_vectors[(j+num_coeff_rhs)*num_pts+i];

    delete [] b_vectors;
    delete [] s_vector;
  }
  delete [] A_matrix;
  delete [] work;

  if (err_flag) { // if numerical problems in LLS, abort with error message
    PCerr << "Error: nonzero return code from LAPACK linear least squares in "
	  << "OrthogPolyApproximation::compute_coefficients()" << std::endl;
    abort_handler(-1);
  }
}


/** The coefficients of the PCE for the response are calculated using a
    spectral projection of the response against each multivariate orthogonal
    polynomial basis fn using the inner product ratio <f,Psi>/<Psi^2>,
    where inner product <a,b> is the n-dimensional integral of a*b*weighting
    over the support range of the n-dimensional (composite) weighting
    function.  When interpreting the weighting function as a probability
    density function, <a,b> = expected value of a*b, which can be evaluated
    by sampling from the probability density function and computing the mean
    statistic.  It is not necessary to compute the mean statistic for the
    denominator, since this is available analytically. */
void OrthogPolyApproximation::expectation()
{
  // "lhs" or "random", no weights needed
  std::vector<SurrogateDataPoint>::iterator it;
  size_t i, j, num_pts = dataPoints.size(),
    num_deriv_vars = expansionCoeffGrads.numRows();
  if (!anchorPoint.is_null())
    ++num_pts;

  /*
  // The following implementation evaluates all PCE coefficients
  // using a consistent expectation formulation
  for (i=0; i<numExpansionTerms; ++i) {
    Real& exp_coeff_i = expansionCoeffs[i];
    exp_coeff_i = (anchorPoint.is_null()) ? 0.0 :
      anchorPoint.response_function() * multivariate_polynomial(
        anchorPoint.continuous_variables(), multiIndex[i]);
    for (it=dataPoints.begin(); it!=dataPoints.end(); ++it)
      exp_coeff_i += it->response_function() * 
        multivariate_polynomial(it->continuous_variables(), multiIndex[i]);
    exp_coeff_i /= num_pts * norm_squared(multiIndex[i]);
#ifdef DEBUG
    PCout << "coeff[" << i << "] = " << exp_coeff_i
	  << " norm squared[" << i <<"] = " << norm_squared(multiIndex[i])
	  << '\n';
#endif // DEBUG
  }
  */

  // This alternate implementation evaluates the first PCE coefficient (the
  // response mean) as an expectation and then removes the mean from the
  // expectation evaluation of all subsequent coefficients.  This approach
  // has been observed to result in better results for small sample sizes.
  Real empty_r;
  Real& mean      = (configOptions.expansionCoeffFlag) ?
    expansionCoeffs[0] : empty_r;
  Real* mean_grad = (configOptions.expansionCoeffGradFlag) ?
    (Real*)expansionCoeffGrads[0] : NULL;
  if (anchorPoint.is_null()) {
    if (configOptions.expansionCoeffFlag)     expansionCoeffs     = 0.;
    if (configOptions.expansionCoeffGradFlag) expansionCoeffGrads = 0.;
  }
  else {
    if (configOptions.expansionCoeffFlag)
      mean = anchorPoint.response_function();
    if (configOptions.expansionCoeffGradFlag)
      mean_grad = anchorPoint.response_gradient().values();
  }
  for (it=dataPoints.begin(); it!=dataPoints.end(); ++it) {
    if (configOptions.expansionCoeffFlag)
      mean += it->response_function();
    if (configOptions.expansionCoeffGradFlag) {
      const RealVector& curr_pt_grad = it->response_gradient();
      for (j=0; j<num_deriv_vars; ++j)
	mean_grad[j] += curr_pt_grad[j];
    }
  }
  if (configOptions.expansionCoeffFlag)
    mean /= num_pts;
  if (configOptions.expansionCoeffGradFlag)
    for (j=0; j<num_deriv_vars; ++j)
      mean_grad[j] /= num_pts;
  Real chaos_sample, resp_fn_minus_mean, term; Real* exp_grad_i;
  RealVector resp_grad_minus_mean;
  if (configOptions.expansionCoeffGradFlag)
    resp_grad_minus_mean.sizeUninitialized(num_deriv_vars);
  if (!anchorPoint.is_null()) {
    if (configOptions.expansionCoeffFlag)
      resp_fn_minus_mean = anchorPoint.response_function() - mean;
    if (configOptions.expansionCoeffGradFlag) {
      const RealVector& anchor_grad = anchorPoint.response_gradient();
      for (j=0; j<num_deriv_vars; ++j)
	resp_grad_minus_mean[j] = anchor_grad[j] - mean_grad[j];
    }
    const RealVector& c_vars = anchorPoint.continuous_variables();
    for (i=1; i<numExpansionTerms; ++i) {
      chaos_sample = multivariate_polynomial(c_vars, multiIndex[i]);
      if (configOptions.expansionCoeffFlag)
	expansionCoeffs[i] = resp_fn_minus_mean * chaos_sample;
      if (configOptions.expansionCoeffGradFlag) {
	exp_grad_i = expansionCoeffGrads[i];
	for (j=0; j<num_deriv_vars; ++j)
	  exp_grad_i[j] = resp_grad_minus_mean[j] * chaos_sample;
      }
    }
  }
  for (it=dataPoints.begin(); it!=dataPoints.end(); ++it) {
    if (configOptions.expansionCoeffFlag)
      resp_fn_minus_mean = it->response_function() - mean;
    if (configOptions.expansionCoeffGradFlag) {
      const RealVector& resp_grad = it->response_gradient();
      for (j=0; j<num_deriv_vars; ++j)
	resp_grad_minus_mean[j] = resp_grad[j] - mean_grad[j];
    }
    const RealVector& c_vars = it->continuous_variables();
    for (i=1; i<numExpansionTerms; ++i) {
      chaos_sample = multivariate_polynomial(c_vars, multiIndex[i]);
      if (configOptions.expansionCoeffFlag)
	expansionCoeffs[i] += resp_fn_minus_mean * chaos_sample;
      if (configOptions.expansionCoeffGradFlag) {
	exp_grad_i = expansionCoeffGrads[i];
	for (j=0; j<num_deriv_vars; ++j)
	  exp_grad_i[j] += resp_grad_minus_mean[j] * chaos_sample;
      }
    }
  }
  for (i=1; i<numExpansionTerms; ++i) {
    term = num_pts * norm_squared(multiIndex[i]);
    if (configOptions.expansionCoeffFlag)
      expansionCoeffs[i] /= term;
    if (configOptions.expansionCoeffGradFlag) {
      exp_grad_i = expansionCoeffGrads[i];
      for (j=0; j<num_deriv_vars; ++j)
	exp_grad_i[j] /= term;
    }
#ifdef DEBUG
    PCout << "coeff[" << i << "] = " << expansionCoeffs[i]
        //<< "coeff_grad[" << i <<"] = " << exp_grad_i
	  << " norm squared[" << i <<"] = " << norm_squared(multiIndex[i])
	  << '\n';
#endif // DEBUG
  }
}


const Real& OrthogPolyApproximation::get_value(const RealVector& x)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response value prediction
  approxValue = 0.0;
  for (size_t i=0; i<numExpansionTerms; ++i)
    approxValue += expansionCoeffs[i]
                *  multivariate_polynomial(x, multiIndex[i]);
  return approxValue;
}


const RealVector& OrthogPolyApproximation::get_gradient(const RealVector& x)
{
  // this could define a default_dvv and call get_gradient(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_gradient()" << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != numVars)
    approxGradient.sizeUninitialized(numVars);
  approxGradient = 0.0;

  // sum expansion to get response gradient prediction
  size_t i, j;
  for (i=0; i<numExpansionTerms; ++i) {
    const RealVector& term_i_grad
      = multivariate_polynomial_gradient(x, multiIndex[i]);
    for (j=0; j<numVars; ++j)
      approxGradient[j] += expansionCoeffs[i]*term_i_grad[j];
  }
  return approxGradient;
}


const RealVector& OrthogPolyApproximation::
get_gradient(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_gradient()" << std::endl;
    abort_handler(-1);
  }

  size_t i, j, num_deriv_vars = dvv.size();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.0;

  // sum expansion to get response gradient prediction
  for (i=0; i<numExpansionTerms; ++i) {
    const RealVector& term_i_grad
      = multivariate_polynomial_gradient(x, multiIndex[i], dvv);
    for (j=0; j<num_deriv_vars; ++j)
      approxGradient[j] += expansionCoeffs[i]*term_i_grad[j];
  }
  return approxGradient;
}


/** In this case, all expansion variables are random variables and the
    mean of the expansion is simply the first chaos coefficient. */
const Real& OrthogPolyApproximation::get_mean()
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_mean()" << std::endl;
    abort_handler(-1);
  }

  //expansionMean = expansionCoeffs[0];
  //return expansionMean;
  return expansionCoeffs[0];
}


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves evaluating the
    expectation over this subset. */
const Real& OrthogPolyApproximation::get_mean(const RealVector& x)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_mean()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  expansionMean = expansionCoeffs[0];

  size_t i;
  SizetList::iterator it;
  for (i=1; i<numExpansionTerms; ++i) {
    bool include = true;
    for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
      if (multiIndex[i][*it]) // expectation of this expansion term is zero
	{ include = false; break; }
    if (include) {
      Real Psi = 1.0;
      for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	size_t index = *it;
	unsigned short order_1d = multiIndex[i][index];
	if (order_1d)
	  Psi *= polynomialBasis[index].get_value(x[index], order_1d);
      }
      expansionMean += Psi*expansionCoeffs[i];
#ifdef DEBUG
      PCout << "Mean estimate inclusion: term index = " << i
	    << " Psi = " << Psi << " PCE coeff = " << expansionCoeffs[i]
	    << " total = " << expansionMean << '\n';
#endif // DEBUG
    }
  }

  return expansionMean;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  In
    this case, the derivative of the expectation is the expectation of
    the derivative.  The mixed derivative case (some design variables
    are inserted and some are augmented) requires no special treatment. */
const RealVector& OrthogPolyApproximation::get_mean_gradient()
{
  // d/ds \mu_R = d/ds \alpha_0 = <dR/ds>

  // Error check for required data
  if (!configOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in "
	  << "OrthogPolyApproximation::get_mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  //expansionMeanGrad = expansionCoeffGrads[0];
  return expansionMeanGrad
    = Teuchos::getCol(Teuchos::Copy, expansionCoeffGrads, 0);
  return expansionMeanGrad;
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
get_mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  size_t i, j, deriv_index, num_deriv_vars = dvv.size();
  if (expansionMeanGrad.length() != num_deriv_vars)
    expansionMeanGrad.sizeUninitialized(num_deriv_vars);
  expansionMeanGrad = 0.0;

  SizetList::iterator it;
  size_t cntr = 0; // insertions carried in order within expansionCoeffGrads
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (randomVarsKey[deriv_index]) { // deriv w.r.t. design var insertion
      // Error check for required data
      if (!configOptions.expansionCoeffGradFlag) {
	PCerr << "Error: expansion coefficient gradients not defined in "
	      << "OrthogPolyApproximation::get_mean_gradient()." << std::endl;
	abort_handler(-1);
      }
      expansionMeanGrad[i] = expansionCoeffGrads[0][cntr];
    }
    else if (!configOptions.expansionCoeffFlag) { // Error check for reqd data
      PCerr << "Error: expansion coefficients not defined in "
	    << "OrthogPolyApproximation::get_mean_gradient()" << std::endl;
      abort_handler(-1);
    }
    for (j=1; j<numExpansionTerms; ++j) {
      bool include = true;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	if (multiIndex[j][*it]) // expectation of this expansion term is zero
	  { include = false; break; }
      if (include) {
	// In both cases below, the term to differentiate is alpha_j(s) Psi_j(s)
	// since <Psi_j>_xi = 1 for included terms.  The difference occurs based
	// on whether a particular s_i dependence appears in alpha (for
	// inserted) or Psi (for augmented).
	if (randomVarsKey[deriv_index]) {
	  // -------------------------------------------
	  // derivative w.r.t. design variable insertion
	  // -------------------------------------------
	  Real Psi = 1.0;
	  for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	    size_t index = *it;
	    unsigned short order_1d = multiIndex[j][index];
	    if (order_1d)
	      Psi *= polynomialBasis[index].get_value(x[index], order_1d);
	  }
	  expansionMeanGrad[i] += Psi*expansionCoeffGrads[j][cntr];
	}
	else {
	  // ----------------------------------------------
	  // derivative w.r.t. design variable augmentation
	  // ----------------------------------------------
	  Real term_j_grad_i = 1.0;
	  for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	    size_t index = *it;
	    term_j_grad_i *= (index == deriv_index) ?
	      polynomialBasis[index].get_gradient(x[index],multiIndex[j][index])
	    : polynomialBasis[index].get_value(x[index], multiIndex[j][index]);
	  }
	  expansionMeanGrad[i] += expansionCoeffs[j]*term_j_grad_i;
	}
      }
    }
    if (randomVarsKey[deriv_index]) // derivative w.r.t. design var insertion
      ++cntr;
  }

  return expansionMeanGrad;
}

// Computes centered (non-normalized) moments numerically
const RealVector& OrthogPolyApproximation::get_numeric_moments()
{
  size_t i, j, num_pts = dataPoints.size();
  size_t num_moments = 4; // may want this to be a user option?
  numericMoments.sizeUninitialized(num_moments);
  if (!anchorPoint.is_null())
    ++num_pts;
  if (!num_pts) {
    PCerr << "Error: nonzero number of sample points required in "
	  << "OrthogPolyApproximation::get_numeric_moments()." << std::endl;
    abort_handler(-1);
  }

  // Assumes coefficients have been calculated
  std::vector<SurrogateDataPoint>::const_iterator cit;
 
  // calculate moments 
  switch (configOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    // Define centered data
    RealVector c_function_vals(dataPoints.size());
    for (i=0,cit=dataPoints.begin(); cit!=dataPoints.end(); ++i,++cit) 
      c_function_vals[i] = cit->response_function() - expansionCoeffs[0];
    // presumably all exceptions are caught in computing coeffs of the expansion
    integrate_moments(multiIndex, dataPoints, c_function_vals,
		      driverRep->weight_sets(),	numericMoments);
    break;
  }
  case SPARSE_GRID:
    switch (sparseGridExpansion) {
    case TENSOR_INT_TENSOR_SUM_EXP: {
      // multiple expansion integration
      // Note: allocate_arrays() calls sparse_grid_multi_index() which uses
      // append_multi_index() to build multiIndex.
      SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
      const IntArray& sm_coeffs = ssg_driver->smolyak_coefficients();
      size_t i, num_tensor_grids = tpMultiIndex.size();
      std::vector<SurrogateDataPoint> tp_data_pts;

      RealVector tp_wts;//, tp_coeffs; RealMatrix tp_coeff_grads;
      RealVector tp_numeric_moments(num_moments);
      numericMoments = 0;
      // loop over tensor-products, forming sub-expansions, and sum them up
      for (i=0; i<num_tensor_grids; ++i) {
	// form tp_data_pts, tp_wts using collocKey et al.
	integration_data(i, tp_data_pts, tp_wts);
        // Define centered data
        RealVector c_function_vals(tp_data_pts.size());
        // center data
        for (j=0, cit=tp_data_pts.begin(); cit!=tp_data_pts.end(); ++j, ++cit)
          c_function_vals[j] = cit->response_function() - expansionCoeffs[0];
        // get moments for constituent tensor grid
        integrate_moments(tpMultiIndex[i], tp_data_pts, c_function_vals,
			  tp_wts, tp_numeric_moments);

	// add contribution to sparse grid
	for (j=0; j<num_moments;++j)
          numericMoments[j] += sm_coeffs[i]* tp_numeric_moments[j];

      }
      break;
    }
    default: // SPARSE_INT_*
      RealVector c_function_vals(dataPoints.size());
      for (i=0,cit=dataPoints.begin(); cit!=dataPoints.end(); ++i,++cit) 
        c_function_vals[i] = cit->response_function() - expansionCoeffs[0];
      integrate_moments(multiIndex, dataPoints, c_function_vals,
			driverRep->weight_sets(), numericMoments);
      break;
    }
    break;
  }
  return numericMoments;
}


/** In this case, all expansion variables are random variables and the
    variance of the expansion is the sum over all but the first term
    of the coefficients squared times the polynomial norms squared. */
const Real& OrthogPolyApproximation::get_variance()
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_variance()" << std::endl;
    abort_handler(-1);
  }

  return get_covariance(expansionCoeffs);
}


const Real& OrthogPolyApproximation::
get_covariance(const RealVector& exp_coeffs_2)
{
  expansionVariance = 0.0;
  for (size_t i=1; i<numExpansionTerms; ++i)
    expansionVariance += expansionCoeffs[i] * exp_coeffs_2[i] *
      norm_squared(multiIndex[i]);
  return expansionVariance;
}


/** In this case, a subset of the expansion variables are random variables
    and the variance of the expansion involves summations over this subset. */
const Real& OrthogPolyApproximation::get_variance(const RealVector& x)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_variance()" << std::endl;
    abort_handler(-1);
  }

  expansionVariance = 0.0;
  size_t i, j;
  SizetList::iterator it;
  for (i=1; i<numExpansionTerms; ++i) {
    // For r = random_vars and nr = non_random_vars,
    // sigma^2_R(nr) = < (R(r,nr) - \mu_R(nr))^2 >_r
    // -> include only those terms from R(r,nr) which do not appear in \mu_R(nr)
    bool include_i = false;
    for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
      if (multiIndex[i][*it]) // term does not appear in mean(nr)
	{ include_i = true; break; }
    if (include_i) {
      Real norm_sq_i = norm_squared_random(multiIndex[i]);
      for (j=1; j<numExpansionTerms; ++j) {

	// random part of polynomial must be identical to contribute to variance
	// (else orthogonality drops term)
	bool include_j = true;
	for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	  if (multiIndex[i][*it] != multiIndex[j][*it])
	    { include_j = false; break; }
	if (include_j) {
	  Real var_ij = expansionCoeffs[i] * expansionCoeffs[j] * norm_sq_i;
	  for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end(); ++it) {
	    size_t index = *it;
	    unsigned short order_1d_i = multiIndex[i][index],
	      order_1d_j = multiIndex[j][index];
	    BasisPolynomial& poly_1d = polynomialBasis[index];
	    if (order_1d_i)
	      var_ij *= poly_1d.get_value(x[index], order_1d_i);
	    if (order_1d_j)
	      var_ij *= poly_1d.get_value(x[index], order_1d_j);
	  }
	  expansionVariance += var_ij;
#ifdef DEBUG
	  PCout << "Variance estimate inclusion: term index = " << i
		<< " variance = " << var_ij << " total = " << expansionVariance
		<<'\n';
#endif // DEBUG
	}
      }
    }
  }

  return expansionVariance;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  The
    mixed derivative case (some design variables are inserted and some
    are augmented) requires no special treatment. */
const RealVector& OrthogPolyApproximation::get_variance_gradient()
{
  // d/ds \sigma^2_R = Sum_{j=1}^P <Psi^2_j> d/ds \alpha^2_j
  //                 = 2 Sum_{j=1}^P \alpha_j <dR/ds, Psi_j>

  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_variance_gradient()" << std::endl;
    abort_handler(-1);
  }
  if (!configOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in "
	  << "OrthogPolyApproximation::get_variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  size_t i, j, num_deriv_vars = expansionCoeffGrads.numRows();
  if (expansionVarianceGrad.length() != num_deriv_vars)
    expansionVarianceGrad.sizeUninitialized(num_deriv_vars);
  expansionVarianceGrad = 0.0;
  for (i=1; i<numExpansionTerms; ++i) {
    Real term_i = 2. * expansionCoeffs[i] * norm_squared(multiIndex[i]);
    for (j=0; j<num_deriv_vars; ++j)
      expansionVarianceGrad[j] += term_i * expansionCoeffGrads[i][j];
  }
  return expansionVarianceGrad;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  This function must handle the mixed case, where
    some design/state variables are augmented (and are part of the
    expansion) and some are inserted (derivatives are obtained from
    expansionCoeffGrads). */
const RealVector& OrthogPolyApproximation::
get_variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!configOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_variance_gradient()" << std::endl;
    abort_handler(-1);
  }

  size_t i, j, k, deriv_index, num_deriv_vars = dvv.size();
  if (expansionVarianceGrad.length() != num_deriv_vars)
    expansionVarianceGrad.sizeUninitialized(num_deriv_vars);
  expansionVarianceGrad = 0.0;

  SizetList::iterator it;
  size_t cntr = 0; // insertions carried in order within expansionCoeffGrads
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (randomVarsKey[deriv_index] && !configOptions.expansionCoeffGradFlag) {
      PCerr << "Error: expansion coefficient gradients not defined in "
	    << "OrthogPolyApproximation::get_variance_gradient()." << std::endl;
      abort_handler(-1);
    }
    for (j=1; j<numExpansionTerms; ++j) {
      bool include_j = false;
      for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	if (multiIndex[j][*it]) // term does not appear in mean(nr)
	  { include_j = true; break; }
      if (include_j) {
	Real norm_sq_j = norm_squared_random(multiIndex[j]);
	for (k=1; k<numExpansionTerms; ++k) {
	  // random part of polynomial must be identical to contribute to
	  // variance (else orthogonality drops term)
	  bool include_k = true;
	  for (it=randomIndices.begin(); it!=randomIndices.end(); ++it)
	    if (multiIndex[j][*it] != multiIndex[k][*it])
	      { include_k = false; break; }
	  if (include_k) {
	    // In both cases below, the term to differentiate is
	    // alpha_j(s) alpha_k(s) <Psi_j^2>_xi Psi_j(s) Psi_k(s) and the
	    // difference occurs based on whether a particular s_i dependence
	    // appears in alpha (for inserted) or Psi (for augmented).
	    if (randomVarsKey[deriv_index]) {
	      // -------------------------------------------
	      // derivative w.r.t. design variable insertion
	      // -------------------------------------------
	      Real var_jk = norm_sq_j *
		(expansionCoeffs[j] * expansionCoeffGrads[k][cntr] +
		 expansionCoeffs[k] * expansionCoeffGrads[j][cntr]);
	      for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end();
		   ++it) {
		size_t index = *it;
		unsigned short order_1d_j = multiIndex[j][index],
		  order_1d_k = multiIndex[k][index];
		BasisPolynomial& poly_1d = polynomialBasis[index];
		if (order_1d_j)
		  var_jk *= poly_1d.get_value(x[index], order_1d_j);
		if (order_1d_k)
		  var_jk *= poly_1d.get_value(x[index], order_1d_k);
	      }
	      expansionVarianceGrad[i] += var_jk;
	    }
	    else {
	      // ----------------------------------------------
	      // derivative w.r.t. design variable augmentation
	      // ----------------------------------------------
	      Real var_jk = expansionCoeffs[j] * expansionCoeffs[k] * norm_sq_j;
	      Real Psi_j = 1., Psi_k = 1., dPsi_j_ds_i = 1., dPsi_k_ds_i = 1.;
	      for (it=nonRandomIndices.begin(); it!=nonRandomIndices.end();
		   ++it) {
		size_t index = *it;
		unsigned short order_1d_j = multiIndex[j][index],
		  order_1d_k = multiIndex[k][index];
		BasisPolynomial& poly_1d = polynomialBasis[index];
		if (order_1d_j)
		  Psi_j *= poly_1d.get_value(x[index], order_1d_j);
		if (order_1d_k)
		  Psi_k *= poly_1d.get_value(x[index], order_1d_k);
		dPsi_j_ds_i *= (index == deriv_index) ?
		  poly_1d.get_gradient(x[index], order_1d_j) :
		  poly_1d.get_value(   x[index], order_1d_j);
		dPsi_k_ds_i *= (index == deriv_index) ?
		  poly_1d.get_gradient(x[index], order_1d_k) :
		  poly_1d.get_value(   x[index], order_1d_k);
	      }
	      expansionVarianceGrad[i] += var_jk * 
		(Psi_j*dPsi_k_ds_i + dPsi_j_ds_i*Psi_k);
	    }
	  }
	}
      }
    }
    if (randomVarsKey[deriv_index]) // derivative w.r.t. design var insertion
      ++cntr;
  }

  return expansionVarianceGrad;
}


/** This test works in combination with DEBUG settings in
    (Legendre,Laguerre,Jacobi,GenLaguerre)OrthogPolynomial::get_gradient(). */
void OrthogPolyApproximation::gradient_check()
{
  BasisPolynomial hermite_poly(HERMITE), legendre_poly(LEGENDRE),
    laguerre_poly(LAGUERRE), jacobi_poly(JACOBI),
    gen_laguerre_poly(GENERALIZED_LAGUERRE), chebyshev_poly(CHEBYSHEV);
  // alpha/beta selections mirror dakota_uq_rosenbrock_pce.in
  jacobi_poly.alpha_stat(1.5);
  jacobi_poly.beta_stat(2.);
  gen_laguerre_poly.alpha_stat(2.5);

  Real x = 0.5; // valid for all support ranges: [-1,1], [0,Inf], [-Inf, Inf]
  PCout << "-------------------------------------------------\n";
  for (size_t n=0; n<=10; n++) {
    PCout << "Gradients at " << x << " for order " << n << '\n';
    hermite_poly.get_gradient(x, n);
    legendre_poly.get_gradient(x, n);
    laguerre_poly.get_gradient(x, n);
    jacobi_poly.get_gradient(x, n);
    gen_laguerre_poly.get_gradient(x, n);
    chebyshev_poly.get_gradient(x, n);
    PCout << "-------------------------------------------------\n";
  }
}


const Real& OrthogPolyApproximation::norm_squared(const UShortArray& indices)
{
  // the norm squared of a particular multivariate polynomial is the product of
  // the norms squared of the numVars univariate polynomials that comprise it.
  multiPolyNormSq = 1.0;
  unsigned short order;
  for (size_t i=0; i<numVars; ++i) {
    order = indices[i];
    if (order)
      multiPolyNormSq *= polynomialBasis[i].norm_squared(order);
  }
  return multiPolyNormSq;
}


const Real& OrthogPolyApproximation::
norm_squared_random(const UShortArray& indices)
{
  // the norm squared of a particular multivariate polynomial is the product of
  // the norms squared of the numVars univariate polynomials that comprise it.
  multiPolyNormSq = 1.0;
  unsigned short order;
  for (size_t i=0; i<numVars; ++i) {
    if (randomVarsKey[i]) {
      order = indices[i];
      if (order)
	multiPolyNormSq *= polynomialBasis[i].norm_squared(order);
    }
  }
  return multiPolyNormSq;
}


void OrthogPolyApproximation::compute_component_effects()
{
  // sobolIndices are index via binary number represenation
  // e.g. if there are 5 variables -> 5 bit represenation
  // if a variable belongs to a subset \alpha, it has a value of 1; otherwise 0
  // | 0 | 0 | 0 | 0 | 1 | represents the subset that only contains variable 1
  // sobolIndices[0] is set to 1; to calculate it is redundant

  // Index the set using binary number representation
  size_t index_bin, i, j;
  // Compute a pseudo-variance that makes no distinction between probabilistic
  // variables and non-probabilistic variables
  Real p_var = 0, p_var_i;
  for (i=1; i<numExpansionTerms; i++) 
    p_var += norm_squared(multiIndex[i])*expansionCoeffs(i)*expansionCoeffs(i);

  // iterate through multiIndex and store sensitivities
  sobolIndices    = 0.; // initialize
  sobolIndices[0] = 1.; // just a place holder; zero index is never invoked
  for (i=1; i<numExpansionTerms; ++i) {
    index_bin = 0;
    p_var_i = norm_squared(multiIndex[i]) * expansionCoeffs(i) *
      expansionCoeffs(i) / p_var;
    for (j=0; j<numVars; ++j)
      if (multiIndex[i][j]) // convert this subset into binary number
	index_bin += (size_t)pow(2.,(int)j);
    // If term is main effect (exists in map), keep; otherwise, discard
    if (sobolIndexMap.find(index_bin) != sobolIndexMap.end())
      sobolIndices[sobolIndexMap[index_bin]] += p_var_i;
  }
#ifdef DEBUG
  PCout << "In OrthogPolyApproximation::compute_component_effects(), "
	<< "sobolIndices =\n"; write_data(PCout, sobolIndices);
#endif // DEBUG
}


void OrthogPolyApproximation::compute_total_effects() 
{
  // iterate through existing indices if all component indices are available
  totalSobolIndices = 0.;
  if (configOptions.vbdControl == ALL_VBD) {
    for (IntIntMIter itr=sobolIndexMap.begin(); itr!=sobolIndexMap.end(); ++itr)
      for (int k=0; k<numVars; k++) 
        if (itr->first & (1 << k))
          totalSobolIndices[k] += sobolIndices[itr->second];
  }
  // otherwise, iterate over the expansion terms as it is more 
  // computationally efficient than performing ANOVA operators
  else {
    // Index the set using binary number representation
    size_t index_bin, i, j;
    Real p_var = 0, p_var_i;
    for (i=1; i<numExpansionTerms; i++) 
      p_var += norm_squared(multiIndex[i]) * expansionCoeffs(i)
	* expansionCoeffs(i);
  
    // Computing total indices by iterating through expansion terms is simpler
    // and more computationally efficient than computing via ANOVA operators 
    for (i=1; i<numExpansionTerms; ++i) {
      index_bin = 0;
      p_var_i = norm_squared(multiIndex[i]) * expansionCoeffs(i) *
        expansionCoeffs(i) / p_var;
      for (j=0; j<numVars; ++j) {
        if (multiIndex[i][j]) {
          // convert this subset multiIndex[i] into binary number
          index_bin += (size_t)pow(2.,(int)j);
          // for any constituent variable j in exansion term i, the expansion
          // term contributes to the total sensitivty of variable j
          totalSobolIndices[j] += p_var_i;
        }
      }
    }
  }
#ifdef DEBUG
  PCout << "In OrthogPolyApproximation::compute_total_effects(), "
	<< "totalSobolIndices =\n"; write_data(PCout, totalSobolIndices);
#endif // DEBUG
}


void OrthogPolyApproximation::print_coefficients(std::ostream& s) const
{
  size_t i, j;
  char tag[10];

  // terms and term identifiers
  for (i=0; i<numExpansionTerms; ++i) {
    s << "\n  " << std::setw(17) << expansionCoeffs[i];
    for (j=0; j<numVars; ++j) {
      switch (basisTypes[j]) {
      case HERMITE:
	std::sprintf(tag,  "He%i", multiIndex[i][j]);
	break;
      case LEGENDRE:
	std::sprintf(tag,   "P%i", multiIndex[i][j]);
	break;
      case LAGUERRE:
	std::sprintf(tag,   "L%i", multiIndex[i][j]);
	break;
      case JACOBI:
	std::sprintf(tag, "Pab%i", multiIndex[i][j]);
	break;
      case GENERALIZED_LAGUERRE:
	std::sprintf(tag,  "La%i", multiIndex[i][j]);
	break;
      case CHEBYSHEV:
	std::sprintf(tag,   "T%i", multiIndex[i][j]);
	break;
      case NUMERICALLY_GENERATED:
	std::sprintf(tag, "Num%i", multiIndex[i][j]);
	break;
      default:
	PCerr << "Error: bad polynomial type = " << basisTypes[j] << " in "
	      << "OrthogPolyApproximation::print_coefficients()." << std::endl;
	abort_handler(-1);
	break;
      }
      s << std::setw(5) << tag;
    }
  }
  s << '\n';
}

} // namespace Pecos
