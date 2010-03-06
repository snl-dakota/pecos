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
        gauss_modes[i] = (short)int_rules[i];//                           break;
      //default:
      //  gauss_modes[i] = 0;                                             break;
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
  if (expansionCoeffFlag || expansionGradFlag)
    switch (expCoeffsSolnApproach) {
    case QUADRATURE: case SPARSE_GRID: case SAMPLING:
      return 1; // quadrature: (int)pow((Real)MIN_GAUSS_PTS, numVars);
      break;
    case REGRESSION:
      // numExpansionTerms is either set from the NonDPolynomialChaos ctor if
      // expansion_terms is specified or computed by the allocate_arrays() call
      // in find_coefficients() if expansion_order is specified.  The latter
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


void OrthogPolyApproximation::find_coefficients()
{
  if (!expansionCoeffFlag && !expansionGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in "
	  << "OrthogPolyApproximation::find_coefficients().\n         "
	  << "Bypassing approximation construction." << std::endl;
    return;
  }

  // For testing of anchorPoint logic:
  //anchorPoint = dataPoints.front();
  //dataPoints.pop_front();

  // anchorPoint, if present, is handled differently for different
  // expCoeffsSolnApproach settings:
  //   SAMPLING:               treat it as another currentPoint
  //   QUADRATURE/SPARSE_GRID: error
  //   REGRESSION:             use equality-constrained least squares
  size_t i, j, num_pts = dataPoints.size();
  if (!anchorPoint.is_null())
    ++num_pts;
  if (!num_pts) {
    PCerr << "Error: nonzero number of sample points required in "
	  << "OrthogPolyApproximation::find_coefficients()." << std::endl;
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
  switch (expCoeffsSolnApproach) {
  case QUADRATURE: {
    if (!driverRep) {
      PCerr << "Error: pointer to integration driver required in "
	    << "OrthogPolyApproximation::find_coefficients()." << std::endl;
      abort_handler(-1);
    }
    // verify quad_order stencil matches num_pts
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray& quad_order = tpq_driver->quadrature_order();
    if (quad_order.size() != numVars) {
      PCerr << "Error: quadrature order array is not consistent with number of "
	    << "variables (" << numVars << ")\n       in "
	    << "OrthogPolyApproximation::find_coefficients()." << std::endl;
      abort_handler(-1);
    }
    size_t num_gauss_pts = 1;
    for (i=0; i<numVars; ++i)
      num_gauss_pts *= quad_order[i];
    if (num_pts != num_gauss_pts ||
	num_pts != driverRep->weight_sets().length()) {
      PCerr << "Error: number of current points (" << num_pts << ") is not "
	    << "consistent with\n       quadrature data in "
	    << "OrthogPolyApproximation::find_coefficients()." << std::endl;
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

    integration();
    break;
  }
  case SPARSE_GRID:
    if (!driverRep) {
      PCerr << "Error: pointer to integration driver required in "
	    << "OrthogPolyApproximation::find_coefficients()." << std::endl;
      abort_handler(-1);
    }
    if (num_pts != driverRep->weight_sets().length()) {
      PCerr << "Error: number of current points (" << num_pts << ") is not "
	    << "consistent with\n       sparse grid data in "
	    << "OrthogPolyApproximation::find_coefficients()." << std::endl;
      abort_handler(-1);
    }

    integration();
    break;
  case REGRESSION:
    regression();
    break;
  case SAMPLING:
    expectation();
    break;
  }
}


void OrthogPolyApproximation::allocate_arrays()
{
  PolynomialApproximation::allocate_arrays();

  // Infer expansion formulation from quadrature_order or sparse_grid_level
  // spec, as in SC.  Preserve previous capability (quadrature_order and
  // sparse_grid_level with total-order expansions) for paper results via
  // quadratureExpansion/sparseGridExpansion (compile-time) switches.

  // update_exp_form controls when to update (refinement) and when not to
  // update (subIterator execution) an expansion's multiIndex definition.
  // Simple logic of updating if previous number of points != current number
  // is not robust enough for anisotropic updates --> track using Prev arrays.
  switch (expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray& quad_order = tpq_driver->quadrature_order();
    bool update_exp_form = (quad_order != quadOrderPrev);
    quadOrderPrev = quad_order;
    // *** TO DO: capture updates to parameterized/numerical polynomials?

    if (update_exp_form) {
      UShortArray integrand_order(numVars);
      quadrature_order_to_integrand_order(quad_order, integrand_order);
      integrand_order_to_expansion_order(integrand_order, approxOrder);
    }
    PCout << "Orthogonal polynomial approximation order = { ";
    for (size_t i=0; i<numVars; ++i)
      PCout << approxOrder[i] << ' ';
    if (quadratureExpansion == TENSOR_PRODUCT) {
      if (update_exp_form)
	tensor_product_multi_index(approxOrder, multiIndex, true);
      PCout << "} using tensor-product expansion of ";
    }
    else if (quadratureExpansion == TOTAL_ORDER) {
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
    break;
  }
  case SPARSE_GRID: {
    const RealVector& aniso_wts  = driverRep->anisotropic_weights();
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    unsigned short    ssg_level  = ssg_driver->level();
    bool update_exp_form
      = (ssg_level != ssgLevelPrev || aniso_wts != ssgAnisoWtsPrev);
    ssgLevelPrev = ssg_level; ssgAnisoWtsPrev = aniso_wts;
    // *** TO DO: capture updates to parameterized/numerical polynomials?

    if (update_exp_form) {
      // set sparseGridExpansion differently for isotropic vs. anisotropic SSG
      sparseGridExpansion = (driverRep->isotropic()) ?
	TOTAL_ORDER : TENSOR_PRODUCT_SUM;
      // compute and output number of terms
      sparse_grid_multi_index(multiIndex);
      numExpansionTerms = multiIndex.size();
    }
    if (sparseGridExpansion == TENSOR_PRODUCT_SUM)
      PCout << "Orthogonal polynomial approximation level = " << ssg_level
	    << " using sparse grid expansion of " << numExpansionTerms
	    << " terms\n";
    else if (sparseGridExpansion == TOTAL_ORDER ||
	     sparseGridExpansion == HEURISTIC_TOTAL_ORDER) {
      PCout << "Orthogonal polynomial approximation order = { ";
      for (size_t i=0; i<numVars; ++i)
	PCout << approxOrder[i] << ' ';
      PCout << "} using total-order expansion of " << numExpansionTerms
	    << " terms\n";
    }
    else {
      PCerr << "Error: unsupported setting for sparseGridExpansion in "
	    << "OrthogPolyApproximation::allocate_arrays()." << std::endl;
      abort_handler(-1);
    }
    break;
  }
  default: { // SAMPLING and REGRESSION
    // Support efficiency in PCBDO for now.
    bool update_exp_form
      = (expansionCoeffs.empty() && expansionCoeffGrads.empty());
    // TO DO: support uniform/adaptive refinement.
    //bool update_exp_form = (approxOrder       != approxOrderPrev ||
    //                        numExpansionTerms != numExpTermsPrev);
    //approxOrderPrev = approxOrder;
    //numExpTermsPrev = numExpansionTerms;

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
    break;
  }
  }

  // now that terms & order are known, shape some arrays.  This is done here,
  // rather than in find coefficients(), in order to support array sizing for
  // the data import case.
  if (expansionCoeffFlag && expansionCoeffs.length() != numExpansionTerms)
    expansionCoeffs.sizeUninitialized(numExpansionTerms);
  if (expansionGradFlag) {
    const SurrogateDataPoint& sdp = (anchorPoint.is_null()) ?
      *dataPoints.begin() : anchorPoint;
    size_t num_deriv_vars = sdp.response_gradient().length();
    if (expansionCoeffGrads.numRows() != num_deriv_vars ||
	expansionCoeffGrads.numCols() != numExpansionTerms)
      expansionCoeffGrads.shapeUninitialized(num_deriv_vars, numExpansionTerms);
  }
}


void OrthogPolyApproximation::
sparse_grid_multi_index(UShort2DArray& multi_index)
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;

  // Smolyak lower bound is w-n+1 --> offset is n-1
  UShort2DArray sm_multi_index; RealArray sm_coeffs;
  smolyak_multi_index(sm_multi_index, sm_coeffs);
  size_t i, num_smolyak_indices = sm_multi_index.size();
#ifdef DEBUG
  PCout << "num_smolyak_indices = " << num_smolyak_indices
	<< "\nsm_multi_index =\n" << sm_multi_index << '\n';
#endif // DEBUG

  UShortArray quad_order(numVars), integrand_order(numVars);
  if (sparseGridExpansion == TENSOR_PRODUCT_SUM) {
    // assemble a complete list of individual polynomial coverage
    // defined from the linear combination of mixed tensor products
    multi_index.clear();
    UShort2DArray tp_multi_index;
    UShortArray expansion_order(numVars);
    for (i=0; i<num_smolyak_indices; ++i) {
      ssg_driver->level_to_order(sm_multi_index[i], quad_order);
      quadrature_order_to_integrand_order(quad_order, integrand_order);
      integrand_order_to_expansion_order(integrand_order, expansion_order);
      tensor_product_multi_index(expansion_order, tp_multi_index, true);
      append_unique(tp_multi_index, multi_index);
#ifdef DEBUG
      PCout << "level =\n" << sm_multi_index[i] << "quadrature_order =\n"
	    << quad_order << "integrand_order =\n" << integrand_order
	    << "expansion_order =\n" << expansion_order << "tp_multi_index =\n"
	    << tp_multi_index << '\n';//<< "multi_index =\n" << multi_index;
#endif // DEBUG
      tp_multi_index.clear();
    }
  }
  else if (sparseGridExpansion == TOTAL_ORDER) {
    // back out approxOrder & use total_order_multi_index()
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
  }
  else if (sparseGridExpansion == HEURISTIC_TOTAL_ORDER) { // previous heuristic
    sparse_grid_level_to_expansion_order(ssg_driver->level(), approxOrder);
    total_order_multi_index(approxOrder, multi_index);
  }
  else { // possible addtnl fallback: revert to [optional] expansion_order spec.
    PCerr << "Error: unsupported sparseGridExpansion option in "
	  << "OrthogPolyApproximation::sparse_grid_multi_index()" << std::endl;
    abort_handler(-1);
  }
}


/* These mappings reflect simplified heuristics (expected to be exact
   for CC, Gauss-Patterson, and linear growth Gaussian, but not for
   exponential Gaussian). */
void OrthogPolyApproximation::
sparse_grid_level_to_expansion_order(unsigned short ssg_level,
				     UShortArray& exp_order)
{
  if (exp_order.size() != numVars)
    exp_order.resize(numVars);
  const IntArray& rules = ((SparseGridDriver*)driverRep)->integration_rules();
  for (size_t i=0; i<numVars; ++i) {
    switch (rules[i]) {
    case CLENSHAW_CURTIS: case CLENSHAW_CURTIS_SLOW:
    case FEJER2:          case FEJER2_SLOW: // TO DO: verify
      // integrand order = 2*level+1 ("sharp" result from Novak & Ritter, 1996)
    case GAUSS_PATTERSON_SLOW: {
      // level is scaled back to be consistent with CC guarantee above:
      // -> consistency is good for mixing & matching anisotropic rules, but...
      // -> seems overly conservative: for isotropic, more desirable to scale
      //    back only where needed to avoid overshooting total-order precision.
      unsigned short integrand = 2*ssg_level + 1;
      exp_order[i] = integrand / 2; // remainder truncated
      // results in exp_order = level
      break;
    }
    /*
    case GAUSS_PATTERSON: { // TO DO
      // Exponential growth: o(l) = 2^{l+1}-1
      // integrand = 2o-1 (Gaussian) - o(l-1) (constraints from previous rule)
      //           + 1 (odd order rule -> symmetry)
      // *** Can be applied on a per-tensor product basis but not to full SSG.
      unsigned short integrand = (ssg_level) ?
	3*(unsigned short)std::pow(2.,(int)ssg_level) - 1 : 1;
      exp_order[i] = integrand / 2; // remainder truncated
      break;
    }
    */
    default: {
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
  else
    for (i=0; i<n; ++i)
      switch (gaussModes[i]) {
      case CLENSHAW_CURTIS:          case FEJER2:
      case CLENSHAW_CURTIS_MODERATE: case FEJER2_MODERATE:
      case CLENSHAW_CURTIS_SLOW:     case FEJER2_SLOW:
	// i = m (odd m), m-1 (even m).  Note that growth rule enforces odd.
	// TO DO: verify FEJER2 same as CC
	int_order[i] = (quad_order[i] % 2) ? quad_order[i] : quad_order[i] - 1;
	break;
      case GAUSS_PATTERSON: case GAUSS_PATTERSON_MODERATE:
      case GAUSS_PATTERSON_SLOW: {
	// for o(l)=2^{l+1}-1, o(l-1) = (o(l)-1)/2
	unsigned short previous = std::max(1,(quad_order[i] - 1)/2);
	int_order[i] = 2*quad_order[i] - previous;
	break;
      }
      default: // standard non-nested Gauss rules
	int_order[i] =  2*quad_order[i] - 1; break; // i = 2m - 1
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
append_unique(const UShort2DArray& tp_multi_index, UShort2DArray& multi_index)
{
  if (multi_index.empty())
    multi_index = tp_multi_index;
  else {
    size_t i, num_tp_mi = tp_multi_index.size();
    for (i=0; i<num_tp_mi; ++i) {
      const UShortArray& search_mi = tp_multi_index[i];
      if (std::find(multi_index.begin(), multi_index.end(), search_mi) ==
	  multi_index.end()) // search_mi does not yet exist in multi_index
	multi_index.push_back(search_mi);
    }
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


/** The coefficients of the PCE for the response are calculated using a
    Galerkin projection of the response against each multivariate orthogonal
    polynomial basis fn using the inner product ratio <f,Psi>/<Psi^2>, where
    inner product <a,b> is the n-dimensional integral of a*b*weighting over
    the support range of the n-dimensional (composite) weighting function.
    1-D quadrature rules are defined for specific 1-D weighting functions
    and support ranges and approximate the integral of f*weighting as the
    Sum_i of w_i f_i.  To extend this to n-dimensions, a tensor product
    quadrature rule or Smolyak sparse grid rule is applied using the product
    of 1-D weightings applied to the n-dimensional stencil of points.  It is
    not necessary to approximate the integral for the denominator numerically,
    since this is available analytically. */
void OrthogPolyApproximation::integration()
{
  if (!anchorPoint.is_null()) { // TO DO: verify this creates a problem
    PCerr << "Error: anchor point not supported for numerical integration in "
	  << "OrthogPolyApproximation::integration()." << std::endl;
    abort_handler(-1);
  }

  // Perform tensor-product quadrature/Smolyak sparse grid using
  // point & weight sets computed in TensorProductDriver/SparseGridDriver
  size_t i, j, k, num_deriv_vars = expansionCoeffGrads.numRows();
  const RealVector& wt_sets = driverRep->weight_sets();
  std::list<SurrogateDataPoint>::iterator it;
  Real empty_r;
  for (i=0; i<numExpansionTerms; ++i) {

    Real& chaos_coeff_i = (expansionCoeffFlag) ? expansionCoeffs[i] : empty_r;
    Real* chaos_grad_i  = (expansionGradFlag)  ? expansionCoeffGrads[i] : NULL;
    if (expansionCoeffFlag) chaos_coeff_i = 0.;
    if (expansionGradFlag)  std::fill_n(chaos_grad_i, num_deriv_vars, 0.);
    for (j=0, it=dataPoints.begin(); it!=dataPoints.end(); ++j, ++it) {

      // sum contribution to inner product <f, Psi_j>
      //Real chaos_sample
      //  = multivariate_polynomial(it->continuous_variables(), i)
      //num += wt_sets[j] * chaos_sample * it->response_function();
      // sum contribution to inner product <Psi_j, Psi_j>
      //den += wt_sets[j] * pow(chaos_sample, 2);

      Real wt_samp_prod = wt_sets[j] *
	multivariate_polynomial(it->continuous_variables(), i);
      if (expansionCoeffFlag)
	chaos_coeff_i += wt_samp_prod * it->response_function();
      if (expansionGradFlag) {
	const RealVector curr_pt_grad = it->response_gradient();
	for (k=0; k<num_deriv_vars; ++k)
	  chaos_grad_i[k] += wt_samp_prod * curr_pt_grad[k];
      }
    }

    //expansionCoeffs[i] = num / den;
    const Real& norm_sq = norm_squared(i);
    if (expansionCoeffFlag)
      chaos_coeff_i /= norm_sq;
    if (expansionGradFlag)
      for (k=0; k<num_deriv_vars; ++k)
	chaos_grad_i[k] /= norm_sq;
#ifdef DEBUG
    PCout << "coeff[" << i <<"] = " << chaos_coeff_i
        //<< "coeff_grad[" << i <<"] = " << chaos_grad_i
	  << " norm squared[" << i <<"] = " << norm_sq << "\n\n";
#endif // DEBUG
  }
}


/** In this case, regression is used in place of Galerkin projection.  That
    is, instead of calculating the PCE coefficients using inner product
    ratios, linear least squares is used to estimate the PCE coefficients
    which best match a set of response samples.  This approach is also known
    as stochastic response surfaces.  The least squares estimation is
    performed using DGELSS (SVD) or DGGLSE (equality-constrained) from
    LAPACK, based on the presence of an anchorPoint. */
void OrthogPolyApproximation::regression()
{
  size_t i, j, num_pts = dataPoints.size(), cntr = 0;
  int info       = 0; // output flag from DGELSS/DGGLSE subroutines
  int num_rows_A = num_pts;           // Number of rows in matrix A
  int num_cols_A = numExpansionTerms; // Number of columns in matrix A
  // Matrix of polynomial terms unrolled into a vector.
  double* A_matrix = new double [num_rows_A*num_cols_A]; // "A" in A*x = b
  std::list<SurrogateDataPoint>::iterator it;
  Teuchos::LAPACK<int, Real> la;
  bool err_flag = false;

  if (anchorPoint.is_null()) {
    // Use DGELSS for LLS soln using SVD.  Solves min ||b - Ax||_2
    // where {b} may have multiple RHS -> multiple {x} solutions.

    size_t num_deriv_vars = (expansionGradFlag)  ?
      expansionCoeffGrads.numRows() : 0;
    size_t num_coeff_rhs  = (expansionCoeffFlag) ? 1 : 0;
    int num_rhs  = num_coeff_rhs + num_deriv_vars; // input: # of RHS vectors
    double rcond = -1.; // input:  use macheps to rank singular vals of A
    int rank     =  0;  // output: effective rank of matrix A
    // input: vectors of response data that correspond to samples in matrix A.
    double* b_vectors = new double [num_rows_A*num_rhs]; // "b" in A*x = b
    // output: vector of singular values, dimensioned for overdetermined system
    double* s_vector  = new double [num_cols_A];

    // Get the optimal work array size
    int lwork    = -1; // special code for workspace query
    double* work = new double [1]; // temporary work array
    la.GELSS(num_rows_A, num_cols_A, num_rhs, A_matrix, num_rows_A, b_vectors,
	     num_rows_A, s_vector, rcond, &rank, work, lwork, &info);
    lwork = (int)work[0]; // optimal work array size returned by query
    delete [] work;
    work = new double [lwork]; // Optimal work array

    // The "A" matrix is a contiguous block of memory packed in column-major
    // ordering as required by F77 for the GELSS subroutine from LAPACK.  For
    // example, the 6 elements of A(2,3) are stored in the order A(1,1), A(2,1),
    // A(1,2), A(2,2), A(1,3), A(2,3).
    for (i=0; i<numExpansionTerms; ++i)
      for (it=dataPoints.begin(); it!=dataPoints.end(); ++it, ++cntr)
	A_matrix[cntr] = multivariate_polynomial(it->continuous_variables(), i);

    // response data (values/gradients) define the multiple RHS which are
    // matched in the LS soln.  b_vectors is num_pts (rows) x num_rhs (cols),
    // arranged in column-major order.
    for (i=0, it=dataPoints.begin(); i<num_rows_A; ++i, ++it) {
      if (expansionCoeffFlag)
	b_vectors[i] = it->response_function(); // i-th point: response value
      if (expansionGradFlag) {
	const RealVector curr_pt_grad = it->response_gradient();
	for (j=0; j<num_deriv_vars; ++j) // i-th point, j-th gradient component
	  b_vectors[(j+num_coeff_rhs)*num_pts+i] = curr_pt_grad[j];
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
    if (expansionCoeffFlag)
      copy_data(b_vectors, numExpansionTerms, expansionCoeffs);
    if (expansionGradFlag)
      for (i=0; i<numExpansionTerms; ++i)
	for (j=0; j<num_deriv_vars; ++j)
	  expansionCoeffGrads(j,i) = b_vectors[(j+num_coeff_rhs)*num_pts+i];

    delete [] b_vectors;
    delete [] s_vector;
    delete [] work;
  }
  else {
    // Use DGGLSE for equality-constrained LLS soln using GRQ factorization.
    // Solves min ||b - Ax||_2 s.t. Cx = d (Note: b,C switched from LAPACK docs)
    // where {b,d} are single vectors (multiple RHS not supported).

    size_t offset = 1;
    num_pts += offset;

    int num_cons = 1; // only 1 anchor point constraint
    // Vector of response values that correspond to the samples in matrix A.
    double* b_vector = new double [num_rows_A]; // "b" in A*x = b
    // Matrix of constraints unrolled into a vector
    double* C_matrix = new double [num_cols_A]; // "C" in C*x = d
    // RHS of constraints
    double* d_vector = new double [num_cons];   // "d" in C*x = d
    // Solution vector
    double* x_vector = new double [num_cols_A]; // "x" in b - A*x, C*x = d

    // Get the optimal work array size
    int lwork    = -1; // special code for workspace query
    double* work = new double [1]; // temporary work array
    la.GGLSE(num_rows_A, num_cols_A, num_cons, A_matrix, num_rows_A, C_matrix,
	     num_cons, b_vector, d_vector, x_vector, work, lwork, &info);
    lwork = (int)work[0]; // optimal work array size returned by query
    delete [] work;
    work = new double [lwork]; // Optimal work array

    if (expansionCoeffFlag) {
      // pack dataPoints in A and anchor point in C.  The "A" matrix is a
      // contiguous block of memory packed in column-major ordering as required
      // by F77 for the GGLSE subroutine from LAPACK.
      for (i=0; i<numExpansionTerms; ++i) {
	for (it=dataPoints.begin(); it!=dataPoints.end(); ++it, ++cntr)
	  A_matrix[cntr]
	    = multivariate_polynomial(it->continuous_variables(), i);
	C_matrix[i]
	  = multivariate_polynomial(anchorPoint.continuous_variables(), i);
      }
      for (it=dataPoints.begin(), i=0; it!=dataPoints.end(); ++it, ++i)
	b_vector[i] = it->response_function();
      d_vector[0] = anchorPoint.response_function();    // anchor data
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

    if (expansionGradFlag) {
      size_t num_deriv_vars = expansionCoeffGrads.numRows();
      for (i=0; i<num_deriv_vars; ++i) {
	// must be recomputed each time since DGGLSE solves in place
	cntr = 0;
	for (i=0; i<numExpansionTerms; ++i) {
	  for (it=dataPoints.begin(); it!=dataPoints.end(); ++it, ++cntr)
	    A_matrix[cntr]
	      = multivariate_polynomial(it->continuous_variables(), i);
	  C_matrix[i]
	    = multivariate_polynomial(anchorPoint.continuous_variables(), i);
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
    delete [] work;
  }

  delete [] A_matrix;

  if (err_flag) { // if numerical problems in LLS, abort with error message
    PCerr << "Error: nonzero return code from LAPACK linear least squares in "
	  << "OrthogPolyApproximation::find_coefficients()" << std::endl;
    abort_handler(-1);
  }
}


/** The coefficients of the PCE for the response are calculated using a
    Galerkin projection of the response against each multivariate orthogonal
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
  std::list<SurrogateDataPoint>::iterator it;
  size_t i, j, k, offset = 0, num_pts = dataPoints.size(),
    num_deriv_vars = expansionCoeffGrads.numRows();
  if (!anchorPoint.is_null()) {
    offset   = 1;
    num_pts += offset;
  }

  Real empty_r;
  /*
  // The following implementation evaluates all PCE coefficients
  // using a consistent expectation formulation
  for (i=0; i<numExpansionTerms; ++i) {
    Real& chaos_coeff_i = expansionCoeffs[i];
    chaos_coeff_i = (anchorPoint.is_null()) ? 0.0 :
      anchorPoint.response_function() *
      multivariate_polynomial(anchorPoint.continuous_variables(), i);
    for (j=offset, it=dataPoints.begin(); j<num_pts; ++j, ++it)
      chaos_coeff_i += it->response_function() * 
        multivariate_polynomial(it->continuous_variables(), i);
    chaos_coeff_i /= num_pts * norm_squared(i);
#ifdef DEBUG
    PCout << "coeff[" << i << "] = " << chaos_coeff_i
	  << " norm squared[" << i <<"] = " << norm_squared(i) << '\n';
#endif // DEBUG
  }
  */

  // This alternate implementation evaluates the first PCE coefficient (the
  // response mean) as an expectation and then removes the mean from the
  // expectation evaluation of all subsequent coefficients.  This approach
  // has been observed to result in better results for small sample sizes.
  Real& mean = (expansionCoeffFlag) ? expansionCoeffs[0] : empty_r;
  Real* mean_grad
    = (expansionGradFlag) ? (Real*)expansionCoeffGrads[0] : NULL;
  if (anchorPoint.is_null()) {
    if (expansionCoeffFlag) mean = 0.0;
    if (expansionGradFlag)  std::fill_n(mean_grad, num_deriv_vars, 0.);
  }
  else {
    if (expansionCoeffFlag)mean      = anchorPoint.response_function();
    if (expansionGradFlag) mean_grad = anchorPoint.response_gradient().values();
  }
  for (it=dataPoints.begin(); it!=dataPoints.end(); ++it) {
    if (expansionCoeffFlag)
      mean += it->response_function();
    if (expansionGradFlag) {
      const RealVector curr_pt_grad = it->response_gradient();
      for (j=0; j<num_deriv_vars; ++j)
	mean_grad[j] += curr_pt_grad[j];
    }
  }
  if (expansionCoeffFlag)
    mean /= num_pts;
  if (expansionGradFlag)
    for (j=0; j<num_deriv_vars; ++j)
      mean_grad[j] /= num_pts;
  Real chaos_sample;
  for (i=1; i<numExpansionTerms; ++i) {
    Real& chaos_coeff_i = (expansionCoeffFlag) ? expansionCoeffs[i] : empty_r;
    Real* chaos_grad_i  = (expansionGradFlag)  ? expansionCoeffGrads[i] : NULL;
    if (anchorPoint.is_null()) {
      if (expansionCoeffFlag) chaos_coeff_i = 0.0;
      if (expansionGradFlag)  std::fill_n(chaos_grad_i, num_deriv_vars, 0.);
    }
    else {
      chaos_sample
	= multivariate_polynomial(anchorPoint.continuous_variables(), i);
      if (expansionCoeffFlag)
	chaos_coeff_i = (anchorPoint.response_function() - mean) * chaos_sample;
      if (expansionGradFlag) {
	const RealVector anchor_grad = anchorPoint.response_gradient();
	for (j=0; j<num_deriv_vars; ++j)
	  chaos_grad_i[j] = (anchor_grad[j] - mean_grad[j]) * chaos_sample;
      }
    }
    for (j=offset, it=dataPoints.begin(); j<num_pts; ++j, ++it) {
      chaos_sample = multivariate_polynomial(it->continuous_variables(), i);
      if (expansionCoeffFlag)
	chaos_coeff_i += (it->response_function() - mean) * chaos_sample;
      if (expansionGradFlag) {
	const RealVector curr_pt_grad = it->response_gradient();
	for (k=0; k<num_deriv_vars; ++k)
	  chaos_grad_i[k] += (curr_pt_grad[k] - mean_grad[k]) * chaos_sample;
      }
    }
    Real term = num_pts * norm_squared(i);
    if (expansionCoeffFlag)
      chaos_coeff_i /= term;
    if (expansionGradFlag)
      for (j=0; j<num_deriv_vars; ++j)
	chaos_grad_i[j] /= term;
#ifdef DEBUG
    PCout << "coeff[" << i << "] = " << chaos_coeff_i
        //<< "coeff_grad[" << i <<"] = " << chaos_grad_i
	  << " norm squared[" << i <<"] = " << norm_squared(i) << '\n';
#endif // DEBUG
  }
}


const Real& OrthogPolyApproximation::get_value(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  approxValue = 0.0;
  Real Psi;
  size_t i, j;
  unsigned short order_1d;
  for (i=0; i<numExpansionTerms; ++i) {
    Psi = 1.0;
    for (j=0; j<numVars; ++j) {
      order_1d = multiIndex[i][j];
      if (order_1d)
	Psi *= polynomialBasis[j].get_value(x[j], order_1d);
    }
    approxValue += Psi*expansionCoeffs[i];
  }
  return approxValue;
}


const RealVector& OrthogPolyApproximation::
get_gradient(const RealVector& x)
{
  // this could define a default_dvv and call get_gradient(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_gradient()" << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != numVars)
    approxGradient.sizeUninitialized(numVars);
  approxGradient = 0.0;

  // sum expansion to get response gradient prediction
  size_t i, j, k;
  for (i=0; i<numExpansionTerms; ++i) {
    for (j=0; j<numVars; ++j) {
      Real term_i_grad_j = 1.0;
      for (k=0; k<numVars; ++k)
	term_i_grad_j *= (k == j) ?
	  polynomialBasis[k].get_gradient(x[k], multiIndex[i][k]) :
	  polynomialBasis[k].get_value(x[k],    multiIndex[i][k]);
      approxGradient[j] += expansionCoeffs[i]*term_i_grad_j;
    }
  }
  return approxGradient;
}


const RealVector& OrthogPolyApproximation::
get_gradient(const RealVector& x, const UIntArray& dvv)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_gradient()" << std::endl;
    abort_handler(-1);
  }

  size_t num_deriv_vars = dvv.size();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.0;

  // sum expansion to get response gradient prediction
  size_t i, j, k, deriv_index;
  for (i=0; i<numExpansionTerms; ++i) {
    for (j=0; j<num_deriv_vars; ++j) {
      deriv_index = dvv[j] - 1; // *** requires an "All" view
      Real term_i_grad_j = 1.0;
      for (k=0; k<numVars; ++k)
	term_i_grad_j *= (k == deriv_index) ?
	  polynomialBasis[k].get_gradient(x[k], multiIndex[i][k]) :
	  polynomialBasis[k].get_value(x[k],    multiIndex[i][k]);
      approxGradient[j] += expansionCoeffs[i]*term_i_grad_j;
    }
  }
  return approxGradient;
}


/** In this case, all expansion variables are random variables and the
    mean of the expansion is simply the first chaos coefficient. */
const Real& OrthogPolyApproximation::get_mean()
{
  // Error check for required data
  if (!expansionCoeffFlag) {
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
  if (!expansionCoeffFlag) {
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
      PCout << "Mean estimate inclusion: term index = " << i << " Psi = " << Psi
	    << " PCE coeff = " << expansionCoeffs[i] << " total = "
	    << expansionMean << '\n';
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
  if (!expansionGradFlag) {
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
get_mean_gradient(const RealVector& x, const UIntArray& dvv)
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
      if (!expansionGradFlag) {
	PCerr << "Error: expansion coefficient gradients not defined in "
	      << "OrthogPolyApproximation::get_mean_gradient()." << std::endl;
	abort_handler(-1);
      }
      expansionMeanGrad[i] = expansionCoeffGrads[0][cntr];
    }
    else if (!expansionCoeffFlag) { // Error check for required data
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


/** In this case, all expansion variables are random variables and the
    variance of the expansion is the sum over all but the first term
    of the coefficients squared times the polynomial norms squared. */
const Real& OrthogPolyApproximation::get_variance()
{
  // Error check for required data
  if (!expansionCoeffFlag) {
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
    expansionVariance += expansionCoeffs[i] * exp_coeffs_2[i] * norm_squared(i);
  return expansionVariance;
}


/** In this case, a subset of the expansion variables are random
    variables and the variance of the expansion involves summations
    over this subset. */
const Real& OrthogPolyApproximation::get_variance(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
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
      Real norm_sq_i = norm_squared_random(i);
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
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::get_variance_gradient()" << std::endl;
    abort_handler(-1);
  }
  if (!expansionGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in "
	  << "OrthogPolyApproximation::get_variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  size_t i, j, num_deriv_vars = expansionCoeffGrads.numRows();
  if (expansionVarianceGrad.length() != num_deriv_vars)
    expansionVarianceGrad.sizeUninitialized(num_deriv_vars);
  expansionVarianceGrad = 0.0;
  for (i=1; i<numExpansionTerms; ++i) {
    Real term_i = 2. * expansionCoeffs[i] * norm_squared(i);
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
get_variance_gradient(const RealVector& x, const UIntArray& dvv)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
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
    if (randomVarsKey[deriv_index] && !expansionGradFlag) { // Error check 
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
	Real norm_sq_j = norm_squared_random(j);
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


const Real& OrthogPolyApproximation::norm_squared(size_t expansion_index)
{
  // the norm squared of a particular multivariate polynomial is the product of
  // the norms squared of the numVars univariate polynomials that comprise it.
  multiPolyNormSq = 1.0;
  unsigned short order;
  for (size_t i=0; i<numVars; ++i) {
    order = multiIndex[expansion_index][i];
    if (order)
      multiPolyNormSq *= polynomialBasis[i].norm_squared(order);
  }
  return multiPolyNormSq;
}


const Real& OrthogPolyApproximation::norm_squared_random(size_t expansion_index)
{
  // the norm squared of a particular multivariate polynomial is the product of
  // the norms squared of the numVars univariate polynomials that comprise it.
  multiPolyNormSq = 1.0;
  unsigned short order;
  for (size_t i=0; i<numVars; ++i) {
    if (randomVarsKey[i]) {
      order = multiIndex[expansion_index][i];
      if (order)
	multiPolyNormSq *= polynomialBasis[i].norm_squared(order);
    }
  }
  return multiPolyNormSq;
}


void OrthogPolyApproximation::compute_global_sensitivity()
{
  // Allocation of memory for sensitivity variables done in 
  // PolynomialApproximation::allocate_arrays()

  // sobolIndices are index via binary number represenation
  // e.g. if there are 5 variables -> 5 bit represenation
  // if a variable belongs to a subset \alpha, it has a value of 1; otherwise 0
  // | 0 | 0 | 0 | 0 | 1 | represents the subset that only contains variable 1
  // sobolIndices[0] is set to 1; to calculate it is redundant

  // Index the set using binary number representation
  size_t index_bin, i, j;
  // Compute a pseudo-variance that makes no distinction between probabilistic
  // variables and non-probabilistic variables
  Real p_var = 0;
  for (i=1; i<numExpansionTerms; i++) 
    p_var += norm_squared(i)*expansionCoeffs(i)*expansionCoeffs(i);

  // iterate through multiIndex and store sensitivities
  sobolIndices      = 0.;
  sobolIndices[0]   = 1.;
  totalSobolIndices = 0.;
  for (i=1; i<numExpansionTerms; ++i) {
    index_bin = 0;
    Real p_var_i = norm_squared(i)*expansionCoeffs(i)*expansionCoeffs(i)/p_var;
    for (j=0; j<numVars; ++j) {
      if (multiIndex[i][j]) {
	// convert this subset multiIndex[i] into binary number
	index_bin += (size_t)pow(2.,(int)j);
	// for any constituent variable j in exansion term i, the expansion
	// term contributes to the total sensitivty of variable j
	totalSobolIndices[j] += p_var_i;
      }
    }
    sobolIndices[index_bin] += p_var_i;
  }
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
