/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 SparseGridDriver
//- Description: Implementation code for SparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "SparseGridDriver.hpp"
#include "sandia_rules.H"
#include "sparse_grid_mixed_growth.H"
#include "sgmga.H"
#include "OrthogPolyApproximation.hpp"
#include "NumericGenOrthogPolynomial.hpp"
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: SparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void SparseGridDriver::dimension_preference(const RealVector& dim_pref)
{
  RealVector aniso_wts;
  if (!dim_pref.empty()) {
    size_t num_pref = dim_pref.length();
    aniso_wts.sizeUninitialized(num_pref);
    webbur::sgmga_importance_to_aniso(num_pref, dim_pref.values(),
				      aniso_wts.values());
#ifdef DEBUG
    PCout << "dimension preference:\n"; write_data(PCout, dim_pref);
    PCout << "anisotropic weights after sgmga_importance_to_aniso():\n";
    write_data(PCout, aniso_wts);
#endif
  }
  anisotropic_weights(aniso_wts);
}


void SparseGridDriver::anisotropic_weights(const RealVector& aniso_wts)
{
  if (aniso_wts.empty())
    dimIsotropic = true;
  else {
    if (aniso_wts.length() != numVars) {
      PCerr << "Error: length of sparse grid anisotropic weights "
	    << "specification is inconsistent with\n       number of variables "
	    << "in SparseGridDriver::anisotropic_weights()." << std::endl;
      abort_handler(-1);
    }

    dimIsotropic = true;
    anisoLevelWts.resize(numVars);
    // truncate any negative values
    size_t i;
    for (i=0; i<numVars; ++i)
      anisoLevelWts[i] = (aniso_wts[i] < 0.) ? 0. : aniso_wts[i];
    // detect anisotropy
    Real wt0 = anisoLevelWts[0];
    for (i=1; i<numVars; ++i)
      if (std::abs(anisoLevelWts[i] - wt0) > DBL_EPSILON)
	{ dimIsotropic = false; break; }
    // normalize and enforce axis lower bounds/weight upper bounds
    if (!dimIsotropic) {
      int option = 1; // weights scaled so that minimum nonzero entry is 1
      webbur::sgmga_aniso_normalize(option, numVars, anisoLevelWts.values());
#ifdef DEBUG
      PCout << "anisoLevelWts after sgmga_aniso_normalize():\n";
      write_data(PCout, anisoLevelWts);
#endif
      // enforce axis lower bounds, if present, for current ssgLevel.  An axis
      // lower bound defines a weight upper bound based on the current ssgLevel:
      // LB_i = level*wt_min/wt_i --> wt_i = level*wt_min/LB_i and wt_min=1.
      // Catch special case of dim_pref_i = 0 --> wt_i = LB_i = 0.
      if (!axisLowerBounds.empty()) {
	for (i=0; i<numVars; ++i)
	  if (axisLowerBounds[i] > 1.e-10) {                 // nonzero LB
	    Real wt_u_bnd = (Real)ssgLevel/axisLowerBounds[i];
	    anisoLevelWts[i] = (anisoLevelWts[i] > 1.e-10) ? // nonzero wt
	      std::min(wt_u_bnd, anisoLevelWts[i]) : wt_u_bnd;
	  }
#ifdef DEBUG
	PCout << "anisoLevelWts after axisLowerBounds enforcement:\n";
	write_data(PCout, anisoLevelWts);
#endif
      }
    }
  }
}


void SparseGridDriver::update_axis_lower_bounds()
{
  if (axisLowerBounds.empty())
    axisLowerBounds.sizeUninitialized(numVars);
  // An axisLowerBound is the maximum index coverage achieved on a coordinate
  // axis (when all other indices are zero); it defines a constraint for
  // minimum coordinate coverage in future refinements.  The linear index set
  // constraint is level*wt_min-|wt| < j.wt <= level*wt_min, which becomes
  // level-|wt| < j_i w_i <= level for wt_min=1 and all other indices=0.
  // The max feasible j_i is then level/w_i (except for special case w_i=0).
  if (dimIsotropic)
    axisLowerBounds = (Real)ssgLevel; // all weights = 1
  else // min nonzero weight scaled to 1 --> just catch special case w_i=0
    for (size_t i=0; i<numVars; ++i)
      axisLowerBounds[i] = (anisoLevelWts[i] > 1.e-10) ? // nonzero wt
	(Real)ssgLevel/anisoLevelWts[i] : 0.;
}


void SparseGridDriver::
initialize_grid(const ShortArray& u_types,  size_t ssg_level,
		const RealVector& dim_pref, const String& ssg_usage,
		short exp_growth, short nested_uniform_rule)
{
  numVars = u_types.size();
  sparseGridUsage = ssg_usage;
  level(ssg_level);
  dimension_preference(dim_pref);

  integrationRules.resize(numVars);
  growthRules.resize(numVars);
  compute1DPoints.resize(numVars);
  compute1DWeights.resize(numVars);

  // For STANDARD exponential growth, use of nested rules is restricted to
  // isotropic uniform in order to enforce consistent growth rates:
  bool nested_rules = true;
  if (exp_growth == UNRESTRICTED_GROWTH)
    for (size_t i=0; i<numVars; ++i)
      if (u_types[i] != STD_UNIFORM)
	{ nested_rules = false; break; }
  // For MODERATE exponential growth, nested rules can be used heterogeneously
  // and synchronized with STANDARD Gaussian linear growth.
  // For SLOW exponential growth, nested rules can be used heterogeneously
  // and synchronized with SLOW Gaussian linear growth (not yet available).
  for (size_t i=0; i<numVars; i++) {
    switch (u_types[i]) {
    case STD_NORMAL:
      compute1DPoints[i]  = webbur::hermite_compute_points_np;
      compute1DWeights[i] = webbur::hermite_compute_weights_np;
      integrationRules[i] = GAUSS_HERMITE;
      growthRules[i]      = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR_ODD : MODERATE_LINEAR; break;
    case BOUNDED_NORMAL:
      compute1DPoints[i]  = bounded_normal_gauss_points;
      compute1DWeights[i] = bounded_normal_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH;
      growthRules[i]      = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    case LOGNORMAL:
      compute1DPoints[i]  = lognormal_gauss_points;
      compute1DWeights[i] = lognormal_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH;
      growthRules[i]      = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    case BOUNDED_LOGNORMAL:
      compute1DPoints[i]  = bounded_lognormal_gauss_points;
      compute1DWeights[i] = bounded_lognormal_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH;
      growthRules[i]      = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    case STD_UNIFORM:
      // For tensor-product quadrature, Gauss-Legendre is used due to greater
      // polynomial exactness since nesting is not a concern.  For nested sparse
      // grids, Clenshaw-Curtis or Gauss-Patterson can be better selections.
      // However, sparse grids that are isotropic in level but anisotropic in
      // rule become skewed when mixing Gauss rules with CC.  For this reason,
      // CC is selected only if isotropic in rule (for now).
      if (nested_rules) {
	integrationRules[i] = nested_uniform_rule;
	switch (nested_uniform_rule) {
	case GAUSS_PATTERSON: // closed fully nested
	  compute1DPoints[i]  = webbur::patterson_lookup_points_np;
	  compute1DWeights[i] = webbur::patterson_lookup_weights_np; break;
	case CLENSHAW_CURTIS:
	  compute1DPoints[i]  = webbur::clenshaw_curtis_compute_points_np;
	  compute1DWeights[i] = webbur::clenshaw_curtis_compute_weights_np;
	  break;
	case FEJER2:
	  compute1DPoints[i]  = webbur::fejer2_compute_points_np;
	  compute1DWeights[i] = webbur::fejer2_compute_weights_np;   break;
	}
	switch (exp_growth) {
	case SLOW_RESTRICTED_GROWTH:
	  growthRules[i] = SLOW_EXPONENTIAL;     break;
	case MODERATE_RESTRICTED_GROWTH:
	  growthRules[i] = MODERATE_EXPONENTIAL; break;
	case UNRESTRICTED_GROWTH:
	  growthRules[i] = FULL_EXPONENTIAL;     break;
	} 
      }
      else {
	compute1DPoints[i]  = webbur::legendre_compute_points_np;
	compute1DWeights[i] = webbur::legendre_compute_weights_np;
	integrationRules[i] = GAUSS_LEGENDRE;
	growthRules[i]      = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	  SLOW_LINEAR_ODD : MODERATE_LINEAR;
      }
      break;
    case LOGUNIFORM:
      compute1DPoints[i]  = loguniform_gauss_points;
      compute1DWeights[i] = loguniform_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH;
      growthRules[i] = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    case TRIANGULAR:
      compute1DPoints[i]  = triangular_gauss_points;
      compute1DWeights[i] = triangular_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH;
      growthRules[i] = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    case STD_EXPONENTIAL:
      compute1DPoints[i]  = webbur::laguerre_compute_points_np;
      compute1DWeights[i] = webbur::laguerre_compute_weights_np;
      integrationRules[i] = GAUSS_LAGUERRE;
      growthRules[i] = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    case STD_BETA:
      compute1DPoints[i]  = webbur::jacobi_compute_points_np;
      compute1DWeights[i] = webbur::jacobi_compute_weights_np;
      integrationRules[i] = GAUSS_JACOBI;
      growthRules[i] = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    case STD_GAMMA:
      compute1DPoints[i]  = webbur::gen_laguerre_compute_points_np;
      compute1DWeights[i] = webbur::gen_laguerre_compute_weights_np;
      integrationRules[i] = GEN_GAUSS_LAGUERRE;
      growthRules[i] = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    case GUMBEL:
      compute1DPoints[i]  = gumbel_gauss_points;
      compute1DWeights[i] = gumbel_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH;
      growthRules[i] = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    case FRECHET:
      compute1DPoints[i]  = frechet_gauss_points;
      compute1DWeights[i] = frechet_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH;
      growthRules[i] = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    case WEIBULL:
      compute1DPoints[i]  = weibull_gauss_points;
      compute1DWeights[i] = weibull_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH;
      growthRules[i] = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    case HISTOGRAM_BIN: {
      compute1DPoints[i]  = histogram_bin_gauss_points;
      compute1DWeights[i] = histogram_bin_gauss_weights;
      integrationRules[i] = GOLUB_WELSCH;
      growthRules[i] = (exp_growth == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    }
    default:
      PCerr << "Error: unsupported distribution type in SparseGridDriver for "
	    << "u_type = " << u_types[i] << std::endl;
      abort_handler(-1);
      break;
    }
  }

  // avoid regenerating polynomialBasis, if up to date
  ShortArray basis_types, gauss_modes;
  OrthogPolyApproximation::distribution_types(u_types, integrationRules,
					      basis_types, gauss_modes);
  OrthogPolyApproximation::distribution_basis(basis_types, gauss_modes,
					      polynomialBasis);
}


void SparseGridDriver::
initialize_grid_parameters(const ShortArray& u_types,
			   const DistributionParams& dp)
{
  numPolyParams.assign(numVars, 0);
  size_t i, num_total_params = 0, hbuv_cntr = 0;
  for (i=0; i<numVars; i++) {
    switch (u_types[i]) {
    case STD_GAMMA:
      numPolyParams[i] = 1; break;
    case STD_BETA: case LOGNORMAL: case LOGUNIFORM:
    case GUMBEL:   case FRECHET:   case WEIBULL:
      numPolyParams[i] = 2; break;
    case TRIANGULAR:
      numPolyParams[i] = 3; break;
    case BOUNDED_NORMAL: case BOUNDED_LOGNORMAL:
      numPolyParams[i] = 4; break;
    case HISTOGRAM_BIN:
      numPolyParams[i] = dp.histogram_bin_pairs(hbuv_cntr).length();
      ++hbuv_cntr;          break;
    }
    num_total_params += numPolyParams[i];
  }
  polyParams.resize(num_total_params); // may be zero size

  size_t nuv_cntr = 0, lnuv_cntr = 0, luuv_cntr = 0, tuv_cntr = 0, buv_cntr = 0,
    gauv_cntr = 0, guuv_cntr = 0, fuv_cntr = 0, wuv_cntr = 0, pp_cntr = 0;
  hbuv_cntr = 0;
  for (i=0; i<numVars; i++) {
    switch (u_types[i]) {
    case STD_NORMAL:
      ++nuv_cntr; break;
    case BOUNDED_NORMAL:
      polyParams[pp_cntr] = dp.normal_mean(nuv_cntr);          ++pp_cntr;
      polyParams[pp_cntr] = dp.normal_std_deviation(nuv_cntr); ++pp_cntr;
      polyParams[pp_cntr] = dp.normal_lower_bound(nuv_cntr);   ++pp_cntr;
      polyParams[pp_cntr] = dp.normal_upper_bound(nuv_cntr);   ++pp_cntr;
      ++nuv_cntr; break;
    case LOGNORMAL:
      moments_from_lognormal_spec(dp.lognormal_means(),
				  dp.lognormal_std_deviations(),
				  dp.lognormal_lambdas(), dp.lognormal_zetas(),
				  dp.lognormal_error_factors(), lnuv_cntr,
				  polyParams[pp_cntr], polyParams[pp_cntr+1]);
      pp_cntr += 2;
      ++lnuv_cntr; break;
    case BOUNDED_LOGNORMAL:
      moments_from_lognormal_spec(dp.lognormal_means(),
				  dp.lognormal_std_deviations(),
				  dp.lognormal_lambdas(), dp.lognormal_zetas(),
				  dp.lognormal_error_factors(), lnuv_cntr,
				  polyParams[pp_cntr], polyParams[pp_cntr+1]);
      pp_cntr += 2;
      polyParams[pp_cntr] = dp.lognormal_lower_bound(lnuv_cntr); ++pp_cntr;
      polyParams[pp_cntr] = dp.lognormal_upper_bound(lnuv_cntr); ++pp_cntr;
      ++lnuv_cntr; break;
    case STD_UNIFORM:
      break;
    case LOGUNIFORM:
      polyParams[pp_cntr] = dp.loguniform_lower_bound(luuv_cntr); ++pp_cntr;
      polyParams[pp_cntr] = dp.loguniform_upper_bound(luuv_cntr); ++pp_cntr;
      ++luuv_cntr; break;
    case TRIANGULAR:
      polyParams[pp_cntr] = dp.triangular_mode(tuv_cntr);        ++pp_cntr;
      polyParams[pp_cntr] = dp.triangular_lower_bound(tuv_cntr); ++pp_cntr;
      polyParams[pp_cntr] = dp.triangular_upper_bound(tuv_cntr); ++pp_cntr;
      ++tuv_cntr; break;
    case STD_EXPONENTIAL:
      break;
    case STD_BETA:
      // convert stat to poly and assume consistent ordering via b_cntr:
      polyParams[pp_cntr] = dp.beta_beta(buv_cntr)  - 1.; ++pp_cntr;
      polyParams[pp_cntr] = dp.beta_alpha(buv_cntr) - 1.; ++pp_cntr;
      ++buv_cntr; break;
    case STD_GAMMA:
      // convert stat to poly and assume consistent ordering via g_cntr:
      polyParams[pp_cntr] = dp.gamma_alpha(gauv_cntr) - 1.; ++pp_cntr;
      ++gauv_cntr; break;
    case GUMBEL:
      polyParams[pp_cntr] = dp.gumbel_alpha(guuv_cntr); ++pp_cntr;
      polyParams[pp_cntr] = dp.gumbel_beta(guuv_cntr);  ++pp_cntr;
      ++guuv_cntr; break;
    case FRECHET:
      polyParams[pp_cntr] = dp.frechet_alpha(fuv_cntr); ++pp_cntr;
      polyParams[pp_cntr] = dp.frechet_beta(fuv_cntr);  ++pp_cntr;
      ++fuv_cntr; break;
    case WEIBULL:
      polyParams[pp_cntr] = dp.weibull_alpha(wuv_cntr); ++pp_cntr;
      polyParams[pp_cntr] = dp.weibull_beta(wuv_cntr);  ++pp_cntr;
      ++wuv_cntr; break;
    case HISTOGRAM_BIN: {
      const RealVector& hbuv_prs_cntr = dp.histogram_bin_pairs(hbuv_cntr);
      size_t j, num_params = hbuv_prs_cntr.length();
      for (j=0; j<num_params; ++j, ++pp_cntr)
	polyParams[pp_cntr] = hbuv_prs_cntr[j];
      ++hbuv_cntr; break;
    }
    default:
      PCerr << "Error: unsupported distribution type in SparseGridDriver for "
	    << "u_type = " << u_types[i] << std::endl;
      abort_handler(-1);
      break;
    }
  }

  OrthogPolyApproximation::distribution_parameters(u_types, dp,
						   polynomialBasis);
}


int SparseGridDriver::grid_size()
{
  // polyParams may be zero size
  double* poly_params = (polyParams.size()) ? &polyParams[0] : NULL;
  return (dimIsotropic) ?
    webbur::sparse_grid_mixed_growth_size(numVars, ssgLevel,
      &integrationRules[0], &numPolyParams[0], poly_params,
      &compute1DPoints[0], duplicateTol, &growthRules[0]) :
    webbur::sgmga_size(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &numPolyParams[0], poly_params,
      &compute1DPoints[0], duplicateTol, &growthRules[0]);
}


int SparseGridDriver::grid_size_total()
{
  return (dimIsotropic) ?
    webbur::sparse_grid_mixed_growth_size_total(numVars, ssgLevel,
      &integrationRules[0], &growthRules[0]) :
    webbur::sgmga_size_total(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &growthRules[0]);
}


void SparseGridDriver::compute_grid()
{
  // --------------------------------
  // Get number of collocation points
  // --------------------------------
  int num_total_pts = grid_size_total(), num_colloc_pts = grid_size();
#ifdef DEBUG
  PCout << "Total number of sparse grid integration points: "
	<< num_colloc_pts << '\n';
#endif // DEBUG

  // ----------------------------------------------
  // Get collocation points and integration weights
  // ----------------------------------------------
  weightSets.sizeUninitialized(num_colloc_pts);
  variableSets.shapeUninitialized(numVars, num_colloc_pts);// Teuchos: col major
  uniqueIndexMapping.resize(num_total_pts);

  int* sparse_order  = new int [num_colloc_pts*numVars];
  int* sparse_index  = new int [num_colloc_pts*numVars];
  // polyParams may be zero size
  double* poly_params = (polyParams.size()) ? &polyParams[0] : NULL;
  if (dimIsotropic) {
    webbur::sparse_grid_mixed_growth_unique_index(numVars, ssgLevel,
      &integrationRules[0], &numPolyParams[0], poly_params,
      &compute1DPoints[0], duplicateTol, num_colloc_pts, num_total_pts,
      &growthRules[0], &uniqueIndexMapping[0]);
    webbur::sparse_grid_mixed_growth_index(numVars, ssgLevel,
      &integrationRules[0], num_colloc_pts, num_total_pts,
      &uniqueIndexMapping[0], &growthRules[0], sparse_order, sparse_index);
    webbur::sparse_grid_mixed_growth_weight(numVars, ssgLevel,
      &integrationRules[0], &numPolyParams[0], poly_params,
      &compute1DWeights[0], num_colloc_pts, num_total_pts,
      &uniqueIndexMapping[0], &growthRules[0], weightSets.values());
    webbur::sparse_grid_mixed_growth_point(numVars, ssgLevel,
      &integrationRules[0], &numPolyParams[0], poly_params,
      &compute1DPoints[0], num_colloc_pts, sparse_order, sparse_index,
      &growthRules[0], variableSets.values());
  }
  else {
    webbur::sgmga_unique_index(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &numPolyParams[0], poly_params,
      &compute1DPoints[0], duplicateTol, num_colloc_pts, num_total_pts,
      &growthRules[0], &uniqueIndexMapping[0]);
    webbur::sgmga_index(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], num_colloc_pts, num_total_pts,
      &uniqueIndexMapping[0], &growthRules[0], sparse_order, sparse_index);
    webbur::sgmga_weight(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &numPolyParams[0], poly_params,
      &compute1DWeights[0], num_colloc_pts, num_total_pts,
      &uniqueIndexMapping[0], &growthRules[0], weightSets.values());
    webbur::sgmga_point(numVars, anisoLevelWts.values(), ssgLevel,
      &integrationRules[0], &numPolyParams[0], poly_params,
      &compute1DPoints[0], num_colloc_pts, sparse_order, sparse_index,
      &growthRules[0], variableSets.values());
  }
  delete [] sparse_order;
  delete [] sparse_index;
#ifdef DEBUG
  PCout << "uniqueIndexMapping:\n" << uniqueIndexMapping << '\n';
#endif

  // compute scale factors for VPISparseGrid outputs
  Real wt_factor = 1.;
  RealVector pt_factor(numVars, false); pt_factor = 1.;
  size_t i, j;
  BasisPolynomial chebyshev_poly;
  for (i=0; i<numVars; ++i) {
    switch (integrationRules[i]) {
    case GAUSS_HERMITE:
      pt_factor[i] = polynomialBasis[i].point_factor();
      wt_factor   *= polynomialBasis[i].weight_factor(); break;
    case GAUSS_LEGENDRE:  case GAUSS_PATTERSON:
    case GAUSS_JACOBI:    case GEN_GAUSS_LAGUERRE:
      wt_factor *= polynomialBasis[i].weight_factor();   break;
    case CLENSHAW_CURTIS: case FEJER2:
      if (chebyshev_poly.is_null())
	chebyshev_poly = BasisPolynomial(CHEBYSHEV);//, mode);
      //chebyshev_poly.gauss_mode(integrationRules[i]);// no effect on wtFactor
      wt_factor *= chebyshev_poly.weight_factor();       break;
    //case GAUSS_LAGUERRE: case GOLUB_WELSCH: // scaling is OK
    }
  }

  // apply point/weight scaling
  weightSets.scale(wt_factor);
  for (i=0; i<num_colloc_pts; ++i)
    for (j=0; j<numVars; ++j)
      variableSets[i][j] *= pt_factor[j]; // ith column, jth row
#ifdef DEBUG
  PCout << "\nPoint factors = " << pt_factor
	<< "\nWeight factor = "	<< wt_factor << std::endl;
#endif // DEBUG

  if (sparseGridUsage == "interpolation") { // 1D arrays not needed for PCE
    // ----------------------------
    // Define 1-D point/weight sets
    // ----------------------------
    if (gaussPts1D.empty())
      gaussPts1D.resize(numVars);
    if (gaussWts1D.empty())
      gaussWts1D.resize(numVars);
    // level_index (j indexing) range is 0:w, level (i indexing) range is 1:w+1
    unsigned short level_index, order;
    for (i=0; i<numVars; i++) {
      gaussPts1D[i].resize(ssgLevel + 1); gaussWts1D[i].resize(ssgLevel + 1);
      switch (integrationRules[i]) {
      case CLENSHAW_CURTIS: case FEJER2:
	chebyshev_poly.gauss_mode(integrationRules[i]); // integration mode
	for (level_index=0; level_index<=ssgLevel; level_index++) {
	  level_to_order(i, level_index, order);
	  gaussPts1D[i][level_index] = chebyshev_poly.gauss_points(order);
	  gaussWts1D[i][level_index] = chebyshev_poly.gauss_weights(order);
	}
	break;
      default: // Gaussian rules
	for (level_index=0; level_index<=ssgLevel; level_index++) {
	  level_to_order(i, level_index, order);
	  gaussPts1D[i][level_index] = polynomialBasis[i].gauss_points(order);
	  gaussWts1D[i][level_index] = polynomialBasis[i].gauss_weights(order);
	}
	break;
      }
    }
  }

  /*
  // -----------------------------------
  // Get sparse grid index/base mappings
  // -----------------------------------
  size_t size = numVars*num_colloc_pts, cntr = 0;
  int* indices = new int [size];
  int* bases   = new int [size];

  webbur::sparse_grid_mixed_growth_index(numVars, ssgLevel,
    integrationRules, num_colloc_pts, num_total_pts, uniqueIndexMapping,
    &growthRules[0], bases, indices);

  IntArray key(2*numVars);
  unsigned short closed_order_max;
  level_to_order_closed_exponential(ssgLevel, closed_order_max);
  for (i=0; i<num_colloc_pts; i++) {
    for (j=0; j<numVars; j++, cntr++) {
      switch (integrationRules[j]) {
      case GAUSS_HERMITE: case GAUSS_LEGENDRE:
	key[j] = 2 * bases[cntr] + 1;                 // map to quad order
	key[j+numVars] = indices[cntr] + bases[cntr]; // 0-based index
	break;
      case CLENSHAW_CURTIS:
	key[j] = closed_order_max;      // promotion to highest grid
	key[j+numVars] = indices[cntr]; // already 0-based
	break;
      case GAUSS_LAGUERRE:
	key[j] = bases[cntr];               // already quad order
	key[j+numVars] = indices[cntr] - 1; // map to 0-based
	break;
      }
    }
    ssgIndexMap[key] = i;
#ifdef DEBUG
    PCout << "i = " << i << " key =\n" << key << std::endl;
#endif // DEBUG
  }
  delete [] indices;
  delete [] bases;
  */
}


void SparseGridDriver::
anisotropic_multi_index(Int2DArray& multi_index, RealArray& coeffs) const
{
  multi_index.clear();
  coeffs.clear();
  // Utilize webbur::sgmga_vcn_{ordered,coef} for 0-based index sets
  // (w*alpha_min-|alpha| < |alpha . j| <= w*alpha_min).
  // With scaling alpha_min = 1: w-|alpha| < |alpha . j| <= w.
  // In the isotropic case, reduces to w-N < |j| <= w, which is the same as
  // w-N+1 <= |j| <= w.
  IntArray x(numVars), x_max(numVars); //x_max = ssgLevel;
  Real wt_sum = 0., q_max = ssgLevel;
  for (size_t i=0; i<numVars; ++i) {
    const Real& wt_i = anisoLevelWts[i];
    wt_sum += wt_i;
    // minimum nonzero weight is scaled to 1, so just catch special case of 0
    x_max[i] = (wt_i > 1.e-10) ? (int)std::ceil(q_max/wt_i) : 0;
  }
  Real q_min = ssgLevel - wt_sum;
#ifdef DEBUG
  PCout << "q_min = " << q_min << " q_max = " << q_max;
#endif // DEBUG

  bool more = false;
  webbur::sgmga_vcn_ordered(numVars, anisoLevelWts.values(), &x_max[0], &x[0],
			    q_min, q_max, &more);
  while (more) {
    Real coeff = webbur::sgmga_vcn_coef(numVars, anisoLevelWts.values(),
					&x_max[0], &x[0], q_min, q_max);
    if (std::abs(coeff) > 1.e-10) {
      multi_index.push_back(x);
      coeffs.push_back(coeff);
    }
    webbur::sgmga_vcn_ordered(numVars, anisoLevelWts.values(), &x_max[0], &x[0],
			      q_min, q_max, &more);
  }
}


// TO DO: avoid recalculating existing Gauss pts within following functions
//
// Burkardt: instead of params/data, pass index and retain array of NGPoly
//
// Then add interface to pass in array of polynomials in initialize() instead
// of current initialize_grid_parameters().
//
// KDE: look at figTree --> if lightweight; else, code simple box kernel and
// Gaussian kernel
//
// add/activate STOCHASTIC_EXPANSION allowing moments or KDE
void SparseGridDriver::
bounded_normal_gauss_points(int order, int num_params, double* params,
			    double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::bounded_normal_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_normal_distribution(params[0], params[1], params[2], params[3]);
                                //(mean,      stdev,     lwr,       upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
bounded_normal_gauss_weights(int order, int num_params, double* params,
			     double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::bounded_normal_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_normal_distribution(params[0], params[1], params[2], params[3]);
                                //(mean,      stdev,     lwr,       upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
lognormal_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::lognormal_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.lognormal_distribution(params[0], params[1]); //(mean, stdev);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
lognormal_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::lognormal_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.lognormal_distribution(params[0], params[1]); //(mean, stdev);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
bounded_lognormal_gauss_points(int order, int num_params, double* params,
			       double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::bounded_lognormal_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_lognormal_distribution(params[0], params[1], params[2],
				      params[3]); //(mean, stdev, lwr, upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
bounded_lognormal_gauss_weights(int order, int num_params, double* params,
				double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::bounded_lognormal_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_lognormal_distribution(params[0], params[1], params[2],
				      params[3]); //(mean, stdev, lwr, upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
loguniform_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::loguniform_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.loguniform_distribution(params[0], params[1]); //(lwr, upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
loguniform_gauss_weights(int order, int num_params, double* params,
			 double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::loguniform_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.loguniform_distribution(params[0], params[1]); //(lwr, upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
triangular_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 3) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::triangular_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.triangular_distribution(params[0], params[1], params[2]);
                            //(mode,      lwr,       upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
triangular_gauss_weights(int order, int num_params, double* params,
			 double* data)
{
  if (num_params != 3) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::triangular_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.triangular_distribution(params[0], params[1], params[2]);
                            //(mode,      lwr,       upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
gumbel_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::gumbel_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.gumbel_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
gumbel_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::gumbel_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.gumbel_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
frechet_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::frechet_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.frechet_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
frechet_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::frechet_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.frechet_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
weibull_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::weibull_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.weibull_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
weibull_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::weibull_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.weibull_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGridDriver::
histogram_bin_gauss_points(int order, int num_params, double* params,
			   double* data)
{
  if (num_params < 4 || num_params % 2) { // need at least 2 (x,c) pairs
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::histogram_bin_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  RealVector param_vec;
  copy_data(params, num_params, param_vec);
  ngop.histogram_bin_distribution(param_vec);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGridDriver::
histogram_bin_gauss_weights(int order, int num_params, double* params,
			    double* data)
{
  if (num_params < 4 || num_params % 2) { // need at least 2 (x,c) pairs
    PCerr << "Error: wrong number of distribution parameters in "
	  << "SparseGridDriver::histogram_bin_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  RealVector param_vec;
  copy_data(params, num_params, param_vec);
  ngop.histogram_bin_distribution(param_vec);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}

} // namespace Pecos
