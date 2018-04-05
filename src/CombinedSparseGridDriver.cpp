/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 CombinedSparseGridDriver
//- Description: Implementation code for CombinedSparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "CombinedSparseGridDriver.hpp"
#include "SharedPolyApproxData.hpp"
#include "sandia_sgmg.hpp"
#include "sandia_sgmga.hpp"
#include "sandia_sgmgg.hpp"
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: CombinedSparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {

/// initialize static member pointer to active driver instance
CombinedSparseGridDriver* CombinedSparseGridDriver::sgdInstance(NULL);


void CombinedSparseGridDriver::
initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
		const ShortArray& u_types,
		const ExpansionConfigOptions& ec_options,
		BasisConfigOptions& bc_options, short growth_rate,
		bool track_colloc, bool track_uniq_prod_wts)
{
  SparseGridDriver::initialize_grid(ssg_level, dim_pref, u_types, ec_options,
				    bc_options, growth_rate);
  trackCollocDetails     = track_colloc;
  trackUniqueProdWeights = track_uniq_prod_wts;

  // set a rule-dependent duplicateTol
  initialize_duplicate_tolerance();
  // set compute1D{Points,Type1Weights,Type2Weights}
  initialize_rule_pointers();
  // set levelGrowthToOrder
  initialize_growth_pointers();
}


void CombinedSparseGridDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis)
{
  IntegrationDriver::initialize_grid(poly_basis);

  // set a rule-dependent duplicateTol
  initialize_duplicate_tolerance();
  // set compute1D{Points,Type1Weights,Type2Weights}
  initialize_rule_pointers();
  // set levelGrowthToOrder
  initialize_growth_pointers();
}


void CombinedSparseGridDriver::clear_inactive()
{
  //SparseGridDriver::clear_inactive();

  std::map<UShortArray, UShort2DArray>::iterator sm_it
    = smolyakMultiIndex.begin();
  std::map<UShortArray, IntArray>::iterator  sc_it = smolyakCoeffs.begin();
  std::map<UShortArray, UShort3DArray>::iterator ck_it = collocKey.begin();
  std::map<UShortArray, Sizet2DArray>::iterator  ci_it = collocIndices.begin();
  std::map<UShortArray, RealVector>::iterator t1_it = type1WeightSets.begin();
  std::map<UShortArray, RealMatrix>::iterator t2_it = type2WeightSets.begin();

  while (sm_it != smolyakMultiIndex.end())
    if (sm_it == smolMIIter) { // preserve active
      ++sm_it; ++sc_it; ++ck_it; ++ci_it;
      if (trackUniqueProdWeights)
	{ ++t1_it; if (computeType2Weights) ++t2_it; }
    }
    else { // clear inactive: postfix increments manage iterator invalidations
      smolyakMultiIndex.erase(sm_it++);   smolyakCoeffs.erase(sc_it++);
      collocKey.erase(ck_it++);           collocIndices.erase(ci_it++);
      if (trackUniqueProdWeights) {
	type1WeightSets.erase(t1_it++);
	if (computeType2Weights)          type2WeightSets.erase(t2_it++);
      }
    }
}


const UShortArray& CombinedSparseGridDriver::maximal_grid() const
{
  std::map<UShortArray, RealVector>::const_iterator
    w_cit = type1WeightSets.begin(), max_cit = w_cit;
  size_t num_wts, max_wts = w_cit->second.length(); ++w_cit;
  for (; w_cit!=type1WeightSets.end(); ++w_cit) {
    num_wts = w_cit->second.length();
    if (num_wts > max_wts)
      { max_wts = num_wts; max_cit = w_cit; }
  }
  return max_cit->first;
}


void CombinedSparseGridDriver::initialize_duplicate_tolerance()
{
  bool parameterized_basis = false, numerical_basis = false;
  for (size_t i=0; i<numVars; ++i) {
    short rule = collocRules[i];
    if (rule == GOLUB_WELSCH)
      { numerical_basis = true; break; }
    else if (rule == GEN_GAUSS_LAGUERRE || rule == GAUSS_JACOBI)
      parameterized_basis = true;
  }

  // Allow a looser duplication tolerance for numerically-generated rules than
  // for lookup tables.  sgmg,sgmga use an absolute tolerance which is fine for
  // standardized probability distributions (scaled to O(1)), but exhibits scale
  // dependence for numerically-generated quadrature rules.

  // rules from eigensolves or parameterized solves:
  if (numerical_basis || parameterized_basis) duplicateTol = 1.e-14;
  // rules mostly from lookup tables:
  else duplicateTol = 1.e-15; // ~= 4.5 * DBL_EPSILON

  // For numerically-generated rules, convert duplicateTol to a relative tol by
  // using a length scale.  sandia_rules.cpp uses "(dist <= tol)", so we could
  // modify sandia_rules (in many places) to use "(dist / length_scale <= tol)"
  // or just modify duplicateTol such that new_tol = tol * length_scale.  The
  // length_scale() function returns max(mean,stdev) in u-space (stdev provides
  // a backup for small/zero mean).
  if (numerical_basis) {
    Real length_scale = 0.;
    for (size_t i=0; i<numVars; ++i)
      length_scale += std::pow(polynomialBasis[i].length_scale(), 2);
                    //std::pow(std::max(ranVarMeansU[i], ranVarStdDevsU[i]), 2);
    if (length_scale > DBL_MIN)
      duplicateTol *= std::sqrt(length_scale);
  }
}


void CombinedSparseGridDriver::initialize_rule_pointers()
{
  size_t i, j;
  // compute1DPoints needed for grid_size() and for sgmg/sgmga
  compute1DPoints.resize(numVars);
  for (i=0; i<numVars; ++i)
    compute1DPoints[i] = basis_collocation_points;
  // compute1D{Type1,Type2}Weights only needed for sgmg/sgmga
  if (!refineControl) {
    compute1DType1Weights.resize(numVars);
    for (i=0; i<numVars; i++)
      compute1DType1Weights[i] = basis_type1_collocation_weights;
    /*
    if (computeType2Weights) {
      compute1DType2Weights.resize(numVars);
      for (i=0; i<numVars; ++i) {
	std::vector<CollocFnPtr>& comp_1d_t2_wts_i = compute1DType2Weights[i];
	comp_1d_t2_wts_i.resize(numVars);
	for (j=0; j<numVars; ++j)
	  comp_1d_t2_wts_i[j] = (j==i) ? basis_type2_collocation_weights :
	                                 basis_type1_collocation_weights;
      }
    }
    */
  }
}


void CombinedSparseGridDriver::initialize_growth_pointers()
{
  levelGrowthToOrder.resize(numVars);

  // if INTEGRATION with restricted growth, sync on integrand precision.
  // if INTERPOLATION with restricted growth, sync on number of points
  // (or Lagrange interpolant order = #pts - 1).
  if (driverMode == INTERPOLATION_MODE)
    for (size_t i=0; i<numVars; ++i)
      switch (collocRules[i]) {
      // nested rules with exponential growth:
      case GENZ_KEISTER:
	levelGrowthToOrder[i] = level_to_order_exp_hgk_interp;    break;
      case CLENSHAW_CURTIS: case NEWTON_COTES:
	levelGrowthToOrder[i] = level_to_order_exp_closed_interp; break;
      case GAUSS_PATTERSON: case FEJER2:
	levelGrowthToOrder[i] = level_to_order_exp_open_interp;   break;
      // non-nested or weakly-nested Gauss rules with linear growth
      case GAUSS_HERMITE: case GAUSS_LEGENDRE: // symmetric Gauss linear growth
	levelGrowthToOrder[i] = webbur::level_to_order_linear_wn; break;
      default: // asymmetric Gauss linear growth
	levelGrowthToOrder[i] = webbur::level_to_order_linear_nn; break;
      }
  else // INTEGRATION_MODE or DEFAULT_MODE
    for (size_t i=0; i<numVars; ++i)
      switch (collocRules[i]) {
      // nested rules with exponential growth:
      case GAUSS_PATTERSON:
	levelGrowthToOrder[i] = webbur::level_to_order_exp_gp;    break;
      case GENZ_KEISTER:
	levelGrowthToOrder[i] = webbur::level_to_order_exp_hgk;   break;
      case CLENSHAW_CURTIS: case NEWTON_COTES:
	levelGrowthToOrder[i] = webbur::level_to_order_exp_cc;    break;
      case FEJER2:
	levelGrowthToOrder[i] = webbur::level_to_order_exp_f2;    break;
      // non-nested or weakly-nested Gauss rules with linear growth
      case GAUSS_HERMITE: case GAUSS_LEGENDRE: // symmetric Gauss linear growth
	levelGrowthToOrder[i] = webbur::level_to_order_linear_wn; break;
      default: // asymmetric Gauss linear growth
	levelGrowthToOrder[i] = webbur::level_to_order_linear_nn; break;
      }
}


void CombinedSparseGridDriver::
assign_smolyak_arrays(UShort2DArray& sm_mi, IntArray& sm_coeffs)
{
  // Populate smolyakMultiIndex and smolyakCoeffs.  Identifies
  // use of polynomialBasis[variable][index] based on index 0:num_levels-1.
  // w = q - N = dimension-independent level.  For isotropic,
  //   w + 1 <= |i| <= w + N for i starts at 1 (used for index set defn.)
  //   w - N + 1 <= |j| <= w for j = i - 1 starts at 0 (used for generation)
  // For anisotropic, a weighted linear index set constraint is used.

  size_t i;
  unsigned short ssg_lev = ssgLevIter->second;
  if (dimIsotropic) { // initialize multi_index
    UShortArray levels(numVars, ssg_lev);
    SharedPolyApproxData::total_order_multi_index(levels, sm_mi, numVars-1);
    size_t num_terms = sm_mi.size();
    // initialize sm_coeffs
    sm_coeffs.resize(num_terms);
    for (i=0; i<num_terms; i++) {
      int wpNmi = ssg_lev - l1_norm(sm_mi[i]); // w+N-|i| = w-|j|
      sm_coeffs[i] = (int)std::pow(-1., wpNmi)
	* (int)std::floor(BasisPolynomial::n_choose_k(numVars - 1, wpNmi)+.5);
    }
  }
  else { // utilize webbur::sgmga_vcn_{ordered,coef}
    sm_mi.clear();
    sm_coeffs.clear();
    // Utilize webbur::sandia_sgmga_vcn_{ordered,coef} for 0-based index sets
    // (w*alpha_min-|alpha| < |alpha . j| <= w*alpha_min).
    // With scaling alpha_min = 1: w-|alpha| < |alpha . j| <= w.
    // In the isotropic case, reduces to w-N < |j| <= w, which is the same as
    // w-N+1 <= |j| <= w.
    IntArray x(numVars), x_max(numVars); //x_max = ssg_lev;
    UShortArray index_set(numVars);
    Real wt_sum = 0., q_max = ssg_lev;
    for (i=0; i<numVars; ++i) {
      const Real& wt_i = anisoLevelWts[i];
      wt_sum += wt_i;
      // minimum nonzero weight is scaled to 1, so just catch special case of 0
      x_max[i] = (wt_i > 1.e-10) ? (int)std::ceil(q_max/wt_i) : 0;
    }
    Real q_min = ssg_lev - wt_sum;
#ifdef DEBUG
    PCout << "q_min = " << q_min << " q_max = " << q_max;
#endif // DEBUG

    bool more = false;
    Real *aniso_wts = anisoLevelWts.values();
    int  *x0 = &x[0], *xm0 = &x_max[0], coeff;
    webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0,
				     q_min, q_max, &more);
    while (more) {
      coeff = (int)webbur::sandia_sgmga_vcn_coef(numVars, aniso_wts, x0, q_max);
      if (coeff) {
	sm_coeffs.push_back(coeff);
	for (i=0; i<numVars; ++i)
	  index_set[i] = (unsigned short)x[i];
	sm_mi.push_back(index_set);
      }
      webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0,
				       q_min, q_max, &more);
    }
  }

#ifdef DEBUG
  size_t num_terms = sm_coeffs.size();
  PCout << "\nnum Smolyak terms = " << num_terms << '\n';
  for (i=0; i<num_terms; i++)
    PCout << "multi_index[" << i << "]:\n" << sm_mi[i]
	  << "coeffs[" << i << "] = " << sm_coeffs[i] << "\n\n";
#endif // DEBUG
}


void CombinedSparseGridDriver::assign_collocation_key()
{
  // define mapping from 1:numCollocPts to set of 1d interpolation indices
  const UShort2DArray& sm_mi = smolMIIter->second;
  UShort3DArray& colloc_key = collocKeyIter->second;
  size_t i, num_smolyak_indices = sm_mi.size();
  colloc_key.resize(num_smolyak_indices);
  UShortArray quad_order(numVars); //, collocation_indices(numVars);
  for (i=0; i<num_smolyak_indices; ++i) {
    level_to_order(sm_mi[i], quad_order);
    SharedPolyApproxData::
      tensor_product_multi_index(quad_order, colloc_key[i], false);
  }
}


void CombinedSparseGridDriver::
assign_collocation_indices(const IntArray& unique_index_map, size_t start_index)
{
  // define mapping from 1:numCollocPts to set of 1d interpolation indices
  const UShort3DArray& colloc_key = collocKeyIter->second;
  Sizet2DArray&        colloc_ind = collocIndIter->second;
  size_t i, j, num_tp_pts, cntr = 0, num_sm_indices = colloc_key.size();
  colloc_ind.resize(num_sm_indices);
  // unique_index_map covers both reference and increment from start_index
  for (i=0; i<start_index; ++i)
    cntr += colloc_key[i].size();
  for (i=start_index; i<num_sm_indices; ++i) {
    num_tp_pts = colloc_key[i].size();
    SizetArray& indices_i = colloc_ind[i];
    indices_i.resize(num_tp_pts);
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      indices_i[j] = unique_index_map[cntr];
#ifdef DEBUG
      PCout << "collocKey[" << i << "][" << j << "]:\n" << colloc_key[i][j]
	    << "collocIndices[" << i << "][" << j << "] = " << indices_i[j]
	    << '\n';
#endif // DEBUG
    }
  }
}


int CombinedSparseGridDriver::grid_size()
{
  if (updateGridSize) {
    sgdInstance = this; // sgdInstance required within compute1DPoints below
    unsigned short ssg_lev = ssgLevIter->second;
    numCollocPts = (dimIsotropic) ?
      webbur::sgmg_size(numVars, ssg_lev, &compute1DPoints[0], duplicateTol,
	growthRate, &levelGrowthToOrder[0]) :
      webbur::sandia_sgmga_size(numVars, anisoLevelWts.values(), ssg_lev,
	&compute1DPoints[0], duplicateTol, growthRate, &levelGrowthToOrder[0]);
    updateGridSize = false;
  }
  return numCollocPts;
}


void CombinedSparseGridDriver::
reinterpolated_tensor_grid(const UShortArray& lev_index,
			   const SizetList& reinterp_indices)
{
  std::map<UShortArray, size_t>::iterator map_it = reinterpMap.find(lev_index);
  if (map_it == reinterpMap.end()) {

    // update arrays in place
    activeReinterpIndex = reinterpLevelIndices.size();
    UShortArray usa; UShort2DArray us2a; RealMatrix rm;
    reinterpLevelIndices.push_back(usa);
    reinterpQuadOrders.push_back(usa);
    reinterpVarSets.push_back(rm);
    reinterpCollocKeys.push_back(us2a);
    UShortArray& reinterp_lev_index  = reinterpLevelIndices.back();
    UShortArray& reinterp_quad_order = reinterpQuadOrders.back();
    reinterp_quad_order.resize(numVars); reinterp_lev_index.resize(numVars);

    // adjust level for reinterpolation of covariance (goal = 2x the interpolant
    // order).  For nested rules, this may not double the integrand exactness,
    // but it is unclear how this affects the accuracy of the interpolant.
    unsigned short l, m, target;
    SizetList::const_iterator cit = reinterp_indices.begin();
    for (size_t i=0; i<numVars; ++i) {
      if (cit != reinterp_indices.end() && i == *cit) { // reinterpolated index
	l = lev_index[i]; level_to_order(i, l, m);
	target = 2*m - 1; // target m doubles the interp order (t = 2(m-1)+1)
	while (m < target)
	  { ++l; level_to_order(i, l, m); }
	reinterp_lev_index[i] = l; reinterp_quad_order[i] = m;
	// advance to the next reinterp index
	++cit;
      }
      else { // not a reinterpolated index --> no change from reference
	reinterp_lev_index[i] = lev_index[i];
	level_to_order(i, reinterp_lev_index[i], reinterp_quad_order[i]);
      }
    }

    // compute the reinterpolation grid
    compute_tensor_grid(reinterp_quad_order, reinterp_lev_index,
			reinterp_indices, reinterpVarSets.back(),
			reinterpCollocKeys.back());

    // update reiterpMap bookkeeping
    reinterpMap[lev_index] = activeReinterpIndex;
  }
  else
    activeReinterpIndex = map_it->second;
}


void CombinedSparseGridDriver::compute_grid(RealMatrix& var_sets)
{
  assign_smolyak_arrays();

  // For efficiency reasons, incremental sparse grid definition uses
  // different point orderings than sgmg/sgmga.  Therefore, the
  // reference grid computations are kept completely separate.

  // ------------------------------------
  // Compute number of collocation points
  // ------------------------------------
  grid_size(); // ensure numCollocPts is up to date

  // ----------------------------------------------
  // Get collocation points and integration weights
  // ----------------------------------------------
  var_sets.shapeUninitialized(numVars, numCollocPts);
  RealVector& t1_wts = type1WeightSets[activeKey];
  RealMatrix& t2_wts = type2WeightSets[activeKey];
  if (trackUniqueProdWeights) {
    t1_wts.sizeUninitialized(numCollocPts);
    if (computeType2Weights)
      t2_wts.shapeUninitialized(numVars, numCollocPts);
  }
  IntArray unique_index_map;
  int* sparse_order = new int [numCollocPts*numVars];
  int* sparse_index = new int [numCollocPts*numVars];
  sgdInstance = this; // sgdInstance required within compute1D fn pointers
  unsigned short ssg_lev = ssgLevIter->second;
  if (dimIsotropic) {
    int num_total_pts = webbur::sgmg_size_total(numVars, ssg_lev,
      growthRate, &levelGrowthToOrder[0]);
    unique_index_map.resize(num_total_pts);
    webbur::sgmg_unique_index(numVars, ssg_lev, &compute1DPoints[0],
      duplicateTol, numCollocPts, num_total_pts, growthRate,
      &levelGrowthToOrder[0], &unique_index_map[0]);
    webbur::sgmg_index(numVars, ssg_lev, numCollocPts, num_total_pts,
      &unique_index_map[0], growthRate, &levelGrowthToOrder[0],
      sparse_order, sparse_index);
    webbur::sgmg_point(numVars, ssg_lev, &compute1DPoints[0], numCollocPts,
      sparse_order, sparse_index, growthRate, &levelGrowthToOrder[0],
      var_sets.values());
    if (trackUniqueProdWeights) {
      webbur::sgmg_weight(numVars, ssg_lev, &compute1DType1Weights[0],
	numCollocPts, num_total_pts, &unique_index_map[0], growthRate,
	&levelGrowthToOrder[0], t1_wts.values());
      if (computeType2Weights) {
	std::vector<CollocFnPtr> comp_1d_t2_wts = compute1DType1Weights;//copy
	RealVector t2_wt_set(numCollocPts);
	for (int i=0; i<numVars; ++i) {
	  comp_1d_t2_wts[i] = basis_type2_collocation_weights;//change ith ptr
	  webbur::sgmg_weight(numVars, ssg_lev, &comp_1d_t2_wts[0],
	    numCollocPts, num_total_pts, &unique_index_map[0], growthRate,
	    &levelGrowthToOrder[0], t2_wt_set.values());
	  copy_row(t2_wt_set, t2_wts, i);
	  comp_1d_t2_wts[i] = basis_type1_collocation_weights; // restore ptr
	}
      }
    }
  }
  else {
    int num_total_pts = webbur::sandia_sgmga_size_total(numVars,
      anisoLevelWts.values(), ssg_lev, growthRate, &levelGrowthToOrder[0]);
    unique_index_map.resize(num_total_pts);
    webbur::sandia_sgmga_unique_index(numVars, anisoLevelWts.values(),
      ssg_lev, &compute1DPoints[0], duplicateTol, numCollocPts,num_total_pts,
      growthRate, &levelGrowthToOrder[0], &unique_index_map[0]);
    webbur::sandia_sgmga_index(numVars, anisoLevelWts.values(), ssg_lev,
      numCollocPts, num_total_pts, &unique_index_map[0], growthRate,
      &levelGrowthToOrder[0], sparse_order, sparse_index);
    webbur::sandia_sgmga_point(numVars, anisoLevelWts.values(), ssg_lev,
      &compute1DPoints[0], numCollocPts, sparse_order, sparse_index,
      growthRate, &levelGrowthToOrder[0], var_sets.values());
    if (trackUniqueProdWeights) {
      webbur::sandia_sgmga_weight(numVars, anisoLevelWts.values(), ssg_lev,
        &compute1DType1Weights[0], numCollocPts, num_total_pts,
	&unique_index_map[0], growthRate, &levelGrowthToOrder[0],
	t1_wts.values());
      if (computeType2Weights) {
	std::vector<CollocFnPtr> comp_1d_t2_wts = compute1DType1Weights;//copy
	RealVector t2_wt_set(numCollocPts);
	for (int i=0; i<numVars; ++i) {
	  comp_1d_t2_wts[i] = basis_type2_collocation_weights;//change ith ptr
	  webbur::sandia_sgmga_weight(numVars, anisoLevelWts.values(),
	    ssg_lev, &comp_1d_t2_wts[0], numCollocPts, num_total_pts,
	    &unique_index_map[0], growthRate, &levelGrowthToOrder[0],
	    t2_wt_set.values());
	  copy_row(t2_wt_set, t2_wts, i);
	  comp_1d_t2_wts[i] = basis_type1_collocation_weights; // restore ptr
	}
      }
    }
  }
  delete [] sparse_order;
  delete [] sparse_index;

  if (trackCollocDetails) {
    assign_collocation_key();                     // define collocKey
    assign_collocation_indices(unique_index_map); // define collocIndices
    assign_1d_collocation_points_weights(); // define 1-D point/weight sets
  }

#ifdef DEBUG
  PCout << "CombinedSparseGridDriver::compute_grid() results:\n"
	<< "unique index mapping:\n" << unique_index_map << "\nvar_sets:\n";
  write_data(PCout, var_sets, false, true, true);
  if (trackUniqueProdWeights) {
    PCout << "\ntype1WeightSets:\n";
    write_data(PCout, type1WeightSets[activeKey]);
    if (computeType2Weights) {
      PCout << "\ntype2WeightSets:\n";
      write_data(PCout, type2WeightSets[activeKey], false, true, true);
    }
  }
#endif
}

} // namespace Pecos
