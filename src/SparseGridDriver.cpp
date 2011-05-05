/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 SparseGridDriver
//- Description: Implementation code for SparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "SparseGridDriver.hpp"
#include "PolynomialApproximation.hpp"
#include "sandia_rules.H"
#include "sandia_sgmg.H"
#include "sandia_sgmga.H"
#include "sandia_sgmgg.H"
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: SparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {

/// initialize static member pointer to active driver instance
SparseGridDriver* SparseGridDriver::sgdInstance(NULL);


void SparseGridDriver::
allocate_smolyak_arrays(UShort2DArray& multi_index, IntArray& coeffs)
{
  // Populate smolyakMultiIndex and smolyakCoeffs.  Identifies
  // use of polynomialBasis[variable][index] based on index 0:num_levels-1.
  // w = q - N = dimension-independent level.  For isotropic,
  //   w + 1 <= |i| <= w + N for i starts at 1 (used for index set defn.)
  //   w - N + 1 <= |j| <= w for j = i - 1 starts at 0 (used for generation)
  // For anisotropic, a weighted linear index set constraint is used.

  size_t i;
  if (dimIsotropic) { // initialize multi_index
    UShortArray levels(numVars, ssgLevel);
    PolynomialApproximation::total_order_multi_index(levels, multi_index,
						     numVars-1);
    size_t num_terms = multi_index.size();
    // initialize coeffs
    coeffs.resize(num_terms);
    for (i=0; i<num_terms; i++) {
      int wpNmi = ssgLevel - index_norm(multi_index[i]); // w+N-|i| = w-|j|
      coeffs[i] = (int)std::pow(-1., wpNmi)
	* (int)BasisPolynomial::n_choose_k(numVars - 1, wpNmi);
    }
  }
  else { // utilize Pecos wrapper to sgmga_vcn_{ordered,coef}
    multi_index.clear();
    coeffs.clear();
    // Utilize webbur::sandia_sgmga_vcn_{ordered,coef} for 0-based index sets
    // (w*alpha_min-|alpha| < |alpha . j| <= w*alpha_min).
    // With scaling alpha_min = 1: w-|alpha| < |alpha . j| <= w.
    // In the isotropic case, reduces to w-N < |j| <= w, which is the same as
    // w-N+1 <= |j| <= w.
    IntArray x(numVars), x_max(numVars); //x_max = ssgLevel;
    UShortArray index_set(numVars);
    Real wt_sum = 0., q_max = ssgLevel;
    for (i=0; i<numVars; ++i) {
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
    Real *aniso_wts = anisoLevelWts.values();
    int  *x0 = &x[0], *xm0 = &x_max[0], coeff;
    webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0,
				     q_min, q_max, &more);
    while (more) {
      coeff = (int)webbur::sandia_sgmga_vcn_coef(numVars, aniso_wts, x0, q_max);
      if (coeff) {
	coeffs.push_back(coeff);
	for (i=0; i<numVars; ++i)
	  index_set[i] = (unsigned short)x[i];
	multi_index.push_back(index_set);
      }
      webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0,
				       q_min, q_max, &more);
    }
  }

#ifdef DEBUG
  size_t num_terms = coeffs.size();
  PCout << "\nnum Smolyak terms = " << num_terms << '\n';
  for (i=0; i<num_terms; i++)
    PCout << "multi_index[" << i << "]:\n" << multi_index[i]
	  << "coeffs[" << i << "] = " << coeffs[i] << "\n\n";
#endif // DEBUG
}


void SparseGridDriver::
allocate_smolyak_coefficients(size_t start_index,
			      const UShort2DArray& multi_index,
			      IntArray& coeffs)
{
  size_t i, j, cntr = 0, num_sets = multi_index.size();
  if (coeffs.size() != num_sets)
    coeffs.resize(num_sets);
  int *s1 = new int [numVars*num_sets], *c1 = new int [num_sets],
      *s2 = new int [numVars];
  // initialize s1 and c1
  for (i=0; i<start_index; ++i) {
    c1[i] = coeffs[i];
    for (j=0; j<numVars; ++j, ++cntr) // no copy_data() since ushort -> int
      s1[cntr] = multi_index[i][j]; // sgmgg packs by variable groups
  }
  // for each s2, update coeffs
  for (i=start_index; i<num_sets; ++i) {
    for (j=0; j<numVars; ++j) // no copy_data() since ushort -> int
      s2[j] = multi_index[i][j];
    webbur::sandia_sgmgg_coef_inc2(numVars, i, s1, c1, s2, &coeffs[0]);
#ifdef DEBUG
    PCout << "allocate_smolyak_coefficients(): updated Smolyak coeffs =\n"
	  << coeffs << '\n';
#endif // DEBUG
    if (i<num_sets-1) { // update s1 and c1 state for next index
      for (j=0; j<=i; ++j) // coeffs updated to len i+1
	c1[j] = coeffs[j];
      for (j=0; j<numVars; ++j, ++cntr)
	s1[cntr] = s2[j];
    }
  }
  delete [] s1; delete [] c1; delete [] s2;
}


void SparseGridDriver::allocate_collocation_key()
{
  // define mapping from 1:numCollocPts to set of 1d interpolation indices
  size_t i, num_smolyak_indices = smolyakMultiIndex.size();
  collocKey.resize(num_smolyak_indices);
  UShortArray quad_order(numVars); //, collocation_indices(numVars);
  for (i=0; i<num_smolyak_indices; ++i) {
    level_to_order(smolyakMultiIndex[i], quad_order);
    PolynomialApproximation::tensor_product_multi_index(quad_order,
							collocKey[i]);
  }
}


void SparseGridDriver::update_collocation_key(size_t start_index)
{
  UShortArray quad_order(numVars);
  size_t i, num_sm_mi = smolyakMultiIndex.size();
  collocKey.resize(num_sm_mi);
  for (i=start_index; i<num_sm_mi; ++i) {
    level_to_order(smolyakMultiIndex[i], quad_order);
    PolynomialApproximation::tensor_product_multi_index(quad_order,
							collocKey[i]);
  }
}


void SparseGridDriver::allocate_collocation_indices()
{
  // define mapping from 1:numCollocPts to set of 1d interpolation indices
  size_t i, j, num_tp_pts, cntr = 0,
    num_smolyak_indices = collocKey.size();
  collocIndices.resize(num_smolyak_indices);
  for (i=0; i<num_smolyak_indices; ++i) {
    num_tp_pts = collocKey[i].size();
    SizetArray& indices_i = collocIndices[i];
    indices_i.resize(num_tp_pts);
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      indices_i[j] = uniqueIndexMapping[cntr];
#ifdef DEBUG
      PCout << "collocKey[" << i << "][" << j << "]:\n" << collocKey[i][j]
	    << "collocIndices[" << i << "][" << j << "] = " << indices_i[j]
	    << '\n';
#endif // DEBUG
    }
  }
}


void SparseGridDriver::allocate_1d_collocation_points_weights()
{
  size_t i, num_levels = ssgLevel + 1;
  if (collocPts1D.size() != num_levels || collocWts1D.size() != num_levels) {
    collocPts1D.resize(num_levels); collocWts1D.resize(num_levels);
    for (i=0; i<num_levels; ++i)
      { collocPts1D[i].resize(numVars); collocWts1D[i].resize(numVars); }
  }
  // level_index (j indexing) range is 0:w, level (i indexing) range is 1:w+1
  unsigned short level_index, order;
  for (i=0; i<numVars; i++) {
    BasisPolynomial& poly_i = polynomialBasis[i];
    for (level_index=0; level_index<num_levels; ++level_index) {
      level_to_order(i, level_index, order);
      collocPts1D[level_index][i] = poly_i.collocation_points(order);
      collocWts1D[level_index][i] = poly_i.collocation_weights(order);
    }
  }
}


void SparseGridDriver::
update_1d_collocation_points_weights(const UShortArray& trial_set,
				     const Real2DArray& pts_1d,
				     const Real2DArray& wts_1d)
{
  size_t i, num_levels = collocPts1D.size(), max_level = 0;
  for (i=0; i<numVars; ++i)
    if (trial_set[i] > max_level)
      max_level = trial_set[i];
  if (max_level >= num_levels) {
    collocPts1D.resize(max_level+1); collocWts1D.resize(max_level+1);
    for (i=num_levels; i<=max_level; ++i)
      { collocPts1D[i].resize(numVars); collocWts1D[i].resize(numVars); }
  }
  for (i=0; i<numVars; ++i) {
    unsigned short trial_index = trial_set[i];
    if (collocPts1D[trial_index][i].empty() ||
	collocWts1D[trial_index][i].empty()) {
      collocPts1D[trial_index][i] = pts_1d[i];
      collocWts1D[trial_index][i] = wts_1d[i];
    }
  }
}


void SparseGridDriver::dimension_preference(const RealVector& dim_pref)
{
  RealVector aniso_wts;
  if (!dim_pref.empty()) {
    size_t num_pref = dim_pref.length();
    aniso_wts.sizeUninitialized(num_pref);
#ifdef DEBUG
    PCout << "dimension preference:\n"; write_data(PCout, dim_pref);
#endif
    webbur::sandia_sgmga_importance_to_aniso(num_pref, dim_pref.values(),
					     aniso_wts.values());
#ifdef DEBUG
    PCout << "anisotropic weights after sandia_sgmga_importance_to_aniso():\n";
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
      webbur::sandia_sgmga_aniso_normalize(option, numVars,
					   anisoLevelWts.values());
#ifdef DEBUG
      PCout << "anisoLevelWts after sandia_sgmga_aniso_normalize():\n";
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
initialize_grid(const ShortArray& u_types,  unsigned short ssg_level,
		const RealVector& dim_pref, //short refine_type,
		short refine_control, bool store_colloc,
		bool track_ensemble_wts, bool nested_rules,
		bool equidistant_rules, short growth_rate,
		short nested_uniform_rule)
{
  numVars = u_types.size();

  //refineType         = refine_type;
  refineControl        = refine_control;
  storeCollocDetails   = store_colloc;
  trackEnsembleWeights = track_ensemble_wts;

  level(ssg_level);
  dimension_preference(dim_pref);

  // For unrestricted exponential growth, use of nested rules is restricted
  // to uniform/normal in order to enforce similar growth rates:
  if (nested_rules && growth_rate == UNRESTRICTED_GROWTH)
    for (size_t i=0; i<numVars; ++i)
      if (u_types[i] != STD_UNIFORM && u_types[i] != PIECEWISE_STD_UNIFORM && 
	  u_types[i] != STD_NORMAL)
	{ nested_rules = false; break; }
  // For MODERATE and SLOW restricted exponential growth, nested rules
  // can be used heterogeneously and synchronized with STANDARD and SLOW
  // linear growth, respectively.

  // define collocRules
  initialize_rules(u_types, nested_rules, equidistant_rules, 
		   nested_uniform_rule);
  // convert collocRules/growth_rate to apiIntegrationRules/apiGrowthRules
  initialize_api_arrays(growth_rate);
  // set compute1DPoints and compute1DWeights
  initialize_rule_pointers();
}


void SparseGridDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis,
		unsigned short ssg_level, const RealVector& dim_pref,
		//short refine_type,
		short refine_control, bool store_colloc,
		bool track_ensemble_wts, short growth_rate)
{
  numVars              = poly_basis.size();
  polynomialBasis      = poly_basis; // shallow copy
  //refineType         = refine_type;
  refineControl        = refine_control;
  storeCollocDetails   = store_colloc;
  trackEnsembleWeights = track_ensemble_wts;

  level(ssg_level);
  dimension_preference(dim_pref);

  // define collocRules
  initialize_rules(poly_basis);
  // convert collocRules/growth_rate to apiIntegrationRules/apiGrowthRules
  initialize_api_arrays(growth_rate);
  // set compute1DPoints and compute1DWeights
  initialize_rule_pointers();
}


/** Convert Pecos rule settings to IntArrays for input to VPISparseGrid. */
void SparseGridDriver::initialize_api_arrays(short growth_rate)
{
  apiIntegrationRules.resize(numVars);
  apiGrowthRules.resize(numVars);

  for (size_t i=0; i<numVars; i++) {
    // convert collocRules to apiIntegrationRules.  NEWTON_COTES is not
    // recognized by sgmg/sgmga and must be converted to the GOLUB_WELSCH
    // user-defined rule type.
    apiIntegrationRules[i] = (collocRules[i] == NEWTON_COTES) ?
      GOLUB_WELSCH : (int)collocRules[i];

    // convert growth_rate to apiGrowthRules
    switch (collocRules[i]) {
    case GAUSS_HERMITE: case GAUSS_LEGENDRE: // symmetric Gaussian linear growth
      apiGrowthRules[i] = (growth_rate == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR_ODD : MODERATE_LINEAR; break;
    case GAUSS_PATTERSON: case GENZ_KEISTER:
    case CLENSHAW_CURTIS: case FEJER2: case NEWTON_COTES:
      // nested rules with exponential growth
      switch (growth_rate) {
      case SLOW_RESTRICTED_GROWTH:
	apiGrowthRules[i] = SLOW_EXPONENTIAL;     break;
      case MODERATE_RESTRICTED_GROWTH:
	apiGrowthRules[i] = MODERATE_EXPONENTIAL; break;
      case UNRESTRICTED_GROWTH:
	apiGrowthRules[i] = FULL_EXPONENTIAL;     break;
      }
      break;
    default: // asymmetric Gaussian linear growth
      apiGrowthRules[i] = (growth_rate == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    }
  }
}


void SparseGridDriver::initialize_rule_pointers()
{
  // compute1DPoints needed for grid_size() and for sgmg/sgmga
  compute1DPoints.resize(numVars);
  for (size_t i=0; i<numVars; i++)
    compute1DPoints[i] = basis_collocation_points;
  // compute1DWeights only needed for sgmg/sgmga
  if (refineControl != DIMENSION_ADAPTIVE_GENERALIZED_SPARSE) {
    compute1DWeights.resize(numVars);
    for (size_t i=0; i<numVars; i++)
      compute1DWeights[i] = basis_collocation_weights;
  }
}


int SparseGridDriver::grid_size()
{
  // do this here (called at beginning of compute_grid()) since sgdInstance
  // required within compute1DPoints below
  sgdInstance = this;

  return (dimIsotropic) ?
    webbur::sgmg_size(numVars, ssgLevel, &apiIntegrationRules[0],
      &compute1DPoints[0], duplicateTol, &apiGrowthRules[0]) :
    webbur::sandia_sgmga_size(numVars, anisoLevelWts.values(), ssgLevel,
      &apiIntegrationRules[0], &compute1DPoints[0], duplicateTol,
      &apiGrowthRules[0]);
}


void SparseGridDriver::compute_grid(RealMatrix& var_sets)
{
  allocate_smolyak_arrays();

  // For efficiency reasons, incremental sparse grid definition uses different
  // point orderings than sgmg/sgmga.  Therefore, the reference grid
  // computations are kept completely separate.

  if (refineControl == DIMENSION_ADAPTIVE_GENERALIZED_SPARSE) {
    // compute reference grid only
    allocate_collocation_key();               // compute collocKey
    allocate_1d_collocation_points_weights(); // define 1-D point/weight sets
    reference_unique(); // updates collocIndices,uniqueIndexMapping,numCollocPts
    update_sparse_points(0, 0, a1Points, isUnique1, uniqueIndex1, var_sets);
#ifdef DEBUG
    PCout << "compute_grid() reference var_sets:\n" << var_sets;
#endif // DEBUG
  }
  else { // compute reference and any refined grids
    // --------------------------------
    // Get number of collocation points
    // --------------------------------
    // Note:  grid_size() sets sgdInstance
    // TO DO: improve efficiency of these calls through data reuse
    numCollocPts = grid_size();

    // ----------------------------------------------
    // Get collocation points and integration weights
    // ----------------------------------------------
    weightSets.sizeUninitialized(numCollocPts);
    var_sets.shapeUninitialized(numVars, numCollocPts);//Teuchos: col major
    int* sparse_order = new int [numCollocPts*numVars];
    int* sparse_index = new int [numCollocPts*numVars];
    if (dimIsotropic) {
      int num_total_pts = webbur::sgmg_size_total(numVars, ssgLevel,
	&apiIntegrationRules[0], &apiGrowthRules[0]);
      uniqueIndexMapping.resize(num_total_pts);
      webbur::sgmg_unique_index(numVars, ssgLevel, &apiIntegrationRules[0],
        &compute1DPoints[0], duplicateTol, numCollocPts, num_total_pts,
        &apiGrowthRules[0], &uniqueIndexMapping[0]);
      webbur::sgmg_index(numVars, ssgLevel, &apiIntegrationRules[0],
	numCollocPts, num_total_pts, &uniqueIndexMapping[0], &apiGrowthRules[0],
	sparse_order, sparse_index);
      webbur::sgmg_weight(numVars, ssgLevel, &apiIntegrationRules[0],
        &compute1DWeights[0], numCollocPts, num_total_pts,
        &uniqueIndexMapping[0], &apiGrowthRules[0], weightSets.values());
      webbur::sgmg_point(numVars, ssgLevel, &apiIntegrationRules[0],
        &compute1DPoints[0], numCollocPts, sparse_order, sparse_index,
        &apiGrowthRules[0], var_sets.values());
    }
    else {
      int num_total_pts = webbur::sandia_sgmga_size_total(numVars,
	anisoLevelWts.values(), ssgLevel, &apiIntegrationRules[0],
	&apiGrowthRules[0]);
      uniqueIndexMapping.resize(num_total_pts);
      webbur::sandia_sgmga_unique_index(numVars, anisoLevelWts.values(),
	ssgLevel, &apiIntegrationRules[0], &compute1DPoints[0], duplicateTol,
	numCollocPts, num_total_pts, &apiGrowthRules[0],&uniqueIndexMapping[0]);
      webbur::sandia_sgmga_index(numVars, anisoLevelWts.values(), ssgLevel,
        &apiIntegrationRules[0], numCollocPts, num_total_pts,
        &uniqueIndexMapping[0], &apiGrowthRules[0], sparse_order, sparse_index);
      webbur::sandia_sgmga_weight(numVars, anisoLevelWts.values(), ssgLevel,
        &apiIntegrationRules[0], &compute1DWeights[0], numCollocPts,
	num_total_pts, &uniqueIndexMapping[0], &apiGrowthRules[0],
	weightSets.values());
      webbur::sandia_sgmga_point(numVars, anisoLevelWts.values(), ssgLevel,
        &apiIntegrationRules[0], &compute1DPoints[0], numCollocPts,
	sparse_order, sparse_index, &apiGrowthRules[0], var_sets.values());
    }
    delete [] sparse_order;
    delete [] sparse_index;

    if (storeCollocDetails) {
      allocate_collocation_key();               // compute collocKey
      allocate_collocation_indices();           // compute collocIndices
      allocate_1d_collocation_points_weights(); // define 1-D point/weight sets
    }
  }

#ifdef DEBUG
  PCout << "uniqueIndexMapping:\n" << uniqueIndexMapping
	<< "\nvar_sets:\n"; write_data(PCout, var_sets, false, true, true);
  PCout << "\nweightSets:\n"; write_data(PCout, weightSets); PCout << '\n';
#endif
}


void SparseGridDriver::initialize_sets()
{
  // define set O (old) from smolyakMultiIndex and smolyakCoeffs:
  //oldMultiIndex = smolyakMultiIndex;
  oldMultiIndex.clear();
  oldMultiIndex.insert(smolyakMultiIndex.begin(), smolyakMultiIndex.end());
  update_reference();

  // compute initial set A (active) by applying add_active_neighbors()
  // to the frontier of smolyakMultiIndex:
  size_t i, num_old_sets = smolyakCoeffs.size();
  // anisotropic test on coeff==1 is necessary but not sufficient for presence
  // on index set frontier, requiring an additional logic test within
  // add_active_neighbors().  This is currently the best we can do since the
  // weighted norm of the index set may differ from the level.
  for (i=0; i<num_old_sets; ++i)
    if ( smolyakCoeffs[i] == 1 && ( !dimIsotropic || // imperfect for aniso
	 ( dimIsotropic && index_norm(smolyakMultiIndex[i]) == ssgLevel ) ) )
	add_active_neighbors(smolyakMultiIndex[i]);

#ifdef DEBUG
  PCout << "SparseGridDriver::initialize_sets():\nold sets:\n" << oldMultiIndex
	<< "active sets:\n" << activeMultiIndex << std::endl;
#endif // DEBUG
}


void SparseGridDriver::push_trial_set(const UShortArray& set)
{
  size_t last_index = smolyakMultiIndex.size();
  smolyakMultiIndex.push_back(set);

  // update smolyakCoeffs from smolyakMultiIndex
  allocate_smolyak_coefficients(last_index);

  // collocKey, collocIndices, and uniqueIndexMapping updated within
  // either restore_set() or compute_trial_grid()
}


void SparseGridDriver::restore_set()
{
  // SparseGridDriver currently retains no memory, so updates are recomputed

  size_t last_index = smolyakMultiIndex.size() - 1;
  // update collocKey
  update_collocation_key(last_index);
  // update collocIndices and uniqueIndexMapping
  increment_unique();
}


void SparseGridDriver::compute_trial_grid(RealMatrix& unique_var_sets)
{
  // compute trial variable/weight sets and update collocKey
  const UShortArray& trial_set = smolyakMultiIndex.back();
  UShortArray quad_order(numVars);
  level_to_order(trial_set, quad_order);
  Real2DArray pts_1d(numVars), wts_1d(numVars);
  UShort2DArray new_key;
  collocKey.push_back(new_key); // empty array updated in place
  compute_tensor_grid(quad_order, a2Points, a2Weights, collocKey.back(),
		      pts_1d, wts_1d);

  // track trial sets that have been evaluated (do here since
  // push_trial_set() used for both new trials and restorations)
  trialSets.insert(trial_set);

  // update 3D with new 2D collocation pts/wts (in correct location)
  update_1d_collocation_points_weights(trial_set, pts_1d, wts_1d);

  // update collocIndices and uniqueIndexMapping
  size_t i, num_tp_pts, index, last_index = smolyakMultiIndex.size() - 1;
  increment_unique(false); // don't recompute a2 data

  // update unique_var_sets
  update_sparse_points(last_index, numUnique1, a2Points, isUnique2,
		       uniqueIndex2, unique_var_sets);
#ifdef DEBUG
  PCout << "compute_trial_grid() increment:\nunique variable sets:\n"
	<< unique_var_sets;
#endif // DEBUG
}


void SparseGridDriver::pop_trial_set()
{
  numCollocPts -= numUnique2; // subtract number of trial points
  uniqueIndexMapping.resize(numCollocPts); // prune trial set from end
  smolyakMultiIndex.pop_back();
  smolyakCoeffs = smolyakCoeffsRef;
  collocKey.pop_back();
  collocIndices.pop_back();
}


void SparseGridDriver::update_sets(const UShortArray& set_star)
{
  // set_star is passed as *cit_star from the best entry in activeMultiIndex.
  // Therefore, we must use caution in updates to activeMultiIndex that can
  // invalidate cit_star.

  // update evaluation set smolyakMultiIndex (permanently, will not be popped)
  push_trial_set(set_star);
  restore_set();  // calls increment_unique() --> INC2
  merge_unique(); // reset a1 --> INC3

  // use smolyakMultiIndex's copy, rather than incoming set_star due to
  // iterator invalidation
  const UShortArray& last_sm_set = smolyakMultiIndex.back();

  // update set O by adding set_star to oldMultiIndex:
  oldMultiIndex.insert(last_sm_set);
  // remove set_star from set A by erasing from activeMultiIndex:
  activeMultiIndex.erase(last_sm_set); // invalidates cit_star -> set_star
  // update subset of A that have been evaluated as trial sets
  trialSets.erase(last_sm_set);

  // update set A (activeMultiIndex) based on neighbors of set_star:
  add_active_neighbors(last_sm_set);

  // Note: pruning irrelevant sets that have Coeff = 0 would be tricky,
  //       since a 0 close to the frontier can become nonzero

#ifdef DEBUG
  PCout << "Sets updated: (Smolyak,Old,Active,Trial) = ("
	<< smolyakMultiIndex.size() << ',' << oldMultiIndex.size() << ','
	<< activeMultiIndex.size() << ',' << trialSets.size() << ')'<<std::endl;
#endif // DEBUG
}


void SparseGridDriver::print_final_sets(bool converged_within_tol)
{
  // this call should precede finalize_sets()
  size_t i, j;
  if (converged_within_tol) {
    size_t last = smolyakMultiIndex.size() - 1;
    PCout << "Above tolerance index sets:\n";
    for (i=0; i<last; ++i) {
      for (j=0; j<numVars; ++j)
	PCout << std::setw(5) << smolyakMultiIndex[i][j];
      PCout << '\n';
    }
    PCout << "Below tolerance index sets:\n";
    for (j=0; j<numVars; ++j)
      PCout << std::setw(5) << smolyakMultiIndex[last][j];
    PCout << '\n';
  }
  else {
    size_t num_sm_mi = smolyakMultiIndex.size();
    PCout << "Final index sets:\n";
    for (i=0; i<num_sm_mi; ++i) {
      for (j=0; j<numVars; ++j)
	PCout << std::setw(5) << smolyakMultiIndex[i][j];
      PCout << '\n';
    }
  }
  std::set<UShortArray>::const_iterator cit;
  for (cit=trialSets.begin(); cit!=trialSets.end(); ++cit) {
    const UShortArray& mi = *cit;
    for (j=0; j<numVars; ++j)
      PCout << std::setw(5) << mi[j];
    PCout << '\n';
  }
}


void SparseGridDriver::finalize_sets()
{
  // For final answer, push all evaluated sets into old and clear active.
  // Multiple trial insertion approach must be compatible with
  // Dakota::Approximation::savedSDPSet behavior (i.e., inc2/inc3 set 
  // insertions must occur one at a time without mixing).

  size_t start_index = smolyakMultiIndex.size();
  // don't insert activeMultiIndex, as this may include sets which have not
  // been evaluated (due to final update_sets() call); instead use trialSets
  smolyakMultiIndex.insert(smolyakMultiIndex.end(), trialSets.begin(),
			   trialSets.end());
  activeMultiIndex.clear(); trialSets.clear();
  // update smolyakCoeffs from smolyakMultiIndex
  allocate_smolyak_coefficients(start_index);
  // update collocKey
  update_collocation_key(start_index);
  // update a2 data, uniqueIndexMapping, collocIndices, numCollocPts
  finalize_unique(start_index);// assure no mixing of discrete a2's
  //merge_unique(); // a1 reference update not needed, no addtnl increments
  //update_reference();
}


void SparseGridDriver::add_active_neighbors(const UShortArray& set)
{
  UShortArray trial_set = set;
  std::set<UShortArray>::const_iterator cit;
  size_t i, j;
  for (i=0; i<numVars; ++i) {
    // i^{th} candidate for set A (active) computed from forward neighbor:
    // increment by 1 in dimension i
    unsigned short& trial_set_i = trial_set[i];
    trial_set_i += 1;
    // anisotropic initialize_sets() candidates could be in oldMultiIndex
    // since smolyakCoeffs[i]==1 test is necessary but not sufficient
    if (dimIsotropic || oldMultiIndex.find(trial_set) == oldMultiIndex.end()) {
      // test all backwards neighbors for membership in set O (old)
      bool backward_old = true;
      for (j=0; j<numVars; ++j) {
	unsigned short& trial_set_j = trial_set[j];
	if (trial_set_j) { // if 0, then admissible by default
	  trial_set_j -= 1;
	  cit = oldMultiIndex.find(trial_set);
	  trial_set_j += 1; // restore
	  if (cit == oldMultiIndex.end())
	    { backward_old = false; break; }
	}
      }
      if (backward_old) // std::set<> will discard any active duplicates
	activeMultiIndex.insert(trial_set);
    }
    trial_set_i -= 1; // restore
  }
}


void SparseGridDriver::reference_unique()
{
  // define a1 pts/wts
  size_t num_sm_mi = smolyakMultiIndex.size();
  compute_tensor_points_weights(0, num_sm_mi, a1Points, a1Weights);

  // ----
  // INC1
  // ----
  int m = numVars, n1 = a1Points.numCols(), seed = 1234567;
  zVec.sizeUninitialized(m);  r1Vec.sizeUninitialized(n1);
  sortIndex1.resize(n1);      uniqueIndex1.resize(n1);
  uniqueSet1.resize(n1); // numUnique1 if count_inc1 used
  bool* is_unique1 = new bool[n1]; // BoolDeque not guaranteed contiguous

  webbur::point_radial_tol_unique_index_inc1(m, n1, a1Points.values(),
    duplicateTol, &seed, zVec.values(), r1Vec.values(), &sortIndex1[0],
    is_unique1, &numUnique1, &uniqueSet1[0], &uniqueIndex1[0]);

  copy_data(is_unique1, n1, isUnique1);
  delete [] is_unique1;

#ifdef DEBUG
  PCout << "Reference unique: numUnique1 = " << numUnique1 << "\na1 =\n"
	<< a1Points << "\n               r1   indx1 unique1   undx1   xdnu1:\n";
  for (size_t i=0; i<n1; ++i)
    std::cout << std::setw(17) << r1Vec[i]     << std::setw(8) << sortIndex1[i]
	      << std::setw(8)  << isUnique1[i] << std::setw(8) << uniqueSet1[i]
	      << std::setw(8)  << uniqueIndex1[i] << '\n';
  PCout << std::endl;
#endif // DEBUG

  uniqueIndexMapping = uniqueIndex1; // copy
  assign_tensor_collocation_indices(0, uniqueIndex1);
  numCollocPts = numUnique1;
  if (trackEnsembleWeights) {
    weightSets = 0.;
    update_sparse_weights(0, a1Weights, uniqueIndex1, weightSets);
#ifdef DEBUG
    PCout << "reference_unique() reference weightSets:\n" << weightSets;
#endif // DEBUG
  }
}


void SparseGridDriver::increment_unique(bool compute_a2)
{
  // increment_unique processes the trailing Smolyak index set
  size_t last_index = smolyakMultiIndex.size() - 1;

  // define a1 pts/wts
  if (compute_a2) // else already computed (e.g., within compute_trial_grid())
    compute_tensor_points_weights(last_index, 1, a2Points, a2Weights);

  // ----
  // INC2
  // ----
  int m = numVars, n1 = a1Points.numCols(), n2 = a2Points.numCols();
  r2Vec.sizeUninitialized(n2); sortIndex2.resize(n2);
  uniqueSet2.resize(n2); // numUnique2 if count_inc2 used
  uniqueIndex2.resize(n2);
  bool *is_unique1 = new bool[n1], *is_unique2 = new bool[n2];
  copy_data(isUnique1, is_unique1, n1);

  webbur::point_radial_tol_unique_index_inc2(m, n1, a1Points.values(),
    n2, a2Points.values(), duplicateTol, zVec.values(), r1Vec.values(),
    &sortIndex1[0], is_unique1, numUnique1, &uniqueSet1[0],
    &uniqueIndex1[0], r2Vec.values(), &sortIndex2[0], is_unique2,
    &numUnique2, &uniqueSet2[0], &uniqueIndex2[0]);

  copy_data(is_unique2, n2, isUnique2);
  delete [] is_unique1;
  delete [] is_unique2;

#ifdef DEBUG
  PCout << "Increment unique: numUnique2 = " << numUnique2 << "\na2 =\n"
	<< a2Points << "\n               r2   indx2 unique2   undx2   xdnu2:\n";
  for (size_t i=0; i<n2; ++i)
    std::cout << std::setw(17) << r2Vec[i]     << std::setw(8) << sortIndex2[i]
	      << std::setw(8)  << isUnique2[i] << std::setw(8) << uniqueSet2[i]
	      << std::setw(8)  << uniqueIndex2[i] << '\n';
  PCout << std::endl;
#endif // DEBUG

  uniqueIndexMapping.insert(uniqueIndexMapping.end(), uniqueIndex2.begin(),
			    uniqueIndex2.end());
  assign_tensor_collocation_indices(last_index, uniqueIndex2);
  numCollocPts = numUnique1 + numUnique2;
  // update weightSets
  if (trackEnsembleWeights) {
    weightSets = weightSetsRef; // to be augmented by last_index data
    update_sparse_weights(last_index, a2Weights, uniqueIndex2, weightSets);
#ifdef DEBUG
    PCout << "\nupdated weight sets:\n" << weightSets;
#endif // DEBUG
  }
}


void SparseGridDriver::merge_unique()
{
  int m = numVars, n1 = a1Points.numCols(), n2 = a2Points.numCols(),
    n1n2 = n1+n2, n3, num_unique3;
  RealVector r3_vec(n1n2, false);
  RealMatrix a3_pts(m, n1n2, false);
  IntArray sort_index3(n1n2), unique_set3(n1n2), unique_index3(n1n2);
  bool *is_unique1 = new bool[n1], *is_unique2 = new bool[n2],
       *is_unique3 = new bool[n1n2];
  copy_data(isUnique1, is_unique1, n1);
  copy_data(isUnique2, is_unique2, n2);

  // ----
  // INC3
  // ----
  webbur::point_radial_tol_unique_index_inc3(m, n1, a1Points.values(),
    r1Vec.values(), &sortIndex1[0], is_unique1, numUnique1, &uniqueSet1[0],
    &uniqueIndex1[0], n2, a2Points.values(), r2Vec.values(), &sortIndex2[0],
    is_unique2, numUnique2, &uniqueSet2[0], &uniqueIndex2[0], &n3,
    a3_pts.values(), r3_vec.values(), &sort_index3[0], is_unique3,
    &num_unique3, &unique_set3[0], &unique_index3[0]);

#ifdef DEBUG
  PCout << "Merge unique: num_unique3 = " << num_unique3 << "\na3 =\n" << a3_pts
	<< "\n               r3   indx3 unique3   undx3   xdnu3:\n";
  for (size_t i=0; i<n1n2; ++i)
    std::cout << std::setw(17)    << r3_vec[i] << std::setw(8) << sort_index3[i]
	      << std::setw(8) << is_unique3[i] << std::setw(8) << unique_set3[i]
	      << std::setw(8) << unique_index3[i] << '\n';
  PCout << std::endl;
#endif // DEBUG

  // update reference points/weights (originally defined by _inc1)
  a1Points = a3_pts;
  if (trackEnsembleWeights) {
    a1Weights.resize(n1n2);
    for (size_t i=0; i<n2; ++i)
      a1Weights[n1+i] = a2Weights[i];
  }
  // update reference indices, counts, radii
  r1Vec        = r3_vec;
  sortIndex1   = sort_index3;
  numUnique1   = num_unique3;
  uniqueSet1   = unique_set3;
  uniqueIndex1 = unique_index3;
  copy_data(is_unique3, n1n2, isUnique1);
  delete [] is_unique1; delete [] is_unique2; delete [] is_unique3;
  // update uniqueIndexMapping, collocIndices, numCollocPts
  uniqueIndexMapping = unique_index3;
  //assign_tensor_collocation_indices(0, unique_index3);
  numCollocPts = num_unique3;
}


void SparseGridDriver::finalize_unique(size_t start_index)
{
  // This fn supports multiple indices and ensures no order mixing among sets
  // (which causes a bookkeeping mismatch with Dakota::Approximation::
  // savedSDPSets) by using inc2/inc3 in careful succession.

  // *** TO DO ***: This doesn't address issue of potential point replication
  // changes between initial trial set status and finalization.  Need an
  // improved mechanism for point restore/finalize in Dakota::Approximation.
  // Could add a virtual fn to interrogate collocation_indices() from the 
  // Approximation level.  Perhaps run some performance tests first to verify
  // that this condition is possible (or does structure of admissible indices
  // prevent replication in trial sets that is not first detected in old sets).

  size_t i, j, num_sm_mi = smolyakMultiIndex.size();
  int m = numVars, n1, n2, n1n2, n3, num_unique3, all_n2 = 0;
  RealVector all_a2_wts, r3_vec; RealMatrix a3_pts;
  IntArray all_unique_index2, sort_index3, unique_set3, unique_index3;
  bool *is_unique1, *is_unique2, *is_unique3;

  for (i=start_index; i<num_sm_mi; ++i) {

    compute_tensor_points_weights(i, 1, a2Points, a2Weights);
    n1 = a1Points.numCols(); n2 = a2Points.numCols();
    all_a2_wts.resize(all_n2+n2);
    for (j=0; j<n2; ++j)
      all_a2_wts[all_n2+j] = a2Weights[j];
    all_n2 += n2;

    // INC2
    r2Vec.sizeUninitialized(n2); sortIndex2.resize(n2);
    uniqueSet2.resize(n2);       uniqueIndex2.resize(n2);
    is_unique1 = new bool[n1];   is_unique2 = new bool[n2];
    copy_data(isUnique1, is_unique1, n1);
    webbur::point_radial_tol_unique_index_inc2(m, n1, a1Points.values(),
      n2, a2Points.values(), duplicateTol, zVec.values(), r1Vec.values(),
      &sortIndex1[0], is_unique1, numUnique1, &uniqueSet1[0],
      &uniqueIndex1[0], r2Vec.values(), &sortIndex2[0], is_unique2,
      &numUnique2, &uniqueSet2[0], &uniqueIndex2[0]);
#ifdef DEBUG
    PCout << "Finalize unique: numUnique2 = " << numUnique2 << "\na2 =\n"
	  << a2Points<<"\n               r2   indx2 unique2   undx2   xdnu2:\n";
    for (j=0; j<n2; ++j)
      std::cout << std::setw(17)    << r2Vec[j] << std::setw(8) << sortIndex2[j]
		<< std::setw(8) << isUnique2[j] << std::setw(8) << uniqueSet2[j]
		<< std::setw(8) << uniqueIndex2[j] << '\n';
    PCout << std::endl;
#endif // DEBUG

    all_unique_index2.insert(all_unique_index2.end(), uniqueIndex2.begin(),
			     uniqueIndex2.end());
    numCollocPts += numUnique2;

    if (i < num_sm_mi - 1) {
      // INC3
      n1n2 = n1+n2;                       r3_vec.sizeUninitialized(n1n2);
      a3_pts.shapeUninitialized(m, n1n2); sort_index3.resize(n1n2);
      unique_set3.resize(n1n2);           unique_index3.resize(n1n2);
      is_unique3 = new bool[n1n2];
      webbur::point_radial_tol_unique_index_inc3(m, n1, a1Points.values(),
        r1Vec.values(), &sortIndex1[0], is_unique1, numUnique1, &uniqueSet1[0],
        &uniqueIndex1[0], n2, a2Points.values(), r2Vec.values(), &sortIndex2[0],
        is_unique2, numUnique2, &uniqueSet2[0], &uniqueIndex2[0], &n3,
        a3_pts.values(), r3_vec.values(), &sort_index3[0], is_unique3,
        &num_unique3, &unique_set3[0], &unique_index3[0]);
#ifdef DEBUG
      PCout << "Finalize unique: num_unique3 = " << num_unique3 << "\na3 =\n"
	    << a3_pts<<"\n               r3   indx3 unique3   undx3   xdnu3:\n";
      for (j=0; j<n1n2; ++j)
	std::cout << std::setw(17) << r3_vec[j] << std::setw(8)
		  << sort_index3[j] << std::setw(8) << is_unique3[j]
		  << std::setw(8) << unique_set3[j] << std::setw(8)
		  << unique_index3[j] << '\n';
      PCout << std::endl;
#endif // DEBUG

      // update reference points, indices, counts, radii
      a1Points     = a3_pts;
      r1Vec        = r3_vec;
      sortIndex1   = sort_index3;
      numUnique1   = num_unique3;
      uniqueSet1   = unique_set3;
      uniqueIndex1 = unique_index3;
      copy_data(is_unique3, n1n2, isUnique1);
      delete [] is_unique3;
    }

    delete [] is_unique1; delete [] is_unique2;
  }

  uniqueIndexMapping.insert(uniqueIndexMapping.end(), all_unique_index2.begin(),
			    all_unique_index2.end());
  assign_tensor_collocation_indices(start_index, all_unique_index2);
  if (trackEnsembleWeights) {
    weightSets = weightSetsRef; // to be augmented
    update_sparse_weights(start_index, all_a2_wts, all_unique_index2,
			  weightSets);
#ifdef DEBUG
    PCout << "weightSets =\n"; write_data(PCout, weightSets);
#endif // DEBUG
  }
}


void SparseGridDriver::
compute_tensor_points_weights(size_t start_index, size_t num_indices,
			      RealMatrix& pts, RealVector& wts)
{
  size_t i, j, k, cntr, num_tp_pts, num_colloc_pts = 0,
    end = start_index + num_indices;
  // define num_colloc_pts
  for (i=start_index; i<end; ++i)
    num_colloc_pts += collocKey[i].size();
  // define pts/wts: wts are raw product weights; Smolyak combinatorial
  // coefficient applied in compute_grid()/compute_trial_grid()
  pts.shapeUninitialized(numVars, num_colloc_pts);
  wts.sizeUninitialized(num_colloc_pts);
  for (i=start_index, cntr=0; i<end; ++i) {
    const UShortArray& sm_index = smolyakMultiIndex[i];
    num_tp_pts = collocKey[i].size();
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      const UShortArray& key_ij = collocKey[i][j];
      Real* pt = pts[cntr]; // column vector of size m
      Real& wt = wts[cntr]; wt = 1.;
      for (k=0; k<numVars; ++k) {
	pt[k] = collocPts1D[sm_index[k]][k][key_ij[k]];
	wt   *= collocWts1D[sm_index[k]][k][key_ij[k]];
      }
    }
  }
}


void SparseGridDriver::
update_sparse_points(size_t start_index, int new_index_offset,
		     const RealMatrix& tensor_pts, const BoolDeque& is_unique,
		     const IntArray& unique_index, RealMatrix& new_sparse_pts)
{
  size_t i, j, cntr, num_sm_mi = smolyakMultiIndex.size(), num_tp_pts,
    num_pts = is_unique.size(), num_unique_pts = 0;
  for (i=0; i<num_pts; ++i)
    if (is_unique[i])
      ++num_unique_pts;

  // update sizes
  new_sparse_pts.shapeUninitialized(numVars, num_unique_pts);

  int index;
  // add contributions for new index sets
  for (i=start_index, cntr=0; i<num_sm_mi; ++i) {
    num_tp_pts = collocKey[i].size();
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      if (is_unique[cntr]) {
	index = unique_index[cntr] - new_index_offset;
	copy_data(tensor_pts[cntr], numVars, new_sparse_pts[index]);
      }
    }
  }
}


void SparseGridDriver::
update_sparse_weights(size_t start_index, const RealVector& tensor_wts,
		      const IntArray& unique_index,
		      RealVector& updated_sparse_wts)
{
  size_t i, j, cntr, num_sm_mi = smolyakMultiIndex.size(), num_tp_pts;

  // update sizes
  updated_sparse_wts.resize(numCollocPts); // new entries initialized to 0

  int index, delta_coeff, sm_coeff;
  // back out changes in Smolyak coeff for existing index sets
  for (i=0, cntr=0; i<start_index; ++i) {
    delta_coeff = smolyakCoeffs[i] - smolyakCoeffsRef[i];
    if (delta_coeff) {
      num_tp_pts = collocKey[i].size();
      for (j=0; j<num_tp_pts; ++j, ++cntr) {
	index = uniqueIndex1[cntr];
	updated_sparse_wts[index] += delta_coeff * a1Weights[cntr];
      }
    }
    else
      cntr += collocKey[i].size();
  }
  // add contributions for new index sets
  for (i=start_index, cntr=0; i<num_sm_mi; ++i) {
    sm_coeff = smolyakCoeffs[i];
    if (sm_coeff) {
      num_tp_pts = collocKey[i].size();
      for (j=0; j<num_tp_pts; ++j, ++cntr) {
	index = unique_index[cntr];
	updated_sparse_wts[index] += sm_coeff * tensor_wts[cntr];
      }
    }
    else
      cntr += collocKey[i].size();
  }
}


void SparseGridDriver::
assign_tensor_collocation_indices(size_t start_index, 
				  const IntArray& unique_index)
{
  size_t i, j, cntr, num_tp_pts, num_sm_mi = smolyakMultiIndex.size();
  if (collocIndices.size() < num_sm_mi)
    collocIndices.resize(num_sm_mi);
  for (i=start_index, cntr=0; i<num_sm_mi; ++i) {
    num_tp_pts = collocKey[i].size();
    SizetArray& indices_i = collocIndices[i];
    indices_i.resize(num_tp_pts);
    for (j=0; j<num_tp_pts; ++j, ++cntr)
      indices_i[j] = unique_index[cntr];
  }
}

void SparseGridDriver::
level_to_order(size_t i, unsigned short level, unsigned short& order)
{
  int ilevel = level, iorder;
  webbur::level_growth_to_order(1, &ilevel, &apiIntegrationRules[i],
				&apiGrowthRules[i], &iorder);
  order = iorder;
}

} // namespace Pecos
