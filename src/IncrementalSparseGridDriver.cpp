/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 IncrementalSparseGridDriver
//- Description: Implementation code for IncrementalSparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "IncrementalSparseGridDriver.hpp"
#include "SharedPolyApproxData.hpp"
#include "sandia_sgmgg.hpp"
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: IncrementalSparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void IncrementalSparseGridDriver::
initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
		const ShortArray& u_types,
		const ExpansionConfigOptions& ec_options,
		BasisConfigOptions& bc_options, short growth_rate,
		bool track_uniq_prod_wts)
{
  SparseGridDriver::initialize_grid(ssg_level, dim_pref, u_types, ec_options,
				    bc_options, growth_rate);
  trackUniqueProdWeights = track_uniq_prod_wts;

  // set a rule-dependent duplicateTol
  initialize_duplicate_tolerance();
  // set compute1D{Points,Type1Weights,Type2Weights}
  //initialize_rule_pointers();
  // set levelGrowthToOrder
  //initialize_growth_pointers();
}


void IncrementalSparseGridDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis)
{
  IntegrationDriver::initialize_grid(poly_basis);

  // set a rule-dependent duplicateTol
  initialize_duplicate_tolerance();
  // set compute1D{Points,Type1Weights,Type2Weights}
  //initialize_rule_pointers();
  // set levelGrowthToOrder
  //initialize_growth_pointers();
}


void IncrementalSparseGridDriver::clear_inactive()
{
  CombinedSparseGridDriver::clear_inactive();

  std::map<UShortArray, int>::iterator nu1_it = numUnique1.begin();
  std::map<UShortArray, int>::iterator nu2_it = numUnique2.begin();
  std::map<UShortArray, RealVector>::iterator z_it = zVec.begin();
  std::map<UShortArray, RealVector>::iterator r1_it = r1Vec.begin();
  std::map<UShortArray, RealVector>::iterator r2_it = r2Vec.begin();
  std::map<UShortArray, RealMatrix>::iterator a1p_it = a1Points.begin();
  std::map<UShortArray, RealVector>::iterator a11w_it = a1Type1Weights.begin();
  std::map<UShortArray, RealMatrix>::iterator a12w_it = a1Type2Weights.begin();
  std::map<UShortArray, RealMatrix>::iterator a2p_it = a2Points.begin();
  std::map<UShortArray, RealVector>::iterator a21w_it = a2Type1Weights.begin();
  std::map<UShortArray, RealMatrix>::iterator a22w_it = a2Type2Weights.begin();
  std::map<UShortArray, IntArray>::iterator si1_it = sortIndex1.begin();
  std::map<UShortArray, IntArray>::iterator si2_it = sortIndex2.begin();
  std::map<UShortArray, IntArray>::iterator us1_it = uniqueSet1.begin();
  std::map<UShortArray, IntArray>::iterator us2_it = uniqueSet2.begin();
  std::map<UShortArray, IntArray>::iterator ui1_it = uniqueIndex1.begin();
  std::map<UShortArray, IntArray>::iterator ui2_it = uniqueIndex2.begin();
  std::map<UShortArray, BitArray>::iterator iu1_it = isUnique1.begin();
  std::map<UShortArray, BitArray>::iterator iu2_it = isUnique2.begin();

  std::map<UShortArray, IntArray>::iterator scr_it = smolyakCoeffsRef.begin();
  std::map<UShortArray, RealVector>::iterator t1r_it
    = type1WeightSetsRef.begin();
  std::map<UShortArray, RealMatrix>::iterator t2r_it
    = type2WeightSetsRef.begin();

  while (a2p_it != a2Points.end())
    if (a2p_it == a2PIter) { // preserve active
      ++nu1_it; ++nu2_it; ++z_it; ++r1_it; ++r2_it; ++a1p_it; ++a11w_it;
      ++a12w_it; ++a2p_it; ++a21w_it; ++a22w_it; ++si1_it; ++si2_it;
      ++us1_it; ++us2_it; ++ui1_it; ++ui2_it; ++iu1_it; ++iu2_it; ++scr_it;
      if (trackUniqueProdWeights)
	{ ++t1r_it; if (computeType2Weights) ++t2r_it; }
    }
    else { // clear inactive: postfix increments manage iterator invalidations
      numUnique1.erase(nu1_it++);       numUnique2.erase(nu2_it++);
      zVec.erase(z_it++); r1Vec.erase(r1_it++); r2Vec.erase(r2_it++);
      a1Points.erase(a1p_it++);         a1Type1Weights.erase(a11w_it++);
      a1Type2Weights.erase(a12w_it++);  a2Points.erase(a2p_it++);
      a2Type1Weights.erase(a21w_it++);  a2Type2Weights.erase(a22w_it++);
      sortIndex1.erase(si1_it++);       sortIndex2.erase(si2_it++);
      uniqueSet1.erase(us1_it++);       uniqueSet2.erase(us2_it++);
      uniqueIndex1.erase(ui1_it++);     uniqueIndex2.erase(ui2_it++);
      isUnique1.erase(iu1_it++);        isUnique2.erase(iu2_it++);
      smolyakCoeffsRef.erase(scr_it++);
      if (trackUniqueProdWeights) {
	type1WeightSetsRef.erase(t1r_it++);
	if (computeType2Weights) type2WeightSetsRef.erase(t2r_it++);
      }
    }
}


void IncrementalSparseGridDriver::
update_smolyak_coefficients(size_t start_index,
			    const UShort2DArray& multi_index, IntArray& coeffs)
{
  size_t j, cntr = 0, num_sets = multi_index.size(), len1 = num_sets-1;
  int i, m = numVars;
  if (coeffs.size() != num_sets)
    coeffs.resize(num_sets);
  int *s1 = new int [numVars*len1], *c1 = new int [len1],
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
    webbur::sandia_sgmgg_coef_inc2(m, i, s1, c1, s2, &coeffs[0]);
#ifdef DEBUG
    PCout << "update_smolyak_coefficients(): updated Smolyak coeffs =\n"
	  << coeffs << '\n';
#endif // DEBUG
    if (i<num_sets-1) { // if not last i, update s1 and c1 state for next pass
      for (j=0; j<=i; ++j) // coeffs updated to len i+1; max len = num_sets-1
	c1[j] = coeffs[j];
      for (j=0; j<numVars; ++j, ++cntr)
	s1[cntr] = s2[j]; // max len = (num_sets-1)*numVars
    }
  }
  delete [] s1; delete [] c1; delete [] s2;
}


void IncrementalSparseGridDriver::update_collocation_key()
{
  const UShort2DArray& sm_mi = smolMIIter->second;
  UShort3DArray& colloc_key = collocKeyIter->second;
  UShortArray quad_order(numVars);
  size_t i, start_index = colloc_key.size(), num_sm_mi = sm_mi.size();
  colloc_key.resize(num_sm_mi);
  for (i=start_index; i<num_sm_mi; ++i) {
    level_to_order(sm_mi[i], quad_order);
    SharedPolyApproxData::tensor_product_multi_index(quad_order, colloc_key[i],
						     false);
  }
}


int IncrementalSparseGridDriver::grid_size()
{
  // Note: additional efficiency could be gained by defining a1Points and
  // detecting up-to-date status in compute_grid() / reference_unique().

  if (updateGridSize) {
    UShort2DArray& sm_mi = smolMIIter->second;
    assign_smolyak_multi_index(sm_mi);

    RealMatrix a1_pts, a1_t2_wts;  RealVector a1_t1_wts;
    compute_tensor_points_weights(0, sm_mi.size(), false, a1_pts,
				  a1_t1_wts, a1_t2_wts);

    int m = numVars, n1 = a1_pts.numCols(), seed = 1234567;
    RealVector zv(m, false), r1v(n1, false);
    IntArray sind1(n1), uind1(n1), uset1(n1);
    bool* is_unique1 = new bool[n1];

    webbur::point_radial_tol_unique_index_inc1(m, n1, a1_pts.values(),
      duplicateTol, &seed, zv.values(), r1v.values(), &sind1[0], is_unique1,
      &numCollocPts, &uset1[0], &uind1[0]);
    delete [] is_unique1;

    updateGridSize = false;
  }
  return numCollocPts;
}


void IncrementalSparseGridDriver::compute_grid(RealMatrix& var_sets)
{
  assign_smolyak_arrays();

  // For efficiency reasons, incremental sparse grid definition uses
  // different point orderings than sgmg/sgmga.  Therefore, the
  // reference grid computations are kept completely separate.

  // compute reference grid only
  assign_collocation_key();               // compute collocKey
  assign_1d_collocation_points_weights(); // define 1-D point/weight sets
  reference_unique(var_sets);             // define reference grid

#ifdef DEBUG
  PCout << "IncrementalSparseGridDriver::compute_grid() results:\n"
	<< "uniqueIndexMapping:\n" << uniqueIndexMapping[activeKey];
  PCout << "\nvar_sets:\n"; write_data(PCout, var_sets, false, true, true);
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


void IncrementalSparseGridDriver::compute_trial_grid(RealMatrix& var_sets)
{
  /* Old approach now embedded within increment_unique():

  // compute trial variable/weight sets and update collocKey
  UShortArray quad_order(numVars);
  level_to_order(trialSet, quad_order);
  UShort2DArray new_key;
  UShort3DArray& colloc_key = collocKeyIter->second;
  colloc_key.push_back(new_key); // empty array updated in place
  compute_tensor_grid(quad_order, trialSet, a2PIter->second,
		      a2T1WIter->second, a2T2WIter->second, colloc_key.back());
  */

  // trialSet already appended, update collocation key
  size_t last_index = collocKeyIter->second.size();
  update_collocation_key(); // needed for compute_tensor_points_weights()
  // compute a2 pts/wts; update collocIndices, uniqueIndexMapping
  increment_unique(last_index);
  // update unique var_sets
  update_sparse_points(last_index, numUniq1Iter->second, a2PIter->second,
		       isUniq2Iter->second, uniqInd2Iter->second, var_sets);

#ifdef DEBUG
  PCout << "compute_trial_grid():\nunique variable sets:\n" << var_sets;
#endif // DEBUG

  // track trial sets that have been evaluated (do here since
  // push_trial_set() used for both new trials and restorations)
  computedTrialSets[activeKey].insert(trialSet);
}


void IncrementalSparseGridDriver::compute_grid_increment(RealMatrix& var_sets)
{
  // assumes smolyak{MultiIndex,Coeffs} have been incremented already

  size_t start_index = collocKeyIter->second.size();
  // synchronize collocKey with smolyakMultiIndex
  update_collocation_key();
  // update var_sets for multiple trial sets
  increment_unique(start_index);
  // update unique var_sets
  update_sparse_points(start_index, numUniq1Iter->second, a2PIter->second,
		       isUniq2Iter->second, uniqInd2Iter->second, var_sets);
}


void IncrementalSparseGridDriver::initialize_sets()
{
  // provide a helpful error message in the case where refineControl is not
  // set for generalized adaptation
  // > this traps the error where the reference grid was computed
  //   inconsistently (see logic in compute_grid())
  if (refineControl != DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
    PCerr << "Error: IncrementalSparseGridDriver::initialize_sets() called for "
	  << "inconsistent refinement control setting (" << refineControl
	  << ")." << std::endl;
    abort_handler(-1);
  }

  // define set O (old) from smolyakMultiIndex and smolyakCoeffs:
  unsigned short     ssg_lev = ssgLevIter->second;
  const UShort2DArray& sm_mi = smolMIIter->second;
  const IntArray&  sm_coeffs = smolCoeffsIter->second;
  UShortArraySet&     old_mi = oldMultiIndex[activeKey];
  //old_mi = sm_mi;
  old_mi.clear(); old_mi.insert(sm_mi.begin(), sm_mi.end());
  update_reference();

  computedTrialSets[activeKey].clear(); // no longer cleared in finalize_sets()
  activeMultiIndex[activeKey].clear();  // also cleared in finalize_sets()

  // compute initial set A (active) by applying add_active_neighbors()
  // to the frontier of smolyakMultiIndex:
  size_t i, num_old_sets = sm_coeffs.size();
  // anisotropic test on coeff==1 is necessary but not sufficient for presence
  // on index set frontier, requiring an additional logic test within
  // add_active_neighbors().  For anisotropic, the weighted norm of the index
  // set may differ from the level --> need to compute Pareto set.
  for (i=0; i<num_old_sets; ++i)
    if ( sm_coeffs[i] == 1 && ( !dimIsotropic || // imperfect for aniso
	 ( dimIsotropic && l1_norm(sm_mi[i]) == ssg_lev ) ) )
      add_active_neighbors(sm_mi[i], dimIsotropic);

#ifdef DEBUG
  PCout << "IncrementalSparseGridDriver::initialize_sets():\n  active key:\n"
	<< activeKey << "\n  sm_mi:\n" << sm_mi << "\n  sm_coeffs:\n"
	<< sm_coeffs << "\n  ssg level = " << ssg_lev << "\n  active sets:\n"
	<< activeMultiIndex[activeKey] << std::endl;
#endif // DEBUG
}


void IncrementalSparseGridDriver::push_trial_set(const UShortArray& set)
{
  trialSet = set;
  UShort2DArray& sm_mi = smolMIIter->second;
  size_t last_index = sm_mi.size();
  sm_mi.push_back(set);

  // update smolyakCoeffs from smolyakMultiIndex
  update_smolyak_coefficients(last_index);

  // collocKey, collocIndices, and uniqueIndexMapping updated within
  // either restore_set() or compute_trial_grid()
}


void IncrementalSparseGridDriver::restore_set()
{
  // SparseGridDriver currently retains no memory, so updates are recomputed

  // synchronize collocKey with smolyakMultiIndex
  update_collocation_key();
  // compute a2; update collocIndices, uniqueIndexMapping
  // no new var_sets and 1D updates have already been performed
  increment_unique(smolMIIter->second.size()-1, false); // no 1D pts/wts update
}


void IncrementalSparseGridDriver::pop_trial_set()
{
  smolMIIter->second.pop_back();
  collocKeyIter->second.pop_back();
  collocIndIter->second.pop_back();
  smolCoeffsIter->second = smolyakCoeffsRef[activeKey];

  numCollocPts -= numUniq2Iter->second; // subtract number of trial points
  uniqueIndexMapping[activeKey].resize(numCollocPts); // prune trial from end
}


void IncrementalSparseGridDriver::
finalize_sets(bool output_sets, bool converged_within_tol)
{
  // For final answer, push all evaluated sets into old and clear active.
  // Multiple trial insertion approach must be compatible with bookkeeping
  // elsewhere (e.g., Dakota::Approximation), i.e., inc2/inc3 set insertions
  // occur one at a time without mixing.

  UShort2DArray&  sm_mi           = smolMIIter->second;
  UShortArraySet& computed_trials = computedTrialSets[activeKey];
  size_t start_index = sm_mi.size();
  // don't insert activeMultiIndex, as this may include sets which have not
  // been evaluated (due to final update_sets() call); use computedTrialSets
  sm_mi.insert(sm_mi.end(), computed_trials.begin(), computed_trials.end());
  activeMultiIndex[activeKey].clear();
  // defer since needed for SharedPolyApproxData::finalization_index()
  //computed_trials.clear();

  // update smolyakCoeffs from smolyakMultiIndex
  update_smolyak_coefficients(start_index);
  // synchronize collocKey with smolyakMultiIndex
  update_collocation_key();
  // update a2 data, uniqueIndexMapping, collocIndices, numCollocPts
  finalize_unique(start_index);// assure no mixing of discrete a2's
  //update_reference(); // not needed, no addtnl increments

  if (output_sets) {
    size_t i, j, num_sm_mi = sm_mi.size();
    if (converged_within_tol) {
      PCout << "Above tolerance index sets:\n";
      size_t last = start_index - 1;
      for (i=0; i<last; ++i)
	print_index_set(PCout, sm_mi[i]);
      PCout << "Below tolerance index sets:\n";
      for (i=last; i<num_sm_mi; ++i)
	print_index_set(PCout, sm_mi[i]);
    }
    else {
      PCout << "Final index sets:\n";
      for (i=0; i<num_sm_mi; ++i)
	print_index_set(PCout, sm_mi[i]);
    }
  }
}


void IncrementalSparseGridDriver::reference_unique(RealMatrix& var_sets)
{
  // define a1 pts/wts
  size_t num_sm_mi = smolMIIter->second.size();
  RealMatrix& a1_pts    = a1PIter->second;
  RealVector& a1_t1_wts = a1T1WIter->second;
  RealMatrix& a1_t2_wts = a1T2WIter->second;
  compute_tensor_points_weights(0, num_sm_mi,
				false, // 1d pts/wts already computed
				a1_pts, a1_t1_wts, a1_t2_wts);

  // ----
  // INC1
  // ----
  int m = numVars, n1 = a1_pts.numCols(), seed = 1234567;
  RealVector& zv  =         zVec[activeKey];    zv.sizeUninitialized(m);
  RealVector& r1v =        r1Vec[activeKey];   r1v.sizeUninitialized(n1);
  IntArray& sind1 =   sortIndex1[activeKey]; sind1.resize(n1);
  IntArray& uind1 = uniqInd1Iter->second;    uind1.resize(n1);
  IntArray& uset1 = uniqueSet1[activeKey]; 
  uset1.resize(n1); // numUnique1 if count_inc1 used
  int&  num_u1    = numUniq1Iter->second;
  bool* is_unique1 = new bool[n1];

  webbur::point_radial_tol_unique_index_inc1(m, n1, a1_pts.values(),
    duplicateTol, &seed, zv.values(), r1v.values(), &sind1[0], is_unique1,
    &num_u1, &uset1[0], &uind1[0]);

  BitArray& isu1 = isUniq1Iter->second;
  copy_data(is_unique1, n1, isu1);
  delete [] is_unique1;

#ifdef DEBUG
  PCout << "Reference unique: numUnique1 = " << num_u1 << "\na1 =\n"
	<< a1_pts << "\n               r1   indx1 unique1   undx1   xdnu1:\n";
  for (size_t i=0; i<n1; ++i)
    std::cout << std::setw(17) << r1v[i]   << std::setw(8) << sind1[i]
	      << std::setw(8)  << isu1[i]  << std::setw(8) << uset1[i]
	      << std::setw(8)  << uind1[i] << '\n';
  PCout << std::endl;
#endif // DEBUG

  uniqueIndexMapping[activeKey] = uind1; // copy
  assign_tensor_collocation_indices(0, uind1);
  numCollocPts = num_u1;
  update_sparse_points(0, 0, a1_pts, isu1, uind1, var_sets);
  if (trackUniqueProdWeights) {
    RealVector& t1_wts = type1WeightSets[activeKey];
    RealMatrix& t2_wts = type2WeightSets[activeKey];
    t1_wts = 0.; if (computeType2Weights) t2_wts = 0.;
    update_sparse_weights(0, a1_t1_wts, a1_t2_wts, uind1, t1_wts, t2_wts);
#ifdef DEBUG
    PCout << "\nreference_unique() reference type1WeightSets:\n";
    write_data(PCout, t1_wts);
    if (computeType2Weights) {
      PCout << "\nreference_unique() reference type2WeightSets:\n";
      write_data(PCout, t2_wts, false, true, true);
    }
#endif // DEBUG
  }
}


void IncrementalSparseGridDriver::
increment_unique(size_t start_index, bool update_1d_pts_wts)
{
  RealMatrix&    a1_pts =   a1PIter->second;
  RealVector& a1_t1_wts = a1T1WIter->second;
  RealMatrix& a1_t2_wts = a1T2WIter->second;

  RealMatrix&    a2_pts =   a2PIter->second;
  RealVector& a2_t1_wts = a2T1WIter->second;
  RealMatrix& a2_t2_wts = a2T2WIter->second;

  RealVector&  zv =       zVec[activeKey];

  RealVector& r1v =      r1Vec[activeKey];
  IntArray& sind1 = sortIndex1[activeKey];
  IntArray& uind1 = uniqInd1Iter->second;
  IntArray& uset1 = uniqueSet1[activeKey]; 
  int&     num_u1 = numUniq1Iter->second;
  BitArray&  isu1 =  isUniq1Iter->second;

  RealVector& r2v =      r2Vec[activeKey];
  IntArray& sind2 = sortIndex2[activeKey];
  IntArray& uind2 = uniqInd2Iter->second;
  IntArray& uset2 = uniqueSet2[activeKey]; 
  int&     num_u2 = numUniq2Iter->second;
  BitArray&  isu2 =  isUniq2Iter->second;

  size_t i, j, num_sm_mi = smolMIIter->second.size();
  int m = numVars, n1 = a1_pts.numCols(), tp_n2, n2 = 0;
  RealVector tp_t1_wts; RealMatrix tp_pts, tp_t2_wts;
  bool *is_unique1, *is_unique2;

  // compute the points/weights for each tensor grid, updating 1D if needed
  for (i=start_index; i<num_sm_mi; ++i) {
    compute_tensor_points_weights(i, 1, update_1d_pts_wts, tp_pts,
				  tp_t1_wts, tp_t2_wts);
    tp_n2 = tp_pts.numCols();
    a2_pts.reshape(numVars, n2+tp_n2);  a2_t1_wts.resize(n2+tp_n2);
    if (computeType2Weights) a2_t2_wts.reshape(numVars, n2+tp_n2);
    for (j=0; j<tp_n2; ++j) {
      copy_data(tp_pts[j], numVars, a2_pts[n2+j]);
      a2_t1_wts[n2+j] = tp_t1_wts[j];
      if (computeType2Weights)
	copy_data(tp_t2_wts[j], numVars, a2_t2_wts[n2+j]);
    }
    n2 += tp_n2;
  }
  
  // ----
  // INC2
  // ----
  // INC2 detects a2 duplication with respect to a1 reference grid as well as
  // internal duplicates (multiple TP grids in a2 have internal duplication).
  // This allows for a single update spanning multiple Smolyak index sets.
  r2v.sizeUninitialized(n2);  sind2.resize(n2);
  uset2.resize(n2);           uind2.resize(n2);
  is_unique1 = new bool[n1];  copy_data(isu1, is_unique1, n1);
  is_unique2 = new bool[n2];  // bridges inc2 to inc3: isUnique2 not needed
  webbur::point_radial_tol_unique_index_inc2(m, n1, a1_pts.values(), n2,
    a2_pts.values(), duplicateTol, zv.values(), r1v.values(), &sind1[0],
    is_unique1,  num_u1, &uset1[0], &uind1[0],  r2v.values(), &sind2[0],
    is_unique2, &num_u2, &uset2[0], &uind2[0]);
#ifdef DEBUG
  PCout << "Increment unique: numUnique2 = " << num_u2 << "\na2 =\n"
	<< a2_pts <<"\n               r2   indx2 unique2   undx2   xdnu2:\n";
  for (j=0; j<n2; ++j)
    std::cout << std::setw(17) << r2v[j]        << std::setw(8) << sind2[j]
	      << std::setw(8)  << is_unique2[j] << std::setw(8) << uset2[j]
	      << std::setw(8)  << uind2[j]      << '\n';
  PCout << std::endl;
#endif // DEBUG
  copy_data(is_unique2, n2, isu2);
  delete [] is_unique1; delete [] is_unique2;
  numCollocPts = num_u1 + num_u2;

  IntArray& uind_map = uniqueIndexMapping[activeKey];
  uind_map.insert(uind_map.end(), uind2.begin(), uind2.end());
  assign_tensor_collocation_indices(start_index, uind2);
  if (trackUniqueProdWeights) {  // update type{1,2}WeightSets
    RealVector& t1_wts = type1WeightSets[activeKey];
    RealMatrix& t2_wts = type2WeightSets[activeKey];
    t1_wts = type1WeightSetsRef[activeKey]; // to be augmented
    if (computeType2Weights)
      t2_wts = type2WeightSetsRef[activeKey]; // to be augmented
    update_sparse_weights(start_index, a2_t1_wts, a2_t2_wts, uind2,
			  t1_wts, t2_wts);
#ifdef DEBUG
    PCout << "type1WeightSets =\n"; write_data(PCout, t1_wts);
    if (computeType2Weights)
      { PCout << "type2WeightSets =\n"; write_data(PCout, t2_wts); }
#endif // DEBUG
  }
}


void IncrementalSparseGridDriver::merge_unique()
{
  RealMatrix& a1_pts    = a1PIter->second;
  RealVector& a1_t1_wts = a1T1WIter->second;
  RealMatrix& a1_t2_wts = a1T2WIter->second;

  RealMatrix& a2_pts    = a2PIter->second;
  RealVector& a2_t1_wts = a2T1WIter->second;
  RealMatrix& a2_t2_wts = a2T2WIter->second;

  RealVector& r1v =      r1Vec[activeKey];
  IntArray& sind1 = sortIndex1[activeKey];
  IntArray& uind1 = uniqInd1Iter->second;
  IntArray& uset1 = uniqueSet1[activeKey]; 
  int&     num_u1 = numUniq1Iter->second;
  BitArray&  isu1 =  isUniq1Iter->second;

  RealVector& r2v =      r2Vec[activeKey];
  IntArray& sind2 = sortIndex2[activeKey];
  IntArray& uind2 = uniqInd2Iter->second;
  IntArray& uset2 = uniqueSet2[activeKey]; 
  int&     num_u2 = numUniq2Iter->second;
  BitArray&  isu2 =  isUniq2Iter->second;

  int m = numVars, n1 = a1_pts.numCols(), n2 = a2_pts.numCols(),
    n1n2 = n1+n2, n3, num_u3;
  RealVector r3v(n1n2, false);
  RealMatrix a3_pts(m, n1n2, false);
  IntArray sind3(n1n2), uset3(n1n2), uind3(n1n2);
  bool *is_unique1 = new bool[n1], *is_unique2 = new bool[n2],
       *is_unique3 = new bool[n1n2];
  copy_data(isu1, is_unique1, n1);
  copy_data(isu2, is_unique2, n2);

  // ----
  // INC3
  // ----
  webbur::point_radial_tol_unique_index_inc3(m, n1, a1_pts.values(),
    r1v.values(), &sind1[0], is_unique1, num_u1, &uset1[0], &uind1[0], n2,
    a2_pts.values(), r2v.values(), &sind2[0], is_unique2, num_u2, &uset2[0],
    &uind2[0], &n3, a3_pts.values(), r3v.values(), &sind3[0], is_unique3,
    &num_u3, &uset3[0], &uind3[0]);

#ifdef DEBUG
  PCout << "Merge unique: num_unique3 = " << num_u3 << "\na3 =\n" << a3_pts
	<< "\n               r3   indx3 unique3   undx3   xdnu3:\n";
  for (size_t i=0; i<n1n2; ++i)
    std::cout << std::setw(17) << r3v[i]     << std::setw(8) << sind3[i]
	      << std::setw(8)  << is_unique3[i] << std::setw(8) << uset3[i]
	      << std::setw(8)  << uind3[i] << '\n';
  PCout << std::endl;
#endif // DEBUG

  // update reference points/weights (originally defined by _inc1)
  a1_pts = a3_pts;
  if (trackUniqueProdWeights) {
    a1_t1_wts.resize(n1n2);
    if (computeType2Weights) a1_t2_wts.reshape(numVars, n1n2);
    for (size_t i=0; i<n2; ++i) {
      a1_t1_wts[n1+i] = a2_t1_wts[i];
      if (computeType2Weights)
	copy_data(a2_t2_wts[i], numVars, a1_t2_wts[n1+i]);
    }
  }
  // update reference indices, counts, radii
  r1v   = r3v;    sind1  = sind3;  num_u1 = num_u3;
  uset1 = uset3;  uind1  = uind3;  copy_data(is_unique3, n1n2, isu1);
  delete [] is_unique1; delete [] is_unique2; delete [] is_unique3;
  // update uniqueIndexMapping, collocIndices, numCollocPts
  uniqueIndexMapping[activeKey] = uind3;
  //assign_tensor_collocation_indices(0, uind3);
  numCollocPts = num_u3;
}


void IncrementalSparseGridDriver::finalize_unique(size_t start_index)
{
  increment_unique(start_index, false);
  merge_unique();

  /*
  // This fn supports multiple indices and ensures no order mixing among sets
  // by using inc2/inc3 in careful succession.

  // *** TO DO ***: This doesn't address issue of potential point replication
  // changes between initial trial set status and finalization.  Need an
  // improved mechanism for point restore/finalize in Dakota::Approximation.
  // Could add a virtual fn to interrogate collocation_indices() from the 
  // Approximation level.  Perhaps run some performance tests first to verify
  // that this condition is possible (or does structure of admissible indices
  // prevent replication in trial sets that is not first detected in old sets).

  RealMatrix& a1_pts    = a1PIter->second;
  RealVector& a1_t1_wts = a1T1WIter->second;
  RealMatrix& a1_t2_wts = a1T2WIter->second;

  RealMatrix& a2_pts    = a2PIter->second;
  RealVector& a2_t1_wts = a2T1WIter->second;
  RealMatrix& a2_t2_wts = a2T2WIter->second;

  RealVector&  zv =         zVec[activeKey];

  RealVector& r1v =      r1Vec[activeKey];
  IntArray& sind1 = sortIndex1[activeKey];
  IntArray& uind1 = uniqInd1Iter->second;
  IntArray& uset1 = uniqueSet1[activeKey]; 
  int&     num_u1 = numUniq1Iter->second;
  BitArray&  isu1 =  isUniq1Iter->second;

  RealVector& r2v =      r2Vec[activeKey];
  IntArray& sind2 = sortIndex2[activeKey];
  IntArray& uind2 = uniqInd2Iter->second;
  IntArray& uset2 = uniqueSet2[activeKey]; 
  int&     num_u2 = numUniq2Iter->second;
  //BitArray&  isu2 =  isUniq2Iter->second;

  numCollocPts = num_u1;

  size_t i, j, num_sm_mi = smolMIIter->second.size();
  int m = numVars, n1, n2, n1n2, n3, num_u3, all_n2 = 0;
  RealVector all_a2t1_wts, r3v; RealMatrix a3_pts, all_a2t2_wts;
  IntArray all_uind2, sind3, uset3, uind3;
  bool *is_unique1, *is_unique2, *is_unique3;

  // Increment one TP grid at a time since we know there are no internal
  // duplicates and INC2 can manage duplicates with latest reference grid
  // (multiple TP grids in a2 would have internal duplicates that would
  // not be managed properly, since routines are focused on a1 duplication).

  for (i=start_index; i<num_sm_mi; ++i) {

    compute_tensor_points_weights(i, 1, false, a2_pts, a2_t1_wts, a2_t2_wts);
    n1 = a1_pts.numCols(); n2 = a2_pts.numCols();
    all_a2t1_wts.resize(all_n2+n2);
    if (computeType2Weights) all_a2t2_wts.reshape(numVars, all_n2+n2);
    for (j=0; j<n2; ++j) {
      all_a2t1_wts[all_n2+j] = a2_t1_wts[j];
      if (computeType2Weights)
	copy_data(a2_t2_wts[j], numVars, all_a2t2_wts[all_n2+j]);
    }
    all_n2 += n2;

    // ----
    // INC2
    // ----
    r2v.sizeUninitialized(n2);  sind2.resize(n2);
    uset2.resize(n2);           uind2.resize(n2);
    is_unique1 = new bool[n1];  copy_data(isu1, is_unique1, n1);
    is_unique2 = new bool[n2];  // bridges inc2 to inc3: isUnique2 not needed
    webbur::point_radial_tol_unique_index_inc2(m, n1, a1_pts.values(), n2,
      a2_pts.values(), duplicateTol, zv.values(), r1v.values(), &sind1[0],
      is_unique1,  num_u1, &uset1[0], &uind1[0],  r2v.values(), &sind2[0],
      is_unique2, &num_u2, &uset2[0], &uind2[0]);
#ifdef DEBUG
    PCout << "Finalize unique: numUnique2 = " << num_u2 << "\na2 =\n"
	  << a2_pts <<"\n               r2   indx2 unique2   undx2   xdnu2:\n";
    for (j=0; j<n2; ++j)
      std::cout << std::setw(17) << r2v[j]        << std::setw(8) << sind2[j]
		<< std::setw(8)  << is_unique2[j] << std::setw(8) << uset2[j]
		<< std::setw(8)  << uind2[j]      << '\n';
    PCout << std::endl;
#endif // DEBUG

    all_uind2.insert(all_uind2.end(), uind2.begin(), uind2.end());
    numCollocPts += num_u2;

    if (i < num_sm_mi - 1) { // *** TO DO ***: WHY?
      // ----
      // INC3
      // ----
      n1n2 = n1+n2;                       r3v.sizeUninitialized(n1n2);
      a3_pts.shapeUninitialized(m, n1n2); sind3.resize(n1n2);
      uset3.resize(n1n2);                 uind3.resize(n1n2);
      is_unique3 = new bool[n1n2];
      webbur::point_radial_tol_unique_index_inc3(m, n1, a1_pts.values(),
        r1v.values(), &sind1[0], is_unique1, num_u1, &uset1[0], &uind1[0], n2,
	a2_pts.values(), r2v.values(), &sind2[0], is_unique2, num_u2, &uset2[0],
	&uind2[0], &n3, a3_pts.values(), r3v.values(), &sind3[0], is_unique3,
        &num_u3, &uset3[0], &uind3[0]);
#ifdef DEBUG
      PCout << "Finalize unique: num_unique3 = " << num_u3 << "\na3 =\n"
	    << a3_pts<<"\n               r3   indx3 unique3   undx3   xdnu3:\n";
      for (j=0; j<n1n2; ++j)
	std::cout << std::setw(17) << r3v[j]        << std::setw(8) << sind3[j]
		  << std::setw(8)  << is_unique3[j] << std::setw(8) << uset3[j]
		  << std::setw(8)  << uind3[j]      << '\n';
      PCout << std::endl;
#endif // DEBUG

      // update reference points, indices, counts, radii
      a1_pts = a3_pts;  r1v   = r3v;    sind1 = sind3;
      num_u1 = num_u3;  uset1 = uset3;  uind1 = uind3;
      copy_data(is_unique3, n1n2, isu1);
      delete [] is_unique3;
    }

    delete [] is_unique1; delete [] is_unique2;
  }

  IntArray& uind_map = uniqueIndexMapping[activeKey];
  uind_map.insert(uind_map.end(), all_uind2.begin(), all_uind2.end());
  assign_tensor_collocation_indices(start_index, all_uind2);
  if (trackUniqueProdWeights) {
    RealVector& t1_wts = type1WeightSets[activeKey];
    RealMatrix& t2_wts = type2WeightSets[activeKey];
    t1_wts = type1WeightSetsRef[activeKey]; // to be augmented
    if (computeType2Weights)
      t2_wts = type2WeightSetsRef[activeKey]; // to be augmented
    update_sparse_weights(start_index, all_a2t1_wts, all_a2t2_wts,
			  all_uind2, t1_wts, t2_wts);
#ifdef DEBUG
    PCout << "type1WeightSets =\n"; write_data(PCout, t1_wts);
#endif // DEBUG
  }
  */
}


void IncrementalSparseGridDriver::
update_sparse_points(size_t start_index, int new_index_offset,
		     const RealMatrix& tensor_pts, const BitArray& is_unique,
		     const IntArray& unique_index, RealMatrix& new_sparse_pts)
{
  size_t i, j, cntr, num_sm_mi = smolMIIter->second.size(),
    num_tp_pts, num_pts = is_unique.size(), num_unique_pts = is_unique.count();

  // update sizes
  new_sparse_pts.shapeUninitialized(numVars, num_unique_pts);

  int index;
  // add contributions for new index sets
  const UShort3DArray& colloc_key = collocKeyIter->second;
  for (i=start_index, cntr=0; i<num_sm_mi; ++i) {
    num_tp_pts = colloc_key[i].size();
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      if (is_unique[cntr]) {
	index = unique_index[cntr] - new_index_offset;
	copy_data(tensor_pts[cntr], numVars, new_sparse_pts[index]);
      }
    }
  }
}


void IncrementalSparseGridDriver::
update_sparse_weights(size_t start_index, const RealVector& tensor_t1_wts,
		      const RealMatrix& tensor_t2_wts,
		      const IntArray& unique_index, RealVector& updated_t1_wts,
		      RealMatrix& updated_t2_wts)
{
  size_t i, j, k, cntr, num_sm_mi = smolMIIter->second.size(),
    num_tp_pts;

  // update sizes
  updated_t1_wts.resize(numCollocPts); // new entries initialized to 0
  if (computeType2Weights)
    updated_t2_wts.reshape(numVars, numCollocPts); // new entries init to 0

  RealVector& a1_t1_wts =    a1T1WIter->second;
  RealMatrix& a1_t2_wts =    a1T2WIter->second;
  IntArray&       uind1 = uniqInd1Iter->second;
  int index, delta_coeff, sm_coeff;
  const UShort3DArray& colloc_key =  collocKeyIter->second;
  const IntArray&       sm_coeffs = smolCoeffsIter->second;
  const IntArray&   sm_coeffs_ref = smolyakCoeffsRef[activeKey];
  // back out changes in Smolyak coeff for existing index sets
  for (i=0, cntr=0; i<start_index; ++i) {
    delta_coeff = sm_coeffs[i] - sm_coeffs_ref[i];
    if (delta_coeff) {
      num_tp_pts = colloc_key[i].size();
      for (j=0; j<num_tp_pts; ++j, ++cntr) {
	index = uind1[cntr];
	updated_t1_wts[index] += delta_coeff * a1_t1_wts[cntr];
	if (computeType2Weights) {
	  Real*       up_t2_wts_j = updated_t2_wts[index];
	  const Real* a1_t2_wts_j = a1_t2_wts[cntr];
	  for (k=0; k<numVars; ++k)
	    up_t2_wts_j[k] += delta_coeff * a1_t2_wts_j[k];
	}
      }
    }
    else
      cntr += colloc_key[i].size();
  }
  // add contributions for new index sets
  for (i=start_index, cntr=0; i<num_sm_mi; ++i) {
    sm_coeff = sm_coeffs[i];
    if (sm_coeff) {
      num_tp_pts = colloc_key[i].size();
      for (j=0; j<num_tp_pts; ++j, ++cntr) {
	index = unique_index[cntr];
	updated_t1_wts[index] += sm_coeff * tensor_t1_wts[cntr];
	if (computeType2Weights) {
	  Real*       up_t2_wts_j = updated_t2_wts[index];
	  const Real* te_t2_wts_j = tensor_t2_wts[cntr];
	  for (k=0; k<numVars; ++k)
	    up_t2_wts_j[k] += sm_coeff * te_t2_wts_j[k];
	}
      }
    }
    else
      cntr += colloc_key[i].size();
  }
}


void IncrementalSparseGridDriver::
compute_tensor_points_weights(size_t start_index, size_t num_indices,
			      bool update_1d_pts_wts, RealMatrix& pts,
			      RealVector& t1_wts, RealMatrix& t2_wts)
{
  // Requirements: updated smolMIIter->second, collocKeyIter->second
  //               for [start_index,start_index+num_indices]
  // 1D Pts/Wts will be updated as indicated by update_1d_pts_wts

  size_t i, j, k, l, cntr, num_tp_pts, num_colloc_pts = 0,
    end = start_index + num_indices;
  const UShort3DArray& colloc_key = collocKeyIter->second;
  const UShort2DArray& sm_mi = smolMIIter->second;
  // define num_colloc_pts
  for (i=start_index; i<end; ++i)
    num_colloc_pts += colloc_key[i].size();
  // define pts/wts: wts are raw product weights; Smolyak combinatorial
  // coefficient applied in compute_grid()/compute_trial_grid()
  pts.shapeUninitialized(numVars, num_colloc_pts);
  t1_wts.sizeUninitialized(num_colloc_pts);
  if (computeType2Weights)
    t2_wts.shapeUninitialized(numVars, num_colloc_pts);
  for (i=start_index, cntr=0; i<end; ++i) {
    const UShortArray& sm_index = sm_mi[i];
    if (update_1d_pts_wts) { // update collocPts1D, {type1,type2}CollocWts1D
      UShortArray quad_order(numVars);
      level_to_order(sm_index, quad_order);
      update_1d_collocation_points_weights(quad_order, sm_index);
    }
    num_tp_pts = colloc_key[i].size();
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      const UShortArray& key_ij = colloc_key[i][j];
      Real* pt    =    pts[cntr]; // column vector
      Real& t1_wt = t1_wts[cntr]; t1_wt = 1.;
      for (k=0; k<numVars; ++k) {
	pt[k]  =      collocPts1D[sm_index[k]][k][key_ij[k]];
	t1_wt *= type1CollocWts1D[sm_index[k]][k][key_ij[k]];
      }
      if (computeType2Weights) {
	Real* t2_wt = t2_wts[cntr]; // column vector
	for (k=0; k<numVars; ++k) {
	  Real& t2_wt_k = t2_wt[k]; t2_wt_k = 1.;
	  for (l=0; l<numVars; ++l)
	    t2_wt_k *= (l==k) ? type2CollocWts1D[sm_index[l]][l][key_ij[l]] :
	                        type1CollocWts1D[sm_index[l]][l][key_ij[l]];
	}
      }
    }
  }
#ifdef DEBUG
    PCout << "Tensor product weights =\ntype1:\n"; write_data(PCout, t1_wts);
    PCout << "type2:\n"; write_data(PCout, t2_wts, false, true, true);
#endif // DEBUG
}


void IncrementalSparseGridDriver::
assign_tensor_collocation_indices(size_t start_index, 
				  const IntArray& unique_index)
{
  size_t i, j, cntr, num_tp_pts,
    num_sm_mi = smolMIIter->second.size();
  const UShort3DArray& colloc_key = collocKeyIter->second;
  Sizet2DArray&        colloc_ind = collocIndIter->second;
  if (colloc_ind.size() < num_sm_mi)
    colloc_ind.resize(num_sm_mi);
  for (i=start_index, cntr=0; i<num_sm_mi; ++i) {
    num_tp_pts = colloc_key[i].size();
    SizetArray& indices_i = colloc_ind[i];
    indices_i.resize(num_tp_pts);
    for (j=0; j<num_tp_pts; ++j, ++cntr)
      indices_i[j] = unique_index[cntr];
  }
}

} // namespace Pecos
