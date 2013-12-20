/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedPolyApproxData
//- Description:  Base Class for Orthogonal/Interpolation Polynomial
//-               Approximations
//-               
//- Owner:        Mike Eldred

#ifndef SHARED_POLY_APPROX_DATA_HPP
#define SHARED_POLY_APPROX_DATA_HPP

#include "SharedBasisApproxData.hpp"
#include "SparseGridDriver.hpp"

namespace Pecos {


/// Container class for various expansion configuration options

/** The ExpansionConfigOptions class provides a simple container class
    for expansion configuration options related to data modes, verbosity,
    refinement, and VBD controls. */

class ExpansionConfigOptions
{
  //
  //- Heading: Friends
  //

  //friend class PolynomialApproximation;

public:

  /// default constructor
  ExpansionConfigOptions();
  /// constructor
  ExpansionConfigOptions(short exp_soln_approach, short output_level,
			 bool vbd_flag, unsigned short vbd_order,
			 //short refine_type,
			 short refine_cntl, int max_iter, Real conv_tol);
  /// copy constructor
  ExpansionConfigOptions(const ExpansionConfigOptions& ec_options);
  /// destructor
  ~ExpansionConfigOptions();

//private:

  /// identifies the approach taken in compute_coefficients(): QUADRATURE,
  /// CUBATURE, COMBINED_SPARSE_GRID, HIERARCHICAL_SPARSE_GRID, REGRESSION,
  /// or SAMPLING
  short expCoeffsSolnApproach;

  /// output verbosity level: {SILENT,QUIET,NORMAL,VERBOSE,DEBUG}_OUTPUT
  short outputLevel;

  /// flag indicated use of variance-based decomposition for computing
  /// Sobol' indices
  bool vbdFlag;
  /// limit for order of interactions computed in variance-based decomposition
  unsigned short vbdOrderLimit;

  // type of refinement: {NO,P,H}_REFINEMENT
  //short refinementType;
  /// approach for control of refinement: {NO,UNIFORM,LOCAL_ADAPTIVE}_CONTROL
  /// or DIMENSION_ADAPTIVE_CONTROL_{SOBOL,DECAY,GENERALIZED}
  short refinementControl;

  /// control for limiting the maximum number of iterations in adapted
  /// or iterated approximation algorithms
  int maxIterations;
  /// convergence tolerance for adapted or iterated approximation algorithms
  Real convergenceTol;
};


inline ExpansionConfigOptions::ExpansionConfigOptions():
  expCoeffsSolnApproach(SAMPLING), outputLevel(NORMAL_OUTPUT),
  vbdFlag(false), vbdOrderLimit(0), //refinementType(NO_REFINEMENT),
  refinementControl(NO_CONTROL), maxIterations(100), convergenceTol(1.e-4)
{ }


inline ExpansionConfigOptions::
ExpansionConfigOptions(short exp_soln_approach, short output_level,
		       bool vbd_flag, unsigned short vbd_order,
		       //short refine_type,
		       short refine_cntl, int max_iter, Real conv_tol):
  expCoeffsSolnApproach(exp_soln_approach), outputLevel(output_level),
  vbdFlag(vbd_flag), vbdOrderLimit(vbd_order), //refinementType(refine_type),
  refinementControl(refine_cntl), maxIterations(max_iter),
  convergenceTol(conv_tol)
{ }


inline ExpansionConfigOptions::
ExpansionConfigOptions(const ExpansionConfigOptions& ec_options):
  expCoeffsSolnApproach(ec_options.expCoeffsSolnApproach),
  outputLevel(ec_options.outputLevel), vbdFlag(ec_options.vbdFlag),
  vbdOrderLimit(ec_options.vbdOrderLimit),
  //refinementType(ec_options.refinementType),
  refinementControl(ec_options.refinementControl),
  maxIterations(ec_options.maxIterations),
  convergenceTol(ec_options.convergenceTol)
{ }


inline ExpansionConfigOptions::~ExpansionConfigOptions()
{ }


/// Container class for various basis configuration options

/** The BasisConfigOptions class provides a simple container class
    for basis configuration options related to rule nesting, 
    piecewise basis polynomials, and derivative enhancement. */

class BasisConfigOptions
{
  //
  //- Heading: Friends
  //

  //friend class PolynomialApproximation;

public:

  /// default constructor
  BasisConfigOptions();
  /// constructor
  BasisConfigOptions(bool nested_rules, bool piecewise_basis,
		     bool equidistant_rules, bool use_derivs);
  /// copy constructor
  BasisConfigOptions(const BasisConfigOptions& bc_options);
  /// destructor
  ~BasisConfigOptions();

//private:

  /// flag for use of nested integration rules
  bool nestedRules;
  /// flag for use of piecewise basis polynomials
  bool piecewiseBasis;
  /// flag for use of equidistant points for forming piecewise basis polynomials
  bool equidistantRules;
  /// flag for utilizing derivatives during formation/calculation of expansions
  bool useDerivs;
};


inline BasisConfigOptions::BasisConfigOptions():
  nestedRules(true), piecewiseBasis(false), equidistantRules(true),
  useDerivs(false)
{ }


inline BasisConfigOptions::
BasisConfigOptions(bool nested_rules, bool piecewise_basis,
		   bool equidistant_rules, bool use_derivs):
  nestedRules(nested_rules), piecewiseBasis(piecewise_basis),
  equidistantRules(equidistant_rules), useDerivs(use_derivs)
{ }


inline BasisConfigOptions::
BasisConfigOptions(const BasisConfigOptions& bc_options):
  nestedRules(bc_options.nestedRules),
  piecewiseBasis(bc_options.piecewiseBasis),
  equidistantRules(bc_options.equidistantRules), useDerivs(bc_options.useDerivs)
{ }


inline BasisConfigOptions::~BasisConfigOptions()
{ }


/// Derived approximation class for global basis polynomials.

/** The SharedPolyApproxData class provides a global approximation
    based on basis polynomials.  This includes orthogonal polynomials
    used for polynomial chaos expansions and interpolation polynomials
    used for stochastic collocation. */

class SharedPolyApproxData: public SharedBasisApproxData
{
  //
  //- Heading: Friends
  //

  friend class PolynomialApproximation;
  friend class InterpPolyApproximation;
  friend class NodalInterpPolyApproximation;
  friend class HierarchInterpPolyApproximation;
  friend class OrthogPolyApproximation;
  friend class ProjectOrthogPolyApproximation;
  friend class RegressOrthogPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  SharedPolyApproxData();
  /// standard constructor
  SharedPolyApproxData(short basis_type, size_t num_vars);
  /// alternate constructor with full configuration options
  SharedPolyApproxData(short basis_type, size_t num_vars,
		       const ExpansionConfigOptions& ec_options,
		       const BasisConfigOptions&     bc_options);
  /// destructorboth
  ~SharedPolyApproxData();

  //
  //- Heading: Virtual member functions
  //

  /// allocate the shared data prior to building the set of approximations
  virtual void allocate_data() = 0;

  /// update the shared data prior to rebuilding the set of approximations
  virtual void increment_data();
  /// decrement the previous increment and store its shared data for
  /// later restoration
  virtual void decrement_data();
  /// restores previously saved approximation data
  virtual void pre_restore_data();
  /// restores previously saved approximation data
  virtual void post_restore_data();
  /// finalizes the shared approximation data following a set of increments
  virtual void pre_finalize_data();
  /// finalizes the shared approximation data following a set of increments
  virtual void post_finalize_data();

  /// stores current approximation data for later combination
  virtual void store_data() = 0;
  /// combines current and stored approximation data
  virtual void pre_combine_data(short combine_type);
  /// combines current and stored approximation data
  virtual void post_combine_data(short combine_type);

  /// increment the Sobol' index bookkeeping due to expansion refinement
  virtual void allocate_component_sobol() = 0;
  /// increment the Sobol' index bookkeeping due to expansion refinement
  virtual void increment_component_sobol();

  //
  //- Heading: Member functions
  //

  // allocate basis_types based on u_types
  //static bool initialize_basis_types(const ShortArray& u_types,
  //				       const BasisConfigOptions& bc_options,
  //				       ShortArray& basis_types);
  /// allocate colloc_rules based on u_types and rule options
  static void initialize_collocation_rules(const ShortArray& u_types,
    const BasisConfigOptions& bc_options, ShortArray& colloc_rules);
  /// allocate poly_basis based on basis_types and colloc_rules
  static void initialize_polynomial_basis(const ShortArray& basis_types,
    const ShortArray& colloc_rules, std::vector<BasisPolynomial>& poly_basis);
  /// pass distribution parameters from dp to poly_basis
  static void update_basis_distribution_parameters(const ShortArray& u_types,
    const AleatoryDistParams& dp, std::vector<BasisPolynomial>& poly_basis);

  /// test for whether current trial set requires a new approximation
  /// increment or can be restored from a previous trial
  bool restore_available();
  /// returns index of the data set to be restored from within saved
  /// bookkeeping (e.g., PolynomialApproximation::savedLevMultiIndex)
  size_t restoration_index();
  /// returns index of the i-th data set to be restored from within saved
  /// bookkeeping (e.g., PolynomialApproximation::savedLevMultiIndex)
  /// during finalization
  size_t finalization_index(size_t i);

  // number of data points to remove in a decrement (implemented at this
  // intermediate level since surrData not defined at base level)
  //size_t pop_count();

  /// set expConfigOptions as a group (instead of per option below)
  void configuration_options(const ExpansionConfigOptions& ec_options);
  /// set basisConfigOptions (instead of per option below)
  void configuration_options(const BasisConfigOptions& bc_options);

  /// set ExpansionConfigOptions::expCoeffsSolnApproach
  void solution_approach(short soln_approach);
  /// get ExpansionConfigOptions::expCoeffsSolnApproach
  short solution_approach() const;

  /// set ExpansionConfigOptions::vbdFlag
  void vbd_flag(bool flag);
  /// get ExpansionConfigOptions::vbdFlag
  bool vbd_flag() const;

  /// set ExpansionConfigOptions::vbdOrderLimit
  void vbd_order_limit(unsigned short vbd_order);
  /// get ExpansionConfigOptions::vbdOrderLimit
  unsigned short vbd_order_limit() const;

  // set ExpansionConfigOptions::refinementType
  //void refinement_type(short refine_type);
  // get ExpansionConfigOptions::refinementType
  //short refinement_type() const;

  /// set ExpansionConfigOptions::refinementControl
  void refinement_control(short refine_cntl);
  /// get ExpansionConfigOptions::refinementControl
  short refinement_control() const;

  /// set ExpansionConfigOptions::maxIterations
  void maximum_iterations(int max_iter);
  /// get ExpansionConfigOptions::maxIterations
  int maximum_iterations() const;

  /// set ExpansionConfigOptions::convergenceTol
  void convergence_tolerance(Real conv_tol);
  /// get ExpansionConfigOptions::convergenceTol
  Real convergence_tolerance() const;

  /// return sobolIndexMap
  const BitArrayULongMap& sobol_index_map() const;

  /// set driverRep
  void integration_driver_rep(IntegrationDriver* driver_rep);

  /// set randomVarsKey
  void random_variables_key(const BitArray& random_vars_key);

  /// return the number of expansion terms for a total order expansion
  /// with the provided (anisotropic) upper_bound array specification
  static size_t total_order_terms(const UShortArray& upper_bound,
				  short lower_bound_offset = -1);
  /// return the number of expansion terms for a tensor-product expansion
  /// with the provided (anisotropic) quadrature orders (default) or
  /// expansion orders (offset = true)
  static size_t tensor_product_terms(const UShortArray& order,
				     bool include_upper_bound = true);

  /// initialize expansion multi_index using a tensor-product expansion
  static void tensor_product_multi_index(const UShortArray& order,
					 UShort2DArray& multi_index,
					 bool include_upper_bound = true);
  /// initialize multi_index using a hierarchical tensor-product expansion
  static void hierarchical_tensor_product_multi_index(
    const UShort2DArray& delta_quad, UShort2DArray& multi_index);
  /// initialize multi_index using a total-order expansion
  /// from a scalar level
  static void total_order_multi_index(unsigned short level, size_t num_vars, 
				      UShort2DArray& multi_index);
  /// initialize expansion multi_index using a total-order expansion
  /// from an upper_bound array specification
  static void total_order_multi_index(const UShortArray& upper_bound,
    UShort2DArray& multi_index, short lower_bound_offset = -1,
    size_t max_terms = _NPOS); // SIZE_MAX is a non-portable extension

  /// utility function for incrementing a set of multidimensional indices
  static void increment_indices(UShortArray& indices, const UShortArray& limits,
				bool include_upper_bound);
  /// utility function for incrementing a set of multidimensional terms
  static void increment_terms(UShortArray& terms, size_t& last_index,
			      size_t& prev_index, size_t  term_limit,
			      bool& order_complete);//,bool unique_terms=false);

protected:

  //
  //- Heading: Member functions
  //

  /// return true if matching key values within random variable subset
  bool match_random_key(const UShortArray& key_1, 
			const UShortArray& key_2) const;
  /// return true if matching variable values within nonrandom variable subset
  bool match_nonrandom_vars(const RealVector& vars_1,
			    const RealVector& vars_2) const;

  /// allocate sobolIndices and sobolIndexMap for main effects only
  void allocate_main_sobol();
  // allocate sobolIndices and sobolIndexMap for main and interaction
  // effects including m-way interactions for m <= max_order
  //void allocate_main_interaction_sobol(unsigned short max_order);

  /// Define the sobolIndexMap (which defines the set of Sobol' indices)
  /// from the incoming multi-index.  sobolIndexMap values are initialized
  /// to interaction orders, prior to updating with multi-index increments
  /// in assign_sobol_index_map_values().
  void multi_index_to_sobol_index_map(const UShort2DArray& mi);
  /// return the sobolIndexMap values to interaction orders, prior to updating
  /// with multi-index increments in assign_sobol_index_map_values().
  void reset_sobol_index_map_values();
  /// Define the mapping from sobolIndexMap into sobolIndices
  void assign_sobol_index_map_values();

  //
  //- Heading: Data
  //

  /// pointer to integration driver instance
  IntegrationDriver* driverRep;

  /// an encapsulation of expansion configuration options
  ExpansionConfigOptions expConfigOptions;
  /// an encapsulation of basis configuration options
  BasisConfigOptions basisConfigOptions;

  /// previous quadrature order;
  /// used for tracking need for expansion form updates
  UShortArray quadOrderPrev;
  /// previous Smolyak sparse grid level;
  /// used for tracking need for expansion form updates
  unsigned short ssgLevelPrev;
  /// previous Smolyak sparse grid anisotropic weighting;
  /// used for tracking need for expansion form updates
  RealVector ssgAnisoWtsPrev;

  /// array of bits identifying the random variable subset within the
  /// active variables (used in all_variables mode)
  BitArray randomVarsKey;
  /// list of indices identifying the random variable subset within the active
  /// variables (used in all_variables mode; defined from randomVarsKey)
  SizetList randomIndices;
  /// list of indices identifying the non-random variable subset within the
  /// active variables (used in all_variables mode; defined from randomVarsKey)
  SizetList nonRandomIndices;

  /// saved trial sets that were computed but not selected
  std::deque<UShortArray> savedLevMultiIndex;

  /// mapping to manage different global sensitivity index options
  /// (e.g. univariate/main effects only vs all effects)
  BitArrayULongMap sobolIndexMap;

private:

  //
  //- Heading: Data
  //
};


inline SharedPolyApproxData::SharedPolyApproxData():
  driverRep(NULL)
{ }


inline SharedPolyApproxData::
SharedPolyApproxData(short basis_type, size_t num_vars):
  SharedBasisApproxData(BaseConstructor(), basis_type, num_vars),
  driverRep(NULL)
{ }


inline SharedPolyApproxData::
SharedPolyApproxData(short basis_type, size_t num_vars,
		     const ExpansionConfigOptions& ec_options,
		     const BasisConfigOptions& bc_options):
  SharedBasisApproxData(BaseConstructor(), basis_type, num_vars),
  expConfigOptions(ec_options), basisConfigOptions(bc_options), driverRep(NULL)
{ }


inline SharedPolyApproxData::~SharedPolyApproxData()
{ }


inline void SharedPolyApproxData::
configuration_options(const ExpansionConfigOptions& ec_options)
{
  expConfigOptions = ec_options;
  // Note: additional settings that are dependent on configuration options
  // (such as SharedInterpPolyApproxData::barycentricFlag and
  // SharedNodalInterpPolyApproxData::momentInterpType) are updated at run
  // time in allocate_data().
}


inline void SharedPolyApproxData::
configuration_options(const BasisConfigOptions& bc_options)
{
  basisConfigOptions = bc_options;
  // Note: additional settings that are dependent on configuration options
  // (such as SharedInterpPolyApproxData::barycentricFlag and
  // SharedNodalInterpPolyApproxData::momentInterpType) are updated at run
  // time in allocate_data().
}


inline void SharedPolyApproxData::solution_approach(short soln_approach)
{ expConfigOptions.expCoeffsSolnApproach = soln_approach; }


inline short SharedPolyApproxData::solution_approach() const
{ return expConfigOptions.expCoeffsSolnApproach; }


inline void SharedPolyApproxData::vbd_flag(bool flag)
{ expConfigOptions.vbdFlag = flag; }


inline bool SharedPolyApproxData::vbd_flag() const
{ return expConfigOptions.vbdFlag; }


inline void SharedPolyApproxData::vbd_order_limit(unsigned short vbd_order)
{ expConfigOptions.vbdOrderLimit = vbd_order; }


inline unsigned short SharedPolyApproxData::vbd_order_limit() const
{ return expConfigOptions.vbdOrderLimit; }


//inline void SharedPolyApproxData::refinement_type(short refine_type)
//{ expConfigOptions.refinementType = refine_type; }


//inline short SharedPolyApproxData::refinement_type() const
//{ return expConfigOptions.refinementType; }


inline void SharedPolyApproxData::refinement_control(short refine_cntl)
{ expConfigOptions.refinementControl = refine_cntl; }


inline short SharedPolyApproxData::refinement_control() const
{ return expConfigOptions.refinementControl; }


inline void SharedPolyApproxData::maximum_iterations(int max_iter)
{ expConfigOptions.maxIterations = max_iter; }


inline int SharedPolyApproxData::maximum_iterations() const
{ return expConfigOptions.maxIterations; }


inline void SharedPolyApproxData::convergence_tolerance(Real conv_tol)
{ expConfigOptions.convergenceTol = conv_tol; }


inline Real SharedPolyApproxData::convergence_tolerance() const
{ return expConfigOptions.convergenceTol; }


inline const BitArrayULongMap& SharedPolyApproxData::sobol_index_map() const
{ return sobolIndexMap; }


inline void SharedPolyApproxData::
integration_driver_rep(IntegrationDriver* driver_rep)
{ driverRep = driver_rep; }


inline void SharedPolyApproxData::
random_variables_key(const BitArray& random_vars_key)
{
  randomVarsKey = random_vars_key;
  randomIndices.clear();
  nonRandomIndices.clear();
  for (size_t i=0; i<numVars; i++)
    if (random_vars_key[i]) randomIndices.push_back(i);
    else                 nonRandomIndices.push_back(i);
}


inline void SharedPolyApproxData::
increment_indices(UShortArray& indices, const UShortArray& limits,
		  bool include_upper_bound)
{
  // perform increment
  size_t n = indices.size(), increment_index = 0;
  ++indices[increment_index];
  // if limit exceeded (including or excluding upper bound value within
  // range of indices), push to next index
  while ( increment_index < n &&
	  ( ( !include_upper_bound && 
	       indices[increment_index] >= limits[increment_index] ) ||
	    (  include_upper_bound && 
	       indices[increment_index] >  limits[increment_index] ) ) ) {
    indices[increment_index] = 0;
    ++increment_index;
    if (increment_index < n)
      ++indices[increment_index];
  }
}


inline void SharedPolyApproxData::
increment_terms(UShortArray& terms, size_t& last_index, size_t& prev_index,
		size_t term_limit, bool& order_complete)//, bool unique_terms)
{
  bool increment_complete = false;
  while (!increment_complete) {
    terms[last_index] = 1; // [1,limit] (not [0,limit-1])
    ++terms[prev_index];
    if (prev_index == 0) {
      increment_complete = true;
      if (terms[prev_index] > term_limit)
	order_complete = true;
    }
    else {
      last_index = prev_index;
      --prev_index;
      if ( //( unique_terms && terms[last_index] <  terms[prev_index]) ||
	   //(!unique_terms &&
	      terms[last_index] <= terms[prev_index]) //)
	increment_complete = true;
    }
  }
}


inline bool SharedPolyApproxData::restore_available()
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  return (std::find(savedLevMultiIndex.begin(), savedLevMultiIndex.end(),
		    ssg_driver->trial_set()) != savedLevMultiIndex.end());
}


inline size_t SharedPolyApproxData::restoration_index()
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  return find_index(savedLevMultiIndex, ssg_driver->trial_set());
}


inline size_t SharedPolyApproxData::finalization_index(size_t i)
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  const std::set<UShortArray>& trial_sets = ssg_driver->computed_trial_sets();
  // {Combined,Hierarch}SparseGridDriver::finalize_sets() updates the grid data
  // with remaining computed trial sets (in sorted order from SparseGridDriver::
  // computedTrialSets).  Below, we determine the order with which these
  // appended trial sets appear in savedLevMultiIndex.
  std::set<UShortArray>::const_iterator cit = trial_sets.begin();
  std::advance(cit, i);
  return find_index(savedLevMultiIndex, *cit);
}


inline bool SharedPolyApproxData::
match_random_key(const UShortArray& key_1, const UShortArray& key_2) const
{
  SizetList::const_iterator cit;
  for (cit=randomIndices.begin(); cit!=randomIndices.end(); ++cit)
    if (key_1[*cit] != key_2[*cit])
      return false;
  return true;
}


inline bool SharedPolyApproxData::
match_nonrandom_vars(const RealVector& vars_1, const RealVector& vars_2) const
{
  SizetList::const_iterator cit;
  for (cit=nonRandomIndices.begin(); cit!=nonRandomIndices.end(); ++cit)
    if (vars_1[*cit] != vars_2[*cit]) // double match w/o tolerance
      return false;
  return true;
}

} // namespace Pecos

#endif
