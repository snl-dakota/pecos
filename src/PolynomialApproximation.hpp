/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        PolynomialApproximation
//- Description:  Base Class for Orthogonal/Interpolation Polynomial
//-               Approximations
//-               
//- Owner:        Mike Eldred

#ifndef POLYNOMIAL_APPROXIMATION_HPP
#define POLYNOMIAL_APPROXIMATION_HPP

#include "BasisApproximation.hpp"
#include "SurrogateData.hpp"
#include "SparseGridDriver.hpp"


namespace Pecos {

class IntegrationDriver;
class SparseGridDriver;
class DistributionParams;


/// Container class for various polynomial approximation configuration options

/** The ConfigurationOptions class provides a simple container class for a 
    polynomial approximation configuration options related to data modes,
    verbosity, and refinement and VBD controls. */

class ConfigurationOptions
{
  //
  //- Heading: Friends
  //

  friend class PolynomialApproximation;
  friend class InterpPolyApproximation;
  friend class NodalInterpPolyApproximation;
  friend class OrthogPolyApproximation;
  friend class HierarchInterpPolyApproximation;
public:

  /// default constructor
  ConfigurationOptions();
  /// constructor
  ConfigurationOptions(short exp_soln_approach, bool exp_coeff_flag,
		       bool exp_grad_flag, bool use_derivs,
		       //short output_level, short refine_type,
		       short refine_cntl, short vbd_type);
  /// destructor
  ~ConfigurationOptions();

private:

  /// identifies the approach taken in compute_coefficients():
  /// QUADRATURE, CUBATURE, SPARSE_GRID, REGRESSION, or SAMPLING
  short expCoeffsSolnApproach;
  /// flag for calculation of expansion coefficients from response values
  bool expansionCoeffFlag;
  /// flag for calculation of gradients of expansion coefficients from
  /// response gradients
  bool expansionCoeffGradFlag;
  /// flag for utilizing derivatives during formation/calculation of expansions
  bool useDerivs;

  // output verbosity level: {SILENT,QUIET,NORMAL,VERBOSE,DEBUG}_OUTPUT
  //short outputLevel;
  // nesting override options: NO_NESTING_OVERRIDE, NESTED, NON_NESTED
  //short nestingOverride;

  // type of refinement: {NO,P,H}_REFINEMENT
  //short refinementType;
  /// approach for control of refinement: {NO,UNIFORM}_CONTROL or
  /// DIMENSION_ADAPTIVE_{TOTAL_SOBOL,SPECTRAL_DECAY,GENERALIZED_SPARSE}
  short refinementControl;

  /// control for amount of data computed in variance-based decomposition:
  /// {NO,UNIVARIATE,ALL}_VBD
  short vbdControl;
};


inline ConfigurationOptions::ConfigurationOptions():
  expCoeffsSolnApproach(SAMPLING), expansionCoeffFlag(true),
  expansionCoeffGradFlag(false), useDerivs(false),
  //outputLevel(NORMAL_OUTPUT), refinementType(NO_REFINEMENT),
  refinementControl(NO_CONTROL), vbdControl(NO_VBD)
{ }


inline ConfigurationOptions::
ConfigurationOptions(short exp_soln_approach, bool exp_coeff_flag,
		     bool exp_grad_flag, bool use_derivs,
		     //short output_level, short refine_type,
		     short refine_cntl, short vbd_type):
  expCoeffsSolnApproach(exp_soln_approach), expansionCoeffFlag(exp_coeff_flag),
  expansionCoeffGradFlag(exp_grad_flag), useDerivs(use_derivs),
  //outputLevel(output_level), refinementType(refine_type),
  refinementControl(refine_cntl), vbdControl(vbd_type)
{ }


inline ConfigurationOptions::~ConfigurationOptions()
{ }


/// Derived approximation class for global basis polynomials.

/** The PolynomialApproximation class provides a global approximation
    based on basis polynomials.  This includes orthogonal polynomials
    used for polynomial chaos expansions and interpolation polynomials
    used for stochastic collocation. */

class PolynomialApproximation: public BasisApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  PolynomialApproximation();
  /// standard constructor
  PolynomialApproximation(size_t num_vars, bool use_derivs);
  /// destructorboth
  ~PolynomialApproximation();

  //
  //- Heading: Virtual member functions
  //

  /// Computes sensitivity indices according to verbosity (from main
  /// to higher order effects)
  virtual void compute_component_effects() = 0;
  /// Computes total sensitivity indices according to verbosity and
  /// existing computations from compute_component_effects()
  virtual void compute_total_effects() = 0;

  /// size derived class data attributes
  virtual void allocate_arrays() = 0;

  /// retrieve the response value for a stored expansion using the
  /// given parameter vector
  virtual Real stored_value(const RealVector& x) = 0;
  /// retrieve the response gradient for a stored expansion using the
  /// given parameter vector and default DVV
  virtual const RealVector& stored_gradient(const RealVector& x) = 0;

  /// return the mean of the expansion, treating all variables as random
  virtual Real mean() = 0;
  /// return the mean of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  virtual Real mean(const RealVector& x) = 0;
  /// return the gradient of the expansion mean for a given parameter
  /// vector, treating all variables as random
  virtual const RealVector& mean_gradient() = 0;
  /// return the gradient of the expansion mean for a given parameter vector
  /// and given DVV, treating a subset of the variables as random
  virtual const RealVector& mean_gradient(const RealVector& x,
					  const SizetArray& dvv) = 0;

  /// return the variance of the expansion, treating all variables as random
  virtual Real variance() = 0;
  /// return the variance of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  virtual Real variance(const RealVector& x) = 0;
  /// return the gradient of the expansion variance for a given parameter
  /// vector, treating all variables as random
  virtual const RealVector& variance_gradient() = 0;
  /// return the gradient of the expansion variance for a given parameter
  /// vector and given DVV, treating a subset of the variables as random
  virtual const RealVector& variance_gradient(const RealVector& x,
					      const SizetArray& dvv) = 0;

  /// return the variance of the expansion, treating all variables as random
  virtual Real covariance(PolynomialApproximation* poly_approx_2) = 0;
  /// return the variance of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  virtual Real covariance(const RealVector& x,
			  PolynomialApproximation* poly_approx_2) = 0;

  /// compute central response moments using some combination of expansion
  /// post-processing and numerical integration
  virtual void compute_moments() = 0;
  /// compute central response moments in all-variables mode using some
  /// combination of expansion post-processing and numerical integration
  virtual void compute_moments(const RealVector& x) = 0;
  /// return central response moments defined from either expansion
  /// post-processing or numerical integration
  virtual const RealVector& moments() const = 0;
  /// compute central moments of response using numerical integration
  virtual void compute_numerical_response_moments(size_t num_moments);

  /// estimate expansion coefficient decay rates for each random
  /// variable dimension (OrthogPolyApproximation only)
  virtual const RealVector& dimension_decay_rates();
  /// increment the approximation order (OrthogPolyApproximation only)
  virtual void increment_order();

  //
  //- Heading: Member functions
  //

  /// allocate basis_types based on u_types
  static bool distribution_types(const ShortArray& u_types,
				 bool piecewise_basis, bool use_derivs,
				 ShortArray& basis_types);
  /// allocate colloc_rules based on u_types and rule options
  static void distribution_rules(const ShortArray& u_types, bool nested_rules,
				 bool  piecewise_basis, bool equidistant_rules,
				 short nested_uniform_rule,
				 ShortArray& colloc_rules);
  /// allocate poly_basis based on basis_types and colloc_rules
  static void distribution_basis(const ShortArray& basis_types,
				 const ShortArray& colloc_rules,
				 std::vector<BasisPolynomial>& poly_basis);
  /// pass distribution parameters from dp to poly_basis
  static void distribution_parameters(const ShortArray& u_types,
				      const DistributionParams& dp,
				      std::vector<BasisPolynomial>& poly_basis);

  /// return expansionMoments
  const RealVector& expansion_moments() const;
  /// return numericalMoments
  const RealVector& numerical_moments() const;

  /// size component Sobol arrays
  void allocate_component_effects();
  /// size total Sobol arrays
  void allocate_total_effects();

  /// set surrData (shared representation)
  void surrogate_data(const SurrogateData& data);
  /// get surrData
  const SurrogateData& surrogate_data() const;

  // number of data points to remove in a decrement (implemented at this
  // intermediate level since surrData not defined at base level)
  //size_t pop_count();

  /// set ConfigurationOptions::expCoeffsSolnApproach
  void solution_approach(short soln_approach);
  /// get ConfigurationOptions::expCoeffsSolnApproach
  short solution_approach() const;

  /// set ConfigurationOptions::expansionCoeffFlag
  void expansion_coefficient_flag(bool coeff_flag);
  /// get ConfigurationOptions::expansionCoeffFlag
  bool expansion_coefficient_flag() const;

  /// set ConfigurationOptions::expansionCoeffGradFlag
  void expansion_coefficient_gradient_flag(bool grad_flag);
  /// get ConfigurationOptions::expansionCoeffGradFlag
  bool expansion_coefficient_gradient_flag() const;

  // set ConfigurationOptions::refinementType
  //void refinement_type(short refine_type);
  // get ConfigurationOptions::refinementType
  //short refinement_type() const;

  /// set ConfigurationOptions::refinementControl
  void refinement_control(short refine_cntl);
  /// get ConfigurationOptions::refinementControl
  short refinement_control() const;

  /// set ConfigurationOptions::vbdControl
  void vbd_control(short vbd_cntl);
  /// get ConfigurationOptions::vbdControl
  short vbd_control() const;

  /// return sobolIndexMap 
  const IntIntMap& sobol_index_map() const;
  /// return sobolIndices
  const RealVector& sobol_indices() const;
  /// return totalSobolIndices
  const RealVector& total_sobol_indices() const;

  /// set driverRep
  void integration_driver_rep(IntegrationDriver* driver_rep);

  /// set randomVarsKey
  void random_variables_key(const BoolDeque& random_vars_key);

  /// return the number of expansion terms for a total order expansion
  /// with the provided (anisotropic) upper_bound array specification
  static size_t total_order_terms(const UShortArray& upper_bound,
				  short lower_bound_offset = -1);
  /// return the number of expansion terms for a tensor-product expansion
  /// with the provided (anisotropic) quadrature orders (default) or
  /// expansion orders (offset = true)
  static size_t tensor_product_terms(const UShortArray& order,
				     bool exp_order_offset = false);

  /// initialize expansion multi_index using a tensor-product expansion
  static void tensor_product_multi_index(const UShortArray& order,
					 UShort2DArray& multi_index,
					 bool exp_order_offset = false);
  /// initialize expansion multi_index using a total-order expansion
  /// from an upper_bound array specification
  static void total_order_multi_index(const UShortArray& upper_bound,
				      UShort2DArray& multi_index,
				      short lower_bound_offset = -1,
				      size_t max_terms = UINT_MAX);

  /// utility function for incrementing a set of multidimensional indices
  static void increment_indices(UShortArray& indices, const UShortArray& limits,
				bool include_limit_equality);
  /// utility function for incrementing a set of multidimensional terms
  static void increment_terms(UShortArray& terms, size_t& last_index,
			      size_t& prev_index, const size_t& term_limit,
			      bool& order_complete);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// query trial index for presence in savedSmolyakMultiIndex, indicating
  /// the ability to restore a previously evaluated increment
  bool restore_available();
  /// returns index of the data set to be restored from within saved
  /// bookkeeping (i.e., savedSmolyakMultiIndex)
  size_t restoration_index();
  /// returns i-th index of the data sets to be restored from within saved
  /// bookkeeping (i.e., savedSmolyakMultiIndex) during finalization
  size_t finalization_index(size_t i);

  //
  //- Heading: Member functions
  //

  /// standardize third and higher central moments and eliminate excess kurtosis
  void standardize_moments(RealVector& moments);

  //
  //- Heading: Data
  //

  /// instance containing the variables/response data arrays for
  /// constructing a surrogate model
  SurrogateData surrData;

  /// pointer to integration driver instance
  IntegrationDriver* driverRep;

  /// a basic encapsulation of configuration options
  ConfigurationOptions configOptions;

  /// previous quadrature order;
  /// used for tracking need for expansion form updates
  UShortArray quadOrderPrev;
  /// previous Smolyak sparse grid level;
  /// used for tracking need for expansion form updates
  unsigned short ssgLevelPrev;
  /// previous Smolyak sparse grid anisotropic weighting;
  /// used for tracking need for expansion form updates
  RealVector ssgAnisoWtsPrev;

  /// array of booleans identifying the random variable subset within
  /// the active variables (used in all_variables mode)
  BoolDeque randomVarsKey;
  /// list of indices identifying the random variable subset within the active
  /// variables (used in all_variables mode; defined from randomVarsKey)
  SizetList randomIndices;
  /// list of indices identifying the non-random variable subset within the
  /// active variables (used in all_variables mode; defined from randomVarsKey)
  SizetList nonRandomIndices;

  /// mean, variance, skewness, and kurtosis (raw, central, standardized,
  /// and offset standardized moments, respectively) computed from the
  /// stochastic expansion form.  For OrthogPolyApproximation, these are
  /// primary, and for InterpPolyApproximation, they are currently inactive.
  RealVector expansionMoments;
  /// mean, variance, skewness, and kurtosis (raw, central, standardized,
  /// and offset standardized moments, respectively) computed via numerical
  /// integration of the response.  For OrthogPolyApproximation, these are
  /// secondary, and for InterpPolyApproximation, they are primary.
  RealVector numericalMoments;

  /// gradient of the primary mean (expansion mean for OrthogPoly,
  /// numerical integration mean for InterpPoly)
  RealVector meanGradient;
  /// gradient of the primary variance (expansion variance for OrthogPoly,
  /// numerical integration variance for InterpPoly)
  RealVector varianceGradient;

  // saved Smolyak coefficients corresponding to savedSmolyakMultiIndex
  //IntArray savedSmolyakCoeffs;
  /// saved trial sets that were computed but not selected
  std::deque<UShortArray> savedSmolyakMultiIndex;

  /// introduce mapping to unify disparate enumeration of sensitivity
  /// indices (e.g. main effects only vs all effects)
  IntIntMap sobolIndexMap;
  /// global sensitivities as given by Sobol'
  RealVector sobolIndices;
  /// total global sensitivities as given by Sobol'
  RealVector totalSobolIndices;

private:

  //
  //- Heading: Data
  //
};


inline PolynomialApproximation::PolynomialApproximation():
  driverRep(NULL), ssgLevelPrev(USHRT_MAX)
{ }


inline PolynomialApproximation::
PolynomialApproximation(size_t num_vars, bool use_derivs):
  BasisApproximation(BaseConstructor(), num_vars), driverRep(NULL),
  ssgLevelPrev(USHRT_MAX)
{ configOptions.useDerivs = use_derivs; }


inline PolynomialApproximation::~PolynomialApproximation()
{ }


inline const RealVector& PolynomialApproximation::expansion_moments() const
{ return expansionMoments; }


inline const RealVector& PolynomialApproximation::numerical_moments() const
{ return numericalMoments; }


inline const SurrogateData& PolynomialApproximation::surrogate_data() const
{ return surrData; }


inline void PolynomialApproximation::surrogate_data(const SurrogateData& data)
{ surrData = data; /* shared representation */ }


//inline size_t PolynomialApproximation::pop_count()
//{
//  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
//  return (size_t)ssg_driver->unique_trial_points();
//}


inline void PolynomialApproximation::solution_approach(short soln_approach)
{ configOptions.expCoeffsSolnApproach = soln_approach; }


inline short PolynomialApproximation::solution_approach() const
{ return configOptions.expCoeffsSolnApproach; }


inline void PolynomialApproximation::expansion_coefficient_flag(bool coeff_flag)
{ configOptions.expansionCoeffFlag = coeff_flag; }


inline bool PolynomialApproximation::expansion_coefficient_flag() const
{ return configOptions.expansionCoeffFlag; }


inline void PolynomialApproximation::
expansion_coefficient_gradient_flag(bool grad_flag)
{ configOptions.expansionCoeffGradFlag = grad_flag; }


inline bool PolynomialApproximation::
expansion_coefficient_gradient_flag() const
{ return configOptions.expansionCoeffGradFlag; }


//inline void PolynomialApproximation::refinement_type(short refine_type)
//{ configOptions.refinementType = refine_type; }


//inline short PolynomialApproximation::refinement_type() const
//{ return configOptions.refinementType; }


inline void PolynomialApproximation::refinement_control(short refine_cntl)
{ configOptions.refinementControl = refine_cntl; }


inline short PolynomialApproximation::refinement_control() const
{ return configOptions.refinementControl; }


inline void PolynomialApproximation::vbd_control(short vbd_cntl)
{ configOptions.vbdControl = vbd_cntl; }


inline short PolynomialApproximation::vbd_control() const
{ return configOptions.vbdControl; }


inline const IntIntMap& PolynomialApproximation::sobol_index_map() const
{ return sobolIndexMap; }


inline const RealVector& PolynomialApproximation::sobol_indices() const
{ return sobolIndices; }


inline const RealVector& PolynomialApproximation::total_sobol_indices() const
{ return totalSobolIndices; }


inline void PolynomialApproximation::
integration_driver_rep(IntegrationDriver* driver_rep)
{ driverRep = driver_rep; }


inline void PolynomialApproximation::
random_variables_key(const BoolDeque& random_vars_key)
{
  randomVarsKey = random_vars_key;
  randomIndices.clear();
  nonRandomIndices.clear();
  for (size_t i=0; i<numVars; i++)
    if (random_vars_key[i])
      randomIndices.push_back(i);
    else
      nonRandomIndices.push_back(i);
}


inline void PolynomialApproximation::
increment_indices(UShortArray& indices, const UShortArray& limits,
		  bool include_limit_equality)
{
  size_t n = indices.size(), increment_index = 0;
  ++indices[increment_index];
  while ( increment_index < n &&
	  ( (  include_limit_equality && 
	       indices[increment_index] >= limits[increment_index] ) ||
	    ( !include_limit_equality && 
	       indices[increment_index] >  limits[increment_index] ) ) ) {
    indices[increment_index] = 0;
    ++increment_index;
    if (increment_index < n)
      ++indices[increment_index];
  }
}


inline void PolynomialApproximation::
increment_terms(UShortArray& terms, size_t& last_index, size_t& prev_index,
		const size_t& term_limit, bool& order_complete)
{
  bool increment_complete = false;
  while (!increment_complete) {
    terms[last_index] = 1;
    ++terms[prev_index];
    if (prev_index == 0) {
      increment_complete = true;
      if (terms[prev_index] > term_limit)
	order_complete = true;
    }
    else {
      last_index = prev_index;
      --prev_index;
      if (terms[last_index] <= terms[prev_index])
	increment_complete = true;
    }
  }
}


inline bool PolynomialApproximation::restore_available()
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  return (std::find(savedSmolyakMultiIndex.begin(),
    savedSmolyakMultiIndex.end(), ssg_driver->trial_index_set()) !=
    savedSmolyakMultiIndex.end());
}


inline size_t PolynomialApproximation::restoration_index()
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  return find_index(savedSmolyakMultiIndex, ssg_driver->trial_index_set());
}


inline size_t PolynomialApproximation::finalization_index(size_t i)
{
  SparseGridDriver*    ssg_driver     = (SparseGridDriver*)driverRep;
  const UShort2DArray& sm_multi_index = ssg_driver->smolyak_multi_index();

  // sm_multi_index is updated in SparseGridDriver::finalize_sets() with all
  // of the remaining trial sets.  Below, we determine the order with which
  // these appended trial sets appear in savedSmolyakMultiIndex
  // (SparseGridDriver::trialSets is a sorted set, but savedSmolyakMultiIndex
  // retains order of insertion).
  size_t num_saved_indices = savedSmolyakMultiIndex.size(),
         start = sm_multi_index.size() - num_saved_indices;
  return find_index(savedSmolyakMultiIndex, sm_multi_index[start+i]);
}

} // namespace Pecos

#endif
