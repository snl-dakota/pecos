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
  //friend class InterpPolyApproximation;
  //friend class NodalInterpPolyApproximation;
  //friend class HierarchInterpPolyApproximation;
  //friend class OrthogPolyApproximation;

public:

  /// default constructor
  ExpansionConfigOptions();
  /// constructor
  ExpansionConfigOptions(short exp_soln_approach, bool exp_coeff_flag,
			 bool exp_grad_flag,
			 //short output_level, short refine_type,
			 short refine_cntl, short vbd_type);
  /// destructor
  ~ExpansionConfigOptions();

//private:

  /// identifies the approach taken in compute_coefficients(): QUADRATURE,
  /// CUBATURE, COMBINED_SPARSE_GRID, HIERARCHICAL_SPARSE_GRID, REGRESSION,
  /// or SAMPLING
  short expCoeffsSolnApproach;

  /// flag for calculation of expansion coefficients from response values
  bool expansionCoeffFlag;
  /// flag for calculation of gradients of expansion coefficients from
  /// response gradients
  bool expansionCoeffGradFlag;

  // output verbosity level: {SILENT,QUIET,NORMAL,VERBOSE,DEBUG}_OUTPUT
  //short outputLevel;

  // type of refinement: {NO,P,H}_REFINEMENT
  //short refinementType;
  /// approach for control of refinement: {NO,UNIFORM,LOCAL_ADAPTIVE}_CONTROL
  /// or DIMENSION_ADAPTIVE_CONTROL_{SOBOL,DECAY,GENERALIZED}
  short refinementControl;

  /// control for amount of data computed in variance-based decomposition:
  /// {NO,UNIVARIATE,ALL}_VBD
  short vbdControl;
};


inline ExpansionConfigOptions::ExpansionConfigOptions():
  expCoeffsSolnApproach(SAMPLING), expansionCoeffFlag(true),
  expansionCoeffGradFlag(false),
  //outputLevel(NORMAL_OUTPUT), refinementType(NO_REFINEMENT),
  refinementControl(NO_CONTROL), vbdControl(NO_VBD)
{ }


inline ExpansionConfigOptions::
ExpansionConfigOptions(short exp_soln_approach, bool exp_coeff_flag,
		       bool exp_grad_flag,
		       //short output_level, short refine_type,
		       short refine_cntl, short vbd_type):
  expCoeffsSolnApproach(exp_soln_approach),
  expansionCoeffFlag(exp_coeff_flag), expansionCoeffGradFlag(exp_grad_flag),
  //outputLevel(output_level), refinementType(refine_type),
  refinementControl(refine_cntl), vbdControl(vbd_type)
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
  //friend class InterpPolyApproximation;
  //friend class NodalInterpPolyApproximation;
  //friend class HierarchInterpPolyApproximation;
  //friend class OrthogPolyApproximation;

public:

  /// default constructor
  BasisConfigOptions();
  /// constructor
  BasisConfigOptions(bool nested_rules, bool piecewise_basis,
		     bool equidistant_rules, bool use_derivs);
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


inline BasisConfigOptions::~BasisConfigOptions()
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

  /// retrieve the gradient for a response expansion with respect to
  /// all variables included in the polynomial bases using the given
  /// parameter vector and default DVV
  virtual const RealVector& gradient_basis_variables(const RealVector& x) = 0;
  /// retrieve the gradient for a response expansion with respect to
  /// variables included in the polynomial basis for a given parameter
  /// vector and a given DVV subset
  virtual const RealVector& gradient_basis_variables(const RealVector& x,
						     const SizetArray& dvv) = 0;
  /// retrieve the gradient for a response expansion with respect to
  /// all variables not included in the polynomial bases using the
  /// given parameter vector and default DVV
  virtual const RealVector&
    gradient_nonbasis_variables(const RealVector& x) = 0;

  /// retrieve the response value for a stored expansion using the
  /// given parameter vector
  virtual Real stored_value(const RealVector& x) = 0;
  /// retrieve the response gradient for a stored expansion with
  /// respect to all variables included in the polynomial bases;
  /// evaluate for the given parameter vector.
  virtual const RealVector&
    stored_gradient_basis_variables(const RealVector& x) = 0;
  /// retrieve the response gradient for a stored expansion with
  /// respect to all variables not included in the polynomial bases;
  /// evaluate for the given parameter vector.
  virtual const RealVector&
    stored_gradient_nonbasis_variables(const RealVector& x) = 0;

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

  /// estimate expansion coefficient decay rates for each random
  /// variable dimension (OrthogPolyApproximation only)
  virtual const RealVector& dimension_decay_rates();
  /// increment the approximation order (OrthogPolyApproximation only)
  virtual void increment_order();

  //
  //- Heading: Member functions
  //

  // allocate basis_types based on u_types
  //static bool initialize_basis_types(const ShortArray& u_types,
  //				       bool piecewise_basis, bool use_derivs,
  //				       ShortArray& basis_types);
  /// allocate colloc_rules based on u_types and rule options
  static void initialize_collocation_rules(const ShortArray& u_types,
    const BasisConfigOptions& bc_options, ShortArray& colloc_rules);
  /// allocate poly_basis based on basis_types and colloc_rules
  static void initialize_polynomial_basis(const ShortArray& basis_types,
    const ShortArray& colloc_rules, std::vector<BasisPolynomial>& poly_basis);
  /// pass distribution parameters from dp to poly_basis
  static void update_basis_distribution_parameters(const ShortArray& u_types,
    const DistributionParams& dp, std::vector<BasisPolynomial>& poly_basis);

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

  /// set ExpansionConfigOptions::expCoeffsSolnApproach
  void solution_approach(short soln_approach);
  /// get ExpansionConfigOptions::expCoeffsSolnApproach
  short solution_approach() const;

  /// set ExpansionConfigOptions::expansionCoeffFlag
  void expansion_coefficient_flag(bool coeff_flag);
  /// get ExpansionConfigOptions::expansionCoeffFlag
  bool expansion_coefficient_flag() const;

  /// set ExpansionConfigOptions::expansionCoeffGradFlag
  void expansion_coefficient_gradient_flag(bool grad_flag);
  /// get ExpansionConfigOptions::expansionCoeffGradFlag
  bool expansion_coefficient_gradient_flag() const;

  // set ExpansionConfigOptions::refinementType
  //void refinement_type(short refine_type);
  // get ExpansionConfigOptions::refinementType
  //short refinement_type() const;

  /// set ExpansionConfigOptions::refinementControl
  void refinement_control(short refine_cntl);
  /// get ExpansionConfigOptions::refinementControl
  short refinement_control() const;

  /// set ExpansionConfigOptions::vbdControl
  void vbd_control(short vbd_cntl);
  /// get ExpansionConfigOptions::vbdControl
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
			      size_t& prev_index, const size_t& term_limit,
			      bool& order_complete);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// generic base class function mapped to gradient_basis_variables(x)
  const RealVector& gradient(const RealVector& x);

  bool restore_available();
  size_t restoration_index();
  size_t finalization_index(size_t i);

  //
  //- Heading: Member functions
  //

  /// compute central moments of response using type1 numerical integration
  void compute_numerical_moments(const RealVector& coeffs,
				 const RealVector& t1_wts, RealVector& moments);
  /// compute central moments of response using type1/2 numerical integration
  void compute_numerical_moments(const RealVector& t1_coeffs,
				 const RealMatrix& t2_coeffs,
				 const RealVector& t1_wts,
				 const RealMatrix& t2_wts, RealVector& moments);
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

  /// saved trial sets that were computed but not selected
  std::deque<UShortArray> savedLevMultiIndex;

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
{ basisConfigOptions.useDerivs = use_derivs; }


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
{ expConfigOptions.expCoeffsSolnApproach = soln_approach; }


inline short PolynomialApproximation::solution_approach() const
{ return expConfigOptions.expCoeffsSolnApproach; }


inline void PolynomialApproximation::expansion_coefficient_flag(bool coeff_flag)
{ expConfigOptions.expansionCoeffFlag = coeff_flag; }


inline bool PolynomialApproximation::expansion_coefficient_flag() const
{ return expConfigOptions.expansionCoeffFlag; }


inline void PolynomialApproximation::
expansion_coefficient_gradient_flag(bool grad_flag)
{ expConfigOptions.expansionCoeffGradFlag = grad_flag; }


inline bool PolynomialApproximation::
expansion_coefficient_gradient_flag() const
{ return expConfigOptions.expansionCoeffGradFlag; }


//inline void PolynomialApproximation::refinement_type(short refine_type)
//{ expConfigOptions.refinementType = refine_type; }


//inline short PolynomialApproximation::refinement_type() const
//{ return expConfigOptions.refinementType; }


inline void PolynomialApproximation::refinement_control(short refine_cntl)
{ expConfigOptions.refinementControl = refine_cntl; }


inline short PolynomialApproximation::refinement_control() const
{ return expConfigOptions.refinementControl; }


inline void PolynomialApproximation::vbd_control(short vbd_cntl)
{ expConfigOptions.vbdControl = vbd_cntl; }


inline short PolynomialApproximation::vbd_control() const
{ return expConfigOptions.vbdControl; }


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


inline const RealVector& PolynomialApproximation::gradient(const RealVector& x)
{ return gradient_basis_variables(x); }


inline bool PolynomialApproximation::restore_available()
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  return (std::find(savedLevMultiIndex.begin(), savedLevMultiIndex.end(),
		    ssg_driver->trial_set()) != savedLevMultiIndex.end());
}


inline size_t PolynomialApproximation::restoration_index()
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  return find_index(savedLevMultiIndex, ssg_driver->trial_set());
}


inline size_t PolynomialApproximation::finalization_index(size_t i)
{
  SparseGridDriver* sg_driver = (SparseGridDriver*)driverRep;
  const std::set<UShortArray>& trial_sets = sg_driver->computed_trial_sets();
  // {Combined,Hierarch}SparseGridDriver::finalize_sets() updates the grid data
  // with remaining computed trial sets (in sorted order from SparseGridDriver::
  // computedTrialSets).  Below, we determine the order with which these
  // appended trial sets appear in savedLevMultiIndex.
  std::set<UShortArray>::const_iterator cit = trial_sets.begin();
  std::advance(cit, i);
  return find_index(savedLevMultiIndex, *cit);
}

} // namespace Pecos

#endif
