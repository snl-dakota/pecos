/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
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
#include "SurrogateDataPoint.hpp"
#include "SparseGridDriver.hpp"


namespace Pecos {

class IntegrationDriver;
class SparseGridDriver;


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
  friend class OrthogPolyApproximation;

public:

  ConfigurationOptions();
  ConfigurationOptions(short exp_soln_approach, bool exp_coeff_flag,
		       bool exp_grad_flag,
		       //short output_level, short refine_type,
		       short refine_cntl, short vbd_type);
  ~ConfigurationOptions();

private:

  /// identifies the approach taken in compute_coefficients():
  /// QUADRATURE, CUBATURE, SPARSE_GRID, REGRESSION, or SAMPLING
  short expCoeffsSolnApproach;
  /// flag for calculation of expansionCoeffs from response values
  bool expansionCoeffFlag;
  /// flag for calculation of expansionCoeffGrads from response gradients
  bool expansionCoeffGradFlag;

  // output verbosity level: {SILENT,QUIET,NORMAL,VERBOSE,DEBUG}_OUTPUT
  //short outputLevel;
  // nesting override options: NO_OVERRIDE, NESTED, NON_NESTED
  //short nestingOverride;
  // type of refinement: {NO,UNIFORM_P,ADAPTIVE_P}_REFINEMENT
  //short refinementType;

  /// approach for control of refinement: DEFAULT_CONTROL,
  /// TOTAL_SOBOL, SPECTRAL_DECAY, or GENERALIZED_SPARSE
  short refinementControl;
  /// control for amount of data computed in variance-based decomposition:
  /// {NO,UNIVARIATE,ALL}_VBD
  short vbdControl;
};


inline ConfigurationOptions::ConfigurationOptions():
  expCoeffsSolnApproach(SAMPLING), expansionCoeffFlag(true),
  expansionCoeffGradFlag(false),
  //outputLevel(NORMAL_OUTPUT), refinementType(NO_REFINEMENT),
  refinementControl(DEFAULT_CONTROL), vbdControl(NO_VBD)
{ }


inline ConfigurationOptions::
ConfigurationOptions(short exp_soln_approach, bool exp_coeff_flag,
		     bool exp_grad_flag,//short output_level, short refine_type,
		     short refine_cntl, short vbd_type):
  expCoeffsSolnApproach(exp_soln_approach), expansionCoeffFlag(exp_coeff_flag),
  expansionCoeffGradFlag(exp_grad_flag),
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
  PolynomialApproximation(size_t num_vars);
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

  /// return the mean of the expansion, treating all variables as random
  virtual const Real& get_mean() = 0;
  /// return the mean of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  virtual const Real& get_mean(const RealVector& x) = 0;
  /// return the gradient of the expansion mean for a given parameter
  /// vector, treating all variables as random
  virtual const RealVector& get_mean_gradient() = 0;
  /// return the gradient of the expansion mean for a given parameter vector
  /// and given DVV, treating a subset of the variables as random
  virtual const RealVector& get_mean_gradient(const RealVector& x,
					      const SizetArray& dvv) = 0;

  /// return moments 1-4
  virtual const RealVector& get_numeric_moments() = 0;

  /// return the variance of the expansion, treating all variables as random
  virtual const Real& get_variance() = 0;
  /// return the variance of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  virtual const Real& get_variance(const RealVector& x) = 0;
  /// return the gradient of the expansion variance for a given parameter
  /// vector, treating all variables as random
  virtual const RealVector& get_variance_gradient() = 0;
  /// return the gradient of the expansion variance for a given parameter
  /// vector and given DVV, treating a subset of the variables as random
  virtual const RealVector& get_variance_gradient(const RealVector& x,
						  const SizetArray& dvv) = 0;

  /// return the variance of the expansion, treating all variables as random
  virtual const Real& get_covariance(const RealVector& exp_coeffs_2) = 0;
  // return the variance of the expansion for a given parameter vector,
  // treating a subset of the variables as random
  //virtual const Real& get_covariance(const RealVector& x,
  //                                   const RealVector& exp_coeffs_2) = 0;

  //
  //- Heading: Member functions
  //

  /// size component Sobol arrays
  void allocate_component_effects();
  /// size total Sobol arrays
  void allocate_total_effects();

  /// set dataPoints
  void data_points(const std::list<SurrogateDataPoint>& pts);
  /// return size of dataPoints
  size_t data_size() const;
  /// set anchorPoint
  void anchor_point(const SurrogateDataPoint& pt);
  /// queries the status of anchorPoint
  bool anchor() const;

  /// push incoming pt onto the end of dataPoints (implemented at this
  /// intermediate level since dataPoints are not defined at base level)
  void push_back(const SurrogateDataPoint& pt);
  /// pop pop_count() instances off the end of dataPoints (implemented at this
  /// intermediate level since dataPoints are not defined at base level)
  void pop();
  /// pop num_pop_pts instances off the end of dataPoints (implemented at this
  /// intermediate level since dataPoints are not defined at base level)
  void pop(size_t num_pop_pts);
  /// number of data points to remove in a decrement (implemented at this
  /// intermediate level since dataPoints are not defined at base level)
  size_t pop_count();

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

  const RealVector& approximation_coefficients() const;
  void approximation_coefficients(const RealVector& approx_coeffs);

  /// query trial index for presence in savedSmolyakMultiIndex, indicating
  /// the ability to restore a previously evaluated increment
  bool restore_available();
  /// index of the data set to be restored from within saved bookkeeping
  /// (i.e.,savedSmolyakMultiIndex)
  size_t restoration_index();

  //
  //- Heading: Member functions
  //

  //
  //- Heading: Data
  //

  /// a special sample (often at the center of the approximation region)
  /// for which exact matching is enforced (e.g., using equality-constrained
  /// least squares regression).
  SurrogateDataPoint anchorPoint;
  /// set of samples used to build the approximation.  These sample points
  /// are fit approximately (e.g., using least squares regression); exact
  /// matching is not enforced.
  std::vector<SurrogateDataPoint> dataPoints;

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

  /// expected value of the expansion
  Real expansionMean;
  /// gradient of the expected value of the expansion
  RealVector expansionMeanGrad;
  /// variance of the expansion
  Real expansionVariance;
  /// gradient of the variance of the expansion
  RealVector expansionVarianceGrad;
  /// approximation to moments (mean, variance, skewness, kurtosis)
  /// via quadrature 
  RealVector numericMoments;

  /// the coefficients of the expansion
  RealVector expansionCoeffs;
  /// the gradients of the expansion coefficients
  /** may be interpreted as either the gradients of the expansion
      coefficients or the coefficients of expansions for the response
      gradients.  This array is used when sensitivities of moments are
      needed with respect to variables that do not appear in the
      expansion (e.g., with respect to design variables for an
      expansion only over the random variables). */
  RealMatrix expansionCoeffGrads;

  /// saved trial sets that were computed but not selected
  UShort2DArray savedSmolyakMultiIndex;

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


inline PolynomialApproximation::PolynomialApproximation(size_t num_vars):
  BasisApproximation(BaseConstructor(), num_vars), driverRep(NULL),
  ssgLevelPrev(USHRT_MAX)
{ }


inline PolynomialApproximation::~PolynomialApproximation()
{ }


inline void PolynomialApproximation::
data_points(const std::list<SurrogateDataPoint>& pts)
{
  size_t i, num_pts = pts.size();
  dataPoints.resize(num_pts);
  std::list<SurrogateDataPoint>::const_iterator cit;
  for (i=0, cit=pts.begin(); i<num_pts; ++i, ++cit)
    dataPoints[i] = *cit; 
}


inline size_t PolynomialApproximation::data_size() const
{ return dataPoints.size(); }


inline void PolynomialApproximation::anchor_point(const SurrogateDataPoint& pt)
{ anchorPoint = pt; }


inline bool PolynomialApproximation::anchor() const
{ return !anchorPoint.is_null(); }


inline void PolynomialApproximation::push_back(const SurrogateDataPoint& pt)
{ dataPoints.push_back(pt); }


inline void PolynomialApproximation::pop()
{ pop(pop_count()); }


inline void PolynomialApproximation::pop(size_t num_pop_pts)
{
  if (num_pop_pts) {
    size_t data_size = dataPoints.size();
    if (data_size >= num_pop_pts)
      dataPoints.resize(data_size - num_pop_pts);
    else {
      PCerr << "Error: pop count (" << num_pop_pts << ") exceeds data size ("
	    << data_size << ") in PolynomialApproximation::pop(size_t)."
	    << std::endl;
      abort_handler(-1);
    }
  }
}


inline size_t PolynomialApproximation::pop_count()
{
  SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
  return (size_t)ssg_driver->unique_trial_points();
}


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


inline const RealVector& PolynomialApproximation::
approximation_coefficients() const
{ return expansionCoeffs; }


inline void PolynomialApproximation::
approximation_coefficients(const RealVector& approx_coeffs)
{ expansionCoeffs = approx_coeffs; }


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
  UShort2DArray::iterator  sit = std::find(savedSmolyakMultiIndex.begin(),
    savedSmolyakMultiIndex.end(), ssg_driver->trial_index_set());
  return (sit == savedSmolyakMultiIndex.end()) ? _NPOS :
    std::distance(savedSmolyakMultiIndex.begin(), sit);
}

} // namespace Pecos

#endif
