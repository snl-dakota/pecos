/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 SparseGridDriver
//- Description: Wrapper class for C++ code from packages/quadrature/sparse_grid
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef SPARSE_GRID_DRIVER_HPP
#define SPARSE_GRID_DRIVER_HPP

#include "IntegrationDriver.hpp"
#include "sandia_rules.H"

namespace Pecos {

/// pointer to a Gauss point or weight evaluation function, matching
/// the prototype required by Pecos/packages/VPISparseGrid
typedef void ( *FPType ) ( int order, int num_params, double* params,
			   double* data );

class DistributionParams;


/// Derived nondeterministic class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis and Gaussian quadrature rules within Smolyak
    sparse grids. */

class SparseGridDriver: public IntegrationDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  SparseGridDriver();  ///< default constructor
  ~SparseGridDriver(); ///< destructor

  //
  //- Heading: Virtual function redefinitions
  //

  /// compute scaled variable and weight sets for the sparse grid
  void compute_grid();
  /// number of collocation points with duplicates removed
  int grid_size();

  /// set ssgAnisoLevelWts
  void dimension_preference(const RealVector& dim_pref);
  /// set ssgAnisoLevelWts
  void anisotropic_weights(const RealVector& aniso_wts);

  //
  //- Heading: Member functions
  //

  /// set ssgLevel, isotropicSSG, and ssgAnisoLevelWts
  void initialize_grid(const ShortArray& u_types, size_t ssg_level,
		       const RealVector& dim_pref, const String& ssg_usage,
		       short exp_growth, short nested_uniform_rule);

  /// set polyParams, integrationRules, and FPType function pointers
  void initialize_grid_parameters(const ShortArray& u_types,
				  const DistributionParams& dp);

  /// Use webbur::sgmga_vcn_* functions to compute index sets satisfying
  /// the anisotropic index set constraint, along with their corresponding
  /// coefficients
  void anisotropic_multi_index(Int2DArray& multi_index,
			       RealArray& coeffs) const;

  /// total number of collocation points including duplicates
  int grid_size_total();

  /// set ssgLevel
  void level(unsigned short ssg_level);
  /// return ssgLevel
  unsigned short level() const;

  /// return integrationRules
  const IntArray& integration_rules() const;
  /// return growthRules
  const IntArray& growth_rules() const;

  /// update axisLowerBounds
  void update_axis_lower_bounds();

  // return duplicateTol
  //const Real& duplicate_tolerance() const;

  /// return uniqueIndexMapping
  const IntArray& unique_index_mapping() const;

  /// converts an array of sparse grid levels to an array of
  /// integration orders based on integrationRules
  void level_to_order(size_t i, unsigned short level,
		      unsigned short& order);
  /// converts an array of sparse grid levels to an array of
  /// integration orders based on integrationRules
  void level_to_order(const UShortArray& levels, UShortArray& orders);

private:

  //
  //- Heading: Convenience functions
  //

  /// function for numerically-generated Gauss points for
  /// BOUNDED_NORMAL distribution
  static void bounded_normal_gauss_points(int order, int num_params,
					  double* params, double* data);
  /// function for numerically-generated Gauss weights for
  /// BOUNDED_NORMAL distribution
  static void bounded_normal_gauss_weights(int order, int num_params,
					   double* params, double* data);
  /// function for numerically-generated Gauss points for LOGNORMAL distribution
  static void lognormal_gauss_points(int order, int num_params, double* params,
				     double* data);
  /// function for numerically-generated Gauss weights for LOGNORMAL
  /// distribution
  static void lognormal_gauss_weights(int order, int num_params, double* params,
				      double* data);
  /// function for numerically-generated Gauss points for
  /// BOUNDED_LOGNORMAL distribution
  static void bounded_lognormal_gauss_points(int order, int num_params,
					     double* params, double* data);
  /// function for numerically-generated Gauss weights for
  /// BOUNDED_LOGNORMAL distribution
  static void bounded_lognormal_gauss_weights(int order, int num_params,
					      double* params, double* data);
  /// function for numerically-generated Gauss points for LOGUNIFORM
  /// distribution
  static void loguniform_gauss_points(int order, int num_params, double* params,
				      double* data);
  /// function for numerically-generated Gauss weights for LOGUNIFORM
  /// distribution
  static void loguniform_gauss_weights(int order, int num_params,
				       double* params, double* data);
  /// function for numerically-generated Gauss points for TRIANGULAR
  /// distribution
  static void triangular_gauss_points(int order, int num_params, double* params,
				      double* data);
  /// function for numerically-generated Gauss weights for TRIANGULAR
  /// distribution
  static void triangular_gauss_weights(int order, int num_params,
				       double* params, double* data);
  /// function for numerically-generated Gauss points for GUMBEL distribution
  static void gumbel_gauss_points(int order, int num_params, double* params,
				  double* data);
  /// function for numerically-generated Gauss weights for GUMBEL distribution
  static void gumbel_gauss_weights(int order, int num_params, double* params,
				   double* data);
  /// function for numerically-generated Gauss points for FRECHET distribution
  static void frechet_gauss_points(int order, int num_params, double* params,
				   double* data);
  /// function for numerically-generated Gauss weights for FRECHET distribution
  static void frechet_gauss_weights(int order, int num_params, double* params,
				    double* data);
  /// function for numerically-generated Gauss points for WEIBULL distribution
  static void weibull_gauss_points(int order, int num_params, double* params,
				   double* data);
  /// function for numerically-generated Gauss weights for WEIBULL distribution
  static void weibull_gauss_weights(int order, int num_params, double* params,
				    double* data);
  /// function for numerically-generated Gauss points for
  /// HISTOGRAM_BIN distribution
  static void histogram_bin_gauss_points(int order, int num_params,
					 double* params, double* data);
  /// function for numerically-generated Gauss weights for
  /// HISTOGRAM_BIN distribution
  static void histogram_bin_gauss_weights(int order, int num_params,
					  double* params, double* data);

  //
  //- Heading: Data
  //

  /// the Smolyak sparse grid level
  unsigned short ssgLevel;

  /// integer codes for integration rule options
  IntArray integrationRules;
  /// integer codes for growth rule options
  IntArray growthRules;

  /// "integration" or "interpolation"
  String sparseGridUsage;

  /// refinement constraints that ensure that level/anisotropic weight updates
  /// contain all previous multi-index sets
  RealVector axisLowerBounds;

  /// number of parameters for each polynomial; input to sgmga routines
  /// (corresponds to set of variables defined by integrationRules)
  IntArray numPolyParams;
  /// concatenated array of polynomial parameters for input to sgmga routines
  /// (corresponds to set of variables defined by integrationRules)
  RealArray polyParams;
  /// duplication tolerance used in sgmga routines
  Real duplicateTol;

  // maps indices and bases from sgmga_index() to collocation point index
  //IntArraySizetMap ssgIndexMap;

  /// output from sgmga_unique_index()
  IntArray uniqueIndexMapping;

  /// array of pointers to Gauss point evaluation functions
  std::vector<FPType> compute1DPoints;
  /// array of pointers to Gauss weight evaluation functions
  std::vector<FPType> compute1DWeights;
};


inline SparseGridDriver::SparseGridDriver(): duplicateTol(1.e-15)
{ }


inline SparseGridDriver::~SparseGridDriver()
{ }


inline unsigned short SparseGridDriver::level() const
{ return ssgLevel; }


inline void SparseGridDriver::level(unsigned short ssg_level)
{ ssgLevel = ssg_level; }


inline const IntArray& SparseGridDriver::integration_rules() const
{ return integrationRules; }


inline const IntArray& SparseGridDriver::growth_rules() const
{ return growthRules; }


//inline const Real& SparseGridDriver::duplicate_tolerance() const
//{ return duplicateTol; }


inline const IntArray& SparseGridDriver::unique_index_mapping() const
{ return uniqueIndexMapping; }


inline void SparseGridDriver::
level_to_order(size_t i, unsigned short level, unsigned short& order)
{
  int ilevel = level, iorder;
  webbur::level_growth_to_order(1, &ilevel, &integrationRules[i],
				&growthRules[i], &iorder);
  order = iorder;
}


inline void SparseGridDriver::
level_to_order(const UShortArray& levels, UShortArray& orders)
{
  size_t i, num_levels = levels.size();
  if (orders.size() != num_levels)
    orders.resize(num_levels);
  for (i=0; i<num_levels; i++)
    level_to_order(i, levels[i], orders[i]);
}

} // namespace Pecos

#endif
