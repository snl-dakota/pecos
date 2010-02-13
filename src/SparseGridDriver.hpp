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

#include "pecos_data_types.hpp"
#include "sandia_rules.H"

namespace Pecos {

/// pointer to a Gauss point or weight evaluation function, matching
/// the prototype required by Dakota/packages/quadrature/sparse_grid
typedef void ( *FPType ) ( int order, int num_params, double* params,
			   double* data );


/// Derived nondeterministic class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis and Gaussian quadrature rules within Smolyak
    sparse grids. */

class SparseGridDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  SparseGridDriver();  ///< default constructor
  ~SparseGridDriver(); ///< destructor

  //
  //- Heading: Member functions
  //

  /// set ssgLevel, isotropicSSG, and ssgAnisoLevelWts
  void initialize_grid_level(size_t num_vars, size_t ssg_level,
			     const RealVector& dimension_pref);

  /// set polyParams, integrationRules, and FPType function pointers
  void initialize_grid_parameters(const ShortArray& u_types,
    const RealVector& nuv_means,      const RealVector& nuv_std_devs,
    const RealVector& nuv_l_bnds,     const RealVector& nuv_u_bnds,
    const RealVector& lnuv_means,     const RealVector& lnuv_std_devs,
    const RealVector& lnuv_lambdas,   const RealVector& lnuv_zetas,
    const RealVector& lnuv_err_facts, const RealVector& lnuv_l_bnds,
    const RealVector& lnuv_u_bnds,    const RealVector& luuv_l_bnds,
    const RealVector& luuv_u_bnds,    const RealVector& tuv_modes,
    const RealVector& tuv_l_bnds,     const RealVector& tuv_u_bnds,
    const RealVector& buv_alphas,     const RealVector& buv_betas,
    const RealVector& gauv_alphas,    const RealVector& guuv_alphas,
    const RealVector& guuv_betas,     const RealVector& fuv_alphas,
    const RealVector& fuv_betas,      const RealVector& wuv_alphas,
    const RealVector& wuv_betas,      const RealVectorArray& hbuv_bin_prs);

  /// compute scaled variable and weight sets for the sparse grid
  void compute_grid();

  /// Use webbur::sgmga_vcn_* functions to compute index sets satisfying the
  /// anisotropic index set constraint, along with their corresponding
  /// coefficients
  void anisotropic_multi_index(Int2DArray& multi_index,
			       RealArray& coeffs) const;

  /// number of collocation points with duplicates removed
  int grid_size();
  /// total number of collocation points including duplicates
  int grid_size_total();

  /// return weightSets
  const RealVector& weight_sets() const;
  /// return variableSets
  const RealMatrix& variable_sets() const;

  /// return isotropicSSG
  bool isotropic() const;
  /// set ssgLevel
  void level(unsigned short ssg_level);
  /// return ssgLevel
  unsigned short level() const;
  /// set ssgAnisoLevelWts
  void dimension_preference(const RealVector& dim_pref);
  /// set ssgAnisoLevelWts
  void anisotropic_weights(const RealVector& aniso_wts);
  /// return ssgAnisoLevelWts
  const RealVector& anisotropic_weights() const;
  /// update axisLowerBounds
  void update_axis_lower_bounds();
  // return duplicateTol
  //const Real& duplicate_tolerance() const;
  /// return integrationRules
  const IntArray& integration_rules() const;
  /// return uniqueIndexMapping
  const IntArray& unique_index_mapping() const;

  /// converts an array of sparse grid levels to an array of
  /// integration orders based on integrationRules
  void level_to_order(size_t i, unsigned short level,
		      unsigned short& order) const;
  /// converts an array of sparse grid levels to an array of
  /// integration orders based on integrationRules
  void level_to_order(const UShortArray& levels, UShortArray& orders) const;

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

private:

  //
  //- Heading: Convenience functions
  //

  //
  //- Heading: Data
  //

  /// number of variables in the sparse grid
  size_t numVars;

  /// the Smolyak sparse grid level
  unsigned short ssgLevel;
  /// flag indicating an isotropic Smolyak sparse grid
  /** sgmga routines are used for anisotropic, but sparse_grid_mixed_growth
      routines are used for isotropic due to reduced computational overhead */
  bool isotropicSSG;
  // vector of dimension preference levels for anisotropic sparse grid
  //RealVector ssgDimPref;
  /// weighting vector for anisotropic sparse grids: higher weight equates
  /// to lower dimension preference due to lb <= |alpha|.|i| <= ub
  RealVector ssgAnisoLevelWts;
  /// refinement constraints that ensure that level/anisotropic weight updates
  /// contain all previous multi-index sets
  RealVector axisLowerBounds;

  /// integer codes for sgmga routine integration rule options
  IntArray integrationRules;
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

  /// the set of weights associated with each point in the sparse grid
  RealVector weightSets;
  /// the set of points in the sparse grid, arranged num points by numVars
  RealMatrix variableSets;

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


inline const RealVector& SparseGridDriver::weight_sets() const
{ return weightSets; }


inline const RealMatrix& SparseGridDriver::variable_sets() const
{ return variableSets; }


inline bool SparseGridDriver::isotropic() const
{ return isotropicSSG; }


inline unsigned short SparseGridDriver::level() const
{ return ssgLevel; }


inline void SparseGridDriver::level(unsigned short ssg_level)
{ ssgLevel = ssg_level; }


inline const RealVector& SparseGridDriver::anisotropic_weights() const
{ return ssgAnisoLevelWts; }


//inline const Real& SparseGridDriver::duplicate_tolerance() const
//{ return duplicateTol; }


inline const IntArray& SparseGridDriver::integration_rules() const
{ return integrationRules; }


inline const IntArray& SparseGridDriver::unique_index_mapping() const
{ return uniqueIndexMapping; }


inline void SparseGridDriver::
initialize_grid_level(size_t num_vars, size_t ssg_level,
		      const RealVector& dim_pref)
{
  numVars = num_vars;
  level(ssg_level);
  dimension_preference(dim_pref);
}


inline void SparseGridDriver::
level_to_order(size_t i, unsigned short level, unsigned short& order) const
{
  int ilevel = level, irule = integrationRules[i], iorder;
  webbur::level_to_order_default(1, &ilevel, &irule, &iorder);
  order = iorder;
}


inline void SparseGridDriver::
level_to_order(const UShortArray& levels, UShortArray& orders) const
{
  size_t i, num_levels = levels.size();
  if (orders.size() != num_levels)
    orders.resize(num_levels);
  for (i=0; i<num_levels; i++)
    level_to_order(i, levels[i], orders[i]);
}

} // namespace Dakota

#endif
