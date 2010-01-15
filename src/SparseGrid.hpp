
//- Class:	 SparseGrid
//- Description: Wrapper class for C++ code from packages/quadrature/sparse_grid
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef SPARSE_GRID_HPP
#define SPARSE_GRID_HPP

#include "pecos_global_defs.hpp"
#include "pecos_data_types.hpp"

namespace Pecos {


/// Derived nondeterministic class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis and Gaussian quadrature rules within Smolyak
    sparse grids. */

class SparseGrid
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  SparseGrid();
  /// destructor
  ~SparseGrid();

  //
  //- Heading: Member functions
  //

  /// converts sparse grid level to integration order for closed rules
  /// with exponential growth
  void level_to_order_closed_exponential(unsigned short level,
					 unsigned short& order) const;
  /// converts sparse grid level to integration order for open rules
  /// with exponential growth
  void level_to_order_open_exponential(unsigned short level,
				       unsigned short& order) const;
  /// converts sparse grid level to integration order for open rules
  /// with linear growth
  void level_to_order_open_linear(unsigned short level,
				  unsigned short& order) const;

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

};


inline SparseGrid::SparseGrid()
{ }


inline SparseGrid::~SparseGrid()
{ }

// NOTE: these mappings are intended for 0-based level indices,
// i.e. the j index set

/** Adapted from webbur::level_to_order_default() for DAKOTA data types. */
inline void SparseGrid::
level_to_order_closed_exponential(unsigned short level,
				  unsigned short& order) const
{ order = (level) ? (unsigned short)std::pow(2., (int)level) + 1 : 1; }


/** Adapted from webbur::level_to_order_default() for DAKOTA data types. */
inline void SparseGrid::
level_to_order_open_exponential(unsigned short level,
				unsigned short& order) const
{ order = (unsigned short)std::pow(2., (int)level + 1) - 1; }


/** Adapted from webbur::level_to_order_default() for DAKOTA data types. */
inline void SparseGrid::
level_to_order_open_linear(unsigned short level, unsigned short& order) const
{ order = 2 * level + 1; }

} // namespace Dakota

#endif
