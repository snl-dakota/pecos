/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 CubatureDriver
//- Description: Wrapper class for cubature components within VPISparseGrid
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef CUBATURE_DRIVER_HPP
#define CUBATURE_DRIVER_HPP

#include "IntegrationDriver.hpp"

namespace Pecos {

class DistributionParams;


/// Generates N-dimensional cubature grids for numerical evaluation of
/// expectation integrals over independent standard random variables.

/** Includes Stroud rules and extensions.  This class is used by
    Dakota::NonDCubature, but could also be used for general numerical
    integration of moments. */

class CubatureDriver: public IntegrationDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  CubatureDriver();  ///< default constructor
  ~CubatureDriver(); ///< destructor

  //
  //- Heading: Member functions
  //

  /// set numVars, integrandOrder, and integrationRule
  void initialize_grid(const ShortArray& u_types, unsigned short order,
		       unsigned short rule);
  /// initialize settings for parameterized cubature rules
  void initialize_grid_parameters(const ShortArray& u_types,
				  const DistributionParams& dist_params);

  // set integrandOrder
  void integrand_order(unsigned short order);
  // get integrandOrder
  unsigned short integrand_order() const;

  // set integrationRule
  void integration_rule(unsigned short rule);
  // get integrationRule
  unsigned short integration_rule() const;

  /// number of collocation points with duplicates removed
  int grid_size();
  /// compute scaled variable and weight sets for the cubature grid
  void compute_grid();

private:

  //
  //- Heading: Convenience functions
  //

  //
  //- Heading: Data
  //

  /// integrand order
  unsigned short integrandOrder;
  /// integer code for integration rule
  unsigned short integrationRule;

  /// alpha value for parameterized polynomials: Jacobi, generalized Laguerre
  Real alphaStat;
  /// beta value for parameterized polynomials: Jacobi, generalized Laguerre
  Real betaStat;
};


inline CubatureDriver::CubatureDriver()
{ }


inline CubatureDriver::~CubatureDriver()
{ }


inline void CubatureDriver::integrand_order(unsigned short order)
{ integrandOrder = order; }


inline unsigned short CubatureDriver::integrand_order() const
{ return integrandOrder; }


inline void CubatureDriver::integration_rule(unsigned short rule)
{ integrationRule = rule; }


inline unsigned short CubatureDriver::integration_rule() const
{ return integrationRule; }

} // namespace Pecos

#endif
