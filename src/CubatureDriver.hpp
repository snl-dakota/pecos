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

#include "pecos_data_types.hpp"

namespace Pecos {


/// Generates N-dimensional cubature grids for numerical evaluation of
/// expectation integrals over independent standard random variables.

/** Includes Stroud rules and extensions.  This class is used by
    Dakota::NonDCubature, but could also be used for general numerical
    integration of moments. */

class CubatureDriver
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

  /// set numVars
  void initialize(size_t num_vars, unsigned short prec, int rule);

  /// number of collocation points with duplicates removed
  int grid_size();
  /// compute scaled variable and weight sets for the cubature grid
  void compute_grid();

  /// return weightSets
  const RealVector& weight_sets() const;
  /// return variableSets
  const RealMatrix& variable_sets() const;

private:

  //
  //- Heading: Convenience functions
  //

  //
  //- Heading: Data
  //

  /// number of variables in the cubature grid
  size_t numVars;
  /// integrand precision
  unsigned short integrandPrec;
  /// integer code for integration rule
  int integrationRule;

  /// the set of weights associated with each point in the cubature grid
  RealVector weightSets;
  /// the set of points in the cubature grid, arranged num points by numVars
  RealMatrix variableSets;
};


inline CubatureDriver::CubatureDriver()
{ }


inline CubatureDriver::~CubatureDriver()
{ }


inline const RealVector& CubatureDriver::weight_sets() const
{ return weightSets; }


inline const RealMatrix& CubatureDriver::variable_sets() const
{ return variableSets; }


inline void CubatureDriver::
initialize(size_t num_vars, unsigned short prec, int rule)
{ numVars = num_vars; integrandPrec = prec; integrationRule = rule; }

} // namespace Dakota

#endif
