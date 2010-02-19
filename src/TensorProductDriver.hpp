/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 TensorProductDriver
//- Description: 
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef TENSOR_PRODUCT_DRIVER_HPP
#define TENSOR_PRODUCT_DRIVER_HPP

#include "pecos_data_types.hpp"

namespace Pecos {


/// generates N-dimensional tensor-product quadrature grids for
/// numerical evaluation of expectation integrals over independent
/// standard random variables.

/** This class is used by Dakota::NonDQuadrature, but could also be
    used for general numerical integration of moments. */

class TensorProductDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  TensorProductDriver();  ///< default constructor
  ~TensorProductDriver(); ///< destructor

  //
  //- Heading: Member functions
  //

  /// set numVars, isotropicTPQ, and tpqAnisoLevelWts
  void initialize_grid_level(size_t num_vars, const RealVector& dim_pref);

  /// set integrationRules
  void initialize_grid_parameters(const ShortArray& u_types);

  /// compute scaled variable and weight sets for the TPQ grid
  void compute_grid();

  /// number of collocation points
  size_t grid_size();

  /// return weightSets
  const RealVector& weight_sets() const;
  /// return variableSets
  const RealMatrix& variable_sets() const;

  /// return isotropicTPQ
  bool isotropic() const;
  /// set tpqAnisoLevelWts
  void dimension_preference(const RealVector& dim_pref);
  /// set tpqAnisoLevelWts
  void anisotropic_weights(const RealVector& aniso_wts);
  /// return tpqAnisoLevelWts
  const RealVector& anisotropic_weights() const;
  /// return integrationRules
  const IntArray& integration_rules() const;

private:

  //
  //- Heading: Convenience functions
  //

  //
  //- Heading: Data
  //

  /// number of variables in the tensor-product grid
  size_t numVars;

  /// the isotropic/anisotropic quadrature order
  UShortArray quadOrder;

  /// flag indicating an isotropic tensor-product grid
  bool isotropicTPQ;
  // vector of dimension preference levels for anisotropic tensor-product grid
  //RealVector tpqDimPref;
  /// weighting vector for anisotropic tensor-product grids
  RealVector tpqAnisoLevelWts;

  /// integer codes for sgmga routine integration rule options
  IntArray integrationRules;

  /// the set of weights associated with each point in the tensor-product grid
  RealVector weightSets;
  /// the set of points in the tensor-product grid,
  /// arranged num points by numVars
  RealMatrix variableSets;
};


inline TensorProductDriver::TensorProductDriver()
{ }


inline TensorProductDriver::~TensorProductDriver()
{ }


inline size_t TensorProductDriver::grid_size()
{
  size_t i, size = 1;
  for (i=0; i<numVars; ++i)
    size *= quadOrder[i];
  return size;
}


inline const RealVector& TensorProductDriver::weight_sets() const
{ return weightSets; }


inline const RealMatrix& TensorProductDriver::variable_sets() const
{ return variableSets; }


inline bool TensorProductDriver::isotropic() const
{ return isotropicTPQ; }


inline const RealVector& TensorProductDriver::anisotropic_weights() const
{ return tpqAnisoLevelWts; }


inline const IntArray& TensorProductDriver::integration_rules() const
{ return integrationRules; }


inline void TensorProductDriver::
initialize_grid_level(size_t num_vars, const RealVector& dim_pref)
{
  numVars = num_vars;
  dimension_preference(dim_pref);
}

} // namespace Dakota

#endif
