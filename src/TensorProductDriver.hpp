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

#include "IntegrationDriver.hpp"

namespace Pecos {


/// generates N-dimensional tensor-product quadrature grids for
/// numerical evaluation of expectation integrals over independent
/// standard random variables.

/** This class is used by Dakota::NonDQuadrature, but could also be
    used for general numerical integration of moments. */

class TensorProductDriver: public IntegrationDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  TensorProductDriver();  ///< default constructor
  ~TensorProductDriver(); ///< destructor

  //
  //- Heading: Virtual function redefinitions
  //

  /// set anisoLevelWts
  void anisotropic_weights(const RealVector& aniso_wts);

  /// compute scaled variable and weight sets for the TPQ grid
  void compute_grid();
  /// number of collocation points
  int grid_size();

  //
  //- Heading: Member functions
  //

  /// set quadOrder
  void quadrature_order(const UShortArray& quad_order);
  /// set ith entry in quadOrder
  void quadrature_order(unsigned short order, size_t i);
  /// return quadOrder
  const UShortArray& quadrature_order() const;
  /// return ith entry in quadOrder
  unsigned short quadrature_order(size_t i) const;

  /// return integrationRules
  const IntArray& integration_rules() const;

  /// invoke initialize_rules() to set integration and growth rules
  void initialize_grid(const ShortArray& u_types, bool nested_rules,
		       short growth_rate, short nested_uniform_rule);

private:

  //
  //- Heading: Convenience functions
  //

  //
  //- Heading: Data
  //

  /// the isotropic/anisotropic quadrature order
  UShortArray quadOrder;

  /// integer codes for integration rule options
  IntArray integrationRules;
  /// integer codes for growth rule options
  IntArray growthRules;
};


inline TensorProductDriver::TensorProductDriver()
{ }


inline TensorProductDriver::~TensorProductDriver()
{ }


inline void TensorProductDriver::quadrature_order(const UShortArray& quad_order)
{ quadOrder = quad_order; }


inline void TensorProductDriver::
quadrature_order(unsigned short order, size_t i)
{ quadOrder[i] = order; }


inline const UShortArray& TensorProductDriver::quadrature_order() const
{ return quadOrder; }


inline unsigned short TensorProductDriver::quadrature_order(size_t i) const
{ return quadOrder[i]; }


inline const IntArray& TensorProductDriver::integration_rules() const
{ return integrationRules; }


inline int TensorProductDriver::grid_size()
{
  int size = 1;
  for (size_t i=0; i<numVars; ++i)
    size *= quadOrder[i];
  return size;
}

} // namespace Pecos

#endif
