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

  /// compute scaled variable and weight sets for the TPQ grid
  void compute_grid(RealMatrix& variable_sets);
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

  /// return collocKey
  const UShort2DArray& collocation_key() const;

  /// return gaussPts1D
  const Real2DArray& gauss_points_array()  const;
  /// return gaussWts1D
  const Real2DArray& gauss_weights_array() const;

  /// invoke initialize_rules() to set integration and growth rules
  void initialize_grid(const ShortArray& u_types, bool nested_rules = false,
		       short growth_rate = MODERATE_RESTRICTED_GROWTH,
		       short nested_uniform_rule = GAUSS_PATTERSON);
  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis,
		       const UShortArray& quad_order,
		       short growth_rate = MODERATE_RESTRICTED_GROWTH);

private:

  //
  //- Heading: Convenience functions
  //

  //
  //- Heading: Data
  //

  /// the isotropic/anisotropic quadrature order
  UShortArray quadOrder;

  /// numCollocPts-by-numVars array for identifying the 1-D point
  /// indices for sets of tensor-product collocation points
  UShort2DArray collocKey;

  /// numContinuousVars sets of 1D Gauss points
  Real2DArray gaussPts1D;
  /// numContinuousVars sets of 1D Gauss weights
  Real2DArray gaussWts1D;
};


inline TensorProductDriver::TensorProductDriver():
  IntegrationDriver(BaseConstructor())
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


inline const UShort2DArray& TensorProductDriver::collocation_key() const
{ return collocKey; }


inline const Real2DArray& TensorProductDriver::gauss_points_array() const
{ return gaussPts1D; }


inline const Real2DArray& TensorProductDriver::gauss_weights_array() const
{ return gaussWts1D; }


inline int TensorProductDriver::grid_size()
{
  int size = 1;
  for (size_t i=0; i<numVars; ++i)
    size *= quadOrder[i];
  return size;
}

} // namespace Pecos

#endif
