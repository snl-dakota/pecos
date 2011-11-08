/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
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

  /// return levelIndex
  const UShortArray& level_index() const;

  /// return collocKey
  const UShort2DArray& collocation_key() const;

  /// invoke initialize_rules() to set collocation rules
  void initialize_grid(const ShortArray& u_types, bool nested_rules = false,
		       bool piecewise_basis = false,
		       bool equidistant_rules = true, bool use_derivs = false);
  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis,
		       const UShortArray& quad_order);

private:

  //
  //- Heading: Convenience functions
  //

  /// update levelIndex from quadOrder
  void update_level_index();
  /// update levelIndex[i] from quadOrder[i]
  void update_level_index(size_t i);

  //
  //- Heading: Data
  //

  /// the isotropic/anisotropic quadrature order
  UShortArray quadOrder;

  /// quadrature order offset by one for use as 0-based indices
  UShortArray levelIndex;

  /// numCollocPts-by-numVars array for identifying the 1-D point
  /// indices for sets of tensor-product collocation points
  UShort2DArray collocKey;
};


inline TensorProductDriver::TensorProductDriver():
  IntegrationDriver(BaseConstructor())
{ }


inline TensorProductDriver::~TensorProductDriver()
{ }


inline void TensorProductDriver::update_level_index()
{
  size_t i, len = quadOrder.size();
  if (levelIndex.size() != len) levelIndex.resize(len);
  for (i=0; i<len; ++i)
    levelIndex[i] = quadOrder[i] - 1;
}


inline void TensorProductDriver::update_level_index(size_t i)
{ levelIndex[i] = quadOrder[i] - 1; }


inline void TensorProductDriver::quadrature_order(const UShortArray& quad_order)
{ quadOrder = quad_order; update_level_index(); }


inline void TensorProductDriver::
quadrature_order(unsigned short order, size_t i)
{ quadOrder[i] = order; update_level_index(i); }


inline const UShortArray& TensorProductDriver::quadrature_order() const
{ return quadOrder; }


inline unsigned short TensorProductDriver::quadrature_order(size_t i) const
{ return quadOrder[i]; }


inline const UShortArray& TensorProductDriver::level_index() const
{ return levelIndex; }


inline const UShort2DArray& TensorProductDriver::collocation_key() const
{ return collocKey; }


inline int TensorProductDriver::grid_size()
{
  int size = 1;
  for (size_t i=0; i<numVars; ++i)
    size *= quadOrder[i];
  return size;
}

} // namespace Pecos

#endif
