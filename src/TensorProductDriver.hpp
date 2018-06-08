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

  /// default constructor
  TensorProductDriver();
  /// constructor
  TensorProductDriver(const UShortArray& quad_order);
  /// constructor
  TensorProductDriver(const UShortArray& quad_order, const UShortArray& key);
  /// destructor
  ~TensorProductDriver();

  //
  //- Heading: Virtual function redefinitions
  //

  void active_key(const UShortArray& key);
  void clear_keys();
  void clear_inactive();

  void compute_grid(RealMatrix& variable_sets);
  int  grid_size();
  void reinterpolated_tensor_grid(const UShortArray& lev_index,
				  const SizetList& reinterp_indices);
  const UShortArray& maximal_grid() const;

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

  /// determine the lowest quadrature order that provides integrand
  /// exactness at least as great as the specified goal, while
  /// satisfying any nestedness constraints
  void integrand_goal_to_nested_quadrature_order(size_t i,
    unsigned short integrand_goal, unsigned short& nested_quad_order);
  /// determine the lowest quadrature order that provides at least as many
  /// points as the specified goal, while satisfying any nestedness constraints
  void quadrature_goal_to_nested_quadrature_order(size_t i,
    unsigned short quad_goal, unsigned short& nested_quad_order);
  /// update quadOrder and levelIndex from ref_quad_order while
  /// satisfying nested rule constraints
  void nested_quadrature_order(const UShortArray& ref_quad_order);

  /// return type1WeightSets
  const RealVector& type1_weight_sets() const;
  /// return type2WeightSets
  const RealMatrix& type2_weight_sets() const;

  /// return levelIndex[activeKey]
  const UShortArray& level_index() const;
  /// return levelIndex[key]
  const UShortArray& level_index(const UShortArray& key) const;
  /// return collocKey[activeKey]
  const UShort2DArray& collocation_key() const;
  /// return collocKey[key]
  const UShort2DArray& collocation_key(const UShortArray& key) const;

  /// stand-alone initializer of tensor grid settings (except for
  /// distribution params)
  void initialize_grid(const ShortArray& u_types,
		       const ExpansionConfigOptions& ec_options,
		       const BasisConfigOptions& bc_options);
  /// helper initializer of tensor grid settings (except distribution params)
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis);

  /// precompute quadrature rules to the maximum current order for each basis
  /// polynomial (efficiency optimization when rules are expensive to compute)
  void precompute_rules();

private:

  //
  //- Heading: Convenience functions
  //

  /// update {levelInd,collocKey}Iter based on activeKey
  void update_active_iterators();
  
  /// update levelIndex from quadOrder
  void update_level_index_from_quadrature_order();
  /// update levelIndex[i] from quadOrder[i]
  void update_level_index_from_quadrature_order(size_t i);

  /// update quadOrder from levelIndex
  void update_quadrature_order_from_level_index();

  //
  //- Heading: Data
  //

  /// the isotropic/anisotropic quadrature order
  UShortArray quadOrder;

  /// quadrature order offset by one for use as 0-based indices
  std::map<UShortArray, UShortArray> levelIndex;
  /// iterator to active entry within levelIndex
  std::map<UShortArray, UShortArray>::iterator levelIndIter;

  /// num points-by-numVars array for identifying the 1-D point
  /// indices for sets of tensor-product collocation points
  std::map<UShortArray, UShort2DArray> collocKey;
  /// iterator to active entry within levelIndex
  std::map<UShortArray, UShort2DArray>::iterator collocKeyIter;

  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each point in the tensor grid
  std::map<UShortArray, RealVector> type1WeightSets;
  /// the set of type2 weights (for integration of gradient interpolants)
  /// for each derivative component and for each point in the tensor grid
  std::map<UShortArray, RealMatrix> type2WeightSets;
  
  /// database key indicating the currently active integration configuration.
  /// the key is a multi-index managing multiple modeling dimensions such as
  /// model form, discretization level, etc.
  UShortArray activeKey;
};


inline TensorProductDriver::TensorProductDriver():
  IntegrationDriver(BaseConstructor()), levelIndIter(levelIndex.end())
{ update_active_iterators(); } // default activeKey is empty array


inline TensorProductDriver::TensorProductDriver(const UShortArray& quad_order):
  IntegrationDriver(BaseConstructor()), levelIndIter(levelIndex.end())
{
  update_active_iterators(); // default activeKey is empty array
  quadrature_order(quad_order);
}


inline TensorProductDriver::
TensorProductDriver(const UShortArray& quad_order, const UShortArray& key):
  IntegrationDriver(BaseConstructor()), levelIndIter(levelIndex.end()),
  activeKey(key)
{
  update_active_iterators(); // activeKey set in initializer list
  quadrature_order(quad_order);
}


inline TensorProductDriver::~TensorProductDriver()
{ }


inline void TensorProductDriver::active_key(const UShortArray& key)
{
  if (activeKey != key) {
    activeKey = key;
    update_active_iterators();
  }
}


inline void TensorProductDriver::clear_keys()
{
  activeKey.clear();
  levelIndex.clear(); levelIndIter  = levelIndex.end();
  collocKey.clear();  collocKeyIter = collocKey.end();
  type1WeightSets.clear(); type2WeightSets.clear();
}


inline void TensorProductDriver::update_active_iterators()
{
  // Test for change
  if (levelIndIter != levelIndex.end() && levelIndIter->first == activeKey)
    return;

  levelIndIter = levelIndex.find(activeKey);
  if (levelIndIter == levelIndex.end()) {
    std::pair<UShortArray, UShortArray> ua_pair(activeKey, UShortArray());
    levelIndIter = levelIndex.insert(ua_pair).first;
  }
  update_quadrature_order_from_level_index(); // empty for new levelIndex

  collocKeyIter = collocKey.find(activeKey);
  if (collocKeyIter == collocKey.end()) {
    std::pair<UShortArray, UShort2DArray> u2a_pair(activeKey, UShort2DArray());
    collocKeyIter = collocKey.insert(u2a_pair).first;
  }
}


inline void TensorProductDriver::update_level_index_from_quadrature_order()
{
  UShortArray& lev_index = levelIndIter->second;
  size_t i, len = quadOrder.size();
  if (lev_index.size() != len) lev_index.resize(len);
  for (i=0; i<len; ++i)
    lev_index[i] = quadOrder[i] - 1;
}


inline void TensorProductDriver::
update_level_index_from_quadrature_order(size_t i)
{ levelIndIter->second[i] = quadOrder[i] - 1; }


inline void TensorProductDriver::update_quadrature_order_from_level_index()
{
  const UShortArray& lev_index = levelIndIter->second;
  size_t i, len = lev_index.size();
  if (quadOrder.size() != len) quadOrder.resize(len);
  for (i=0; i<len; ++i)
    quadOrder[i] = lev_index[i] + 1;
}


inline void TensorProductDriver::quadrature_order(const UShortArray& quad_order)
{ quadOrder = quad_order; update_level_index_from_quadrature_order(); }


inline void TensorProductDriver::
quadrature_order(unsigned short order, size_t i)
{ quadOrder[i] = order; update_level_index_from_quadrature_order(i); }


inline const UShortArray& TensorProductDriver::quadrature_order() const
{ return quadOrder; }


inline unsigned short TensorProductDriver::quadrature_order(size_t i) const
{ return quadOrder[i]; }


inline void TensorProductDriver::
nested_quadrature_order(const UShortArray& ref_quad_order)
{
  size_t i, len = ref_quad_order.size();
  if (quadOrder.size()            != len)            quadOrder.resize(len);
  if (levelIndIter->second.size() != len) levelIndIter->second.resize(len);
  unsigned short nested_order;
  for (i=0; i<len; ++i) {
    // synchronize on number of points: Lagrange poly order = #pts - 1
    if (driverMode == INTERPOLATION_MODE)
      quadrature_goal_to_nested_quadrature_order(i, ref_quad_order[i],
						 nested_order);
    else // {INTEGRATION,DEFAULT}_MODE: synchronize on integrand prec 2m-1
      integrand_goal_to_nested_quadrature_order(i, 2 * ref_quad_order[i] - 1,
						nested_order);
    quadrature_order(nested_order, i); // sets quadOrder and levelIndex
  }
}


inline const RealVector& TensorProductDriver::type1_weight_sets() const
{
  std::map<UShortArray, RealVector>::const_iterator cit
    = type1WeightSets.find(activeKey);
  if (cit == type1WeightSets.end()) {
    PCerr << "Error: active key not found in TensorProductDriver::"
	  << "()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const RealMatrix& TensorProductDriver::type2_weight_sets() const
{
  std::map<UShortArray, RealMatrix>::const_iterator cit
    = type2WeightSets.find(activeKey);
  if (cit == type2WeightSets.end()) {
    PCerr << "Error: key not found in TensorProductDriver::"
	  << "type2_weight_sets()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const UShortArray& TensorProductDriver::level_index() const
{ return levelIndIter->second; }


inline const UShortArray& TensorProductDriver::
level_index(const UShortArray& key) const
{
  std::map<UShortArray, UShortArray>::const_iterator cit = levelIndex.find(key);
  if (cit == levelIndex.end()) {
    PCerr << "Error: key not found in TensorProductDriver::level_index()."
	  << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const UShort2DArray& TensorProductDriver::collocation_key() const
{ return collocKeyIter->second; }


inline const UShort2DArray& TensorProductDriver::
collocation_key(const UShortArray& key) const
{
  std::map<UShortArray, UShort2DArray>::const_iterator cit
    = collocKey.find(key);
  if (cit == collocKey.end()) {
    PCerr << "Error: key not found in TensorProductDriver::"
	  << "collocation_key()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline int TensorProductDriver::grid_size()
{
  int size = 1;
  for (size_t i=0; i<numVars; ++i)
    size *= quadOrder[i];
  return size;
}

} // namespace Pecos

#endif
