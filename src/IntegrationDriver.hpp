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

#ifndef INTEGRATION_DRIVER_HPP
#define INTEGRATION_DRIVER_HPP

#include "pecos_data_types.hpp"
#include "BasisPolynomial.hpp"

namespace Pecos {

class DistributionParams;


/// base class for generating N-dimensional grids for numerical evaluation
/// of expectation integrals over independent standard random variables.

/** This class enables Dakota::NonD{Quadrature,Cubature,SparseGrid}. */

class IntegrationDriver
{
public:

  //
  //- Heading: Constructors, destructor, and operator=
  //

  /// default constructor
  IntegrationDriver();
  /// standard constructor for envelope
  IntegrationDriver(short driver_type);
  /// copy constructor
  IntegrationDriver(const IntegrationDriver& driver);

  /// destructor
  virtual ~IntegrationDriver();

  /// assignment operator
  IntegrationDriver operator=(const IntegrationDriver& driver);

  //
  //- Heading: Virtual functions
  //

  /// set anisoLevelWts
  virtual void dimension_preference(const RealVector& dim_pref);
  /// set anisoLevelWts
  virtual void anisotropic_weights(const RealVector& aniso_wts);

  /// update polynomialBasis with data from dist_params
  virtual void initialize_grid_parameters(const ShortArray& u_types,
					  const DistributionParams& dp);

  /// compute scaled variable and weight sets for the TPQ grid
  virtual void compute_grid();
  /// compute number of collocation points
  virtual int grid_size();

  //
  //- Heading: Member functions
  //

  /// return isotropicTPQ
  bool isotropic() const;
  /// return anisoLevelWts
  const RealVector& anisotropic_weights() const;
  /// return polynomialBasis
  const std::vector<BasisPolynomial>& polynomial_basis() const;

  /// return weightSets
  const RealVector& weight_sets() const;
  /// return variableSets
  const RealMatrix& variable_sets() const;

  /// return orderGenzKeister
  const UShortArray& genz_keister_order()     const;
  /// return precGenzKeister
  const UShortArray& genz_keister_precision() const;

  /// returns approxRep for access to derived class member functions
  /// that are not mapped to the base level
  IntegrationDriver* driver_rep() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  IntegrationDriver(BaseConstructor);

  //
  //- Heading: Member functions
  //

  /// set int_rules and growth_rules from u_types, nested_rules, growth_rate,
  /// and nested_uniform_rule
  void initialize_rules(const ShortArray& u_types, bool nested_rules,
			short growth_rate, short nested_uniform_rule,
			IntArray& int_rules, IntArray& growth_rules);
  /// set int_rules and growth_rules from poly_basis and growth_rate
  void initialize_rules(const std::vector<BasisPolynomial>& poly_basis,
			short growth_rate, IntArray& int_rules,
			IntArray& growth_rules);

  /// compute variable and weight sets for a tensor-product grid
  void compute_tensor_grid(const UShortArray& order,
			   UShort2DArray& colloc_key,
			   Real2DArray& pts_1d, Real2DArray& wts_1d);

  //
  //- Heading: Data
  //

  /// number of variables in the tensor-product grid
  size_t numVars;

  /// flag indicating a dimension isotropic grid
  bool dimIsotropic;
  // vector of dimension preference levels for dimension anisotropic grids
  //RealVector dimPref;
  /// weighting vector for dimension anisotropic grids
  RealVector anisoLevelWts;

  /// array of one-dimensional orthogonal polynomials used in
  /// computing Gaussian quadrature points and weights
  std::vector<BasisPolynomial> polynomialBasis;

  /// the set of weights associated with each point in the tensor-product grid
  RealVector weightSets;
  /// the set of points in the tensor-product grid,
  /// arranged num points by numVars
  RealMatrix variableSets;

  /// lookup for set of 1-D Genz-Keister quadrature orders
  static UShortArray orderGenzKeister;
  /// lookup for set of 1-D Genz-Keister integrand precisions
  static UShortArray precGenzKeister;

private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// basisApproxRep to the appropriate derived type.
  IntegrationDriver* get_driver(short driver_type);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  IntegrationDriver* driverRep;
  /// number of objects sharing driverRep
  int referenceCount;
};


inline bool IntegrationDriver::isotropic() const
{ return (driverRep) ? driverRep->dimIsotropic : dimIsotropic; }


inline const RealVector& IntegrationDriver::anisotropic_weights() const
{ return (driverRep) ? driverRep->anisoLevelWts : anisoLevelWts; }


inline const std::vector<BasisPolynomial>& 
IntegrationDriver::polynomial_basis() const
{ return (driverRep) ? driverRep->polynomialBasis : polynomialBasis; }


inline const RealVector& IntegrationDriver::weight_sets() const
{ return (driverRep) ? driverRep->weightSets : weightSets; }


inline const RealMatrix& IntegrationDriver::variable_sets() const
{ return (driverRep) ? driverRep->variableSets : variableSets; }


inline const UShortArray& IntegrationDriver::genz_keister_order() const
{ return (driverRep) ? driverRep->orderGenzKeister : orderGenzKeister; }


inline const UShortArray& IntegrationDriver::genz_keister_precision() const
{ return (driverRep) ? driverRep->precGenzKeister : precGenzKeister; }


inline IntegrationDriver* IntegrationDriver::driver_rep() const
{ return driverRep; }

} // namespace Pecos

#endif
