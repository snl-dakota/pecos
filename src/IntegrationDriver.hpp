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

  /// return polynomialBasis
  const std::vector<BasisPolynomial>& polynomial_basis() const;

  /// return weightSets
  const RealVector& weight_sets() const;

  /// return collocRules
  const ShortArray& collocation_rules() const;

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
			bool  equidistant_rules,  short nested_uniform_rule);
  /// set int_rules and growth_rules from poly_basis and growth_rate
  void initialize_rules(const std::vector<BasisPolynomial>& poly_basis);

  /// compute variable and weight sets for a tensor-product grid
  void compute_tensor_grid(const UShortArray& order, RealMatrix& variable_sets,
			   RealVector& weight_sets, UShort2DArray& colloc_key,
			   Real2DArray& pts_1d, Real2DArray& wts_1d);

  //
  //- Heading: Data
  //

  /// number of variables in the tensor-product grid
  size_t numVars;

  /// enumeration codes for integration rule options
  ShortArray collocRules;

  /// array of one-dimensional orthogonal polynomials used in
  /// computing Gaussian quadrature points and weights
  std::vector<BasisPolynomial> polynomialBasis;

  /// the set of weights associated with each point in the {TPQ,SSG,Cub} grid
  RealVector weightSets;

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


inline const std::vector<BasisPolynomial>& 
IntegrationDriver::polynomial_basis() const
{ return (driverRep) ? driverRep->polynomialBasis : polynomialBasis; }


inline const RealVector& IntegrationDriver::weight_sets() const
{ return (driverRep) ? driverRep->weightSets : weightSets; }


inline const ShortArray& IntegrationDriver::collocation_rules() const
{ return collocRules; }


inline const UShortArray& IntegrationDriver::genz_keister_order() const
{ return (driverRep) ? driverRep->orderGenzKeister : orderGenzKeister; }


inline const UShortArray& IntegrationDriver::genz_keister_precision() const
{ return (driverRep) ? driverRep->precGenzKeister : precGenzKeister; }


inline IntegrationDriver* IntegrationDriver::driver_rep() const
{ return driverRep; }

} // namespace Pecos

#endif
