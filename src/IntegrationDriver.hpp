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

  /// compute scaled variable and weight sets for the TPQ grid
  virtual void compute_grid();
  /// compute number of collocation points
  virtual int grid_size();

  //
  //- Heading: Member functions
  //

  /// return isotropicTPQ
  bool isotropic() const;
  /// return integrationRules
  const IntArray& integration_rules() const;
  /// return anisoLevelWts
  const RealVector& anisotropic_weights() const;
  /// return polynomialBasis
  const std::vector<BasisPolynomial>& polynomial_basis() const;

  /// return weightSets
  const RealVector& weight_sets() const;
  /// return variableSets
  const RealMatrix& variable_sets() const;

  /// return gaussPts1D
  const Real3DArray& gauss_points_array()  const;
  /// return gaussWts1D
  const Real3DArray& gauss_weights_array() const;

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

  /// integer codes for sgmga routine integration rule options
  IntArray integrationRules;

  /// array of one-dimensional orthogonal polynomials used in
  /// computing Gaussian quadrature points and weights
  std::vector<BasisPolynomial> polynomialBasis;

  /// the set of weights associated with each point in the tensor-product grid
  RealVector weightSets;
  /// the set of points in the tensor-product grid,
  /// arranged num points by numVars
  RealMatrix variableSets;

  /// numContinuousVars x num_levels_per_var sets of 1D Gauss points
  Real3DArray gaussPts1D;
  /// numContinuousVars x num_levels_per_var sets of 1D Gauss weights
  Real3DArray gaussWts1D;

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


inline const IntArray& IntegrationDriver::integration_rules() const
{ return (driverRep) ? driverRep->integrationRules : integrationRules; }


inline const std::vector<BasisPolynomial>& IntegrationDriver::
polynomial_basis() const
{ return (driverRep) ? driverRep->polynomialBasis : polynomialBasis; }


inline const RealVector& IntegrationDriver::weight_sets() const
{ return (driverRep) ? driverRep->weightSets : weightSets; }


inline const RealMatrix& IntegrationDriver::variable_sets() const
{ return (driverRep) ? driverRep->variableSets : variableSets; }


inline const Real3DArray& IntegrationDriver::gauss_points_array() const
{ return (driverRep) ? driverRep->gaussPts1D : gaussPts1D; }


inline const Real3DArray& IntegrationDriver::gauss_weights_array() const
{ return (driverRep) ? driverRep->gaussWts1D : gaussWts1D; }


inline IntegrationDriver* IntegrationDriver::driver_rep() const
{ return driverRep; }

} // namespace Pecos

#endif
