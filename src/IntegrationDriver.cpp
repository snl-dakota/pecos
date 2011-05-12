/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 IntegrationDriver
//- Description: Implementation code for IntegrationDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "IntegrationDriver.hpp"
#include "CubatureDriver.hpp"
#include "SparseGridDriver.hpp"
#include "TensorProductDriver.hpp"
#include "PolynomialApproximation.hpp"

static const char rcsId[]="@(#) $Id: IntegrationDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {

UShortArray IntegrationDriver::orderGenzKeister;
UShortArray IntegrationDriver::precGenzKeister;


/** This constructor is the one which must build the base class data
    for all derived classes.  get_driver() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_driver() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~IntegrationDriver). */
IntegrationDriver::IntegrationDriver(BaseConstructor):
  computeType2Weights(false), driverRep(NULL), referenceCount(1)
{
  if (orderGenzKeister.empty()) {
    orderGenzKeister.resize(5); //orderGenzKeister = { 1, 3, 9, 19, 35 };
    orderGenzKeister[0] =  1; orderGenzKeister[1] =  3; orderGenzKeister[2] = 9;
    orderGenzKeister[3] = 19; orderGenzKeister[4] = 35;
  }
  if (precGenzKeister.empty()) {
    precGenzKeister.resize(5); //precGenzKeister = { 1, 5, 15, 29, 51 }; 
    precGenzKeister[0] =  1; precGenzKeister[1] =  5; precGenzKeister[2] = 15;
    precGenzKeister[3] = 29; precGenzKeister[4] = 51;
  }

#ifdef REFCOUNT_DEBUG
  PCout << "IntegrationDriver::IntegrationDriver(BaseConstructor) called to "
        << "build base class for letter." << std::endl;
#endif
}


/** The default constructor: driverRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
IntegrationDriver::IntegrationDriver(): driverRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "IntegrationDriver::IntegrationDriver() called to build empty "
        << "envelope." << std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_driver, since IntegrationDriver(BaseConstructor)
    builds the actual base class data for the derived basis functions. */
IntegrationDriver::IntegrationDriver(short driver_type):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "IntegrationDriver::IntegrationDriver(short) called to "
        << "instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  driverRep = get_driver(driver_type);
  if ( !driverRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize driverRep to the 
    appropriate derived type. */
IntegrationDriver* IntegrationDriver::get_driver(short driver_type)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_driver(short)." << std::endl;
#endif

  if (driver_type == QUADRATURE)
    return new TensorProductDriver();
  else if (driver_type == CUBATURE)
    return new CubatureDriver();
  else if (driver_type == SPARSE_GRID)
    return new SparseGridDriver();
  else {
    PCerr << "Error: IntegrationDriver type " << driver_type
	  << " not available." << std::endl;
    return NULL;
  }
}


/** Copy constructor manages sharing of driverRep and incrementing
    of referenceCount. */
IntegrationDriver::IntegrationDriver(const IntegrationDriver& driver)
{
  // Increment new (no old to decrement)
  driverRep = driver.driverRep;
  if (driverRep) // Check for an assignment of NULL
    driverRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  PCout << "IntegrationDriver::IntegrationDriver(IntegrationDriver&)"
	<< std::endl;
  if (driverRep)
    PCout << "driverRep referenceCount = " << driverRep->referenceCount
	  << std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old driverRep,
    assigns new driverRep, and increments referenceCount for new
    driverRep. */
IntegrationDriver IntegrationDriver::operator=(const IntegrationDriver& driver)
{
  if (driverRep != driver.driverRep) { // std case: old != new
    // Decrement old
    if (driverRep) // Check for null pointer
      if (--driverRep->referenceCount == 0) 
	delete driverRep;
    // Assign and increment new
    driverRep = driver.driverRep;
    if (driverRep) // Check for an assignment of NULL
      driverRep->referenceCount++;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "IntegrationDriver::operator=(IntegrationDriver&)" << std::endl;
  if (driverRep)
    PCout << "driverRep referenceCount = " << driverRep->referenceCount
	  << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes driverRep
    when referenceCount reaches zero. */
IntegrationDriver::~IntegrationDriver()
{ 
  // Check for NULL pointer 
  if (driverRep) {
    --driverRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "driverRep referenceCount decremented to "
	  << driverRep->referenceCount << std::endl;
#endif
    if (driverRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting driverRep" << std::endl;
#endif
      delete driverRep;
    }
  }
}


void IntegrationDriver::
initialize_grid_parameters(const ShortArray& u_types, 
			   const DistributionParams& dp)
{
  if (driverRep)
    driverRep->initialize_grid_parameters(u_types, dp); // forward to letter
  else // default implementation
    PolynomialApproximation::distribution_parameters(u_types, dp,
						     polynomialBasis);
}


void IntegrationDriver::compute_grid()
{
  if (driverRep)
    driverRep->compute_grid(); // forward to letter
  else {
    PCerr << "Error: compute_grid() not available for this driver type."
	  << std::endl;
    abort_handler(-1);
  }
}


int IntegrationDriver::grid_size()
{
  if (!driverRep) {
    PCerr << "Error: grid_size() not available for this driver type."
	  << std::endl;
    abort_handler(-1);
  }
  return driverRep->grid_size(); // forward to letter
}


/** protected function called only from derived class letters. */
void IntegrationDriver::
initialize_rules(const ShortArray& u_types, bool nested_rules,
		 bool piecewise_basis,      bool equidistant_rules, 
		 bool use_derivs,          short nested_uniform_rule)
{
  numVars = u_types.size();
  ShortArray basis_types;
  PolynomialApproximation::distribution_types(u_types, piecewise_basis,
					      use_derivs, basis_types);
  PolynomialApproximation::distribution_rules(u_types, nested_rules,
					      piecewise_basis,equidistant_rules,
					      nested_uniform_rule, collocRules);
  PolynomialApproximation::distribution_basis(basis_types, collocRules,
					      polynomialBasis);
  for (size_t i=0; i<numVars; i++)
    if (basis_types[i] == HERMITE_INTERP ||
	basis_types[i] == PIECEWISE_CUBIC_INTERP)
      computeType2Weights = true;
}


/** protected function called only from derived class letters. */
void IntegrationDriver::
initialize_rules(const std::vector<BasisPolynomial>& poly_basis)
{
  numVars         = poly_basis.size();
  polynomialBasis = poly_basis; // shallow copy
  collocRules.resize(numVars);
  for (size_t i=0; i<numVars; i++) {
    // update collocRules
    collocRules[i] = poly_basis[i].collocation_rule();
    // define computeType2Weights
    short basis_type = poly_basis[i].basis_type();
    if (basis_type == HERMITE_INTERP || basis_type == PIECEWISE_CUBIC_INTERP)
      computeType2Weights = true;
  }
}


void IntegrationDriver::
compute_tensor_grid(const UShortArray& quad_order, RealMatrix&  variable_sets,
		    RealVector&    t1_weight_sets, RealMatrix&  t2_weight_sets,
		    UShort2DArray& colloc_key,     Real2DArray& pts_1d,
		    Real2DArray&   t1_wts_1d,      Real2DArray& t2_wts_1d)
{
  size_t i, j, k, num_colloc_pts = 1;
  for (i=0; i<numVars; ++i)
    num_colloc_pts *= quad_order[i];
  if (pts_1d.empty())
    pts_1d.resize(numVars);
  if (t1_wts_1d.empty())
    t1_wts_1d.resize(numVars);
  if (computeType2Weights && t2_wts_1d.empty())
    t2_wts_1d.resize(numVars);
  for (i=0; i<numVars; ++i) {
    pts_1d[i]    = polynomialBasis[i].collocation_points(quad_order[i]);
    t1_wts_1d[i] = polynomialBasis[i].type1_collocation_weights(quad_order[i]);
    if (computeType2Weights)
      t2_wts_1d[i]
	= polynomialBasis[i].type2_collocation_weights(quad_order[i]);
  }
  // Tensor-product quadrature: Integral of f approximated by
  // Sum_i1 Sum_i2 ... Sum_in (w_i1 w_i2 ... w_in) f(x_i1, x_i2, ..., x_in)
  // > project 1-D colloc point arrays (of potentially different type and order)
  //   into an n-dimensional stencil
  // > compute and store products of 1-D colloc weights at each point in stencil
  t1_weight_sets.sizeUninitialized(num_colloc_pts);
  if (computeType2Weights)
    t2_weight_sets.shapeUninitialized(numVars, num_colloc_pts);
  variable_sets.shapeUninitialized(numVars, num_colloc_pts);//Teuchos: col major
  colloc_key.resize(num_colloc_pts);
  UShortArray colloc_indices(numVars, 0);
  for (i=0; i<num_colloc_pts; ++i) {
    Real& t1_wt_i = t1_weight_sets[i]; t1_wt_i = 1.;
    Real*    pt_i = variable_sets[i]; // column vector i
    for (j=0; j<numVars; ++j) {
      pt_i[j]  =    pts_1d[j][colloc_indices[j]];
      t1_wt_i *= t1_wts_1d[j][colloc_indices[j]];
    }
    if (computeType2Weights) {
      Real* t2_wt_i = t2_weight_sets[i]; // column vector i
      for (j=0; j<numVars; ++j) {
	Real& t2_wt_ij = t2_wt_i[j]; t2_wt_ij = 1.;
	for (k=0; k<numVars; ++k)
	  t2_wt_ij *= (k==j) ? t2_wts_1d[k][colloc_indices[k]] :
	                       t1_wts_1d[k][colloc_indices[k]];
      }
    }
    colloc_key[i] = colloc_indices;
    // increment the n-dimensional collocation point index set
    if (i != num_colloc_pts-1)
      PolynomialApproximation::increment_indices(colloc_indices,
						 quad_order, true);
  }
}

} // namespace Pecos
