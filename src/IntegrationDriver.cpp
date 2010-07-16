/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
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
#include "OrthogPolyApproximation.hpp"

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
  driverRep(NULL), referenceCount(1)
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


void IntegrationDriver::dimension_preference(const RealVector& dim_pref)
{
  if (driverRep)
    driverRep->dimension_preference(dim_pref); // forward to letter
  else { // default implementation
    RealVector aniso_wts;
    if (!dim_pref.empty()) {
      size_t num_pref = dim_pref.length();
      aniso_wts.sizeUninitialized(num_pref);
      for (size_t i=0; i<num_pref; ++i)
	aniso_wts[i] = (dim_pref[i] == 0.) ? 0. : 1./dim_pref[i];
      // scaling occurs in anisotropic_weights() below
    }
    anisotropic_weights(aniso_wts);
  }
}


void IntegrationDriver::anisotropic_weights(const RealVector& aniso_wts)
{
  if (driverRep)
    driverRep->anisotropic_weights(aniso_wts); // forward to letter
  else
    anisoLevelWts = aniso_wts; // default implementation
}


/** protected function called only from derived class letters. */
void IntegrationDriver::
initialize_rules(const ShortArray& u_types, bool nested_rules,
		 short growth_rate, short nested_uniform_rule,
		 IntArray& int_rules, IntArray& growth_rules)
{
  int_rules.resize(numVars);
  growth_rules.resize(numVars);

  for (size_t i=0; i<numVars; i++) {
    // set int_rules
    switch (u_types[i]) {
    case STD_NORMAL:
      int_rules[i] = (nested_rules) ? GENZ_KEISTER : GAUSS_HERMITE; break;
    case STD_UNIFORM:
      // For tensor-product quadrature, Gauss-Legendre is used due to greater
      // polynomial exactness since nesting is not a concern.  For nested sparse
      // grids, Clenshaw-Curtis or Gauss-Patterson can be better selections.
      // However, sparse grids that are isotropic in level but anisotropic in
      // rule become skewed when mixing Gauss rules with CC.  For this reason,
      // CC is selected only if isotropic in rule (for now).
      int_rules[i] = (nested_rules) ? nested_uniform_rule : GAUSS_LEGENDRE;
      break;
    case STD_EXPONENTIAL: int_rules[i] = GAUSS_LAGUERRE;     break;
    case STD_BETA:        int_rules[i] = GAUSS_JACOBI;       break;
    case STD_GAMMA:       int_rules[i] = GEN_GAUSS_LAGUERRE; break;
    default:              int_rules[i] = GOLUB_WELSCH;       break;
    }

    // set growth_rules
    switch (u_types[i]) {
    case STD_NORMAL: case STD_UNIFORM:
      if (nested_rules) // symmetric exponential growth
	switch (growth_rate) {
	case SLOW_RESTRICTED_GROWTH:
	  growth_rules[i] = SLOW_EXPONENTIAL;     break;
	case MODERATE_RESTRICTED_GROWTH:
	  growth_rules[i] = MODERATE_EXPONENTIAL; break;
	case UNRESTRICTED_GROWTH:
	  growth_rules[i] = FULL_EXPONENTIAL;     break;
	}
      else // symmetric Gaussian linear growth
	growth_rules[i] = (growth_rate == SLOW_RESTRICTED_GROWTH) ?
	  SLOW_LINEAR_ODD : MODERATE_LINEAR;
      break;
    default: // asymmetric Gaussian linear growth
      growth_rules[i] = (growth_rate == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
    }
  }

  ShortArray basis_types, gauss_modes;
  OrthogPolyApproximation::distribution_types(u_types, int_rules, basis_types,
					      gauss_modes);
  OrthogPolyApproximation::distribution_basis(basis_types, gauss_modes,
					      polynomialBasis);
}


/** protected function called only from derived class letters. */
void IntegrationDriver::
initialize_rules(const std::vector<BasisPolynomial>& poly_basis,
		 short growth_rate, IntArray& int_rules, IntArray& growth_rules)
{
  int_rules.resize(numVars);
  growth_rules.resize(numVars);

  for (size_t i=0; i<numVars; i++) {
    int_rules[i] = poly_basis[i].gauss_mode();
    switch (int_rules[i]) {
    case GAUSS_HERMITE: case GAUSS_LEGENDRE: // symmetric Gaussian linear growth
      growth_rules[i] = (growth_rate == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR_ODD : MODERATE_LINEAR; break;
    case GAUSS_PATTERSON: case CLENSHAW_CURTIS: case FEJER2: case GENZ_KEISTER:
      // nested rules with exponential growth
      switch (growth_rate) {
      case SLOW_RESTRICTED_GROWTH:
	growth_rules[i] = SLOW_EXPONENTIAL;     break;
      case MODERATE_RESTRICTED_GROWTH:
	growth_rules[i] = MODERATE_EXPONENTIAL; break;
      case UNRESTRICTED_GROWTH:
	growth_rules[i] = FULL_EXPONENTIAL;     break;
      }
      break;
    default: // asymmetric Gaussian linear growth
      growth_rules[i] = (growth_rate == SLOW_RESTRICTED_GROWTH) ?
	SLOW_LINEAR : MODERATE_LINEAR; break;
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
    OrthogPolyApproximation::distribution_parameters(u_types, dp,
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

} // namespace Pecos
