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

static const char rcsId[]="@(#) $Id: IntegrationDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


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
