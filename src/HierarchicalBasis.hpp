/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchicalBasis
//- Description:  Class for 1-D HierarchicalLinearBasis
//-               
//- Owner:        Christopher Miller: University of Maryland at College Park
//- Contact:      cmiller@math.umd.edu

#ifndef HIERARCHICABASIS_HPP
#define HIERARCHICABASIS_HPP

#include "InterpolationPolynomial.hpp"
#include "pecos_data_types.hpp"
#include "RefinablePointSet.hpp"

namespace Pecos {

/// Derived basis for 1-D piecewise linear interpolants expressed in a hierarchical basis. 

/** The HierarchicalBasis class evaluates a univariate piecewise
    polynomial interpolation function defined on a hierarchical grid.  The grid is
    provided by the RefinablePointSet passed to the constructor. 

    For additional information see:
    X. Ma and N. Zabras, "An adaptive hierarchical sparse grid 
    collocation algorithm for the solution of stochastic 
    differential equations", J. Comput. Phys. 228:3084--3113,
    2009.
*/

class HierarchicalBasis: public InterpolationPolynomial
{
public:

  //
  //- Heading: Constructor and Destructor
  //

  /// Standard constructor
  HierarchicalBasis(RefinablePointSet& pointSet_, 
		    short interpType_=PIECEWISE_LINEAR_INTERP);
    
  /// Destructor
  virtual ~HierarchicalBasis();

  //
  //- Heading: Virtual function redefinitions
  //

  const Real& get_type1_value(const Real& x, const unsigned int i);
  const Real& get_type2_value(const Real& x, const unsigned int i);

  const Real& get_type1_gradient(const Real& x, const unsigned int i);
  const Real& get_type2_gradient(const Real& x, const unsigned int i);
       
protected:
    
  //
  //- Heading: Data
  //

  /// Grid of definition
  RefinablePointSet& pointSet;

private:

  /// type of polynomial interpolant: PIECEWISE_LINEAR_INTERP,
  /// PIECEWISE_QUADRATIC_INTERP, or PIECEWISE_CUBIC_INTERP
  short interpType;

}; //End class definition

} //End namespace Pecos

#endif
