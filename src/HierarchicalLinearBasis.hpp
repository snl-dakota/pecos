/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchicalLinearBasis
//- Description:  Class for 1-D HierarchicalLinearBasis
//-               
//- Owner:        Christopher Miller: University of Maryland at College Park
//- Contact:      cmiller@math.umd.edu

#ifndef HIERARCHICALLINEARBASIS_HPP
#define HIERARCHICALLINEARBASIS_HPP

#include "InterpolationPolynomial.hpp"
#include "pecos_data_types.hpp"

namespace Pecos {

/// Derived basis for 1-D piecewise linear interpolants expressed in a hierarchical basis. 

/** The HierarchicalLinearBasis class evaluates a univariate piecewise
    linear interpolation functions defined on a hierarchical grid of
    equidistant spaced points on an interval [a,b].  The first level
    of the grid consists of the midpoint between a and b.  The second
    level contains the endpoints. The third level contains the midpoints
    between the endpoints and the average of a and b and so on.  This
    induces a binary tree structure on the set of interpolation nodes.

    Individual elements within the basis can be accessed by a size
    2 IntArray [level,index], where index runs left to right over a
    given level.  Or basis elements can be accessed by a single unsigned
    int where the ordering is lexiographic from low level to high and
    then from left to right inside a given level.

    For additional information see:
    X. Ma and N. Zabras, "An adaptive hierarchical sparse grid 
    collocation algorithm for the solution of stochastic 
    differential equations", J. Comput. Phys. 228:3084--3113,
    2009.
*/

class HierarchicalLinearBasis: public InterpolationPolynomial
{
public:

  //
  //- Heading: Constructor and Destructor
  //

  /// Standard constructor
  HierarchicalLinearBasis(const Real& left_end_ = 0.0,
			  const Real& right_end_ = 1.0);
    
  /// Destructor
  virtual ~HierarchicalLinearBasis();

  //
  //- Heading: Virtual function redefinitions
  //

  virtual const Real& get_type1_value(const Real& x, unsigned int i); 
  virtual const Real& get_type1_gradient(const Real& x, unsigned int i);

  //
  //- Heading: Function definitions
  //

  virtual const Real& get_type1_value(const Real& x, 
				      const IntArray& basis_index);
  virtual const Real& get_type1_gradient(const Real& x, 
					 const IntArray& basis_index);
    
  virtual const Real get_interp_point(const IntArray& basis_index) const;
  virtual const Real get_interp_point(unsigned int i) const;
    
  virtual const Real& get_lower_interp_bound() const;
  virtual const Real& get_upper_interp_bound() const;

  //
  //- Heading: Static helper function definitions
  //

    
  static const IntArray int_to_intArray(unsigned int i);
  static const unsigned int intArray_to_int(const IntArray& basis_index);

  static void check_valid_index(const IntArray& basis_index);

protected:
    
  //
  //- Heading: Data
  //

  /// Grid Lower Bound
  Real left_end;
  /// Grid Upper Bound
  Real right_end;

private:

}; //End class definition

} //End namespace Pecos

#endif
