/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchicalLinearBasis
//- Description:  Implementation code for HierarchicalLinearBasis class
//-               
//- Owner:        Christopher Miller, University of Maryland
//- Contact:      cmiller@math.umd.edu

#include "HierarchicalLinearBasis.hpp"

namespace Pecos{

/// Default constructor throws an exception if left_end_ >= right_end_.
HierarchicalLinearBasis::HierarchicalLinearBasis(const Real& left_end_,
						 const Real& right_end_):
  InterpolationPolynomial(),
  left_end(left_end_),
  right_end(right_end_)
{
  if (left_end_ >= right_end_) {
    PCerr << "Error: The left interpolation endpoint should be less than the right." 
	  << " Got left = " << left_end_ << " right = " << right_end_ 
	  << ". Throwing an exception." << std::endl;
    throw std::invalid_argument("");
  }
}

HierarchicalLinearBasis::~HierarchicalLinearBasis()
{}

/** Point evaluation of a basis element.  get_value indexes into the basis
  via the integer tuple basis_index whose first element specifies the level
  and whos second element specifies the index in that level.*/
const Real& HierarchicalLinearBasis::
get_type1_value(const Real& x, const IntArray& basis_index)
{
  check_valid_index(basis_index);
  const unsigned int level = basis_index[0];
  const unsigned int index = basis_index[1];
  const Real middle = (left_end + right_end) / 2.0;

  //basis at level 1 is 1 on [left,right] and zero off.
  if( level == 1 ) { 
    if( x < left_end || x > right_end ) basisPolyValue = 0.0;
    else basisPolyValue = 1.0;
    return basisPolyValue;
  } else if( level == 2 ){ 
    if( index == 0 ){ // linear connecting (left,1) and (middle,0)
      if( x < left_end || x > middle ) basisPolyValue = 0.0;
      else basisPolyValue = 
	     1 - (x - left_end) / (middle - left_end);
    } else { //linear connecting (middle,0) and (right,1)
      if( x < middle || x > right_end ) basisPolyValue = 0.0;
      else basisPolyValue = 
	     (x - middle) / (right_end - middle);
    }
    return basisPolyValue;
  } else { // General case.
    const Real offset = ( (left_end + right_end) / 2.0 - left_end) 
      / pow(2.0,level-2);
    const Real node_point = get_interp_point(basis_index);
    const Real left_support = node_point - offset;
    const Real right_support = node_point + offset;
    if( x < left_support || x > right_support ) basisPolyValue = 0.0;
    else{
      if( x <= node_point ){
	basisPolyValue = 
	  (x - left_support)/(node_point - left_support);
      } else {
	basisPolyValue = 
	  1 - (node_point - x)/(node_point - right_support);
      }
    }
    return basisPolyValue;
  }
} 

/** Point evaluation of the gradient of a  basis element.  get_value indexes
  into the basis via the integer tuple basis_index whose first element specifies
  the level and whos second element specifies the index in that level.*/
const Real& HierarchicalLinearBasis::
get_type1_gradient(const Real& x, const IntArray& basis_index)
{
  check_valid_index(basis_index);
  const unsigned int level = basis_index[0];
  const unsigned int index = basis_index[1];
  const Real middle = (left_end + right_end) / 2.0;

  //Gradient is defined to be the slope of the linear function where that makes
  //sense and zero off the support and at the transition point.

  if( level == 1 ){
    basisPolyGradient = 0.0;
    return basisPolyGradient;
  } else if( level == 2 ){
    if( index == 0 ){
      if( x <= left_end || x >= middle ) basisPolyGradient = 0.0;
      else basisPolyGradient = 
	     - 1 / (middle - left_end);
    } else {
      if( x <= middle || x >= right_end ) basisPolyGradient = 0.0;
      else basisPolyGradient = 
	     1 / (right_end - middle);
    }
    return basisPolyGradient;
  } else {
    const Real offset = ( (left_end + right_end) / 2.0 - left_end) 
      / pow(2.0,level-2);
    const Real node_point = get_interp_point(basis_index);
    const Real left_support = node_point - offset;
    const Real right_support = node_point + offset;
    if( x <= left_support ||
	x >= right_support ||
	x == node_point) basisPolyGradient = 0.0;
    else{
      if( x < node_point ) basisPolyGradient = 1/offset;
      else basisPolyGradient = -1/offset; //x > node_point since == is above
    }
    return basisPolyGradient;
  }
}

/** Point evaluation of a basis element.  The index i indexes into the basis
  by a lexiographic ordering with lower levels before higher levels and then
  left to right within a given level.*/
const Real& HierarchicalLinearBasis::
get_type1_value(const Real& x, unsigned int i) 
{
  const IntArray basis_index = int_to_intArray(i);
  get_type1_value(x,basis_index);
  return basisPolyValue;
}

/** Point evaluation of the gradient of a basis element.  The index i 
  indexes into the basis by a lexiographic ordering with lower levels 
  before higher levels and then left to right within a given level.*/
const Real& HierarchicalLinearBasis::
get_type1_gradient(const Real& x, unsigned int i) 
{
  const IntArray basis_index = int_to_intArray(i);
  get_type1_gradient(x, basis_index);
  return basisPolyGradient;
}

/// Returns the interpolation point associated with a given basis function.  Address by a level/index pair
const Real HierarchicalLinearBasis::
get_interp_point(const IntArray& basis_index) const
{
  //This seems bad but it's not so rough.  I handle
  //the level 1 and 2 cases manually.  The first
  //point in level 3 is half way between left_end
  //and the middle.  The first point in level 4
  //is a quarter way, level 5 an eighth and so on.
  //This is the offset of the first point from the
  //left end.  Subsequent points on the same level
  //are spaced in increments 2*offset.
  check_valid_index(basis_index);
  const unsigned int level = basis_index[0];
  const unsigned int index = basis_index[1];

  if( level == 1 ) return (left_end + right_end) / 2.0;
  else if( level == 2 ){
    if( index == 0 ) return left_end;
    else return right_end;
  } else {
    Real offset = ( (left_end + right_end) / 2.0 - left_end) 
      / pow(2.0,level-2);
    return left_end + (1 + 2*index)*offset;
  }
}

/// Returns the interpolation point associated with a given basis function. Address by a integer index
const Real HierarchicalLinearBasis::
get_interp_point(unsigned int i) const
{
  const IntArray levelOffset = int_to_intArray( i );
  return this->get_interp_point( levelOffset );
}

/// Returns the left endpoint of interpolation.
const Real& HierarchicalLinearBasis::
get_lower_interp_bound() const
{
  return left_end;
}

/// Returns the right endpoint of interpolation.
const Real& HierarchicalLinearBasis::
get_upper_interp_bound() const
{
  return right_end;
}


/// Converts an integer index for a basis element into a level/index pair.
const IntArray HierarchicalLinearBasis::
int_to_intArray(unsigned int i)
{
  IntArray basis_index(2);
  basis_index[0] = 1;
  basis_index[1] = i;
  if( i == 0 ) return basis_index;
  else if( i == 1 ){
    basis_index[0] = 2;
    basis_index[1] = 0;
    return basis_index;
  }else{
    basis_index[0] = ceil(log2(i)) + 1;  // Works to get the level
    // Subtract all the indices in levels above you.  The
    // remainder is your index in the current level
    for( unsigned int idx = basis_index[0]-1; idx>=1; idx-- ){
      if( idx == 1 ) basis_index[1] -= 1;
      else if( idx == 2 || idx == 3 ) basis_index[1] -= 2;
      else basis_index[1] -= pow(2,idx-2);
    }
    return basis_index;
  }
}

/// Converts a level/index pair into an integer index for a basis element. 
const unsigned int HierarchicalLinearBasis::
intArray_to_int(const IntArray& basis_index)
{
    
  check_valid_index(basis_index);
  const unsigned int level = basis_index[0];
  const unsigned int index = basis_index[1];
  // Start with your index in the current level
  unsigned int integer_index = index;

  // Then add in all of the basis functions from
  // levels above you.
  for( unsigned int i = level-1; i>=1; i-- ){
    if( i == 1 ) integer_index += 1;
    else if( i == 2 || i == 3 ) integer_index += 2;
    else integer_index += pow(2,i-2);
  }
    
  return integer_index;
}

/** Checks that a given IntArray represents a vailid
    level/index pair.  Throws a std::length_error if
    the vector is not size 2 and a std::invalid_argument
    if the pair is not valid.  Also writes error to PCerr.
*/
void HierarchicalLinearBasis::
check_valid_index(const IntArray& basis_index)
{
  // Test for valid index length
  if( basis_index.size() != 2 ){
    PCerr << "IntArray basis_index should have size == 2." 
	  << "  Got basis_index.size() == " << basis_index.size()
	  << std::endl;
    throw( std::length_error("") );
  }
  // Test for valid level index
  else if ( basis_index[0] < 1 ){
    PCerr << "IntArray basis_index[0] must be >=1."
	  << "  Got basis_index[0] == " << basis_index[0]
	  << std::endl;
    throw( std::invalid_argument("") );
  }
  // Test for valid index at a given level.  Note the logic
  // here is odd due to the boundary points only having one
  // child.  C'est la vie.
  else if ( basis_index[0] == 1 && basis_index[1] != 0){
    PCerr << "basis_index[1] must be 0<=basis_index[1]<= 0.  Got"
	  << "  Got basis_index[1] == " << basis_index[1] 
	  << std::endl;
    throw( std::invalid_argument("") );
  }
  else if ( ( basis_index[0] == 2 || basis_index[0] == 3 ) && 
	    (basis_index[1] < 0 || basis_index[1] > 1) ){
    PCerr << "basis_index[1] must be 0<=basis_index[1]<= 1."
          << "  Got basis_index[1] == " << basis_index[1] 
          << std::endl;
    throw( std::invalid_argument("") );
  }
  else if ( basis_index[0] > 3 && 
	    ( basis_index[1] < 0 || 
	      basis_index[1] > pow(2,basis_index[0]-2)-1) ){
    PCerr << "basis_index[1] must be 0<=basis_index[1]<= "
          << pow(2,basis_index[0]-2)-1 
          << "  Got basis_index[1] == " << basis_index[1] 
          << std::endl;
    throw( std::invalid_argument("") );
  }
}

} // End namespace Pecos
