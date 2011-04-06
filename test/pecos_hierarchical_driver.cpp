/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_hierarchical_driver.cpp
    \brief A test program for HierarchicalLinearBasis class. */


#include "HierarchicalBasis.hpp"
#include "UniformRefinablePointSet.hpp"
using namespace Pecos;
int main(int argc, char** argv)
{
  
  RefinablePointSet *pointSet = new UniformRefinablePointSet(-1,1);
  HierarchicalBasis *a = new HierarchicalBasis(*pointSet);

  if ( ( a->get_type1_value(0.1, 0) != 1 ) ||
       ( a->get_type1_value(-1.1,0) != 0 ) ){
    std::cout << "test failure:  basis[0] did not evaluate correctly."
              << " Got f(0.1) = " << a->get_type1_value(0.1, 0)
              << " and f(-1.1) = " << a->get_type1_value(-1.1, 0) 
              << std::endl;
    return EXIT_FAILURE;
  }

  pointSet->refine_all();
  if ( ( a->get_type1_value(-1, 1) != 1.0 ) ||
       ( a->get_type1_value(-.5,1) != 0.5 ) ||
       ( a->get_type1_value(-1.1,1) != 0.0 ) ||
       ( a->get_type1_value(0.1,1) != 0.0 ) ){
    std::cout << "test failure:  basis[1] did not evaluate correctly."
              << " Got f(-1) = " << a->get_type1_value(-1, 1)
              << " and f(-.5) = " << a->get_type1_value(-.5, 1) 
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->get_type1_value(1, 2) != 1.0 ) ||
       ( a->get_type1_value(.5,2) != 0.5 ) ||
       ( a->get_type1_value(-.1,2) != 0.0 ) ||
       ( a->get_type1_value(1.1,2) != 0.0 ) ){
    std::cout << "test failure:  basis[2] did not evaluate correctly."
              << " Got f(1) = " << a->get_type1_value(1, 2)
              << " and f(.5) = " << a->get_type1_value(.5, 2) 
              << std::endl;
    return EXIT_FAILURE;
  }

  pointSet->refine_all();
  BoolDeque bitset(2,false);
  bitset[1] = true;
  pointSet->refine(bitset);

  if ( ( a->get_type1_value(.25, 5) != 1.0 ) ||
       ( a->get_type1_value(.5,5) != 0.0 ) ||
       ( a->get_type1_value(.75,5) != 0.0 ) ||
       ( a->get_type1_value(.125,5) != 0.5 ) ||
       ( a->get_type1_value(.25 + .125,5) != .5 ) ||
       ( a->get_type1_value(-.1,5) != 0 ) ||
       ( a->get_type1_value(.51,5) != 0 ) ){
    std::cout << "test failure:  basis[5] did not evaluate correctly."
              << std::endl;
    return EXIT_FAILURE;
  } 

  return EXIT_SUCCESS;
  
}
