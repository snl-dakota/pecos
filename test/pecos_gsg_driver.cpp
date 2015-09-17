/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_int_driver.cpp
    \brief A driver program for PECOS */

#include "CombinedSparseGridDriver.hpp"
#include "TensorProductDriver.hpp"
#include "CubatureDriver.hpp"
#include "pecos_data_types.hpp"
//#include "LocalRefinableDriver.hpp"

#define NUMVARS  3
#define STARTLEV 1
#define NITER    10

/// A driver program for PECOS.
void write_USAS(std::ostream& s, const Pecos::UShortArraySet &a);

int main(int argc, char* argv[])
{

  std::cout << "Instantiating CombinedSparseGridDriver:\n";
  unsigned short level = STARTLEV;  // reference grid level
  Pecos::RealVector dimension_pref; // empty -> isotropic
  Pecos::CombinedSparseGridDriver csg_driver(level, dimension_pref);

  std::cout << "Instantiating basis:\n";
  size_t num_vars = NUMVARS;
  std::vector<Pecos::BasisPolynomial> poly_basis(num_vars);
  for (int i=0; i<num_vars; ++i)
    poly_basis[i] = Pecos::BasisPolynomial(Pecos::LEGENDRE_ORTHOG);
  csg_driver.initialize_grid(poly_basis);

  // initial grid
  Pecos::RealMatrix variable_sets;
  csg_driver.refinement_control(Pecos::DIMENSION_ADAPTIVE_CONTROL_GENERALIZED);
  csg_driver.compute_grid(variable_sets);

#define DEBUG
#ifdef DEBUG
  Pecos::write_data(std::cout, variable_sets, true, true, true);
#endif

  // start refinement
  csg_driver.initialize_sets();

#ifdef DEBUG
  Pecos::UShortArraySet aold = csg_driver.old_multi_index();
  write_USAS(std::cout, aold) ;
#endif

  Pecos::UShortArraySet a;
  for ( int iter = 0; iter<NITER; iter++) {
    std::cout<<"Refine, iteration:"<<iter+1<<std::endl;
    a = csg_driver.active_multi_index();
    write_USAS(std::cout, a) ;
    //std::cout<<a<<std::endl;
    //std::cout<<a.size()<<std::endl;
  
    std::vector<short unsigned int> asave;
    Pecos::RealMatrix vsets1;
    int choose = 0;
    for (Pecos::UShortArraySet::iterator it=a.begin(); it!=a.end(); ++it) {
      int pick = std::rand();
      if ( pick > choose) {
        asave = *it;
        choose = pick;
      }
      csg_driver.push_trial_set(*it);
      csg_driver.compute_trial_grid(vsets1); 
      //std::cout << "Sparse grid points:\n";
      //Pecos::write_data(std::cout, vsets1, false, true, true);
      csg_driver.pop_trial_set();
    }
    csg_driver.update_sets(asave);
    csg_driver.update_reference();
  }

  // Print final set
  aold = csg_driver.old_multi_index();
  std::cout<<"Final set:"<<std::endl;
  write_USAS(std::cout, aold) ;
  
  return (0);

}

void write_USAS(std::ostream& s, const Pecos::UShortArraySet &a) {

  s << "-----------------------------------------\n";
  Pecos::UShortArraySet::iterator it;
  for (it=a.begin(); it!=a.end(); it++) {
    Pecos::UShortArray aa = *it;
    for (int i=0; i < aa.size(); i++ )
      s<<aa[i]<<" ";
    s<<"\n";
  }
  s << "-----------------------------------------\n";
  return ;
}
