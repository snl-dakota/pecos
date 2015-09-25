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

void write_USAS(std::ostream& s, const Pecos::UShortArraySet &a);
void write_US2A(std::ostream& s, const Pecos::UShort2DArray  &a);

/// A driver program for PECOS.
int main(int argc, char* argv[])
{
  using namespace Pecos;

  std::cout << "Instantiating CombinedSparseGridDriver:\n";
  unsigned short level = STARTLEV;  // reference grid level
  RealVector dimension_pref; // empty -> isotropic
  short refine_cntl = DIMENSION_ADAPTIVE_CONTROL_GENERALIZED;
  CombinedSparseGridDriver csg_driver(level, dimension_pref, refine_cntl);

  std::cout << "Instantiating basis:\n";
  size_t num_vars = NUMVARS;
  std::vector<BasisPolynomial> poly_basis(num_vars);
  for (int i=0; i<num_vars; ++i)
    poly_basis[i] = BasisPolynomial(LEGENDRE_ORTHOG);
  csg_driver.initialize_grid(poly_basis);

  // initial grid
  RealMatrix variable_sets;
  csg_driver.compute_grid(variable_sets);

#define DEBUG
#ifdef DEBUG
  write_data(std::cout, variable_sets, true, true, true);
  write_US2A(std::cout, csg_driver.smolyak_multi_index());
#endif

  // start refinement
  csg_driver.initialize_sets();

  UShortArraySet a;
  for ( int iter = 0; iter<NITER; iter++) {
    a = csg_driver.active_multi_index();
    std::cout<<"Refine, iteration: "<<iter+1<<'\n';
    write_USAS(std::cout, a) ;
  
    std::vector<short unsigned int> asave;
    RealMatrix vsets1;
    int choose = 0;
    for (UShortArraySet::iterator it=a.begin(); it!=a.end(); ++it) {
      int pick = std::rand();
      if ( pick > choose) {
        asave = *it;
        choose = pick;
      }
      csg_driver.push_trial_set(*it);
      csg_driver.compute_trial_grid(vsets1); 
      //std::cout << "Sparse grid points:\n";
      //write_data(std::cout, vsets1, false, true, true);
      csg_driver.pop_trial_set();
    }
    csg_driver.update_sets(asave);
    csg_driver.update_reference();
  }
  csg_driver.finalize_sets(true, false); // use embedded output option

  // Print final sets
  //std::cout<<"Final set:\n";
  //write_US2A(std::cout, csg_driver.smolyak_multi_index());
  
  return (0);
}

void write_US2A(std::ostream& s, const Pecos::UShort2DArray &a)
{
  s << "-----------------------------------------\n";
  size_t i, j, num_a = a.size();
  for (i=0; i<num_a; ++i) {
    const Pecos::UShortArray& aa = a[i];
    for (j=0; j < aa.size(); ++j)
      s<<std::setw(5)<<aa[j];
    s<<"\n";
  }
  s << "-----------------------------------------\n";
  return ;
}

void write_USAS(std::ostream& s, const Pecos::UShortArraySet &a)
{
  s << "-----------------------------------------\n";
  Pecos::UShortArraySet::const_iterator cit;
  for (cit=a.begin(); cit!=a.end(); ++cit) {
    const Pecos::UShortArray& aa = *cit;
    for (size_t j=0; j < aa.size(); ++j)
      s<<std::setw(5)<<aa[j];
    s<<'\n';
  }
  s << "-----------------------------------------\n";
  return ;
}
