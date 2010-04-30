/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_int_driver.cpp
    \brief A driver program for PECOS */

#include "SparseGridDriver.hpp"
#include "TensorProductDriver.hpp"
#include "CubatureDriver.hpp"


/// A driver program for PECOS.

/** Generates sparse, tensor, and cubature grids for numerical integration. */

int main(int argc, char* argv[])
{
  std::cout << "Instantiating basis:\n";
  size_t i, num_vars = 4;
  std::vector<Pecos::BasisPolynomial> poly_basis(num_vars);
  for (i=0; i<num_vars; ++i)
    poly_basis[i] = Pecos::BasisPolynomial(Pecos::HERMITE);

  // Smolyak sparse grids
  std::cout << "Instantiating SparseGridDriver:\n";
  Pecos::SparseGridDriver sg_driver;
  unsigned short level = 3; Pecos::RealVector dimension_pref; // isotropic
  sg_driver.initialize_grid(poly_basis, level, dimension_pref);
                          //store_1d_gauss, growth_rate);
  sg_driver.compute_grid();
  std::cout << "Sparse grid points:\n";
  Pecos::write_data(std::cout, sg_driver.variable_sets(), false, true, true);
  std::cout << "Sparse grid weights:\n";
  Pecos::write_data(std::cout, sg_driver.weight_sets());

  // Tensor-product quadrature
  std::cout << "Instantiating TensorProductDriver:\n";
  Pecos::TensorProductDriver tp_driver;
  Pecos::UShortArray quad_order(num_vars, 3);
  tp_driver.initialize_grid(poly_basis, quad_order);//, growth_rate);
  tp_driver.compute_grid();
  std::cout << "Tensor grid points:\n";
  Pecos::write_data(std::cout, tp_driver.variable_sets(), false, true, true);
  std::cout << "Tensor grid weights:\n";
  Pecos::write_data(std::cout, tp_driver.weight_sets());

  // Cubature
  std::cout << "Instantiating CubatureDriver:\n";
  Pecos::CubatureDriver c_driver;
  unsigned short integrand_order = 5;
  c_driver.initialize_grid(poly_basis, integrand_order);
  c_driver.compute_grid();
  std::cout << "Cubature points:\n";
  Pecos::write_data(std::cout, c_driver.variable_sets(), false, true, true);
  std::cout << "Cubature weights:\n";
  Pecos::write_data(std::cout, c_driver.weight_sets());

  return 0;
}
