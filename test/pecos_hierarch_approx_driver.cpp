/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_hierarch_approx_driver.cpp
    \brief A test program for HierarchInterpPolyApproximation class. */

#include "HierarchInterpPolyApproximation.hpp"
#include "LocalRefinableDriver.hpp"

using namespace Pecos;

Real test_function(const RealVector& x);
RealVector test_function_grad(const RealVector& x);


int main(int argc, char** argv)
{

  HierarchInterpPolyApproximation *a = 
    new HierarchInterpPolyApproximation(0,1,false);
  LocalRefinableDriver l_driver;
  l_driver.initialize_grid(RealArray(1,0),RealArray(1,1),1);
  const std::vector<CollocationPoint>& col_pts = 
    l_driver.get_collocation_points();
  a->integration_driver_rep(&l_driver);

  RealVector x(1);
  Real fn_val;
  RealVector fn_grad;
  RealSymMatrix fn_hess;
  fn_val = 2.718;
  x[0] = col_pts[0].get_point()[0];

  SurrogateData points;
  points.push_back( SurrogateDataVars(x),
		    SurrogateDataResp(fn_val,fn_grad,fn_hess,1) );

  a->surrogate_data(points);

  a->compute_coefficients();

  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient(x)[0] << std::endl;

  assert(a->value(x) == 2.718);
  assert(a->gradient(x)[0] == 0);

  l_driver.refine_globally();

  x[0] = col_pts[1].get_point()[0];
  fn_val = 0.0;
  points.push_back(SurrogateDataVars(x),
		   SurrogateDataResp(fn_val,fn_grad,fn_hess,1) );

  x[0] = col_pts[2].get_point()[0];
  fn_val = 5.6;
  points.push_back(SurrogateDataVars(x),
		   SurrogateDataResp(fn_val,fn_grad,fn_hess,1) );

  a->surrogate_data(points);
  a->compute_coefficients();
  
  x[0] = col_pts[1].get_point()[0];

  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient(x)[0] << std::endl;

  assert(a->value(x) == 0);
  assert(a->gradient(x)[0] == 0);

  x[0] = col_pts[2].get_point()[0];

  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient(x)[0] << std::endl;

  assert(a->value(x) == 5.6);
  assert(a->gradient(x)[0] == 0);

  x[0] = .25;

  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient(x)[0] << std::endl;

  assert(std::abs(a->value(x) - 1.359) < 1e-10);
  assert(std::abs(a->gradient(x)[0] - 5.436) < 1e-10);

  x[0] = .75;

  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient(x)[0] << std::endl;

  assert(std::abs(a->value(x) - 4.159) < 1e-10);
  assert(std::abs(a->gradient(x)[0] - 5.764) < 1e-10);

  BoolDeque refinementSelector(2,false);
  refinementSelector[0] = true;
  l_driver.refine_locally(refinementSelector);
  x[0] = (col_pts[3].get_point())[0];
  fn_val = -1;
  
  points.push_back( SurrogateDataVars(x),
		    SurrogateDataResp(fn_val,fn_grad,fn_hess,1) );

  a->surrogate_data(points);
  a->increment_coefficients();
  
  x[0] = .25;
  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient(x)[0] << std::endl;
  
  //2D example.
  //delete a;
  //a = new HierarchInterpPolyApproximation(0,1,false);
  points.clear_data();
  l_driver.initialize_grid(RealArray(2,-1),RealArray(2,1),11);

  RealMatrix col_pts_mat;
  l_driver.compute_grid(col_pts_mat);
  
  for (unsigned int idx = 0; idx < l_driver.grid_size() ; ++idx) {
    RealVector x(Teuchos::Copy,col_pts_mat[idx],2);
    Real fn_val = test_function(x);
    points.push_back( SurrogateDataVars(x),
		      SurrogateDataResp(fn_val,fn_grad,fn_hess,1) );
  }
  a->surrogate_data(points);
  std::cout << "Grid size is: " << l_driver.grid_size() << std::endl;
  a->compute_coefficients();

  for (unsigned int idx = 0; idx < l_driver.grid_size(); ++idx ) {
    assert( std::abs( a->value( points.continuous_variables(idx) ) -  
		      test_function(points.continuous_variables(idx)) < 1e-9) );
  }
  x.size(2);
  x[0] = .69324;
  x[1] = .84529;
  std::cout << "Value =  " << test_function(x) << std::endl;
  std::cout << "Approximate = " << a->value(x) << std::endl;
  std::cout << "Error = " << std::abs(a->value(x) - test_function(x)) << std::endl;

  std::cout << "Mean = " << a->mean() << std::endl;

  //2D example with gradients.  Not fully implemented yet.
  /*
  delete a;
  a = new HierarchInterpPolyApproximation(0,1,true);
  a->integration_driver_rep(&l_driver);
  points.clear_data();
  l_driver.initialize_grid(RealArray(2,-1),RealArray(2,1),9);

  l_driver.compute_grid(col_pts_mat);
  
  fn_grad.size(2);
  for (unsigned int idx = 0; idx < l_driver.grid_size() ; ++idx) {
    RealVector x(Teuchos::Copy,col_pts_mat[idx],2);
    Real fn_val = test_function(x);
    fn_grad = test_function_grad(x);
    points.push_back( SurrogateDataVars(x),
		      SurrogateDataResp(fn_val,fn_grad,fn_hess) );
  }
  
  a->surrogate_data(points);
  std::cout << "Grid size is: " << l_driver.grid_size() << std::endl;
  a->compute_coefficients();

  for (unsigned int idx = 0; idx < l_driver.grid_size(); ++idx ) {
    std::cout << "Value = " << test_function(points.continuous_variables(idx)) 
	      << " approx = " << a->value( points.continuous_variables(idx) ) << std::endl;
    assert( std::abs( a->value( points.continuous_variables(idx) ) -  
		      test_function(points.continuous_variables(idx)) < 1e-9) );
  }
  x.size(2);
  x[0] = .69324;
  x[1] = .84529;
  std::cout << "Value =  " << test_function(x) << std::endl;
  std::cout << "Approximate = " << a->value(x) << std::endl;
  std::cout << "Error = " << std::abs(a->value(x) - test_function(x)) << std::endl;

  std::cout << "Mean = " << a->mean() << std::endl;
  */

  return EXIT_SUCCESS;

}

Real test_function(const RealVector& x) {

  unsigned int dimension = x.length();
  Real return_val = 1;
  for ( unsigned int idx = 0; idx < dimension ; ++idx ) {
    if ( x[idx] == 0 ) return 0;
    else return_val *= std::abs(x[idx])*std::sin(1/x[idx]);
  }
  return return_val;
}

RealVector test_function_grad(const RealVector& x){
  
  unsigned int dimension = x.length();
  RealVector grad(dimension,1);
  for ( unsigned int idx = 0; idx < dimension; ++idx ){
    for ( unsigned int idx2 = 0; idx2 < dimension; ++idx2) {
      if ( idx != idx2 ) {
	if ( x[idx2] != 0 ) {
	  grad[idx] *= std::abs(x[idx2])*std::sin(1/x[idx2]);
	} else  grad[idx] *= 0;  
      }
    }
  }
  for ( unsigned int idx = 0; idx< dimension; ++idx ){
    if ( x[idx] < 0 ) {
      grad[idx] *= std::abs(x[idx])*cos(1/x[idx])*(-x[idx]*x[idx]) - sin(1/x[idx]);
    } else if (x[idx] > 0) {
      grad[idx] *= std::abs(x[idx])*cos(1/x[idx])*(-x[idx]*x[idx]) + sin(1/x[idx]);
    } else {
      grad[idx] *= 0;
    }
  }
  return grad;
}
