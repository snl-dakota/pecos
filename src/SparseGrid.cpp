
//- Class:	 SparseGrid
//- Description: Implementation code for SparseGrid class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "SparseGrid.hpp"
#include "sandia_rules.H"
#include "sparse_grid_mixed_growth.H"
#include "sgmga.H"
#include "NumericGenOrthogPolynomial.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: SparseGrid.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


// TO DO: avoid recalculating existing Gauss pts within following functions
void SparseGrid::
bounded_normal_gauss_points(int order, int num_params, double* params,
			    double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "bounded_normal_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_normal_distribution(params[0], params[1], params[2], params[3]);
                                //(mean,      stdev,     lwr,       upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGrid::
bounded_normal_gauss_weights(int order, int num_params, double* params,
			     double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "bounded_normal_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_normal_distribution(params[0], params[1], params[2], params[3]);
                                //(mean,      stdev,     lwr,       upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGrid::
lognormal_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "lognormal_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.lognormal_distribution(params[0], params[1]); //(mean, stdev);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGrid::
lognormal_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "lognormal_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.lognormal_distribution(params[0], params[1]); //(mean, stdev);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGrid::
bounded_lognormal_gauss_points(int order, int num_params, double* params,
			       double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "bounded_lognormal_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_lognormal_distribution(params[0], params[1], params[2],
				      params[3]); //(mean, stdev, lwr, upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGrid::
bounded_lognormal_gauss_weights(int order, int num_params, double* params,
				double* data)
{
  if (num_params != 4) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "bounded_lognormal_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.bounded_lognormal_distribution(params[0], params[1], params[2],
				      params[3]); //(mean, stdev, lwr, upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGrid::
loguniform_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "loguniform_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.loguniform_distribution(params[0], params[1]); //(lwr, upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGrid::
loguniform_gauss_weights(int order, int num_params, double* params,
			 double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "loguniform_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.loguniform_distribution(params[0], params[1]); //(lwr, upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGrid::
triangular_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 3) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "triangular_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.triangular_distribution(params[0], params[1], params[2]);
                            //(mode,      lwr,       upr);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGrid::
triangular_gauss_weights(int order, int num_params, double* params,
			 double* data)
{
  if (num_params != 3) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "triangular_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.triangular_distribution(params[0], params[1], params[2]);
                            //(mode,      lwr,       upr);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGrid::
gumbel_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "gumbel_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.gumbel_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGrid::
gumbel_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "gumbel_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.gumbel_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGrid::
frechet_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "frechet_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.frechet_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGrid::
frechet_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "frechet_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.frechet_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGrid::
weibull_gauss_points(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "weibull_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.weibull_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGrid::
weibull_gauss_weights(int order, int num_params, double* params, double* data)
{
  if (num_params != 2) {
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "weibull_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  ngop.weibull_distribution(params[0], params[1]); //(alpha, beta);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}


void SparseGrid::
histogram_bin_gauss_points(int order, int num_params, double* params,
			   double* data)
{
  if (num_params < 4 || num_params % 2) { // need at least 2 (x,c) pairs
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "histogram_bin_gauss_points()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  RealVector param_vec;
  copy_data(params, num_params, param_vec);
  ngop.histogram_bin_distribution(param_vec);
  const RealArray& gauss_pts = ngop.gauss_points(order);
  std::copy(gauss_pts.begin(), gauss_pts.begin()+order, data);
}


void SparseGrid::
histogram_bin_gauss_weights(int order, int num_params, double* params,
			    double* data)
{
  if (num_params < 4 || num_params % 2) { // need at least 2 (x,c) pairs
    PCerr << "Error: wrong number of distribution parameters in SparseGrid::"
	  << "histogram_bin_gauss_weights()" << std::endl;
    abort_handler(-1);
  }
  NumericGenOrthogPolynomial ngop;
  RealVector param_vec;
  copy_data(params, num_params, param_vec);
  ngop.histogram_bin_distribution(param_vec);
  const RealArray& gauss_wts = ngop.gauss_weights(order);
  std::copy(gauss_wts.begin(), gauss_wts.begin()+order, data);
}

} // namespace Pecos
