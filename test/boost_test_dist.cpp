/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "boost_test_dist.hpp"
#ifdef HAVE_GSL
#include "gsl/gsl_sf_gamma.h"
#endif
#ifdef HAVE_BOOST
#include <boost/math/distributions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/exponential.hpp>
#endif
#include <algorithm>


boost_test_dist::boost_test_dist()
{
}

boost_test_dist::~boost_test_dist()
{ 
}

  void boost_test_dist::print_comparison()
{
  double cdf_value,quantile_value,this_gamma;
  int i,j;
  
  PCout << "Quantiles for normal" << '\n';
  for (i=1; i<10; i++){
    cdf_value = 0.1*i;
#ifdef HAVE_BOOST
    boost::math::normal norm(0,1);
    PCout << "Boost " << quantile(norm,cdf_value);
#endif
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_cdf_ugaussian_Pinv(cdf_value) << '\n';
#endif
  }

  PCout << "CDF values for normal" << '\n';
  for (i=-4; i<5; i++){
    quantile_value = 1.0*i;
#ifdef HAVE_BOOST
    boost::math::normal norm(0,1);
    PCout << "Boost " << cdf(norm,quantile_value);
#endif
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_cdf_ugaussian_P(quantile_value) << '\n';
#endif // HAVE_GSL
  }
  
  PCout << "PDF values for normal" << '\n';
  for (i=-4; i<5; i++){
    quantile_value = 1.0*i;
#ifdef HAVE_BOOST
    boost::math::normal norm(0,1);
    PCout << "Boost " << pdf(norm,quantile_value);
#endif
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_ran_ugaussian_pdf(quantile_value) << '\n';
#endif // HAVE_GSL
  }

  PCout << "Gamma function values" << '\n';
  for (i=2; i<6; i++){
    quantile_value = 1.0*i+0.5;
#ifdef HAVE_BOOST
    this_gamma = boost::math::tgamma(quantile_value);
    PCout << "Boost " << this_gamma ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    this_gamma = gsl_sf_gamma(quantile_value);
    PCout << "  GSL " << this_gamma << '\n';
#endif // HAVE_GSL
  }

  PCout << "Quantiles for gamma distribution" << '\n';
  for (i=1; i<10; i++){
    cdf_value = 0.1*i;
#ifdef HAVE_BOOST
    boost::math::gamma_distribution<> gamma1(3.,1.);
    PCout << "Boost " << quantile(gamma1,cdf_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_cdf_gamma_Pinv(cdf_value, 3.,1.) << '\n';
#endif // HAVE_GSL
  }

  PCout << "CDF values for gamma distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 1.0*i; 
#ifdef HAVE_BOOST
    boost::math::gamma_distribution<> gamma1(3.,1.);
    PCout << "Boost " << cdf(gamma1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_cdf_gamma_P(quantile_value, 3.,1.) << '\n';
#endif // HAVE_GSL
  }

  PCout << "PDF values for gamma distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 1.0*i; 
#ifdef HAVE_BOOST
    boost::math::gamma_distribution<> gamma1(3.,1.);
    PCout << "Boost " << pdf(gamma1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_ran_gamma_pdf(quantile_value, 3.,1.) << '\n';
#endif // HAVE_GSL
  }
  
  PCout << "PDF values for weibull distribution" << '\n';
  for (i=1; i<5; i++){
    quantile_value = 0.1*i; 
#ifdef HAVE_BOOST
    boost::math::weibull_distribution<> weibull1(3.,1.);
    PCout << "Boost " << pdf(weibull1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_ran_weibull_pdf(quantile_value, 1.,3.) << '\n';
#endif // HAVE_GSL
  }

  PCout << "CDF values for weibull distribution" << '\n';
  for (i=1; i<5; i++){
    quantile_value = 0.1*i; 
#ifdef HAVE_BOOST
    boost::math::weibull_distribution<> weibull1(3.,1.);
    PCout << "Boost " << cdf(weibull1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_cdf_weibull_P(quantile_value, 1.,3.) << '\n';
#endif // HAVE_GSL
  }

  PCout << "CDF values for beta distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 0.1*i; 
#ifdef HAVE_BOOST
    boost::math::beta_distribution<> beta1(3.,1.);
    PCout << "Boost " << cdf(beta1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_cdf_beta_P(quantile_value, 3.,1.) << '\n';
#endif // HAVE_GSL
  }

  PCout << "PDF values for beta distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 0.1*i; 
#ifdef HAVE_BOOST
    boost::math::beta_distribution<> beta1(3.,1.);
    PCout << "Boost " << pdf(beta1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_ran_beta_pdf(quantile_value, 3., 1.) << '\n';
#endif // HAVE_GSL
  }

  PCout << "PDF values for F distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 0.1*i; 
#ifdef HAVE_BOOST
    boost::math::fisher_f_distribution<> fisher(10.,4.);
    PCout << "Boost " << pdf(fisher,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_ran_fdist_pdf(quantile_value, 10., 4.) << '\n';
#endif // HAVE_GSL
  }

  PCout << "Quantiles for t-distribution" << '\n';
  for (i=1; i<10; i++){
    cdf_value = 0.1*i;
#ifdef HAVE_BOOST
    boost::math::students_t_distribution<> tdist(10);
    PCout << "Boost " << quantile(tdist,cdf_value);
#endif
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_cdf_tdist_Pinv(cdf_value,10) << '\n';
#endif
  }

  PCout << "Quantiles for Chi-squared distribution" << '\n';
  for (i=1; i<10; i++){
    cdf_value = 0.1*i;
#ifdef HAVE_BOOST
    boost::math::chi_squared_distribution<> chisq(10);
    PCout << "Boost " << quantile(chisq,cdf_value);
#endif
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_cdf_chisq_Pinv(cdf_value,10) << '\n';
#endif
  }

 PCout << "CDF values for exponential" << '\n';
  for (i=1; i<5; i++){
    quantile_value = 0.1*i;
#ifdef HAVE_BOOST
    boost::math::exponential_distribution<> expon1(1.5);
    PCout << "Boost " << cdf(expon1,quantile_value);
#endif
#ifdef HAVE_GSL
    PCout << "  GSL " << gsl_cdf_exponential_P(quantile_value,(1/1.5)) << '\n';
#endif // HAVE_GSL
  }
}

int main(int argc, char* argv[])
{

  boost_test_dist bt;
  bt.print_comparison();
  return 0;
}
