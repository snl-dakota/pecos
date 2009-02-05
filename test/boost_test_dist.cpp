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
#include <iomanip>


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
  PCout.precision(16); 
  PCout.setf(std::ios::fixed);
  PCout.setf(std::ios::showpoint);
  
  PCout << "Quantiles for normal" << '\n';
  for (i=0; i<11; i++){
    cdf_value = 0.1*i;
#ifdef HAVE_BOOST
    Pecos:: normal_dist norm(0,1);
    // Pecos::Real q = bmth::quantile(norm, cdf_value);
    PCout << std::setw(16) << "Boost " << bmth::quantile(norm,cdf_value);
#endif
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_cdf_ugaussian_Pinv(cdf_value) << '\n';
#endif
  }

  PCout << std::setw(16) << "CDF values for normal" << '\n';
  for (i=-4; i<5; i++){
    quantile_value = 1.0*i;
#ifdef HAVE_BOOST
    Pecos::normal_dist norm(0,1);
    PCout << std::setw(16) << "Boost " << bmth::cdf(norm,quantile_value);
#endif
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_cdf_ugaussian_P(quantile_value) << '\n';
#endif // HAVE_GSL
  }
  
  PCout << std::setw(16) << "PDF values for normal" << '\n';
  for (i=-4; i<5; i++){
    quantile_value = 1.0*i;
#ifdef HAVE_BOOST
    Pecos::normal_dist norm(0,1);
    PCout << std::setw(16) << "Boost " << bmth::pdf(norm,quantile_value);
#endif
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_ran_ugaussian_pdf(quantile_value) << '\n';
#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "Gamma function values" << '\n';
  for (i=2; i<6; i++){
    quantile_value = 1.0*i+0.5;
#ifdef HAVE_BOOST
    this_gamma = bmth::tgamma(quantile_value);
    PCout << std::setw(16) << "Boost " << this_gamma ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    this_gamma = gsl_sf_gamma(quantile_value);
    PCout << std::setw(16) << "  GSL " << this_gamma << '\n';
#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "Beta function values" << '\n';
  for (i=1; i<5; i++){
    double alpha1 = 1.0*i;
#ifdef HAVE_BOOST
    this_gamma = bmth::beta(alpha1, 1.0);
    PCout << std::setw(16) << "Boost " << this_gamma ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    this_gamma = gsl_sf_beta(alpha1, 1.0);
    PCout << std::setw(16) << "  GSL " << this_gamma << '\n';
#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "Quantiles for gamma distribution" << '\n';
  for (i=1; i<10; i++){
    cdf_value = 0.1*i;
#ifdef HAVE_BOOST
    Pecos::gamma_dist gamma1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::quantile(gamma1,cdf_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_cdf_gamma_Pinv(cdf_value, 3.,1.) << '\n';
#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "CDF values for gamma distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 1.0*i; 
#ifdef HAVE_BOOST
    Pecos::gamma_dist gamma1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::cdf(gamma1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_cdf_gamma_P(quantile_value, 3.,1.) << '\n';
#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "PDF values for gamma distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 1.0*i; 
#ifdef HAVE_BOOST
    Pecos::gamma_dist gamma1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::pdf(gamma1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_ran_gamma_pdf(quantile_value, 3.,1.) << '\n';
#endif // HAVE_GSL
  }
  
  PCout << std::setw(16) << "PDF values for weibull distribution" << '\n';
  for (i=1; i<5; i++){
    quantile_value = 0.1*i; 
#ifdef HAVE_BOOST
    Pecos::weibull_dist weibull1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::pdf(weibull1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_ran_weibull_pdf(quantile_value, 1.,3.) << '\n';
#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "CDF values for weibull distribution" << '\n';
  for (i=1; i<5; i++){
    quantile_value = 0.1*i; 
#ifdef HAVE_BOOST
    Pecos::weibull_dist weibull1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::cdf(weibull1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_cdf_weibull_P(quantile_value, 1.,3.) << '\n';
#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "CDF values for beta distribution" << '\n';
  for (i=0; i<11; i++){
    quantile_value = i*0.1; 
#ifdef HAVE_BOOST
    Pecos::beta_dist beta1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::cdf(beta1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_cdf_beta_P(quantile_value, 3.,1.) << '\n';
#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "Quantile values for beta distribution" << '\n';
  for (i=0; i<11; i++){
    cdf_value = 0.1*i; 
#ifdef HAVE_BOOST
    Pecos::beta_dist beta1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::quantile(beta1,cdf_value) <<'\n';
#endif //HAVE_BOOST
    //#ifdef HAVE_GSL
      // Pecos::ProbabilityTransformation this_trans;
    //Pecos::Real q = this_trans.cdf_beta_Pinv(cdf_value,3.0,1.0); 
    //PCout << std::setw(16) << "  GSL " << q  << '\n';
    //#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "PDF values for beta distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 0.1*i; 
#ifdef HAVE_BOOST
    Pecos::beta_dist beta1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::pdf(beta1,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_ran_beta_pdf(quantile_value, 3., 1.) << '\n';
#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "PDF values for F distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 0.1*i; 
#ifdef HAVE_BOOST
    Pecos::fisher_f_dist fisher(10.,4.);
    PCout << std::setw(16) << "Boost " << bmth::pdf(fisher,quantile_value) ;
#endif //HAVE_BOOST
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_ran_fdist_pdf(quantile_value, 10., 4.) << '\n';
#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "Quantiles for t-distribution" << '\n';
  for (i=1; i<10; i++){
    cdf_value = 0.1*i;
#ifdef HAVE_BOOST
    Pecos::students_t_dist tdist(10);
    PCout << std::setw(16) << "Boost " << bmth::quantile(tdist,cdf_value);
#endif
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_cdf_tdist_Pinv(cdf_value,10) << '\n';
#endif
  }

  PCout << std::setw(16) << "Quantiles for Chi-squared distribution" << '\n';
  for (i=1; i<10; i++){
    cdf_value = 0.1*i;
#ifdef HAVE_BOOST
    Pecos::chi_squared_dist chisq(10);
    PCout << std::setw(16) << "Boost " << bmth::quantile(chisq,cdf_value);
#endif
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_cdf_chisq_Pinv(cdf_value,10) << '\n';
#endif
  }

 PCout << std::setw(16) << "CDF values for exponential" << '\n';
  for (i=1; i<5; i++){
    quantile_value = 0.1*i;
#ifdef HAVE_BOOST
    Pecos::exponential_dist expon1(1.5);
    PCout << std::setw(16) << "Boost " << bmth::cdf(expon1,quantile_value);
#endif
#ifdef HAVE_GSL
    PCout << std::setw(16) << "  GSL " << gsl_cdf_exponential_P(quantile_value,(1/1.5)) << '\n';
#endif // HAVE_GSL
  }
}


// WJB: RNG testing -- eventually split out into its own rng_test.cpp file
#ifdef HAVE_BOOST
#include <boost/random/linear_congruential.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

// This is a reproducible simulation experiment.  See main().
void experiment(boost::mt19937& generator)
{
  using namespace boost;
  // Define a uniform random number distribution of integer values between
  // 1 and 6 inclusive.
  //typedef variate_generator<mt19937&, uniform_int<> > gen_type;
  variate_generator<mt19937&, uniform_int<> >
    die_gen(generator, uniform_int<>(1, 6));

  // If you want to use an STL iterator interface, use iterator_adaptors.hpp.
  generator_iterator<variate_generator<mt19937&, uniform_int<> > >
    die(&die_gen);
  for(int i = 0; i < 10; ++i)
    PCout << *die++ << " ";
  PCout << '\n';
}
#endif // HAVE_BOOST


int main(int argc, char* argv[])
{

  boost_test_dist bt;
  bt.print_comparison();

// WJB: now some RNG testing -- eventually split out into its own main/.cpp file
#ifdef HAVE_BOOST
  using namespace boost;
  PCout.setf(std::ios::fixed);

  // Define a RNG and initialize it with a reproducible seed.
  // (The seed is unsigned, otherwise the wrong overload may be selected
  // when using mt19937 as the base_generator_type.)
  mt19937 generator(42u);

  PCout << "10 samples of a uniform distribution in [0..1):\n";

  // Define a uniform random number distribution which produces "double"
  // values between 0 and 1 (0 inclusive, 1 exclusive).
  uniform_real<> uni_dist(0,1);
  variate_generator<mt19937&, uniform_real<> > uni(generator, uni_dist);

  // You can now retrieve random numbers from that distribution by means
  // of a STL Generator interface, i.e. calling the generator like a C-function
  // with no arguments.
  for(int i = 0; i < 10; ++i)
    PCout << uni() << '\n';

  generator.seed(static_cast<unsigned int>(std::time(0)));

  PCout << "\nexperiment: roll a die 10 times:\n";

  // You can save a generator's state by copy construction/assignment.
  mt19937 saved_generator = generator;

  // When calling other functions which take a generator or distribution
  // as a parameter, make sure to always call by reference (or pointer).
  // Calling by value invokes the copy constructor... (yadda-yadda -- UNDERSTOOD/BAD IDEA)
  experiment(generator);

  PCout << "redo the experiment to verify it:\n";
  experiment(saved_generator);

  // After that, both generators are equivalent
  //assert(generator == saved_generator);

  // as a degenerate case, you can set min = max for uniform_int
  uniform_int<> degen_dist(4,4);
  variate_generator<mt19937&, uniform_int<> > deg(generator, degen_dist);
  PCout << deg() << " " << deg() << " " << deg() << std::endl;
#endif // HAVE_BOOST

  return 0;
}
