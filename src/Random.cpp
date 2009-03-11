/*  ______________________________________________________________________
    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */
 
//- Class:       Random
//- Description: Wrapper for various implementations of Random Number Generators 
//- Owner:       Laura Swiler and Bill Bohnhoff 
//- Checked by:
//- Version: $Id: Random.C 5721 2009-03-03 23:51:34Z wjbohnh $
 
// Are we heading toward the need for a "lhs_config.h"??
#include "pecos_config.h"

#ifdef HAVE_BOOST
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
//#include <boost/generator_iterator.hpp> // WJB: probably not needed
 

class BoostRNG_Monostate
{
 static unsigned int rngSeed;
 static boost::mt19937 rnumGenerator;
 static boost::uniform_real<> uniDist;
 static boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uniMT;

public:

 // WJB: consider adding a counter to ensure set_seed is called only ONCE!
 static void seed(unsigned int rng_seed)
 { rngSeed = rng_seed; rnumGenerator.seed(rng_seed); };

 static unsigned int seed()
 { return rngSeed; }

 static double random_num1()
 { return uniMT(); };
};

unsigned int BoostRNG_Monostate::rngSeed(41u); // WJB: is 41 OK for starters??
boost::mt19937 BoostRNG_Monostate::rnumGenerator( BoostRNG_Monostate::seed() );
boost::uniform_real<> BoostRNG_Monostate::uniDist(0, 1);

boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
  BoostRNG_Monostate::uniMT(BoostRNG_Monostate::rnumGenerator,
                            BoostRNG_Monostate::uniDist);
 
extern "C" void set_boost_rng_seed(unsigned int rng_seed)
{
  BoostRNG_Monostate::seed(rng_seed);
}


#define rnum1 FC_FUNC(rnumlhs1,RNUMLHS1)
#define rnum2 FC_FUNC(rnumlhs2,RNUMLHS2)
extern "C" double rnum1(void), rnum2(void);
 
double rnum1(void)
{ std::cout << "running Boost MT" << "\n";
  return BoostRNG_Monostate::random_num1();
}
 
double rnum2(void)
{
// clone of rnum1
 std::cout << "running Boost MT" << "\n";
 return BoostRNG_Monostate::random_num1();
}

#endif // HAVE_BOOST
 
