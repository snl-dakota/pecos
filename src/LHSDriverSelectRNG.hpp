/*  ______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */
 
//- Class:       BoostRNG_Monostate
//- Description: Wrapper for various implementations of Random Number Generators
//- Owner:       Laura Swiler, Dave Gay, and Bill Bohnhoff 
//- Checked by:
//- Version: $Id: Random.cpp 5721 2009-03-03 23:51:34Z wjbohnh $

#include "pecos_data_types.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

// WJB - 5/3/09 - ToDo: address proper namespace utilization next iteration
using Pecos::Real;
typedef Real (*Rfunc)();


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

 static Real random_num1()
 { return uniMT(); };

 static Rfunc random_num;
 static Rfunc random_num2;
};


unsigned int BoostRNG_Monostate::rngSeed(41u); // 41 used in the Boost examples
boost::mt19937 BoostRNG_Monostate::rnumGenerator( BoostRNG_Monostate::seed() );
boost::uniform_real<> BoostRNG_Monostate::uniDist(0, 1);

#define rnum1 FC_FUNC(rnumlhs1,RNUMLHS1)
#define rnum2 FC_FUNC(rnumlhs2,RNUMLHS2)
#define rnumlhs10 FC_FUNC(rnumlhs10,RNUMLHS10)
#define rnumlhs20 FC_FUNC(rnumlhs20,RNUMLHS20)
#define lhs_setseed FC_FUNC(lhssetseed,LHSSETSEED)

extern "C" Real rnum1(void), rnum2(void), rnumlhs10(void), rnumlhs20(void);
extern "C" void lhs_setseed(int*);

Real (*BoostRNG_Monostate::random_num)()  = BoostRNG_Monostate::random_num1;
Real (*BoostRNG_Monostate::random_num2)() = BoostRNG_Monostate::random_num1;

boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
  BoostRNG_Monostate::uniMT(BoostRNG_Monostate::rnumGenerator,
                            BoostRNG_Monostate::uniDist);
 
extern "C" void set_boost_rng_seed(unsigned int rng_seed)
{
  BoostRNG_Monostate::seed(rng_seed);
}
 
Real rnum1(void)
{
  //std::cout << "running Boost MT" << "\n";
  return BoostRNG_Monostate::random_num();
}
 
Real rnum2(void)
{
  // clone of rnum1
  //std::cout << "running Boost MT" << "\n";
  return BoostRNG_Monostate::random_num2();
}
