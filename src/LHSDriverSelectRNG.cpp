/*  ______________________________________________________________________
    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */
 
//- Class:       BoostRNG_Monostate
//- Description: Wrapper for various implementations of Random Number Generators
//- Owner:       Laura Swiler, Dave Gay, and Bill Bohnhoff 
//- Checked by:
//- Version: $Id: Random.cpp 5721 2009-03-03 23:51:34Z wjbohnh $

#include "pecos_data_types.hpp"
#include "LHSDriver.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <cstdio>
#include <cstring>

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
{ //std::cout << "running Boost MT" << "\n";
  return BoostRNG_Monostate::random_num();
}
 
Real rnum2(void)
{
// clone of rnum1
 //std::cout << "running Boost MT" << "\n";
 return BoostRNG_Monostate::random_num2();
}

void Pecos::LHSDriver::seed(int seed)
{
	randomSeed = seed;
	if (BoostRNG_Monostate::random_num == BoostRNG_Monostate::random_num1)
		BoostRNG_Monostate::seed(seed);
	else
		lhs_setseed(&seed);
	}

void Pecos::LHSDriver::seed(int seed, const Pecos::String &unifGen)
{
	static int first = 1;
	static const char *s;
	randomSeed = seed;
	if (first) {
		s = std::getenv("DAKOTA_LHS_UNIFGEN");
		first = 0;
		}
	if (s) {
		if (!std::strcmp(s,"rnumlhs1"))
			goto use_rnum;
		else if (!std::strcmp(s, "mt19937"))
			goto use_mt;
		else if (*s) {
			std::fprintf(stderr, "Expected $DAKOTA_LHS_UNIFGEN"
				" to be \"rnumlhs1\" or \"mt19937\","
				" not \"%s\"\n", s);
			std::exit(1);
			}
		}
	if (unifGen == "mt19937" || unifGen == "") {
 use_mt:
		BoostRNG_Monostate::random_num  = BoostRNG_Monostate::random_num1;
		BoostRNG_Monostate::random_num2 = BoostRNG_Monostate::random_num1;
		BoostRNG_Monostate::seed(seed);
		}
	else {
 use_rnum:
		BoostRNG_Monostate::random_num  = (Rfunc)rnumlhs10;
		BoostRNG_Monostate::random_num2 = (Rfunc)rnumlhs20;
		lhs_setseed(&seed);
		}
	}

const char* Pecos::LHSDriver::unifGen()
{
	if (BoostRNG_Monostate::random_num == BoostRNG_Monostate::random_num1)
		return "mt19937";
	if (BoostRNG_Monostate::random_num == (Rfunc)rnumlhs10)
		return "mndp";
	return "unknown (bug?)";
	}

