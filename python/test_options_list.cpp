#include "test_options_list.hpp"

#define TEST_FUNCS(TYPE, SNAME)					\
								\
  OptionsList SNAME ## set_entry(OptionsList &opts,	\
		    const std::string &name, const TYPE &item){	\
    opts.set(name,item);					\
    return opts;						\
  }								

TEST_FUNCS(int, int)
TEST_FUNCS(double, double)
TEST_FUNCS(std::string, string)
TEST_FUNCS(OptionsList, optionslist)
