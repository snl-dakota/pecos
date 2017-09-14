#ifndef TEST_OPTIONS_LIST_HPP
#define TEST_OPTIONS_LIST_HPP

#include "OptionsList.hpp"
#include <boost/shared_ptr.hpp>

#define TEST_FUNC_PROTOS(TYPE, SNAME)				  \
								  \
  OptionsList SNAME ## set_entry(OptionsList &opts,	  \
		    const std::string &name, const TYPE & item);  \


TEST_FUNC_PROTOS(int, int)
TEST_FUNC_PROTOS(double, double)
TEST_FUNC_PROTOS(std::string, string)
TEST_FUNC_PROTOS(OptionsList, optionslist)

#endif // TEST_OPTIONS_LIST_HPP

// based upon Vector.h Vector.cxx testVector.py at
// https://github.com/numpy/numpy/tree/master/tools/swig/test
