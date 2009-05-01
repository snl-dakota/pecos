/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef BOOST_TEST_DIST_HPP
#define BOOST_TEST_DIST_HPP

#include "pecos_stat_util.hpp"


/// Class for testing Boost vs. GSL statistical functions


class boost_test_dist
{
public:

  /// default constructor
  boost_test_dist();
  /// destructor
  ~boost_test_dist();

  void print_comparison();
};

#endif
