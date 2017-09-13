#!/usr/bin/env python

import sys
# This driver requires the path to the compiled SWIG module for now:
sys.path.append(sys.argv[1])

import swig_example

four_factorial = swig_example.fact(4)
five_mod_two = swig_example.my_mod(5,2)

print "4! = ", four_factorial
print "5%2 = ", five_mod_two

assert(four_factorial == 24)
assert(five_mod_two == 1)
