from PyDakota.swig_examples import dot
import numpy
vec1=[1,2,3]
vec2=[4,5,6]
assert dot.dot(vec1,vec2)==numpy.dot(vec1,vec2)
