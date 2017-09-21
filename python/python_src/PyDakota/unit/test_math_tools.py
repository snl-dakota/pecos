import unittest, numpy
from PyDakota.math_tools import *
class TestMathTools(unittest.TestCase):
    def setUp(self):
        pass

    def test_cartesian_product(self):
        """
        """
                  
        # test when num elems = 1
        s1 = numpy.arange( 0, 3 )
        s2 = numpy.arange( 3, 5 )

        sets = numpy.array( [[0,3], [1,3], [2,3], [0,4], 
                                [1,4], [2,4]], numpy.int )
        output_sets = cartesian_product( [s1,s2], 1 )
        assert numpy.array_equal( output_sets.T, sets )
          
        # test when num elems > 1
        s1 = numpy.arange( 0, 6 )
        s2 = numpy.arange( 6, 10 )
        
        sets = numpy.array( [[ 0, 1, 6, 7], [ 2, 3, 6, 7],
                                [ 4, 5, 6, 7], [ 0, 1, 8, 9],
                                [ 2, 3, 8, 9], [ 4, 5, 8, 9]], numpy.int )
        output_sets = cartesian_product( [s1,s2], 2 )
        assert numpy.array_equal( output_sets.T, sets )
        
if __name__ == '__main__':
    unittest.main()
