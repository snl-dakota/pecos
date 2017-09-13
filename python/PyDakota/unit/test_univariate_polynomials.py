import unittest, numpy, os
from PyDakota.univariate_polynomials import *
class TestUnivariatePolynomials(unittest.TestCase):
    def setUp( self ):
        pass

    def test_legendre_polynomial(self):
        poly = LegendreOrthogPolynomial()
        poly.type1_value(0.,2)

if __name__ == '__main__':
    unittest.main()
