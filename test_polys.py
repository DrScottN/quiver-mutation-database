from polynomial import *
import unittest
import random

class PolynomialTestCase(unittest.TestCase):
    def setUp(self):
        self.poly1 = polynomial([1,0,1,0,1,0])
        self.poly2 = polynomial([0,0,0,0,0,1])
        self.poly3 = polynomial([1,4,6,4,1])

    def testEval(self):
        assert self.poly1.eval([0])==1
        assert self.poly2.eval([0])==0
        assert self.poly3.eval([0])==1

        assert self.poly1.eval([2])==1 + 2**2 + 2**4
        assert self.poly2.eval([2])==2**5
        assert self.poly3.eval([2])==(3)**4

        assert self.poly1.eval([-2])==1 + 2**2 + 2**4
        assert self.poly2.eval([-2])==-2**5
        assert self.poly3.eval([-2])==(-1)**4

    def testSum(self):
        assert self.poly1 + self.poly2 == self.poly2 + self.poly1, "addition not commutative"
        assert self.poly1 + self.poly1 == 2*self.poly1
        assert (self.poly1 + self.poly2).eval([0]) == 1
        assert (self.poly1 + self.poly2).eval([2]) == 1 + 2**2 + 2**4 + 2**5
        assert (self.poly3 + self.poly2).eval([2]) == 3**4 + 2**5
        assert self.poly1 + 2 == polynomial([3,0,1,0,1,0]), f"Incorrectly add constants."
        assert self.poly2 + 2 == polynomial([2,0,0,0,0,1]), f"Incorrectly add constants."
        assert self.poly3 + 2 == polynomial([3,4,6,4,1]), f"Incorrectly add constants."
        assert self.poly1 + self.poly2 + self.poly3 == self.poly3 + self.poly1 + self.poly2
        assert self.poly1 + self.poly2 + self.poly3 == self.poly1 + (self.poly2 + self.poly3)

    def testProd(self):
        assert (self.poly1*self.poly2).eval([0]) == 0
        assert self.poly3*self.poly2 == self.poly2*self.poly3, "mult not commutative"
        assert self.poly1.eval([3])*self.poly2.eval([3]) == (self.poly1 * self.poly2).eval([3])
        assert self.poly1 + self.poly1 == self.poly1*2
        assert self.poly1 * self.poly2 * self.poly3 == self.poly3 * self.poly1 * self.poly2
        assert self.poly1 * self.poly2 * self.poly3 == self.poly1 * (self.poly2*self.poly3)

    def testEq(self):
        assert self.poly1 != self.poly2, "Incorrectly identifies different polynomials"
        assert self.poly1 != self.poly3, "Incorrectly identifies different polynomials"
        assert self.poly1 == polynomial([1]) + polynomial([0,0,1]) + polynomial([0,0,0,0,1])
        self.poly1 + 2
        self.poly1 * 2
        assert self.poly1 == polynomial([1,0,1,0,1])
        assert self.poly1 +2 != self.poly1
        assert self.poly3.eval([polynomial([-1,1])]) == polynomial([0,0,0,0,1]), "Incorrectly composes polynomials"

    def testSub(self):
        assert -self.poly1 == polynomial([-1,0,-1,0,-1]), "Incorrectly negates polynomials"
        assert self.poly1 - self.poly1 == polynomial([0])
        assert self.poly3 - self.poly2 != polynomial([0])
        assert -self.poly3 == (-1)*self.poly3
        assert 3*self.poly3 - self.poly3 == self.poly3*2


class PolynomialVarnamesTestCase(unittest.TestCase):
    def setUp(self):
        self.oney = polynomial([1], ('y',))
        self.onez = polynomial([1], ('z',))
        
        self.fy = polynomial([0,1], ('y',))
        self.fz = polynomial([0,1], ('z',))

        self.gz = polynomial([0,1,2,1], ('z',))
        self.gw = polynomial([0,1,2,1], ('w',))

    def testEqVarnames(self):
        assert self.oney == self.onez, "Incorrectly distinguishes constant polys"
        assert self.oney == 1
        assert 10*self.oney == 10
        assert self.fy*self.oney == self.fy
        assert self.fy != self.fz
        assert self.gz != self.gw
        assert self.fy*0 +self.fz == self.fz
    
    def testSumVarnames(self):
        assert self.fy + self.fz == self.fz + self.fy
        assert self.gw + self.fz == self.fz + self.gw
        assert self.fy + self.oney == self.fy + self.onez
        assert 3*self.fy - 2*self.fy == self.fy
    
    def testProdVarnames(self):
        assert self.fy * self.fy * self.fy == polynomial([0,0,0,1], ('y',))
        assert self.fy * self.fz == self.fz * self.fy
        assert self.fy*(self.gw + self.gz) == self.fy*self.gw + self.gz*self.fy
        assert 3*self.gw == self.gw + self.gw + self.gw


class PolynomialMultivarTestCase(unittest.TestCase):
    def setUp(self):
        self.onexyz = polynomial([(1,0,0,0)], ('x','y','z'))
        self.onezyx = polynomial([(1,0,0,0)], ('z','y','x'))
        
        self.fxy = polynomial([(0,0,0),(2,0,1), (1,1,0)], ('x','y'))
        self.fyx = polynomial([(0,0,0),(2,1,0), (1,0,1)], ('y','x'))

    def testEqMulti(self):
        assert self.onexyz == self.onezyx
        assert self.fxy == self.fyx, f"Incorrectly distinguishes {str(self.fxy), str(self.fyx)}"
        assert self.fxy + self.fyx == 2*self.fxy


if __name__ == "__main__":
    unittest.main()