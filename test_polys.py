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
        #Unsupported: assert self.poly3.eval(polynomial([-1,1])) == polynomial([0,0,0,0,1]), "Incorrectly composes polynomials"

    def testSub(self):
        assert -self.poly1 == polynomial([-1,0,-1,0,-1]), "Incorrectly negates polynomials"
        assert self.poly1 - self.poly1 == polynomial([0])
        assert self.poly3 - self.poly2 != polynomial([0])
        assert -self.poly3 == (-1)*self.poly3
        assert 3*self.poly3 - self.poly3 == self.poly3*2

if __name__ == "__main__":
    unittest.main()