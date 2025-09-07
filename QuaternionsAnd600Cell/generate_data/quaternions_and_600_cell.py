# this file parallels the computations of the corresponding mathematica file
# for the exact generation of the c600 we need to define quaternions over field extensions
from fractions import Fraction
import numpy as np


# let's start with field extension Q[r5], which are the rational numbers and the  root of five adjoined

class elem:
    def __init__(self,x:Fraction,y:Fraction):
        """
        >>> elem(Fraction(1, 2), Fraction(3, 4))
        (1/2+3/4*r5)

        :param x:
        :param y:
        """
        self.x=x
        self.y=y
        #self.zero = elem(Fraction(0,1),Fraction(0,1))
        #self.one = elem(Fraction(1,0),Fraction(0,1))

    @classmethod
    def from_integers(cls, a:int, b:int, c:int, d:int):
        """
        >>> elem.from_integers(1,2,3,4)
        (1/2+3/4*r5)
        """

        return cls(Fraction(a,b), Fraction(c,d))

    def __str__(self):
        if self.y>0:
            return "("+str(self.x)+"+"+str(self.y)+"*r5)"
        elif self.y<0:
            return "("+str(self.x)+"-"+str(Fraction(np.abs(self.y.numerator),np.abs(self.y.denominator)))+"*r5)"
        else:
            return str(self.x)

    def __repr__(self):
        return self.__str__()

    def __mul__(self, other):
        """
        >>> elem(Fraction(1, 2), Fraction(3, 4))*elem(Fraction(2,1), Fraction(4,3))
        (6+13/6*r5)

        >>> np.random.seed(1234)
        >>> z = elem.random()
        >>> w = elem.random()
        >>> z*w
        (70-26*r5)

        :param other:
        :return:
        """
        return elem(self.x*other.x+5*self.y*other.y, self.x*other.y+self.y*other.x)

    def __add__(self,other):
        """
        >>> np.random.seed(1234)
        >>> z = elem.random()
        >>> w = elem.random()
        >>> z+w
        (1+11*r5)

        :param other:
        :return:
        """
        return elem(self.x+other.x,self.y+other.y)

    def __sub__(self,other):
        """
        >>> np.random.seed(1234)
        >>> z = elem.random()
        >>> w = elem.random()
        >>> z-w
        (9+7*r5)

        :param other:
        :return:
        """
        return elem(self.x-other.x,self.y-other.y)

    def conj(self):
        return elem(self.x,-self.y)

    def norm(self):
        return self.x**2-5*self.y**2

    def __neg__(self):
        return elem(-self.x,-self.y)

    def __truediv__(self,other):
        """
        >>> np.random.seed(1234)
        >>> z = elem.random()
        >>> w = elem.random()
        >>> z/w
        (55/2+23/2*r5)

        >>> z = elem.random()
        >>> w = elem.random()
        >>> z/w*w,z

        :param other:
        :return:
        """
        return elem(1/other.norm(),Fraction(0,1))*(self*other.conj())

    @classmethod
    def random(cls,range=10):
        x=Fraction(np.random.randint(-range,range),1)
        y=Fraction(np.random.randint(-range,range),1)
        return elem(x,y)

class FComplex:
    def __init__(self,re:elem,im:elem):
        self.re=re
        self.im=im
    def __str__(self):
        return str(self.re)+"+"+str(self.im)+"*i"
    def __repr__(self):
        return self.__str__()
    def __mul__(self, other):
        return FComplex(self.re*other.re-self.im*other.im,self.re*other.im+self.im*other.re)
    def __add__(self,other):
        return FComplex(self.re+other.re,self.im+other.im)
    def __sub__(self,other):
        return FComplex(self.re-other.re,self.im-other.im)
    def conj(self):
        return FComplex(self.re,self.im.neg)
    def norm(self):
        return self*self.conj()
    def __neg__(self):
        return FComplex(-self.re,-self.im)
    def __truediv__(self,other):
        return 1/other.norm()*(self*other.conj())
    @classmethod
    def random(cls,range=10):
        x=elem.random(range)
        y=elem.random(range)
        return FComplex(x,y)

np.random.seed(0)
for i in range(10):
    z =FComplex.random(range=10)
    print(z)