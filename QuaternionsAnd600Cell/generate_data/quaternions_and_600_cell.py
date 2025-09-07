# this file parallels the computations of the corresponding mathematica file
# for the exact generation of the c600 we need to define quaternions over field extensions
from fractions import Fraction

# let's start with field extension Q[r5], which are the rational numbers and the  root of five adjoined

class elem:
    def __init__(self,x:Fraction,y:Fraction):
        self.x=x
        self.y=y

    @classmethod
    def from_integers(cls, a:int, b:int, c:int, d:int):
        return cls(Fraction(a,b), Fraction(c,d))

    def __str__(self):
        return "("+str(self.x)+"+"+str(self.y)+"*r5)"

    def __repr__(self):
        return self.__str__()

    def __mul__(self, other):
        return elem(self.x*other.x+5*self.y*other.y, self.x*other.y+self.y*other.x)


print(elem.from_integers(1,1,1,1)*elem.from_integers(1,2,0,1))