# this file parallels the computations of the corresponding mathematica file
# for the exact generation of the c600 we need to define quaternions over field extensions
from fractions import Fraction
import numpy as np


# let's start with field extension Q[r5], which are the rational numbers and the  root of five adjoined

class QR5:


    def __init__(self,x:Fraction,y:Fraction):
        """
        >>> QR5(Fraction(1, 2), Fraction(3, 4))
        (1/2+3/4*r5)

        :param x:
        :param y:
        """
        self.x=x
        self.y=y



    @classmethod
    def from_integers(cls, a:int, b:int, c:int, d:int):
        """
        >>> QR5.from_integers(1,2,3,4)
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
        >>> QR5(Fraction(1, 2), Fraction(3, 4))*QR5(Fraction(2,1), Fraction(4,3))
        (6+13/6*r5)

        >>> np.random.seed(1234)
        >>> z = QR5.random()
        >>> w = QR5.random()
        >>> z*w
        (70-26*r5)

        :param other:
        :return:
        """
        return QR5(self.x * other.x + 5 * self.y * other.y, self.x * other.y + self.y * other.x)

    def __add__(self,other):
        """
        >>> np.random.seed(1234)
        >>> z = QR5.random()
        >>> w = QR5.random()
        >>> z+w
        (1+11*r5)

        :param other:
        :return:
        """
        return QR5(self.x + other.x, self.y + other.y)

    def __sub__(self,other):
        """
        >>> np.random.seed(1234)
        >>> z = QR5.random()
        >>> w = QR5.random()
        >>> z-w
        (9+7*r5)

        :param other:
        :return:
        """
        return QR5(self.x - other.x, self.y - other.y)

    def conj(self):
        """
        >>> QR5(Fraction(1, 2), Fraction(3, 4)).conj()
        (1/2-3/4*r5)

        :return:
        """
        return QR5(self.x, -self.y)

    def norm(self):
        """
        >>> QR5(Fraction(1, 2), Fraction(3, 4)).norm()
        Fraction(-41, 16)

        :return:
        """
        return self.x**2-5*self.y**2

    def __neg__(self):
        """
        >>> -QR5.from_integers(-1,1,2,-3)
        (1+2/3*r5)

        :return:
        """
        return QR5(-self.x, -self.y)

    def __truediv__(self,other):
        """
        >>> np.random.seed(1234)
        >>> z = QR5.random()
        >>> w = QR5.random()
        >>> z/w
        (55/2+23/2*r5)

        >>> z = QR5.random()
        >>> w = QR5.random()
        >>> z/w*w,z
        ((5+7*r5), (5+7*r5))

        :param other:
        :return:
        """
        return QR5(1 / other.norm(), Fraction(0, 1))*(self * other.conj())

    @classmethod
    def random(cls,range=10):
        x=Fraction(np.random.randint(-range,range),1)
        y=Fraction(np.random.randint(-range,range),1)
        return QR5(x, y)

    def __eq__(self, other):
        return self.x==other.x and self.y==other.y

    def __hash__(self):
        return hash((self.x,self.y))

    def to_float(self):
        return float(self.x)+float(self.y)*np.sqrt(5)

class FComplex:
    def __init__(self, re:QR5, im:QR5):
        self.re=re
        self.im=im

    def __str__(self):
        """
        >>> str(FComplex(QR5(Fraction(1, 2), Fraction(3, 4)), QR5(Fraction(2, 1), Fraction(4, 3))))
        '(1/2+3/4*r5)+(2+4/3*r5)*i'

        :return:
        """
        return str(self.re)+"+"+str(self.im)+"*i"

    def __repr__(self):
        return self.__str__()
    def __mul__(self, other):
        """
        >>> FComplex(QR5(Fraction(1, 2), Fraction(3, 4)), QR5(Fraction(5, 6), Fraction(7, 8))) * FComplex(QR5(Fraction(9, 8), Fraction(7, 6)), QR5(Fraction(5, 4), Fraction(3, 2)))
        (-8/3-11/12*r5)+(295/24+2099/576*r5)*i

        >>> u = FComplex.random(10)
        >>> v = FComplex.random(10)
        >>> w = FComplex.random(10)
        >>> w*(u+v)-w*u-w*v
        0+0*i

        :param other:
        :return:
        """
        return FComplex(self.re*other.re-self.im*other.im,self.re*other.im+self.im*other.re)

    def __add__(self,other):
        """
        >>> FComplex(QR5(Fraction(1, 2), Fraction(3, 4)), QR5(Fraction(5, 6), Fraction(7, 8))) + FComplex(QR5(Fraction(9, 8), Fraction(7, 6)), QR5(Fraction(5, 4), Fraction(3, 2)))
        (13/8+23/12*r5)+(25/12+19/8*r5)*i

        :param other:
        :return:
        """
        return FComplex(self.re+other.re,self.im+other.im)

    def __sub__(self,other):
        """
        >>> FComplex(QR5(Fraction(1, 2), Fraction(3, 4)), QR5(Fraction(5, 6), Fraction(7, 8))) - FComplex(QR5(Fraction(9, 8), Fraction(7, 6)), QR5(Fraction(5, 4), Fraction(3, 2)))
        (-5/8-5/12*r5)+(-5/12-5/8*r5)*i

        :param other:
        :return:
        """
        return FComplex(self.re-other.re,self.im-other.im)


    def conj(self):
        return FComplex(self.re,-self.im)


    def norm(self):
        """
        >>> FComplex(QR5(Fraction(1, 2), Fraction(3, 4)), QR5(Fraction(5, 6), Fraction(7, 8))).norm()
        Fraction(10998241, 331776)

        :return:
        """
        return (self*self.conj()).re.norm()

    def __neg__(self):
        return FComplex(-self.re,-self.im)

    def __truediv__(self,other):
        """
        >>> u = FComplex.random(10)
        >>> v = FComplex.random(10)
        >>> w = FComplex.random(10)
        >>> (u-v)/w-u/w+v/w
        0+0*i

        :param other:
        :return:
        """
        return FComplex(QR5(1/other.norm(),Fraction(0,1)),QR5.from_integers(0,1,0,1))*(self*other.conj())

    def __eq__(self, other):
        return self.re==other.re and self.im==other.im

    @classmethod
    def random(cls,range=10):
        x=QR5.random(range)
        y=QR5.random(range)
        return FComplex(x,y)

    def __hash__(self):
        return hash((self.re,self.im))

class FQuaternion:
    def __init__(self, a:FComplex, b:FComplex):
        self.a = a
        self.b = b

    def __str__(self):
        return str(self.a.re)+"+"+str(self.a.im)+"*i+"+str(self.b.re)+"*j"+"+"+str(self.b.im)+"*k"
    def __repr__(self):
        return self.__str__()

    @classmethod
    def from_vector(cls,a:QR5,b:QR5,c:QR5,d:QR5):
        return FQuaternion(FComplex(a,b),FComplex(c,d))

    def __add__(self,other):
        return FQuaternion(self.a+other.a,self.b+other.b)

    def __sub__(self,other):
        return FQuaternion(self.a-other.a,self.b-other.b)

    def __mul__(self,other):
        """
        Caley-Dickson construction

        :param other:
        :return:
        """
        return FQuaternion(self.a*other.a-other.b.conj()*self.b,other.b*self.a+self.b*other.a.conj())

    def __neg__(self):
        return FQuaternion(-self.a,-self.b)

    def __eq__(self, other):
        """
        Check algebra of quaternions

        >>> one = QR5.from_integers(1,1,0,1)
        >>> zero = QR5.from_integers(0,1,0,1)
        >>> q_one=FQuaternion.from_vector(one,zero,zero,zero)
        >>> q_i=FQuaternion.from_vector(zero,one,zero,zero)
        >>> q_j=FQuaternion.from_vector(zero,zero,one,zero)
        >>> q_k=FQuaternion.from_vector(zero,zero,zero,one)
        >>> assert(q_i*q_j==q_j*q_i.conj())
        >>> assert(q_i*q_k==q_k*q_i.conj())
        >>> assert(q_j*q_k==q_k*q_j.conj())
        >>> assert(q_i*q_i==-q_one)
        >>> assert(q_j*q_j==-q_one)
        >>> assert(q_k*q_k==-q_one)

        :param other:
        :return:
        """
        return self.a==other.a and self.b==other.b

    def conj(self):
        return FQuaternion(self.a.conj(),-self.b)

    def __hash__(self):
        return hash((self.a,self.b))

    def norm(self):
        return (self*self.conj()).a.re

    def to_vector(self):
        return [self.a.re,self.a.im,self.b.re,self.b.im]


def generate_group():
    # generate 120 elements of the 600 cell
    omega=gen_a = FQuaternion.from_vector(QR5.from_integers(-1,2,0,1),QR5.from_integers(1,2,0,1),QR5.from_integers(1,2,0,1),QR5.from_integers(1,2,0,1))
    q=gen_b = FQuaternion.from_vector(QR5.from_integers(0,1,0,1),QR5.from_integers(1,2,0,1),QR5.from_integers(1,4,1,4),QR5.from_integers(-1,4,1,4))
    print(gen_a)
    print(gen_b)
    # check group constraints
    one = QR5.from_integers(1, 1, 0, 1)
    zero = QR5.from_integers(0, 1, 0, 1)
    q_one = FQuaternion.from_vector(one, zero, zero, zero)
    assert(q*q*q*q==q_one)
    assert(omega*omega*omega==q_one)
    assert(q*omega*q*omega*q*omega*q*omega*q*omega==q_one)

    elements = {gen_a, gen_b}
    new_elements = {gen_a, gen_b}
    generators = {gen_a, gen_b}

    while len(new_elements)>0:
        next_elements = set()
        for e in new_elements:
            for g in generators:
                element = e*g
                if element not in elements:
                    next_elements.add(element)
                    elements.add(element)
        new_elements = next_elements

    assert(len(elements)==120)
    return elements

def get_4D_vectors():
    elements = generate_group()
    vectors = []
    for element in elements:
        vectors.append(element.to_vector())
    return vectors

def detect_edges(elements):
    edges = []
    min_dist = np.inf
    min = Fraction(0,1)
    for i in range(len(elements)):
        for j in range(i+1,len(elements)):
            norm = (elements[i]-elements[j]).norm()
            dist = norm.to_float()
            if dist<min_dist:
                min_dist = dist
                min = norm

    print("minimal edge length: ",min)

    for i in range(len(elements)):
        for j in range(i+1,len(elements)):
            dist = (elements[i]-elements[j]).norm().to_float()
            if dist==min_dist:
                edges.append((i,j))

    return edges

if __name__ == '__main__':
    elements = list(generate_group())
    vectors = get_4D_vectors()
    edges = detect_edges(elements)
    print(len(edges))
    print(edges)