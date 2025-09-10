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

class FVector:
    def __init__(self, components:list):
        self.dim = len(components)
        self.components = components

    def dot(self,other):
        return sum([a*b for a,b in zip(self.components,other.components)],start=QR5(Fraction(0,1),Fraction(0,1)))

class FTensor:
    def __init__(self,components):
        self.components = components
        self.rank = self.get_rank()
        self.dims = self.get_dimensions()


    def get_dimensions(self):
        """
        >>> FTensor([1,2,3]).dims
        [3]

        >>> FTensor([[[1,2],[3,4]],[[5,6],[7,8]]]).dims
        [2, 2, 2]

        >>> FTensor(5).dims
        []

        :return:
        """
        return self._get_dimensions(self.components)

    def get_rank(self):
        """
        >>> FTensor([[1,2],[3,4]]).rank
        2
        >>> FTensor([[[1,2],[3,4]],[[5,6],[7,8]]]).rank
        3
        >>> FTensor(5).rank
        0

        :param components:
        :return:
        """
        return self._rank_rec(self.components)

    def _rank_rec(self,components):
        if isinstance(components,list):
            return self._rank_rec(components[0])+1
        else:
            return 0

    def _get_dimensions(self,components):
        if isinstance(components,list):
            result = self._get_dimensions(components[0])
            result.append(len(components))
            return result
        else:
            return list()

    def scale(self,alpha:QR5):
        """
        >>> np.random.seed(1234)
        >>> comps = np.array([QR5.random() for i in range(12)])
        >>> components = comps.reshape(3,4)

        >>> tensor = FTensor(components)
        >>> tensor.components.tolist()
        [[(5+9*r5), (-4+2*r5), (5+7*r5), (-1+1*r5)], [(2+6*r5), (-5+6*r5), (-1+5*r5), (8+6*r5)], [(2-5*r5), (-8-4*r5), (-7-3*r5), (1-10*r5)]]

        >>> tensor = tensor.scale(QR5(Fraction(2,1),Fraction(0,1)))
        >>> tensor.components.tolist()
        [[(10+18*r5), (-8+4*r5), (10+14*r5), (-2+2*r5)], [(4+12*r5), (-10+12*r5), (-2+10*r5), (16+12*r5)], [(4-10*r5), (-16-8*r5), (-14-6*r5), (2-20*r5)]]

        :param alpha:
        :return:
        """
        return FTensor(self.components*alpha)

    def __mul__(self,other):
        """
        tensor product between two tensors
        >>> np.random.seed(1234)
        >>> comps = np.array([QR5.random() for i in range(12)])
        >>> components = comps.reshape(3,4)
        >>> comps2 = np.array([QR5.random() for i in range(12)])
        >>> components2 = comps2.reshape(3,2,2)
        >>> product = FTensor(components)*FTensor(components2)
        >>> product.components.tolist()
        [[[[[(40-4*r5), (-285+19*r5)], [(365-27*r5), (-395-27*r5)]], [[(410+54*r5), (320+44*r5)], [(425+81*r5), (-15-27*r5)]], [[(185+29*r5), (170+78*r5)], [(40-80*r5), (290-10*r5)]]], [[[(14-6*r5), (-94+40*r5)], [(122-52*r5), (-98+40*r5)]], [[(86-34*r5), (66-26*r5)], [(74-28*r5), (12-6*r5)]], [[(36-14*r5), (2+2*r5)], [(60-28*r5), (90-38*r5)]]], [[[(30-2*r5), (-215+7*r5)], [(275-11*r5), (-305-31*r5)]], [[(320+52*r5), (250+42*r5)], [(335+73*r5), (-15-21*r5)]], [[(145+27*r5), (140+64*r5)], [(20-60*r5), 220]]], [[[(6-2*r5), (-41+13*r5)], [(53-17*r5), (-47+11*r5)]], [[(44-8*r5), (34-6*r5)], [(41-5*r5), (3-3*r5)]], [[(19-3*r5), (8+4*r5)], [(20-12*r5), (40-12*r5)]]]], [[[[(28-4*r5), (-198+22*r5)], [(254-30*r5), (-266-6*r5)]], [[(272+24*r5), (212+20*r5)], [(278+42*r5), (-6-18*r5)]], [[(122+14*r5), (104+48*r5)], [(40-56*r5), (200-16*r5)]]], [[[(35-11*r5), (-240+71*r5)], [(310-93*r5), (-280+57*r5)]], [[(265-39*r5), (205-29*r5)], [(250-21*r5), (15-18*r5)]], [[(115-14*r5), (55+27*r5)], [(110-70*r5), (235-65*r5)]]], [[[(26-6*r5), (-181+37*r5)], [(233-49*r5), (-227+19*r5)]], [[(224-4*r5), (174-2*r5)], [(221+11*r5), (3-15*r5)]], [[(99+1*r5), (68+32*r5)], [(60-52*r5), (180-32*r5)]]], [[[(22+2*r5), (-162-20*r5)], [(206+24*r5), (-254-60*r5)]], [[(278+78*r5), (218+62*r5)], [(302+96*r5), (-24-18*r5)]], [[(128+38*r5), (146+66*r5)], [(-20-44*r5), (170+26*r5)]]]], [[[[(-27+7*r5), (187-44*r5)], [(-241+58*r5), (229-28*r5)]], [[(-223+13*r5), (-173+9*r5)], [(-217-2*r5), (-6+15*r5)]], [[(-98+3*r5), (-61-29*r5)], [(-70+54*r5), (-185+39*r5)]]], [[[(-12-4*r5), (92+32*r5)], [(-116-40*r5), (164+64*r5)]], [[(-188-76*r5), (-148-60*r5)], [(-212-88*r5), (24+12*r5)]], [[(-88-36*r5), (-116-52*r5)], [(40+24*r5), (-100-36*r5)]]], [[[(-8-4*r5), (63+31*r5)], [(-79-39*r5), (121+57*r5)]], [[(-142-66*r5), (-112-52*r5)], [(-163-75*r5), (21+9*r5)]], [[(-67-31*r5), (-94-42*r5)], [(40+16*r5), (-70-34*r5)]]], [[[(-51+11*r5), (356-67*r5)], [(-458+89*r5), (452-29*r5)]], [[(-449-1*r5), (-349-3*r5)], [(-446-31*r5), (-3+30*r5)]], [[(-199-6*r5), (-143-67*r5)], [(-110+102*r5), (-355+57*r5)]]]]]
        >>> product.dims

        :param other:
        :return:
        """
        return FTensor(np.tensordot(self.components,other.components,axes=0))

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

def detect_cells(faces,edges):
    # find to faces with a common edge
    pairs = []
    for i in range(len(faces)):
        for j in range(i+1,len(faces)):
            common = set(faces[i]).intersection(set(faces[j]))
            if len(common)==1:
                pairs.append((i,j))

    print(len(pairs),"pairs")

    # merge to pairs of faces into a cell
    cells = set()
    for i in range(len(pairs)):
        for j in range(i+1,len(pairs)):
            pair_faces = set(pairs[i]+pairs[j])
            if len(pair_faces)==4:
                cell_edges = set()
                for face_index in pair_faces:
                     cell_edges= cell_edges.union(set([i for i in faces[face_index]]))

                if len(cell_edges)==6:
                    cells.add(tuple(sorted(pair_faces)))

    print(len(cells),"cells")
    return list(cells)

def detect_faces(edges):
    faces = []
    for i in range(len(edges)):
        for j in range(i+1,len(edges)):
            for k in range(j+1,len(edges)):
                edge1 = set(edges[i])
                edge2 = set(edges[j])
                edge3 = set(edges[k])
                all = edge1.union(edge2,edge3)
                if len(all)==3:
                    faces.append([i,j,k])
    return faces

def save(data,filename):
    with open(filename,"w") as f:
        for d in data:
            f.write(str(d)+"\n")


def read(filename):
    with open(filename,"r") as f:
        data = []
        for line in f:
            data.append(eval(line))
    return data

def compute_equation(cell,faces,edges,vectors):
    # grap the four vertices of the cell
    cell_faces = [faces[i] for i in cell]
    cell_edges = [edges[j] for c in cell_faces for j in c]
    cell_vertices = set([tuple(vectors[j]) for e in cell_edges for j in e])
    matrix = [list(v) for v in cell_vertices]

    center = [sum(e,QR5.from_integers(0,1,0,1))*QR5.from_integers(1,4,0,1) for e in zip(*matrix)]
    basis = []




if __name__ == '__main__':
    elements = list(generate_group())
    vectors = get_4D_vectors()

    # edges in terms of vertices

    # edges = detect_edges(elements)
    # save(edges,"edges.data")
    edges = read("edges.data")
    print("number of edges ",len(edges),edges)

    # faces in terms of edges

    #faces = detect_faces(edges)
    #save(faces,"faces.data")
    faces = read("faces.data")
    print("number of faces ",len(faces),faces)

    # cells in terms of faces
    # cells = detect_cells(faces,edges)
    # save(cells,"cells.data")
    cells = read("cells.data")

    print(len(cells))
    print(compute_equation(cells[0],faces,edges,vectors))
