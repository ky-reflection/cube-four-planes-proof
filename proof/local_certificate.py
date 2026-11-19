"""Exact data for the local minimum of the four-plane cube problem.

Run with Python 3. The cube is [-1/2,1/2]^3. The four moving planes, in
sign-mask order, have equations

    (1+u1)x + (1+u2)y + (1-u1-u2)z = d0
    sqrt(2)z - (1+v31)x - (1+v32)y = d3
    sqrt(2)y - (1+v21)x - (1+v23)z = d2
    sqrt(2)x - (1+v12)y - (1+v13)z = d1

The parameter order is (u1,u2,v31,v32,v21,v23,v12,v13,d0,d3,d2,d1).
"""
from fractions import Fraction
from itertools import combinations
from functools import cmp_to_key

class Quadratic:
    """Number r+s*sqrt(2) with rational r,s."""
    __slots__=('r','s')
    def __init__(self,r=0,s=0):
        self.r,self.s=Fraction(r),Fraction(s)
    def __add__(self,other):
        other=number(other)
        return Quadratic(self.r+other.r,self.s+other.s)
    __radd__=__add__
    def __neg__(self):
        return Quadratic(-self.r,-self.s)
    def __sub__(self,other):
        return self+-number(other)
    def __rsub__(self,other):
        return number(other)+-self
    def __mul__(self,other):
        other=number(other)
        return Quadratic(self.r*other.r+2*self.s*other.s,self.r*other.s+self.s*other.r)
    __rmul__=__mul__
    def __truediv__(self,other):
        other=number(other)
        norm=other.r*other.r-2*other.s*other.s
        return self*Quadratic(other.r/norm,-other.s/norm)
    def __bool__(self):
        return bool(self.r or self.s)
    def __repr__(self):
        return f'({self.r},{self.s})'
    def sign(self):
        if not self.r:
            return (self.s>0)-(self.s<0)
        if not self.s or (self.r>0)==(self.s>0):
            return 1 if self.r>0 else -1
        return (1 if self.r>0 else -1)*(1 if self.r*self.r>2*self.s*self.s else -1)
    def __abs__(self):
        return self*self.sign()

def number(x):
    return x if isinstance(x,Quadratic) else Quadratic(x)

def add(x,y):
    return tuple(a+b for a,b in zip(x,y))

def sub(x,y):
    return tuple(a-b for a,b in zip(x,y))

def mul(x,s):
    return tuple(a*s for a in x)

def dot(x,y):
    return sum(a*b for a,b in zip(x,y))

def cross(x,y):
    return tuple(x[i-2]*y[i-1]-x[i-1]*y[i-2] for i in range(3))

def determinant2(x,y):
    return x[0]*y[1]-x[1]*y[0]

def perimeter(vertices):
    """Convex hull of exact projected vertices, returned counterclockwise."""
    def compare(p,q):
        return (p[0]-q[0]).sign() or (p[1]-q[1]).sign()
    vertices=sorted(vertices,key=cmp_to_key(compare))
    lower=[]
    upper=[]
    for sequence,chain in [(vertices,lower),(reversed(vertices),upper)]:
        for p in sequence:
            while len(chain)>=2 and determinant2(sub(chain[-1],chain[-2]),sub(p,chain[-1])).sign()<=0:
                chain.pop()
            chain.append(p)
    return lower[:-1]+upper[:-1]

a=Quadratic(0,1)
zero=number(0)
normals=[tuple(map(number,row)) for row in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),
         (1,1,1),(-1,-1,a),(-1,a,-1),(a,-1,-1)]]
offsets=[number(Fraction(1,2))]*6+[zero]*4
A_masks=(3,5,9,12,10,6)
B_masks=(7,11,13,8,4,2)
masks=A_masks+B_masks
dimension=12
normal_variations={
    6: [(0,(1,0,-1)),(1,(0,1,-1))],
    7: [(2,(-1,0,0)),(3,(0,-1,0))],
    8: [(4,(-1,0,0)),(5,(0,0,-1))],
    9: [(6,(0,-1,0)),(7,(0,0,-1))],
}

def linear_form(i,point):
    """Coefficients of n_i' * point - b_i' along a parameter direction."""
    out=[zero]*dimension
    if i>=6:
        for j,direction in normal_variations[i]:
            out[j]=dot(direction,point)
        out[i+2]=number(-1)
    return out

def all_vertices():
    vertices={}
    singular_without_opposite_faces=0
    for indices in combinations(range(10),3):
        i,j,k=indices
        columns=(cross(normals[j],normals[k]),cross(normals[k],normals[i]),cross(normals[i],normals[j]))
        det=dot(normals[i],columns[0])
        numer=add(add(mul(columns[0],offsets[i]),mul(columns[1],offsets[j])),mul(columns[2],offsets[k]))
        if not det:
            if not any(p<6 and p%2==0 and p+1 in indices for p in indices):
                assert any((abs(q)-1).sign()>0 for q in numer)
                singular_without_opposite_faces+=1
            continue
        assert (abs(det)-1).sign()>=0
        inverse_rows = [tuple(x/det for x in row) for row in zip(*columns)]
        for row in inverse_rows:
            assert (sum(map(abs,row))-5).sign()<0
        p=mul(numer,number(1)/det)
        sides=tuple((dot(n,p)-h).sign() for n,h in zip(normals,offsets))
        if 1 in sides[:6]:
            assert any((abs(q)-Fraction(7,10)).sign()>0 for q in p)
            continue

        # If the system of three vertex equations is A(t)x(t)=b(t),
        # these bound A(0)^(-1)A'(t) and x'(0), per unit parameter radius.
        for row in inverse_rows:
            matrix_row_bound=zero
            for axis in range(3):
                coeffs=[zero]*dimension
                for i, coefficient in zip(indices,row):
                    for parameter,direction in normal_variations.get(i,()):
                        coeffs[parameter]+=coefficient*direction[axis]
                matrix_row_bound+=sum(map(abs,coeffs))
            velocity=[zero]*dimension
            for i,coefficient in zip(indices,row):
                add_to(velocity,linear_form(i,p),-coefficient)
            assert (matrix_row_bound-4).sign()<=0
            assert (sum(map(abs,velocity))-Fraction(5,2)).sign()<=0
        vertices[tuple((q.r,q.s) for q in p)]=(p,sides)
    assert singular_without_opposite_faces==6
    return tuple(vertices.values())

def hull_on_plane(vertices,mask,i,axis):
    # A point on the boundary is admitted by either sign.
    pts={}
    for p,sides in vertices:
        if sides[i] or any(sides[j+6]*(1 if mask&(1<<j) else -1)<0 for j in range(4)):
            continue
        projected=tuple(p[j] for j in range(3) if j!=axis)
        pts[tuple((q.r,q.s) for q in projected)]=(p,sides)
    ordered=perimeter([tuple(Quadratic(*q) for q in key) for key in pts])
    return ordered,[pts[tuple((q.r,q.s) for q in y)] for y in ordered]

def add_to(acc,form,coefficient):
    for i,x in enumerate(form):
        if x:
            acc[i]+=x*coefficient

def add_product(matrix,left,right,coefficient):
    # Quadratic form in the parameters, evaluated as p^T matrix p.
    if not coefficient:
        return
    for i,x in enumerate(left):
        if x:
            ax=x*coefficient
            for j,y in enumerate(right):
                if y:
                    matrix[i][j]+=ax*y

def derivatives(vertices,mask):
    gradient=[zero]*dimension
    hessian=[[zero]*dimension for _ in range(dimension)]
    for i in range(6,10):
        # The chosen eliminated component is nonzero for all nearby cuts.
        axis=9-i
        if i==6:
            axis=2
        coordinate=normals[i][axis]
        y,vertices_on_facet=hull_on_plane(vertices,mask,i,axis)
        if len(y)<3:
            continue
        signs=1 if mask&(1<<(i-6)) else -1
        at_origin=[zero]*3
        at_origin[axis]=offsets[i]/coordinate
        W_origin=linear_form(i,at_origin)
        W=[linear_form(i,p) for p,_ in vertices_on_facet]
        integral=[zero]*dimension
        for k,p in enumerate(y):
            l=(k+1)%len(y)
            weight=determinant2(p,y[l])/6
            add_to(integral,add(add(W[k],W[l]),W_origin),weight)
        add_to(gradient,integral,number(signs)/abs(coordinate))

        # Differentiate W_i/|n_{i,axis}| inside the projected polygon.
        moving_component=[zero]*dimension
        for index,direction in normal_variations[i]:
            moving_component[index]=number(direction[axis])
        add_product(hessian,integral,moving_component,number(-2*signs)/(coordinate*abs(coordinate)))

        # The polygon boundary moves according to each other plane.
        for k,(p,sp) in enumerate(vertices_on_facet):
            l=(k+1)%len(y)
            q,sq=vertices_on_facet[l]
            shared=[j for j in range(10) if j!=i and not sp[j] and not sq[j]]
            assert len(shared)==1, (mask,i,k,shared)
            j=shared[0]
            tau=1 if j>=6 and mask&(1<<(j-6)) else -1
            tangential=tuple(normals[j][b]-normals[j][axis]*normals[i][b]/coordinate
                             for b in range(3) if b!=axis)
            free=0 if tangential[1] else 1
            factor=abs(y[l][free]-y[k][free])/abs(tangential[1-free])
            ratio=normals[j][axis]/coordinate
            p_motion=mul(sub(linear_form(j,p),mul(W[k],ratio)),tau)
            q_motion=mul(sub(linear_form(j,q),mul(W[l],ratio)),tau)
            weight=number(signs)*factor/(abs(coordinate)*6)
            add_product(hessian,W[k],p_motion,2*weight)
            add_product(hessian,W[k],q_motion,weight)
            add_product(hessian,W[l],p_motion,weight)
            add_product(hessian,W[l],q_motion,2*weight)
    return gradient,hessian

def inverse(matrix):
    d=len(matrix)
    rows=[list(row)+[number(i==j) for j in range(d)] for i,row in enumerate(matrix)]
    for col in range(d):
        pivot=next(i for i in range(col,d) if rows[i][col])
        rows[col],rows[pivot]=rows[pivot],rows[col]
        rows[col]=[q/rows[col][col] for q in rows[col]]
        for i in range(d):
            if i!=col:
                rows[i]=list(sub(rows[i],mul(rows[col],rows[i][col])))
    return [row[d:] for row in rows]

def matrix_norm_sum(H):
    return sum(abs((H[i][j]+H[j][i])/2) for i in range(dimension) for j in range(dimension))

def form(H,v):
    return dot(v,[dot(row,v) for row in H])

def centered_volume(vertices,mask):
    # Each face polygon subtends a pyramid of height 1/2 at the center.
    twice_area=zero
    for face in range(6):
        y,_=hull_on_plane(vertices,mask,face,face//2)
        twice_area+=sum(determinant2(p,y[(i+1)%len(y)]) for i,p in enumerate(y))
    return twice_area/12

if __name__=='__main__':
    vertices=all_vertices()
    m=(1+2*a)/48
    expected=[(3-2*a)/8]*16
    expected[0]=expected[15]=zero
    for mask in masks:
        expected[mask]=m
    volumes=[centered_volume(vertices,mask) for mask in range(16)]
    assert all(not (v-w) for v,w in zip(volumes,expected))
    assert not (sum(volumes)-1)

    J,H=map(tuple,zip(*(derivatives(vertices,mask) for mask in masks)))
    weights=[number(1)/(6*a)]*6+[(1-number(1)/a)/6]*6
    assert all((w-Fraction(1,21)).sign()>0 for w in weights)
    assert not (sum(weights)-1)
    for col in zip(*J):
        assert not dot(weights,col)
    F=[[dot(weights,[cell[i][j] for cell in H]) for j in range(dimension)] for i in range(dimension)]
    excluded={0,1,2}
    chosen_rows=(0,2,3,5,6,7,8,9,10)
    inv=inverse([[J[i][j] for j in range(dimension) if j not in excluded] for i in chosen_rows])
    for row in inv:
        assert (sum(map(abs,row))-70).sign()<0
    b=(8+9*a)/14
    c=number(2)/7
    r=list(map(number,[0,0,1,-1,-1,1,1,-1,0,0,0,0]))
    u=list(map(number,[1,-1,c,-c,-b+c,-b,b-c,b,0,0,0,0]))
    v=list(map(number,[0,1,-b,-b+c,b,b-c,c,-c,0,0,0,0]))
    # Two coordinate directions and the alternating direction span the kernel.
    for w in [r,u,v]:
        assert all(not dot(row,w) for row in J),w
    alternating_coefficient=3-101*a/48
    coordinate_coefficient=-number(137)/784+53*a/336
    assert not (form(F,r)/2-alternating_coefficient)
    assert not (form(F,u)/2-coordinate_coefficient)
    assert not (form(F,v)-form(F,u))
    assert not dot(r,[dot(add(x,y),u) for x,y in zip(F,zip(*F))])
    assert not dot(r,[dot(add(x,y),v) for x,y in zip(F,zip(*F))])
    assert not (dot(u,[dot(add(x,y),v) for x,y in zip(F,zip(*F))])+form(F,u))
    assert (alternating_coefficient-Fraction(1,50)).sign()>0
    assert (coordinate_coefficient/2-Fraction(1,50)).sign()>0
    assert (b+c-2).sign()<0
    assert all((matrix_norm_sum(M)-6).sign()<0 for M in H)
    assert (matrix_norm_sum(F)-Fraction(5,4)).sign()<0

    radius=Fraction(1,2_500_000)
    remainder=Fraction(26_000,6)
    assert 70*20*(3+remainder*radius)<4300
    assert 4300*radius<Fraction(1,500)
    assert ((Fraction(5,8)*4300*Fraction(1001,500)+remainder)*radius
            < Fraction(499,500)**2/250)
    print("Certified exact local derivative data and radius",radius)
