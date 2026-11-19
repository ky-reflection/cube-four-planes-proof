"""Exact first-order separation near the candidate for column weights.

Parameters are (c0,c1,c2,H0,H1,H2,H3), with c3=1 and cube [-1,1]^3.
The volume comparison and remainder bounds are in column-weight-affine-proof.md.
"""

from fractions import Fraction as F
from itertools import combinations

from local_certificate import Quadratic as Q


def volume_gradients():
    c = (Q(-1,1), Q(1), Q(1), Q(1))
    total = sum(c)
    size = 7
    d_c = tuple(tuple(Q(int(i == j)) for j in range(size)) for i in range(4))
    d_h = tuple(tuple(Q(int(j == 3+i)) for j in range(size)) for i in range(4))
    faces = []
    for pair in combinations(range(4),2):
        mask = sum(1 << i for i in pair)
        signs = [1 if i in pair else -1 for i in range(4)]
        shift = pair[0] ^ pair[1]
        area = [(Q(), [Q()]*size) for _ in range(16)]
        d0 = [Q()] * size
        area0 = Q(1)
        for i in range(4):
            j = i ^ shift
            shares = []
            for k in range(4):
                if k in (i,j):
                    continue
                denominator = c[j] + c[k]
                value = c[j] / denominator
                derivative = tuple((a + signs[i]*b - value*(a+e))/denominator
                                   for a,b,e in zip(d_c[j],d_h[i],d_c[k]))
                shares.append((value,derivative))
            (a,da),(b,db) = shares
            triangle = a*b/2
            derivative = [(a*v+b*u)/2 for u,v in zip(da,db)]
            area[mask ^ (1<<i)] = triangle,derivative
            area0 -= triangle
            d0 = [u-v for u,v in zip(d0,derivative)]
        area[mask] = area0,d0
        height = tuple(tuple(sum((signs[i]-signs[j]*c[i]/c[j])*d_h[i][k]
                                 for i in range(4))/total for k in range(size)) for j in range(4))
        faces.append((area,height))

    values = [Q() for _ in range(16)]
    derivatives = [[Q()]*size for _ in range(16)]
    for cardinality in (3,2,1):
        for mask in range(1,15):
            if bin(mask).count('1') != cardinality:
                continue
            j = next(i for i in range(4) if not mask & (1<<i))
            larger = mask | (1<<j)
            value = Q()
            gradient = [Q()]*size
            for area,height in faces:
                (a,da),(b,db) = area[mask],area[larger]
                value += (a+b)/6
                for k in range(size):
                    gradient[k] += (da[k]+db[k]-(a+b)*height[j][k])/6
            values[mask] = value-values[larger]
            derivatives[mask] = tuple(a-b for a,b in zip(gradient,derivatives[larger]))
    return values,derivatives


def verify():
    values, gradients = volume_gradients()
    m = Q(1,2)/48
    active = tuple(i for i in range(1,15) if not (values[i]-m))
    assert active == tuple(range(2,14))
    assert all((value-m).sign() <= 0 for value in values)
    weights = (
        (142662,0,142662,0,0,272512,102310,131645,0,104104,0,104104),
        (0,254775,40782,127482,79168,0,0,68516,163851,0,265427,0),
        (0,21785,104746,75753,174559,188147,0,314989,120021,0,0,0),
        (157733,269491,30414,25089,21785,0,0,0,0,104746,390742,0),
        (104746,75753,0,21785,174559,188147,0,314989,0,0,120021,0),
        (27659,294580,160488,0,21785,0,0,0,56936,104746,333806,0),
        (0,189226,0,189226,0,157215,0,207334,0,128499,0,128499),
        (128499,0,128499,0,207334,0,157215,0,189226,0,189226,0),
        (28595,235702,0,0,235702,97631,0,0,235702,97631,0,69036),
        (0,0,89817,279558,0,0,165126,232067,0,0,195482,37951),
        (0,0,28595,235702,235702,97631,0,0,0,69036,235702,97631),
        (89817,279558,0,0,0,0,165126,232067,195482,37951,0,0),
        (0,0,0,0,0,0,37951,195482,279558,89817,232067,165126),
        (97631,235702,97631,235702,235702,28595,69036,0,0,0,0,0),
    )
    for axis in range(7):
        for sign in (-1,1):
            w = weights[2*axis+(sign+1)//2]
            assert all(q >= 0 for q in w) and sum(w) > 0
            direction = [sum(q*gradients[i][j] for i,q in zip(active,w)) for j in range(7)]
            lower = sign*direction[axis] - sum(abs(g) for j,g in enumerate(direction) if j != axis)
            assert (lower/sum(w)-F(1,80)).sign() > 0

    r = F(1,1600)
    geometric_radius = F(1,200)
    root = Q(0,1)
    assert (root-1-geometric_radius-F(2,5)).sign() > 0
    assert (3+root-1-3*geometric_radius-F(27,8)).sign() > 0
    assert (root-2*geometric_radius-F(7,5)).sign() > 0
    assert 1-geometric_radius > F(99,100)
    assert F(2,1)/(F(27,8)*F(7,5)) < F(1,2)
    assert 4*F(150,49)+2+4*F(1,2) < 17
    assert 3*(F(150,49)+1) < 17
    assert F(40,7)+17*geometric_radius < 6
    assert (F(1,2)-6*m).sign() > 0
    assert F(5,2)+F(2,1)/F(27,8) < F(13,4)
    assert 12/F(27,8)**2+6/(F(27,8)*F(2,5))+1/F(2,5)**2 < 12
    assert F(5,2)+F(100,99) < F(18,5)
    assert 6*F(18,5)/F(27,8)+F(5,2)**2+F(100,99)**2 < 15
    assert (17+F(12,2)+F(13,4)*6)/6 + 2*(F(15,2)+F(18,5)*6)/6 < 18
    assert 18*r < F(1,80)
    return r


if __name__ == '__main__':
    print('certified the independent-offset candidate neighborhood of radius', verify())
