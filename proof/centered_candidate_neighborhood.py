"""Exact constants for the centered candidate neighborhood.

The proof is in centered-bernstein-completion.md. No numerical linear
program or approximate derivative is needed for this neighborhood.
"""

from fractions import Fraction as F
from itertools import combinations

from centered_bernstein_certificate import plus, scale, times
from local_certificate import Quadratic as Q


def less(a, b):
    return (Q() + b - a).sign() > 0


def coefficient(poly, degree):
    return poly.get((degree,), Q())


def verify_row(c, weights, x, w, symmetric):
    directions = ((F(1), F(-1,2), F(-1,2)), (F(0), F(1), F(-1)),
                  (F(1), F(1,2), F(-3,2)))
    for h in directions:
        a, b, d = ({(0,): u, (1,): u*v} for u,v in zip(c,h))
        ab, ad, bd = plus(a,b), plus(a,d), plus(b,d)
        denominator = times(times(ab,ad),bd)
        terms = (times(times(a,a),bd), times(times(b,b),ad), times(times(d,d),ab))
        numerator = {}
        for term, weight in zip(terms,weights):
            numerator = plus(numerator,scale(term,weight))
        base = coefficient(numerator,0) / coefficient(denominator,0)
        residue = plus(numerator,scale(denominator,-base))
        assert not coefficient(residue,0) and not coefficient(residue,1)
        quadratic = coefficient(residue,2) / coefficient(denominator,0)
        if symmetric:
            expected = w*sum((u-v)**2 for u,v in combinations(h,2))/32
        else:
            first = h[0]-(h[1]+h[2])/2
            second = (h[1]-h[2])/2
            expected = x*x*(first**2+(1+2*x)*second**2)/((1+x)*(1+x)*(1+x)*(1+x))
        assert not (quadratic-expected)


def verify():
    root = Q(0,1)
    x = root-1
    w = 2*x/(1+x)
    b = Q(1)/(3+x)
    radius, e, norm = F(1,200), F(1,80), F(1,40)
    raw = Q(radius)/x
    normalization = 1-2*raw/3
    assert less(raw,e*normalization)
    assert less(radius**2*(Q(3)/(x*x)+5),norm**2*normalization*normalization)

    verify_row((Q(1),)*3,(w,)*3,x,w,True)
    verify_row((x,Q(1),Q(1)),(Q(1),w,w),x,w,False)
    fourth = (1+x)*(1+x)*(1+x)*(1+x)
    assert less(3*w/16,F(129,1000))
    assert less(3*x*x/fourth,F(129,1000))
    high_triangle = Q(1)/(2*(1+x))
    low_triangle = x*x/((1+x)*(1+x))
    for triangle, weighted in ((Q(F(1,4)),w/4),(high_triangle,w*high_triangle),(low_triangle,low_triangle)):
        assert less(triangle,F(354,1000)) and less(weighted,F(21,100))
    assert F(21,100)+F(354,1000)+F(1,12) < F(65,100)

    assert less(F(28,1000),(1+w)*x*x*b/((1+x)*(1+x)*(1+x)))
    assert less(F(28,1000),w*b/4)
    penalty = F(28,1000)/(1+e)**4
    assert F(15,4)*penalty > F(99,1000)
    assert (F(129,2000)+F(65,100)*e)/(1-e)**3 + 30*penalty*e**2 < F(76,1000)

    assert less(root/(4*(1-e)),F(9,25))
    assert less(Q(1)/(1+x),F(3,4))
    assert 2*F(3,4)*F(9,25)+F(9,25)**2 < F(7,10)
    assert F(1,4)/(1-e)**2 < F(26,100)
    assert 3*(F(7,10)+2*F(26,100))**2 < F(11,5)**2

    a = (1+2*x)/(2*(1+x)*(1+x))
    d = (1-x)/(4*(1+x)*(1+x))
    g = F(3,8)+(5*a-d)/3
    k = (-8*a+16*d)/3
    diff = F(3,8)-a
    diff_sq = diff*diff+d*d
    assert less(-k/2,g) and less(0,-k)
    assert less(F(8,3)*diff_sq,g*g/36)
    assert less(F(8,3),4*k*k)
    assert less(diff_sq+F(3,8)*g*g,F(7,10)**2)
    assert 4*F(11,5) < 9

    contraction = F(6,5)*(F(1,6)+9*norm)**2
    assert F(1,1)/(1-contraction) < F(5,4)
    assert 96/(1-contraction) < 120
    assert F(76,1000)*F(5,4) < F(99,1000)
    assert F(7,10)*norm+F(11,5)*norm**2 < F(19,1000)
    assert less(120*F(76,1000)*F(19,1000),1-w)
    offset_fraction = F(3,25)
    assert less((1+2*root)/12, F(1,3)-offset_fraction**2+F(2,3)*offset_fraction**3)
    return radius


if __name__ == "__main__":
    print("certified the centered candidate neighborhood of radius", verify())
