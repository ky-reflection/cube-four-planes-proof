"""Prove the affine minimum for the edge-separated tetrahedral family.

The geometry and exact cell formulas are in edge-separated-family-proof.md.
Write the weights as (x,y,z,1)/(1+x+y+z), with 0<x<=y<=z<=1. In the
scaled cube [-1,1]^3 the common threshold is 4*tau. Here s=tau/min(w)
is in [0,1]. All volumes are normalized by the cube volume.

For x<=1/20, the four cases in check_small_weights bound the largest
singleton below by S_d minus the tetrahedron volume. The latter is at
most (1+x+y+z)*min(1,x/z,x^2/(y*z))/12.

In the analytic neighborhood, use the centered proof from the note.
If x<=sqrt(2)-1, then S_d>=m. Since a<36*b*c*S_d and tau<=a, the
tetrahedron tau^3/(12*a*b*c*d) is less than 3*tau*S_d/d whenever
tau>0. Thus the affine singleton S_d*(1+tau/d)^3 minus the tetrahedron
is larger than S_d. If x>=sqrt(2)-1, the centered pair P_ad is at least
m. The larger of the two opposite affine pair cells is at least
P_ad+tau*(abs(ad-bc)-tau)/(4*K), where K is the positive product of
the four cross-pair weight sums. Here tau<abs(ad-bc), so any positive
tau again gives a strict bound. At tau=0 the centered equality
characterization applies.

All other parameter boxes use rational interval lower bounds on the
central tetrahedron, the largest-weight singleton, and the six pair
cells. A box is rejected only when a lower bound exceeds m by more
than the requested rational margin.
"""
from fractions import Fraction as F

from fourth_plane_certificate import exceeds_candidate
from seven_box_certificate import (
    interval_abs, interval_add, interval_multiply, interval_negate,
)


def point(value):
    return F(value), F(value)


def subtract(a, b):
    return interval_add(a, interval_negate(b))


def divide(a, b):
    assert b[0] > 0
    return interval_multiply(a, (1/b[1], 1/b[0]))


def power(interval, exponent):
    assert interval[0] >= 0
    return tuple(value**exponent for value in interval)


def scale(value, interval):
    return interval_multiply(point(value), interval)


def singleton_largest(x, y, z):
    return (3+x+y+z)/(12*(1+x)*(1+y)*(1+z))


def check_analytic_regions(margin):
    # All x <= 1/20. The four cases use z<=1/10; 1/10<=z<=1/4;
    # z>=1/4,y<=1/4; and z>=1/4,y>=1/4.
    small = F(1,20)
    for singleton_lower, central_upper in (
        (F(1,4)/F(11,10)**2, F(13,120)),
        (F(4,25), F(1,2)*(1+small+F(1,2))/12),
        (singleton_largest(small,F(1,4),F(1)), F(1,5)*(2+small+F(1,4))/12),
        (singleton_largest(small,F(1),F(1)), F(1,25)*(3+small)/12),
    ):
        assert exceeds_candidate(singleton_lower-central_upper-margin)

    largest_x, smallest_yz, largest_s = F(23,50), F(9,10), F(1,5)
    lower_m = F(79,1000)
    assert exceeds_candidate(F(2,25))
    assert not exceeds_candidate(lower_m)
    smallest_large_weight = smallest_yz/(3+largest_x)
    largest_small_weight = largest_x/(1+largest_x+2*smallest_yz)
    assert largest_small_weight < 36*smallest_large_weight**2*lower_m
    assert largest_s*largest_x*(3+largest_x) < smallest_yz**2-largest_x


def proves_pair(weights, total, threshold, margin):
    # These are the unnormalized weights with largest weight one, so
    # the formula's linear, quadratic and cubic terms use total=sum(w).
    a, b, c, d = weights
    A, B = interval_add(a,b), interval_multiply(a,b)
    C, E = interval_add(c,d), interval_multiply(c,d)
    denominator = interval_multiply(
        interval_multiply(interval_add(a,c), interval_add(a,d)),
        interval_multiply(interval_add(b,c), interval_add(b,d)),
    )
    difference = subtract(E,B)
    gap = interval_abs(difference)
    center = (1-gap[1]**2/denominator[0])/12
    quadratic = interval_multiply(total, power(threshold,2))
    gain = divide(interval_multiply(total, subtract(
        interval_multiply(threshold,gap), quadratic,
    )), scale(4,denominator))[0]
    if exceeds_candidate(center+gain-margin):
        return True
    for high_sum, high_product, low_sum, low_product, delta in (
        (C,E,A,B,difference), (A,B,C,E,interval_negate(difference)),
    ):
        cubic_factor = divide(interval_add(
            interval_add(power(high_sum,2), subtract(low_product,high_product)),
            interval_multiply(low_sum,high_sum),
        ), high_product)
        correction = interval_add(
            scale(3, subtract(interval_multiply(threshold,delta),quadratic)),
            interval_multiply(power(threshold,3),cubic_factor),
        )
        lower = divide(interval_multiply(total,correction),scale(12,denominator))[0]
        if exceeds_candidate(center+lower-margin):
            return True
    return False


def verify(max_nodes=100000, margin=F(1, 10**14)):
    check_analytic_regions(margin)
    start = ((F(1,20),F(1)),)*3 + ((F(0),F(1)),)
    pending = [(start,0)]
    visited = analytic = deepest = 0
    while pending:
        box, depth = pending.pop()
        visited += 1
        if visited > max_nodes or depth > 52:
            raise RuntimeError(f"unresolved affine box at depth {depth}: {box}")
        x, y, z, s = box
        x = x[0], min(x[1],y[1],z[1])
        z = max(x[0],y[0],z[0]), z[1]
        y = max(y[0],x[0]), min(y[1],z[1])
        box = x, y, z, s
        if any(lo > hi for lo,hi in box):
            continue
        deepest = max(deepest,depth)
        if x[1]<=F(23,50) and y[0]>=F(9,10) and z[0]>=F(9,10) and s[1]<=F(1,5):
            analytic += 1
            continue

        weights = x, y, z, point(1)
        total = interval_add(interval_add(x,y),interval_add(z,point(1)))
        threshold = interval_multiply(s,x)
        central_upper = s[1]**3*total[1]*min(x[1]**2/(y[0]*z[0]),x[1]/z[0],F(1))/12
        central_lower = s[0]**3*x[0]**2*total[0]/(12*y[1]*z[1])
        singleton_lower = singleton_largest(x[1],y[1],z[1])*(1+threshold[0])**3-central_upper
        if exceeds_candidate(central_lower-margin) or exceeds_candidate(singleton_lower-margin):
            continue
        if any(proves_pair(tuple(weights[i] for i in order),total,threshold,margin)
               for order in ((0,3,1,2),(0,2,1,3),(0,1,2,3))):
            continue

        split = max(range(4),key=lambda i:(box[i][1]-box[i][0])/(2 if i==3 else 1))
        lo, hi = box[split]
        middle = (lo+hi)/2
        for interval in ((lo,middle),(middle,hi)):
            child = list(box)
            child[split] = interval
            pending.append((tuple(child),depth+1))
    return visited, analytic, deepest


if __name__ == "__main__":
    visited, analytic, deepest = verify()
    print(f"certified gap 10^-14 outside the analytic region: visited={visited}, "
          f"analytic boxes={analytic}, depth={deepest}")
