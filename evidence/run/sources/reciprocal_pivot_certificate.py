"""Exact agreement and tetrahedron checks for reciprocal maximum scores."""

from fractions import Fraction as F
from itertools import combinations
from math import comb

from affine_pyramid_matching_certificate import CAP, G, M, round_up
from centered_bernstein_certificate import bernstein


AGREEMENT = F(3619,10000)


def verify_cubic(coefficients, lo, hi):
    """Check nonnegativity on the closed interval by rational Bernstein bisection."""
    stack = [(lo,hi,0)]
    while stack:
        left,right,depth = stack.pop()
        width = right-left
        local = {(j,): width**j*sum(c*comb(k,j)*left**(k-j)
                                  for k,c in enumerate(coefficients) if k >= j)
                 for j in range(4)}
        if min(bernstein(local,(3,)).values()) >= 0:
            continue
        assert depth < 20, (coefficients,left,right)
        midpoint = (left+right)/2
        stack.extend(((left,midpoint,depth+1),(midpoint,right,depth+1)))


def square(line):
    a,b = line
    return a*a,2*a*b,b*b


def verify_cap(b, radius, cap=CAP):
    """If both nonpivots are at most b and |e| <= cap, then |r| <= radius."""
    s_max = b/(1+b)
    if radius >= s_max:
        return
    switch = (b+radius)/(2+b)
    for lo,hi,line in (
        (radius,min(s_max,switch),(-radius,F(1))),
        (max(radius,switch),s_max,(b,-(1+b))),
    ):
        if lo > hi:
            continue
        c0,c1,c2 = square(line)
        # radius*(2s-s^2-t^2-radius^2/3) - 2s*(1-t^2)*cap
        polynomial = (-radius*(c0+radius**2/3),
                      2*radius-radius*c1-2*cap*(1-c0),
                      -radius*(1+c2)+2*cap*c1,
                      2*cap*c2)
        verify_cubic(polynomial,lo,hi)


def section_vertices(constraints, upper_s):
    """Return affine vertices (t(s),r(s)) and their feasible s intervals."""
    for (a,b,c,d),(e,f,g,h) in combinations(constraints,2):
        determinant = a*f-b*e
        if not determinant:
            continue
        t = ((f*c-b*g)/determinant,(f*d-b*h)/determinant)
        r = ((a*g-e*c)/determinant,(a*h-e*d)/determinant)
        lo,hi = F(0),upper_s
        for i,j,k,l in constraints:
            value = i*t[0]+j*r[0]-k
            slope = i*t[1]+j*r[1]-l
            if slope > 0:
                hi = min(hi,-value/slope)
            elif slope < 0:
                lo = max(lo,-value/slope)
            elif value > 0:
                break
        else:
            if lo <= hi:
                yield t,r,lo,hi


def verify_agreement(a,b,radius,upper_s=None):
    """Check the bound for averaged pair coefficients in the stated domain."""
    assert 0 < a <= b <= 1
    average = (a+b)/2
    a_s = a/(1+a)
    geometric_limit = (a_s+b/(1+b))/2
    upper_s = geometric_limit if upper_s is None else min(upper_s,geometric_limit)
    zero = F(0)
    # Rows encode i*t+j*r <= k+l*s.
    constraints = (
        (-1,0,zero,zero), (0,-1,zero,zero), (1,1,zero,F(1)),
        (0,1,radius,zero), (1,0,average,-(1+a)),
        (1,0,average+(b-a)*a_s/2,-(1+b)),
    )
    checked = 0
    for t,r,lo,hi in section_vertices(constraints,upper_s):
        t0,t1,t2 = square(t)
        r0,r1,r2 = square(r)
        # s*(AGREEMENT - (s-s^2/3+(1/s-1)*r^2+t^2/(3s))).
        polynomial = (-r0-t0/3, AGREEMENT-r1+r0-t1/3,
                      -1-r2+r1-t2/3, F(1,3)+r2)
        verify_cubic(polynomial,lo,hi)
        checked += 1
    assert checked
    return checked


def verify_reciprocal_matching():
    b = F(347,500)
    offset_upper = F(353,2000)
    negative_upper = F(46,625)
    minimum_s = F(303,1000)
    minimum_r = F(301,4000)
    verify_cap(b,offset_upper)
    verify_cap(b,negative_upper,G)
    branches_r = verify_agreement(b,b,minimum_r)
    branches_s = verify_agreement(b,b,offset_upper,minimum_s)
    assert negative_upper <= minimum_r

    # The two pair sums of s_i exceed 2*minimum_s, so each pair's
    # sum of |t_i| is at most this number.
    s_max = b/(1+b)
    t_sum = 2*b-2*(1+b)*minimum_s
    determinant_slope = s_max*(1+t_sum)
    d_low,d_high = F(101,250),F(41,50)
    assert determinant_slope*(1+b) < 1
    assert d_low < 2*minimum_s-determinant_slope*t_sum
    assert 2*s_max < d_high

    slope_lower = 1-s_max/2-offset_upper/(6*(1-s_max**2))
    max_negative_total = 2*max(F(0),offset_upper-2*minimum_r)
    volume_bounds = []
    for count in (4,3):
        # Each pair has positive offset sum. With a negative offset,
        # replace one r_i by zero when taking the upper bound.
        tetrahedron_upper = max(
            (offset_upper*(max(0,count-k)*d_low+min(k,count)*d_high))**3
            /(24*d_low**(4-k)*d_high**k)
            for k in range(5)
        )
        negative = F(0) if count == 4 else max_negative_total
        tetrahedron_lower = (4*minimum_r*slope_lower-(1-slope_lower)*negative-2*G)/3
        assert tetrahedron_upper < tetrahedron_lower
        volume_bounds.append((tetrahedron_upper,tetrahedron_lower))
    return branches_r+branches_s,volume_bounds


def main():
    assert 48*M > 1 and (48*M-1)**2 > 8
    assert AGREEMENT < 1-8*M
    cases = (
        (F(8,15),F(8,15),F(8501,50000),F(8501,50000)),
        (F(23,100),F(1),F(157,1000),F(47,250)),
        (F(1,4),F(19,20),F(1571,10000),F(371,2000)),
        (F(3,10),F(17,20),F(319,2000),F(911,5000)),
    )
    for a,b,first,second in cases:
        verify_cap(a,first)
        verify_cap(b,second)
        checked = verify_agreement(a,b,(first+second)/2)
        print(f"reciprocal row limits ({a}, {b}): {checked} vertex branches verified")
    checked,bounds = verify_reciprocal_matching()
    print(f"two reciprocal pairs, all nonpivots <= 347/500: {checked} vertex branches verified")
    for case,(upper,lower) in zip(("all offsets nonnegative","some offset negative"),bounds):
        print(f"{case}: tetrahedron volume <= {float(round_up(upper)):.10f}, "
              f"requires > {float(-round_up(-lower)):.10f}")


if __name__ == "__main__":
    main()
