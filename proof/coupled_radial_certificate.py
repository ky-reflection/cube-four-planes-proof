"""Rational radial improvement using each plane's three linked face depths.

The geometry is in coupled-radial-proof.md. The main function returns a
lower bound on the excess over the six axial disks at a specified radius.
"""

import argparse
from bisect import bisect_left, bisect_right
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, lcm

from four_line_radial_certificate import PI_LOWER


RADIAL_GRID = 2000
DEPTH_GRID = 2048
COMPONENT_GRID = 512
ANGLE_GRID = 4000
COS_SCALE = 10**12
HORIZONS = (900, 1000, 1100, 1200)


def squared_cosine_tables():
    """Negated lower and upper bounds for the squared positive cosine."""
    lower, upper = [], []
    denominators = tuple(factorial(2*k) for k in range(8))
    for i in range(6285):
        theta = F(i, ANGLE_GRID)
        top = sum((-1)**k * theta**(2*k) / denominators[k] for k in range(7))
        bottom = top - theta**14 / denominators[7]
        a, b = (COS_SCALE * max(F(0), v)**2 for v in (bottom, top))
        lower.append(-(a.numerator // a.denominator))
        upper.append((-b.numerator) // b.denominator)
    assert lower == sorted(lower) and upper == sorted(upper)
    assert upper[-1] == 0 and all(a >= b for a, b in zip(lower, upper))
    return lower, upper


def integrated_weights():
    kappa = F(8, 25)
    coefficients = (F(1), F(5,4)*kappa, F(35,24)*kappa**2,
                    F(105,64)*kappa**3, F(231,128)*kappa**4)
    antiderivative = tuple(c / (k+2) / RADIAL_GRID**(k+2)
                           for k, c in enumerate(coefficients))
    scale = lcm(*(v.denominator for v in antiderivative))
    integer_coefficients = tuple(int(scale*v) for v in antiderivative)
    sums = tuple(sum(c*i**(k+2) for k,c in enumerate(integer_coefficients))
                 for i in range(HORIZONS[-1]+1))
    return scale, sums


@lru_cache(maxsize=1)
def line_moments():
    lower, upper = squared_cosine_tables()
    weight_scale, cumulative = integrated_weights()
    # Radius squared of the annulus endpoints, relative to the disk radius.
    radial_bins = tuple((5*RADIAL_GRID-2*i, 5*RADIAL_GRID+3*i,
                         RADIAL_GRID-i, cumulative[i]-cumulative[i-1])
                        for i in range(1,HORIZONS[-1]+1))
    values = {h: [] for h in HORIZONS}
    for depth in range(2*DEPTH_GRID):
        moment = 0
        for i, (outer_den, outer_num, inner_num, weight) in enumerate(radial_bins, 1):
            if depth*outer_den < DEPTH_GRID*outer_num:
                outer_target = COS_SCALE*depth*outer_den // (DEPTH_GRID*outer_num)
                outer_angle = bisect_left(upper, -outer_target)
                if (depth+1)*RADIAL_GRID >= DEPTH_GRID*inner_num:
                    inner_angle = 0
                else:
                    inner_target = -(-COS_SCALE*(depth+1)*RADIAL_GRID // (DEPTH_GRID*inner_num))
                    inner_angle = bisect_right(lower, -inner_target)-1
                assert 0 <= inner_angle <= outer_angle < len(upper)
                moment += 2*weight*(outer_angle-inner_angle)
            if i in values:
                values[i].append(moment)
    return values, weight_scale, cumulative


def component_bins(values, radius_squared):
    """Enclose line moments as a function of one squared normal component."""
    bins = []
    numerator, denominator = radius_squared.numerator, radius_squared.denominator
    for i in range(COMPONENT_GRID):
        first = i*denominator*DEPTH_GRID // ((COMPONENT_GRID-i)*numerator)
        if i+1 == COMPONENT_GRID:
            last = len(values)-1
        else:
            last = (i+1)*denominator*DEPTH_GRID // ((COMPONENT_GRID-i-1)*numerator)
        bins.append(max(values[first:min(last+1,len(values))],default=0))
    return bins


def three_face_moment(bins):
    """Maximize the enclosing sum over bins meeting p+q+r=1."""
    best = 0
    for i in range(COMPONENT_GRID//3+1):
        for j in range(i,min(COMPONENT_GRID,(COMPONENT_GRID-i)//2+1)):
            first = max(j,COMPONENT_GRID-i-j-3)
            last = min(COMPONENT_GRID-1,COMPONENT_GRID-i-j)
            if first <= last:
                best = max(best,bins[i]+bins[j]+max(bins[first:last+1]))
    return best


def excess_lower_bound(radius_squared):
    """Certified excess for 8/25 <= r^2/(1+r^2) and r <= 4/5."""
    assert F(8,17) <= radius_squared <= F(16,25)
    values, weight_scale, cumulative = line_moments()
    best = F(0)
    for h in HORIZONS:
        outer_squared = F(5*RADIAL_GRID+3*h,5*RADIAL_GRID-2*h)
        if radius_squared*outer_squared > 1:
            continue
        moment = F(three_face_moment(component_bins(values[h],radius_squared)),
                   weight_scale*ANGLE_GRID)
        bracket = 2*PI_LOWER*F(cumulative[h],weight_scale)-F(4,3)*moment
        best = max(best,bracket)
    return F(3,16)*radius_squared**2/(1+radius_squared)*best


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("radius_squared",type=F,nargs="*",default=[F(9,16),F(5929,10000),F(16,25)])
    for radius_squared in parser.parse_args().radius_squared:
        bound = excess_lower_bound(radius_squared)
        print(f"r^2={float(radius_squared):.9f}: coupled radial excess >= {float(bound):.12f}",flush=True)
