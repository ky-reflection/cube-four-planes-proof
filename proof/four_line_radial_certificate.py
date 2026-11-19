"""Exact quadrature for the four-line improvement to cube radial rearrangement.

See four-line-radial-proof.md for the geometric reduction.
"""
from bisect import bisect_left, bisect_right
from fractions import Fraction
from math import factorial

F = Fraction
GRID = 1024
ANGLE_GRID = 2000
COS_SCALE = 10**12
PI_LOWER = F(6283, 2000)
BRACKET_LOWER = F(619, 10000)
EXCESS_LOWER = F(47, 10000)


def cosine_tables():
    """Negated bounds on COS_SCALE*cos(i/ANGLE_GRID), in increasing order."""
    lower = []
    upper = []
    denominators = [factorial(2 * k) for k in range(8)]
    for i in range(3143):
        angle = F(i, ANGLE_GRID)
        cosine_upper = sum(
            (-1)**k * angle**(2 * k) / denominators[k] for k in range(7)
        )
        cosine_lower = cosine_upper - angle**14 / denominators[7]
        rounded_lower = cosine_lower * COS_SCALE
        rounded_upper = cosine_upper * COS_SCALE
        lower.append(-(rounded_lower.numerator // rounded_lower.denominator))
        upper.append((-rounded_upper.numerator) // rounded_upper.denominator)
    assert all(a >= b for a, b in zip(lower, upper))
    assert lower == sorted(lower) and upper == sorted(upper)
    return lower, upper


def weight_numerator(i):
    """80*GRID**5 times the integral of y(1-y)(2-y)(1+3y/4) over [0,i/GRID]."""
    return 80 * i**2 * GRID**3 - 40 * i**3 * GRID**2 - 25 * i**4 * GRID + 12 * i**5


def excess_lower_bound(radius_squared):
    """Radius-dependent bound once verify() has certified BRACKET_LOWER."""
    assert F(9, 16) <= radius_squared <= F(16, 25)
    return F(3, 8) * radius_squared**2 / (1 + radius_squared) * BRACKET_LOWER


def verify(verbose=True):
    lower_cosines, upper_cosines = cosine_tables()
    horizon = GRID // 4
    weights = [
        weight_numerator(i) - weight_numerator(i - 1)
        for i in range(1, horizon + 1)
    ]

    max_weighted_angle = 0
    for distance in range(GRID + horizon):
        weighted_angle = 0
        for radius, weight in enumerate(weights, 1):
            if distance >= GRID + radius:
                continue
            # For this rectangle in (distance, radius), bound the angle
            # on the outer circle at the smallest distance.
            outer = bisect_left(
                upper_cosines, -(COS_SCALE * distance // (GRID + radius))
            )
            assert outer < len(upper_cosines)
            # On the inner circle use the largest distance and round the
            # subtracted angle down.
            if distance + 1 >= GRID - radius:
                inner = 0
            else:
                target = -(-COS_SCALE * (distance + 1) // (GRID - radius))
                inner = bisect_right(lower_cosines, -target) - 1
                assert inner >= 0
            assert outer >= inner
            weighted_angle += 2 * weight * (outer - inner)
        max_weighted_angle = max(max_weighted_angle, weighted_angle)

    weight_integral = F(weight_numerator(horizon), 80 * GRID**5)
    angle_moment = F(max_weighted_angle, 80 * GRID**5 * ANGLE_GRID)
    bracket = 2 * PI_LOWER * weight_integral - 4 * angle_moment
    assert bracket > BRACKET_LOWER
    minimum_radius = F(3, 4)
    excess = F(3, 8) * minimum_radius**4 / (1 + minimum_radius**2) * bracket
    assert excess_lower_bound(minimum_radius**2) > EXCESS_LOWER
    if verbose:
        print(f"one-line angular moment <= {float(angle_moment):.12f}")
        print(f"radial bracket >= {float(bracket):.12f} > {float(BRACKET_LOWER)}")
        print(f"excess cube volume >= {float(excess):.12f} > {float(EXCESS_LOWER)}")
    return EXCESS_LOWER


if __name__ == "__main__":
    verify()
