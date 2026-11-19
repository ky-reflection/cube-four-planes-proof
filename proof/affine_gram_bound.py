"""Rational arithmetic for the affine Gram-matrix lower bound.

The geometric implications used here are in affine-gram-bound-proof.md.
The angular anchor and derivative are verified by wedge_certificate.cpp
and angular_moment_certificate.cpp, respectively.
"""
from fractions import Fraction
from math import factorial, isqrt

F = Fraction
GRID = 10**6
DELTA = F(283, 100000)
M = F(1, 14) + DELTA
ANCHORS = (
    (F(1, 3), F(143, 250)),
    (F(200**2 - 131**2, 200**2 + 131**2), F(596, 1000)),
    (F(100**2 - 61**2, 100**2 + 61**2), F(616, 1000)),
)
ANGULAR_SLOPE = F(23, 100)
ONE_THIRD = F(1, 3)
PI_LOWER = F(6283, 2000)
PI_UPPER = F(3927, 1250)
SQRT_TWO_UPPER = F(1414214, GRID)
ANCHOR_ANGLE_LOWER = F(1910633, GRID)


def round_up(value):
    scaled = value * GRID
    return F(-(-scaled.numerator // scaled.denominator), GRID)


def sqrt_down(value):
    scaled = value * GRID**2
    return F(isqrt(scaled.numerator // scaled.denominator), GRID)


def sqrt_up(value):
    lower = sqrt_down(value)
    return lower if lower**2 == value else lower + F(1, GRID)


def triangle_area_ratio(offset):
    return (1 + F(7, 2) * offset) / (2 * (1 - F(7, 8) * offset))


def negative_dot_bound(disagreement):
    """Upper bound on the negative dot product from a disagreement upper bound."""
    bounds = []
    for initial_dot, initial_mass in ANCHORS:
        if disagreement <= initial_mass:
            bounds.append(initial_dot)
            continue

        def excludes(grid_index):
            candidate = F(grid_index, GRID)
            midpoint = (candidate + initial_dot) / 2
            return (ANGULAR_SLOPE * (candidate - initial_dot))**2 > (
                (disagreement - initial_mass)**2 * (1 - midpoint**2)
            )

        lower = int(GRID * initial_dot) + 1
        upper = GRID
        if not excludes(upper):
            bounds.append(F(1))
            continue
        while lower < upper:
            middle = (lower + upper) // 2
            if excludes(middle):
                upper = middle
            else:
                lower = middle + 1
        bounds.append(F(lower, GRID))
    return min(bounds)

def weight_bounds(negative_dot):
    upper = negative_dot / (1 + negative_dot)
    lower = max(F(0), 1 - 3 * upper)
    if negative_dot < F(1, 2):
        ratio = sqrt_down((1 - 2 * negative_dot) / 3)
        lower = max(lower, ratio / (1 + ratio))
    return lower, upper


def mean_norm_squared_bound(negative_dot, positive_dot):
    lower_weight, upper_weight = weight_bounds(negative_dot)
    total = 4 * upper_weight - 1
    if total == 0:
        return F(0)
    cap = upper_weight - lower_weight
    assert cap > 0 and total <= 4 * cap
    full = total // cap
    sum_squares = full * cap**2 + (total - full * cap)**2
    return ((1 - positive_dot) * sum_squares + positive_dot * total**2) / (
        16 * upper_weight**2
    )


def weight_product_lower_bound(negative_dot):
    lower, upper = weight_bounds(negative_dot)
    total = 1 - upper
    assert 0 < lower <= upper and 3 * lower <= total <= 3 * upper
    return min(
        first * second * (total - first - second)
        for first in (lower, upper)
        for second in (lower, upper)
        if lower <= total - first - second <= upper
    )


def disagreement_error_constants(absolute_dot):
    assert 0 <= absolute_dot < 1
    mixed = sqrt_up(3 / (1 - absolute_dot**2))
    square = 2 + mixed * absolute_dot
    return square, mixed


def verify():
    tetrahedron = 28 * DELTA
    assert tetrahedron + 7 * DELTA < F(2, 7)
    for _ in range(80):
        offset = tetrahedron + 7 * DELTA
        tetrahedron = min(tetrahedron, round_up(triangle_area_ratio(offset) * offset))
    offset = tetrahedron + 7 * DELTA
    negative_dot = negative_dot_bound(8 * M + 2 * offset)
    positive_dot = negative_dot_bound(7 * M + 2 * offset)

    for _ in range(20):
        square, mixed = disagreement_error_constants(max(negative_dot, positive_dot))
        error = 2 * (square + mixed) * offset**2
        negative_dot = min(negative_dot, negative_dot_bound(8 * M + error))
        positive_dot = min(positive_dot, negative_dot_bound(7 * M + error))
        mean_squared = mean_norm_squared_bound(negative_dot, positive_dot)
        positive_dot = min(
            positive_dot,
            round_up(-ONE_THIRD + 5 * (negative_dot - ONE_THIRD) + 8 * mean_squared),
        )
        lower_weight, upper_weight = weight_bounds(negative_dot)
        absolute_dot = max(negative_dot, positive_dot)
        determinant_squared = 1 - 3 * absolute_dot**2 - 2 * absolute_dot**3
        if determinant_squared > 0 and lower_weight > 0:
            positive_offsets = 3 * tetrahedron + 21 * DELTA
            denominator = (
                6 * weight_product_lower_bound(negative_dot) * sqrt_down(determinant_squared)
            )
            tetrahedron = min(
                tetrahedron, round_up((upper_weight * positive_offsets)**3 / denominator)
            )
            offset = tetrahedron + 7 * DELTA

    assert positive_dot <= ONE_THIRD
    mean_squared = mean_norm_squared_bound(negative_dot, positive_dot)
    # The Taylor series through degree 18 is a lower bound for cos here.
    cosine_lower = sum(
        (-1)**index * ANCHOR_ANGLE_LOWER**(2 * index) / factorial(2 * index)
        for index in range(10)
    )
    assert cosine_lower > -ONE_THIRD and SQRT_TWO_UPPER**2 > 2
    spherical_fraction = (
        6 * ANCHOR_ANGLE_LOWER - 6 * SQRT_TWO_UPPER * mean_squared
    ) / PI_UPPER - 3
    cap_cosine = 1 - spherical_fraction / 3
    assert 0 < spherical_fraction < 1
    assert cap_cosine > F(71, 100) and 2 * F(71, 100)**2 > 1
    centered_mass = PI_LOWER / 4 * (1 / cap_cosine**2 - 1)

    square, mixed = disagreement_error_constants(max(negative_dot, positive_dot))
    total_offsets = 3 * tetrahedron + 28 * DELTA
    positive_offsets = 3 * tetrahedron + 21 * DELTA
    # After crediting the +3t in the affine 2+2 mass, the loss
    # upper bound is convex in t and largest at one endpoint.
    loss = max(
        (3 * square - mixed) * (7 * DELTA) * (28 * DELTA)
        + mixed * ((21 * DELTA)**2 + (7 * DELTA)**2),
        (3 * square - mixed) * offset * total_offsets
        + mixed * (positive_offsets**2 + (7 * DELTA)**2) - 3 * tetrahedron,
    )
    margin = centered_mass - loss - 6 * M
    assert margin > 0
    print(f"t <= {float(tetrahedron):.6f}; Gram dots in [{-float(negative_dot):.6f}, {float(positive_dot):.6f}]")
    print(f"centered 2+2 mass >= {float(centered_mass):.12f}")
    print(f"affine 2+2 mass >= {float(centered_mass - loss):.12f}")
    print(f"affine 2+2 mass <= {float(6 * M):.12f} under M <= {float(M):.12f}")
    print(f"contradiction margin >= {float(margin):.12f}")


if __name__ == "__main__":
    verify()
