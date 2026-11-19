"""Verify the affine Gram lower bound using exact rational arithmetic.

See affine-gram-bm-proof.md for the geometry and required C++ anchor checks.
Run with Python 3. No third-party Python packages are required.
"""

import argparse
from fractions import Fraction
from itertools import combinations
from math import factorial, isqrt

from four_line_radial_certificate import excess_lower_bound, verify as verify_radial

F = Fraction
GRID = 10**6
DEFAULT_DELTA = F(428, 100000)
LEGACY_PAIRS = ((200, 131, 404), (100, 61, 384))
SHARP_PAIRS = (
    (14, 9, 398), (25, 16, 397), (5, 4, 469), (20, 16, 468), (17, 14, 478),
    (16, 13, 475), (100, 81, 473), (100, 82, 477), (100, 63, 392), (43, 27, 391),
    (19, 12, 393), (8, 5, 390), (30, 19, 394), (76, 49, 399),
)
ANCHORS = ((F(1, 3), F(143, 250)),) + tuple(
    (F(a * a - b * b, a * a + b * b), 1 - F(upper, 1000))
    for a, b, upper in LEGACY_PAIRS + SHARP_PAIRS
)
LOWER_ANGULAR_SLOPE = F(23, 100)
UPPER_ANGULAR_SLOPE = F(51, 100)
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


def triangle_area_ratio(offset):
    eta = F(7, 4) * offset
    assert 0 <= eta <= F(1, 2)
    u = F(2, 3) * eta * (1 - eta)
    return (1 + 2 * u) / (2 - eta + u)


def negative_dot_bound(disagreement):
    """Upper bound on the negative dot product from a disagreement upper bound."""
    bounds = []
    for initial_dot, initial_mass in ANCHORS:
        if disagreement <= initial_mass:
            backwards = (initial_mass - disagreement) * sqrt_down(1 - initial_dot**2) / UPPER_ANGULAR_SLOPE
            bounds.append(max(F(0), round_up(initial_dot - backwards)))
            continue

        def excludes(grid_index):
            candidate = F(grid_index, GRID)
            midpoint = (candidate + initial_dot) / 2
            return (LOWER_ANGULAR_SLOPE * (candidate - initial_dot))**2 > (
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


def mean_norm_squared_bound(negative_dot, positive_dot, alpha_floor = F(0)):
    lower_weight, upper_weight = weight_bounds(negative_dot)
    original_total = 4 * upper_weight - 1
    total = original_total - 4 * alpha_floor
    cap = upper_weight - lower_weight - alpha_floor
    assert total >= 0 and cap >= 0 and total <= 4 * cap
    if total == 0:
        return F(0)
    full = total // cap
    sum_squares = full * cap**2 + (total - full * cap)**2
    if lower_weight >= upper_weight / 2:
        sum_squares = min(sum_squares, original_total**2 / 2 - 2 * alpha_floor * original_total + 4 * alpha_floor**2)
    return ((1 - positive_dot) * sum_squares + positive_dot * total**2) / (16 * (upper_weight - alpha_floor)**2)


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


def pair_error(sum_deficit, non_disagreement):
    assert 0 <= sum_deficit <= 3 * non_disagreement
    return sum_deficit**2 / (3 * non_disagreement) - sum_deficit**3 / (27 * non_disagreement**2)


def centered_disagreement_upper(affine, H, G, reverse):
    bound = affine + 2 * H
    if bound >= 1 or 2 * H > 3 * (1 - bound):
        return min(bound, F(1))
    for _ in range(20):
        q = 1 - bound
        if reverse:
            err = max(pair_error(H + G, q), 2 * pair_error(H, q))
        else:
            err = pair_error(2 * H, q)
        nxt = min(bound, round_up(affine + err))
        if nxt == bound:
            break
        bound = nxt
    return bound


ALL_EDGES = tuple(combinations(range(4), 2))
STAR = ((0, 1), (0, 2), (0, 3))
FIVE_EDGES = ALL_EDGES[1:]


def graph_pair_bound(maximum, delta, t_max, disagree_upper, edges, tetra_credit, use_cut):
    G = 7 * delta
    uniform_error = pair_error(2 * (t_max + G), 1 - disagree_upper)
    same_q = 1 - 8 * maximum - uniform_error
    candidates = []
    for t in (F(0), t_max):
        H = t + G
        deficit = 14 * delta + t
        opposite_q = same_q + deficit
        plus = pair_error(2 * H, same_q)
        minus = pair_error(2 * G, same_q)
        opposite = pair_error(H, opposite_q) + pair_error(G, opposite_q)
        if use_cut:
            opposite -= deficit
        for mask in range(1, 15):
            total = 0
            for i, j in edges:
                first = bool(mask & (1<<i))
                second = bool(mask & (1<<j))
                total += opposite if first != second else (plus if first else minus)
            candidates.append(total - tetra_credit * t)
    bound = max(candidates)
    return 8 * len(edges) * maximum + bound if use_cut else bound


def graph_deficit_lower(edge_count, centered_sum_upper, negative_dot, positive_dot):
    bounds = [F(0)]
    for c0, D0 in ANCHORS:
        max_abs = max(c0, positive_dot)
        slope = UPPER_ANGULAR_SLOPE / sqrt_down(1 - max_abs**2)
        gain = (edge_count * D0 - centered_sum_upper) / slope - edge_count * max(F(0), c0 - negative_dot)
        bounds.append(gain)
    return max(bounds)


def positive_dot_and_weight_floor(maximum, delta, t_max, disagree_upper, negative_dot, positive_dot):
    if (positive_dot >= 1 or disagree_upper >= 1 or 2 * (t_max + 7 * delta) > 3 * (1 - disagree_upper)):
        return positive_dot, F(0)
    uniform_error = pair_error(2 * (t_max + 7 * delta), 1 - disagree_upper)
    if 2 * (t_max + 7 * delta) > 3 * (1 - 8 * maximum - uniform_error):
        return positive_dot, F(0)
    star_centered = graph_pair_bound(maximum, delta, t_max, disagree_upper, STAR, 0, True)
    star_lower = graph_deficit_lower(3, star_centered, negative_dot, positive_dot)
    weight_lower, _ = weight_bounds(negative_dot)
    alpha_floor = weight_lower * star_lower / (1 + negative_dot)
    mean_squared = mean_norm_squared_bound(negative_dot, positive_dot, alpha_floor)
    five_centered = graph_pair_bound(maximum, delta, t_max, disagree_upper, FIVE_EDGES, 0, True)
    five_lower = graph_deficit_lower(5, five_centered, negative_dot, positive_dot)
    positive_dot = min(positive_dot, round_up(-2 + 5 * negative_dot + 8 * mean_squared - five_lower))
    return positive_dot, alpha_floor


def angular_quadratic(upper_dot):
    theta_upper = F(1910634, GRID)
    tangent_lower = F(1060659, GRID)
    q_max = F(19, 100)
    assert cosine_taylor(theta_upper, 8) < -ONE_THIRD and tangent_lower**2 < F(9, 8)
    assert (2 * q_max)**2 * F(8, 9)**3 < ONE_THIRD**2
    if upper_dot <= -ONE_THIRD:
        return q_max
    lower, upper = 0, 2 * GRID
    while upper - lower > 1:
        midpoint = (lower + upper)//2
        if cosine_taylor(F(midpoint, GRID), 9) > upper_dot:
            lower = midpoint
        else:
            upper = midpoint
    x = upper_dot + ONE_THIRD
    bound = (F(lower, GRID) - theta_upper + tangent_lower * x) / (x * x)
    scaled = bound * GRID
    return min(q_max, max(F(0), F(scaled.numerator//scaled.denominator, GRID)))


def cosine_taylor(x, terms):
    return sum((-1)**i * x**(2 * i) / factorial(2 * i) for i in range(terms + 1))


def verify(delta = DEFAULT_DELTA):
    M = F(1, 14) + delta
    tetrahedron = 28 * delta
    assert tetrahedron + 7 * delta < F(2, 7)
    for _ in range(80):
        offset = tetrahedron + 7 * delta
        tetrahedron = min(tetrahedron, round_up(triangle_area_ratio(offset) * offset))
    offset = tetrahedron + 7 * delta
    negative_dot = F(1)
    positive_dot = F(1)
    disagree_upper = F(1)
    for _ in range(20):
        disagree_upper = min(disagree_upper, centered_disagreement_upper(8 * M, offset, 7 * delta, False))
        agree_upper = centered_disagreement_upper(7 * M, offset, 7 * delta, True)
        negative_dot = min(negative_dot, negative_dot_bound(disagree_upper))
        positive_dot = min(positive_dot, negative_dot_bound(agree_upper))
        mean_squared = mean_norm_squared_bound(negative_dot, positive_dot)
        positive_dot = min(
            positive_dot,
            round_up(-ONE_THIRD + 5 * (negative_dot - ONE_THIRD) + 8 * mean_squared),
        )
        positive_dot, alpha_floor = positive_dot_and_weight_floor(
            M, delta, tetrahedron, disagree_upper, negative_dot, positive_dot
        )
        lower_weight, upper_weight = weight_bounds(negative_dot)
        absolute_dot = max(negative_dot, positive_dot)
        determinant_squared = 1 - 3 * absolute_dot**2 - 2 * absolute_dot**3
        if determinant_squared > 0 and lower_weight > 0:
            positive_offsets = 3 * tetrahedron + 21 * delta
            denominator = (
                6 * weight_product_lower_bound(negative_dot) * sqrt_down(determinant_squared)
            )
            tetrahedron = min(
                tetrahedron, round_up((upper_weight * positive_offsets)**3 / denominator)
            )
            offset = tetrahedron + 7 * delta

    assert negative_dot >= ONE_THIRD and positive_dot >= -negative_dot
    assert positive_dot <= ONE_THIRD
    mean_squared = mean_norm_squared_bound(negative_dot, positive_dot, alpha_floor)
    # The Taylor series through degree 18 is a lower bound for cos here.
    cosine_lower = cosine_taylor(ANCHOR_ANGLE_LOWER, 9)
    assert cosine_lower > -ONE_THIRD and SQRT_TWO_UPPER**2 > 2
    q = angular_quadratic(positive_dot)
    spherical_fraction = (
        6 * ANCHOR_ANGLE_LOWER - (6 * SQRT_TWO_UPPER - F(16, 3) * q) * mean_squared
    ) / PI_UPPER - 3
    cap_cosine = 1 - spherical_fraction / 3
    assert 0 < spherical_fraction < 1
    assert cap_cosine > F(71, 100) and 2 * F(71, 100)**2 > 1
    radius_squared = 1 / cap_cosine**2 - 1
    centered_mass = PI_LOWER / 4 * radius_squared
    if F(9, 16) <= radius_squared <= F(16, 25):
        verify_radial(verbose = False)
        centered_mass += excess_lower_bound(radius_squared)

    loss = graph_pair_bound(M, delta, tetrahedron, disagree_upper, ALL_EDGES, 3, False)
    margin = centered_mass - loss - 6 * M
    assert margin > 0, f"The comparison did not close for delta = {delta}."
    print(f"t <= {float(tetrahedron):.6f}")
    print(f"Gram dots in [{-float(negative_dot):.6f}, {float(positive_dot):.6f}]")
    print(f"centered 2+2 mass >= {float(centered_mass):.12f}")
    print(f"affine 2+2 mass >= {float(centered_mass - loss):.12f}")
    print(f"affine 2+2 mass <= {float(6 * M):.12f} under M <= {float(M):.12f}")
    print(f"contradiction margin >= {float(margin):.12f}")
    return margin


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("delta", nargs = "?", type = F, default = DEFAULT_DELTA)
    verify(parser.parse_args().delta)
