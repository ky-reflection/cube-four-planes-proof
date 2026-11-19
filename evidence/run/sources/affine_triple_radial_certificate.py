"""Exact certificate for C > 383/5000. See affine-triple-radial-proof.md."""

from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

from affine_frame_height_certificate import verify as height_bounds
from affine_gram_bm_certificate import PI_LOWER, PI_UPPER, cosine_taylor
from affine_weight_graph_certificate import EDGES, intercept
from four_line_radial_certificate import excess_lower_bound, verify as check_four_lines


M = F(383, 5000)
DELTA = M - F(1, 14)
PAIR_INTERCEPT, PAIR_SLOPE = F(".44470"), F(".37")
FACE_SCALE, TAIL = F("1.52"), F(".0095")
CAP_COSINE = F(".79")
# Each triple omits one of the four labels.
TRIANGLES = tuple(tuple(e for e, pair in enumerate(EDGES) if i not in pair)
                  for i in range(4))
STAR_LINE = (F("-.19"), F("1.56843"), F("1.018554"))
TRIANGLE_LINE = (F("-.425"), F("1.54024"), F("1.104736"))


def corner_tail_upper():
    """Upper Riemann sum on one face quadrant, including the area normalization."""
    grid, scale = 512, 10**9
    numerator = FACE_SCALE.numerator**3 * grid**3 * scale
    total = 0
    for i in range(1, grid + 1):
        for j in range(1, i + 1):
            q = grid**2 + i*i + j*j
            root_upper = isqrt(q * scale**2) + 1
            denominator = FACE_SCALE.denominator**3 * q * root_upper
            # Round 1 - s^3 / (1+x^2+y^2)^(3/2) upwards to 1/scale.
            upper = max(0, scale - numerator * scale // denominator)
            total += upper if i == j else 2 * upper
    return F(total, grid**2 * scale)


def check_angle_line(line, c, u):
    center, intercept_value, slope = line
    tolerance = F(1, 10**6)
    assert -c < center < 0 < u and max(center+c, u-center) < 1
    assert (slope-tolerance)**2 * (1-center**2) < 1
    assert (slope+tolerance)**2 * (1-center**2) > 1
    # The tangent at center supports arccos on the negative half. On the
    # positive half the difference is concave, so only its far endpoint remains.
    for x, allowance in ((center, tolerance), (u, 0)):
        angle = intercept_value - slope*x + allowance
        assert 0 < angle < PI_LOWER and cosine_taylor(angle, 9) > x


def deficit_states(t_max, upper_dot, upper_disagreement):
    states = []
    _, angle_intercept, angle_slope = TRIANGLE_LINE
    for t in (0, t_max):
        deficit = 14*DELTA + t
        for mask in range(1, 15):
            signs = [bool(mask & (1 << i)) for i in range(4)]
            caps = [deficit*sign - 7*DELTA for sign in signs]
            centered_upper, errors = [], []
            for i, j in EDGES:
                affine = 8*M - deficit*(signs[i] != signs[j])
                error = (caps[i]+caps[j])**2 / (3*(1-upper_disagreement))
                centered_upper.append(affine+error)
                errors.append(error)
            lower = tuple((PAIR_INTERCEPT-d)/PAIR_SLOPE for d in centered_upper)
            assert max(lower) <= -F(".18") < upper_dot
            triple_lower = []
            for triple in TRIANGLES:
                mono = 1 - sum(centered_upper[e] for e in triple)/2
                triple_lower.append((12*(mono-TAIL)/FACE_SCALE**3 - 2*PI_UPPER
                                     + 3*angle_intercept)/angle_slope)
            loss = sum(errors) - 3*t - deficit*(sum(signs) == 2)
            states.append((lower, triple_lower, loss))
    return states


# Rational multipliers in units of 1/100000, found by linear programming.
# Four coefficients for G*lambda=0, then four nonnegative triangle multipliers.
WITNESSES = (
    (218161, -98761, 1283896, 5000000, 232676, 0, 0, 0),
    (104143, 47481, 2081351, 5000000, 506657, 0, 0, 0),
    (78679, 79270, 78009, 93468, 4063, 0, 0, 0),
    (50169, 123276, 123276, 155920, 27186, 0, 0, 0),
    (32293, 151658, 157741, 197610, 44669, 0, 0, 0),
    (-2608763, 4344957, 4344957, 5000000, 2303632, 0, 0, 0),
    (-2531297, 4240669, 4221961, 5000000, 2242405, 0, 0, 0),
    (-1726304, 3482782, 5000000, 5000000, 1982486, 0, 0, 0),
    (74280, 84565, 74249, 105145, 0, 0, 0, 0),
    (-273755, 500680, 3257655, 5000000, 841829, 0, 0, 0),
    (-164534, 372532, 1875562, 5000000, 484095, 0, 0, 0),
    (74280, 84565, 67898, 127079, 0, 0, 0, 0),
    (-424172, 627217, 3244654, 5000000, 859707, 0, 0, 0),
    (77028, 85243, 64001, 102565, 466, 0, 0, 0),
    (73045, 85800, 114947, 136861, 10639, 9003, 0, 0),
    (71761, 90938, 95721, 115984, 8903, 2643, 0, 0),
    (71296, 87548, 85751, 106831, 3951, 0, 0, 0),
    (68340, 99307, 1165020, 5000000, 253819, 0, 0, 0),
    (69121, 93792, 121555, 223349, 16464, 0, 0, 0),
    (82425, 76250, 78828, 102365, 3100, 0, 0, 0),
    (74652, 87739, 162576, 174381, 28038, 0, 0, 0),
    (46825, 128869, 143946, 143440, 34890, 0, 0, 0),
    (69579, 97906, 97147, 95631, 14397, 0, 0, 0),
    (76900, 87831, 96882, 96882, 9857, 0, 0, 0),
    (73228, 87860, 127823, 127823, 17040, 0, 0, 0),
    (75491, 87952, 93045, 91291, 8845, 0, 0, 0),
    (71017, 91742, 112851, 121760, 13952, 0, 0, 0),
    (75971, 82874, 79433, 87651, 378, 0, 0, 0),
    (77586, 81925, 77250, 68486, 705, 0, 0, 0),
    (80696, 74806, 74806, 71381, 1848, 0, 0, 0),
    (67815, 101299, 101299, 101299, 16672, 0, 0, 0),
)


def cover_weights(c, upper_dot, states, coefficients, threshold):
    cap = c/(1+c)
    assert F(1, 4) < cap < F(1, 3)
    root = ((F(1, 4),)*4,) + tuple(
        ((1-(4-i)*cap)/i,)*i + (cap,)*(4-i) for i in (1, 2, 3))
    witnesses = [tuple(F(x, 100000) for x in row) for row in WITNESSES]
    assert all(min(w[4:]) >= 0 for w in witnesses)

    @lru_cache(maxsize=None)
    def margin(weights, witness):
        k, y = witness[:4], witness[4:]
        residual = [a-k[i]*weights[j]-k[j]*weights[i]
                    + sum(y[h] for h in range(4) if h not in (i, j))
                    for a, (i, j) in zip(coefficients, EDGES)]
        constant = threshold + sum(a*b for a, b in zip(k, weights))
        return min(constant-loss + sum(a*b for a, b in zip(y, triple_lower))
                   - sum(r*(upper_dot if r >= 0 else bound)
                         for r, bound in zip(residual, lower))
                   for lower, triple_lower, loss in states)

    stack = [(root, 0)]
    nodes = leaves = deepest = 0
    least = None
    while stack:
        simplex, depth = stack.pop()
        nodes += 1
        deepest = max(deepest, depth)
        for witness in witnesses:
            if all(margin(v, witness) > 0 for v in simplex):
                bound = min(margin(v, witness) for v in simplex)
                least = bound if least is None else min(least, bound)
                leaves += 1
                break
        else:
            assert depth < 22, f"unresolved weight simplex at depth {depth}: {simplex}"
            i, j = max(EDGES, key=lambda edge: sum((a-b)**2 for a, b in
                       zip(simplex[edge[0]], simplex[edge[1]])))
            midpoint = tuple((a+b)/2 for a, b in zip(simplex[i], simplex[j]))
            for endpoint in (i, j):
                stack.append((simplex[:endpoint]+(midpoint,)+simplex[endpoint+1:], depth+1))
    assert nodes == 2*leaves-1 and least > 0
    return nodes, leaves, deepest, least


def verify():
    t_max, c, u, _, disagreement = height_bounds(DELTA)
    assert 0 < disagreement < 1 and 14*DELTA+t_max < M
    assert intercept(PAIR_SLOPE, x_lo=F(".18"), x_hi=c,
                     upper_slope=F(".59")) >= PAIR_INTERCEPT
    check_angle_line(STAR_LINE, c, u)
    check_angle_line(TRIANGLE_LINE, c, u)
    tail = corner_tail_upper()
    assert tail < TAIL
    check_four_lines(verbose=False)
    radius_squared = 1/CAP_COSINE**2 - 1
    lines = [STAR_LINE if 0 in edge else TRIANGLE_LINE for edge in EDGES]
    factor = 1/(6*CAP_COSINE**3)
    coefficients = [factor*slope for _, _, slope in lines]
    threshold = (PI_LOWER*radius_squared/4 + excess_lower_bound(radius_squared) - 6*M
                 + factor*(sum(a for _, a, _ in lines) - (6-3*CAP_COSINE)*PI_UPPER))
    states = deficit_states(t_max, u, disagreement)
    nodes, leaves, depth, least = cover_weights(c, u, states, coefficients, threshold)
    assert least > F(18569, 10**9)
    print(f"candidate cap {M}; t <= {t_max}, dots in [{-c}, {u}], D0 <= {disagreement}")
    print(f"corner-tail upper bound {float(tail):.10f}; {len(states)} deficit vertices")
    print(f"{nodes} weight simplices, {leaves} certified leaves, depth {depth}, no unresolved leaves")
    print(f"minimum certified margin > {int(least*10**9)}/1000000000")


if __name__ == "__main__":
    verify()
