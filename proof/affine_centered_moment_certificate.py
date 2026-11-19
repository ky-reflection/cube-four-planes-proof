"""Upper bounds on centered cell combinations under the affine volume cap.

Print the 128 integer upper bounds (in the centered checker's units) for
standard input to its --affine-moments mode.
"""
from fractions import Fraction as F
from itertools import combinations

from affine_frame_height_certificate import capped_deficit_vertices, vertex_caps

Q = F(79759, 10**6)
T = F(8452, 10**6)
G = 7*Q - F(1, 2)
GRID = 2**40
EDGES = tuple(combinations(range(4), 2))
REPRESENTATIVES = (1, 2, 4, 8, 3, 5, 9)
CENTERED_UPPER = F(654538, 10**6)
CENTERED_LOWER = 1 - 6*Q - T - 2*(G+T)


def centered_pair_bounds():
    assert Q < 2*G <= 2*G+T < 2*Q
    assert CENTERED_LOWER >= F(379464, 10**6)
    for t in (F(0), T):
        for support in capped_deficit_vertices(Q, 2*G+t):
            cap = vertex_caps(support, G)
            bounds = []
            for i, j in EDGES:
                disagreement = 8*Q - sum(d for m, d in support if ((m>>i) ^ (m>>j)) & 1)
                lo = disagreement - (cap[i]-cap[j])**2/(3*F(379464, 10**6))
                hi = disagreement + (cap[i]+cap[j])**2/(3*(1-CENTERED_UPPER))
                bounds.append(((6*GRID*lo)//1, -((-6*GRID*hi)//1)))
            yield bounds


def thresholds():
    vertices = list(centered_pair_bounds())
    result = []
    for mask in range(128):
        chosen = [r for k, r in enumerate(REPRESENTATIVES) if (mask>>k) & 1]
        singleton = [r for r in chosen if r in REPRESENTATIVES[:4]]
        pairs = [r for r in chosen if r in REPRESENTATIVES[4:]]
        constant = 6 * GRID * (len(singleton)-len(pairs))
        coefficients = [
            sum(((r>>i) ^ (r>>j)) & 1 for r in pairs)
            - sum(not (r & ((1<<i) | (1<<j))) for r in singleton)
            for i, j in EDGES
        ]
        largest = max(sum(c*(hi if c>=0 else lo) for c, (lo, hi) in zip(coefficients, vertex))
                      for vertex in vertices)
        result.append(constant - ((-largest)//2))
    assert result[0] == 0 and result[-1] == 6*GRID
    return result


if __name__ == '__main__':
    print(*thresholds())
