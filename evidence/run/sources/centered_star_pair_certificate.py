"""Check the centered two-cut exclusion in the README exact-threshold roadmap."""

import argparse
from fractions import Fraction as F
from itertools import permutations

from seven_box_certificate import CUBE, volume


def verify():
    divisions = 32
    normals = [tuple(map(F, (divisions-i, divisions-j, i+j)))
               for i in range(divisions+1) for j in range(divisions-i+1)]
    best = F(0)
    checks = 0
    for p in normals:
        # Simultaneous coordinate permutation may sort the first cut.
        if not (p[0] <= p[1] <= p[2]):
            continue
        for q in normals:
            if p == q:
                continue
            disagreement = 2*volume(CUBE, [(*p, F(0)), (*(-x for x in q), F(0))])
            best = max(best, disagreement)
            checks += 1
    assert best == F(1, 3)
    # 0.51*sqrt(2)/divisions is the error for moving both normals.
    error = F(23, 1000)
    assert 2*F(51, 100)**2 < (divisions*error)**2
    cap = F(2, 25)
    assert best + error < 1 - 8*cap
    print(f"{len(normals)} normals, {checks} exact pairs, grid maximum {best}")
    print(f"continuous disagreement < {best+error} < {1-8*cap}; largest cell > {cap}")


def verify_recentered():
    # Extend the segment between the two normals to the boundary of the triangle.
    divisions = 64
    first = [(divisions, k, divisions-k) for k in range(divisions//2+1)]
    boundary = {v for k in range(divisions+1) for v in permutations((divisions, k, divisions-k))}
    largest = F(0)
    for p in first:
        for q in boundary:
            if p != q:
                pp, qq = tuple(map(F,p)), tuple(map(F,q))
                largest = max(largest, 2*volume(CUBE, [(*pp,F(0)), (*(-x for x in qq),F(0))]))
    assert largest == F(1,3)
    error = F(1,100)
    assert 6*F(51,100)**2 < (2*divisions*error)**2
    assert largest+error < 1-F(654538,10**6)
    print(f"{len(first)} x {len(boundary)} boundary grid normals, largest {largest}")
    print(f"same recentered star disagreement < {largest+error}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recentered", action="store_true")
    if parser.parse_args().recentered:
        verify_recentered()
    else:
        verify()
