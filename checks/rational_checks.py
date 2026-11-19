"""Cross-check two geometric kernels with exact rational polytope clipping."""

import argparse
from fractions import Fraction
from itertools import combinations
from math import prod
from pathlib import Path
import random
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
PROOF = ROOT / "proof"
sys.path.insert(0, str(PROOF))
from affine_vertex_volumes import normals_from_scores
from seven_box_certificate import CUBE, cap, cap_quantile, volume


def run_rows(binary, inputs):
    text = "\n".join(" ".join(map(str, row)) for row in inputs) + "\n"
    rows = subprocess.check_output([binary], text=True, input=text).splitlines()
    if len(rows) != len(inputs):
        raise ValueError("the C++ probe returned the wrong number of results")
    return [list(map(int, row.split())) for row in rows]


def check_symmetric_family():
    """Check the hand-derived formula on both sides of the optimal slope."""
    slopes = tuple(map(Fraction, ("7/6", "5/4", "3/2", "17/10", "19/10")))
    for t in slopes:
        normals = [(1, 1, 1), (t, -1, -1), (-1, t, -1), (-1, -1, t)]
        pair = (3 + 2*t - t*t) / 48
        small = (1 - t/2)**2 / 4
        other = (1 - 6*pair - 2*small) / 6
        for mask in range(16):
            expected = [0, other, pair, other, 0][bin(mask).count("1")]
            if mask in (1, 14):
                expected = small
            cuts = []
            for i, normal in enumerate(normals):
                sign = 1 if mask & (1 << i) else -1
                cuts.append((*(sign*x for x in normal), 0))
            actual = volume(CUBE, cuts)
            if not isinstance(actual, Fraction) or actual != expected:
                raise ValueError((t, mask, actual, expected))
    print("All 80 exact cell volumes agree with the article's symmetric-family formulas.")


def check_slices(binary):
    rng = random.Random(20261119)
    bound = 65536
    special = (0, 0, 1, -1, 2, -2, 32768, -32768, bound, -bound)
    cases = []
    for case in range(4096):
        cuts = []
        for k in range(4):
            if case % 5 == 0 and k and rng.randrange(2):
                previous = rng.choice(cuts)
                sign = rng.choice((-1, 1))
                offset = sign * previous[3] + rng.choice((0, 0, 0, 1, -1))
                cut = (*(sign * x for x in previous[:3]), max(-bound, min(bound, offset)))
            elif case % 3 == 0:
                cut = tuple(rng.choice(special) for _ in range(4))
            else:
                cut = tuple(rng.randint(-bound, bound) for _ in range(4))
            cuts.append(cut)
        cases.append(cuts)
    rows = run_rows(binary, [[x for cut in cuts for x in cut] for cuts in cases])
    for index, (cuts, bounds) in enumerate(zip(cases, rows)):
        exact = volume(CUBE, [tuple(map(Fraction, cut)) for cut in cuts]) * (1 << 24)
        low, high, face_low, face_high = bounds
        if not (low <= exact <= high and face_low <= exact <= face_high):
            raise ValueError(f"slice {index}: {cuts}, {bounds}, exact {exact}")
    print(f"{len(cases)} exact rational polytope volumes enclosed by both C++ integrators")


def check_interpolation(binary):
    rng = random.Random(6911392)
    grid = 1 << 15
    inputs, references = [], []
    for word in ("1000", "1001", "1002", "1032", "1200", "1230"):
        pivots = list(map(int, word))
        for pair, (i, j) in enumerate(combinations(range(4), 2)):
            for _ in range(32):
                scores = [[Fraction(0)] * 4 for _ in range(4)]
                row, sample = [word, pair], []
                for a in range(4):
                    for b in range(4):
                        if a == b:
                            continue
                        if b == pivots[a]:
                            scores[a][b] = Fraction(1)
                            continue
                        denominator = rng.choice((32, 64, 128, 256, 2, 4, 8, 16))
                        low = rng.randrange(denominator) * grid // denominator
                        high = low + grid // denominator
                        t = rng.choice(tuple(map(Fraction, ("1/2", "1/3", "1/16", "1", "0", "7/9"))))
                        if a in (i, j):
                            sample.append(t)
                        row.extend((low, high))
                        scores[a][b] = ((1-t) * low + t * high) / grid
                normals = normals_from_scores(scores)
                exact = 2 * volume(CUBE, [(*normals[i], 0), (*(-x for x in normals[j]), 0)])
                inputs.append(row)
                references.append((sample, exact))

    for index, (values, (sample, exact)) in enumerate(zip(run_rows(binary, inputs), references)):
        if len(values) != 40:
            raise ValueError(f"expected 16 vertex intervals and 4 error intervals in row {index}")
        weights = [prod(t if mask & (1 << k) else 1-t for k, t in enumerate(sample)) for mask in range(16)]
        limits = []
        for side in (0, 1):
            interpolation = sum(w * values[2*m + side] for m, w in enumerate(weights))
            error = sum(values[32 + 2*k + side] for k in range(4))
            limits.append((interpolation + error) / (1 << 40))
        if not limits[0] <= exact <= limits[1]:
            raise ValueError(f"direction box {index}: {inputs[index]}, {limits}, exact {exact}")
    print(f"{len(inputs)} exact centered pair volumes enclosed by the directional interpolant plus its certified error")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cxx", default="c++")
    args = parser.parse_args()
    # Integer coefficients giving a non-dyadic result catch accidental float division.
    exact = volume([(0, 1)] * 3, [(1, 0, 3, 1)])
    corner = cap((1, 1, 1), Fraction(1, 2))
    limits = cap_quantile((1, 1, 1), Fraction(1, 6))
    if not (isinstance(exact, Fraction) and exact == Fraction(5, 6)
            and isinstance(corner, Fraction) and corner == Fraction(1, 6)
            and all(isinstance(v, Fraction) for v in limits) and limits[0] <= Fraction(1, 2) <= limits[1]):
        raise ValueError((exact, corner, limits))
    print("Integer inputs give exact 5/6 and 1/6 volumes and a rational cap enclosure.")
    check_symmetric_family()
    with tempfile.TemporaryDirectory() as temp:
        for mode, check in (("slice", check_slices), ("interpolation", check_interpolation)):
            binary = Path(temp) / mode
            source = Path(__file__).with_name(f"cube-proof-{mode}-audit.cpp")
            subprocess.run([args.cxx, "-std=c++20", "-O2", "-UNDEBUG", "-fsanitize=undefined",
                            "-fno-sanitize-recover=all", "-I", str(PROOF), str(source), "-o", str(binary)],
                           check=True)
            check(binary)


if __name__ == "__main__":
    main()
