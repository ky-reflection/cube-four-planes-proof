"""Eliminate the fourth affine plane using exact comparisons of coefficients.

Each trial cell that is too large excludes a cone of possible fourth
planes. Checking coverage of their union uses rational ratio constraints.

For fixed independent first-three planes this form of certificate exists
whenever every fourth plane gives a strict gap. Given any candidate w,
take a cell whose volume exceeds the threshold. Lower its favorable
coefficients slightly (keeping favorable zeros zero) and increase all
unfavorable coefficients to obtain a nearby rational trial v. The ratio
tests are then strict at w except for the automatic zero comparisons,
so they cover a neighborhood of w in the closed weight simplex. Cell
volume is continuous there, including at the constant positive form.
Compactness gives a finite cover. The search below is a bounded way
to find such certificates for a box of first-three planes.
"""
from fractions import Fraction as F
from itertools import combinations, product
import sys

from fourth_plane_certificate import (
    combine_rows, common_volume, diagonal_box, exceeds_candidate, signed_row,
)
from seven_box_certificate import ONE, ZERO, axis_data, interval_negate

SIGNS = tuple(product((1, -1), repeat=3))
SUPPORTS = tuple(support for size in (4, 3, 2, 1)
                 for support in combinations(range(4), size))


def oversized_cells(rows, offsets, trial, precise, requested=None):
    """Yield trial cells whose robust volume strictly exceeds the candidate."""
    fourth = combine_rows(rows, trial[:3])
    fourth_offset = tuple(sum(t * offset[end] for t, offset in zip(trial, offsets))
                          - trial[3] for end in (0, 1))
    patterns = product(SIGNS, (1, -1)) if requested is None else (requested,)
    for signs, side in patterns:
        # The coefficient of the constant trial term is positive in f.
        good = tuple(i for i, sign in enumerate((*signs, 1)) if sign == side)
        if not any(trial[i] for i in good):
            continue
        oriented = tuple(signed_row(row, sign) for row, sign in zip(rows, signs))
        upper = tuple((hi if sign == 1 else -lo) for (lo, hi), sign in zip(offsets, signs))
        last = fourth_offset[1] if side == 1 else -fourth_offset[0]
        lower = common_volume((*oriented, signed_row(fourth, side)), (*upper, last), precise)
        if exceeds_candidate(lower):
            yield signs, side, good, lower


def support_comparisons(good, trial, support):
    """Write a witness on a fixed support as a list of lower ratio bounds.

    If all good coordinates are scaled at least as much as all bad ones,
    the witness halfspace grows inside this cell. Zero coordinates follow
    from the cross-multiplied comparisons, before taking any ratios.
    """
    comparisons = []
    for i in good:
        for j in support:
            if j in good or trial[i] == 0:
                continue
            if i not in support or trial[j] == 0:
                return None
            comparisons.append((support.index(i), support.index(j), trial[i] / trial[j]))
    return tuple(comparisons)


def extend_ratio_bounds(bounds, i, j, upper):
    """Add a strict upper bound to the multiplicative shortest paths."""
    size = len(bounds)
    result = [list(row) for row in bounds]
    for start in range(size):
        if bounds[start][i] is None:
            continue
        for end in range(size):
            if bounds[j][end] is None:
                continue
            length = bounds[start][i] * upper * bounds[j][end]
            if result[start][end] is None or length < result[start][end]:
                result[start][end] = length
    return tuple(tuple(row) for row in result)


def uncovered_ratios(tests, bounds):
    """Return a consistent set of strict upper bounds, if any."""
    choices = None
    remaining = []
    for test in tests:
        if any(bounds[i][j] is not None and bounds[i][j] <= upper for i, j, upper in test):
            continue
        # Staying outside this witness means violating at least one of its
        # closed lower bounds. A strict upper bound closing a cycle with
        # product at most one is impossible.
        options = tuple((i, j, upper) for i, j, upper in test
                        if bounds[j][i] is None or upper * bounds[j][i] > ONE)
        if not options:
            return None
        if choices is None or len(options) < len(choices):
            if choices is not None:
                remaining.append(choices)
            choices = options
        else:
            remaining.append(options)
    if choices is None:
        return bounds
    for i, j, upper in choices:
        uncovered = uncovered_ratios(remaining, extend_ratio_bounds(bounds, i, j, upper))
        if uncovered is not None:
            return uncovered
    return None


def realize_ratios(bounds, support):
    """Produce rational weights satisfying a feasible system of strict ratios."""
    size = len(bounds)
    denominator = 2
    while True:
        shrink = ONE - F(1, denominator)
        point = [ONE] * size
        for _ in range(size):
            next_point = point.copy()
            for i, j in product(range(size), repeat=2):
                if i != j and bounds[i][j] is not None:
                    next_point[i] = min(next_point[i], point[j] * shrink * bounds[i][j])
            if next_point == point:
                total = sum(point)
                return tuple(point[support.index(i)] / total if i in support else ZERO
                             for i in range(4))
            point = next_point
        denominator *= 2


def uncovered_weights(witnesses):
    """Return rational weights outside the witness union, or None if covered."""
    for support in SUPPORTS:
        tests = []
        for good, trial in witnesses:
            comparisons = support_comparisons(good, trial, support)
            if comparisons is not None:
                tests.append(comparisons)
        bounds = tuple(tuple(ONE if i == j else None for j in support) for i in support)
        uncovered = uncovered_ratios(tests, bounds)
        if uncovered is not None:
            return realize_ratios(uncovered, support)
    return None


def round_away(point, good, denominator):
    """Make a trial halfspace smaller, so its witness excludes a neighborhood."""
    trial = []
    for i, weight in enumerate(point):
        floor = weight.numerator * denominator // weight.denominator
        component = max(0, floor - 1) if i in good else floor + 1
        trial.append(F(component, denominator))
    return tuple(trial)


def find_affine_witnesses(rows, offsets, precise=True, max_rounds=32):
    """Find a finite rational witness list, or return None on exhaustion."""
    witnesses = []
    ratios = []
    for c in (ZERO, F(1, 64), F(1, 16), F(1, 4), ONE):
        trial = (ONE, ONE, ONE, c)
        for signs, side, good, lower in oversized_cells(rows, offsets, trial, precise):
            witnesses.append((signs, side, trial, lower))
            ratios.append((good, trial))
    for _ in range(max_rounds):
        point = uncovered_weights(ratios)
        if point is None:
            return tuple(witnesses)
        added = False
        for signs, side, good, _ in oversized_cells(rows, offsets, point, precise):
            for power in (6, 8, 10, 12, 16, 24):
                trial = round_away(point, good, 2**power)
                for _, _, _, lower in oversized_cells(rows, offsets, trial, precise,
                                                       requested=(signs, side)):
                    witnesses.append((signs, side, trial, lower))
                    ratios.append((good, trial))
                    added = True
                if added:
                    break
            if added:
                break
        if not added:
            return None
    return None


def verify_affine_witnesses(rows, offsets, witnesses, precise=True):
    """Verify the given trial cells without running witness discovery."""
    ratios = []
    for signs, side, trial, expected in witnesses:
        if (len(trial) != 4 or any(not isinstance(w, (int, F)) or w < 0 for w in trial)
                or not any(trial) or side not in (1, -1) or signs not in SIGNS):
            return False
        matches = tuple(oversized_cells(rows, offsets, tuple(map(F, trial)), precise,
                                       requested=(signs, side)))
        if len(matches) != 1 or matches[0][-1] != expected:
            return False
        ratios.append((matches[0][2], tuple(map(F, trial))))
    return uncovered_weights(ratios) is None


def certify_affine_box(normal_box, offset_intervals, precise=True, max_rounds=32):
    """Certify all orientations and offsets of plane four for the given box."""
    values = [v for row in normal_box for pair in row for v in pair]
    values.extend(v for pair in offset_intervals for v in pair)
    if not all(isinstance(v, (int, F)) for v in values):
        raise TypeError("coefficient and offset endpoints must be rational")
    box = tuple(tuple((F(lo), F(hi)) for lo, hi in row) for row in normal_box)
    offsets = tuple((F(lo), F(hi)) for lo, hi in offset_intervals)
    if len(box) != 3 or any(len(row) != 3 for row in box) or len(offsets) != 3:
        raise ValueError("expected three rows of three coefficients and one offset")
    if any(lo > hi for row in (*box, offsets) for lo, hi in row):
        raise ValueError("interval endpoints must be ordered")
    if all(axis_data(box, axis)[1][0] == 0 for axis in range(3)):
        raise ValueError("box does not certify independence of the first three normals")

    proofs = []
    for signs in SIGNS:
        rows = tuple(signed_row(row, sign) for row, sign in zip(box, signs))
        oriented_offsets = tuple(pair if sign == 1 else interval_negate(pair)
                                 for pair, sign in zip(offsets, signs))
        lower = common_volume(rows, tuple(hi for _, hi in oriented_offsets), precise)
        if exceeds_candidate(lower):
            proofs.append((signs, lower))
            continue
        proof = find_affine_witnesses(rows, oriented_offsets, precise, max_rounds)
        if proof is None:
            return None
        if not verify_affine_witnesses(rows, oriented_offsets, proof, precise):
            raise RuntimeError("generated witnesses did not verify")
        proofs.append((signs, proof))
    return tuple(proofs)


def certify_affine_band():
    """Subdivide a slope and an offset scale, allowing independent errors."""
    pending = [(F(-17, 20), F(-41, 50), F(1, 200), F(1, 100), 0)]
    normal_error, offset_error = F(1, 2000), F(1, 2500)
    certified = deepest = 0
    while pending:
        lo, hi, a, b, depth = pending.pop()
        box = diagonal_box(lo - normal_error, hi + normal_error)
        offsets = ((a - offset_error, b + offset_error),
                   (-b - offset_error, -a + offset_error),
                   (a/2 - offset_error, b/2 + offset_error))
        if certify_affine_box(box, offsets, max_rounds=24) is not None:
            certified += 1
            deepest = max(deepest, depth)
            continue
        if depth == 11:
            raise RuntimeError(f"unresolved band over [{lo}, {hi}] and [{a}, {b}]")
        if hi - lo > 2*(b - a):
            mid = (lo + hi)/2
            pending.extend(((lo, mid, a, b, depth + 1), (mid, hi, a, b, depth + 1)))
        else:
            mid = (a + b)/2
            pending.extend(((lo, hi, a, mid, depth + 1), (lo, hi, mid, b, depth + 1)))
    return certified, deepest


if __name__ == "__main__":
    radius = F(1, 1000)
    box = diagonal_box(F(-5, 6) - radius, F(-5, 6) + radius)
    noise = F(1, 2000)
    offsets = tuple((center - noise, center + noise)
                    for center in (F(1, 100), F(-1, 100), F(1, 200)))
    proof = certify_affine_box(box, offsets, max_rounds=40)
    assert proof is not None
    trials = tuple(witness for _, value in proof if isinstance(value, tuple) for witness in value)
    lower = min((value for _, value in proof if isinstance(value, F)),
                default=ONE)
    lower = min((lower, *(witness[-1] for witness in trials)))
    print(f"certified nine-parameter affine box: trial cells={len(trials)}, "
          f"maximum volume at least {float(lower):.9f}")
    if sys.argv[1:] == ["--band"]:
        leaves, depth = certify_affine_band()
        print(f"certified affine slope and offset band: leaves={leaves}, depth={depth}")
