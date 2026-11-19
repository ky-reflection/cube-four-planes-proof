"""Exclude all fourth planes for boxes of the first three planes.

Normals are scaled to have one coefficient fixed at +1 or -1. The first
three offsets have symmetric bounds about zero in the centered test.
The fourth plane has arbitrary direction and offset.
"""
from fractions import Fraction as F
from itertools import product

from seven_box_certificate import CUBE, OCTANTS, ONE, ZERO, interval_negate, volume

SIGN_REPRESENTATIVES = tuple((1, b, c) for b, c in product((1, -1), repeat=2))
SINGLE_POSITIVE = tuple(tuple(1 if i == j else -1 for i in range(3)) for j in range(3))
ZERO_OFFSETS = (ZERO,) * 3
EQUAL_WEIGHTS = (ONE,) * 3


def exceeds_candidate(value, copies=1):
    """Test exactly whether value is above copies * (1 + 2 sqrt(2)) / 48."""
    remainder = 48 * value - copies
    return remainder > 0 and remainder**2 > 8 * copies**2


def signed_row(row, sign):
    return row if sign == 1 else tuple(interval_negate(interval) for interval in row)


def combine_rows(rows, weights):
    """Enclose the coefficients of a positive linear combination of rows."""
    return tuple((sum(weight * row[j][0] for weight, row in zip(weights, rows)),
                  sum(weight * row[j][1] for weight, row in zip(weights, rows)))
                 for j in range(3))


def common_volume(rows, offsets, precise):
    """Volume common to the indicated positive sides across the whole box."""
    if not precise:
        halfspaces = []
        for row, offset in zip(rows, offsets):
            middle = tuple((lo + hi) / 2 for lo, hi in row)
            error = sum(hi - lo for lo, hi in row) / 4
            halfspaces.append((*middle, offset + error))
        return volume(CUBE, halfspaces)
    result = ZERO
    for octant in OCTANTS:
        halfspaces = []
        for row, offset in zip(rows, offsets):
            endpoint = tuple(lo if bounds[0] == 0 else hi
                             for (lo, hi), bounds in zip(row, octant))
            halfspaces.append((*endpoint, offset))
        result += volume(octant, halfspaces)
    return result


def normalize_box(normal_box, offset_bounds):
    if len(normal_box) != 3 or any(len(row) != 3 for row in normal_box):
        raise ValueError("expected three rows of three coefficient intervals")
    if len(offset_bounds) != 3 or any(not isinstance(v, (int, F)) or v < 0
                                      for v in offset_bounds):
        raise ValueError("expected three nonnegative rational offset bounds")
    if any(not isinstance(v, (int, F)) for row in normal_box for pair in row for v in pair):
        raise TypeError("coefficient endpoints must be rational")
    box = tuple(tuple((F(lo), F(hi)) for lo, hi in row) for row in normal_box)
    for row in box:
        if any(lo > hi or lo < -1 or hi > 1 for lo, hi in row):
            raise ValueError("expected ordered coefficient intervals in [-1,1]")
        if not any(lo == hi and abs(lo) == 1 for lo, hi in row):
            raise ValueError("each normal needs a coefficient fixed at +1 or -1")
    offsets = tuple(F(v) for v in offset_bounds)
    if not exceeds_candidate(F(1, 12) - sum(offsets)):
        raise ValueError("offsets too large for the rank-deficient counting bound")
    return box, offsets


def certify_box(normal_box, offset_bounds=ZERO_OFFSETS, trial_weights=(EQUAL_WEIGHTS,),
                precise=True):
    """Return exclusion witnesses, or None if some fourth directions remain.

    Three independent central normals span the fourth normal. The signs of
    its three coordinates choose the unsplit pair of opposite cones. Once
    that pair is oriented as +++/---, these coordinates are nonnegative.
    Comparing the extreme ratios to a positive trial vector gives a cell
    at least as large as one of the three certified cells. For noncentral
    cuts use symmetric eroded cones, contained in both allowed affine cells
    after reflection, so the fourth offset can still be arbitrary.
    """
    box, offsets = normalize_box(normal_box, offset_bounds)
    weights = tuple(tuple(v) for v in trial_weights)
    if any(len(v) != 3 or any(not isinstance(x, (F, int)) or x <= 0 for x in v)
           for v in weights):
        raise ValueError("each trial needs three positive rational weights")
    volumes = {}
    witnesses = []
    for signs in SIGN_REPRESENTATIVES:
        rows = tuple(signed_row(row, sign) for row, sign in zip(box, signs))
        lower = common_volume(rows, offsets, precise)
        if exceeds_candidate(lower, copies=2):
            return (("cone must be split", signs, lower),)
        volumes[signs] = lower
    for signs, lower in volumes.items():
        if exceeds_candidate(lower):
            witnesses.append(("unsplit cone", signs, lower))
            continue
        rows = tuple(signed_row(row, sign) for row, sign in zip(box, signs))
        for trial, side in product(weights, (1, -1)):
            fourth = signed_row(combine_rows(rows, trial), side)
            lower_bounds = []
            for pattern in SINGLE_POSITIVE:
                cone = tuple(signed_row(row, sign) for row, sign in zip(rows, pattern))
                lower = common_volume((*cone, fourth), (*offsets, ZERO), precise)
                if not exceeds_candidate(lower):
                    break
                lower_bounds.append(lower)
            if len(lower_bounds) == 3:
                witnesses.append(("trial direction", signs, trial, side, min(lower_bounds)))
                break
        else:
            return None
    return tuple(witnesses)


def diagonal_box(lo, hi):
    return tuple(tuple((ONE, ONE) if i == j else (lo, hi) for j in range(3))
                 for i in range(3))


def certify_band(lower, upper, radius, offset, max_depth=18):
    """Check a common off-diagonal coefficient plus six independent errors."""
    pending = [(lower, upper, 0)]
    closed = trial_count = deepest = 0
    while pending:
        lo, hi, depth = pending.pop()
        box = diagonal_box(lo - radius, hi + radius)
        proof = certify_box(box, (offset,) * 3, precise=depth >= 3)
        if proof is not None:
            closed += 1
            trial_count += any(witness[0] == "trial direction" for witness in proof)
            deepest = max(deepest, depth)
            continue
        if depth == max_depth:
            raise RuntimeError(f"unresolved band over [{lo}, {hi}]")
        mid = (lo + hi) / 2
        pending.extend(((lo, mid, depth + 1), (mid, hi, depth + 1)))
    return closed, trial_count, deepest


if __name__ == "__main__":
    count, trials, depth = certify_band(
        F(-10, 11), F(-10, 13), F(1, 500), F(1, 2000),
    )
    print(f"certified affine band: leaves={count}, using a trial={trials}, depth={depth}")
    radius = F(1, 200)
    proof = certify_box(diagonal_box(F(-5, 6) - radius, F(-5, 6) + radius),
                        (F(1, 1000),) * 3)
    assert proof is not None and any(witness[0] == "trial direction" for witness in proof)
    print(f"certified independent radius {radius} and offsets at most 1/1000")
