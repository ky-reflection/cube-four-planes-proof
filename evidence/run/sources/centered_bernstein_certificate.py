"""Exact Bernstein exclusion of parameter boxes in the centered vertex chart.

Every row has one fixed entry 1. The other eight off-diagonal entries
vary independently in rational intervals, in row order. A floating-point
linear program proposes weights; only exact Fraction arithmetic accepts a
weighted volume inequality. See centered-bernstein-proof.md.
"""

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, product
from math import comb

from local_certificate import Quadratic as Q


ZERO = F(0)
TARGET = F(957106781186548, 10**15)
LABELS = range(4)
TRIANGLES = tuple((i, j) for i in LABELS for j in LABELS if i != j)
EDGES = tuple(combinations(LABELS, 2))
ATOM_COUNT = len(TRIANGLES) + len(EDGES)


def plus(first, second):
    result = first.copy()
    for powers, coefficient in second.items():
        result[powers] = result.get(powers, ZERO) + coefficient
    return {powers: coefficient for powers, coefficient in result.items() if coefficient}


def scale(poly, value):
    return {powers: coefficient * value for powers, coefficient in poly.items() if coefficient * value}


def times(first, second):
    result = {}
    for a, x in first.items():
        for b, y in second.items():
            powers = tuple(i + j for i, j in zip(a, b))
            result[powers] = result.get(powers, ZERO) + x * y
    return {powers: coefficient for powers, coefficient in result.items() if coefficient}


def midpoint_jet(poly, dimension):
    value = ZERO
    slope = [ZERO] * dimension
    for powers, coefficient in poly.items():
        half_term = coefficient / 2**sum(powers)
        value += half_term
        for axis, power in enumerate(powers):
            slope[axis] += 2 * power * half_term
    return value, slope


def bernstein(poly, degrees):
    """Convert the monomial coefficients on [0,1]^d by tensor products."""
    for axis, degree in enumerate(degrees):
        converted = {}
        for powers, value in poly.items():
            old = powers[axis]
            for new in range(old, degree + 1):
                index = powers[:axis] + (new,) + powers[axis + 1:]
                factor = F(comb(new, old), comb(degree, old))
                converted[index] = converted.get(index, ZERO) + factor * value
        poly = converted
    return poly


@dataclass(frozen=True)
class Atom:
    origin: F
    slope: tuple
    error_low: F
    error_high: F


def rational_atom(arguments, box, formula):
    """Linearize one rational atom and verify the remainder with Bernstein."""
    axes = tuple(sorted(x for x in arguments if x is not None))
    dimension = len(axes)
    constant = (0,) * dimension
    units = tuple(tuple(int(i == j) for i in range(dimension)) for j in range(dimension))
    arguments = tuple(
        {constant: F(1)} if axis is None else
        {constant: box[axis][0], units[axes.index(axis)]: box[axis][1] - box[axis][0]}
        for axis in arguments
    )
    numerator, denominators, selected = formula(arguments, dimension)
    center_num, gradient_num = midpoint_jet(numerator, dimension)
    center_den, gradient_den = midpoint_jet(denominators[selected], dimension)
    value = center_num / center_den
    local_slope = tuple((a - value * b) / center_den for a, b in zip(gradient_num, gradient_den))
    origin = value - sum(local_slope) / 2
    linear = {constant: origin, **dict(zip(units, local_slope))}
    lower, upper = ZERO, ZERO
    for denominator in denominators:
        residue = plus(numerator, scale(times(denominator, linear), -1))
        degrees = tuple(max(p[i] for p in (*residue, *denominator)) for i in range(dimension))
        num_control = bernstein(residue, degrees)
        den_control = bernstein(denominator, degrees)
        assert all(den_control.get(p, ZERO) > 0 for p in product(*(range(d+1) for d in degrees)))
        values = [num_control.get(p, ZERO) / b for p, b in den_control.items()]
        lower, upper = min(lower, min(values)), max(upper, max(values))
    slope = tuple(local_slope[axes.index(i)] if i in axes else ZERO for i in range(8))
    return Atom(origin, slope, lower, upper)


def triangle(arguments, dimension):
    a, b, c = arguments
    return times(a, a), (times(plus(a, b), plus(a, c)),), 0


def overlap(arguments, dimension):
    a, b, e, c, d, f = arguments
    delta = plus(times(a, d), scale(times(b, c), -1))
    u, v = plus(a, b), plus(c, d)
    common = times(u, v)
    forward = times(common, plus(times(plus(a, e), v), times(plus(d, f), u)))
    backward = times(common, plus(times(plus(c, f), u), times(plus(b, e), v)))
    control = tuple(bernstein(delta, (1,) * dimension).values()) or (ZERO,)
    if min(control) >= 0:
        return times(delta, delta), (forward,), 0
    if max(control) <= 0:
        return times(delta, delta), (backward,), 0
    selected = int(midpoint_jet(delta, dimension)[0] < 0)
    return times(delta, delta), (forward, backward), selected


def make_atoms(pivots, box):
    axes = iter(range(8))
    rows = tuple({j: None if pivots[i] == j else next(axes) for j in LABELS if j != i}
                 for i in LABELS)
    atoms = []
    for i, j in TRIANGLES:
        arguments = (rows[i][j], *(rows[i][k] for k in LABELS if k not in (i, j)))
        atoms.append(rational_atom(arguments, box, triangle))
    for i, j in EDGES:
        k, l = (a for a in LABELS if a not in (i, j))
        arguments = rows[k][i], rows[k][j], rows[k][l], rows[l][i], rows[l][j], rows[l][k]
        atoms.append(rational_atom(arguments, box, overlap))
    return atoms


def cell_rows():
    indices = {ij: i for i, ij in enumerate(TRIANGLES)}
    overlaps = {ij: i + len(TRIANGLES) for i, ij in enumerate(EDGES)}
    rows = [[0] * ATOM_COUNT for _ in range(7)]
    for i in LABELS:
        for j in LABELS:
            if j != i:
                rows[i][indices[j, i]] = 1
                rows[i][overlaps[tuple(sorted((i, j)))]] = -1
    for row, j in zip(rows[4:], (1, 2, 3)):
        k, l = (a for a in LABELS if a not in (0, j))
        for ij in ((0, j), (j, 0), (k, l), (l, k)):
            row[indices[ij]] = -1
        for a, b in product((0, j), (k, l)):
            row[overlaps[tuple(sorted((a, b)))]] = 1
    return rows


CELLS = cell_rows()


def weighted_margin(atoms, weights):
    assert len(weights) == 7 and all(w >= 0 for w in weights) and sum(weights) > 0
    coefficients = tuple(sum(c * w for c, w in zip(col, weights)) for col in zip(*CELLS))
    constant = 2 * sum(weights[4:]) - TARGET * sum(weights)
    slope = [ZERO] * 8
    for atom, coefficient in zip(atoms, coefficients):
        constant += coefficient * atom.origin + min(coefficient * atom.error_low, coefficient * atom.error_high)
        for i, s in enumerate(atom.slope):
            slope[i] += coefficient * s
    return (constant + sum(min(ZERO, s) for s in slope)) / sum(weights)


def simplex(rows, rhs, objective):
    """Floating-point witness discovery. The zero vector is feasible."""
    variables, constraints = len(objective), len(rows)
    tableau = [list(row) + [float(i == j) for j in range(constraints)] + [bound]
               for i, (row, bound) in enumerate(zip(rows, rhs))]
    tableau.append([-c for c in objective] + [0.] * (constraints + 1))
    basis = list(range(variables, variables + constraints))
    for _ in range(3000):
        enter = next((j for j, c in enumerate(tableau[-1][:-1]) if c < -1e-10), None)
        if enter is None:
            answer = [0.] * variables
            for i, j in enumerate(basis):
                if j < variables:
                    answer[j] = tableau[i][-1]
            return answer
        candidates = [(max(0., row[-1]) / row[enter], basis[i], i)
                      for i, row in enumerate(tableau[:-1]) if row[enter] > 1e-10]
        if not candidates:
            return None
        least = min(candidates)[0]
        _, _, leave = min((entry for entry in candidates if entry[0] <= least + 1e-11),
                          key=lambda entry: entry[1])
        pivot = tableau[leave][enter]
        tableau[leave] = row = [entry / pivot for entry in tableau[leave]]
        basis[leave] = enter
        for i, current in enumerate(tableau):
            if i != leave and current[enter] != 0:
                factor = current[enter]
                tableau[i] = [a - factor * b for a, b in zip(current, row)]
    return None


def discover_weights(atoms):
    # One nonnegative penalty for each signed remainder and each affine slope.
    columns = tuple(zip(*CELLS))
    penalties = ATOM_COUNT + 8
    constraints = [[1.] * 7 + [0.] * penalties]
    objective = [float(2 * (i >= 4) - TARGET + sum(c * a.origin for c, a in zip(row, atoms)))
                 for i, row in enumerate(CELLS)] + [-1.] * penalties

    def constrain(coefficients, penalty):
        row = [float(x) for x in coefficients] + [0.] * penalties
        row[7 + penalty] = -1.
        constraints.append(row)

    for i, (column, atom) in enumerate(zip(columns, atoms)):
        for error in (atom.error_low, atom.error_high):
            constrain((-error * c for c in column), i)
    for j in range(8):
        constrain((-sum(c * a.slope[j] for c, a in zip(row, atoms)) for row in CELLS), ATOM_COUNT + j)
    answer = simplex(constraints, [1.] + [0.] * (len(constraints) - 1), objective)
    if answer is None:
        return None
    weights = tuple(max(0, round(w * 10**8)) for w in answer[:7])
    return weights if sum(weights) > 0 else None


def check(pivots, box, weights=None):
    assert (4 * TARGET - 1)**2 > 8 and 4 * TARGET > 1
    assert len(pivots) == 4 and all(p in LABELS and p != i for i, p in enumerate(pivots))
    assert len(box) == 8 and all(0 < lo <= hi <= 1 for lo, hi in box)
    atoms = make_atoms(pivots, box)
    single = max((weighted_margin(atoms, tuple(int(i == j) for i in range(7))), j) for j in range(7))
    if weights is None:
        if single[0] > 0:
            return single[0], tuple(int(i == single[1]) for i in range(7)), single[0]
        weights = discover_weights(atoms)
        if weights is None:
            return single[0], None, single[0]
    return weighted_margin(atoms, weights), weights, single[0]


def verify_local_1000():
    """Check the linear separation of the candidate cone in chart 1000."""
    pivots = (1, 0, 0, 0)
    small = 2  # Interchanging labels 2 and 3 preserves this chart.
    axes = tuple((i, j) for i, j in TRIANGLES if j != pivots[i])
    entries = {(i, j): Q(-1, 1) if j == small else Q(1) for i, j in TRIANGLES}
    gradients = []
    for i, j in TRIANGLES:
        a = entries[i, j]
        b, c = (entries[i, k] for k in LABELS if k not in (i, j))
        volume = a*a / ((a+b)*(a+c))
        gradients.append(tuple(
            Q(0) if r != i else
            volume * (Q(2)/a - Q(1)/(a+b) - Q(1)/(a+c)) if s == j else
            -volume / (a+entries[i,s])
            for r, s in axes
        ))
    active = (0, 1, 3, 4, 5, 6)
    matrix = tuple(tuple(sum(c*g[j] for c, g in zip(CELLS[i], gradients)) for j in range(8))
                   for i in active)
    witnesses = (
        (0, -1, (409137, 90863, 0, 0, 458158, 41842)),
        (0, 1, (138071, 0, 0, 333333, 195262, 333333)),
        (1, -1, (146447, 146447, 0, 190744, 258182, 258182)),
        (2, -1, (381021, 0, 237959, 0, 65879, 315141)),
        (2, 1, (0, 138071, 0, 333333, 333333, 195262)),
        (3, -1, (292893, 0, 0, 190744, 325620, 190744)),
        (4, -1, (292893, 0, 0, 190744, 325620, 190744)),
        (5, -1, (195262, 97631, 0, 190744, 280661, 235702)),
        (6, -1, (292893, 0, 0, 190744, 325620, 190744)),
        (7, -1, (490956, 254522, 0, 254522, 0, 0)),
        (7, 1, (0, 0, 0, 79009, 460496, 460496)),
    )
    assert {(j, s) for j, s, _ in witnesses} == {
        (j, s) for j, (_, col) in enumerate(axes) for s in ((-1, 1) if col == small else (-1,))
    }
    for axis, sign, weights in witnesses:
        assert all(w >= 0 for w in weights) and sum(weights) > 0
        linear = tuple(sum(w*g for w, g in zip(weights, column)) for column in zip(*matrix))
        lower = sign * linear[axis]
        for j, (_, col) in enumerate(axes):
            if j != axis:
                if col == small:
                    lower -= abs(linear[j])
                elif linear[j].sign() > 0:
                    lower -= linear[j]
        assert (lower / sum(weights) - F(19,1000)).sign() > 0
    radius = F(1,1000)
    assert (Q(0,1)-2*radius-F(7,5)).sign() > 0
    assert 3*(F(125,49)+F(100,49)) < 14 and 4*F(125,49) < 14
    assert 14*radius < F(19,1000)
    return radius


def parse_interval(text):
    lo, hi = text.split(":")
    return F(lo), F(hi)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pivots", nargs="?", help="four fixed-column indices, for example 1230")
    parser.add_argument("intervals", nargs="*", type=parse_interval, metavar="LOW:HIGH")
    parser.add_argument("--weights", nargs=7, type=F, help="verify these weights without calling the LP")
    parser.add_argument("--local-1000", action="store_true", help="verify the local cone in chart 1000")
    arguments = parser.parse_args()
    if arguments.local_1000:
        if arguments.pivots or arguments.intervals or arguments.weights:
            parser.error("the local chart check takes no box or weights")
        print("certified chart 1000 candidate neighborhoods of radius", verify_local_1000())
    else:
        if arguments.pivots is None or len(arguments.intervals) != 8:
            parser.error("a box check needs four fixed columns and eight intervals")
        margin, weights, single = check(tuple(map(int, arguments.pivots)), arguments.intervals, arguments.weights)
        print("excluded" if margin > 0 else "unresolved", f"margin={float(margin):.12g}",
              f"best individual bound={float(single):.12g}", "weights=", weights)
