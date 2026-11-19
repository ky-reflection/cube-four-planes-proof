"""Count the initial four-cut orbits independently by Burnside's lemma."""

from collections import Counter
from itertools import combinations, permutations, product

VERTICES = list(product((0, 1), repeat=3))


def canonical_side(points):
    mask = sum(1 << VERTICES.index(p) for p in points)
    return min(mask, 255 ^ mask)


def one_cut_partitions():
    sides = set()
    for size in (2, 3, 4):
        for points in combinations(VERTICES, size):
            degrees = sorted(sum(sum(a != b for a, b in zip(x, y)) == 1 for y in points) for x in points)
            if degrees in ([1, 1], [1, 1, 2], [2, 2, 2, 2], [1, 1, 1, 3]):
                sides.add(canonical_side(points))
    return sides


def fixed_multisets(sides, axes, flips):
    def transform(mask):
        points = [tuple(p[a] ^ flip for a, flip in zip(axes, flips))
                  for i, p in enumerate(VERTICES) if mask & (1 << i)]
        return canonical_side(points)

    unseen = sides.copy()
    coefficients = [1, 0, 0, 0, 0]
    while unseen:
        first = unseen.pop()
        value, length = transform(first), 1
        while value != first:
            unseen.remove(value)
            value, length = transform(value), length + 1
        # Every orbit of this symmetry has a common multiplicity, contributing
        # the generating function 1 + t^length + t^(2*length) + ...
        for degree in range(length, 5):
            coefficients[degree] += coefficients[degree-length]
    return coefficients[4]


def main():
    sides = one_cut_partitions()
    counts = Counter(fixed_multisets(sides, axes, flips)
                     for axes in permutations(range(3)) for flips in VERTICES)
    numerator = sum(count * fixed for fixed, count in counts.items())
    if len(sides) != 43 or sum(counts.values()) != 48 or numerator != 48 * 3746:
        raise ValueError((len(sides), counts))
    print("Four-cut multisets fixed by a cube symmetry | number of symmetries")
    for fixed, count in sorted(counts.items()):
        print(f"{fixed:6d} | {count}")
    print(f"Burnside orbit count: {numerator}/48 = {numerator//48}")


if __name__ == "__main__":
    main()
