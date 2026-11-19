"""Enumerate affine vertex splits and verify the centered and affine classifications."""

from collections import Counter
from itertools import combinations, combinations_with_replacement, permutations, product
import sys


def mask(vertices):
    return sum(1 << v for v in vertices)


def unordered_side(bits):
    return min(bits, 255 ^ bits)


def cube_actions():
    for order in permutations(range(3)):
        for flip in range(8):
            image = [sum(((v >> order[k]) & 1) << k for k in range(3)) ^ flip
                     for v in range(8)]
            yield [unordered_side(mask(image[v] for v in range(8) if s & (1 << v)))
                   for s in range(256)]


def representatives(splits, count, actions):
    seen = set()
    for row in combinations_with_replacement(sorted(splits), count):
        if row not in seen:
            seen.update(tuple(sorted(action[s] for s in row)) for action in actions)
            yield row


def passes_vertex_test(row):
    cells = {}
    for v in range(8):
        pattern = tuple((s >> v) & 1 for s in row)
        same = cells.setdefault(pattern, [])
        if len(same) >= 3 or (v ^ 7) in same:
            return False
        same.append(v)
    return True


def passes_two_vertex_cap_tests(row, edge_splits):
    caps = [s if bin(s).count("1") == 2 else 255 ^ s for s in row if s in edge_splits]
    if len(set(caps)) != len(caps):
        return False
    for s in caps:
        opposite = mask(v ^ 7 for v in range(8) if s & (1 << v))
        if opposite in caps:
            return False
    return len(caps) < 4 or (caps[0] | caps[1] | caps[2] | caps[3]) != 255


def determinant_volume(vertices):
    """Six times the tetrahedron volume for four cube vertices."""
    a, b, c, d = [tuple((v >> j) & 1 for j in range(3)) for v in vertices]
    u, v, w = [tuple(x-y for x, y in zip(point, a)) for point in (b, c, d)]
    return abs(sum(u[j]*(v[(j+1) % 3]*w[(j+2) % 3] - v[(j+2) % 3]*w[(j+1) % 3])
                   for j in range(3)))


def common_hull_tests():
    # A pair cell cannot contain 1/3 and a triple cell cannot contain 1/6.
    forbidden = {2: set(), 3: set()}
    for bits in range(256):
        vertices = [v for v in range(8) if bits & (1 << v)]
        maximum = max((determinant_volume(t) for t in combinations(vertices, 4)), default=0)
        if len(vertices) >= 5 or maximum >= 2:
            forbidden[2].add(bits)
        if maximum >= 1:
            forbidden[3].add(bits)
    return forbidden


def passes_hull_tests(row, forbidden):
    for k in (2, 3):
        for cuts in combinations(row, k):
            for sides in product(*[(s, 255 ^ s) for s in cuts]):
                common = 255
                for side in sides:
                    common &= side
                if common in forbidden[k]:
                    return False
    return True


def mixed_chart_orbits(row):
    """Orbits of oriented dominant-axis choices, in the mixed cover's input order."""
    small = [s if bin(s).count("1") <= 4 else 255 ^ s for s in row]
    fixed = []
    for side in small:
        vertices = [v for v in range(8) if side & (1 << v)]
        fixed.append({k for k in range(3) if len({v >> k & 1 for v in vertices}) == 1})
    charts = []
    for index in range(81):
        axes = [(index // 3**i) % 3 for i in range(4)]
        if any(choices and k not in choices for k, choices in zip(axes, fixed)):
            continue
        for g in range(16):
            sides = [s if g >> i & 1 else s ^ 255 for i, s in enumerate(small)]
            if any(len({s >> v & 1 for s in sides}) == 1 for v in range(8)):
                continue
            directions = []
            for side, k in zip(sides, axes):
                steps = {(side >> (v ^ (1 << k)) & 1) - (side >> v & 1)
                         for v in range(8) if not (v >> k & 1)}
                steps.discard(0)
                assert len(steps) == 1
                directions.append(steps.pop())
            if any(axes[i] == axes[j] and directions[i] == directions[j]
                   and fixed[i] and fixed[j] for i, j in combinations(range(4), 2)):
                continue
            charts.append(tuple(sorted(zip(sides, axes))))
    lookup = {chart: i for i, chart in enumerate(charts)}
    assert len(lookup) == len(charts)
    actions = []
    for perm in permutations(range(3)):
        for flip in range(8):
            image = [flip ^ sum((v >> j & 1) << perm[j] for j in range(3)) for v in range(8)]
            masks = {side: mask(image[v] for v in range(8) if side >> v & 1)
                     for side in {*small, *(s ^ 255 for s in small)}}
            actions.append((perm, masks))
    seen = set()
    for index, chart in enumerate(charts):
        if index in seen:
            continue
        orbit = set()
        for perm, masks in actions:
            mapped = tuple(sorted((masks[side], perm[k]) for side, k in chart))
            if mapped in lookup:
                orbit.add(lookup[mapped])
        assert index == min(orbit) and orbit.isdisjoint(seen)
        seen.update(orbit)
        yield sorted(orbit)
    assert len(seen) == len(charts)


def main():
    if len(sys.argv) == 6 and sys.argv[1] == "--mixed-charts":
        for orbit in mixed_chart_orbits(tuple(map(int, sys.argv[2:]))):
            print(*orbit)
        return
    actions = list(cube_actions())
    templates = ((0, 1), (0, 1, 2), (0, 1, 2, 3), (0, 1, 2, 4))
    classes = [{a[mask(template)] for a in actions} for template in templates]
    edges, triples, faces, stars = classes
    splits = edges | triples | faces | stars
    assert [len(c) for c in classes] == [12, 24, 3, 4] and len(splits) == 43
    three = list(representatives(splits, 3, actions))
    four = list(representatives(splits, 4, actions))
    vertex_ok = [row for row in four if passes_vertex_test(row)]
    cap_ok = [row for row in vertex_ok if passes_two_vertex_cap_tests(row, edges)]
    balanced = [row for row in four if set(row) <= faces | stars]
    centered = [row for row in balanced if passes_vertex_test(row)]
    repeats = [row for row in centered if any(row.count(s) > 1 for s in stars)]
    alternating = [row for row in centered if set(row) == stars]
    one_face = [row for row in centered if sum(s in faces for s in row) == 1
                and len(set(row) & stars) == 3]
    multi_face = set(centered) - set(alternating) - set(repeats) - set(one_face)
    counts = (len(three), len(four), len(vertex_ok), len(cap_ok),
              len(balanced), len(centered), len(repeats), len(alternating),
              len(one_face), len(multi_face))
    assert counts == (377, 3746, 3496, 3177, 22, 20, 8, 1, 1, 10)
    edge_triples = {row for row in three if set(row) <= edges
                    and passes_two_vertex_cap_tests(row, edges)}
    assert edge_triples == {(3, 5, 10), (3, 5, 17), (3, 5, 34), (3, 5, 48),
                            (3, 5, 119), (3, 12, 80), (3, 68, 95)}
    def edge_count(row):
        return sum(s in edges for s in row)
    at_most_two = [row for row in cap_ok if edge_count(row) <= 2]
    forbidden = common_hull_tests()
    hull_ok = [row for row in at_most_two if passes_hull_tests(row, forbidden)]
    star_ok = [row for row in hull_ok if not any(row.count(s) > 1 for s in stars)]
    stages = []
    for rows in (cap_ok, at_most_two, hull_ok, star_ok):
        counts_by_edges = Counter(map(edge_count, rows))
        stages.append([counts_by_edges[k] for k in range(5)])
    assert stages == [[1071, 1374, 618, 108, 6], [1071, 1374, 618, 0, 0],
                      [791, 785, 207, 0, 0], [757, 774, 206, 0, 0]]
    if sys.argv[1:] == ["--one-edge"]:
        one = []
        for row in star_ok:
            if edge_count(row) != 1:
                continue
            edge, = (v for v in row if v in edges)
            others = [v for v in row if v not in edges]
            one.append((3, *min(tuple(sorted(g[v] for v in others)) for g in actions if g[edge] == 3)))
        assert len(one) == len(set(one)) == 774
        for row in sorted(one):
            print(*row)
        return
    if sys.argv[1:] == ["--no-edges"]:
        zero = [row for row in star_ok if edge_count(row) == 0]
        assert len(zero) == 757
        for row in zero:
            print(*row)
        return
    two_edges = Counter()
    anchored = {"parallel": [], "skew": []}
    for row in star_ok:
        if edge_count(row) != 2:
            continue
        axes = []
        for split in (s for s in row if s in edges):
            side = split if bin(split).count("1") == 2 else 255 ^ split
            a, b = [v for v in range(8) if (side >> v) & 1]
            axes.append(a ^ b)
        kind = "parallel" if axes[0] == axes[1] else "skew"
        two_edges[kind] += 1
        target = (3, 12) if kind == "parallel" else (3, 68)
        pair = [s for s in row if s in edges]
        rest = [s for s in row if s not in edges]
        images = [tuple(sorted(g[s] for s in rest)) for g in actions
                  if tuple(sorted(g[s] for s in pair)) == target]
        anchored[kind].append((*target, *min(images)))
    assert two_edges == {"parallel": 74, "skew": 132}
    if sys.argv[1:] == ["--two-edges"]:
        for kind in ("parallel", "skew"):
            pair = (3, 12) if kind == "parallel" else (3, 68)
            stabilizer = [g for g in actions if tuple(sorted(g[s] for s in pair)) == pair]
            direct = []
            for a, b in combinations_with_replacement(sorted(triples | faces | stars), 2):
                if a == b and a in stars:
                    continue
                if min(tuple(sorted((g[a], g[b]))) for g in stabilizer) != (a, b):
                    continue
                row = (*pair, a, b)
                if passes_vertex_test(row) and passes_hull_tests(row, forbidden):
                    direct.append(row)
            assert direct == sorted(anchored[kind]) and len(direct) == two_edges[kind]
            for row in direct:
                print(*row)
        return
    normals = dict(X=(1, 0, 0), Y=(0, 1, 0), Z=(0, 0, 1),
                   A=(1, 1, 1), B=(1, 1, -1), C=(1, -1, 1), D=(1, -1, -1))
    codes = {label: unordered_side(mask(v for v in range(8)
                            if sum(a * (2 * ((v >> i) & 1) - 1)
                                   for i, a in enumerate(normal)) > 0))
             for label, normal in normals.items()}
    assert {codes[c] for c in "XYZ"} == faces and {codes[c] for c in "ABCD"} == stars
    groups = (("ZZAD", "ZZAB", "ZYAD", "ZXAD"), ("ZZZA", "ZZYA", "ZYXA"),
              ("ZZZY", "ZZYY", "ZZYX"))
    for f, group in enumerate(groups, start=2):
        assert all(sum(c in "XYZ" for c in word) == f for word in group)
        represented = {tuple(sorted(codes[c] for c in word)) for word in group}
        assert len(represented) == len(group)
        assert represented == {row for row in multi_face if sum(s in faces for s in row) == f}
    assert set(balanced) - set(centered) == {(min(faces),) * 4, (min(stars),) * 4}
    print("three-plane types: 377; four-plane types: 3746 -> 3496 -> 3177")
    print("affine four-plane types by edge count (0 through 4):")
    for label, stage in zip(("previous tests", "three-edge theorem", "common hulls",
                             "affine star pair"), stages):
        print(f"{label}: {stage}; total {sum(stage)}")
    print(f"two edges: {two_edges['parallel']} parallel, {two_edges['skew']} skew")
    print("centered: 22 types, including two constant types; among the other 20:")
    print("repeated star: 8; alternating: 1; one face / three distinct stars: 1")
    for f, group in enumerate(groups, start=2):
        print(f"{f} faces: {' '.join(group)}")


if __name__ == "__main__":
    main()
