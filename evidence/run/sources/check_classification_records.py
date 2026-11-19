"""Match the finite classification records to the independently enumerated domains."""
from collections import Counter, defaultdict
from itertools import permutations, product
from pathlib import Path
import re
import subprocess
import sys

from vertex_partition_census import mixed_chart_orbits


ROOT = Path(__file__).resolve().parent
ALTERNATING = (23, 43, 77, 113)


def census(option):
    output = subprocess.check_output(
        [sys.executable, "-B", str(ROOT / "vertex_partition_census.py"), option], text=True
    )
    rows = [tuple(map(int, line.split())) for line in output.splitlines()]
    assert len(rows) == len(set(rows))
    return set(rows)


def records(name):
    for line in (ROOT / name).read_text().splitlines():
        if line and not line.startswith("#"):
            yield line.split()


def check_vertices():
    two = {}
    for line in (ROOT / "affine-two-edge-proof.md").read_text().splitlines():
        if line.startswith("<summary>"):
            match = re.search(r"固定棱对为 (\d+),(\d+)", line)
            pair = tuple(map(int, match.groups()))
        match = re.fullmatch(r"\| (\d+), (\d+) \| [TFS]{2} \| ([\d,]+) \| ([\d,]+) \| (\d+) \|", line)
        if match:
            a, b, n, leaves, depth = (int(x.replace(",", "")) for x in match.groups())
            row = (*pair, a, b)
            assert row not in two and n + 1 == 2 * leaves and depth >= 0
            two[row] = n
    assert set(two) == census("--two-edges") and len(two) == 206

    one, zero = census("--one-edge"), census("--no-edges")
    assert len(one) == 774 and len(zero) == 757 and one.isdisjoint(zero) and ALTERNATING in zero
    certified = {}
    for fields in records("affine-vertex-type-records.txt"):
        row = tuple(map(int, fields[:4]))
        n, leaves, depth = map(int, fields[5:])
        assert row not in certified and n + 1 == 2 * leaves and depth >= 0
        certified[row] = fields[4], n
    assert set(certified) == one | (zero - {ALTERNATING})

    charts = defaultdict(dict)
    for fields in records("affine-mixed-chart-records.txt"):
        row = tuple(map(int, fields[:4]))
        i, size, n = map(int, fields[5:8])
        assert i not in charts[row] and n > 0 and n % 2 == 1
        charts[row][i] = size, n
    assert set(charts) == {row for row, (method, _) in certified.items() if method == "mixed_faces"}
    for row, values in charts.items():
        assert {i: size for i, (size, _) in values.items()} == {
            orbit[0]: len(orbit) for orbit in mixed_chart_orbits(row)
        }
        assert sum(n for _, n in values.values()) + len(values) - 1 == certified[row][1]
    print("Vertex records: 206 two-edge, 774 one-edge, 756 nonalternating zero-edge.")


def check_alternating():
    labels = range(4)
    symmetries = list(permutations(labels))
    types = Counter()
    for f in product(labels, repeat=4):
        if all(f[i] != i for i in labels):
            rep = min(tuple(s[f[s.index(i)]] for i in labels) for s in symmetries)
            types["".join(map(str, rep))] += 1
    assert types == {"1000": 12, "1001": 12, "1002": 24, "1032": 3, "1200": 24, "1230": 6}

    domains = set()
    totals = Counter()
    for fields in records("affine-alternating-sector-records.txt"):
        word, mode, sector, patch = fields[:4]
        n, cap, bm, moment, depth, game, cuts, deficit, pending = map(int, fields[4:])
        key = mode, int(sector), None if patch == "-" else int(patch)
        assert word == "1200" and key not in domains and pending == 0
        assert n + 1 == 2 * (bm + moment) + cuts and 0 <= game <= moment
        assert cap >= 0 and 0 <= deficit <= moment and depth >= (6 if patch == "-" else 12)
        domains.add(key)
        if mode == "project":
            totals.update(angular=n, cap=cap, leaves=bm + moment, contractions=cuts)
    required = {("project", s, None) for s in range(64)}
    required.update(("inverse", s, p) for s in (47, 63) for p in range(64))
    assert domains == required
    assert totals["angular"] + 63 == 2 * totals["leaves"] - 1 + totals["contractions"]
    print("All 81 pivot choices give the six listed orbits.")
    print("1200: 64 sectors, 128 cross-check patches,", dict(totals))


if __name__ == "__main__":
    if not __debug__:
        sys.exit("Run without -O so that all record checks are enabled.")
    check_vertices()
    check_alternating()
