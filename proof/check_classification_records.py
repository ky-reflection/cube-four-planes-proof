"""Match the finite classification records to the independently enumerated domains."""
from collections import Counter, defaultdict
from itertools import permutations, product
from hashlib import sha256
import json
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


def portable_commands(items):
    """Make archived commands comparable across Python installs and checkout paths."""
    commands = []
    for job in items:
        argv = list(job["argv"])
        if argv[1:2] == ["-B"]:
            argv[0] = "python3"
        for i, value in enumerate(argv):
            path = Path(value)
            if path.is_absolute() and path.parent.name == "bin":
                argv[i] = "bin/" + path.name
        commands.append({**job, "argv": argv})
    return commands


def check_snapshot(snapshot, hashes):
    """Check both the coverage of the source manifest and its file digests."""
    input_suffixes = {".py", ".cpp", ".hpp", ".txt"}
    required = {p.name for p in ROOT.iterdir() if p.is_file() and p.suffix in input_suffixes}
    required.add("affine-two-edge-proof.md")  # The verifier parses its two-edge domain table.
    archived = {p.name for p in snapshot.iterdir()
                if p.is_file() and p.suffix in input_suffixes | {".md"}}
    listed = set(hashes)
    missing = (required | archived) - listed
    absent = listed - archived
    if missing or absent:
        raise ValueError(f"source manifest mismatch: unlisted={sorted(missing)}, absent={sorted(absent)}")
    for name, expected in hashes.items():
        if sha256((snapshot / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"source snapshot changed: {name}")


def check_replay(directory):
    """Check the complete fresh run against its frozen plan and reference records."""
    snapshot = directory / "sources"
    metadata = json.loads((directory / "inputs.json").read_text())
    check_snapshot(snapshot, metadata["sources"])
    log_dir = directory / "logs"
    def regenerate(script):
        return subprocess.check_output([sys.executable, "-B", script], cwd=snapshot, text=True)
    thresholds = regenerate("affine_centered_moment_certificate.py")
    assert len(list(map(int, thresholds.split()))) == 128
    assert thresholds == (log_dir / "thresholds.log").read_text()
    assert regenerate("check_classification_records.py") == (log_dir / "records.log").read_text()
    # Generate the domain plan using the snapshotted enumerator and inputs.
    # A separate interpreter avoids reusing the current directory's imports.
    import_plan = (
        "import json,sys; from dataclasses import asdict; from pathlib import Path; "
        "from verify_global_proof import jobs; "
        "print(json.dumps([asdict(j) for j in jobs(Path(sys.argv[1]),sys.stdin.read())]))"
    )
    plan = json.loads(subprocess.check_output(
        [sys.executable, "-B", "-c", import_plan, str(directory / "bin")],
        cwd=snapshot, text=True, input=thresholds
    ))
    archived = json.loads((directory / "plan.json").read_text())
    from dataclasses import asdict
    from verify_global_proof import jobs

    current_plan = [asdict(job) for job in jobs(directory / "bin", thresholds)]
    assert portable_commands(archived) == portable_commands(plan) == portable_commands(current_plan)
    names = [job["name"] for job in plan]
    assert len(names) == len(set(names)) == 2005
    for prefix, count in (("vertex-", 1731), ("mixed-", 162), ("pivot-1200-", 64)):
        assert sum(name.startswith(prefix) for name in names) == count, prefix
    rows = [json.loads(s) for s in (directory / "results.jsonl").read_text().splitlines()]
    by_name = {r["name"]: r for r in rows}
    assert len(by_name) == len(rows) == len(plan)
    summary = json.loads((directory / "summary.json").read_text())
    assert summary == {"total": len(plan), "failures": []}
    logs = {}
    for job in plan:
        name = job["name"]
        row = by_name[name]
        raw = (directory / "logs" / (name + ".log")).read_bytes()
        assert row["passed"] and row["exit"] == 0 and sha256(raw).hexdigest() == row["log_sha256"], name
        logs[name] = raw.decode()
        assert re.search(job["success"], logs[name], re.M), name

    for name in logs:
        if name.startswith("wedge") or name in ("angular-lower", "angular-upper"):
            assert re.search(r"^certified leaves=[1-9]\d*", logs[name]), name
    assert "all 256 Taylor lower forms are positive definite" in logs["angular-Hessian"]
    assert "Certified exact local derivative data" in logs["local_certificate"]
    assert "unresolved 0" in logs["affine_star_pair"]
    for bound in ("tetrahedron volume <= 0.008452", "Gram dots in [-0.573370, 0.143814]",
                  "absolute signed offsets <= 0.066765"):
        assert bound in logs["height"], bound
    assert "four-plane types: 3746" in logs["vertex_partition_census"]

    differences = []
    def compare(name, actual, expected):
        actual, expected = tuple(actual), tuple(expected)
        if actual != expected:
            differences.append((name, actual, expected))

    for fields in records("affine-vertex-type-records.txt"):
        codes = "-".join(fields[:4])
        if fields[4] == "mixed_faces":
            continue  # Each orbit is checked separately below.
        if fields[4] == "whole_faces":
            m = re.search(r"proved nodes (\d+) caps \d+ depth (\d+)", logs["whole-faces"])
            compare(codes, map(int, m.groups()), map(int, (fields[5], fields[7])))
        else:
            compare(codes, map(int, logs["vertex-" + codes].split()[-3:]), map(int, fields[5:]))
    for line in (ROOT / "affine-two-edge-proof.md").read_text().splitlines():
        m = re.search(r"<summary>.*固定棱对为 (\d+),(\d+)", line)
        if m:
            pair = m.groups()
        m = re.fullmatch(r"\| (\d+), (\d+) \| [TFS]{2} \| ([\d,]+) \| ([\d,]+) \| (\d+) \|", line)
        if m:
            name = "vertex-" + "-".join((*pair, m[1], m[2]))
            compare(name, map(int, logs[name].split()[-3:]), (int(x.replace(",", "")) for x in m.groups()[2:]))
    for fields in records("affine-mixed-chart-records.txt"):
        name = "mixed-" + "-".join([*fields[:4], fields[5]])
        m = re.search(r"proved nodes (\d+) depth (\d+)(?: leaves \d+ contractions (\d+))?", logs[name])
        n, depth, shrink = map(int, m.groups("0"))
        compare(name, (n-shrink, depth, shrink), map(int, [*fields[7:9], fields[9] if len(fields)>9 else 0]))
    for fields in records("affine-alternating-sector-records.txt"):
        if fields[1] != "project":
            continue
        name = f"pivot-1200-{int(fields[2]):02}"
        counters = dict(re.findall(r"([a-z]+) (\d+)", logs[name]))
        keys = "angular cap bm moment deepest game contractions deficit unresolved".split()
        compare(name, (int(counters[k]) for k in keys), map(int, fields[4:]))
    # Every binary split adds two children and every contraction adds one.
    for name, body in logs.items():
        if name.startswith("pivot-") and name != "pivot-1032":
            counts = {k: int(v) for k, v in re.findall(r"([a-z]+) (\d+)", body)}
            assert counts["angular"] + 1 == 2 * (counts["bm"] + counts["moment"]) + counts.get("contractions", 0), name
    # Alternate valid trees for 1001 are both documented in reciprocal-pivot-proof.md.
    reference = {
        "1000": {(1932811, 6394866, 50)},
        "1001": {(3195195, 8938300, 64), (3194576, 8938300, 64)},
        "1002": {(2033224, 6929879, 64)},
        "1230": {(13625249, 77662732, 44)},
    }
    for word, trees in reference.items():
        name = "pivot-" + word
        counts = dict(re.findall(r"([a-z]+) (\d+)", logs[name]))
        actual = tuple(int(counts[k]) for k in ("angular", "cap", "deepest"))
        if actual not in trees:
            differences.append((name, actual, sorted(trees)))
    counts = dict(re.findall(r"([a-z]+)=(\d+)", logs["pivot-1032"]))
    assert int(counts["visited"]) + 1 == 2 * sum(int(counts[k]) for k in ("individual", "combined", "regular"))
    compare("pivot-1032", (int(counts[k]) for k in ("visited", "depth")), (683, 28))
    print(f"{len(plan)} replay jobs passed and all log digests match.")
    print(f"{len(differences)} differences from historical tree statistics.")
    for name, actual, expected in differences:
        print(f"{name}: replay {actual}, record {expected}")


if __name__ == "__main__":
    if not __debug__:
        sys.exit("Run without -O so that all record checks are enabled.")
    check_vertices()
    check_alternating()
    if len(sys.argv) == 3 and sys.argv[1] == "--replay":
        check_replay(Path(sys.argv[2]).resolve())
    elif sys.argv[1:]:
        sys.exit("usage: check_classification_records.py [--replay OUTPUT_DIRECTORY]")
