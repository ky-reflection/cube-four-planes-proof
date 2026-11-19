"""Replay the computer-assisted steps of the global minimum and equality proof."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import fcntl
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time

from affine_frame_height_certificate import WIDE_SHARP_PAIRS
from affine_gram_bm_certificate import LEGACY_PAIRS, SHARP_PAIRS
from check_classification_records import census, portable_commands, records

ROOT = Path(__file__).resolve().parent
PYTHON = (sys.executable, "-B")
FLAGS = ("-O3", "-std=c++20", "-UNDEBUG", "-fsanitize=undefined", "-fno-sanitize-recover=all")


@dataclass(frozen=True)
class Job:
    name: str
    argv: tuple
    stdin: str = ""
    success: str = ""


def jobs(bins, thresholds):
    work = []

    def cpp(name, program, *args, stdin="", success=""):
        work.append(Job(name, (str(bins / program), *map(str, args)), stdin, success))

    def py(name, program, *args):
        work.append(Job(name, (*PYTHON, program + ".py", *map(str, args))))

    # Start long independent single-chart jobs early. All six orbits are required.
    cpp("pivot-1230", "affine_alternating_chart_cover", stdin=thresholds, success=r"1230 complete .*unresolved 0")
    for word, mode in (("1000", "plain"), ("1001", "columns"), ("1002", "flow")):
        cpp("pivot-" + word, "affine_alternating_candidate_cover", word, 96, 200000000, mode,
            stdin=thresholds, success=word + r" complete .*unresolved 0")
    cpp("pivot-1032", "centered_chart_cover", "--affine-moments", "1032", 32, 1000000,
        stdin=thresholds, success=r"chart=1032 .*unresolved=0 pending=0 .*complete")
    cpp("whole-faces", "affine_whole_face_cover", success=r"proved nodes ")
    py("height", "affine_frame_height_certificate", "--candidate", "--edges")
    py("angular-Hessian", "affine_local_merit_certificate", "--angular", bins / "affine_angular_merit", 35)
    cpp("angular-offsets", "affine_angular_offset_transfer", 35, 20, success=r"35 proved angular ")
    for mode, depth in (("horizontal", 30), ("vertical", 36)):
        cpp("seven-" + mode, "seven_scaled_charts_certificate", "--" + mode, depth,
            success=r"certified visited=\d+ unresolved=0 pending=0")
    for program in ("local_certificate", "antipodal_vertex_certificate", "affine_pyramid_row_certificate",
                    "affine_pyramid_matching_certificate", "affine_alternating_pair", "vertex_partition_census"):
        py(program, program)
    for program in ("affine_edge_triples", "affine_star_pair"):
        cpp(program, program)
    cpp("angular-lower", "angular_moment_certificate")
    cpp("angular-upper", "angular_moment_certificate", "--upper")
    cpp("wedge-default", "wedge_certificate")
    for option, anchors in (("--pair", LEGACY_PAIRS), ("--sharp-pair", SHARP_PAIRS + WIDE_SHARP_PAIRS)):
        for anchor in dict.fromkeys(anchors):
            cpp("wedge" + option + "-" + "-".join(map(str, anchor)), "wedge_certificate", option, *anchor)

    # Each mixed chart has its own enumeration index and method. The independent
    # record checker verifies that these representatives exhaust the symmetry orbits.
    mixed = list(records("affine-mixed-chart-records.txt"))
    for row in sorted(mixed, key=lambda r: -int(r[7])):
        ids, mode, index = row[:4], row[4], row[5]
        cpp("mixed-" + "-".join([*ids, index]), "affine_mixed_face_cover", *ids, index, 0, mode,
            success=rf"chart {index} axes \d+ g \d+ proved nodes ")
    for sector in range(64):
        cpp(f"pivot-1200-{sector:02}", "affine_alternating_candidate_cover",
            1200, 90, 120000000, "project", sector, stdin=thresholds,
            success=rf"sector {sector} 1200 complete .*unresolved 0")

    ordinary = []
    for row in records("affine-vertex-type-records.txt"):
        ids, mode = row[:4], row[4]
        if mode == "mixed_faces":
            if ids not in [r[:4] for r in mixed]:
                raise ValueError(("missing charts", ids))
        elif mode == "whole_faces":
            if ids != ["15", "15", "51", "85"]:
                raise ValueError(("unexpected whole-face chart", ids))
        else:
            ordinary.append((int(row[5]), ids, ["--" + mode]))
    two = census("--two-edges")
    # Historical node counts only set the scheduling order.
    for line in (ROOT / "affine-two-edge-proof.md").read_text().splitlines():
        m = re.search(r"<summary>.*固定棱对为 (\d+),(\d+)", line)
        if m:
            pair = m.groups()
        m = re.fullmatch(r"\| (\d+), (\d+) \| [TFS]{2} \| ([\d,]+) \| [\d,]+ \| \d+ \|", line)
        if m:
            ids = [*pair, m[1], m[2]]
            two.remove(tuple(map(int, ids)))
            ordinary.append((int(m[3].replace(",", "")), ids, []))
    if two:
        raise ValueError(("missing two-edge domains", two))
    for _, ids, args in sorted(ordinary, reverse=True):
        row = " ".join(ids)
        cpp("vertex-" + "-".join(ids), "affine_vertex_type_cover", *args, stdin=row + "\n",
            success=rf"^{row}\s+proved \d+ \d+ \d+$")
    if len({j.name for j in work}) != len(work):
        raise ValueError("duplicate domain in replay plan")
    return work


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new directory for the immutable input snapshot and results")
    parser.add_argument("--jobs", type=int, default=4, help="simultaneous certificate processes")
    parser.add_argument("--resume", action="store_true", help="continue an interrupted run in its original directory")
    args = parser.parse_args()
    if args.jobs < 1 or not __debug__:
        parser.error("positive worker count and enabled Python assertions required")
    output = args.output.resolve()
    if args.resume:
        if not output.is_dir():
            parser.error("--resume requires the existing run directory")
        with run_lock(output):
            return resume(output, args.jobs)
    output.mkdir(parents=True, exist_ok=False)
    sources = output / "sources"
    sources.mkdir()
    files = sorted(p for p in ROOT.iterdir() if p.suffix in (".cpp", ".hpp", ".py", ".md", ".txt"))
    hashes = {}
    for p in files:
        copy = sources / p.name
        shutil.copyfile(p, copy)
        hashes[p.name] = sha256(copy.read_bytes()).hexdigest()
    env = dict(os.environ)
    env.pop("PYTHONOPTIMIZE", None)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    compiler = shutil.which(os.environ.get("CXX", "c++"))
    if not compiler:
        raise RuntimeError("missing C++ compiler")
    def checked(cmd):
        return subprocess.check_output(cmd, cwd=sources, env=env, text=True)
    metadata = dict(sources=hashes, compiler=checked([compiler, "--version"]), flags=FLAGS,
                    python=sys.version, platform=platform.platform(), workers=args.jobs)
    (output / "inputs.json").write_text(json.dumps(metadata, indent=2) + "\n")
    # Construct and validate the domain plan with imports from the frozen source tree.
    subprocess.run([*PYTHON, str(sources / Path(__file__).name), "--snapshot-run",
                    str(output), str(args.jobs), compiler], check=True, env=env, cwd=sources)


def run_lock(output):
    handle = (output / ".run.lock").open("a+")
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        handle.close()
        raise RuntimeError("another verifier is writing to this run")
    return handle


def resume(output, workers):
    """Rebuild the frozen inputs and reuse only complete, previously successful jobs."""
    source = output / "sources"
    log_dir = output / "logs"
    metadata = json.loads((output / "inputs.json").read_text())
    for name, expected in metadata["sources"].items():
        if sha256((source / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"source snapshot changed: {name}")
    def checked(cmd, stdin=None):
        return subprocess.check_output(cmd, cwd=source, input=stdin, text=True)
    compiler = shutil.which(os.environ.get("CXX", "c++"))
    if (not compiler or checked([compiler, "--version"]) != metadata["compiler"]
            or tuple(metadata["flags"]) != FLAGS or sys.version != metadata["python"]):
        raise ValueError("compiler, flags, or Python changed; start a fresh run directory")
    thresholds = checked([*PYTHON, "affine_centered_moment_certificate.py"])
    regenerated = checked([*PYTHON, "check_classification_records.py"])
    if thresholds != (log_dir / "thresholds.log").read_text() or regenerated != (log_dir / "records.log").read_text():
        raise ValueError("the verification inputs disagree with the archived run")
    script = (
        "import sys,json; from dataclasses import asdict; from pathlib import Path; "
        "from verify_global_proof import jobs; "
        "print(json.dumps([asdict(j) for j in jobs(Path(sys.argv[1]),sys.stdin.read())]))"
    )
    plan = json.loads(checked([*PYTHON, "-c", script, str(output / "bin")], thresholds))
    if portable_commands(plan) != portable_commands(json.loads((output / "plan.json").read_text())):
        raise ValueError("the archived plan disagrees with the frozen verifier")
    work = [Job(**j) for j in plan]
    by_name = {j.name: j for j in work}
    results = {}
    journal = output / "results.jsonl"
    rows = journal.read_bytes().splitlines(keepends=True) if journal.exists() else []
    for index, row in enumerate(rows):
        try:
            record = json.loads(row)
        except (ValueError, UnicodeDecodeError):
            if index == len(rows) - 1 and not row.endswith(b"\n"):
                break
            raise
        name = record["name"]
        if not record["passed"] or record["exit"] != 0:
            continue
        job = by_name[name]
        raw = (log_dir / (name + ".log")).read_bytes()
        if name in results or sha256(raw).hexdigest() != record["log_sha256"] or not re.search(job.success, raw.decode(), re.M):
            raise ValueError(f"saved result failed integrity verification: {name}")
        results[name] = record
    if len(results) == len(work) and (output / "summary.json").exists():
        if json.loads((output / "summary.json").read_text()) == {"total": len(work), "failures": []}:
            print(f"All {len(work)} archived jobs passed.")
            return 0
    staging = journal.with_suffix(".new")
    staging.write_text("".join(json.dumps(r) + "\n" for r in results.values()))
    staging.replace(journal)
    restarts = output / "restarts"
    restarts.mkdir(exist_ok=True)
    stamp = str(time.time_ns())
    shutil.copyfile(__file__, restarts / (stamp + ".py"))
    (restarts / (stamp + ".json")).write_text(json.dumps({
        "workers": workers, "reused": len(results), "scheduled": len(work)-len(results),
        "driver_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }, indent=2) + "\n")
    print(f"Verified {len(results)} saved jobs; continuing {len(work) - len(results)} jobs.", flush=True)
    return execute_jobs(output, workers, compiler, work, results)


def run_snapshot(output, workers, compiler):
    bins = output / "bin"
    bins.mkdir()
    logs = output / "logs"
    logs.mkdir()
    def checked(cmd):
        return subprocess.check_output(cmd, cwd=ROOT, text=True)
    (logs / "records.log").write_text(checked([*PYTHON, "check_classification_records.py"]))
    thresholds = checked([*PYTHON, "affine_centered_moment_certificate.py"])
    if len(list(map(int, thresholds.split()))) != 128:
        raise ValueError("expected 128 centered-volume thresholds")
    (logs / "thresholds.log").write_text(thresholds)
    work = jobs(bins, thresholds)
    (output / "plan.json").write_text(json.dumps(portable_commands([asdict(j) for j in work]), indent=2) + "\n")
    return execute_jobs(output, workers, compiler, work, {})


def execute_jobs(output, workers, compiler, work, completed):
    bins, logs, sources = output / "bin", output / "logs", output / "sources"
    pending = [job for job in work if job.name not in completed]
    programs = {Path(a).name for j in work for a in j.argv if Path(a).parent == bins}
    def build(name):
        subprocess.run([compiler, *FLAGS, str(sources / (name + ".cpp")), "-o", str(bins / name)], check=True)
    with ThreadPoolExecutor(max_workers=min(workers, 6)) as pool:
        list(pool.map(build, sorted(programs)))

    def run(job):
        start = time.monotonic()
        log = logs / (job.name + ".log")
        with log.open("w") as stream:
            proc = subprocess.run(job.argv, cwd=sources, input=job.stdin, text=True,
                                  stdout=stream, stderr=subprocess.STDOUT)
        body = log.read_text()
        passed = proc.returncode == 0 and (not job.success or bool(re.search(job.success, body, re.M)))
        return dict(name=job.name, passed=passed, exit=proc.returncode,
                    seconds=round(time.monotonic()-start, 3), log_sha256=sha256(log.read_bytes()).hexdigest())

    failures = []
    with (output / "results.jsonl").open("a", buffering=1) as journal, ThreadPoolExecutor(workers) as pool:
        futures = [pool.submit(run, job) for job in pending]
        for n, future in enumerate(as_completed(futures), 1 + len(completed)):
            result = future.result()
            journal.write(json.dumps(result) + "\n")
            if not result["passed"]:
                failures.append(result["name"])
                print("FAILED", result["name"], flush=True)
            if n % 50 == 0 or n == len(work):
                print(f"{n}/{len(work)} completed, {len(failures)} failures", flush=True)
    (output / "summary.json").write_text(json.dumps(dict(total=len(work), failures=failures), indent=2) + "\n")
    return int(bool(failures))


if __name__ == "__main__":
    if sys.argv[1:2] == ["--snapshot-run"]:
        with run_lock(Path(sys.argv[2])):
            sys.exit(run_snapshot(Path(sys.argv[2]), int(sys.argv[3]), sys.argv[4]))
    sys.exit(main())
