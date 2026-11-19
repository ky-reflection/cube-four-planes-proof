"""Compare the release proof code with the source used for the recorded replay."""

import argparse
import ast
from collections import Counter
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def cpp_tokens(path, clang):
    source = path.read_text()
    # These files use only include and guard directives. Check them, including
    # directive boundaries, before comparing C++ tokens independently of layout.
    directives = [line.strip() for line in source.splitlines() if line.lstrip().startswith("#")]
    for line in directives:
        if not re.fullmatch(r"#(?:include\s+[<\"][^>\"]+[>\"]|(?:ifndef|define)\s+\w+|endif)", line):
            raise ValueError(f"unhandled directive in {path.name}: {line}")
    if re.search(r"\\\r?\n|\b__(?:LINE|FILE|COUNTER)__\b", source):
        raise ValueError(f"source-location dependent tokens in {path.name}")

    raw = subprocess.run(
        [clang, "-cc1", "-dump-raw-tokens", str(path)],
        check=True, capture_output=True, text=True,
    ).stderr
    locations = re.compile(r"\tLoc=<[^\n]+:\d+:\d+>\n")
    tokens = []
    for record in locations.split(raw):
        if not record:
            continue
        match = re.fullmatch(r"(\w+) '(.*)'\t[ \t\[\]A-Za-z]*", record, re.S)
        if not match:
            raise ValueError(f"unrecognized Clang record in {path.name}: {record!r}")
        kind, spelling = match.groups()
        if kind == "comment":
            continue
        if kind == "unknown":
            if not spelling.isspace():
                raise ValueError(f"unrecognized lexical token in {path.name}: {spelling!r}")
            continue
        tokens.append((kind, spelling))
    if not tokens:
        raise ValueError(f"empty Clang token stream for {path.name}")
    return directives, tokens


def expected_python_tree(path):
    tree = ast.parse(path.read_text())
    if path.name == "seven_box_certificate.py":
        # The only mathematical-code change is rational coercion at the input
        # boundaries of this auxiliary Python integrator. Match the insertions
        # exactly. Every original arithmetic statement must still match.
        additions = {
            "cap": "weights = tuple(map(F, weights))",
            "cap_quantile": "weights = tuple(map(F, weights))",
            "volume": ("bounds = tuple(tuple(map(F, axis)) for axis in bounds)\n"
                       "halfspaces = (tuple(map(F, row)) for row in halfspaces)"),
        }
        for function in tree.body:
            if isinstance(function, ast.FunctionDef) and function.name in additions:
                statements = ast.parse(additions[function.name]).body
                at = int(ast.get_docstring(function) is not None)
                function.body[at:at] = statements
    return ast.dump(tree)


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("snapshot", type=Path, help="sources/ from a verification run")
    args.add_argument("--proof", type=Path, default=ROOT / "proof")
    args.add_argument("--clang", default="clang")
    opts = args.parse_args()
    totals = Counter()
    normalized_inputs = []
    # The two orchestration scripts have a separate, documented version history.
    drivers = {"check_classification_records.py", "verify_global_proof.py"}
    for suffix in (".cpp", ".hpp", ".py", ".txt"):
        files = sorted(p.name for p in opts.snapshot.glob("*" + suffix) if p.name not in drivers)
        current = sorted(p.name for p in opts.proof.glob("*" + suffix) if p.name not in drivers)
        if not files or files != current:
            sys.exit(f"different {suffix} file sets")
        for name in files:
            old, new = opts.snapshot / name, opts.proof / name
            if suffix in {".cpp", ".hpp"}:
                equal = cpp_tokens(old, opts.clang) == cpp_tokens(new, opts.clang)
            elif suffix == ".py":
                new_tree = ast.dump(ast.parse(new.read_text()))
                old_tree = ast.dump(ast.parse(old.read_text()))
                equal = old_tree == new_tree
                if not equal and new_tree == expected_python_tree(old):
                    equal = True
                    normalized_inputs.append(name)
            else:
                equal = old.read_bytes() == new.read_bytes()
            if not equal:
                sys.exit(f"verification source differs: {name}")
            totals[suffix] += 1
    print("Compared source files:", dict(totals))
    if normalized_inputs:
        print("Checked rational input normalization:", ", ".join(normalized_inputs))


if __name__ == "__main__":
    main()
