#!/usr/bin/env python3
from __future__ import annotations

import ast
import difflib
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "demo"
INP = DEMO / "input"
EXP = DEMO / "expected"
OUT = DEMO / "output"
PY = ROOT / "2_python_script"


def run(cmd, cwd=None):
    p = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed ({p.returncode}): {' '.join(map(str, cmd))}\n"
            f"STDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
        )
    return p


def normalized_text(path: Path) -> str:
    return path.read_text().replace("\r\n", "\n").rstrip() + "\n"


def compare(name: str, got: Path, expected: Path):
    g = normalized_text(got)
    e = normalized_text(expected)
    if g != e:
        diff = "".join(
            difflib.unified_diff(
                e.splitlines(True), g.splitlines(True),
                fromfile=str(expected), tofile=str(got)
            )
        )
        raise AssertionError(f"{name} output mismatch:\n{diff}")
    print(f"PASS  {name}")


def demo_local_filter():
    p = run([
        sys.executable, str(PY / "localseqfilter_py3.py"),
        str(INP / "local_variants.vcf"), str(INP / "local_sequences.txt")
    ])
    got = OUT / "local_filtered.vcf"
    got.write_text(p.stdout)
    compare("local-sequence filter", got, EXP / "local_filtered.vcf")


def demo_denovo_loci():
    got = OUT / "de_novo_loci.txt"
    run([
        sys.executable, str(PY / "de_novo_loci.py"),
        str(INP / "denovo"), str(INP / "sample_list.txt"),
        ".calls.vcf", str(got)
    ])
    compare("de-novo loci aggregation", got, EXP / "de_novo_loci.txt")


def demo_denovo_estimation():
    got = OUT / "de_novo_estimation.txt"
    run([
        sys.executable, str(PY / "De_novo_estimation_by_variant_position.py"),
        str(INP),
        "hetero_mut.tsv", "hetero_ss.tsv", "denovo_ss.tsv", "denovo_mut.tsv",
        "100", "100", "10", str(got)
    ])
    compare("de-novo position estimation", got, EXP / "de_novo_estimation.txt")


def load_assignment_function():
    src = (PY / "DuplexSeq_Mutation_call_aMsN_for_hg19.py").read_text()
    tree = ast.parse(src)
    fn = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "scSNV_assignment"
    )
    module = ast.Module(body=[fn], type_ignores=[])
    ast.fix_missing_locations(module)
    ns = {}
    exec(compile(module, str(PY / "DuplexSeq_Mutation_call_aMsN_for_hg19.py"), "exec"), ns)
    return ns["scSNV_assignment"]


def demo_assignment_logic():
    f = load_assignment_function()
    cases = [
        ((0, 2, 0, 2), 4, 2, "Mutation"),
        ((4, 0, 0, 2), 4, 2, "Damage_0"),
        ((2, 2, 0, 2), 4, 2, "Damage_swap_strand_filtered"),
        ((0, 4, 0, 0), 4, 2, "Single_strand_variant_allele_filtered"),
    ]
    for counts, allele, strand, expected in cases:
        got = f(counts, allele, strand)
        if got != expected:
            raise AssertionError(f"scSNV_assignment{counts} -> {got!r}; expected {expected!r}")
    print("PASS  scSNV assignment logic")


def main():
    start = time.perf_counter()
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    demo_local_filter()
    demo_denovo_loci()
    demo_denovo_estimation()
    demo_assignment_logic()
    elapsed = time.perf_counter() - start
    print(f"DEMO PASS  ({elapsed:.3f} s)")


if __name__ == "__main__":
    main()
