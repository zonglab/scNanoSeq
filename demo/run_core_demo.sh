#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
WORK="$ROOT/demo/work"
INPUT="$ROOT/demo/input"
EXPECTED="$ROOT/demo/expected/toy_call.expected.txt"

command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools is required." >&2; exit 1; }
python3 - <<'PY' >/dev/null 2>&1 || { echo "ERROR: Python package 'pysam' is required." >&2; exit 1; }
import pysam
PY

rm -rf "$WORK"
mkdir -p "$WORK/split"

# Materialize a tiny BAM from a human-readable synthetic SAM fixture.
samtools view -bS "$INPUT/toy_duplex.sam" | samtools sort -o "$WORK/toy_duplex.sorted.bam"
samtools index "$WORK/toy_duplex.sorted.bam"

# Exercise the repository-specific barcode normalization and MR-tag construction.
python3 "$ROOT/2_python_script/bam_split_by_allele_barcode.py" "$WORK" toy_duplex.sorted.bam "$WORK/split/toy"
samtools index "$WORK/split/toy_AAA.bam"

gzip -c "$INPUT/toy_candidates.vcf" > "$WORK/split/toy_candidates.vcf.gz"

# Exercise the repository-specific duplex SNV classifier with the production a4s2 thresholds.
python3 "$ROOT/2_python_script/DuplexSeq_Mutation_call_aMsN_for_hg19.py" "$WORK/split" toy_candidates.vcf.gz toy_AAA.bam 4 2 toy_call.txt

if diff -u "$EXPECTED" "$WORK/split/toy_call.txt"; then
    echo "PASS: synthetic barcode split + duplex SNV scoring matched expected output."
else
    echo "FAIL: demo output differs from expected output." >&2
    exit 1
fi
