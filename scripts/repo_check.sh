#!/usr/bin/env bash
set -u
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

fail=0
warn=0
pass() { printf '[PASS] %s\n' "$1"; }
warning() { printf '[WARN] %s\n' "$1"; warn=$((warn+1)); }
failure() { printf '[FAIL] %s\n' "$1"; fail=$((fail+1)); }

required=(
  README.md VERSION.md KNOWN_LIMITATIONS.md RELEASE_CHECKLIST.md
  0_ref/barcode_3digit.list 0_ref/hg19_cento.bed 0_ref/hg19_tandem.bed
  2_python_script/DuplexSeq_Mutation_call_aMsN_for_hg19.py
  2_python_script/bam_split_by_allele_barcode.py
  demo/run_core_demo.sh demo/input/toy_duplex.sam demo/input/toy_candidates.vcf
)
for f in "${required[@]}"; do
    [[ -f "$f" ]] && pass "found $f" || failure "missing $f"
done

for f in 1_Bulk_WGS/*.sh 1_Bulk_WGS/02_Bulk_variant_calling 2_Somatic_mutation_calling/*.sh 3_Phylogeny/*.sh demo/*.sh scripts/*.sh; do
    [[ -f "$f" ]] || continue
    if bash -n "$f"; then pass "bash syntax: $f"; else failure "bash syntax: $f"; fi
done

for f in 2_python_script/*.py; do
    if python3 -m py_compile "$f" 2>/dev/null; then pass "python syntax: $f"; else failure "python syntax: $f"; fi
done

if [[ -f 2_python_script/gtbulksearch_py3.py ]]; then
    pass "gtbulksearch_py3.py present"
else
    warning "gtbulksearch_py3.py is missing; matched-bulk filtering cannot be reproduced faithfully"
fi

if grep -RInE 'xxxxxx|xxxxxxxx' 1_Bulk_WGS 2_Somatic_mutation_calling 3_Phylogeny >/dev/null 2>&1; then
    warning "site-specific xxxxxx path placeholders remain in production templates"
else
    pass "no xxxxxx path placeholders"
fi

for tool in samtools bcftools bedtools bwa cutadapt seqtk; do
    command -v "$tool" >/dev/null 2>&1 && pass "tool available: $tool" || warning "tool not found on PATH: $tool"
done

python3 - <<'PY' >/dev/null 2>&1
import pysam
PY
[[ $? -eq 0 ]] && pass "Python package available: pysam" || warning "Python package not available: pysam"

printf '\nSummary: %d failure(s), %d warning(s)\n' "$fail" "$warn"
[[ "$fail" -eq 0 ]]
