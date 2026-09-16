# Reviewer smoke demo

This demo is intentionally tiny and deterministic. It tests custom code logic without requiring a reference genome, BAM files, BWA, SAMtools, CellPhy, Ginkgo, or a cluster scheduler.

## Requirements

- Python 3
- NumPy (used by `De_novo_estimation_by_variant_position.py`)

## Run

From the repository root:

```bash
python3 demo/run_demo.py
```

Expected result:

```text
PASS  local-sequence filter
PASS  de-novo loci aggregation
PASS  de-novo position estimation
PASS  scSNV assignment logic
DEMO PASS
```

The script also reports elapsed wall-clock time and writes generated files under `demo/output/`.

## What is tested

1. `localseqfilter_py3.py`: a normal-complexity SNV context passes and a homopolymer-rich context is removed.
2. `de_novo_loci.py`: loci are aggregated across two toy cells.
3. `De_novo_estimation_by_variant_position.py`: position-dependent de-novo estimates are generated from tiny VCF-like call tables.
4. `scSNV_assignment` inside `DuplexSeq_Mutation_call_aMsN_for_hg19.py`: the exact function body is extracted with Python's AST and tested on representative count patterns (`Mutation`, `Damage_0`, `Damage_swap_strand_filtered`, and `Single_strand_variant_allele_filtered`). This avoids importing Pysam or needing a BAM while still executing the actual decision function from the production source file.

## Scope limitation

This is a code smoke test, not an end-to-end biological validation. A full FASTQ/BAM-to-final-call demo would additionally require large third-party tools/resources and the currently missing `gtbulksearch_py3.py`. If the journal/editor treats the complete pipeline as the central new software contribution, provide a compact biologically representative end-to-end example as an additional demo.
