# scNanoSeq / scDuplexSeq

Somatic mutation calling and lineage construction from single-cell Duplex-seq data.

This repository contains the analysis scripts used for bulk-WGS processing, single-cell Duplex-seq somatic mutation calling, CellPhy-based phylogeny construction, and Ginkgo-based CNV analysis. The original analysis was designed for a SLURM/HPC environment and several scripts are historical templates with study-specific paths. The repository now also includes a small deterministic reviewer demo, dependency/environment capture utilities, and documentation intended to satisfy a Nature-style Code and Software Submission Checklist.

> **Important:** this review-preparation snapshot is not yet publication-complete. The author must still provide the items marked **AUTHOR ACTION REQUIRED** in [`REVIEWER_CHECKLIST.md`](REVIEWER_CHECKLIST.md), especially the missing `gtbulksearch_py3.py`, the actual tested software versions/OS, the license, and a release tag/DOI. Do not claim full reproducibility until those are resolved.

## Repository map

| Path | Purpose |
|---|---|
| `1_Bulk_WGS/` | Bulk WGS mapping and heterozygous SNP calling/filtering |
| `2_Somatic_mutation_calling/` | Single-cell preprocessing, mapping, barcode splitting, mutation calling, filtering, and summary |
| `2_python_script/` | Custom Python scripts used by the somatic-calling and phylogeny stages |
| `3_Phylogeny/` | Construction of a multi-cell VCF and CellPhy invocation |
| `4_CNV/` | Ginkgo-based CNV analysis notes and modified processing script |
| `0_ref/` | Small reference/support files distributed with this repository; the human reference FASTA is not bundled |
| `demo/` | Tiny deterministic demo of custom code paths that do not require large sequencing files |
| `docs/` | Installation, usage, workflow, provenance, and manuscript-facing templates |
| `scripts/` | Reviewer preflight and environment/version capture utilities |

## What the reviewer demo does

The demo intentionally **does not rerun BWA, SAMtools, Picard, CellPhy, or Ginkgo**. Those are third-party tools, not the novel software being reviewed. Instead it exercises custom/repository logic on tiny synthetic inputs: local-sequence filtering, de-novo locus aggregation, de-novo burden estimation, and the exact `scSNV_assignment` decision function extracted from the mutation caller.

Run:

```bash
python3 demo/run_demo.py
```

The command writes outputs to `demo/output/`, compares them byte-for-byte (after newline normalization) with `demo/expected/`, and exits non-zero if any result differs. On a normal modern computer this tiny smoke demo should complete in seconds. See [`demo/README.md`](demo/README.md).

This smoke demo is useful for code review, but it is **not a substitute for a biologically representative end-to-end example**. If the paper presents the full FASTQ/BAM-to-mutation workflow as a central contribution, provide a small de-identified real or simulated sequencing example that reaches the central mutation-calling output. That item is marked in the reviewer checklist.

## Installation and system requirements

Start with [`docs/INSTALLATION.md`](docs/INSTALLATION.md). The historical scripts require a Unix-like/HPC environment and external command-line tools including BWA, SAMtools/BCFtools, BEDTools, Seqtk, Cutadapt, Picard, Java/SnpSift, Python 3/Pysam, SLURM, CellPhy, and Ginkgo/R for the CNV stage. The exact versions actually used in the manuscript must be supplied by the authors rather than guessed from the scripts.

Capture the real versions from the analysis environment with:

```bash
bash scripts/collect_environment.sh > docs/TESTED_ENVIRONMENT.txt
```

Then review and commit that file. The script records OS/kernel information and versions of tools it can find.

## Running on your own data

The original shell scripts contain study-specific placeholders such as `xxxxxx`, `xxxxxxxx`, `AYY`, `BYY`, and `AXXX`. They are retained so the published analysis logic is not silently changed. Before reuse, follow [`docs/RUN_ON_YOUR_DATA.md`](docs/RUN_ON_YOUR_DATA.md), replace the path/sample placeholders, provide the required hg19 reference/resources, and verify every external dependency.

The central single-cell pipeline also calls `gtbulksearch_py3.py`, which is **not present in the supplied repository**. This is a hard reproducibility blocker until the original script is added or its exact public source is identified.

## Third-party code and provenance

`2_python_script/extract_tags.py` in the original repository is only a provenance note pointing to NanoSeq. The required upstream implementation is pinned in that note to NanoSeq commit `fd22dec943a3a9cc70643f079c40a8baf13a8e62`. See [`docs/THIRD_PARTY.md`](docs/THIRD_PARTY.md). A pinned fetch helper is provided as `scripts/fetch_nanoseq_extract_tags.sh`. Do not replace the upstream implementation with a newly written look-alike without validating that the scientific behavior is identical.

## Citation, license, and archived release

The public repository is <https://github.com/zonglab/scNanoSeq>. Before journal submission/acceptance, the authors should:

1. choose and add the intended software license;
2. create an immutable GitHub release/tag for the manuscript version;
3. archive that release in a DOI-minting repository such as Zenodo or an institutional repository; and
4. place the exact URL, release/tag, DOI, license, and any restrictions in the manuscript Code Availability statement.

A fill-in template is provided in [`docs/CODE_AVAILABILITY_TEMPLATE.md`](docs/CODE_AVAILABILITY_TEMPLATE.md).

## Reviewer-facing status

See [`REVIEWER_CHECKLIST.md`](REVIEWER_CHECKLIST.md) for a requirement-by-requirement status table. Items labeled **COMPLETED HERE** were added or verified in this review-preparation snapshot. Items labeled **AUTHOR ACTION REQUIRED** depend on facts or code that cannot be responsibly invented.

## Updating an existing clone

See [`LOCAL_UPDATE.md`](LOCAL_UPDATE.md) for a safe branch/backup/`rsync` procedure that preserves your existing `.git` history.
