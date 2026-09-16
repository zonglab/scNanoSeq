# Minimal no-FASTQ demo

This directory contains a **synthetic, non-biological** fixture designed to test the custom scNanoSeq logic without distributing raw sequencing data.

## Why the demo starts from SAM/BAM rather than FASTQ

The repository-specific method starts after NanoSeq tags have been extracted and carried into the alignment. Shipping FASTQ would primarily retest `cutadapt` and `bwa`, which are mature external tools, while making the example much larger and potentially raising unnecessary data-sharing issues.

The smallest useful demo therefore starts from a transparent SAM file that already contains the alignment and tags expected by the custom code. `run_core_demo.sh` converts it to BAM locally with samtools.

## Files

```text
input/toy_reference.fa
    300-bp synthetic chr1 sequence. Position chr1:125 is C.

input/toy_duplex.sam
    Coordinate-sorted synthetic paired alignments. Four reads overlap chr1:125:
    two flag-99 reads and two flag-163 reads. All support T and carry compatible
    rb/mb tags. Their mates are present at the distal end of the fragment.

input/toy_candidates.vcf
    One candidate SNV: chr1:125 C>T.

expected/toy_call.expected.txt
    Expected output from the custom DuplexSeq caller using thresholds a=4, s=2.
```

The overlapping reads are constructed so that `bam_split_by_allele_barcode.py` places them in the `AAA` partition and assigns the shared molecule tag:

```text
MR = AAACCC
```

At the candidate site, the custom caller should count:

```text
(Fwt, Fmut, Rwt, Rmut) = (0, 2, 0, 2)
```

The site is 25 bp from the proximal read end and 76 bp from the distal molecule end, so it passes the >8-bp end-distance filters and is classified as:

```text
Mutation
```

## Run

Requirements:

- `samtools`
- Python 3
- `pysam`

From the repository root:

```bash
bash demo/run_core_demo.sh
```

A successful run ends with:

```text
PASS: synthetic barcode split + duplex SNV scoring matched expected output.
```

## What this demo validates

It validates the custom core path:

```text
SAM/BAM with rb/mb tags
        ↓
bam_split_by_allele_barcode.py
  - orientation normalization
  - MR construction
  - 3-nt barcode partition
        ↓
DuplexSeq_Mutation_call_aMsN_for_hg19.py
  - molecule grouping
  - strand counts
  - base/end-distance filters
  - a4s2 classification
```

## What it intentionally does not validate

This minimal fixture does not run:

- FASTQ tag extraction or adapter trimming;
- BWA mapping;
- Picard optical-duplicate processing;
- raw candidate discovery with samtools/bcftools;
- tandem/centromere/low-complexity filtering;
- matched-bulk filtering via the missing `gtbulksearch_py3.py`;
- SnpSift/dbSNP filtering;
- CellPhy; or
- Ginkgo CNV analysis.

For a public methods repository, that is a feature rather than a defect: the mandatory demo should exercise the code unique to this repository. External packages can be tested separately, and the missing matched-bulk helper must be recovered before a faithful full end-to-end demo can be claimed.

## Minimal publication toy-data recommendation

If you want one slightly larger example beyond this core fixture, the next useful level would add:

1. a tiny matched-bulk BAM covering the same synthetic locus;
2. a few candidate classes (`Mutation`, single-strand, damage/swap) instead of one site;
3. the recovered `gtbulksearch_py3.py` and an expected bulk-filtered output; and
4. two or three synthetic cells plus a prepared multi-sample VCF for an **optional** CellPhy smoke test.

There is no need to add FASTQ unless a reviewer explicitly requires validation of the laboratory barcode-extraction/mapping boundary.
