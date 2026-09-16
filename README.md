# scNanoSeq

**Single-cell Duplex-seq somatic mutation calling, mutation-burden estimation, phylogeny reconstruction, and CNV analysis.**

This repository contains the analysis scripts used to process matched bulk WGS and single-cell Duplex-seq data, classify duplex-supported SNVs, estimate de novo mutation burden, construct multi-cell phylogenies with CellPhy, and run a modified Ginkgo CNV workflow.

> **Repository naming.** The GitHub repository is named `scNanoSeq`. The core single-cell caller was renamed from `NanoSeq_*` to `DuplexSeq_*` in the 2025 code history. This README therefore uses **scNanoSeq** for the repository and **Duplex-seq** for the assay/analysis logic.

## 1. Version and provenance

The supplied repository snapshot is based on:

- Repository: `https://github.com/zonglab/scNanoSeq`
- Base branch: `main`
- Base Git commit: `2f1c8cb4f131a57c7c8e86fc1f2d699a968545d7`
- Base commit date: 2025-08-20
- Documentation/demo refresh: 2026-09-16
- The uploaded working tree had already been reorganized into numbered directories (`0_ref`, `1_Bulk_WGS`, ...); that reorganization was not committed in the embedded Git history.

There is no semantic release tag in the supplied Git history, so the commit hash is the authoritative code version for this snapshot. See [`VERSION.md`](VERSION.md) for details.

The preprocessing step refers to NanoSeq `extract_tags.py` from the fixed upstream commit `fd22dec943a3a9cc70643f079c40a8baf13a8e62`. The upstream NanoSeq project documents the `rb`/`mb` tag extraction and use of `bwa mem -C` for carrying those tags into the BAM. The helper is **not vendored here**; use [`scripts/fetch_nanoseq_extract_tags.sh`](scripts/fetch_nanoseq_extract_tags.sh) to retrieve that exact upstream revision.

## 2. Scope of the repository

The repository is organized as follows:

```text
scNanoSeq/
├── 0_ref/                         # small bundled masks / barcode list
│   ├── barcode_3digit.list        # all 64 canonical 3-nt barcodes
│   ├── hg19_cento.bed             # hg19 centromere/telomere exclusion intervals
│   └── hg19_tandem.bed            # hg19 tandem-repeat exclusion intervals
├── 1_Bulk_WGS/                    # matched-bulk mapping and germline-site calling
│   ├── 01_Bulk_mapping.sh
│   ├── 02_Bulk_variant_calling
│   └── 03_Call_summary.sh
├── 2_Somatic_mutation_calling/    # single-cell preprocessing, mapping and calling
│   ├── 01_PreProcessing.sh
│   ├── 02_Mapping.sh
│   ├── 03_Split_bam.sh
│   ├── 04_Variant_call.sh
│   └── 05_Call_summary.sh
├── 2_python_script/               # custom Python helpers / core caller
│   ├── De_novo_estimation_by_variant_position.py
│   ├── DuplexSeq_Mutation_call_aMsN_for_hg19.py
│   ├── bam_split_by_allele_barcode.py
│   ├── de_novo_loci.py
│   ├── extract_tags.py            # fail-fast placeholder; fetched from pinned NanoSeq revision
│   └── localseqfilter_py3.py
├── 3_Phylogeny/
│   └── 01_CellPhy.sh
├── 4_CNV/
│   ├── README.md
│   └── modified_process.R         # modified Ginkgo process.R
├── demo/                          # minimal synthetic no-FASTQ demo fixture
├── scripts/                       # repository setup / audit helpers
├── KNOWN_LIMITATIONS.md
├── RELEASE_CHECKLIST.md
├── LOCAL_UPDATE.md
├── VERSION.md
└── requirements.txt
```

The historical Git tree also contained `Simulation/` and `VAF_fingerprint/` directories. They were already removed from the uploaded working tree and are **not restored here**, because their removal may be intentional and they are not required for the core scNanoSeq workflow described below.

## 3. What the pipeline expects

### 3.1 Reference build

The current scripts are written for **hg19 / GRCh37-style `chr1`...`chr22`, `chrX`, `chrY` contig names**. They refer to the reference as `hg19.fa` and use the bundled hg19 tandem-repeat and centromere/telomere masks.

The full reference FASTA is intentionally not distributed in this repository. Before production use, provide:

```text
hg19.fa
hg19.fa.fai
BWA index files for hg19.fa
```

The reference build and contig naming must match all BAM/VCF/BED inputs.

### 3.2 Matched bulk sample

The somatic workflow uses a matched bulk sample for two purposes:

1. identify high-confidence heterozygous germline sites, which are used to estimate detection sensitivity; and
2. filter single-cell candidate variants against evidence in bulk.

The original code expects the final matched-bulk BAM to be named approximately:

```text
<PATIENT>.RG.markdup.bam
```

and the high-confidence heterozygous SNP set to include:

```text
snp.total.hetero.filtered.vcf
snp.total.hetero.filtered_with_header.vcf
```

### 3.3 Single-cell Duplex-seq BAM tags

The single-cell logic depends on the NanoSeq-style duplex tags:

- `rb`: tag extracted from the current read;
- `mb`: tag extracted from the mate;
- `MR`: a six-base molecule tag created by `bam_split_by_allele_barcode.py` after strand-orientation normalization.

For reverse-orientation flags `83` and `147`, `bam_split_by_allele_barcode.py` swaps `rb` and `mb`, then sets:

```text
MR = mb + rb
```

The resulting reads are split into the 64 possible 3-base `mb` partitions listed in `0_ref/barcode_3digit.list`.

## 4. Complete analysis workflow

### Stage 1. Matched bulk WGS

#### 1.1 Adapter trimming and alignment

`1_Bulk_WGS/01_Bulk_mapping.sh`

The script:

1. trims Nextera adapter sequence with `cutadapt`;
2. aligns paired-end reads to hg19 using `bwa mem -M`;
3. runs `samtools fixmate` and coordinate sorting;
4. removes marked duplicates with `samtools markdup -r`;
5. adds a read group with Picard; and
6. indexes the final BAM.

Main output:

```text
<SAMPLE>.RG.markdup.bam
<SAMPLE>.RG.markdup.bam.bai
```

#### 1.2 Germline variant calling

`1_Bulk_WGS/02_Bulk_variant_calling`

The legacy command calls variants chromosome-by-chromosome with `samtools mpileup` + `bcftools call`. SNPs are retained with high QUAL and alternate-support criteria; indel coordinates are also written separately.

This script uses legacy `samtools mpileup` options. If you move to a modern samtools/bcftools stack, validate the equivalent `bcftools mpileup | bcftools call` command rather than silently assuming CLI compatibility.

#### 1.3 Heterozygous SNP filtering

`1_Bulk_WGS/03_Call_summary.sh`

The script:

1. merges per-chromosome SNP calls;
2. retains `0/1` heterozygous calls;
3. excludes tandem-repeat intervals;
4. excludes centromere/telomere intervals;
5. extracts ±10-bp local sequence; and
6. removes low-complexity local sequence with `localseqfilter_py3.py`.

Key outputs:

```text
snp.total.hetero.filtered.vcf
snp.total.hetero.filtered_with_header.vcf
snp.total.all.vcf
```

### Stage 2. Single-cell Duplex-seq somatic mutation calling

#### 2.1 FASTQ preprocessing

`2_Somatic_mutation_calling/01_PreProcessing.sh`

For the production workflow, paired FASTQs are concatenated and the fixed NanoSeq `extract_tags.py` helper extracts the 3-nt duplex tags. The original parameters are:

```text
-m 3 -s 4 -l 150
```

The resulting tag-bearing reads are adapter-trimmed with `cutadapt`.

**No FASTQ is required for the bundled demo.** The demo starts from a transparent SAM fixture representing the post-alignment/tagged state, then materializes a BAM locally.

#### 2.2 Mapping and optical-duplicate handling

`2_Somatic_mutation_calling/02_Mapping.sh`

The script:

1. maps reads with `bwa mem -C`, which allows the `rb`/`mb` fields from the FASTQ header to become SAM/BAM tags;
2. keeps properly paired reads with MAPQ ≥ 50;
3. coordinate-sorts and indexes the BAM;
4. uses Picard `MarkDuplicates` with `REMOVE_SEQUENCING_DUPLICATES=true` and a large optical-duplicate pixel distance; and
5. indexes the resulting single-cell BAM.

Main output:

```text
<SAMPLE>_DuplexSeq_sorted.bam
<SAMPLE>_DuplexSeq_sorted.bam.bai
```

#### 2.3 Split reads by 3-nt allele barcode

`2_Somatic_mutation_calling/03_Split_bam.sh` calls `2_python_script/bam_split_by_allele_barcode.py`.

For each read, the helper:

- removes the duplicate flag if present in the expected legacy flag representation;
- normalizes `rb`/`mb` orientation for flags 83/147;
- adds the six-base `MR` tag; and
- writes the read into one of 64 barcode-specific BAMs.

Outputs follow:

```text
<SAMPLE>_AAA.bam
<SAMPLE>_AAC.bam
...
<SAMPLE>_TTT.bam
```

#### 2.4 Per-barcode raw candidate discovery and filtering

`2_Somatic_mutation_calling/04_Variant_call.sh`

For each barcode BAM, the script:

1. runs `samtools mpileup | bcftools call`;
2. removes calls in tandem-repeat intervals;
3. removes calls in centromere/telomere intervals;
4. removes non-SNV calls and low-complexity local sequence with `localseqfilter_py3.py`; and
5. passes the remaining SNVs to `DuplexSeq_Mutation_call_aMsN_for_hg19.py`.

Although `bcftools` can emit indels, the present custom caller explicitly requires `len(REF) == len(ALT) == 1`. Therefore the **custom somatic scoring implemented here is SNV-only**.

#### 2.5 Duplex SNV scoring

`2_python_script/DuplexSeq_Mutation_call_aMsN_for_hg19.py`

For each candidate SNV, the caller:

- processes only canonical hg19 chromosomes (`chr1`-`chr22`, `chrX`, `chrY`);
- fetches overlapping reads from the barcode BAM;
- chooses the more represented read-orientation group (`83/147` versus `99/163`);
- groups reads by mapping coordinates plus `MR` molecule tag;
- requires base quality ≥ 30;
- rejects reads whose CIGAR contains an insertion or deletion;
- rejects reads with >3 lowercase bases in `get_reference_sequence()` (the code uses lowercase mismatch representation as an alignment-quality proxy);
- requires the candidate to lie >8 bp from both molecule ends; and
- classifies the highest-coverage molecule group using total-allele and per-strand read thresholds.

The production call uses:

```text
allele_read threshold = 4
strand_read threshold = 2
```

The principal output categories are:

| Category | Code interpretation |
| --- | --- |
| `Mutation` | alternate allele supported on both strand directions with no reference-supporting reads |
| `Single_strand_variant_allele_filtered` | alternate-only support but only one strand direction represented |
| `Damage_0` | one transformed orientation dominates with sufficient reads |
| `Damage_swap_strand_filtered` | swap-like pattern with strand support |
| other damage/swap labels | diagnostic intermediate classes not propagated as final mutations |

The output line contains the original VCF record plus:

```text
(Fwt, Fmut, Rwt, Rmut)
(position_from_proximal_end, position_from_distal_end)
MR
```

#### 2.6 Matched-bulk filtering

`04_Variant_call.sh` next compares candidate classes against matched bulk evidence. This step calls:

```text
gtbulksearch_py3.py
```

**That helper is missing from both the uploaded snapshot and all reachable Git objects in the embedded repository history.** Its scientific behavior cannot be reconstructed reliably from the filename and command line alone, so it has not been fabricated in this refresh. Production execution of Stage 2.6 remains blocked until the original helper is recovered. See [`KNOWN_LIMITATIONS.md`](KNOWN_LIMITATIONS.md).

#### 2.7 Merge calls, dbSNP filtering, and mutation-burden estimation

`2_Somatic_mutation_calling/05_Call_summary.sh`

After all 64 barcode jobs produce completion flags, the script:

1. concatenates per-barcode results;
2. sorts mutation and single-strand variant VCFs;
3. annotates dbSNP using SnpSift;
4. removes records carrying `rs` identifiers; and
5. estimates the de novo mutation burden as a function of distance from molecule ends.

`De_novo_estimation_by_variant_position.py` computes cumulative counts for each minimum end distance `p`. In the code, the uncorrected estimate is:

```text
N_de_novo(p) = C_de_novo(p) / C_bulk_het_detected(p) × N_bulk_heterozygous
```

The strand-swapping rate is estimated as:

```text
r_swap = N_Damage_swap / N_Damage_0
```

and the swap-corrected estimate is:

```text
N_corrected(p) =
    [C_de_novo(p) - r_swap × C_single_strand_de_novo(p)]
    ---------------------------------------------------------------- × N_bulk_heterozygous
    [C_bulk_het_mut(p) - r_swap × C_bulk_het_single_strand(p)]
```

Main output:

```text
<SAMPLE>_a4s2_de_novo_estimation.txt
```

### Stage 3. Multi-cell phylogeny with CellPhy

`3_Phylogeny/01_CellPhy.sh`

The script:

1. lists all single cells for a patient;
2. restores VCF headers on each final mutation file;
3. uses `de_novo_loci.py` to create the union of observed de novo loci;
4. jointly genotypes those loci across matched bulk + all single-cell BAMs using `samtools mpileup | bcftools call`;
5. retains SNPs; and
6. runs `cellphy.sh -t 8` on the combined multi-sample VCF.

CellPhy is an external program and is not vendored here.

### Stage 4. CNV analysis with modified Ginkgo

`4_CNV/modified_process.R` is a modified form of Ginkgo's CNV processing step. The accompanying `4_CNV/README.md` describes installation into a standalone Ginkgo workspace.

The script expects the Ginkgo genome-bin files and a binned read-count matrix produced by the Ginkgo workflow. It performs GC correction, segmentation, ploidy/copy-number inference, and profile plotting. The modified script writes ordered outputs including:

```text
zout_01_summary_stats.tsv
zout_02_normalized_rr.tsv
zout_03_break__points.tsv
zout_04_smoothened_rr.tsv
zout_05_CN__local_min.tsv
zout_06_CN__close_min.tsv
```

This repository does not reimplement Ginkgo's BAM-to-bin preprocessing. Ginkgo should be installed separately and `modified_process.R` used as documented in `4_CNV/README.md`.

## 5. Software requirements

### Custom Python code

The current Python scripts require Python 3 and the packages listed in [`requirements.txt`](requirements.txt):

```text
numpy
pandas
pysam
```

The exact historical package versions were not recorded in the repository and are therefore not invented here.

### Command-line software used by the shell pipeline

The scripts reference:

- `cutadapt`
- `bwa`
- `samtools`
- `bcftools`
- `bedtools`
- `seqtk`
- Picard (`picard` command)
- Java + `SnpSift.jar`
- CellPhy (`cellphy.sh`) for phylogeny
- Ginkgo + R for CNV analysis
- SLURM (`sbatch`) in the provided job templates

The bulk and single-cell calling scripts use legacy samtools/bcftools syntax. For strict reproducibility, record the versions from the environment in which the original analysis was run. You can capture the currently installed versions with:

```bash
bash scripts/capture_versions.sh > software_versions.txt
```

## 6. Required external resources

The following large or third-party resources are intentionally not bundled:

```text
hg19.fa and indices
SnpSift.jar
dbSNP/All_20151104.vcf
CellPhy installation
Ginkgo installation and genome-bin resources
NanoSeq extract_tags.py (retrievable from a pinned upstream commit)
```

To retrieve the exact `extract_tags.py` revision referenced by the original repository:

```bash
bash scripts/fetch_nanoseq_extract_tags.sh
```

## 7. Minimal demo: what should be distributed instead of FASTQ?

For this repository, shipping raw FASTQ would be large, unnecessary, and mostly demonstrate mature external tools (`cutadapt` and `bwa`) rather than the custom method. The most informative minimal test should begin at the first method-specific intermediate.

The bundled [`demo/`](demo/) therefore contains a **synthetic, non-biological** fixture with:

1. a tiny coordinate-sorted-compatible SAM containing paired reads with `rb`, `mb`, `MD`, and alignment fields needed by the custom logic;
2. one SNV candidate VCF;
3. a tiny matching reference FASTA for transparency; and
4. an expected `Mutation` classification.

The demo tests the two repository-specific operations that matter most:

```text
barcode/orientation normalization + MR creation
        ↓
DuplexSeq_Mutation_call_aMsN_for_hg19.py scoring
```

It deliberately does **not** attempt to prove that `bwa`, Picard, CellPhy, or Ginkgo themselves work. Those are external packages with their own test suites.

Run:

```bash
bash demo/run_core_demo.sh
```

The demo requires `samtools`, Python 3, and `pysam`. It converts the transparent SAM fixture to a sorted/indexed BAM locally, so no FASTQ and no binary BAM need to be stored in Git.

See [`demo/README.md`](demo/README.md) for the exact fixture contract and what a larger publication demo would add.

## 8. Configuring the production scripts

The shell scripts are **analysis templates**, not a portable workflow manager. They contain placeholders such as:

```text
xxxxxx/
xxxxxxxx/
AYY
BYY
AXX
AXXX
```

Before production use, replace these with your actual project paths/sample identifiers or refactor them into your site's configuration system.

At minimum configure:

- reference FASTA path;
- bulk data/output root;
- single-cell data/output root;
- script/helper paths;
- dbSNP and SnpSift paths;
- CellPhy path;
- Ginkgo path; and
- SLURM resources appropriate to your cluster.

The existing job chaining uses `sed` to substitute sample/barcode placeholders and submits downstream jobs with `sbatch`. Preserve that behavior only if it matches your scheduler layout.

## 9. Known limitations before a public release

The most important unresolved items are summarized in [`KNOWN_LIMITATIONS.md`](KNOWN_LIMITATIONS.md). In particular:

- `gtbulksearch_py3.py` is missing and blocks faithful matched-bulk filtering;
- the exact original software environment is not recorded;
- the repository does not currently contain an explicit top-level software license;
- author/citation metadata cannot be filled in without project-owner information;
- several absolute/path placeholders remain in production shell templates;
- CNV execution depends on an external Ginkgo installation; and
- the core custom mutation scorer is SNV-only.

These are documented rather than silently guessed.

## 10. Repository validation

Basic static checks can be run with:

```bash
bash scripts/repo_check.sh
```

This checks shell syntax, Python syntax, expected repository files, unresolved path placeholders, and whether key external executables are available. It does not claim scientific validation of external tools or the missing matched-bulk helper.

## 11. Recommended release contents

A practical public release should contain:

```text
source scripts
small hg19 masks/barcode list already present
top-level README
minimal synthetic no-FASTQ demo
expected demo output
software requirements / version provenance
top-level LICENSE chosen by the project owner
CITATION.cff with correct authors and publication metadata
the original gtbulksearch_py3.py (or a documented replacement validated against the original analysis)
```

Do **not** include raw human FASTQ, production BAMs, large reference genomes, dbSNP databases, or patient-derived example files merely to make the repository look substantial. A tiny synthetic fixture that exercises the custom logic is more useful and vastly less troublesome.

## 12. References and third-party software

Please cite the underlying methods/software used in your analysis as appropriate, including NanoSeq/Duplex Sequencing, CellPhy, Ginkgo, BWA, samtools/bcftools, Picard, bedtools, cutadapt, and SnpSift.

The NanoSeq helper referenced by this repository is maintained by the `cancerit/NanoSeq` project and is licensed separately by its upstream authors. CellPhy and Ginkgo are also separate projects with their own licenses and citation requirements.
