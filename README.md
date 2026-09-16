# scNanoSeq

**Single-cell Duplex-seq somatic mutation calling, mutation-burden estimation, phylogeny reconstruction, and CNV analysis.**

This repository contains the analysis scripts used to process matched bulk WGS and single-cell Duplex-seq data, classify duplex-supported SNVs, estimate de novo mutation burden, reconstruct multi-cell phylogenies with CellPhy, and perform CNV analysis with a modified Ginkgo workflow.

## 1. Repository structure

```text
scNanoSeq/
├── 0_ref/                         # bundled masks and barcode list
│   ├── barcode_3digit.list        # 64 canonical 3-nt barcodes
│   ├── hg19_cento.bed             # hg19 centromere/telomere exclusion intervals
│   └── hg19_tandem.bed            # hg19 tandem-repeat exclusion intervals
├── 1_Bulk_WGS/                    # matched-bulk mapping and germline-site calling
│   ├── 01_Bulk_mapping.sh
│   ├── 02_Bulk_variant_calling
│   └── 03_Call_summary.sh
├── 2_Somatic_mutation_calling/    # single-cell preprocessing, mapping, and calling
│   ├── 01_PreProcessing.sh
│   ├── 02_Mapping.sh
│   ├── 03_Split_bam.sh
│   ├── 04_Variant_call.sh
│   └── 05_Call_summary.sh
├── 2_python_script/               # custom Python helpers and core caller
│   ├── De_novo_estimation_by_variant_position.py
│   ├── DuplexSeq_Mutation_call_aMsN_for_hg19.py
│   ├── bam_split_by_allele_barcode.py
│   ├── de_novo_loci.py
│   ├── extract_tags.py            # placeholder; fetched from pinned NanoSeq revision
│   └── localseqfilter_py3.py
├── 3_Phylogeny/
│   └── 01_CellPhy.sh
├── 4_CNV/
│   ├── README.md
│   └── modified_process.R         # modified Ginkgo process.R
├── demo/                          # chr20 MCF10A tagged BAM demo data
├── RELEASE_CHECKLIST.md
└── requirements.txt
```

## 2. Software and external resources

### Recorded core environment

The original analysis used:

```text
cutadapt  5.0
bwa       0.7.13-r1126
samtools  1.12
bcftools  1.12
Picard    3.4.0.0
bedtools  2.31.1
seqtk     1.4
SnpSift   4.1k
Python    3.10
R         4.4.3
numpy     2.2.6
pandas    2.2.3
pysam     0.23.3
```

The shell workflows also use CellPhy and Ginkgo. 

### External resources not bundled

The following large or third-party resources must be provided separately:

```text
hg19.fa and indices
SnpSift.jar
dbSNP/All_20151104.vcf for hg19
CellPhy installation
Ginkgo installation and genome-bin resources
NanoSeq extract_tags.py
```

The pinned NanoSeq helper can be retrieved with:

```bash
bash 2_python_script/fetch_nanoseq_extract_tags.sh
```

## 3. Input conventions

### 3.1 Reference genome

The workflow is written for **hg19 / GRCh37-style `chr1`...`chr22`, `chrX`, `chrY` contig names** and uses the bundled hg19 tandem-repeat and centromere/telomere masks.

Required reference files:

```text
hg19.fa
hg19.fa.fai
BWA index files for hg19.fa
```

The reference build and contig naming must match all BAM, VCF, and BED inputs.

### 3.2 Matched bulk sample

The matched bulk sample is used to:

1. identify high-confidence heterozygous germline sites for detection-sensitivity estimation; and
2. filter candidate somatic variants against bulk evidence.

Expected outputs from the bulk workflow include:

```text
<PATIENT>.RG.markdup.bam
snp.total.hetero.filtered.vcf
snp.total.hetero.filtered_with_header.vcf
```

### 3.3 Single-cell Duplex-seq BAM tags

The single-cell workflow uses NanoSeq-style duplex tags:

- `rb`: tag extracted from the current read;
- `mb`: tag extracted from the mate;
- `MR`: six-base molecule tag created by `bam_split_by_allele_barcode.py` after strand-orientation normalization.

For flags `83` and `147`, `rb` and `mb` are swapped before setting:

```text
MR = mb + rb
```

Reads are then split into the 64 possible 3-base `mb` partitions listed in `0_ref/barcode_3digit.list`.

## 4. Analysis workflow

### Stage 1. Matched bulk WGS

#### 1.1 Mapping

**Script:** `1_Bulk_WGS/01_Bulk_mapping.sh`

The workflow trims Nextera adapters, aligns paired-end reads with `bwa mem -M`, runs `samtools fixmate` and coordinate sorting, removes duplicates with `samtools markdup -r`, adds read groups with Picard, and indexes the final BAM.

Output:

```text
<SAMPLE>.RG.markdup.bam
<SAMPLE>.RG.markdup.bam.bai
```

#### 1.2 Germline variant calling

**Script:** `1_Bulk_WGS/02_Bulk_variant_calling`

Variants are called chromosome-by-chromosome with `samtools mpileup` + `bcftools call`. High-quality SNPs are retained and indel coordinates are written separately.

This step uses legacy `samtools mpileup` options. When using a newer samtools/bcftools stack, validate an equivalent `bcftools mpileup | bcftools call` workflow.

#### 1.3 Heterozygous SNP filtering

**Script:** `1_Bulk_WGS/03_Call_summary.sh`

Per-chromosome calls are merged, `0/1` heterozygous SNPs are retained, tandem-repeat and centromere/telomere regions are excluded, and low-complexity sequence is removed using `localseqfilter_py3.py`.

Key outputs:

```text
snp.total.hetero.filtered.vcf
snp.total.hetero.filtered_with_header.vcf
snp.total.all.vcf
```

### Stage 2. Single-cell Duplex-seq somatic mutation calling

#### 2.1 FASTQ preprocessing

**Script:** `2_Somatic_mutation_calling/01_PreProcessing.sh`

Paired FASTQs are concatenated, NanoSeq `extract_tags.py` extracts the 3-nt duplex tags, and the resulting reads are adapter-trimmed with `cutadapt`.

Original tag-extraction parameters:

```text
-m 3 -s 4 -l 150
```

#### 2.2 Mapping and optical-duplicate handling

**Script:** `2_Somatic_mutation_calling/02_Mapping.sh`

Reads are mapped with `bwa mem -C`, preserving `rb`/`mb` header fields as BAM tags. Properly paired reads with MAPQ ≥ 50 are retained, coordinate-sorted, and processed with Picard `MarkDuplicates` using `REMOVE_SEQUENCING_DUPLICATES=true` and a large optical-duplicate pixel distance.

Output:

```text
<SAMPLE>_DuplexSeq_sorted.bam
<SAMPLE>_DuplexSeq_sorted.bam.bai
```

#### 2.3 Split reads by 3-nt barcode

**Script:** `2_Somatic_mutation_calling/03_Split_bam.sh`  
**Helper:** `2_python_script/bam_split_by_allele_barcode.py`

The helper normalizes `rb`/`mb` orientation, adds the six-base `MR` tag, and writes reads into one of 64 barcode-specific BAMs.

```text
<SAMPLE>_AAA.bam
<SAMPLE>_AAC.bam
...
<SAMPLE>_TTT.bam
```

#### 2.4 Candidate discovery and filtering

**Script:** `2_Somatic_mutation_calling/04_Variant_call.sh`

For each barcode BAM, the workflow:

1. calls variants with `samtools mpileup | bcftools call`;
2. excludes tandem-repeat and centromere/telomere intervals;
3. removes non-SNV calls and low-complexity local sequence; and
4. passes the remaining SNVs to `DuplexSeq_Mutation_call_aMsN_for_hg19.py`.

#### 2.5 Duplex SNV scoring

**Script:** `2_python_script/DuplexSeq_Mutation_call_aMsN_for_hg19.py`

For each candidate SNV, the caller:

- processes canonical hg19 chromosomes only;
- fetches overlapping reads from the barcode BAM;
- selects the more represented orientation group (`83/147` or `99/163`);
- groups reads by mapping coordinates and `MR` tag;
- requires base quality ≥ 30;
- rejects reads containing CIGAR insertions or deletions;
- rejects reads with >3 lowercase bases in `get_reference_sequence()`;
- requires the candidate to lie >8 bp from both molecule ends; and
- classifies the highest-coverage molecule group using total-allele and per-strand read thresholds.

Production thresholds:

```text
allele_read threshold = 4
strand_read threshold = 2
```

Principal output categories:

| Category | Interpretation |
| --- | --- |
| `Mutation` | alternate allele supported on both strand directions with no reference-supporting reads |
| `Single_strand_variant_allele_filtered` | alternate-only support but only one strand direction represented |
| `Damage_0` | one transformed orientation dominates with sufficient reads |
| `Damage_swap_strand_filtered` | swap-like pattern with strand support |
| other damage/swap labels | diagnostic intermediate classes not propagated as final mutations |

Each output record contains the original VCF fields plus:

```text
(Fwt, Fmut, Rwt, Rmut)
(position_from_proximal_end, position_from_distal_end)
MR
```

#### 2.6 Matched-bulk filtering

`04_Variant_call.sh` compares candidate classes against matched bulk evidence using:

```text
gtbulksearch_py3.py
```

#### 2.7 Final call set and mutation-burden estimation

**Script:** `2_Somatic_mutation_calling/05_Call_summary.sh`

After all 64 barcode jobs finish, the workflow merges the results, sorts mutation and single-strand VCFs, annotates dbSNP with SnpSift, removes records carrying `rs` identifiers, and estimates de novo mutation burden as a function of distance from molecule ends.

`De_novo_estimation_by_variant_position.py` computes cumulative counts for each minimum end distance `p`.

Uncorrected estimate:

```text
N_de_novo(p) = C_de_novo(p) / C_bulk_het_detected(p) × N_bulk_heterozygous
```

Strand-swapping rate:

```text
r_swap = N_Damage_swap / N_Damage_0
```

Swap-corrected estimate:

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

### Stage 3. Multi-cell phylogeny

**Script:** `3_Phylogeny/01_CellPhy.sh`

For each patient, the workflow combines de novo loci observed across cells, jointly genotypes those loci in matched bulk and all single-cell BAMs, retains SNPs, and runs:

```bash
cellphy.sh -t 8
```

CellPhy is installed separately and is not vendored in this repository.

### Stage 4. CNV analysis

**Script:** `4_CNV/modified_process.R`

This is a modified Ginkgo CNV-processing step. It expects Ginkgo genome-bin files and a binned read-count matrix, then performs GC correction, segmentation, ploidy/copy-number inference, and profile plotting.

Representative outputs:

```text
zout_01_summary_stats.tsv
zout_02_normalized_rr.tsv
zout_03_break__points.tsv
zout_04_smoothened_rr.tsv
zout_05_CN__local_min.tsv
zout_06_CN__close_min.tsv
```

Ginkgo's BAM-to-bin preprocessing is not reimplemented here. See `4_CNV/README.md` for setup details.

## 5. Third-party software and citation

Please cite the underlying methods and software as appropriate, including NanoSeq/Duplex Sequencing, CellPhy, Ginkgo, BWA, samtools/bcftools, Picard, bedtools, cutadapt, and SnpSift.

The NanoSeq helper used by this repository is maintained by the `cancerit/NanoSeq` project and is distributed under its upstream license. CellPhy and Ginkgo are separate projects with their own licenses and citation requirements.
