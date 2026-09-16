# Running the historical workflow on your data

The shell scripts in this repository are analysis scripts/templates rather than a fully parameterized workflow engine. To preserve provenance, study-specific placeholder strings have not been silently rewritten. A user must therefore configure paths and sample IDs before execution.

## 1. Prepare external resources

Provide:

- `hg19.fa` plus BWA/SAMtools index files;
- `hg19_tandem.bed` and `hg19_cento.bed` (copies are distributed under `0_ref/`);
- `barcode_3digit.list` (`0_ref/`);
- `SnpSift.jar` and the dbSNP VCF used by the analysis;
- all programs listed in `INSTALLATION.md` on `PATH`;
- the original `gtbulksearch_py3.py` (currently missing from this repository).

## 2. Bulk WGS

Scripts: `1_Bulk_WGS/01_Bulk_mapping.sh`, `02_Bulk_variant_calling`, `03_Call_summary.sh`.

Historical inputs include paired-end bulk FASTQ files and the hg19 reference. `01_Bulk_mapping.sh` trims Nextera adapters, maps reads, removes duplicates, and adds read groups. The variant-calling step expects shell variables including `sample`, `data_dir`, and `chr`; it is therefore not a standalone command until those variables are defined. The summary step combines chromosome calls and removes tandem-repeat/centromeric/local-sequence artifacts.

Expected major outputs include:

- `<sample>.RG.markdup.bam` and BAM index;
- `output.<chr>.bcf`;
- `snp.total.hetero.filtered_with_header.vcf`;
- `snp.total.all.vcf`.

## 3. Single-cell Duplex-seq somatic mutation calling

Scripts are intended to run in numeric order:

1. `01_PreProcessing.sh`: concatenate FASTQ lanes, extract molecular/barcode tags, trim adapters.
2. `02_Mapping.sh`: BWA mapping, mapping-quality filtering, coordinate sort, optical-duplicate handling.
3. `03_Split_bam.sh`: split BAM into the 64 three-base barcode groups.
4. `04_Variant_call.sh`: per-barcode candidate calling, sequence-context filters, Duplex-seq mutation classification, matched-bulk filtering.
5. `05_Call_summary.sh`: wait for per-barcode jobs, merge results, dbSNP filtering, estimate de-novo burden as a function of read position.

Study-specific placeholders that must be replaced include:

- `xxxxxx` / `xxxxxxxx`: project/data/script roots;
- `AYY`, `AXX`, `AXXX`: sample or patient IDs;
- `BYY`: three-base barcode;
- `AAA`: a bulk VCF-header source placeholder in the CellPhy script.

Do not globally replace these strings without reviewing each occurrence: the historical scripts use different placeholder lengths in different directories.

Expected major per-cell outputs include:

- `<sample>_DuplexSeq_sorted.bam`;
- `<sample>_a4s2_mutation_dbSNP_filtered.vcf`;
- `<sample>_a4s2_ss_var_dbSNP_filtered.vcf`;
- `<sample>_a4s2_de_novo_estimation.txt`.

### Critical missing dependency

`2_Somatic_mutation_calling/04_Variant_call.sh` calls `gtbulksearch_py3.py` four times, but that file is absent from the repository and no public provenance is recorded. The full pipeline is therefore **not reproducible yet**. Add the original script and document its origin/license before claiming end-to-end reproducibility.

## 4. Phylogeny with CellPhy

`3_Phylogeny/01_CellPhy.sh` builds a list of bulk/single-cell BAMs, identifies de-novo loci, creates a combined multi-sample VCF with SAMtools/BCFtools, restricts to SNPs, and invokes:

```bash
cellphy.sh -t 8 <combined_SNP.vcf>
```

CellPhy is a third-party tool. You do not need to re-demonstrate the correctness of CellPhy itself. For manuscript reproducibility, record the exact CellPhy version, cite it, preserve the exact command/parameters, and ideally provide the small combined VCF used for a toy/reviewer example plus the expected tree files.

## 5. CNV with Ginkgo

See `4_CNV/README.md`. This stage depends on a separately installed Ginkgo environment and a modified processing script. Record the exact Ginkgo revision/version, R version, and the relationship between `modified_process.R` and the upstream file used as its base.
