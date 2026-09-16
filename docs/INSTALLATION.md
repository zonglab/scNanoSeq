# Installation and system requirements

## Intended environment

The supplied workflow scripts were written for a Unix-like SLURM/HPC environment. They use Bash, standard Unix utilities, and a number of external bioinformatics programs. The repository does not currently provide a fully pinned container/Conda environment because the exact versions used for the manuscript have not been supplied and should not be guessed.

### AUTHOR ACTION REQUIRED: tested environment

Run the following command in the actual environment used to generate the manuscript results:

```bash
bash scripts/collect_environment.sh > docs/TESTED_ENVIRONMENT.txt
```

Review the output and commit it. If the analysis was run on more than one environment, capture each one separately (for example `TESTED_ENVIRONMENT_bulk.txt` and `TESTED_ENVIRONMENT_cellphy.txt`).

## Required software by workflow stage

| Component | Used for | Version status |
|---|---|---|
| Bash | workflow scripts | AUTHOR ACTION REQUIRED: record actual version |
| SLURM / `sbatch` | job submission | AUTHOR ACTION REQUIRED: record cluster/SLURM version if relevant |
| BWA | read mapping | AUTHOR ACTION REQUIRED: record actual version |
| SAMtools | BAM operations / pileup | AUTHOR ACTION REQUIRED: record actual version |
| BCFtools | variant calling/filtering | AUTHOR ACTION REQUIRED: record actual version; command syntax is version-sensitive |
| BEDTools | interval exclusion/intersection | AUTHOR ACTION REQUIRED |
| Seqtk | sequence extraction | AUTHOR ACTION REQUIRED |
| Cutadapt | adapter trimming | AUTHOR ACTION REQUIRED |
| Picard | read groups / duplicate processing | AUTHOR ACTION REQUIRED |
| Java | SnpSift/Picard runtime | AUTHOR ACTION REQUIRED |
| SnpSift | dbSNP annotation | AUTHOR ACTION REQUIRED; record jar/version and dbSNP resource provenance |
| Python 3 | custom scripts | some shebangs mention Python 3.5/3.6 historically; record the version actually used |
| Pysam | BAM access | required by `bam_split_by_allele_barcode.py` and `DuplexSeq_Mutation_call_aMsN_for_hg19.py`; record actual version |
| NumPy | de-novo estimation | record actual version |
| NanoSeq `extract_tags.py` | barcode extraction | upstream source is pinned to commit `fd22dec943a3a9cc70643f079c40a8baf13a8e62`; see `THIRD_PARTY.md` |
| CellPhy | phylogeny | AUTHOR ACTION REQUIRED: record exact version/build and citation |
| Ginkgo + R | CNV stage | AUTHOR ACTION REQUIRED: record exact versions and any modified Ginkgo base revision |

## Required reference/resources

The repository does **not** bundle all large or redistributability-sensitive resources. The following must be obtained separately and documented with build/version/checksum where feasible:

- hg19 reference FASTA (`hg19.fa`) and index files required by BWA/SAMtools;
- dbSNP VCF used by SnpSift (`dbSNP/All_20151104.vcf` in the historical script);
- `SnpSift.jar`;
- any CellPhy reference/installation assets;
- the Ginkgo installation used for CNV analysis.

The repository does include `0_ref/hg19_tandem.bed`, `0_ref/hg19_cento.bed`, and `0_ref/barcode_3digit.list`. The provenance and any redistribution rights for the BED resources should be documented by the authors before publication.

## Install-time statement

A reliable "typical install time" cannot be reconstructed from the repository alone because the original environment and installation method are not recorded. **AUTHOR ACTION REQUIRED:** install the documented environment on a clean machine/container, time it, and report the result in this file and the main README.

The included Python-only reviewer smoke demo requires only Python 3 plus NumPy and should take seconds to run once those are available.
