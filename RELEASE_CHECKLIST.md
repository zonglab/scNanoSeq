# scNanoSeq release checklist

Status legend:

- `[x]` completed in this refresh
- `[ ]` action still required
- `[!]` cannot be completed faithfully without project-owner/original-analysis information

## Documentation and repository structure

- [x] Replace the two-line top-level README with a complete pipeline description.
- [x] Record the exact base Git commit and repository provenance.
- [x] Explain the numbered directory structure in the uploaded working tree.
- [x] Document the matched-bulk, single-cell, phylogeny, and CNV stages.
- [x] Document the custom duplex SNV classifier logic and `a4s2` thresholds.
- [x] Document the de novo burden / strand-swap correction equations implemented in code.
- [x] Document external resources that are intentionally not bundled.
- [x] Add `.gitignore` entries for notebook checkpoints, Python cache, SLURM logs, and demo work products.
- [x] Add a static repository check script.

## Minimal demonstration data

- [x] Do not include FASTQ.
- [x] Add a small, synthetic, human-readable SAM fixture with required `rb`/`mb`/`MD` fields.
- [x] Add a one-site candidate VCF and tiny reference FASTA.
- [x] Add expected custom caller output.
- [x] Add a demo runner that materializes a BAM locally and exercises barcode splitting + duplex scoring.
- [x] Clearly label all demo data as synthetic/non-biological.
- [x] Keep CellPhy/Ginkgo out of the mandatory minimal demo; they are external downstream tools rather than the custom core method.

## Dependencies and third-party code

- [x] Add Python package requirements (`numpy`, `pandas`, `pysam`).
- [x] Add a version-capture helper rather than inventing historical dependency versions.
- [x] Pin the NanoSeq `extract_tags.py` source to the exact upstream commit referenced by the original repository.
- [x] Make the local `extract_tags.py` placeholder fail loudly instead of silently doing nothing.
- [x] Update the preprocessing template from `python2` to `python3`, consistent with the pinned upstream helper.
- [x] Add a helper script to fetch that exact upstream file.
- [!] Recover `gtbulksearch_py3.py` from the original analysis environment, author archive, or lab storage. **Required for faithful Stage 2 matched-bulk filtering.**
- [!] Recover the original versions/modules for samtools, bcftools, BWA, Picard, cutadapt, bedtools, seqtk, SnpSift, Python/pysam, CellPhy, R, and Ginkgo if strict computational reproducibility is required.

## Scientific resources

- [x] Retain the small bundled hg19 tandem-repeat and centromere/telomere masks.
- [ ] Verify and document the exact source/version of `hg19_tandem.bed` and `hg19_cento.bed`.
- [!] Provide the exact hg19 FASTA accession/source used in the published/original analysis.
- [!] Provide the exact dbSNP resource provenance for `dbSNP/All_20151104.vcf`.

## Licensing and citation

- [!] Add a top-level `LICENSE` selected/approved by the project owner.
- [!] Add `CITATION.cff` only after confirming authors, preferred citation, DOI/preprint/publication, and contact metadata.
- [ ] Confirm that any vendored third-party code complies with its upstream license; the current refresh fetches NanoSeq code from upstream instead of redistributing it.
- [ ] Add acknowledgements/citations for NanoSeq, CellPhy, Ginkgo, and other external software in the manuscript/repository release.

## Production script cleanup

- [ ] Replace all `xxxxxx` / `xxxxxxxx` path placeholders with a documented configuration strategy before claiming turnkey execution.
- [ ] Decide whether SLURM job submission should remain part of the public scripts or be converted into scheduler-neutral commands/wrappers.
- [ ] Replace relative helper-script assumptions with robust repository-root paths or an installation step.
- [x] Fix the obvious bulk summary mask filename mismatch from `cento.bed` to the bundled `hg19_cento.bed`.
- [ ] Validate legacy `samtools mpileup` commands under the intended software version; modern samtools may require `bcftools mpileup` equivalents.
- [ ] Add at least one small validated matched-bulk filtering fixture after `gtbulksearch_py3.py` is recovered.

## Final pre-release verification

- [ ] Run `bash scripts/repo_check.sh` in the intended analysis environment.
- [ ] Run `bash demo/run_core_demo.sh` and confirm it reports `PASS`.
- [ ] Run one real, de-identified/controlled sample through the complete production workflow and compare key outputs with the original analysis.
- [ ] Create a clean release commit and tag after all release-blocking items are resolved.
