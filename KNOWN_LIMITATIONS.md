# Known limitations and unresolved provenance

This file intentionally separates verified repository facts from items that cannot be reconstructed safely.

## Release-blocking scientific/reproducibility items

1. **`gtbulksearch_py3.py` is missing.**
   - It is called four times by `2_Somatic_mutation_calling/04_Variant_call.sh`.
   - It is absent from the uploaded working tree.
   - It is absent from all reachable Git objects/history in the supplied repository.
   - Its exact behavior cannot be inferred safely from the filename and three command-line arguments.
   - Do not replace it with a guessed pileup filter without validating against the original analysis/output.

2. **The exact original software versions are not recorded.**
   - The scripts use legacy `samtools mpileup` syntax.
   - Reproducibility therefore requires either recovering the original environment/module list or validating equivalent modern commands.

3. **The top-level project license is missing.**
   - A license cannot be selected on behalf of the project owner.
   - The upstream NanoSeq helper has its own AGPL-3.0-or-later licensing terms.

4. **Author/citation metadata are missing.**
   - `CITATION.cff` should not be fabricated without the correct project authors, preferred citation, DOI/preprint/publication status, and contact information.

## Pipeline portability limitations

- Production shell scripts contain site-specific placeholders (`xxxxxx`, `xxxxxxxx`) and sample placeholders (`AYY`, `BYY`, `AXX`, `AXXX`).
- Job chaining assumes SLURM and uses `sbatch`.
- Helper script calls are written as local filenames and therefore assume those helpers have been copied/symlinked into the run directory or the scripts have been edited to absolute paths.
- The full hg19 FASTA/index, dbSNP VCF, SnpSift JAR, CellPhy, and Ginkgo resources are external.
- The core `DuplexSeq_Mutation_call_aMsN_for_hg19.py` scorer is SNV-only.
- The caller only accepts canonical `chr1`-`chr22`, `chrX`, and `chrY` names.
- `De_novo_estimation_by_variant_position.py` assumes non-empty input classes; empty classes can fail before output is produced.
- The original scripts use broad `except:` clauses in several places, which can hide malformed reads/records.

## Demo limitation

The bundled demo is a synthetic unit/integration fixture for the custom barcode split + duplex SNV scoring path. It is not a biological benchmark and does not exercise the missing `gtbulksearch_py3.py`, SnpSift/dbSNP filtering, CellPhy, or Ginkgo.
