# Author action required before claiming checklist compliance

These items depend on facts, files, legal choices, or validation that cannot be reconstructed responsibly from the repository alone.

1. **Add `gtbulksearch_py3.py` or identify its exact public source.** It is called by the central mutation workflow but is absent.
2. **Record exact tested versions and OS/hardware.** Run `scripts/collect_environment.sh` in the actual analysis environment(s).
3. **State typical install time.** Measure a clean install/recreation of the documented environment.
4. **Choose and add the repository license.** This is a legal/author/institutional decision. Do not infer one from third-party dependencies.
5. **Create a manuscript release/tag and DOI/permanent identifier.** The current public GitHub repository has no release/DOI recorded in the supplied material.
6. **Provide a biologically representative end-to-end test dataset if the full FASTQ/BAM-to-call workflow is central to the paper.** The included demo is a deterministic code smoke test, not a biological validation dataset.
7. **Record expected runtime for the representative end-to-end demo** on a stated machine.
8. **Document reference/database provenance and checksums** for hg19, dbSNP, the tandem-repeat BED, and centromere BED where redistribution permits.
9. **Record CellPhy version and expected demo tree output** if CellPhy-derived lineage results are a central paper result.
10. **Record Ginkgo/R versions and upstream base revision** for `modified_process.R` if the CNV analysis is part of the reported results.
11. **Identify the manuscript location of the detailed code/pseudocode description** (Methods/main text/elsewhere).
12. **Have an unfamiliar colleague run the installation/demo instructions** and incorporate their feedback, as Nature recommends.
13. **Confirm whether deletion of the historical `Simulation/` and `VAF_fingerprint/` directories is intentional.** They exist on the current public GitHub `main` but are absent from the uploaded working tree used to prepare this bundle.
