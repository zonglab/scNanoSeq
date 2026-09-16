# Nature-style Code and Software Submission Checklist: scNanoSeq status

Status legend:

- **COMPLETED HERE**: added/verified in this review-preparation bundle.
- **PARTIAL**: documentation/support exists, but an author-supplied fact/file is still needed.
- **AUTHOR ACTION REQUIRED**: cannot be truthfully invented or recovered from the supplied repository.

| Checklist requirement | Status | Evidence / remaining action |
|---|---|---|
| Source code supplied | **PARTIAL** | Main shell/Python/R scripts are present. `gtbulksearch_py3.py` is called by `04_Variant_call.sh` but absent and is a hard blocker. |
| Small real/simulated demo dataset | **COMPLETED HERE (smoke demo)** | `demo/` contains tiny deterministic synthetic inputs and expected outputs. For full biological capability, see the end-to-end item below. |
| System requirements documented | **PARTIAL** | `docs/INSTALLATION.md` lists software/resources and HPC assumptions; actual tested OS/tool versions must be captured by authors. |
| Dependencies + version numbers | **PARTIAL** | Dependency list is complete enough for the visible scripts; exact historically tested versions are not inferable. Use `scripts/collect_environment.sh`. |
| Versions software tested on | **AUTHOR ACTION REQUIRED** | Must come from the real analysis/reviewer test environment. |
| Non-standard hardware/resources | **PARTIAL** | SLURM/HPC use and large genomic resources documented. Actual cluster/hardware requirements must be confirmed. |
| Installation instructions | **PARTIAL** | Installation/resource requirements documented, but a fully pinned environment cannot be supplied without real versions. |
| Typical install time | **AUTHOR ACTION REQUIRED** | Must be measured on a stated clean environment. |
| Demo instructions | **COMPLETED HERE** | `python3 demo/run_demo.py`. |
| Expected demo output | **COMPLETED HERE** | `demo/expected/` plus automatic comparisons. |
| Expected demo runtime | **PARTIAL** | Smoke demo is seconds-scale; run-time is printed when executed. Authors should report a measured typical value on their stated machine. |
| Instructions for use on own data | **COMPLETED HERE (historical-template level)** | `docs/RUN_ON_YOUR_DATA.md` maps stages, placeholders, inputs, and outputs. Full execution still requires missing `gtbulksearch_py3.py` and external resources. |
| Reproduction instructions for manuscript results (optional in attached form, strongly useful) | **PARTIAL** | Pipeline overview and use docs added; exact manuscript datasets/commands/figure reproduction are not present. |
| License | **AUTHOR ACTION REQUIRED** | No project license is present. Authors/institution must choose it. NanoSeq's upstream component has its own AGPL terms; see `docs/THIRD_PARTY.md`. |
| Open-source repository link | **COMPLETED HERE** | Public repository: `https://github.com/zonglab/scNanoSeq`. |
| DOI/permanent identifier when available / frozen paper version | **AUTHOR ACTION REQUIRED** | Create a release/tag and archive it (e.g. Zenodo/institutional repository). |
| Detailed code functionality/pseudocode in manuscript | **PARTIAL** | `docs/PIPELINE_OVERVIEW.md` provides repository-side pseudocode; authors must identify the final manuscript section/page/line. |
| External dependencies specified for test dataset | **COMPLETED HERE for smoke demo** | Smoke demo uses Python 3 + NumPy only; broader pipeline dependencies are listed separately. |
| Test dataset demonstrates complete software capability | **AUTHOR ACTION REQUIRED if full pipeline is central** | Current demo validates custom logic but deliberately bypasses BWA/Picard/CellPhy/Ginkgo and does not prove the full biological workflow. Add a small real/simulated end-to-end dataset if reviewers must execute the complete workflow. |
| Independent colleague test | **AUTHOR ACTION REQUIRED** | Nature recommends a colleague unfamiliar with the tool perform installation + demo and provide feedback. |

## Static audit findings

The original supplied scripts pass Bash syntax checks (for `.sh` files) and Python compilation checks, but syntax validity is not the same thing as reproducibility. The highest-priority functional/documentation issues are:

- `gtbulksearch_py3.py` is missing;
- `extract_tags.py` is only a provenance note, not executable implementation;
- multiple scripts contain `xxxxxx`/`xxxxxxxx` and sample/barcode placeholders;
- `1_Bulk_WGS/02_Bulk_variant_calling` depends on undefined shell variables (`chr`, `data_dir`, `sample`);
- the uploaded `1_Bulk_WGS/03_Call_summary.sh` referred to `cento.bed`, while the distributed file is `hg19_cento.bed`; this obvious filename mismatch is corrected in the review-ready bundle;
- the scripts assume `hg19.fa`, SnpSift/dbSNP, and custom scripts are available from the current working directory or `PATH` without documenting that layout;
- the uploaded `4_CNV/README.md` referred to `optimized_process.R`, while the supplied file is `modified_process.R`; this documentation mismatch is corrected here;
- the public GitHub README is only a title/one-line description and currently lacks the checklist documentation, demo, license, and release metadata;
- the uploaded working tree renames several public-repository directories with numeric prefixes, but those renames were not committed in the included `.git` state; `Simulation/` and `VAF_fingerprint/` are absent from the uploaded working tree.

## BWA / CellPhy demo decision

You do **not** need to create a tutorial/demo proving that BWA, SAMtools, Picard, CellPhy, or Ginkgo themselves work. They are third-party dependencies. The hard requirement is to demonstrate **your custom code/software** on example data and to specify external dependencies. If a third-party step is necessary to reproduce a manuscript output, preserve its exact invocation, version, parameters, and expected file interface. For CellPhy, a tiny combined VCF plus the exact CellPhy command and expected tree output is a useful integration example, but not a re-validation of CellPhy itself.
