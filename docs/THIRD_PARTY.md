# Third-party software and code provenance

This repository orchestrates several established third-party tools. The journal demo should focus on the new/custom code and workflow behavior, while third-party tools should be documented with exact versions, citations, parameters, and inputs/outputs.

## NanoSeq tag extraction

The historical `2_python_script/extract_tags.py` file is not an implementation; it contains a note directing users to NanoSeq:

- repository: `cancerit/NanoSeq`
- pinned commit recorded by the authors: `fd22dec943a3a9cc70643f079c40a8baf13a8e62`
- upstream path: `python/extract_tags.py`
- upstream license at that commit: GNU Affero General Public License v3 or later.

For a publication release, choose one of these reproducible approaches:

1. declare NanoSeq as an external dependency and provide exact install/fetch commands pinned to the commit above; or
2. vendor the exact upstream file while preserving its copyright/license notices and satisfying the AGPL redistribution requirements.

Do **not** silently substitute a newly written implementation without validating scientific equivalence.

A helper is included to fetch that exact pinned upstream file without redistributing it in this bundle:

```bash
bash scripts/fetch_nanoseq_extract_tags.sh
```

The helper stores it under `.external/` and prints the source/license information.

## CellPhy

CellPhy is an external phylogenetic inference program. The repository should record the exact version/build, citation, and the command-line parameters used. The software itself does not need to be reimplemented or independently validated in this repository.

## BWA/SAMtools/BCFtools/BEDTools/Seqtk/Picard/Cutadapt/SnpSift/Ginkgo

These are external tools. A journal-facing demo does not need to teach reviewers how BWA or SAMtools work. What is required for reproducibility is:

- exact version/build;
- exact invocation/parameters as used here;
- required reference/database versions;
- clear inputs and expected output files;
- installation instructions or links sufficient to recreate the environment.
