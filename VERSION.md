# Version provenance

This repository did not contain a semantic version tag in the supplied Git history.

| Field | Value |
| --- | --- |
| Repository | `zonglab/scNanoSeq` |
| Base branch | `main` |
| Base commit | `2f1c8cb4f131a57c7c8e86fc1f2d699a968545d7` |
| Base commit date | 2025-08-20 |
| Base `git describe` | `2f1c8cb-dirty` in the uploaded snapshot |
| Documentation/demo refresh | 2026-09-16 |
| Reference build assumed by scripts | hg19 / GRCh37 with `chr` contig prefix |
| Core single-cell thresholds | `a4s2` = total allele reads ≥4 and per-strand reads ≥2 |

The uploaded working tree had been reorganized from the historical directory names (`ref`, `Bulk_WGS`, `Somatic_mutation_calling`, etc.) into numbered directories (`0_ref`, `1_Bulk_WGS`, `2_Somatic_mutation_calling`, etc.) without a corresponding commit in the embedded Git history.

For a publication/release, use the final release commit hash as the authoritative code version and create an explicit Git tag (for example, `v1.0.0`) only after the unresolved release-blocking items in `RELEASE_CHECKLIST.md` are addressed.
