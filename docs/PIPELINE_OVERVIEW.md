# Pipeline overview / pseudocode

This document provides a compact algorithmic description suitable for cross-referencing from a manuscript Methods section. It is descriptive documentation, not a replacement for the manuscript's own detailed method description.

## Bulk WGS branch

1. Trim adapters from paired-end bulk reads.
2. Align reads to hg19 with BWA-MEM.
3. Fix mate information, coordinate-sort, remove duplicates, add read groups, and index the BAM.
4. Generate chromosome-wise pileups and call variants with SAMtools/BCFtools.
5. Retain high-quality heterozygous SNPs and remove tandem-repeat, centromeric, and low-complexity local-sequence loci.
6. Use the resulting bulk calls as a matched-normal/germline filter for single-cell candidates.

## Single-cell Duplex-seq branch

1. Concatenate sequencing lanes for each cell.
2. Extract three-base molecular/barcode tags from paired FASTQ reads using the pinned NanoSeq tag-extraction implementation; trim adapters.
3. Map with BWA-MEM and retain properly paired, high-mapping-quality reads.
4. Mark/remove sequencing/optical duplicates as configured in the historical script.
5. Split the BAM by the three-base allele/molecular barcode, preserving an `MR` tag used by downstream grouping.
6. For each barcode:
   - call candidate variants;
   - exclude tandem-repeat and centromeric regions;
   - exclude low-complexity local sequence;
   - examine high-quality reads supporting each candidate;
   - group reads by mapping coordinates and molecular barcode;
   - classify the locus using forward/reverse reference/alternate counts (`scSNV_assignment`);
   - compare candidates with matched bulk data and remove likely germline/bulk-supported loci.
7. Merge barcode-level outputs, annotate/filter dbSNP loci, and estimate de-novo burden versus distance from the read/fragment end.

## Phylogeny branch

1. Collect per-cell de-novo mutation loci.
2. Re-genotype the union of loci across bulk and single-cell BAMs to form a multi-sample VCF.
3. Retain SNPs.
4. Run CellPhy with the recorded thread count/parameters to infer the lineage tree.

## CNV branch

1. Convert per-cell alignment data to the BED/count representation required by Ginkgo.
2. Run the documented Ginkgo standalone workflow using the modified process script.
3. Use the resulting segmented copy-number profiles for downstream CNV interpretation.

## Manuscript placement

The Nature checklist asks the manuscript to provide a complete, detailed description of code functionality/pseudocode. **AUTHOR ACTION REQUIRED:** point the checklist to the exact manuscript section/page/line where these operations are described. This repository cannot truthfully supply that manuscript location without the final manuscript.
