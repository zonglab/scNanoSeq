#!/bin/bash

ls -1v snp.chr*.txt |xargs cat |grep ^chr | grep "0/1" > snp.total.hetero.vcf

bcftools view output.chr1.bcf | grep ^# > vcf.header

cat vcf.header snp.total.hetero.vcf > snp.total.hetero_header.vcf

bedtools intersect -a snp.total.hetero_header.vcf  -b hg19_tandem.bed -v > 1.body
cat vcf.header 1.body > 1.vcf


bedtools intersect -a 1.vcf  -b cento.bed -v > snp.total.hetero.filtered1.vcf





rm 1.vcf
rm 1.body


cat snp.total.hetero.filtered1.vcf | awk '{print $1"\t"$2-10"\t"$2+10}' > var.localseq.bed


seqtk subseq hg19.fa   var.localseq.bed -t  > var.localseq.txt

python3 localseqfilter_py3.py snp.total.hetero.filtered1.vcf var.localseq.txt > snp.total.hetero.filtered.vcf

cat vcf.header snp.total.hetero.filtered.vcf > snp.total.hetero.filtered_with_header.vcf


rm var.localseq.bed
rm var.localseq.txt
rm snp.total.hetero_header.vcf
rm snp.total.hetero.filtered1.vcf




ls -1v snp.chr*.txt |xargs cat |grep ^chr > snp.total.all.body
cat vcf.header snp.total.all.body > snp.total.all.vcf
rm snp.total.all.body



