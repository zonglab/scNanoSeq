#!/bin/bash
#SBATCH --job-name=04_BYY_raw_calling
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=2gb # memory pool for all cores
#SBATCH -t 0-12:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR



sample=AYY
barcode=BYY

patient_ID=`echo ${sample} | awk 'BEGIN { FS = "_" } ; { print $1}'`
seq_dir=xxxxxxxx/${sample}
data_dir=xxxxxxxx/${patient_ID}/${sample}
split_dir=${data_dir}/split

completion_flag_dir=${split_dir}/flag

mkdir -p ${data_dir}
mkdir -p ${split_dir}
cd ${split_dir}


samtools index ${sample}_${barcode}.bam

samtools mpileup -uf hg19.fa ${sample}_${barcode}.bam | bcftools call -mv | gzip  > ${barcode}_var.raw.vcf.gz


zcat ${barcode}_var.raw.vcf.gz > ${barcode}_var.raw.vcf
bedtools intersect -a ${barcode}_var.raw.vcf  -b hg19_tandem.bed -v > ${barcode}_var.notandem.body

grep ^# ${barcode}_var.raw.vcf > ${barcode}_vcf.header
cat ${barcode}_vcf.header ${barcode}_var.notandem.body> ${barcode}_var.notandem.vcf
rm ${barcode}_var.notandem.body

bedtools intersect -a ${barcode}_var.notandem.vcf  -b cento.bed -v > ${barcode}_var.cleaned.body

rm ${barcode}_var.notandem.vcf
rm ${barcode}_var.raw.vcf

cat ${barcode}_var.cleaned.body | awk '{print $1"\t"$2-10"\t"$2+10}' > ${barcode}_var.localseq.bed
seqtk subseq hg19.fa   ${barcode}_var.localseq.bed -t  > ${barcode}_var.localseq.txt

python3 localseqfilter_py3.py ${barcode}_var.cleaned.body ${barcode}_var.localseq.txt | gzip > ${barcode}_var.filtered.vcf.gz

rm ${barcode}_var.cleaned.body
rm ${barcode}_var.localseq.txt
rm ${barcode}_var.localseq.bed


zcat ${barcode}_var.raw.vcf.gz | grep ^# > ${sample}_${barcode}_vcf.header

python3 NanoSeq_Mutation_call_aMsN_for_hg19.py ${split_dir} ${barcode}_var.filtered.vcf.gz ${sample}_${barcode}.bam 4 2 ${sample}_${barcode}_a4s2_call.txt




grep Mutation ${sample}_${barcode}_a4s2_call.txt | cut -f2- - | cat ${sample}_${barcode}_vcf.header - > ${sample}_${barcode}_a4s2_mutation.vcf
grep Single_strand_variant_allele_filtered ${sample}_${barcode}_a4s2_call.txt | cut -f2- - | cat ${sample}_${barcode}_vcf.header - > ${sample}_${barcode}_a4s2_ss_var.vcf
grep Damage_0 ${sample}_${barcode}_a4s2_call.txt | cut  -f2-  -   | cat ${sample}_${barcode}_vcf.header - > ${sample}_${barcode}_a4s2_Damage_0.vcf 
grep Damage_swap_strand_filtered ${sample}_${barcode}_a4s2_call.txt | cut  -f2-  -   | cat ${sample}_${barcode}_vcf.header - > ${sample}_${barcode}_a4s2_Damage_swap.vcf 


bedtools intersect -a ${sample}_${barcode}_a4s2_mutation.vcf -b xxxxxx/bulk/${patient_ID}/samcall/snp.total.hetero.filtered_with_header.vcf  > ${sample}_${barcode}_a4s2_bulk_hetero_detected_double_strand.vcf

bedtools intersect -a ${sample}_${barcode}_a4s2_ss_var.vcf  -b xxxxxx/bulk/${patient_ID}/samcall/snp.total.hetero.filtered_with_header.vcf  > ${sample}_${barcode}_a4s2_bulk_hetero_detected_single_strand.vcf


python3 gtbulksearch_py3.py ${sample}_${barcode}_a4s2_Damage_0.vcf xxxxxx/bulk/${patient_ID}/${patient_ID}.RG.markdup.bam 0 > ${sample}_${barcode}_a4s2_Damage_0_bulk_filtered.vcf



python3 gtbulksearch_py3.py ${sample}_${barcode}_a4s2_Damage_swap.vcf  xxxxxx/bulk/${patient_ID}/${patient_ID}.RG.markdup.bam 0 > ${sample}_${barcode}_a4s2_Damage_swap_bulk_filtered.vcf

python3 gtbulksearch_py3.py ${sample}_${barcode}_a4s2_mutation.vcf  xxxxxx/bulk/${patient_ID}/${patient_ID}.RG.markdup.bam 0 > ${sample}_${barcode}_a4s2_mutation_bulk_filtered.vcf

python3 gtbulksearch_py3.py ${sample}_${barcode}_a4s2_ss_var.vcf  xxxxxx/bulk/${patient_ID}/${patient_ID}.RG.markdup.bam 0 > ${sample}_${barcode}_a4s2_ss_var_bulk_filtered.vcf



mkdir -p ${completion_flag_dir}
echo "Complete" > ${completion_flag_dir}/${sample}_${barcode}.txt

