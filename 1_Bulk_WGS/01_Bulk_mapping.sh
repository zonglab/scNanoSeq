#!/bin/bash
#SBATCH --job-name=AXX_bulk_mapping
#SBATCH -N 1 # number of nodes
#SBATCH -n 32 # number of cores
#SBATCH --mem=15gb # memory pool for all cores
#SBATCH -t 0-72:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR


sample=BTT1

seq_dir=xxxxxx/BTT1_normal_bulk
data_dir=xxxxxx/bulk/${sample}

mkdir -p ${data_dir}
cd ${data_dir}


cutadapt -j 32 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -m 50 -o ${sample}_trimmed_R1.fastq.gz -p ${sample}_trimmed_R2.fastq.gz ${seq_dir}/${sample}_R1.fq.gz ${seq_dir}/${sample}_R2.fq.gz > Nextera_trimming_report.txt


bwa mem -M -t 32 hg19.fa ${sample}_trimmed_R1.fastq.gz ${sample}_trimmed_R2.fastq.gz | samtools view -bS - > ${sample}.bam

rm ${sample}_trimmed_R1.fastq.gz
rm ${sample}_trimmed_R2.fastq.gz

samtools fixmate -m -@ 32 -O bam ${sample}.bam ${sample}.fixmate.bam

if [ -f ${sample}.fixmate.bam ]; then
   rm -rf ${sample}.bam
fi

samtools sort -@ 32 ${sample}.fixmate.bam -o ${sample}.fixmate.sorted.bam 


if [ -f ${sample}.fixmate.sorted.bam ]; then
   rm -rf ${sample}.fixmate.bam
fi

samtools markdup -@ 32 -r ${sample}.fixmate.sorted.bam ${sample}.fixmate.markdup.bam

samtools index -@ 32 ${sample}.fixmate.markdup.bam

if [ -f ${sample}.fixmate.markdup.bam ]; then
   rm -rf ${sample}.fixmate.sorted.bam
fi






picard  AddOrReplaceReadGroups \
       I=${sample}.fixmate.markdup.bam \
       O=${sample}.RG.markdup.bam \
       RGID=1 \
       RGLB=M \
       RGPL=X10 \
       RGPU=1 \
       RGSM=Bulk

samtools index -@ 8 ${sample}.RG.markdup.bam


if [ -f ${sample}.RG.markdup.bam ]; then
   rm -rf ${sample}.fixmate.markdup.bam ${sample}.fixmate.markdup.bam.bai
fi



