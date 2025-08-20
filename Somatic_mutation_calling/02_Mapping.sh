#!/bin/bash
#SBATCH --job-name=02_mapping
#SBATCH -N 1 # number of nodes
#SBATCH -n 6 # number of cores
#SBATCH --mem=15gb # memory pool for all cores
#SBATCH -t 0-120:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

sample=AYY


patient_ID=`echo ${sample} | awk 'BEGIN { FS = "_" } ; { print $1}'`
seq_dir=xxxxxx/${sample}
data_dir=xxxxxx/${patient_ID}/${sample}
split_dir=${data_dir}/split

mkdir -p ${data_dir}
mkdir -p ${split_dir}
cd ${data_dir}


bwa mem -t 6 -C hg19.fa extrR1_trimmed.fastq.gz extrR2_trimmed.fastq.gz | samtools view -bS - > output.bwa.bam

rm extrR1_trimmed.fastq.gz
rm extrR2_trimmed.fastq.gz


samtools view -@ 6 -b -q 50 -f 3 output.bwa.bam | samtools sort -@ 6 - > ${sample}_DuplexSeq_sorted0.bam

samtools index -@ 6 ${sample}_DuplexSeq_sorted0.bam

rm output.bwa.bam

picard MarkDuplicates \
    I=${sample}_DuplexSeq_sorted0.bam \
    O=${sample}_DuplexSeq_sorted.bam \
    M=${sample}_optical_duplicate.txt \
    TAGGING_POLICY=All \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 \
    REMOVE_SEQUENCING_DUPLICATES=true

samtools index -@ 6 ${sample}_DuplexSeq_sorted.bam

rm ${sample}_DuplexSeq_sorted0.bam
rm ${sample}_DuplexSeq_sorted0.bam.bai



#### sbatch 03
mkdir -p xxxxxx/${patient_ID}/run_log
cd xxxxxx/${patient_ID}/run_log

sed "s/AXX/${sample}/g" xxxxxx/${patient_ID}/script/03_Split_bam.sh > ./03_Split_bam_${sample}.sh
sbatch  ./03_Split_bam_${sample}.sh




