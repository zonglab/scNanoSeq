#!/bin/bash
#SBATCH --job-name=01_PreProcessing
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=3gb # memory pool for all cores
#SBATCH -t 0-12:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

sample=BTT1_sc001



patient_ID=`echo ${sample} | awk 'BEGIN { FS = "_" } ; { print $1}'`
seq_dir=xxxxxx/${sample}
data_dir=xxxxxx/${patient_ID}/${sample}
split_dir=${data_dir}/split

mkdir -p ${data_dir}
mkdir -p ${split_dir}
cd ${data_dir}


zcat ${seq_dir}/*1.fq.gz> R1.fastq
zcat ${seq_dir}/*2.fq.gz> R2.fastq
python3 extract_tags.py -a R1.fastq -b R2.fastq -c extrR1.fastq -d extrR2.fastq -m 3 -s 4 -l 150

rm R1.fastq
rm R2.fastq

gzip extrR1.fastq
gzip extrR2.fastq

cutadapt -j 2 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 50 -o extrR1_trimmed.fastq.gz -p extrR2_trimmed.fastq.gz extrR1.fastq.gz extrR2.fastq.gz > trimming_report.txt

rm extrR1.fastq.gz
rm extrR2.fastq.gz


#### sbatch 02
mkdir -p xxxxxx/${patient_ID}/run_log
cd xxxxxx/${patient_ID}/run_log

sed "s/AYY/${sample}/g" xxxxxx/${patient_ID}/script/02_Mapping.sh > ./02_Mapping_${sample}.sh
sbatch  ./02_Mapping_${sample}.sh




