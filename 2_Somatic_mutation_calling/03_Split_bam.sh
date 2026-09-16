#!/bin/bash
#SBATCH --job-name=03_Split_bam
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=2gb # memory pool for all cores
#SBATCH -t 0-72:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR


sample=AXX


patient_ID=`echo ${sample} | awk 'BEGIN { FS = "_" } ; { print $1}'`
seq_dir=xxxxxxxx/${sample}
data_dir=xxxxxxxx/${patient_ID}/${sample}
split_dir=${data_dir}/split

mkdir -p ${data_dir}
mkdir -p ${split_dir}
cd ${data_dir}


mkdir -p ${split_dir}
python3 bam_split_by_allele_barcode.py ${data_dir} ${sample}_DuplexSeq_sorted.bam ${split_dir}/${sample}







#### sbatch 04
mkdir -p xxxxxxxx/${patient_ID}/run_log
cd xxxxxxxx/${patient_ID}/run_log

LINE1=${sample}
cat barcode_3digit.list | while read LINE
do
 echo $LINE
 sed "s/AYY/${LINE1}/g" xxxxxxxx/${patient_ID}/script/04_Variant_call.sh  | sed "s/BYY/${LINE}/g" - > ./04_Variant_call_${LINE1}_${LINE}.sh
 sbatch  ./04_Variant_call_${LINE1}_${LINE}.sh
 #sleep 5m
done

#### sbatch 05
mkdir -p xxxxxxxx/${patient_ID}/run_log
cd xxxxxxxx/${patient_ID}/run_log

sed "s/AYY/${sample}/g" xxxxxxxx/${patient_ID}/script/05_Call_summary.sh > ./05_Call_summary_${sample}.sh
sbatch  ./05_Call_summary_${sample}.sh




