#!/bin/bash
#SBATCH --job-name=05_call_summary
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=2gb # memory pool for all cores
#SBATCH -t 0-12:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

sample=AYY


patient_ID=`echo ${sample} | awk 'BEGIN { FS = "_" } ; { print $1}'`
seq_dir=xxxxxxxx/${sample}
data_dir=xxxxxxxx/${patient_ID}/${sample}
split_dir=${data_dir}/split
completion_flag_dir=${split_dir}/flag

mkdir -p ${data_dir}
mkdir -p ${split_dir}


cd ${completion_flag_dir}


## Wait for bam_add_META_tag to finish
cat barcode_3digit.list | while read chr
do
    FILE=${completion_flag_dir}/${sample}_${chr}.txt
    while true; do
        if [[ -f "$FILE" ]]; then
           echo "$FILE exists."
           break
          else
           echo $(date +"%T")
           echo "waiting for $FILE"
           sleep 1m
         fi
    done
done

cd ${split_dir}
rm -rf ${completion_flag_dir}

cat ${sample}_[A-Z][A-Z][A-Z]_a4s2_bulk_hetero_detected_double_strand.vcf > ${data_dir}/${sample}_a4s2_bulk_hetero_detected_double_strand.vcf
cat ${sample}_[A-Z][A-Z][A-Z]_a4s2_bulk_hetero_detected_single_strand.vcf > ${data_dir}/${sample}_a4s2_bulk_hetero_detected_single_strand.vcf


cat ${sample}_[A-Z][A-Z][A-Z]_a4s2_Damage_0_bulk_filtered.vcf > ${data_dir}/${sample}_a4s2_Damage_0_bulk_filtered.vcf
cat ${sample}_[A-Z][A-Z][A-Z]_a4s2_Damage_swap_bulk_filtered.vcf > ${data_dir}/${sample}_a4s2_Damage_swap_bulk_filtered.vcf

cat ${sample}_[A-Z][A-Z][A-Z]_a4s2_ss_var_bulk_filtered.vcf > ${data_dir}/${sample}_a4s2_ss_var_bulk_filtered.vcf
cat ${sample}_[A-Z][A-Z][A-Z]_a4s2_mutation_bulk_filtered.vcf  > ${data_dir}/${sample}_a4s2_mutation_bulk_filtered.vcf 


cd ${data_dir}

cat ${sample}_a4s2_mutation_bulk_filtered.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > ${sample}_a4s2_mutation_bulk_filtered_sorted.vcf
java -jar SnpSift.jar annotate All_20151104.vcf ${sample}_a4s2_mutation_bulk_filtered_sorted.vcf > ${sample}_a4s2_mutation_bulk_filtered_sifted.vcf
grep -v rs ${sample}_a4s2_mutation_bulk_filtered_sifted.vcf | grep -v "^#" > ${sample}_a4s2_mutation_dbSNP_filtered.vcf
rm ${sample}_a4s2_mutation_bulk_filtered_sifted.vcf
rm ${sample}_a4s2_mutation_bulk_filtered_sorted.vcf


cat ${sample}_a4s2_ss_var_bulk_filtered.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > ${sample}_a4s2_ss_var_bulk_filtered_sorted.vcf
java -jar SnpSift.jar annotate All_20151104.vcf ${sample}_a4s2_ss_var_bulk_filtered_sorted.vcf > ${sample}_a4s2_ss_var_bulk_filtered_sifted.vcf
grep -v rs ${sample}_a4s2_ss_var_bulk_filtered_sifted.vcf | grep -v "^#" > ${sample}_a4s2_ss_var_dbSNP_filtered.vcf
rm ${sample}_a4s2_ss_var_bulk_filtered_sifted.vcf
rm ${sample}_a4s2_ss_var_bulk_filtered_sorted.vcf




python3 De_novo_estimation_by_variant_position.py ${data_dir} ${sample}_a4s2_bulk_hetero_detected_double_strand.vcf ${sample}_a4s2_bulk_hetero_detected_single_strand.vcf ${sample}_a4s2_ss_var_dbSNP_filtered.vcf ${sample}_a4s2_mutation_dbSNP_filtered.vcf `wc -l xxxxxx/bulk/${patient_ID}/samcall/snp.total.hetero.filtered.vcf | awk '{print $1}'` `wc -l ${sample}_a4s2_Damage_0_bulk_filtered.vcf | awk '{print $1}'` `wc -l ${sample}_a4s2_Damage_swap_bulk_filtered.vcf | awk '{print $1}'` ${sample}_a4s2_de_novo_estimation.txt


rm -rf ${split_dir}


