#!/bin/bash
#SBATCH --job-name=BTXX_AF
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=10gb # memory pool for all cores
#SBATCH -t 1-12:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

source ~/.bashrc

patient_ID=BTT4

data_dir=xxxxxx/${patient_ID}
UGcall_dir=xxxxxx/UG_loci_call/${patient_ID}
UGcall_output_dir=xxxxxx/UG_loci_call/summary

rm -rf ${UGcall_dir}
mkdir -p ${UGcall_dir}
mkdir -p ${UGcall_output_dir}
cd ${UGcall_dir}

ls ${data_dir} | grep ${patient_ID}_sc > totalcellslist.txt


############
rm -rf zz_${patient_ID}_bam_list.txt
echo xxxxxx/bulk/${patient_ID}/${patient_ID}.RG.markdup.bam >> zz_${patient_ID}_bam_list.txt
cat totalcellslist.txt| while read line 
do
  cat xxxxxx/bulk/AAA_vcf.header ${data_dir}/${line}/${line}_a4s2_mutation_dbSNP_filtered.vcf > ${line}_a4s2_mutation_dbSNP_filtered_header.vcf
done


python3 de_novo_loci.py ${UGcall_dir} totalcellslist.txt _a4s2_mutation_dbSNP_filtered_header.vcf zz_${patient_ID}_de_novo_loci.txt

cat zz_${patient_ID}_de_novo_loci.txt | sort -k1,1 -k2,2n > zz_${patient_ID}_de_novo_loci_sorted.bed
rm zz_${patient_ID}_de_novo_loci.txt


rm -rf zzc_${patient_ID}_de_novo_loci_summary.txt
for cn in {1..32}
do
 lofreq call -l zz_${patient_ID}_de_novo_loci_sorted.bed -a 1.000 --min-mq 60 --min-alt-bq 30 -f hg19.fa  --plp-summary-only  xxxxxx/Ultradeep/interval.${cn}.combined.bam >> zzc_${patient_ID}_de_novo_loci_summary.txt
done

cat zzc_${patient_ID}_de_novo_loci_summary.txt | grep ^chr > zzd_${patient_ID}_de_novo_loci_summary.txt

cp zzd_${patient_ID}_de_novo_loci_summary.txt ${UGcall_output_dir}





