#!/bin/bash
#SBATCH --job-name=phy_AXXX
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=10gb # memory pool for all cores
#SBATCH -t 1-12:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

patient_ID=AXXX

data_dir=xxxxxx/${patient_ID}
cellphy_dir=xxxxxx/cellphy_0/${patient_ID}
cellphy_output_dir=xxxxxx/cellphy_0/cellphy_matrix

rm -rf ${cellphy_dir}
mkdir -p ${cellphy_dir}
mkdir -p ${cellphy_output_dir}
cd ${cellphy_dir}

ls ${data_dir} | grep ${patient_ID}_sc > totalcellslist.txt


############
rm -rf zz_${patient_ID}_bam_list.txt
echo xxxxxx/bulk/${patient_ID}/${patient_ID}.RG.markdup.bam >> zz_${patient_ID}_bam_list.txt
cat totalcellslist.txt| while read line 
do
  cat xxxxxx/bulk/AAA_vcf.header ${data_dir}/${line}/${line}_a4s2_mutation_dbSNP_filtered.vcf > ${line}_a4s2_mutation_dbSNP_filtered_header.vcf
  echo xxxxxx/${patient_ID}/${line}/${line}_DuplexSeq_sorted.bam >> zz_${patient_ID}_bam_list.txt

done

python3 de_novo_loci.py ${cellphy_dir} totalcellslist.txt _a4s2_mutation_dbSNP_filtered_header.vcf zz_${patient_ID}_de_novo_loci.txt

cat zz_${patient_ID}_de_novo_loci.txt | sort -k1,1 -k2,2n > zz_${patient_ID}_de_novo_loci_sorted.txt
rm zz_${patient_ID}_de_novo_loci.txt

rm ${patient_ID}_sc*_a4s2_mutation_dbSNP_filtered_header.vcf
######
samtools mpileup -uf hg19.fa -l zz_${patient_ID}_de_novo_loci_sorted.txt -b zz_${patient_ID}_bam_list.txt | bcftools call -m > zz_${patient_ID}_combined.vcf

bcftools view -i 'TYPE="SNP"' zz_${patient_ID}_combined.vcf  > zz_${patient_ID}_combined_SNP.vcf

cellphy.sh -t 8 zz_${patient_ID}_combined_SNP.vcf 



