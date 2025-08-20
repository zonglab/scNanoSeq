#!/usr/local/bin/python3.6

import os
import sys
import gzip
import pysam

##zcat sample.vcf.gz | bgzip > sample_bg.vcf.gz
##bcftools index sample_bg.vcf.gz



wd,Fvcf, Fbam, Out_file_name = sys.argv[1],sys.argv[2],sys.argv[3], sys.argv[6]
a_threshold, s_threshold = int(sys.argv[4]), int(sys.argv[5])



os.chdir(wd)
samfile = pysam.AlignmentFile(Fbam, "rb")

def n_lower_chars(string):
    return sum(1 for c in string if c.islower())

def bulk_count(bulk_stat_in_list):
    Sum, Flag, Homo = False, 0, ''
    stat=bulk_stat_in_list.split(":")
    Ref = int(stat[1].split(",")[0])+int(stat[2].split(",")[0])
    Alt = int(stat[1].split(",")[1])+int(stat[2].split(",")[1])
    Sum = Ref + Alt
    if (Alt == 0) and (Sum >= 20):
        Flag = True
        Homo = "Ref"
    if (Ref == 0) and (Sum >= 20):
        Flag = True
        Homo = "Alt"
    return(Sum, Flag, Homo)


def scSNV_assignment(TUPLE_Fwt_Fmut_Rwt_Rmut, allele_read, strand_read):
    Fwt,Fmut,Rwt,Rmut = TUPLE_Fwt_Fmut_Rwt_Rmut[0], TUPLE_Fwt_Fmut_Rwt_Rmut[1], TUPLE_Fwt_Fmut_Rwt_Rmut[2], TUPLE_Fwt_Fmut_Rwt_Rmut[3]
    Fsum, Rsum = Fwt+Fmut, Rwt+Rmut
    assignment = "*"
    if (Fsum+Rsum>=allele_read) and (min(Fsum,Rsum)>=strand_read):
        if Fwt==0 and Rwt==0:
            assignment="Mutation"
        elif Fwt*Rwt*Fmut*Rmut ==0:
            if min(Rwt+Fmut, Fwt+Rmut)==0:
                assignment="Damage_all_transformed"
                if max(Fwt, Fmut, Rwt, Rmut)>=allele_read:
                    assignment="Damage_0"
            elif (Rmut+Fmut ==0):
                pass
            elif (Rmut==0) or (Fmut==0):
                assignment="Damage_mis"
            else:
                assignment="Damage_swap"
                if min(Fmut,Rmut)>=strand_read:
                    if max(Fwt, Rwt)>=strand_read:
                        assignment="Damage_swap_strand_filtered"
        else:
            assignment="Swap&Mis" 
    elif (Fsum+Rsum>=allele_read) and (min(Fsum,Rsum)==0) and (Fwt==0) and (Rwt==0):
        assignment="Single_strand_variant_allele_filtered"

    return assignment
            
   
    
with gzip.open(Fvcf,mode="rt") as f, open(Out_file_name,"w") as out_file:
    for line in f:
        if line.startswith("#"):
            if not line.startswith("##"):
                line = line.replace("#","#TAG\t")
            out_file.write(line)
        else:
            x1=line.strip().split("\t")
            chrnum, pos, snpref, snpvar = x1[0], int(x1[1]), x1[3], x1[4]
            if chrnum not in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']:
                continue           
            if (len(snpref)!=1) or (len(snpvar)!=1):
                continue

            flag_dict={}
            for read in samfile.fetch(chrnum, pos-1, pos):
                try:
                    if int(read.cigarstring.split("S")[0]) >= 10:
                        continue
                except:
                    pass
                flag_dict[read]=read.flag            

            read_83 = [key for (key,values)in flag_dict.items() if values==83]
            read_147 = [key for (key,values)in flag_dict.items() if values==147]
            read_99 = [key for (key,values)in flag_dict.items() if values==99]
            read_163 = [key for (key,values)in flag_dict.items() if values==163]

            if len(read_83)+len(read_147)> len(read_99)+len(read_163):
                read_use=read_83+read_147
            else:
                read_use=read_99+read_163

            count_per_MapC={}
            position_per_MapC={}
            barcode_per_MapC={}

            for read in read_use:
                # read.query_name
                positions = read.get_reference_positions()
                start, end, read_name = positions[0],positions[-1], read.query_name
                if read.flag==83 or read.flag==99:
                    direction = "F"
                elif read.flag==147 or read.flag==163:
                    direction = "R"

                try:
                    base = read.query_alignment_sequence[positions.index(pos-1)]
                    quality = read.query_alignment_qualities[positions.index(pos-1)]
                except ValueError:
                    base, quality = "", 0


                # by mapping cooridinate, allow 1 base shift
                Map_C = (start,read.next_reference_start,dict(read.tags)['MR'])
                if Map_C in count_per_MapC:
                    pass
                #elif (Map_C-1) in count_per_MapC:
                #    Map_C=Map_C-1
                else:
                    count_per_MapC[Map_C]=(0,0,0,0)
                    
                (Fwt,Fmut,Rwt,Rmut)=count_per_MapC[Map_C]  


#dict(read.get_aligned_pairs())

                if quality<30:
                    continue
                if ("I" in read.cigarstring) or ("D" in read.cigarstring):
                    continue
                if n_lower_chars(read.get_reference_sequence())>3:
                    continue


                position_in_readX_from_proximal,position_in_readX_from_distal= 0,0
                position_to_partial_read_end = 0 

                if (read.flag in [99,163]):
                    position_in_readX_from_proximal = positions.index(pos-1)+1
                    position_in_readX_from_distal = abs(read.isize) +1 - position_in_readX_from_proximal
                    position_to_partial_read_end = read.alen+1 - position_in_readX_from_proximal
                if (read.flag in [83,147]):
                    position_in_readX_from_proximal = read.alen - positions.index(pos-1)
                    position_in_readX_from_distal = abs(read.isize) +1 - position_in_readX_from_proximal
                    position_to_partial_read_end = read.alen+1 - position_in_readX_from_proximal

                if position_to_partial_read_end <= 0:  
                    continue
                if position_in_readX_from_proximal <= 8:
                    continue
                if position_in_readX_from_distal <= 8:
                    continue
                
                if   direction =="F" and base ==snpref:
                      Fwt+=1
                if   direction =="F" and base ==snpvar:
                      Fmut+=1
                if   direction =="R" and base ==snpref:
                      Rwt+=1                
                if   direction =="R" and base ==snpvar:
                      Rmut+=1
                count_per_MapC[Map_C]=(Fwt,Fmut,Rwt,Rmut)

                position_per_MapC[Map_C]=(position_in_readX_from_proximal,position_in_readX_from_distal)
                barcode_per_MapC[Map_C]=dict(read.tags)['MR']
                
            sum_count_per_MapC={}
            for (key,value) in count_per_MapC.items():
                sum_count_per_MapC[key] = value[0]+value[1]+value[2]+value[3]
                
            try:
                Key_max = max(sum_count_per_MapC, key=sum_count_per_MapC.get)
                assignment = scSNV_assignment(count_per_MapC[Key_max], a_threshold, s_threshold)
                if assignment != "*":
                    line = assignment + "\t" + line.strip() + "\t" + str(count_per_MapC[Key_max]) + "\t" + str(position_per_MapC[Key_max])+ "\t"+ str(barcode_per_MapC[Key_max]) +"\n"
                    out_file.write(line)

            except:
                pass


samfile.close()
      



