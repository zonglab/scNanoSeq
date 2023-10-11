import os
import pandas as pd
import numpy as np

sample_name='BTT4'

wd='UG_bulk_loci_call/'+sample_name
cell_list_file='totalcellslist.txt'
cell_ref_file='01_breast_cancer_mutation_all_summary_842.txt'
bulk_af_file='zzc_'+sample_name+'_de_novo_loci_summary.txt'
output_file='zze_'+sample_name+"_AF_q30.txt"
ks_test_file='zzf_'+sample_name+"_ks_test_result.txt"


df_cell_ref=pd.read_csv(cell_ref_file, sep="\t", index_col=0)
cell_list_normal=df_cell_ref[df_cell_ref['Tumor_cell_tree']!=1].index.to_list()
cell_list_tumor=df_cell_ref[df_cell_ref['Tumor_cell_tree']==1].index.to_list()


cell_list=[i for i in df_cell_ref.index.to_list() if sample_name in i]


os.chdir(wd)

count_dict={}   #A,C,G,T
with open(bulk_af_file,'r') as f:
  current_variant=0
  for line in f:
    lineX=line.strip().split("\t")
    if line.startswith("chr"):
      current_variant+=1
      locus=lineX[0]+':'+lineX[1]
      count_A=0
      count_C=0
      count_G=0
      count_T=0

    if line.startswith("  A\tBQ ="):
      qA=lineX[2:]
      qA=qA[0][1:].split(' ')
      qA_int=[int(i.replace(" ",'')) for i in qA]
      count_A=sum(1 if i>=30 else 0 for i in qA_int)

    if line.startswith("  C\tBQ ="):
      qC=lineX[2:]
      qC=qC[0][1:].split(' ')
      qC_int=[int(i.replace(" ",'')) for i in qC]
      count_C=sum(1 if i>=30 else 0 for i in qC_int)

    if line.startswith("  G\tBQ ="):
      qG=lineX[2:]
      qG=qG[0][1:].split(' ')
      qG_int=[int(i.replace(" ",'')) for i in qG]
      count_G=sum(1 if i>=30 else 0 for i in qG_int)

    if line.startswith("  T\tBQ ="):
      qT=lineX[2:]
      qT=qT[0][1:].split(' ')
      qT_int=[int(i.replace(" ",'')) for i in qT]
      count_T=sum(1 if i>=30 else 0 for i in qT_int)

    if line=="\n":
      count_dict[locus]=(int(count_A),int(count_C),int(count_G),int(count_T))




with open(output_file,'w') as out:
  line_out="\t".join(['cell_ID','Tumor_cell_tree','clade','locus','base_ref','base_var','depth','AF_var','AF_other','AF_other_major'])+"\n"
  out.write(line_out)
  for cellX in cell_list:
    vcfX=cellX+'_a4s2_mutation_dbSNP_filtered_header.vcf'
    
    try:
      with open(vcfX,'r') as f:
        for line in f:
          lineX=line.strip().split("\t")
          if line.startswith("#"):
            continue
          locus=lineX[0]+':'+lineX[1]
          base_ref=lineX[3]
          base_var=lineX[4]

          try:
            countA, countC, countG, countT = count_dict[locus]
          except KeyError:
            continue

          count_ref = globals()['count'+base_ref]
          count_var = globals()['count'+base_var]
          list_ATCG=['A','C','G','T']
          list_ATCG.remove(base_ref)
          list_ATCG.remove(base_var)
          count_other_major=max(globals()['count'+list_ATCG[0]],globals()['count'+list_ATCG[1]])


          count_total=countA+countC+countG+countT

          if count_total<100:
            continue

          count_other = count_total-count_ref-count_var

          AF_ref=count_ref/count_total
          AF_var=count_var/count_total
          AF_other=count_other/count_total
          AF_other_major=count_other_major/count_total

          tumor_cell_tree=str(df_cell_ref['Tumor_cell_tree'][cellX])
          clade=str(df_cell_ref['Clade'][cellX])

          line_out = "\t".join([cellX,tumor_cell_tree,clade,locus,base_ref,base_var,str(count_total),str(AF_var),str(AF_other),str(AF_other_major)])+"\n"
          out.write(line_out)


    except FileNotFoundError:
      continue




######ks test

df=pd.read_csv(output_file,sep="\t",index_col=None)
cell_list=list(set(df['cell_ID'].to_list()))
cell_list.sort()
cell_num=len(cell_list)


dg=pd.DataFrame(columns=['Cell_A','Tumor_cell_tree_A','Clade_A','Cell_B','Tumor_cell_tree_B','Clade_B','Distance','p_value'])

from scipy.stats import ks_2samp

df['AF_clip'] = df['AF_var'].clip(lower=0.001)
df['AF_log2']= np.log2(df['AF_clip'])


for i in range(0,cell_num):
  for j in range(i+1,cell_num):
    celli=cell_list[i]
    cellj=cell_list[j]
    AFi=df[df['cell_ID']==celli]['AF_log2'].to_numpy()
    AFj=df[df['cell_ID']==cellj]['AF_log2'].to_numpy()
    ks_stat=ks_2samp(AFi, AFj, alternative='two-sided')
    Distance=ks_stat[0]
    p_value=ks_stat[1]

    dg0 = {'Cell_A':celli,
    'Tumor_cell_tree_A':df_cell_ref['Tumor_cell_tree'][celli],
    'Clade_A':df_cell_ref['Clade'][celli],
    'Cell_B':cellj,
    'Tumor_cell_tree_B':df_cell_ref['Tumor_cell_tree'][cellj],
    'Clade_B':df_cell_ref['Clade'][cellj],
    'Distance':Distance,
    'p_value':p_value}
    dg = dg.append(dg0, ignore_index = True)


def BH_FDR(p_value_in_Series):
    pvals=p_value_in_Series.array.fillna(1)
    import statsmodels.api
    fdr=statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05, method='i', is_sorted=False)[1] 
    FDR_in_Series = pd.Series(fdr.tolist(), index=p_value_in_Series.index.to_list())
    FDR_in_Series.name="FDR_BH"
    return FDR_in_Series #this is FDR now

dg=pd.concat([dg,BH_FDR(dg["p_value"])],axis=1)

dg.to_csv(ks_test_file,sep='\t')



















