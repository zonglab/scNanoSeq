import os
import sys
import numpy as np



wd=sys.argv[1]
File_hetero_detected_mut,File_hetero_detected_ss_var = sys.argv[2],sys.argv[3]
File_ss_variant_de_novo,File_de_novo = sys.argv[4],sys.argv[5]

bulk_hetero_number,damage_0_number,damage_swap = int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8])
Output_file_name=sys.argv[9]


def Union(lst1, lst2): 
    final_list = lst1 + lst2 
    return final_list 


def cumulative_count(vcf_file_name):
	position_count_XXX={} 

	with open(vcf_file_name,"r") as f1:
		for line in f1:
			x=line.strip().split("\t")
			position_str = x[11].replace("(","").replace(")","").replace(" ","").split(",")
			positionX= min(int(position_str[0]),int(position_str[1]))
			try:
				position_count_XXX[positionX]=position_count_XXX[positionX]+1
			except KeyError:
				position_count_XXX[positionX]=1

	XXX_keys=list(position_count_XXX.keys())
	XXX_keys.sort(reverse=True)

	cumulative_count_XXX={}
	cumulative_count_XXX[XXX_keys[0]]=position_count_XXX[XXX_keys[0]]

	for i in range(1,len(XXX_keys)):
		current_position=XXX_keys[i]
		previous_position=XXX_keys[i-1]
		cumulative_count_XXX[current_position]=cumulative_count_XXX[previous_position]+position_count_XXX[current_position]

	return cumulative_count_XXX


def dict_add_missing_key_values(dictXXX,startX,endX):
	for i in range(startX,endX+1,1):
		j=i
		while True:
			try:
				dictXXX[i]=dictXXX[j]
				break
			except KeyError:
				j+=1
	return dictXXX



os.chdir(wd)

cumulative_de_novo = cumulative_count(File_de_novo)
cumulative_ss_de_novo = cumulative_count(File_ss_variant_de_novo)
cumulative_hetero_mut = cumulative_count(File_hetero_detected_mut)
cumulative_hetero_ss = cumulative_count(File_hetero_detected_ss_var)

all_position_list=Union(Union(Union(list(cumulative_de_novo.keys()), list(cumulative_ss_de_novo.keys())),list(cumulative_hetero_mut.keys())),list(cumulative_hetero_ss.keys()))
min_pos_shared, max_pos_shared = min(all_position_list), max(all_position_list)

cumulative_de_novo[max_pos_shared+1]=0
cumulative_ss_de_novo[max_pos_shared+1]=0
cumulative_hetero_mut[max_pos_shared+1]=0
cumulative_hetero_ss[max_pos_shared+1]=0


cumulative_de_novo=dict_add_missing_key_values(cumulative_de_novo,min_pos_shared,max_pos_shared)
cumulative_ss_de_novo=dict_add_missing_key_values(cumulative_ss_de_novo,min_pos_shared,max_pos_shared)
cumulative_hetero_mut=dict_add_missing_key_values(cumulative_hetero_mut,min_pos_shared,max_pos_shared)
cumulative_hetero_ss=dict_add_missing_key_values(cumulative_hetero_ss,min_pos_shared,max_pos_shared)

swapping_rate=damage_swap/damage_0_number




with open(Output_file_name,"w") as f:
	line="\t".join(["Position","De_novo","De_novo_swap_corrected","De_novo_mutation_count","Bulk_mutation_detected"])+"\n"
	f.write(line)
	for i in range(min_pos_shared,max_pos_shared+1,1):
		positionX=i
		try:
			De_novo_number = cumulative_de_novo[i]/cumulative_hetero_mut[i]*bulk_hetero_number
			De_novo_swapping_corrected = (cumulative_de_novo[i]-swapping_rate*cumulative_ss_de_novo[i])/(cumulative_hetero_mut[i]-swapping_rate*cumulative_hetero_ss[i])*bulk_hetero_number
		except ZeroDivisionError:
			De_novo_number=np.nan
			De_novo_swapping_corrected=np.nan
		line="\t".join([str(positionX),str(De_novo_number),str(De_novo_swapping_corrected),str(cumulative_de_novo[i]),str(cumulative_hetero_mut[i])])+"\n"
		f.write(line)

