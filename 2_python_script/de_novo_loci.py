import sys
import os
import pandas as pd

wd=sys.argv[1]
sample_list_file=sys.argv[2]
file_input_houzhui=sys.argv[3]
loci_output=sys.argv[4]


os.chdir(wd)


sample_list=[]
with open(sample_list_file,'r') as f:
  for line in f:
    sample_list.append(line.strip())


de_novo={}

for cellX in sample_list:
  fileX=cellX+file_input_houzhui
  with open(fileX,"r") as f:
    for line in f:
      if line.startswith("#"):
        continue
      lineX=line.strip().split("\t")
      chrX,posX,ref,alt=lineX[0],lineX[1],lineX[3],lineX[4]
      try:
        de_novo[("_").join([chrX,posX,ref,alt])]+=1
      except KeyError:
        de_novo[("_").join([chrX,posX,ref,alt])]=1



with open(loci_output,"w") as f:
  for tx in list(de_novo.keys()):
    x=tx.split("_")
    f.write(x[0]+"\t"+x[1]+"\n")
