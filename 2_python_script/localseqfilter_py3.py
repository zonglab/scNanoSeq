#!/usr/bin/env python

import sys

x1=sys.argv[1]
x2=sys.argv[2]

#filename1="denovo."+x1+".vcf"
#filename2="denovo."+x1+".localseq.txt"


filename1=x1
filename2=x2


with open(filename1) as f1, open(filename2) as f2: 
  for line, line2 in zip(f1, f2):

       rows = line.strip().split("\t")
       rows2 = line2.strip().split("\t")
       #snp.append(int(rows[1]))
       #chrnum = int(rows[0][3:])
       chrnum = rows[0]
       pos =  int(rows[1])
       snpref= rows[3]
       snpvar= rows[4]

       snplen = len(snpvar)

       if snplen > 1 :
           continue

       localseq = rows2[2][:10]
       localseq2 = rows2[2][10:]

       if ( localseq.count('A')+localseq.count('a') >=8 or  \
            localseq.count('T')+localseq.count('t') >=8 or  \
            localseq.count('G')+localseq.count('g') >=8 or  \
            localseq.count('C')+localseq.count('c') >=8 ) :
          continue

       if ( localseq2.count('A')+localseq2.count('a') >=8 or  \
            localseq2.count('T')+localseq2.count('t') >=8 or  \
            localseq2.count('G')+localseq2.count('g') >=8 or  \
            localseq2.count('C')+localseq2.count('c') >=8 ) :
          continue



       print (line.strip())
           #print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(snpvcf[pos][0], snpvcf[pos][1],snpvcf[pos][2],snpvcf[pos][3], \
	   #           snpvcf[pos][4],snpvcf[pos][5],snpvcf[pos][6],snpvcf[pos][7]))

