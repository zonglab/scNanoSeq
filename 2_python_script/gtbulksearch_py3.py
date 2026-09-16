#! /usr/local/bin/python2.7

import csv
import sys
import pysam
import gzip

f1=sys.argv[1]
f2=sys.argv[2]
mxreads=sys.argv[3]

#samfile = pysam.AlignmentFile("/data1/czong/genome/X10_031617_analysis/Sample_HHHJ3ALXX-1-IDN702/mergedclone1bulk_RG.bam", "rb")
samfile = pysam.AlignmentFile(f2, "rb")

#print "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"

with open(f1, 'r') as m1:
#with gzip.open(f1, 'rb') as m1:
#    data1 = csv.reader(m1,delimiter='	')
#    mydict1 = {rows[1]:rows[4] for rows in data1}
  for line  in m1:
    if not line.startswith("#"):
       rows = line.strip().split("\t")
       #snp.append(int(rows[1]))
       #chrnum = int(rows[0][3:])
       chrnum = rows[0]
       pos =  int(rows[1])
       snpref= rows[3]
       snpvar= rows[4]

       snplen = len(snpvar)

       if snplen > 1:
         continue

       denovo = 0

       readnum = 0 
       for pileupcolumn in samfile.pileup(chrnum , pos-1, pos):
          for pileupread in pileupcolumn.pileups:
              if pileupcolumn.pos == pos-1:
                   #print "position", pileupread.query_position
                   if pileupread.query_position != None : 
                      if pileupread.alignment.query_qualities[pileupread.query_position] > 0 :                    
                       base = pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+snplen]
                       readnum += 1
                       #print pos, base, snpvar, snpref
                       if base == snpvar :
                          denovo += 1 
                     

       #print chrnum, pos, snpref, snpvar, denovo

       #for read in samfile.fetch( chrnum , pos-1, pos):
       #   if len(read.seq) > pos-1 - read.pos :
       #      base =  read.seq[pos-1 - read.pos]
             #print base 
       #      if base == snpvar :
       #            denovo = 0 
       #            break 
       if denovo <= int(mxreads) and readnum>=10 :
           #print  "denovo"
           print (line.strip())
           #print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(snpvcf[pos][0], snpvcf[pos][1],snpvcf[pos][2],snpvcf[pos][3], \
	   #           snpvcf[pos][4],snpvcf[pos][5],snpvcf[pos][6],snpvcf[pos][7]))

