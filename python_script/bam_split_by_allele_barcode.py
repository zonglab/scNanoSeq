#!/usr/local/bin/python3.5

import sys
import os
import pysam



wd = sys.argv[1]
bam_in, bam_out_prefix = sys.argv[2], sys.argv[3]


os.chdir(wd)

samfile = pysam.AlignmentFile(bam_in, "rb")

out_dict={}
barcode_list=[]
for i in ["A",'C','G','T']:
	for j in ["A",'C','G','T']:
		for k in ["A",'C','G','T']:
			out_dict[i+j+k]=pysam.AlignmentFile(bam_out_prefix+"_"+i+j+k+".bam", template= samfile, mode= 'wb')
			




for read in samfile:
	try:
		mb = dict(read.tags)['mb']
		rb = dict(read.tags)['rb']
		if read.flag >=1024:
			read.flag=read.flag-1024
		if read.flag in [83,147]:
			t=rb; rb=mb; mb=t; del t
		new_tag = read.tags
		new_tag.append(('MR',mb+rb))
		read.set_tags(new_tag)


		out_dict[mb].write(read)

	except:
		pass
		
samfile.close()
for i in ["A",'C','G','T']:
	for j in ["A",'C','G','T']:
		for k in ["A",'C','G','T']:
			barcode_list.append(i+j+k)
			out_dict[i+j+k].close()
			




