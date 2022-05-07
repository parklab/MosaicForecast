#!/usr/bin/env python3
#import commands
import os
import sys
import uuid
import re
import subprocess
from subprocess import *
import gzip

def run_cmd(cmd):
	Popen(cmd, shell=True, stdout=PIPE).communicate()

program=sys.argv[0]
arguments=sys.argv[1:]
count=len(arguments)

if count !=3:
	print ("Usage: python(v3) MT2_filter.py sample_id input_vcf(Mutect2-PoN) repeat_bed(SegDup+clutered regions)\n")
	sys.exit(1)
elif count==3:
	sample=sys.argv[1]
	input_vcf=sys.argv[2]
	repeat_file=sys.argv[3]


if input_vcf.endswith('.vcf.gz'):
	file=gzip.open(input_vcf,'rt')
elif input_vcf.endswith('.vcf'):
	file=open(input_vcf)
else:
	print ("please check if the input is a Mutect2 vcf file")

tmp_filename=str(uuid.uuid4())
fo=open(tmp_filename,'w')
for line in file:
	line=line.rstrip()
	fields=line.split('\t')
	if not re.search("^#",line):
		chr=fields[0]
		pos=int(fields[1])
		rs=fields[2]
		ref=fields[3]
		alt=fields[4]
		filters=fields[6]
		INFO_name=fields[8].split(':')
		INFOs=fields[9].split(':')
		AF=float(INFOs[2])
		ref_count=int(INFOs[1].split(',')[0])
		alt_count=int(INFOs[1].split(',')[1])
		depth=ref_count+alt_count
		if not re.search('str_contraction', filters) and not re.search('t_lod_fstar', filters) and not re.search('triallelic_site', filters) and not re.search('panel_of_normals', filters) and not re.search('multiallelic',filters):
			if re.search(":0|1:", line) or re.search(":1|0:", line):
				if AF>=0.02 and AF<0.4 and alt_count>=2:
					print(chr,pos-1,pos,ref,alt,sample,depth, AF, file=fo,sep="\t")
			elif not re.search(":0|1:", line):
				if AF>=0.03 and AF<0.4 and alt_count>=2:
					print(chr,pos-1,pos,ref,alt,sample,depth, AF, file=fo,sep="\t")

fo.close()
my_cmd="subtractBed -a "+tmp_filename+" -b "+repeat_file+" > "+sample+".bed"
Popen(my_cmd, shell=True, stdout=PIPE).communicate()

subprocess.run("rm "+tmp_filename, shell=True)
