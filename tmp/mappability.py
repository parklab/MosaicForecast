#!/usr/bin/env python
#import commands
import os
import sys
import subprocess
import uuid

program=sys.argv[0]
arguments=sys.argv[1:]
count=len(arguments)
if count !=5:
	print ("Usage: python(v3) ReadLevel_Features_extraction.py input_bed(file_format: chr pos-1 pos ref alt sample, sep=\"\\t\") output_features bam_dir reference_fasta Umap_mappability(bigWig file,k=24) \n\nNote: 1. name of bam files are \"sample.bam\" by default. \n2.there should be a fai file under the same dir of the fasta file (samtools faidx input.fa) \n3. we did not use dbSNP AF as an feature, but you can use it to train your model if you have interest in common variants\n4.The program to extract mappability score: \"bigWigAverageOverBed\" could be downloaded here at http://hgdownload.soe.ucsc.edu/admin/exe/, the program to convert wiggle file to BigWig file \"wigToBigWig\", and the \"fetchChromSizes\" script to create the chrom.sizes file for the UCSC database with which you are working (e.g., hg19) could be downloaded from the same directory. The wiggle file containing mappability score (Umap,k=24) could be downloaded here: https://bismap.hoffmanlab.org/")
	sys.exit(1)
elif count==5:
	program_name = sys.argv[0]
	input_pos=sys.argv[1] #walsh.nocluster.noalt_allele_in_normal.norepeat.bed 1       1015256 1015257 A       G       Walsh
	output=sys.argv[2]
	bam_dir_tmp=sys.argv[3]
	reference_fasta=sys.argv[4]
	unimap_mappability_BigWigfile=sys.argv[5]
	#sequencing_type=sys.argv[5]

import numpy as np
import regex as re
from collections import defaultdict
from pyfaidx import Fasta
import pysam
from scipy.stats import mannwhitneyu


file=open(input_pos)
tmp_filename=str(uuid.uuid4())
input_mappability=open(tmp_filename,'w')
for line in file:
	line=line.rstrip()
	fields=line.split('\t')
	chr=fields[0]
	pos=int(fields[2])
	ref=fields[3]
	alt=fields[4]
	sample=fields[5]
	ID=sample+'~'+chr+"~"+str(pos)+"~"+ref+"~"+alt
	print("chr"+chr,pos,pos+1,ID,file=input_mappability,sep="\t")
file.close()
input_mappability.close()

#bigWigAverageOverBed /n/data1/hms/dbmi/park/yanmei/resources/hg19/k24.umap.wg.bw test.bed test.tab
subprocess.run("bigWigAverageOverBed "+unimap_mappability_BigWigfile+" "+tmp_filename+" "+tmp_filename+".2", shell=True, check=True)

mappability=dict()
file=open(tmp_filename+".2")
#4638~1~86002185~G~A     1       0       0       0       0
#4638~10~39135403~A~T    1       1       0.5     0.5     0.5
#4638~16~21579992~C~T    1       1       0.958333        0.958333        0.958333
for line in file:
	line=line.rstrip()
	fields=line.split('\t')
	ID=fields[0]
	score=float(fields[5])
	mappability[ID]=score
file.close()
subprocess.run("rm "+tmp_filename, shell=True)
subprocess.run("rm "+tmp_filename+".2", shell=True)

	
