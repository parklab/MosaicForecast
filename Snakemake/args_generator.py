#!/usr/bin/env python3
#import commands
import os
import sys
import glob
import re

program=sys.argv[0]
arguments=sys.argv[1:]
count=len(arguments)

if count <3:
	print ("python argv_generator.py interval_dir bam_links_dir output_dir")
	sys.exit(1)

interval_dir_tmp=sys.argv[1]
bam_links_dir_tmp=sys.argv[2]
output_dir_tmp=sys.argv[3]

if interval_dir_tmp.endswith('/'):
	interval_dir=interval_dir_tmp[:-1]
else:
	interval_dir=interval_dir_tmp

if bam_links_dir_tmp.endswith('/'):
	bam_links_dir=bam_links_dir_tmp[:-1]
else:
	bam_links_dir=bam_links_dir_tmp

if output_dir_tmp.endswith('/'):
	output_dir=output_dir_tmp[:-1]
else:
	output_dir=output_dir_tmp


for interval in glob.glob(os.path.join(interval_dir+"/","*.list")):
	# pattern=re.compile(r'([0-9]+)\.intervals\.list')
	pattern=re.compile(r'([0-9XY]+)\.intervals\.list')
	num=re.search(pattern, interval).group(1)
	fo=open(output_dir+"/"+str(num)+".args","w")
#	num=interval.str.extract(pattern)
	print(num)
	for bam in glob.glob(os.path.join(bam_links_dir+"/","*bam")):
		pattern=re.compile(r'\/([\w]+)\.bam')
		sample=re.search(pattern,bam).group(1)
		print(output_dir+"/"+sample+"_"+str(num)+"_pon.vcf.gz",file=fo)
	fo.close()

