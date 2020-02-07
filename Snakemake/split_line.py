#!/usr/bin/env python3
#import commands
import os
import sys
import uuid
import subprocess


program=sys.argv[0]
arguments=sys.argv[1:]
count=len(arguments)

if count !=3:
	print ("Usage: python split_line.py fasta.fai interval_size output_dir(intervals/)")
	sys.exit(1)
elif count==3:
	program_name = sys.argv[0]
	input_fai=sys.argv[1]
	interval_size=int(sys.argv[2])
	output_dir_tmp=sys.argv[3]	

if output_dir_tmp.endswith('/'):
	output_dir=output_dir_tmp[:-1]
else:
	output_dir=output_dir_tmp


subprocess.run("vcfutils.pl splitchr -l "+str(interval_size)+" "+ input_fai+"| xargs -i echo \"{}\" > GRCh37.intervals", shell=True, check=True)

file=open("GRCh37.intervals")
i=1
for line in file:
	line=line.rstrip()
	fo=open(output_dir+"/"+str(i)+".intervals.list","w")
	print(line,file=fo)
	i=i+1
	fo.close()
file.close()

