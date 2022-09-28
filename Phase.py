#!/usr/bin/env python3
import os
import sys
import uuid
import numpy as np
program=sys.argv[0]
arguments=sys.argv[1:]
count=len(arguments)
if count !=8:
	print ("Usage: python Phase.py bam_dir output_dir ref_fasta input_positions(file format:chr pos-1 pos ref alt sample, sep=\\t) min_dp_inforSNPs(int) Umap_mappability(bigWig file,k=24) n_threads_parallel sequencing_file_format(bam/cram)\n\nNote:\n1. Name of bam files should be \"sample.bam\" under the bam_dir, and there should be corresponding index files.\n2. There should be a fai file under the same dir of the fasta file (samtools faidx input.fa).\n3. The \"min_dp_inforSNPs\" is the minimum depth of coverage of trustworthy neaby het SNPs.\n4. Bam file is preferred than cram file, as the program would run much more slowly if using cram format.\n")
#	1       1004864 1004865 G       C       test
#1       13799   T       G       0/1     236     8       0.033   5       3       0.625   8447    282     146     90
	sys.exit(1)
elif count==8:
	#sample=sys.argv[1]
	bam_dir_tmp=sys.argv[1]
	output_dir_tmp=sys.argv[2]
	ref_fasta=sys.argv[3]
	input_pos=sys.argv[4]
	min_dp_inforSNPs=int(sys.argv[5])
	unimap_mappability_BigWigfile=sys.argv[6]	
	n_threads=int(sys.argv[7])
	seq_file_format=sys.argv[8]
	log_file='multiple_inforSNPs.log'

import regex as re
import pysam
from pyfaidx import Fasta
import scipy.stats
from multiprocessing import Pool
from subprocess import *
import subprocess
import pysamstats
#import os
from pandas import *
import itertools
from collections import OrderedDict, defaultdict



def run_cmd(cmd):
	Popen(cmd, shell=True, stdout=PIPE).communicate()

if bam_dir_tmp.endswith('/'):
	bam_dir=bam_dir_tmp[:-1]
else:
	bam_dir=bam_dir_tmp

if output_dir_tmp.endswith('/'):
	output_dir=output_dir_tmp[:-1]
else:
	output_dir=output_dir_tmp

run_cmd("mkdir -p "+output_dir)
run_cmd("mkdir -p "+output_dir+"/tmp")
#os.system("mkdir -p "+output_dir)
#os.system("mkdir -p "+output_dir+"/tmp")
homopolymers=list()
homopolymers=["AAAAA","TTTTT","GGGGG","CCCCC","ATATAT","TATATA","AGAGAG","GAGAGA","ACACAC","CACACA","TGTGTG","GTGTGT","GCGCGC","CGCGCG","CTCTCT","TCTCTC","ATTATT","TAATAA","AATAAT","GCCGCC","CGGCGG","CCGCCG","ATTTATTT","TAAATAAA","GCCCGCCCC","CGGGCGGG","CCGGCCGG","GGCCGGCC","TTTATTT","ATTTATT","TAAATAA","AAATAAT","GCCCGCC","CCCGCCC","GGCGGC","GAAAGAAA","AAAGAAAG","TTTCTTTC","TTCTTTCT","CCCTCCCT","CTTTCTTT"]

reference = Fasta(ref_fasta)
genome=ref_fasta+".fai"
file0=open(genome)
chr_sizes=dict()
for line in file0:
	line=line.rstrip()
	fields=line.split('\t')
	chr_sizes[fields[0]]=chr_sizes.get(fields[0],fields[1])
file0.close()


sites_chr_dict=dict()
sites_pos_dict=dict()
file=open(input_pos)

for line in file:
	line=line.rstrip()
	fields=line.split('\t')
	chr=fields[0]
	pos=int(fields[2])
	ref=fields[3]
	alt=fields[4]
	sample=fields[5]
	ID=sample+';'+chr+";"+str(pos)+";"+ref+";"+alt
	sites_chr_dict[ID]=chr
	sites_pos_dict[ID]=pos
#	print("chr"+chr,pos,pos+1,ID,file=input_mappability,sep="\t")
file.close()

tmp_filename=str(uuid.uuid4())
input_mappability=open(tmp_filename,'w')
for k,v in sorted(sites_chr_dict.items()):
	if not re.search('^chr',v):
		print("chr"+v,sites_pos_dict[k],int(sites_pos_dict[k])+1,k,file=input_mappability,sep="\t")
	else:
		print(v,sites_pos_dict[k],int(sites_pos_dict[k])+1,k,file=input_mappability,sep="\t")
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
##input_bam=bam_dir+"/"+sample+".bam"
##a=pysam.AlignmentFile(input_bam, "rb")
	
##head Vaccarino_brain.mutect2/final.passed.tsv
##1       2533604 G       A       Vaccarino_brain
cmd_list1=list()
cmd_list2=list()

##1st step: generate seperate bams supporting ref allale and alt alleles:
#file=open(input_pos)
#fo1=open(output_dir+"/all_candidates","w")
def process_line0(line):
#for line in file:
	#1       1004864 1004865 G       C       test
	line=line.rstrip()
	fields=line.split('\t')
	sample=fields[5]
	chr=fields[0]
	chrom=str(chr)
	major_allele=fields[3]
	minor_allele=fields[4]
	pos=int(fields[2])
	start=int(pos)-1
	end=int(pos)
	length=len(major_allele)-len(minor_allele)
	if length<0:
		length=0-length
	if not (re.search("MT",chrom) or re.search(",",minor_allele) or major_allele==minor_allele):
		try:
			if seq_file_format=="bam":
				#print ("ok here")
				input_bam=bam_dir+"/"+sample+".bam"
				bai_file=bam_dir+"/"+str(sample)+".bai"
				bai_file2=bam_dir+"/"+str(sample)+".bam.bai"
				if not os.path.exists(input_bam):
					print("no sample.bam under the bam_dir")
				if not os.path.exists(bai_file) and not os.path.exists(bai_file2):
					print("no bam index files under the bam_dir")
				a=pysam.AlignmentFile(input_bam, "rb", reference_filename=ref_fasta)
				f1=pysam.AlignmentFile(output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.bam","wb",template=a)
				#print("ok here")
				f2=pysam.AlignmentFile(output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.bam","wb",template=a)
				f3=pysam.AlignmentFile(output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.bam","wb",template=a)
			elif seq_file_format=="cram":
				input_cram=bam_dir+"/"+sample+".cram"
				crai_file=bam_dir+"/"+str(sample)+".crai"
				crai_file2=bam_dir+"/"+str(sample)+".cram.crai"
				if not os.path.exists(input_cram):
					print("no sample.cram under the cram_dir")
				if not os.path.exists(crai_file) and not os.path.exists(crai_file2):
					print("no cram index files under the cram_dir")
				a=pysam.AlignmentFile(input_cram, "rc", reference_filename=ref_fasta)
				f1=pysam.AlignmentFile(output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.cram","wc",template=a, reference_filename=ref_fasta)
				f2=pysam.AlignmentFile(output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.cram","wc",template=a, reference_filename=ref_fasta)
				f3=pysam.AlignmentFile(output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.cram","wc",template=a, reference_filename=ref_fasta)


				tmp1_localcram_filename=output_dir+"/tmp/"+sample+"_"+str(chr)+"_"+str(pos)+"_"+str(uuid.uuid4())+".cram"
				a_local=pysam.AlignmentFile(tmp1_localcram_filename,'wc',template=a,reference_filename=ref_fasta)
				for read in a.fetch(chrom,start-2001,end+2001):
					a_local.write(read)
				a_local.close()
				pysam.index(tmp1_localcram_filename,tmp1_localcram_filename+".crai")

			name=sample+'_'+chr+'_'+str(pos)
			major_ids=list()
			minor_ids=list()
			major_num=int(0)
			minor_num=int(0)
			
			if seq_file_format=="cram":
				a=pysam.AlignmentFile(tmp1_localcram_filename, "rc",reference_filename=ref_fasta)

			if len(major_allele)==1 and len(minor_allele)==1:# and minor_allele!=".":
				state="SNP"
				for pileupcolumn in a.pileup(chrom, start, end, max_depth=8000):
					for pileupread in pileupcolumn.pileups:
						if pileupread.indel!=0:
							continue
						try:
							querybase=pileupread.alignment.query_sequence[pileupread.query_position]
							if int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(major_allele): #and pileuperead.alignment.mapping_quality>=10:
								major_ids.append(pileupread.alignment.query_name)
								major_num+=1
							elif int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(minor_allele): #and pileupread.alignment.mapping_quality>=10:
								minor_ids.append(pileupread.alignment.query_name)
								minor_num+=1
						except:
							continue
	
			elif len(major_allele)>1 and len(major_allele)==len(minor_allele):
				state="MNP"
				for pileupcolumn in a.pileup(chrom, start, end, max_depth=8000):
					for pileupread in pileupcolumn.pileups:
						if pileupread.indel!=0:
							continue
						try:
							querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(minor_allele)]
							if int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(major_allele): #and pileuperead.alignment.mapping_quality>=10:
								major_ids.append(pileupread.alignment.query_name)
								major_num+=1
							elif int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(minor_allele): #and pileupread.alignment.mapping_quality>=10:
								minor_ids.append(pileupread.alignment.query_name)
								minor_num+=1
						except:
							continue
#			elif len(major_allele)>1 and len(minor_allele)==1:
			elif len(major_allele)> len(minor_allele):
				state="DEL"
				#context1[name]=reference[chrom][int(pos)-2:int(pos)+1]
				context1=reference[chrom][max(1,int(pos)-11):min(int(pos)+1,int(chr_sizes[chrom]))]
				context2=reference[chrom][max(1,int(pos)-1):min(int(pos)+10,int(chr_sizes[chrom]))]
				context=reference[chrom][max(1,int(pos)-11):min(int(pos)+10,int(chr_sizes[chrom]))]
	
				if_homopolymer="No"
				for item in homopolymers:
					if re.search(str(item), str(context1)) or re.search(str(item),str(context2)):
						if_homopolymer="Yes"
						break
				if if_homopolymer=="No":
					for read in a.fetch(chrom,start-length, end+length):
						try:
							#if read.cigar[0][0]==4 and read.cigar[0][1]<=length and read.reference_start>= pos-1 and read.reference_start-read.query_alignment_start< pos-1:
							if (read.cigar[0][0]==4 or read.cigar[0][0]==5) and read.reference_start>= pos-2 and read.reference_start-read.query_alignment_start< pos-1:
								query_clipped = read.query_sequence[:read.query_alignment_start][:length]
								if re.search(query_clipped, major_allele):
									minor_ids.append(read.query_name)
									minor_num+=1
							#elif read.cigar[-1][0]==4 and read.cigar[-1][1]<=length and read.reference_end <= pos-1 and (read.reference_end + read.query_length-read.query_alignment_end>pos-1):
							elif (read.cigar[-1][0]==4 or read.cigar[-1][0]==5) and read.reference_end <= pos and (read.reference_end + read.query_length-read.query_alignment_end>pos-1):
								query_clipped = read.query_sequence[read.query_alignment_end:][-length:]
								if re.search(query_clipped, major_allele):
									minor_ids.append(read.query_name)
									minor_num+=1
						except:
							continue	
							#print (chrom, pos, read.query_name, read.cigar)
					if seq_file_format=="cram":
						a=pysam.AlignmentFile(tmp1_localcram_filename, "rc",reference_filename=ref_fasta)
#					print("DEL here:ok")
					for pileupcolumn in a.pileup(chrom, start, end, max_depth=8000):
						for pileupread in pileupcolumn.pileups:
							try:
								if pileupread.indel<0: #del <0
									querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(minor_allele)]
									#print(querybase,major_allele, minor_allele)
									if int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(minor_allele):
										minor_ids.append(pileupread.alignment.query_name)
										minor_num+=1
								elif pileupread.indel==0:
									#if int(pileupcolumn.pos)==int(pos)-1 and (not (pileupread.alignment.query_alignment_end==int(pos)-1 and pileupread.alignment.cigar[-1][0]==4 and pileupread.alignment.cigar[-1][1]<=length)) and (not (pileupread.alignment.query_alignment_start==int(pos)-1 and pileupread.alignment.cigar[0][0]==4 and pileupread.alignment.cigar[0][1]<=length)): #and pileuperead.alignment.mapping_quality>=10:
									if int(pileupcolumn.pos)==int(pos)-1 and (not (pileupread.alignment.query_alignment_end<=int(pos) and (pileupread.alignment.cigar[-1][0]==4 or pileupread.alignment.cigar[-1][0]==5))) and (not (pileupread.alignment.query_alignment_start>=int(pos)-2 and (pileupread.alignment.cigar[0][0]==4 or pileupread.alignment.cigar[0][0]==5))): #and pileuperead.alignment.mapping_quality>=10:
										major_ids.append(pileupread.alignment.query_name)
										major_num+=1
							except:
								continue						
	
			#elif len(major_allele)==1 and len(minor_allele)>1:
			elif len(major_allele) < len(minor_allele):
				state="INS"
				context1=reference[chrom][max(1,int(pos)-11):min(int(pos)+1,int(chr_sizes[chrom]))]
				context2=reference[chrom][max(1,int(pos)-1):min(int(pos)+10,int(chr_sizes[chrom]))]
				context=reference[chrom][max(1,int(pos)-11):min(int(pos)+10,int(chr_sizes[chrom]))]
	
				if_homopolymer="No"
				for item in homopolymers:
					if re.search(str(item), str(context1)) or re.search(str(item),str(context2) or re.search(str(item),minor_allele)):
						if_homopolymer="Yes"
						break
				if if_homopolymer=="No":
					for read in a.fetch(chrom,start-length, end+length):
						try:
							#if read.cigar[0][0]==4 and read.cigar[0][1]<=length and read.reference_start>= pos-1 and read.reference_start-read.query_alignment_start< pos-1:
							if (read.cigar[0][0]==4 or read.cigar[0][0]==5) and read.reference_start>= pos-2 and read.reference_start-read.query_alignment_start< pos-1:
								query_clipped = read.query_sequence[:read.query_alignment_start][:length]
								if re.search(query_clipped, minor_allele):
									minor_ids.append(read.query_name)
									minor_num+=1
							#elif read.cigar[-1][0]==4 and read.cigar[-1][1]<=length and read.reference_end <= pos-1 and (read.reference_end + read.query_length-read.query_alignment_end>pos-1):
							elif (read.cigar[-1][0]==4 or read.cigar[-1][0]==5) and read.reference_end <= pos and (read.reference_end + read.query_length-read.query_alignment_end>pos-1):
								query_clipped=read.query_sequence[read.query_alignment_end:][-length:]
								if re.search(query_clipped, minor_allele):
									minor_ids.append(read.query_name)
									minor_num+=1
						except:
							continue
	#						print (chrom, pos, read.query_name, read.cigar)
					if seq_file_format=="cram":
						a=pysam.AlignmentFile(tmp1_localcram_filename, "rc",reference_filename=ref_fasta)
#					print("INS here:ok")
					for pileupcolumn in a.pileup(chrom, start, end, max_depth=8000):
						for pileupread in pileupcolumn.pileups:
							try:
								if pileupread.indel > 0: ###ins>0
									#querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+(len(minor_allele)-len(major_allele))+1]
#									querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(minor_allele)]
									#if int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(minor_allele):
									if int(pileupcolumn.pos)==int(pos)-1:
										minor_ids.append(pileupread.alignment.query_name)
										minor_num+=1
								elif pileupread.indel==0:
									#querybase=pileupread.alignment.query_sequence[pileupread.query_position]
									querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(major_allele)]
									if int(pileupcolumn.pos)==int(pos)-1:
										#if str(querybase)==str(major_allele) and (not (pileupread.alignment.query_alignment_end==int(pos)-1 and pileupread.alignment.cigar[-1][0]==4 and pileupread.alignment.cigar[-1][1]<=length)) and (not (pileupread.alignment.query_alignment_start==int(pos)-1 and pileupread.alignment.cigar[0][0]==4 and pileupread.alignment.cigar[0][1]<=length)):
										if str(querybase)==str(major_allele) and (not (pileupread.alignment.query_alignment_end>=int(pos)-2 and (pileupread.alignment.cigar[-1][0]==4 or pileupread.alignment.cigar[-1][0]==5))) and (not (pileupread.alignment.query_alignment_start<=int(pos) and (pileupread.alignment.cigar[0][0]==4 or pileupread.alignment.cigar[0][0]==5))):
											major_ids.append(pileupread.alignment.query_name)
											major_num+=1
							except:
								continue						
					#print(major_num,minor_num)	
	
								
		
			start=max(int(pos)-1000,1)
			end=min(int(pos)+1000,int(chr_sizes[chrom]))
			conflictnum=0
			major_ids=list(set(major_ids))
			minor_ids=list(set(minor_ids))
			if len(major_ids)>=2 and len(minor_ids)>=2:
				conflict_reads=set(major_ids) & set(minor_ids)
				conflictnum=len(conflict_reads)
					
				for read in a.fetch(chrom, start, end):
					if ((read.query_name in major_ids) or (read.query_name in minor_ids)) and (read.query_name not in list(conflict_reads)):
						f3.write(read)
					if (read.query_name in major_ids) and (read.query_name not in list(conflict_reads)):
						f1.write(read)
					if (read.query_name in minor_ids) and (read.query_name not in list(conflict_reads)):
						f2.write(read)
			f1.close()
			f2.close() 
			f3.close()
			if seq_file_format=="bam":	
				f1_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.bam"
				f2_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.bam"
				f3_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.bam"
				
				f1_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.sorted.bam"
				f2_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.sorted.bam"
				f3_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.sorted.bam"
			
				run_cmd("samtools sort "+f1_name+" -o "+f1_sorted_name)
				run_cmd("samtools sort "+f2_name+" -o "+f2_sorted_name)
				run_cmd("samtools sort "+f3_name+" -o "+f3_sorted_name)
			
				run_cmd("samtools index "+f1_sorted_name)	
				run_cmd("samtools index "+f2_sorted_name)	
				run_cmd("samtools index "+f3_sorted_name)
			elif seq_file_format=="cram":
				f1_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.cram"
				f2_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.cram"
				f3_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.cram"
				
				f1_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.sorted.cram"
				f2_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.sorted.cram"
				f3_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.sorted.cram"
			
				run_cmd("samtools sort "+f1_name+" -o "+f1_sorted_name)
				run_cmd("samtools sort "+f2_name+" -o "+f2_sorted_name)
				run_cmd("samtools sort "+f3_name+" -o "+f3_sorted_name)
			
				run_cmd("samtools index "+f1_sorted_name)	
				run_cmd("samtools index "+f2_sorted_name)	
				run_cmd("samtools index "+f3_sorted_name)
					
		
		#	cmd_list1.append("samtools sort "+f1_name+" -o "+f1_sorted_name)
		#	cmd_list1.append("samtools sort "+f2_name+" -o "+f2_sorted_name)
		#	cmd_list1.append("samtools sort "+f3_name+" -o "+f3_sorted_name)
		#	cmd_list2.append("samtools index "+f1_sorted_name)
		#	cmd_list2.append("samtools index "+f2_sorted_name)
		#	cmd_list2.append("samtools index "+f3_sorted_name)
		#	print (sample,chr,pos,major_allele,minor_allele,major_num,minor_num,conflictnum,file=fo1)
			return sample,chr,pos,major_allele,minor_allele,str(major_num),str(minor_num),str(conflictnum),state
		except:
			print("check the format of your input file:\nchr\tpos-1\tpos\tref\talt\tsample")

fo1=open(output_dir+"/all_candidates","w")

if __name__ == "__main__":
	pool = Pool(processes=int(n_threads))
	with open(input_pos) as source_file:
	# chunk the work into batches of 4 lines at a time
		#pool.map(process_line, source_file,1)
		result=pool.map(process_line0, source_file,1)
		for atuple in result:
			try:
				print (' '.join(str(x) for x in atuple),file=fo1)
			except:
				continue

fo1.close()
	
#pool=Pool(processes=int(n_threads))
#pool.map(run_cmd,cmd_list1,1)
#pool.close()
#pool.join()
#
#pool=Pool(processes=int(n_threads))
#pool.map(run_cmd,cmd_list2,1)
#pool.close()
#pool.join()


#2nd step: extract candidate nearby inforSNPs (germline het):
conflict_mosaic=dict()
variant_state=dict()
file=open(output_dir+"/all_candidates")
for line in file:
	line=line.rstrip()
	fields=line.split(' ')
	sample=fields[0]
	chr=fields[1]
	chrom=str(chr)
	pos=int(fields[2])
	major_allele=fields[3]
	minor_allele=fields[4]
	major_num=int(fields[5])
	minor_num=int(fields[6])
	conflictnum=int(fields[7])
	variant_type=fields[8]
	mosaic_name=sample+';'+chr+';'+str(pos)+";"+major_allele+";"+minor_allele
	conflict_mosaic[mosaic_name]=conflictnum
	variant_state[mosaic_name]=variant_type
file.close()

#==> /n/data1/hms/dbmi/park/yanmei/simulated_bams_na12878/12878-50x.merged.inforSNPs.pos <==
#1 18089 G T 1 17385 G A 0
def process_line(line):
#for line in file:
	line=line.rstrip()
	fields=line.split(' ')
	sample=fields[0]
	chr=fields[1]
	chrom=str(chr)
	pos=int(fields[2])
	major_allele=fields[3]
	minor_allele=fields[4]
	major_num=int(fields[5])
	minor_num=int(fields[6])
	conflictnum=int(fields[7])
	variant_type=fields[8]
	mosaic_name=sample+';'+chr+';'+str(pos)+";"+major_allele+";"+minor_allele
	conflict_mosaic[mosaic_name]=conflictnum
	start=int(pos)-1
	end=int(pos)
	if seq_file_format=="bam":
		f3_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.sorted.bam"
	elif seq_file_format=="cram":
		f3_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.sorted.cram"

	x=list()
	if major_num>=3 and minor_num>=3:
		print (chr,pos, major_allele,minor_allele,major_num, minor_num)
		##print (f3_sorted_name)
		if seq_file_format=="bam":
			f3_alignment_file =pysam.AlignmentFile(f3_sorted_name,'rb',reference_filename=ref_fasta)
		elif seq_file_format=="cram":
			f3_alignment_file =pysam.AlignmentFile(f3_sorted_name,'rc',reference_filename=ref_fasta)
		for rec in pysamstats.stat_variation(f3_alignment_file, fafile=ref_fasta, min_mapq=20, min_baseq=20):
			if  ([rec['A'],rec['C'],rec['G'],rec['T']].count(0)<=2) and (rec['reads_pp']>10):
				bases=['A','C','G','T']
				counts=[rec['A'],rec['C'],rec['G'],rec['T']]
				max_index=counts.index(max(counts))
				max_base=bases[max_index]
				max_num=int(rec[max_base])
				subset=list(set(['A','C','G','T'])-set(max_base))
				counts=[rec[subset[0]],rec[subset[1]],rec[subset[2]]]
				max_index_2nd = counts.index(max(counts))
				alt_base=subset[max_index_2nd]
				alt_num=int(rec[alt_base])
				#print (max_base, max_num, alt_base, alt_num)
				if alt_num>=int(3) and max_num>=int(3) and max_num+alt_num>=int(min_dp_inforSNPs) and scipy.stats.binom_test(alt_num, max_num+alt_num)>0.05 and rec['pos']!=int(pos):
					x.append((sample,chr,pos,major_allele,minor_allele,chr,int(rec['pos'])+1,max_base,alt_base,conflictnum,variant_type))
					#return sample,chr,pos,major_allele,minor_allele,chr,str(int(rec['pos'])+1),max_base,alt_base,str(conflictnum)
		print(x)
		return x 	
#file.close()
#merged_inforSNPs.close()
#file=open(output_dir+"/all_candidates")
#data = Parallel(n_jobs=int(n_threads))(delayed(process_line)(line)
#                           for line in open(output_dir+"/all_candidates"))
#https://github.com/bioconda/bioconda-recipes/issues/12100
merged_inforSNPs=open(output_dir+"/all.merged.inforSNPs.pos","w")

if __name__ == "__main__":
	pool = Pool(processes=int(n_threads))
	with open(output_dir+"/all_candidates") as source_file:
	# chunk the work into batches of 4 lines at a time
		#pool.map(process_line, source_file,1)
		result=pool.map(process_line, source_file,1)
		if len(result)>0:
			for item in result:
				if not item is None:
					for atuple in item:
						try:
							print (' '.join(str(x) for x in atuple),file=merged_inforSNPs)
						except:
							continue

merged_inforSNPs.close()


##3rd step: generate the 2X2 table:
n_major_het1=dict()
n_major_het2=dict()
n_minor_het1=dict()
n_minor_het2=dict()



def process_line2(line):
#input_pos=output_dir+"/all.merged.inforSNPs.pos"
#1 1661246 T C 1 1660812 A G 0
#1 2585130 C G 1 2584706 A G 4
#file=open(input_pos)
#for line in file:
	line=line.rstrip()
	fields=line.split(' ')
	sample=fields[0]
	chr=fields[1]
	chrom=str(chr)
	pos=fields[2]
	major_allele=fields[3]
	minor_allele=fields[4]
	inforSNP_pos=fields[6]
	inforSNP_ref=fields[7]
	inforSNP_alt=fields[8]
	conflict=fields[9]
	variant_type=fields[10]
	if pos != inforSNP_pos:
		if seq_file_format=="bam":
			f1_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.sorted.bam"
			f2_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.sorted.bam"
			a1=pysam.AlignmentFile(f1_sorted_name,"rb", reference_filename=ref_fasta)
			a2=pysam.AlignmentFile(f2_sorted_name,"rb", reference_filename=ref_fasta)
		elif seq_file_format=="cram":
			f1_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.sorted.cram"
			f2_sorted_name=output_dir+"/tmp/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.sorted.cram"
			a1=pysam.AlignmentFile(f1_sorted_name,"rc", reference_filename=ref_fasta)
			a2=pysam.AlignmentFile(f2_sorted_name,"rc", reference_filename=ref_fasta)
		start_pos=max(int(inforSNP_pos)-1000,0)
		end_pos=min(int(inforSNP_pos)+1000,int(chr_sizes[chrom]))
	
		name=sample+";"+chr+';'+str(pos)+';'+major_allele+";"+minor_allele+";"+str(inforSNP_pos)+";"+str(inforSNP_ref)+";"+str(inforSNP_alt)+";"+str(conflict)
		n_major_het1[name]=0
		n_major_het2[name]=0
		n_minor_het1[name]=0
		n_minor_het2[name]=0
	
		for rec in pysamstats.stat_variation(a1,fafile=ref_fasta,chrom=chr,start=start_pos,end=end_pos):
			if rec['pos']==int(inforSNP_pos)-1:
				n_major_het1[name]=rec[inforSNP_ref]
				n_major_het2[name]=rec[inforSNP_alt]
		for rec in pysamstats.stat_variation(a2,fafile=ref_fasta,chrom=chr,start=start_pos,end=end_pos):
			if rec['pos']==int(inforSNP_pos)-1:
				n_minor_het1[name]=rec[inforSNP_ref]
				n_minor_het2[name]=rec[inforSNP_alt]
#file.close()
		return sample, chr,str(pos),major_allele, minor_allele, str(inforSNP_pos), str(inforSNP_ref), str(inforSNP_alt), str(conflict),  str(n_major_het1[name]), str(n_major_het2[name]),str(n_minor_het1[name]),str(n_minor_het2[name]), variant_type


fo2=open(output_dir+"/all_2x2table","w")
#print ("sample","chr","pos","ref","alt","pos_inforSNP","het1","het2","conflicted_reads","ref_het1_count","ref_het2_count","alt_het1_count","alt_het2_count",file=fo2)
#data = Parallel(n_jobs=int(n_threads))(delayed(process_line2)(line)
#                           for line in open(output_dir+"/all.merged.inforSNPs.pos"))
#fo2.close()
#12878-200x 10_10002280_A_G_10002816_G_A_0 14 0 0 8
#for k,v in sorted(n_major_het1.items()):
#	#print >>fo, sample,k,n_major_het1[k],n_major_het2[k],n_minor_het1[k],n_minor_het2[k]
#	print (' '.join(str(x) for x in k.split(";")), n_major_het1[k],n_major_het2[k],n_minor_het1[k],n_minor_het2[k],file=fo2)
#fo2.close()

if __name__ == "__main__":
        pool = Pool(processes=int(n_threads))
        with open(output_dir+"/all.merged.inforSNPs.pos") as source_file:
        # chunk the work into batches of 4 lines at a time
                #pool.map(process_line, source_file,1)
                result=pool.map(process_line2, source_file,1)
                for atuple in result:
                        try:
                                print (' '.join(str(x) for x in atuple),file=fo2)
                        except:
                                continue
fo2.close()

##last step: assign phasing state to each site
phase=defaultdict(dict)
inforSNPs=defaultdict(list)
input_file=output_dir+"/all_2x2table"
#12878-250x 10 100022092 A G 100021983 G A 0 21 19 0 6
input_table =open(input_file)
for line in input_table:
	line=line.rstrip()
	if not re.search('sample',line):
		fields=line.split(' ')
		sample=fields[0]
		chr=fields[1]
		pos=fields[2]
		ref=fields[3]
		alt=fields[4]
		pos_inforSNP=fields[5]
		het1_inforSNP=fields[6]
		het2_inforSNP=fields[7]
		conflict=int(fields[8])
		C1=float(fields[9])
		C2=float(fields[10])
		C3=float(fields[11])
		C4=float(fields[12])
		variant_type=fields[13]
		name=sample+";"+chr+";"+pos+";"+ref+";"+alt
		inforSNPs[name].append(sample+";"+chr+";"+pos_inforSNP+";"+het1_inforSNP+";"+het2_inforSNP)
	
		phase[name]['hap=2']=phase[name].get("hap=2",0)+0
		phase[name]['hap=3']=phase[name].get("hap=3",0)+0
		phase[name]['hap>3']=phase[name].get("hap>3",0)+0
		#phase[name]['NA']=phase[name].get("NA",0)+0
		if C1+C2+C3+C4>=10:
			if variant_type=="SNP" or variant_type=="MNP":
			#UMB1349 18 30348609 C T 30348301 G A 0 34 29 0 2
				#if not ( ((C1>C2*10) and (C4>C3*10)) or ((C1<C2/10) and (C4<C3/10)) or (((C1>C2/10) and (C1<C2*10)) and (C3<=C4/10 or C4<=C3/10))  ):
				if not ( ( (C1>C2*10) and (C4>C3*10) and (C1<C4*5 and C1>C4/5) ) or ((C1<C2/10) and (C4<C3/10) and (C2<C3*5 and C2>C3/5)) or (((C1>C2/10) and (C1<C2*10)) and (C3<=C4/10 or C4<=C3/10))  ):
					phase[name]['hap>3']=phase[name].get("hap>3",0)+1
				elif (((C1>C2/10) and (C1<C2*10)) and (C3<C4/10 or C4<C3/10) and (C3+C4>1)):
					phase[name]['hap=3']=phase[name].get("hap=3",0)+1
				#elif ((C1>C2*10) and (C4>C3*10)) or ((C1<C2/10) and (C4<C3/10)):
				#	phase[name]['hap=2']=phase[name].get("hap=2",0)+1
				elif ( (C1>C2*10) and (C4>C3*10) and (C1<C4*5 and C1>C4/5) ) or ((C1<C2/10) and (C4<C3/10) and (C2<C3*5 and C2>C3/5)):
					phase[name]['hap=2']=phase[name].get("hap=2",0)+1
		#	elif C3+C4==0:
		#		phase[name]['NA']=phase[name].get("NA",0)+1
			elif variant_type!="SNP":
				if not ( ((C1>C2*5) and (C4>C3*5) and (C1<C4*5 and C1>C4/5)) or ((C1<C2/5) and (C4<C3/5) and (C2<C3*5 and C2>C3/5)) or (((C1>C2/5) and (C1<C2*5)) and (C3<=C4/5 or C4<=C3/5))  ):
					phase[name]['hap>3']=phase[name].get("hap>3",0)+1
				elif (((C1>C2/5) and (C1<C2*5)) and (C3<C4/5 or C4<C3/5) and (C3+C4>1)):
					phase[name]['hap=3']=phase[name].get("hap=3",0)+1
				elif ((C1>C2*5) and (C4>C3*5) and (C1<C4*5 and C1>C4/5)) or ((C1<C2/5) and (C4<C3/5) and (C2<C3*5 and C2>C3/5)):
					phase[name]['hap=2']=phase[name].get("hap=2",0)+1

phasing_2by2 = dict()
fo=open(output_dir+"/all.phasing_2by2","w")
for k,v in sorted(phase.items()):
	if max(phase[k]['hap=2'],phase[k]['hap=3'],phase[k]['hap>3'])>0:
		phasing_2by2[k]=max(v,key=v.get)
		if phase[k]['hap>3']>=2:
			phasing_2by2[k]='hap>3'
		if phasing_2by2[k]=="hap=3" and phase[k]['hap>3']==phase[k]['hap=3']:
			phasing_2by2[k]='hap>3'
		print (' '.join(str(x) for x in k.split(";")), v,max(v, key=v.get), file=fo)
	elif max(phase[k]['hap=2'],phase[k]['hap=3'],phase[k]['hap>3'])==0:
		phasing_2by2[k]="UNKNOWN"
fo.close()

#phasing_2by2 = dict()
#for k,v in sorted(phase.items()):
#	phasing_2by2[k]=max(v,key=v.get)

#12878-250x 10 10002280 A G 10002816 G A 0 18 0 0 1


##one additional step: using multiple germline het sites to do further phasing
MT2_inforSNPs_phasing=defaultdict(list)
MT2_phasing_num=defaultdict(dict)
for k,v in sorted(inforSNPs.items()):
	reads_type1=list()
	reads_type2=list()
	reads_type3=list()
	reads_type4=list()
	if len(v)>1:
		sample=k.split(';')[0]
		chr=k.split(';')[1]
		mosaic_pos=k.split(';')[2]
		inforSNPs_pos_list=list()
		inforSNPs_alleles_list=list()
		mosaic_name=k
		MT2_phasing_num[mosaic_name]['correct']=MT2_phasing_num[mosaic_name].get('correct',0)
		MT2_phasing_num[mosaic_name]['wrong']=MT2_phasing_num[mosaic_name].get('wrong',0)
		MT2_phasing_num[mosaic_name]['doubt']=MT2_phasing_num[mosaic_name].get('doubt',0)
		for i in range(0,len(v)):
			inforSNPs_pos_list.append(int(inforSNPs[k][i].split(';')[2]))
			inforSNPs_alleles_list.append([inforSNPs[k][i].split(';')[3],inforSNPs[k][i].split(';')[4]])
		if seq_file_format=="bam":
			samfile=pysam.AlignmentFile(bam_dir+"/"+sample+".bam", "rb",reference_filename=ref_fasta)
		elif seq_file_format=="cram":
			samfile=pysam.AlignmentFile(bam_dir+"/"+sample+".bam", "rc",reference_filename=ref_fasta)
		M=defaultdict(dict)
		for read in samfile.fetch(chr, min(inforSNPs_pos_list),max(inforSNPs_pos_list)):
			readID=read.query_name
			for i in range(0,len(v)):
				M[readID][str(i)]=M[readID].get(str(i),".")
				try:
					if (int(inforSNPs[k][i].split(';')[2])-1 > int(read.reference_start)) and (int(inforSNPs[k][i].split(';')[2])-1 <int(read.reference_end)):
						distance=int(inforSNPs[k][i].split(';')[2])-int(read.reference_start)
						cigar=str(read.cigarstring)
						cigar_num=re.split('M|N|D|S|I|H',cigar)
						cigar_patt=re.split('[0-9]*',cigar)
						sum=int(0)
						offset=int(0)
						for l in range(1,len(cigar_patt)):
							if cigar_patt[l]=="M" or cigar_patt[l]=="S" or cigar_patt[l]=="H":
								##print l
								sum+=int(cigar_num[l-1])
			                                        ##print cigar_patt[l],sum
							elif cigar_patt[l]=="I":
								sum+=int(cigar_num[l-1])
								offset=offset-int(cigar_num[l-1])
			                                        #print cigar_patt[l],sum
							elif cigar_patt[l]=="N" or cigar_patt[l]=="D":
								offset+=int(cigar_num[l-1])
			                                        #print cigar_patt[l],sum	
							if sum>=distance:
								break
						distance=distance-int(offset)
						querybase=read.query_sequence[distance-1]
						baseQ=read.query_qualities[distance-1]
						if baseQ>=20:
							M[readID][str(i)]=str(querybase)
				except:
					continue
		df=DataFrame(M).T
		##sort the column index according to values:
		try:
			df=df.reindex_axis(sorted(df.columns,key=lambda x:int(x[0:])), axis=1)
		except:
			print ("here1")
		#df=df.reindex_axis(sorted(df.columns,key=lambda x:int(x[0:])), axis=1)
		##print df
		try:
			##python2:
			##df2 =DataFrame(df.apply(lambda x: ''.join(x[range(0,len(df.columns))]),axis=1))
			##python3:	
			df2 =DataFrame(df.apply(lambda x: ''.join(x),axis=1))
		except:
			print (df)
		df2.columns = ['p0']
		str_lenth=len(df.columns)	
		all_haps=list()
		phasing_probs=dict()
		for m in itertools.combinations(range(0,str_lenth),2):
			if m[1]-m[0]>3:
				continue
			pos1_alleles=inforSNPs_alleles_list[m[0]]
			pos2_alleles=inforSNPs_alleles_list[m[1]]
			string_for_search="^"
			if m[0] >=1:
				for n in range(0,m[0]):
					string_for_search=string_for_search+"."

			string1_for_search=string_for_search+pos1_alleles[0]
			string2_for_search=string_for_search+pos1_alleles[0]
			string3_for_search=string_for_search+pos1_alleles[1]
			string4_for_search=string_for_search+pos1_alleles[1]
			
			string_for_search=""
			if m[1]-m[0]>1 and m[1]-m[0]<=3:
				for n in range(m[0],m[1]):
					string_for_search=string_for_search+"."

			string1_for_search=string1_for_search+string_for_search+pos2_alleles[0]
			string2_for_search=string2_for_search+string_for_search+pos2_alleles[1]
			string3_for_search=string3_for_search+string_for_search+pos2_alleles[0]
			string4_for_search=string4_for_search+string_for_search+pos2_alleles[1]
			
##			print mosaic_name, string1_for_search,string2_for_search,string3_for_search, string4_for_search, np.sum(df2.p0.str.contains(string1_for_search)), np.sum(df2.p0.str.contains(string2_for_search)), np.sum(df2.p0.str.contains(string3_for_search)), np.sum(df2.p0.str.contains(string4_for_search))

			if np.sum(df2.p0.str.contains(string1_for_search)) + np.sum(df2.p0.str.contains(string2_for_search))+ np.sum(df2.p0.str.contains(string3_for_search)) + np.sum(df2.p0.str.contains(string4_for_search)) < 10:
				MT2_inforSNPs_phasing[mosaic_name].append("lack"+"_"+str(np.sum(df2.p0.str.contains(string1_for_search)))+"_"+str(np.sum(df2.p0.str.contains(string2_for_search)))+"_" + str(np.sum(df2.p0.str.contains(string3_for_search))) + "_" + str(np.sum(df2.p0.str.contains(string4_for_search))) )


			elif ((np.sum(df2.p0.str.contains(string1_for_search)) >0 and np.sum(df2.p0.str.contains(string4_for_search))>0) and (np.sum(df2.p0.str.contains(string2_for_search))==0 and np.sum(df2.p0.str.contains(string3_for_search))==0)) or ((np.sum(df2.p0.str.contains(string1_for_search)) ==0 and np.sum(df2.p0.str.contains(string4_for_search))==0) and (np.sum(df2.p0.str.contains(string2_for_search))>0 and np.sum(df2.p0.str.contains(string3_for_search))>0)):
				MT2_inforSNPs_phasing[mosaic_name].append("correct"+"_"+str(np.sum(df2.p0.str.contains(string1_for_search)))+"_"+str(np.sum(df2.p0.str.contains(string2_for_search)))+"_" + str(np.sum(df2.p0.str.contains(string3_for_search))) + "_" + str(np.sum(df2.p0.str.contains(string4_for_search))) )
				MT2_phasing_num[mosaic_name]['correct']=MT2_phasing_num[mosaic_name].get('correct',0)+1

			elif ((np.sum(df2.p0.str.contains(string1_for_search)) >0 and np.sum(df2.p0.str.contains(string4_for_search))>0) and (np.sum(df2.p0.str.contains(string2_for_search))< float(np.sum(df2.p0.str.contains(string1_for_search)))/float(10) and np.sum(df2.p0.str.contains(string3_for_search)) < float(np.sum(df2.p0.str.contains(string4_for_search)))/float(10) )) or ( (np.sum(df2.p0.str.contains(string2_for_search)) >0 and np.sum(df2.p0.str.contains(string3_for_search))>0) and ( np.sum(df2.p0.str.contains(string1_for_search))< float(np.sum(df2.p0.str.contains(string2_for_search)))/float(10) and np.sum(df2.p0.str.contains(string4_for_search)) < float(np.sum(df2.p0.str.contains(string3_for_search)))/float(10) ) ):
				MT2_inforSNPs_phasing[mosaic_name].append("doubt"+"_"+str(np.sum(df2.p0.str.contains(string1_for_search)))+"_"+str(np.sum(df2.p0.str.contains(string2_for_search)))+"_" + str(np.sum(df2.p0.str.contains(string3_for_search))) + "_" + str(np.sum(df2.p0.str.contains(string4_for_search))) )
				MT2_phasing_num[mosaic_name]['doubt']=MT2_phasing_num[mosaic_name].get('doubt',0)+1
			else:
				MT2_inforSNPs_phasing[mosaic_name].append("wrong"+"_"+str(np.sum(df2.p0.str.contains(string1_for_search)))+"_"+str(np.sum(df2.p0.str.contains(string2_for_search)))+"_" + str(np.sum(df2.p0.str.contains(string3_for_search))) + "_" + str(np.sum(df2.p0.str.contains(string4_for_search))) )
				MT2_phasing_num[mosaic_name]['wrong']=MT2_phasing_num[mosaic_name].get('wrong',0)+1

fo3=open(output_dir+"/"+log_file,'w')
for k,v in sorted(MT2_inforSNPs_phasing.items()):
	print (k,v,file=fo3)
fo3.close()

fo4=open(output_dir+"/all.phasing","w")
print ("sample","chr","pos","ref","alt","phasing","conflicting_reads","mappability","variant_type",file=fo4, sep="\t")
for k,v in sorted(phasing_2by2.items()):
	if not v == "UNKNOWN":
		if k in MT2_phasing_num:
			if MT2_phasing_num[k]['wrong']>=MT2_phasing_num[k]['correct'] and MT2_phasing_num[k]['wrong']>0:
	#		if max(MT2_phasing_num[k], key=MT2_phasing_num[k].get)=="wrong":
				v="hap>3"
		print ('\t'.join(k.split(';')), v, conflict_mosaic[k], mappability[k], variant_state[k], sep="\t", file=fo4)
fo4.close()


