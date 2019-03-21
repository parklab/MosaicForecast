#!/usr/bin/env python
#import commands
import os
import sys
import uuid

program=sys.argv[0]
arguments=sys.argv[1:]
count=len(arguments)

if count !=6:
	print ("Usage: python(v3) ReadLevel_Features_extraction.py input_bed(file_format: chr pos-1 pos ref alt sample, sep=\"\\t\") output_features bam_dir reference_fasta Umap_mappability(bigWig file,k=24) num_threads_parallel\n\nNote:\n1. Names of bam files should be \"sample.bam\" under the bam_dir. \n\n2. There should be a fai file under the same dir of the fasta file (samtools faidx input.fa) \n\n3. We did not use gnomad population AF as an feature (instead we use it to filter), but you can use it to train your model if you have interest in common variants\n\n4. The program to extract mappability score: \"bigWigAverageOverBed\" could be downloaded here at http://hgdownload.soe.ucsc.edu/admin/exe/, the program to convert wiggle file to BigWig file \"wigToBigWig\", and the \"fetchChromSizes\" script to create the chrom.sizes file for the UCSC database with which you are working (e.g., hg19) could be downloaded from the same directory. The wiggle file containing mappability score (Umap,k=24) could be downloaded here: https://bismap.hoffmanlab.org/\n")
	sys.exit(1)
elif count==6:
	program_name = sys.argv[0]
	input_pos=sys.argv[1] #walsh.nocluster.noalt_allele_in_normal.norepeat.bed 1       1015256 1015257 A       G       Walsh
	output=sys.argv[2]
	bam_dir_tmp=sys.argv[3]
	reference_fasta=sys.argv[4]
	unimap_mappability_BigWigfile=sys.argv[5]
	n_jobs=sys.argv[6]
	#sequencing_type=sys.argv[5]

import numpy as np
import regex as re
from collections import defaultdict
from pyfaidx import Fasta
import pysam
from scipy.stats import mannwhitneyu
from subprocess import *
from multiprocessing import Pool
import subprocess
base=dict()
base['A']='T'
base['T']='A'
base['G']='C'
base['C']='G'
#base['N']='N'

homopolymers=list()
homopolymers=["AAAAA","TTTTT","GGGGG","CCCCC","ATATAT","TATATA","AGAGAG","GAGAGA","ACACAC","CACACA","TGTGTG","GTGTGT","GCGCGC","CGCGCG","CTCTCT","TCTCTC","ATTATT","TAATAA","AATAAT","GCCGCC","CGGCGG","CCGCCG","ATTTATTT","TAAATAAA","GCCCGCCCC","CGGGCGGG","CCGGCCGG","GGCCGGCC","TTTATTT","ATTTATT","TAAATAA","AAATAAT","GCCCGCC","CCCGCCC","GGCGGC","GAAAGAAA","AAAGAAAG","TTTCTTTC","TTCTTTCT","CCCTCCCT","CTTTCTTT"]

#file=open(input_pos)
#tmp_filename=str(uuid.uuid4())
#input_mappability=open(tmp_filename,'w')
#for line in file:
#	line=line.rstrip()
#	fields=line.split('\t')
#	chr=fields[0]
#	pos=int(fields[2])
#	ref=fields[3]
#	alt=fields[4]
#	sample=fields[5]
#	ID=sample+'~'+chr+"~"+str(pos)+"~"+ref+"~"+alt
#	print("chr"+chr,pos,pos+1,ID,file=input_mappability,sep="\t")
#file.close()
#input_mappability.close()

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
	ID=sample+'~'+chr+"~"+str(pos)+"~"+ref+"~"+alt
	sites_chr_dict[ID]=chr
	sites_pos_dict[ID]=pos
#	print("chr"+chr,pos,pos+1,ID,file=input_mappability,sep="\t")
file.close()

tmp_filename=str(uuid.uuid4())
input_mappability=open(tmp_filename,'w')
for k,v in sorted(sites_chr_dict.items()):
	print("chr"+v,sites_pos_dict[k],int(sites_pos_dict[k])+1,k,file=input_mappability,sep="\t")
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


fo=open(output,"w")
header='id querypos_major querypos_minor leftpos_major leftpos_minor seqpos_major seqpos_minor mapq_major mapq_minor baseq_major baseq_minor baseq_major_near1b baseq_minor_near1b major_plus major_minus minor_plus minor_minus context1 context2 context1_count context2_count mismatches_major mismatches_minor major_read1 major_read2 minor_read1 minor_read2 dp_near dp_far dp_p conflict_num mappability type length GCcontent ref_softclip alt_softclip'.split()
print (' '.join(header),file=fo)
#print (' '.join(header))

if bam_dir_tmp.endswith('/'):
	bam_dir=bam_dir_tmp[:-1]
else:
	bam_dir=bam_dir_tmp

reference = Fasta(reference_fasta)
reference_fai =reference_fasta+".fai"
genome=reference_fai

querypos_major=defaultdict(list)
mapq_major=defaultdict(list)
baseq_major=defaultdict(list)
leftpos_major=defaultdict(list)
mismatches_major=defaultdict(list)
major_plus=dict()
major_minus=dict()
major_read1=dict()
major_read2=dict()
major_ids=defaultdict(list)
minor_ids=defaultdict(list)
conflict_num=dict()
baseq_major_near1b=defaultdict(list)
seqpos_major=defaultdict(list)
context1=dict()
context2=dict()
context1_count=dict()
context2_count=dict()
querypos_minor=defaultdict(list)
dp_far=defaultdict(list)
dp_near=defaultdict(list)
mapq_minor=defaultdict(list)
baseq_minor=defaultdict(list)
leftpos_minor=defaultdict(list)
mismatches_minor=defaultdict(list)
minor_plus=dict()
minor_minus=dict()
minor_read1=dict()
minor_read2=dict()
baseq_minor_near1b=defaultdict(list)
seqpos_minor=defaultdict(list)

#genome="/home/yd65/tools/MosaicHunter/resources/human_g1k_v37_decoy.fasta.fai"
file0=open(genome)
chr_sizes=dict()
for line in file0:
	line=line.rstrip()
	fields=line.split('\t')
	chr_sizes[fields[0]]=chr_sizes.get(fields[0],fields[1])
file0.close()

def process_line(line):
#file=open(input_pos)
#1       1072410 1072411 C       A       Walsh
#for line in file:
	line=line.rstrip()
	fields=line.split('\t')
	sample=fields[5]
	chr=fields[0]
	if not re.search("MT",chr):
		pos=int(fields[2])
		pos1=max(0,int(pos)-1)
		pos2=min(int(chr_sizes[chr]),int(pos)+1)
		major_allele=fields[3]
		minor_allele=fields[4]
		name=str(sample)+'~'+str(chr)+'~'+str(pos)+"~"+str(major_allele)+"~"+str(minor_allele)
		input_bam=bam_dir+"/"+str(sample)+".bam"
		a=pysam.AlignmentFile(input_bam, "rb")
		chrom=str(chr)
		start=int(pos)-1
		end=int(pos)
		major_plus[name]=0
		minor_plus[name]=0
		major_minus[name]=0
		minor_minus[name]=0
		major_read1[name]=0
		minor_read1[name]=0
		major_read2[name]=0
		minor_read2[name]=0
		major_ids[name]=list()
		minor_ids[name]=list()
		conflict_num[name]=0
		context1_count[name]=context1_count.get(name,0)
		context2_count[name]=context2_count.get(name,0)
		mismatches_major[name]=list()
		mismatches_minor[name]=list()
		try:
			if len(major_allele)==len(minor_allele) and not(major_allele==minor_allele):
				if len(major_allele)==1:
					state="SNP"
				elif len(major_allele)>1:
					state="MNP"
				length=0
				major_softclippedreads=0
				minor_softclippedreads=0
				major_num=0
				minor_num=0
				context_20bp=str(reference[chrom][max(1,int(pos)-11):min(int(pos)+10,int(chr_sizes[chrom]))])
				GCcontent=(context_20bp.count('G')+context_20bp.count('C'))/len(context_20bp)
				context1[name]=reference[chrom][int(pos)-2:int(pos)+1]
				context2[name]=(base[str(reference[chrom][int(pos)-2:int(pos)-1])]+base[str(reference[chrom][int(pos)-1:int(pos)])]+base[str(reference[chrom][int(pos):int(pos)+1])])[::-1]
				for pileupcolumn in a.pileup(chrom, start, end):
					for pileupread in pileupcolumn.pileups:
						if pileupread.indel !=0:
							continue
						try:
							querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(major_allele)]
							if pileupcolumn.pos==pos-1 and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024):
								#if sequencing_type="PE":
								#print pileupcolumn.pos
								if querybase==major_allele:
									major_num=major_num+1
									if (pileupread.alignment.get_cigar_stats()[0][4])>=10:
										major_softclippedreads=major_softclippedreads+1
									major_ids[name].append(pileupread.alignment.query_name)
									querypos_major[name].append(pileupread.query_position)
									#print querypos_major[name]
									mapq_major[name].append(pileupread.alignment.mapping_quality)
									baseq_major[name].append(pileupread.alignment.query_qualities[pileupread.query_position])
									leftpos_major[name].append(pileupread.alignment.reference_start)
									#mismatches_major[name].append(filterPick(pileupread.alignment.tags,'NM'))
									mismatches_major[name].append(pileupread.alignment.get_tag('NM'))
									#print mismatches_major[name]
									#print mismatches_major[name]
									if not pileupread.alignment.is_reverse:
										major_plus[name]=major_plus.get(name,0)+1
									elif pileupread.alignment.is_reverse:
										major_minus[name]=major_minus.get(name,0)+1
				
									if pileupread.alignment.flag & 64:
										major_read1[name]=major_read1.get(name,0)+1
									elif pileupread.alignment.flag &128:
										major_read2[name]=major_read2.get(name,0)+1
								
									
									if pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
									##if pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
										seqpos_major[name].append(pileupread.query_position)
										#print seqpos_major[name]
										if pileupread.query_position < len(pileupread.alignment.query_sequence)-1:
											baseq_major_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+1])
										#elif pileupread.query_position == len(pileupread.alignment.query_sequence)-1:
										#	baseq_major_near1b[name].append("end")
									elif pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
									##elif pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
										seqpos_major[name].append(len(pileupread.alignment.query_sequence)-pileupread.query_position)
										if pileupread.query_position >=1 :
											baseq_major_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position-1])
											#print baseq_major_near1b[name]
						
							##elif pileupcolumn.pos==pos-1 and querybase==minor_allele and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024):
								elif querybase==minor_allele:
									minor_num=minor_num+1
									if int(pileupread.alignment.get_cigar_stats()[0][4])>10:
										minor_softclippedreads=minor_softclippedreads+1
									minor_ids[name].append(pileupread.alignment.query_name)
									querypos_minor[name].append(pileupread.query_position)
									mapq_minor[name].append(pileupread.alignment.mapping_quality)
									baseq_minor[name].append(pileupread.alignment.query_qualities[pileupread.query_position])
									leftpos_minor[name].append(pileupread.alignment.reference_start)
									#mismatches_minor[name].append(filterPick(pileupread.alignment.tags,'NM'))
									mismatches_minor[name].append(pileupread.alignment.get_tag('NM'))
									if not pileupread.alignment.is_reverse:
										minor_plus[name]=minor_plus.get(name,0)+1
									elif pileupread.alignment.is_reverse:
										minor_minus[name]=minor_minus.get(name,0)+1
									
									if pileupread.alignment.flag & 64:
										minor_read1[name]=minor_read1.get(name,0)+1
									elif pileupread.alignment.flag &128:
										minor_read2[name]=minor_read2.get(name,0)+1
									
									if pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
									#if pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
										seqpos_minor[name].append(pileupread.query_position)
										if pileupread.query_position < len(pileupread.alignment.query_sequence)-1:
											baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+1])
				#      elif pileupread.query_position == len(pileupread.alignment.query_sequence)-1:
				#       baseq_minor_near1b[name].append("end")
									elif pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
									#elif pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
										seqpos_minor[name].append(len(pileupread.alignment.query_sequence)-pileupread.query_position)
										if pileupread.query_position >=1 :
											baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position-1])
			
									if pileupread.alignment.is_proper_pair:
									##http://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/
										if pileupread.alignment.flag & 64 and (not pileupread.alignment.is_reverse):
											context1_count[name]=context1_count.get(name,int(0))+int(1)
										elif pileupread.alignment.flag & 128 and (pileupread.alignment.is_reverse):
											context2_count[name]=context2_count.get(name,int(0))+int(1)
										elif pileupread.alignment.flag & 64 and (pileupread.alignment.is_reverse):
											context2_count[name]=context2_count.get(name,int(0))+int(1)
										elif pileupread.alignment.flag & 128 and (not pileupread.alignment.is_reverse):
											context1_count[name]=context1_count.get(name,int(0))+int(1)
			#      elif pileupread.query_position <1:
			#       baseq_minor_near1b[name].append("end")
						except:
							#print ("error: ", chrom,pos)
							continue   
				conflict_reads=set(major_ids[name]) & set(minor_ids[name])
				conflict_num[name]=len(conflict_reads)
						
				for pileupcolumn in a.pileup(str(chrom), max(0,int(start)-2000), min(int(end)+2000,int(chr_sizes[str(chr)]))):
					if pileupcolumn.pos==pos-2000:
						dp_far[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos-1500:
						dp_far[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos-1000:
						dp_far[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos-500:
						dp_far[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos+500:
						dp_far[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos+1000:
						dp_far[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos+1500:
						dp_far[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos+2000:
						dp_far[name].append(pileupcolumn.n)
					
					elif pileupcolumn.pos==pos-200:
						dp_near[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos-100:
						dp_near[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos-50:
						dp_near[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos-1:
						dp_near[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos+50:
						dp_near[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos+100:
						dp_near[name].append(pileupcolumn.n)
					elif pileupcolumn.pos==pos+200:
						dp_near[name].append(pileupcolumn.n)
			
				u, p_value = mannwhitneyu(dp_far[name], dp_near[name])
				return name,','.join(str(x) for x in querypos_major[name])+",",  ','.join(str(x) for x in querypos_minor[name])+",",  ','.join(str(x) for x in leftpos_major[name])+",",  ','.join(str(x) for x in leftpos_minor[name])+",",  ','.join(str(x) for x in seqpos_major[name])+",",  ','.join(str(x) for x in seqpos_minor[name])+",",  ','.join(str(x) for x in mapq_major[name])+",",  ','.join(str(x) for x in mapq_minor[name])+",", ','.join(str(x) for x in baseq_major[name])+",",  ','.join(str(x) for x in baseq_minor[name])+",",  ','.join(str(x) for x in baseq_major_near1b[name])+",", ','.join(str(x) for x in baseq_minor_near1b[name])+",", major_plus[name],major_minus[name],minor_plus[name],minor_minus[name], str(context1[name]), str(context2[name]), context1_count[name],context2_count[name],','.join(str(x) for x in mismatches_major[name])+",",','.join(str(x) for x in mismatches_minor[name])+",", major_read1[name],major_read2[name],minor_read1[name],minor_read2[name],np.mean(dp_near[name]),np.mean(dp_far[name]),p_value,conflict_num[name], mappability[name], state, str(length),str(GCcontent), str(major_softclippedreads/major_num), str(minor_softclippedreads/minor_num)
	
			#elif len(major_allele)>1 and len(minor_allele)==1:
			elif len(major_allele) > len(minor_allele):
				major_num=0
				minor_num=0
				major_softclippedreads=0
				minor_softclippedreads=0
				state="DEL"
				length=len(major_allele)-len(minor_allele)
				#context1[name]=reference[chrom][int(pos)-2:int(pos)+1]
				context1_10bp=reference[chrom][max(1,int(pos)-11):min(int(pos)+1,int(chr_sizes[chrom]))]
				context2_10bp=reference[chrom][max(1,int(pos)-1):min(int(pos)+10,int(chr_sizes[chrom]))]
				context_20bp=str(reference[chrom][max(1,int(pos)-11):min(int(pos)+10,int(chr_sizes[chrom]))])
				GCcontent=(context_20bp.count('G')+context_20bp.count('C'))/len(context_20bp)
				context1[name]=reference[chrom][int(pos)-2:int(pos)+1]
				context2[name]=(base[str(reference[chrom][int(pos)-2:int(pos)-1])]+base[str(reference[chrom][int(pos)-1:int(pos)])]+base[str(reference[chrom][int(pos):int(pos)+1])])[::-1]
	
				if_homopolymer="No"
				for item in homopolymers:
					if re.search(str(item), str(context1_10bp)) or re.search(str(item),str(context2_10bp)):
						if_homopolymer="Yes"
						break
				if if_homopolymer=="No":			
					for pileupcolumn in a.pileup(chrom, start, end):
						for pileupread in pileupcolumn.pileups:
							try:
								if pileupread.indel==0:
									#querybase=pileupread.alignment.query_sequence[pileupread.query_position]
									querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(minor_allele)]
									#print(querybase,major_allele, minor_allele)
									#if int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(major_allele)[1]: #and pileuperead.alignment.mapping_quality>=10:
									#if int(pileupcolumn.pos)==int(pos)-1 and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024) and (not (pileupread.alignment.query_alignment_end<=int(pos) and (pileupread.alignment.cigar[-1][0]==4 or pileupread.alignment.cigar[-1][0]==5))) and (not (pileupread.alignment.query_alignment_start>=int(pos)-2 and (pileupread.alignment.cigar[0][0]==4 or pileupread.alignment.cigar[0][0]==5))): #and pileuperead.alignment.mapping_quality>=10:
									if int(pileupcolumn.pos)==int(pos)-1 and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024) and (not (pileupread.alignment.query_alignment_end<=int(pos) and (pileupread.alignment.cigar[-1][0]==4 or pileupread.alignment.cigar[-1][0]==5) and re.search(pileupread.alignment.query_sequence[pileupread.alignment.query_alignment_end:][:length],major_allele))) and (not (pileupread.alignment.query_alignment_start>=int(pos)-2 and (pileupread.alignment.cigar[0][0]==4 or pileupread.alignment.cigar[0][0]==5) and re.search(pileupread.alignment.query_sequence[:pileupread.alignment.query_alignment_start][-length:],major_allele) )): 
#and pileuperead.alignment.mapping_quality>=10:
										major_num=major_num+1
										if pileupread.alignment.get_cigar_stats()[0][4]>10:
											major_softclippedreads=major_softclippedreads+1
										major_ids[name].append(pileupread.alignment.query_name)
										querypos_major[name].append(pileupread.query_position)
										mapq_major[name].append(pileupread.alignment.mapping_quality)
										baseq_major[name].append(pileupread.alignment.query_qualities[pileupread.query_position])
										leftpos_major[name].append(pileupread.alignment.reference_start)
										mismatches_major[name].append(pileupread.alignment.get_tag('NM'))
	
										if not pileupread.alignment.is_reverse:
											major_plus[name]=major_plus.get(name,0)+1
										elif pileupread.alignment.is_reverse:
											major_minus[name]=major_minus.get(name,0)+1
					
										if pileupread.alignment.flag & 64:
											major_read1[name]=major_read1.get(name,0)+1
										elif pileupread.alignment.flag &128:
											major_read2[name]=major_read2.get(name,0)+1
										
										if pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
											seqpos_major[name].append(pileupread.query_position)
											#print seqpos_major[name]
											#if pileupread.query_position < len(pileupread.alignment.query_sequence)-1:
											if pileupread.query_position < len(pileupread.alignment.query_sequence)-len(major_allele)-1:
												qualities=pileupread.alignment.query_qualities[pileupread.query_position+len(major_allele)+1:len(pileupread.alignment.query_sequence)]
												baseq_average_leftbases=sum(qualities)/len(qualities)
												baseq_major_near1b[name].append(baseq_average_leftbases)
												#baseq_major_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+len(major_allele)+1])
										elif pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
											seqpos_major[name].append(len(pileupread.alignment.query_sequence)-pileupread.query_position)
											if pileupread.query_position >=len(major_allele) :
												qualities=pileupread.alignment.query_qualities[0:pileupread.query_position-len(major_allele)]
												baseq_average_leftbases=sum(qualities)/len(qualities)
												baseq_major_near1b[name].append(baseq_average_leftbases)
							
#									elif int(pileupcolumn.pos)==int(pos)-1 and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024) and ((pileupread.alignment.query_alignment_end<=int(pos) and (pileupread.alignment.cigar[-1][0]==4 or pileupread.alignment.cigar[-1][0]==5)) or (pileupread.alignment.query_alignment_start>=int(pos)-2 and (pileupread.alignment.cigar[0][0]==4 or pileupread.alignment.cigar[0][0]==5))): #and pileuperead.alignment.mapping_quality>=10:

								##elif pileupcolumn.pos==pos-1 and querybase==minor_allele and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024):
								elif (pileupread.indel<0 and int(pileupcolumn.pos)==int(pos)-1 and pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(minor_allele)]==minor_allele and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024)) or (pileupread.indel==0 and (int(pileupcolumn.pos)==int(pos)-1 and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024) and ((pileupread.alignment.query_alignment_end<=int(pos) and (pileupread.alignment.cigar[-1][0]==4 or pileupread.alignment.cigar[-1][0]==5) and re.search(pileupread.alignment.query_sequence[pileupread.alignment.query_alignment_end:][:length],major_allele)) or (pileupread.alignment.query_alignment_start>=int(pos)-2 and (pileupread.alignment.cigar[0][0]==4 or pileupread.alignment.cigar[0][0]==5) and re.search(pileupread.alignment.query_sequence[:pileupread.alignment.query_alignment_start][-length:],major_allele))))): 


#del <0 or clipped reads
								#querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+1]
#									querybase=pileupread.alignment.query_sequence[pileupread.query_position]
								#print(querybase,major_allele, minor_allele)
#									if (int(pileupcolumn.pos)==int(pos)-1) and (str(querybase)==str(minor_allele)) and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024):
									minor_num=minor_num+1
									if pileupread.alignment.get_cigar_stats()[0][4]>10:
										minor_softclippedreads=minor_softclippedreads+1
									minor_ids[name].append(pileupread.alignment.query_name)
									querypos_minor[name].append(pileupread.query_position)
									mapq_minor[name].append(pileupread.alignment.mapping_quality)
									baseq_minor[name].append(pileupread.alignment.query_qualities[pileupread.query_position])
									leftpos_minor[name].append(pileupread.alignment.reference_start)
									#mismatches_minor[name].append(filterPick(pileupread.alignment.tags,'NM'))
									mismatches_minor[name].append(pileupread.alignment.get_tag('NM'))
									if not pileupread.alignment.is_reverse:
										minor_plus[name]=minor_plus.get(name,0)+1
									elif pileupread.alignment.is_reverse:
										minor_minus[name]=minor_minus.get(name,0)+1
									
									if pileupread.alignment.flag & 64:
										minor_read1[name]=minor_read1.get(name,0)+1
									elif pileupread.alignment.flag &128:
										minor_read2[name]=minor_read2.get(name,0)+1
									
									if pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
									#if pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
										seqpos_minor[name].append(pileupread.query_position)
										#if pileupread.query_position < len(pileupread.alignment.query_sequence)-1:
										#	baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+1])
										if pileupread.query_position < len(pileupread.alignment.query_sequence)-len(major_allele)-1:
#											baseq_average=sum(pileupread.alignment.query_qualities[pileupread.query_position+len(major_allele):])
											
#											baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+len(major_allele)+1])


											qualities=pileupread.alignment.query_qualities[pileupread.query_position+len(major_allele)+1:len(pileupread.alignment.query_sequence)]
											baseq_average_leftbases=sum(qualities)/len(qualities)
											baseq_minor_near1b[name].append(baseq_average_leftbases)
				#      elif pileupread.query_position == len(pileupread.alignment.query_sequence)-1:
				#       baseq_minor_near1b[name].append("end")
									elif pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
									#elif pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
										seqpos_minor[name].append(len(pileupread.alignment.query_sequence)-pileupread.query_position)
										#if pileupread.query_position >=1 :
										#	baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position-1])
										if pileupread.query_position >=len(major_allele) :
											#baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position-len(major_allele)])
											qualities=pileupread.alignment.query_qualities[0:pileupread.query_position-len(major_allele)]
											baseq_average_leftbases=sum(qualities)/len(qualities)
											baseq_minor_near1b[name].append(baseq_average_leftbases)
			
									if pileupread.alignment.is_proper_pair:
									##http://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/
										if pileupread.alignment.flag & 64 and (not pileupread.alignment.is_reverse):
											context1_count[name]=context1_count.get(name,int(0))+int(1)
										elif pileupread.alignment.flag & 128 and (pileupread.alignment.is_reverse):
											context2_count[name]=context2_count.get(name,int(0))+int(1)
										elif pileupread.alignment.flag & 64 and (pileupread.alignment.is_reverse):
											context2_count[name]=context2_count.get(name,int(0))+int(1)
										elif pileupread.alignment.flag & 128 and (not pileupread.alignment.is_reverse):
											context1_count[name]=context1_count.get(name,int(0))+int(1)
				#      elif pileupread.query_position <1:
				#       baseq_minor_near1b[name].append("end")
							except:
								#print ("error: ", chrom,pos)
								continue   
					conflict_reads=set(major_ids[name]) & set(minor_ids[name])
					conflict_num[name]=len(conflict_reads)
							
					for pileupcolumn in a.pileup(str(chrom), max(0,int(start)-2000), min(int(end)+2000,int(chr_sizes[str(chr)]))):
						if pileupcolumn.pos==pos-2000:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-1500:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-1000:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-500:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+500:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+1000:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+1500:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+2000:
							dp_far[name].append(pileupcolumn.n)
						
						elif pileupcolumn.pos==pos-200:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-100:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-50:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-1:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+50:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+100:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+200:
							dp_near[name].append(pileupcolumn.n)
				
					u, p_value = mannwhitneyu(dp_far[name], dp_near[name])
					return name,','.join(str(x) for x in querypos_major[name])+",",  ','.join(str(x) for x in querypos_minor[name])+",",  ','.join(str(x) for x in leftpos_major[name])+",",  ','.join(str(x) for x in leftpos_minor[name])+",",  ','.join(str(x) for x in seqpos_major[name])+",",  ','.join(str(x) for x in seqpos_minor[name])+",",  ','.join(str(x) for x in mapq_major[name])+",",  ','.join(str(x) for x in mapq_minor[name])+",", ','.join(str(x) for x in baseq_major[name])+",",  ','.join(str(x) for x in baseq_minor[name])+",",  ','.join(str(x) for x in baseq_major_near1b[name])+",", ','.join(str(x) for x in baseq_minor_near1b[name])+",", major_plus[name],major_minus[name],minor_plus[name],minor_minus[name], str(context1[name]), str(context2[name]), context1_count[name],context2_count[name],','.join(str(x) for x in mismatches_major[name])+",",','.join(str(x) for x in mismatches_minor[name])+",", major_read1[name],major_read2[name],minor_read1[name],minor_read2[name],np.mean(dp_near[name]),np.mean(dp_far[name]),p_value,conflict_num[name], mappability[name], state,str(length),str(GCcontent), str(major_softclippedreads/major_num), str(minor_softclippedreads/minor_num)
			#print (fo)
	
	
	##end: paste from phase
	
			#elif len(major_allele)==1 and len(minor_allele)>1:
			elif len(major_allele) < len(minor_allele):
				state="INS"
				length=len(minor_allele)-len(major_allele)
				major_num=0
				minor_num=0
				major_softclippedreads=0
				minor_softclippedreads=0
				context1[name]=reference[chrom][int(pos)-2:int(pos)+1]
				context1_10bp=reference[chrom][max(1,int(pos)-11):min(int(pos)+1,int(chr_sizes[chrom]))]
				context2_10bp=reference[chrom][max(1,int(pos)-1):min(int(pos)+10,int(chr_sizes[chrom]))]
				context_20bp=str(reference[chrom][max(1,int(pos)-11):min(int(pos)+10,int(chr_sizes[chrom]))])
				GCcontent=(context_20bp.count('G')+context_20bp.count('C'))/len(context_20bp)
				context1[name]=reference[chrom][int(pos)-2:int(pos)+1]
				context2[name]=(base[str(reference[chrom][int(pos)-2:int(pos)-1])]+base[str(reference[chrom][int(pos)-1:int(pos)])]+base[str(reference[chrom][int(pos):int(pos)+1])])[::-1]
	
				if_homopolymer="No"
				for item in homopolymers:
					if re.search(str(item), str(context1_10bp)) or re.search(str(item),str(context2_10bp) or re.search(str(item),minor_allele)):
						if_homopolymer="Yes"
						break
				if if_homopolymer=="No":			
					for pileupcolumn in a.pileup(chrom, start, end):
						for pileupread in pileupcolumn.pileups:
							try:
								if pileupread.indel==0:
									querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(major_allele)]
									#print(querybase,major_allele, minor_allele)
									#if int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(major_allele)[1]: #and pileuperead.alignment.mapping_quality>=10:
#									if int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(major_allele) and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024): #and pileuperead.alignment.mapping_quality>=10:

#									if int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(major_allele) and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024) and (not (pileupread.alignment.query_alignment_end<=int(pos) and (pileupread.alignment.cigar[-1][0]==4 or pileupread.alignment.cigar[-1][0]==5))) and (not (pileupread.alignment.query_alignment_start>=int(pos)-2 and (pileupread.alignment.cigar[0][0]==4 or pileupread.alignment.cigar[0][0]==5))): #and pileuperead.alignment.mapping_quality>=10:

									if int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(major_allele) and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024) and (not (pileupread.alignment.query_alignment_end<=int(pos) and (pileupread.alignment.cigar[-1][0]==4 or pileupread.alignment.cigar[-1][0]==5) and re.search(pileupread.alignment.query_sequence[pileupread.alignment.query_alignment_end:][:length],major_allele))) and (not (pileupread.alignment.query_alignment_start>=int(pos)-2 and (pileupread.alignment.cigar[0][0]==4 or pileupread.alignment.cigar[0][0]==5) and re.search(pileupread.alignment.query_sequence[:pileupread.alignment.query_alignment_start][-length:],major_allele) )): 


										major_num=major_num+1
										if pileupread.alignment.get_cigar_stats()[0][4]>10:
											major_softclippedreads=major_softclippedreads+1
										major_ids[name].append(pileupread.alignment.query_name)
										querypos_major[name].append(pileupread.query_position)
										mapq_major[name].append(pileupread.alignment.mapping_quality)
										baseq_major[name].append(pileupread.alignment.query_qualities[pileupread.query_position])
										leftpos_major[name].append(pileupread.alignment.reference_start)
										mismatches_major[name].append(pileupread.alignment.get_tag('NM'))
	
										if not pileupread.alignment.is_reverse:
											major_plus[name]=major_plus.get(name,0)+1
										elif pileupread.alignment.is_reverse:
											major_minus[name]=major_minus.get(name,0)+1
					
										if pileupread.alignment.flag & 64:
											major_read1[name]=major_read1.get(name,0)+1
										elif pileupread.alignment.flag &128:
											major_read2[name]=major_read2.get(name,0)+1
									
										
										if pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
										##if pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
											seqpos_major[name].append(pileupread.query_position)
											#print seqpos_major[name]
											#if pileupread.query_position < len(pileupread.alignment.query_sequence)-1:
											if pileupread.query_position < len(pileupread.alignment.query_sequence)-len(minor_allele)-1:
												#baseq_major_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+1])
	#											baseq_major_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+len(minor_allele)+1])
												qualities=pileupread.alignment.query_qualities[pileupread.query_position+len(minor_allele)+1:len(pileupread.alignment.query_sequence)]
												baseq_average_leftbases=sum(qualities)/len(qualities)
												baseq_major_near1b[name].append(baseq_average_leftbases)
											#elif pileupread.query_position == len(pileupread.alignment.query_sequence)-1:
											#	baseq_major_near1b[name].append("end")
	
		
	
										elif pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
										##elif pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
											seqpos_major[name].append(len(pileupread.alignment.query_sequence)-pileupread.query_position)
											#if pileupread.query_position >=1 :
											if pileupread.query_position >=len(minor_allele) :
	#											baseq_major_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position-len(minor_allele)])
												qualities=pileupread.alignment.query_qualities[0:len(pileupread.alignment.query_sequence)-len(minor_allele)]
												baseq_average_leftbases=sum(qualities)/len(qualities)
												baseq_major_near1b[name].append(baseq_average_leftbases)
												#print baseq_major_near1b[name]
							
#								if pileupread.indel > 0: ###ins>0
#								elif pileupread.indel>0 or (pileupread.indel==0 and (int(pileupcolumn.pos)==int(pos)-1 and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024) and ((pileupread.alignment.query_alignment_end<=int(pos) and (pileupread.alignment.cigar[-1][0]==4 or pileupread.alignment.cigar[-1][0]==5)) or (pileupread.alignment.query_alignment_start>=int(pos)-2 and (pileupread.alignment.cigar[0][0]==4 or pileupread.alignment.cigar[0][0]==5))))): #del <0 or clipped reads

								elif (pileupread.indel>0 and int(pileupcolumn.pos)==int(pos)-1 and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024)) or (pileupread.indel==0 and (int(pileupcolumn.pos)==int(pos)-1 and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024) and ((pileupread.alignment.query_alignment_end<=int(pos) and (pileupread.alignment.cigar[-1][0]==4 or pileupread.alignment.cigar[-1][0]==5) and re.search(pileupread.alignment.query_sequence[pileupread.alignment.query_alignment_end:][:length],minor_allele)) or (pileupread.alignment.query_alignment_start>=int(pos)-2 and (pileupread.alignment.cigar[0][0]==4 or pileupread.alignment.cigar[0][0]==5) and re.search(pileupread.alignment.query_sequence[:pileupread.alignment.query_alignment_start][-length:],minor_allele))))): 

								#querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+(len(minor_allele)-len(major_allele))+1]
								#if int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(minor_allele):
#									if (int(pileupcolumn.pos)==int(pos)-1) and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024):
									minor_num=minor_num+1
									if pileupread.alignment.get_cigar_stats()[0][4]>10:
										minor_softclippedreads=minor_softclippedreads+1
									#print(querybase,major_allele, minor_allele)
									minor_ids[name].append(pileupread.alignment.query_name)
									querypos_minor[name].append(pileupread.query_position)
									mapq_minor[name].append(pileupread.alignment.mapping_quality)
									queryscore=pileupread.alignment.query_qualities[pileupread.query_position:pileupread.query_position+(len(minor_allele)-len(major_allele))+1]
									baseq_mean =sum(queryscore)/len(queryscore)
									#baseq_minor[name].append(pileupread.alignment.query_qualities[pileupread.query_position])
									baseq_minor[name].append(baseq_mean)
									leftpos_minor[name].append(pileupread.alignment.reference_start)
									#mismatches_minor[name].append(filterPick(pileupread.alignment.tags,'NM'))
									mismatches_minor[name].append(pileupread.alignment.get_tag('NM'))
									if not pileupread.alignment.is_reverse:
										minor_plus[name]=minor_plus.get(name,0)+1
									elif pileupread.alignment.is_reverse:
										minor_minus[name]=minor_minus.get(name,0)+1
									
									if pileupread.alignment.flag & 64:
										minor_read1[name]=minor_read1.get(name,0)+1
									elif pileupread.alignment.flag &128:
										minor_read2[name]=minor_read2.get(name,0)+1
									
									if pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
									#if pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
										seqpos_minor[name].append(pileupread.query_position)
										#if pileupread.query_position < len(pileupread.alignment.query_sequence)-1:
										#	baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+1])
										if pileupread.query_position < len(pileupread.alignment.query_sequence)-len(minor_allele)-1:
											#baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+len(minor_allele)+1])
											qualities=pileupread.alignment.query_qualities[pileupread.query_position+len(minor_allele)+1:len(pileupread.alignment.query_sequence)]
											baseq_average_leftbases=sum(qualities)/len(qualities)
											baseq_minor_near1b[name].append(baseq_average_leftbases)
				#      elif pileupread.query_position == len(pileupread.alignment.query_sequence)-1:
				#       baseq_minor_near1b[name].append("end")
									elif pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
									#elif pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
										seqpos_minor[name].append(len(pileupread.alignment.query_sequence)-pileupread.query_position)
										#if pileupread.query_position >=1 :
										#	baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position-1])
										if pileupread.query_position >=len(minor_allele) :
											#baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position-len(minor_allele)])
											qualities=pileupread.alignment.query_qualities[0:len(pileupread.alignment.query_sequence)-len(minor_allele)]
											baseq_average_leftbases=sum(qualities)/len(qualities)
											baseq_minor_near1b[name].append(baseq_average_leftbases)
			
									if pileupread.alignment.is_proper_pair:
									##http://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/
										if pileupread.alignment.flag & 64 and (not pileupread.alignment.is_reverse):
											context1_count[name]=context1_count.get(name,int(0))+int(1)
										elif pileupread.alignment.flag & 128 and (pileupread.alignment.is_reverse):
											context2_count[name]=context2_count.get(name,int(0))+int(1)
										elif pileupread.alignment.flag & 64 and (pileupread.alignment.is_reverse):
											context2_count[name]=context2_count.get(name,int(0))+int(1)
										elif pileupread.alignment.flag & 128 and (not pileupread.alignment.is_reverse):
											context1_count[name]=context1_count.get(name,int(0))+int(1)
				#      elif pileupread.query_position <1:
				#       baseq_minor_near1b[name].append("end")
							except:
								#print ("error: ", chrom,pos)
								continue   
					conflict_reads=set(major_ids[name]) & set(minor_ids[name])
					conflict_num[name]=len(conflict_reads)
							
					for pileupcolumn in a.pileup(str(chrom), max(0,int(start)-2000), min(int(end)+2000,int(chr_sizes[str(chr)]))):
						if pileupcolumn.pos==pos-2000:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-1500:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-1000:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-500:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+500:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+1000:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+1500:
							dp_far[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+2000:
							dp_far[name].append(pileupcolumn.n)
						
						elif pileupcolumn.pos==pos-200:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-100:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-50:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos-1:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+50:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+100:
							dp_near[name].append(pileupcolumn.n)
						elif pileupcolumn.pos==pos+200:
							dp_near[name].append(pileupcolumn.n)
				
					u, p_value = mannwhitneyu(dp_far[name], dp_near[name])
					return name,','.join(str(x) for x in querypos_major[name])+",",  ','.join(str(x) for x in querypos_minor[name])+",",  ','.join(str(x) for x in leftpos_major[name])+",",  ','.join(str(x) for x in leftpos_minor[name])+",",  ','.join(str(x) for x in seqpos_major[name])+",",  ','.join(str(x) for x in seqpos_minor[name])+",",  ','.join(str(x) for x in mapq_major[name])+",",  ','.join(str(x) for x in mapq_minor[name])+",", ','.join(str(x) for x in baseq_major[name])+",",  ','.join(str(x) for x in baseq_minor[name])+",",  ','.join(str(x) for x in baseq_major_near1b[name])+",", ','.join(str(x) for x in baseq_minor_near1b[name])+",", major_plus[name],major_minus[name],minor_plus[name],minor_minus[name], str(context1[name]), str(context2[name]), context1_count[name],context2_count[name],','.join(str(x) for x in mismatches_major[name])+",",','.join(str(x) for x in mismatches_minor[name])+",", major_read1[name],major_read2[name],minor_read1[name],minor_read2[name],np.mean(dp_near[name]),np.mean(dp_far[name]),p_value,conflict_num[name], mappability[name], state,str(length),str(GCcontent),  str(major_softclippedreads/major_num), str(minor_softclippedreads/minor_num)
			#print (fo)
		except:
			print ("not enough alt reads: ",chrom, start,end)
	
	
	
#pool=Pool(processes=int(n_jobs))
#pool.map(run_cmd,cmd_list1,1)
	
if __name__ == "__main__":
	pool = Pool(processes=int(n_jobs))
	with open(input_pos) as source_file:
	# chunk the work into batches of 4 lines at a time
		#pool.map(process_line, source_file,1)
		result=pool.map(process_line, source_file,1)
		for atuple in result:
			try:
				print (' '.join(str(x) for x in atuple),file=fo)
			except:
				continue





