#!/usr/bin/env python
import sys
program=sys.argv[0]
arguments=sys.argv[1:]
count=len(arguments)
if count !=6:
	print ("Usage: python phase_WGS.py bam_dir data_dir ref_fasta n_jobs input_positions(file format:sample chr pos ref alt, sep=\\t) min_dp_inforSNPs(int)")
	sys.exit(1)
elif count==6:
	#sample=sys.argv[1]
	bam_dir_tmp=sys.argv[1]
	data_dir_tmp=sys.argv[2]
	ref_fasta=sys.argv[3]
	n_jobs=sys.argv[4]
	input_pos=sys.argv[5]
	min_dp_inforSNPs=int(sys.argv[6])
	log_file='multiple_inforSNPs.log'

import regex as re
import pysam
from pyfaidx import Fasta
import scipy.stats
from multiprocessing import Pool
from subprocess import *
import subprocess
import pysamstats
import os
from pandas import *
import itertools
from collections import OrderedDict, defaultdict

def run_cmd(cmd):
	Popen(cmd, shell=True, stdout=PIPE).communicate()

if bam_dir_tmp.endswith('/'):
	bam_dir=bam_dir_tmp[:-1]
else:
	bam_dir=bam_dir_tmp

if data_dir_tmp.endswith('/'):
	data_dir=data_dir_tmp[:-1]
else:
	data_dir=data_dir_tmp

output_dir=data_dir+'/phasing'
os.system("mkdir -p "+output_dir)


reference = Fasta(ref_fasta)
genome=ref_fasta+".fai"
file0=open(genome)
chr_sizes=dict()
for line in file0:
	line=line.rstrip()
	fields=line.split('\t')
	chr_sizes[fields[0]]=chr_sizes.get(fields[0],fields[1])
file0.close()

##input_bam=bam_dir+"/"+sample+".bam"
##a=pysam.AlignmentFile(input_bam, "rb")
	
##head Vaccarino_brain.mutect2/final.passed.tsv
##1       2533604 G       A       Vaccarino_brain
cmd_list1=list()
cmd_list2=list()

##1st step: generate seperate bams supporting ref allale and alt alleles:
file=open(input_pos)
fo1=open(output_dir+"/all_candidates","w")
for line in file:
	line=line.rstrip()
	fields=line.split('\t')
	sample=fields[0]
	chr=fields[1]
	chrom=str(chr)
	pos=int(fields[2])
	major_allele=fields[3]
	minor_allele=fields[4]
	input_bam=bam_dir+"/"+sample+".bam"
	a=pysam.AlignmentFile(input_bam, "rb")
	f1=pysam.AlignmentFile(output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.bam","wb",template=a)
	f2=pysam.AlignmentFile(output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.bam","wb",template=a)
	f3=pysam.AlignmentFile(output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.bam","wb",template=a)
	name=sample+'_'+chr+'_'+str(pos)
	major_ids=list()
	minor_ids=list()
	start=int(pos)-1
	end=int(pos)
	major_num=int(0)
	minor_num=int(0)

	for pileupcolumn in a.pileup(chrom, start, end):
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

	start=int(pos)-1000
	end=int(pos)+1000
	conflictnum=0
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

	f1_name=output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.bam"
	f2_name=output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.bam"
	f3_name=output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.bam"
	
	f1_sorted_name=output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.sorted.bam"
	f2_sorted_name=output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.sorted.bam"
	f3_sorted_name=output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.sorted.bam"
	cmd_list1.append("samtools sort "+f1_name+" -o "+f1_sorted_name)
	cmd_list1.append("samtools sort "+f2_name+" -o "+f2_sorted_name)
	cmd_list1.append("samtools sort "+f3_name+" -o "+f3_sorted_name)
	cmd_list2.append("samtools index "+f1_sorted_name)
	cmd_list2.append("samtools index "+f2_sorted_name)
	cmd_list2.append("samtools index "+f3_sorted_name)
	print (sample,chr,pos,major_allele,minor_allele,major_num,minor_num,conflictnum,file=fo1)
file.close()
fo1.close()
	
pool=Pool(processes=int(n_jobs))
pool.map(run_cmd,cmd_list1,1)
pool.close()
pool.join()

pool=Pool(processes=int(n_jobs))
pool.map(run_cmd,cmd_list2,1)
pool.close()
pool.join()

##2nd step: extract candidate nearby inforSNPs (germline het):
file=open(output_dir+"/all_candidates")
merged_inforSNPs=open(output_dir+"/all.merged.inforSNPs.pos","w")
conflict_mosaic=dict()
#==> /n/data1/hms/dbmi/park/yanmei/simulated_bams_na12878/12878-50x.merged.inforSNPs.pos <==
#1 18089 G T 1 17385 G A 0
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
	mosaic_name=sample+';'+chr+';'+str(pos)+";"+major_allele+";"+minor_allele
	conflict_mosaic[mosaic_name]=conflictnum
	start=int(pos)-1
	end=int(pos)
	f3_sorted_name=output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.sorted.bam"

	if major_num>=3 and minor_num>=3:
		##print (chr,pos, major_num, minor_num)
		##print (f3_sorted_name)
		f3_alignment_file =pysam.AlignmentFile(f3_sorted_name,'rb')
		for rec in pysamstats.stat_variation(f3_alignment_file, fafile=ref_fasta, min_mapq=0, min_baseq=0):
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
					print (sample,chr,pos,major_allele,minor_allele,chr,int(rec['pos'])+1,max_base,alt_base,conflictnum, file=merged_inforSNPs)
		
file.close()
merged_inforSNPs.close()



##3rd step: generate the 2X2 table:
n_major_het1=dict()
n_major_het2=dict()
n_minor_het1=dict()
n_minor_het2=dict()

input_pos=output_dir+"/all.merged.inforSNPs.pos"
#1 1661246 T C 1 1660812 A G 0
#1 2585130 C G 1 2584706 A G 4
file=open(input_pos)
for line in file:
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
	if pos != inforSNP_pos:
		f1_sorted_name=output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.sorted.bam"
		f2_sorted_name=output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.sorted.bam"
		a1=pysam.AlignmentFile(f1_sorted_name,"rb")
		a2=pysam.AlignmentFile(f2_sorted_name,"rb")
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
file.close()

fo2=open(output_dir+"/all_2x2table","w")
print ("sample","chr","pos","ref","alt","pos_inforSNP","het1","het2","conflicted_reads","ref_het1_count","ref_het2_count","alt_het1_count","alt_het2_count",file=fo2)
#12878-200x 10_10002280_A_G_10002816_G_A_0 14 0 0 8

for k,v in sorted(n_major_het1.items()):
	#print >>fo, sample,k,n_major_het1[k],n_major_het2[k],n_minor_het1[k],n_minor_het2[k]
	print (' '.join(str(x) for x in k.split(";")), n_major_het1[k],n_major_het2[k],n_minor_het1[k],n_minor_het2[k],file=fo2)
fo2.close()


##last step: assign phasing state to each site
phase=defaultdict(dict)
inforSNPs=defaultdict(list)
input_file="all_2x2table"
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
		name=sample+";"+chr+";"+pos+";"+ref+";"+alt
		inforSNPs[name].append(sample+";"+chr+";"+pos_inforSNP+";"+het1_inforSNP+";"+het2_inforSNP)
	
		phase[name]['hap=2']=phase[name].get("hap=2",0)+0
		phase[name]['hap=3']=phase[name].get("hap=3",0)+0
		phase[name]['hap>3']=phase[name].get("hap>3",0)+0
		if C1+C2+C3+C4>=10:
			if ((((C2==float(0) or C1>C2*float(10) ) and ( C3==float(0) or C4>C3*float(10) )) or ((C1==float(0) or C2>C1*float(10)) and (C4==float(0) or C3>C4*float(10)))) and (C1+C2>float(3) and C3+C4>float(3))):
				phase[name]['hap=2']=phase[name].get("hap=2",0)+1
			elif (C1*C2==float(0) and C1+C2>=float(3) and C3*C4>float(0) and (C3>C4/float(10) or C3<C4*float(10))) or (C3*C4==float(0) and C3+C4>=float(3) and C1*C2>float(0) and (C1>C2/float(10) or C1<C2*float(10))):
				phase[name]['hap=3']=phase[name].get("hap=3",0)+1
			elif C1*C2>float(0) and C3*C4>float(0):
				phase[name]['hap>3']=phase[name].get("hap>3",0)+1

phasing_2by2 = dict()
fo=open(output_dir+"/all.phasing_2by2","w")
for k,v in sorted(phase.items()):
	phasing_2by2[k]=max(v,key=v.get)
	print (' '.join(str(x) for x in k.split(";")), v,max(v, key=v.get), file=fo)
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
		samfile=pysam.AlignmentFile(bam_dir+"/"+sample+".bam", "rb")
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
print ("sample","chr","pos","ref","alt","phasing","conflicting_reads",file=fo4)
for k,v in sorted(phasing_2by2.items()):
	if k in MT2_phasing_num:
		if MT2_phasing_num[k]=="wrong":
			phasing_2by2[k]="hap>3"
	print (k,v,conflict_mosaic[k], file=fo4)
#	print (' '.join(str(x) for x in k.split(";")), v, conflict_mosaic[k], file=fo4)
fo4.close()


