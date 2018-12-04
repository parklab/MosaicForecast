import sys
from pyfaidx import Fasta
from collections import defaultdict
import pysam
from scipy.stats import mannwhitneyu
import numpy as np
import time
from . import mosaicpc_util as mutil
from . import extract_feature_R

def run(opt):
    print ("extract_feature>run")
    if not option_check(opt):
        sys.exit(1)
    step1_extract_features(opt['ref_fasta'],opt['input_pos'], opt['bam_dir'], opt['output_feature'])
    step2_test_features(opt['output_feature'])

def step2_test_features(output_feature):
    print (output_feature)
    outRscript = output_feature + "_extract_feature.R"
    mutil.fileSave(outRscript,extract_feature_R.RSCRIPT ,'w')
    time.sleep(3)
    print ("Rscript "+outRscript+" " + output_feature)
    mutil.run_cmd("Rscript mosaicpc/extract_feature.R " + output_feature)

def step1_extract_features(ref_fasta, input_pos, bam_dir, output_feature):
    base={'A':'T','T':'A','G':'C','C':'G'}
    reference = Fasta(ref_fasta)
    chr_sizes=mutil.get_chr_sizes(ref_fasta+".fai")

    querypos_major=defaultdict(list)
    mapq_major=defaultdict(list)
    baseq_major=defaultdict(list)
    leftpos_major=defaultdict(list)
    mismatches_major=defaultdict(list)
    major_plus=dict()
    major_minus=dict()
    major_read1=dict()
    major_read2=dict()
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
    
    file=open(input_pos)
    #1       1072410 1072411 C       A       Walsh
    for line in file:
        line=line.rstrip()
        fields=line.split('\t')
        sample=fields[5]
        chr=fields[0]
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
        context1_count[name]=context1_count.get(name,0)
        context2_count[name]=context2_count.get(name,0)
        mismatches_major[name]=list()
        mismatches_minor[name]=list()
        context1[name]=reference[chrom][int(pos)-2:int(pos)+1]
        context2[name]=(base[str(reference[chrom][int(pos)-2:int(pos)-1])]+base[str(reference[chrom][int(pos)-1:int(pos)])]+base[str(reference[chrom][int(pos):int(pos)+1])])[::-1]

        for pileupcolumn in a.pileup(chrom, start, end):
            for pileupread in pileupcolumn.pileups:
                if pileupread.indel !=0:
                    continue
                try:
                    querybase=pileupread.alignment.query_sequence[pileupread.query_position]
                    if pileupcolumn.pos==pos-1 and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024):
                        #if sequencing_type="PE":
                        #print pileupcolumn.pos
                        if querybase==major_allele:
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
                                minor_read1[name]=major_read1.get(name,0)+1
                            elif pileupread.alignment.flag &128:
                                minor_read2[name]=major_read2.get(name,0)+1
                            
                            #if pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
                            if pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
                                seqpos_minor[name].append(pileupread.query_position)
                                if pileupread.query_position < len(pileupread.alignment.query_sequence)-1:
                                    baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+1])
        #      elif pileupread.query_position == len(pileupread.alignment.query_sequence)-1:
        #       baseq_minor_near1b[name].append("end")
                            #elif pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
                            elif pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
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
                    continue   
                
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
                
    fo=open(output_feature,"w")

    #header='id querypos_major querypos_minor leftpos_major leftpos_minor seqpos_major seqpos_minor mapq_major mapq_minor baseq_major baseq_minor baseq_major_near1b baseq_minor+near1b major_plus major_minus minor_plus minor_minus context_reference_forward context_reference_reverse context_antireference_forward context_antireference_reverse context_reference_forward_count context_reference_reverse_count context_antireference_forward_count context_antireference_reverse_count mismatches_major mismatches_minor major_read1 major_read2 minor_read1 minor_read2 dp_near dp_far dp_p'.split()
    #print (' '.join(header),file=fo)
    ##print >>fo,'\t'.join(header)
    #for k,v in sorted(querypos_major.items()):
    #	u, p_value = mannwhitneyu(dp_far[k], dp_near[k])
    #	print (k,','.join(str(x) for x in v)+",",  ','.join(str(x) for x in querypos_minor[k])+",",  ','.join(str(x) for x in leftpos_major[k])+",",  ','.join(str(x) for x in leftpos_minor[k])+",",  ','.join(str(x) for x in seqpos_major[k])+",",  ','.join(str(x) for x in seqpos_minor[k])+",",  ','.join(str(x) for x in mapq_major[k])+",",  ','.join(str(x) for x in mapq_minor[k])+",", ','.join(str(x) for x in baseq_major[k])+",",  ','.join(str(x) for x in baseq_minor[k])+",",  ','.join(str(x) for x in baseq_major_near1b[k])+",", ','.join(str(x) for x in baseq_minor_near1b[k])+",", major_plus[k],major_minus[k],minor_plus[k],minor_minus[k],context_reference_forward[k],context_reference_reverse[k],context_antireference_forward[k],context_antireference_reverse[k],context_reference_forward_count[k],context_reference_reverse_count[k],context_antireference_forward_count[k],context_antireference_reverse_count[k],','.join(str(x) for x in mismatches_major[k])+",",','.join(str(x) for x in mismatches_minor[k])+",", major_read1[k],major_read2[k],minor_read1[k],minor_read2[k],np.mean(dp_near[k]),np.mean(dp_far[k]),p_value,file=fo)
    #fo.close()

    header='id querypos_major querypos_minor leftpos_major leftpos_minor seqpos_major seqpos_minor mapq_major mapq_minor baseq_major baseq_minor baseq_major_near1b baseq_minor_near1b major_plus major_minus minor_plus minor_minus context1 context2 context1_count context2_count mismatches_major mismatches_minor major_read1 major_read2 minor_read1 minor_read2 dp_near dp_far dp_p'.split()
    print (' '.join(header),file=fo)
    #print >>fo,'\t'.join(header)
    for k,v in sorted(querypos_major.items()):
        u, p_value = mannwhitneyu(dp_far[k], dp_near[k])
        print (k,','.join(str(x) for x in v)+",",  ','.join(str(x) for x in querypos_minor[k])+",",  ','.join(str(x) for x in leftpos_major[k])+",",  ','.join(str(x) for x in leftpos_minor[k])+",",  ','.join(str(x) for x in seqpos_major[k])+",",  ','.join(str(x) for x in seqpos_minor[k])+",",  ','.join(str(x) for x in mapq_major[k])+",",  ','.join(str(x) for x in mapq_minor[k])+",", ','.join(str(x) for x in baseq_major[k])+",",  ','.join(str(x) for x in baseq_minor[k])+",",  ','.join(str(x) for x in baseq_major_near1b[k])+",", ','.join(str(x) for x in baseq_minor_near1b[k])+",", major_plus[k],major_minus[k],minor_plus[k],minor_minus[k], context1[k], context2[k], context1_count[k],context2_count[k],','.join(str(x) for x in mismatches_major[k])+",",','.join(str(x) for x in mismatches_minor[k])+",", major_read1[k],major_read2[k],minor_read1[k],minor_read2[k],np.mean(dp_near[k]),np.mean(dp_far[k]),p_value,file=fo)

    fo.close()







def option_check(opt):
    flag = True
    mutil.file_check(opt['input_pos'])
    mutil.file_check(opt['ref_fasta'])
    mutil.file_check(opt['ref_fasta']+".fai")
    return flag
