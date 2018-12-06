import sys
import pysam
from multiprocessing import Pool
from . import mosaicpc_util as mutil
import pysamstats
import scipy.stats
from collections import OrderedDict, defaultdict
import regex as re
from pandas import *
import itertools
opt_silence = False

def pl(cont):
    mutil.print_log({'silence':opt_silence},cont)

def run(opt):
    global opt_silence

    opt_silence = opt['silence']
    pl("checking option")
    if not option_check(opt):
        sys.exit(1)

    pl("running phase > step1 make seperate bams")
    step1_mk_seperate_bams(opt['output_dir'], opt['input_pos'], opt['bam_dir'], opt['n_jobs'], opt['mapq_threshold'], opt['bqse_quality_threshold'])
    
    pl("running phase > step2 extract inforSNPs")
    conflict_mosaic = step2_extract_inforSNPs(opt['output_dir'], opt['ref_fasta'], opt['min_dp_inforSNPs'])

    pl("running phase > step3 make 2x2table")
    step3_mk_2x2table(opt['output_dir'], opt['ref_fasta'])

    pl("running phase > step4 assign phasing state")
    (inforSNPs,phasing_2by2) = step4_assign_phasing_state(opt['output_dir'])

    pl("running phase > step5 further phasing")
    step5_further_phasing(inforSNPs,phasing_2by2, conflict_mosaic, opt['output_dir'], opt['log_file'],opt['bam_dir'])

def option_check(opt):
    flag = True
    opt['bam_dir'] = mutil.strip_endslash(opt['bam_dir'])
    opt['output_dir'] = mutil.strip_endslash(opt['output_dir'])
    mutil.mkdir(opt['output_dir'])
    mutil.file_check(opt['ref_fasta'])
    mutil.file_check(opt['ref_fasta']+".fai")
    mutil.file_check(opt['input_pos'])
    return flag

##input_bam=bam_dir+"/"+sample+".bam"
##a=pysam.AlignmentFile(input_bam, "rb")
##head Vaccarino_brain.mutect2/final.passed.tsv
##1       2533604 G       A       Vaccarino_brain

##1st step: generate seperate bams supporting ref allale and alt alleles:
def step1_mk_seperate_bams(output_dir, input_pos, bam_dir, n_jobs, mapq_threshold, bqse_quality_threshold):
    cmd_list1=list()
    cmd_list2=list()
    fo1=open(output_dir+"/all_candidates","w")
    for line in open(input_pos):
        #1       1004864 1004865 G       C       test
        line=line.rstrip()
        fields=line.split('\t')
        sample=fields[5]
        chr=fields[0]
        chrom=str(chr)
        pos=int(fields[1])+1
        start=int(fields[1])
        end=int(fields[2])
        major_allele=fields[3]
        minor_allele=fields[4]
        input_bam=bam_dir+"/"+sample+".bam"
        a=pysam.AlignmentFile(input_bam, "rb")
        name=sample+'_'+chr+'_'+str(pos)
        major_ids=list()
        minor_ids=list()
        major_num=int(0)
        minor_num=int(0)

        pl ("\tprocessing... "+line)
        for pileupcolumn in a.pileup(chrom, start, end):
            for pileupread in pileupcolumn.pileups:
                if pileupread.alignment.mapping_quality>=mapq_threshold:
                    # print ("MQ:",pileupread.alignment.mapping_quality)
                    try:
                        flag_queryscore_th = True
                        if pileupread.indel==len(minor_allele)-len(major_allele) and pileupread.indel > 0: ### In>0
                            querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+(len(minor_allele)-len(major_allele))+1]
                            queryscore=pileupread.alignment.query_qualities[pileupread.query_position:pileupread.query_position+(len(minor_allele)-len(major_allele))+1]
                            # print (pileupread.indel, pileupread.query_position, querybase,major_allele, minor_allele, pileupcolumn.pos, pos-1,queryscore)
                        elif pileupread.indel==len(minor_allele)-len(major_allele) and pileupread.indel<0: ### Del<0
                            querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(minor_allele)]
                            queryscore=pileupread.alignment.query_qualities[pileupread.query_position:pileupread.query_position+len(minor_allele)]
                            # print (pileupread.indel, pileupread.query_position, querybase,major_allele, minor_allele, pileupcolumn.pos, pos-1,queryscore)
                            # print (pileupread.alignment.query_sequence, pileupread.alignment.query_name, pileupread.alignment.cigar)
                            
                            # print (start, end)
                        elif pileupread.indel == 0: ### No InDel read
                            if len(minor_allele)>=len(major_allele): ### SNV & insertion region
                                querybase=pileupread.alignment.query_sequence[pileupread.query_position]
                                queryscore=pileupread.alignment.query_qualities[pileupread.query_position]
                                if queryscore < bqse_quality_threshold:
                                    flag_queryscore_th = False
                            elif len(minor_allele)<len(major_allele): ### deletion region
                                querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+(end-start)]
                                queryscore=pileupread.alignment.query_qualities[pileupread.query_position:pileupread.query_position+(end-start)]
                            if int(pileupcolumn.pos)==int(pos)-1:
                                pass
                                # print (pileupread.indel, pileupread.query_position, querybase,major_allele, minor_allele, pileupcolumn.pos, pos-1,queryscore)
                        if flag_queryscore_th:
                            if int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(major_allele): 
                                major_ids.append(pileupread.alignment.query_name)
                                major_num+=1
                            elif int(pileupcolumn.pos)==int(pos)-1 and str(querybase)==str(minor_allele):
                                minor_ids.append(pileupread.alignment.query_name)
                                minor_num+=1
                    except:
                        continue

        start=int(pos)-1000
        end=int(pos)+1000
        conflictnum=0
        # if len(major_ids)>=2 and len(minor_ids)>=2:
        if True:
            outbam_major = output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.major.bam"
            outbam_minor = output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.minor.bam"
            outbam_merged = output_dir+"/"+sample+"."+str(chr)+"_"+str(pos)+".mosaic.merged.bam"
            f1=pysam.AlignmentFile(outbam_major,"wb",template=a)
            f2=pysam.AlignmentFile(outbam_minor,"wb",template=a)
            f3=pysam.AlignmentFile(outbam_merged,"wb",template=a)

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
            pl("\tsaved "+outbam_major)
            pl("\tsaved "+outbam_minor)
            pl("\tsaved "+outbam_merged)

            outbam_major_sort = outbam_major.replace(".mosaic.major.bam", ".mosaic.major.sorted.bam")
            outbam_minor_sort = outbam_minor.replace(".mosaic.minor.bam", ".mosaic.minor.sorted.bam")
            outbam_merged_sort = outbam_merged.replace(".mosaic.merged.bam", ".mosaic.merged.sorted.bam")
            cmd_list1.append("samtools sort "+outbam_major+" -o "+outbam_major_sort)
            cmd_list1.append("samtools sort "+outbam_minor+" -o "+outbam_minor_sort)
            cmd_list1.append("samtools sort "+outbam_merged+" -o "+outbam_merged_sort)
            cmd_list2.append("samtools index "+outbam_major_sort)
            cmd_list2.append("samtools index "+outbam_minor_sort)
            cmd_list2.append("samtools index "+outbam_merged_sort)

        print (sample,chr,pos,major_allele,minor_allele,major_num,minor_num,conflictnum,file=fo1)
        pl('\t> ' + mutil.arr_join([sample,chr,pos,major_allele,minor_allele,major_num,minor_num,conflictnum]))
    fo1.close()
    pl("\tsaved "+output_dir+"/all_candidates")

    pl("\tsorting... bam files")
    pool=Pool(processes=n_jobs)
    pool.map(mutil.run_cmd,cmd_list1,1)
    pool.close()
    pool.join()

    pl("\tindexing... sorted_bam files")
    pool=Pool(processes=n_jobs)
    pool.map(mutil.run_cmd,cmd_list2,1)
    pool.close()
    pool.join()


## 2nd step: extract candidate nearby inforSNPs (germline het):
## output : all.merged.inforSNPs.pos
def step2_extract_inforSNPs(output_dir, ref_fasta, min_dp_inforSNPs):
    merged_inforSNPs=open(output_dir+"/all.merged.inforSNPs.pos","w")
    conflict_mosaic=dict()
    #==> /n/data1/hms/dbmi/park/yanmei/simulated_bams_na12878/12878-50x.merged.inforSNPs.pos <==
    #1 18089 G T 1 17385 G A 0
    for line in open(output_dir+"/all_candidates"):
        line=line.rstrip()
        # pl ("\tprocessing... "+line)

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
            f3_alignment_file =pysam.AlignmentFile(f3_sorted_name,'rb')
            for rec in pysamstats.stat_variation(f3_alignment_file, fafile=ref_fasta, min_mapq=20, min_baseq=20):
                if  ([rec['A'],rec['C'],rec['G'],rec['T']].count(0)<=2) and (rec['reads_pp']>10):
                    # print (rec) 
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
                    # print (max_base, max_num, alt_base, alt_num)
                    # print (max_base, max_num, alt_base, alt_num, scipy.stats.binom_test(alt_num, max_num+alt_num), rec['pos'], pos)
                    if alt_num>=3 and max_num>=3 and max_num+alt_num>=min_dp_inforSNPs and scipy.stats.binom_test(alt_num, max_num+alt_num)>0.05 and rec['pos']!=int(pos):
                        pl ("\t\t>infoSNP: " + mutil.arr_join([rec['pos'],max_base, max_num, alt_base, alt_num, scipy.stats.binom_test(alt_num, max_num+alt_num)]) )
                        print (sample,chr,pos,major_allele,minor_allele,chr,int(rec['pos'])+1,max_base,alt_base,conflictnum, file=merged_inforSNPs)
    merged_inforSNPs.close()
    pl("\tsaved " + output_dir+"/all.merged.inforSNPs.pos")
    return conflict_mosaic

## 3rd step: generate the 2X2 table:
def step3_mk_2x2table(output_dir, ref_fasta):
    n_major_het1=dict()
    n_major_het2=dict()
    n_minor_het1=dict()
    n_minor_het2=dict()

    chr_sizes = mutil.get_chr_sizes(ref_fasta+'.fai')

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

    pl("\tsaved " + output_dir+"/all_2x2table")


## last step: assign phasing state to each site
## output : all_2x2table
def step4_assign_phasing_state(output_dir):
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
            name=sample+";"+chr+";"+pos+";"+ref+";"+alt
            inforSNPs[name].append(sample+";"+chr+";"+pos_inforSNP+";"+het1_inforSNP+";"+het2_inforSNP)
        
            phase[name]['hap=2']=phase[name].get("hap=2",0)+0
            phase[name]['hap=3']=phase[name].get("hap=3",0)+0
            phase[name]['hap>3']=phase[name].get("hap>3",0)+0
            if C1+C2+C3+C4>=10:
                #UMB1349 18 30348609 C T 30348301 G A 0 34 29 0 2

                if not ( ((C1>C2*10) and (C4>C3*10)) or ((C1<C2/10) and (C4<C3/10)) or (((C1>C2/10) and (C1<C2*10)) and (C3<C4/10 or C4<C3/10))  ):
                    phase[name]['hap>3']=phase[name].get("hap>3",0)+1
                elif (((C1>C2/10) and (C1<C2*10)) and (C3<C4/10 or C4<C3/10)):
                    phase[name]['hap=3']=phase[name].get("hap=3",0)+1
                elif ((C1>C2*10) and (C4>C3*10)) or ((C1<C2/10) and (C4<C3/10)):
                    phase[name]['hap=2']=phase[name].get("hap=2",0)+1

    phasing_2by2 = dict()
    fo=open(output_dir+"/all.phasing_2by2","w")
    for k,v in sorted(phase.items()):
        phasing_2by2[k]=max(v,key=v.get)
        print (' '.join(str(x) for x in k.split(";")), v,max(v, key=v.get), file=fo)
    fo.close()

    pl("\tsaved " + output_dir+"/all.phasing_2by2")

    return (inforSNPs,phasing_2by2)

#phasing_2by2 = dict()
#for k,v in sorted(phase.items()):
#	phasing_2by2[k]=max(v,key=v.get)

#12878-250x 10 10002280 A G 10002816 G A 0 18 0 0 1

## one additional step: using multiple germline het sites to do further phasing
## output : all.phasing
def step5_further_phasing(inforSNPs,phasing_2by2,conflict_mosaic, output_dir, log_file, bam_dir):
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
                # df=df.reindex_axis(sorted(df.columns,key=lambda x:int(x[0:])), axis=1)
                df=df.reindex(sorted(df.columns,key=lambda x:int(x[0:])), axis=1)
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
    pl("\tsaved " + output_dir+"/"+log_file)

    fo4=open(output_dir+"/all.phasing","w")
    print ("sample","chr","pos","ref","alt","phasing","conflicting_reads",file=fo4)
    for k,v in sorted(phasing_2by2.items()):
        if k in MT2_phasing_num:
            if MT2_phasing_num[k]['wrong']>=MT2_phasing_num[k]['correct']:
    #		if max(MT2_phasing_num[k], key=MT2_phasing_num[k].get)=="wrong":
                v="hap>3"
        print (' '.join(k.split(';')), v, conflict_mosaic[k], file=fo4)
    fo4.close()
    pl("\tsaved " + output_dir+"/all.phasing")


if __name__ == "__main__":
    pass


