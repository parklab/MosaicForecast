#! /usr/bin/env python3

__format_spec__ = '.4g'

__info__ = {\
        'mappability':
        {'ID': 'MAPPABILITY', 'Number': 1, 'Type': 'Float', 'Description': '"UMAP mappability score at the candidate site (k=24)"'},
        'type':
        {'ID': 'TYPE', 'Number': 1, 'Type': 'String', 'Description': '"type of the candidate mutation (SNP, MNP, INS or DEL)"'},
        'length':
        {'ID': 'LENGTH', 'Number': 1, 'Type': 'Integer', 'Description': '"difference of base pair lengh of ref and alt allele for candidate sites"'},
        'GCcontent':
        {'ID': 'GCCONTENT', 'Number': 1, 'Type': 'Float', 'Description': '"20-bp local GCcontent"'},
        'querypos_p':
        {'ID': 'QUERYPOS_P', 'Number': 1, 'Type': 'Float', 'Description': "\"p-value or effect size by wilcoxon's rank sum test of base query positions of ref and alt alleles\""},
        'seqpos_p':
        {'ID': 'SEQPOS_P', 'Number': 1, 'Type': 'Float', 'Description': "\"p-value or effect size by wilcoxon's rank sum test of base sequencing cycles of ref and alt alleles\""},
        'baseq_p':
        {'ID': 'BASEQ_P', 'Number': 1, 'Type': 'Float', 'Description': "\"p-value or effect size by Wilcoxon's rank sum test of base qualities of ref and alt alleles\""},
        'baseq_t':
        {'ID': 'BASEQ_T', 'Number': 1, 'Type': 'Float', 'Description': "\"The test statistic under the large-sample approximation that the rank sum statistic is normally distributed (wilcox rank sum test of base qualites of alt alleles vs. ref alleles)\""},
        'context':
        {'ID': 'CONTEXT', 'Number': 1, 'Type': 'String', 'Description': '"three-nucleotide base context on the reads surrounding the mutant position"'},
        'sb_p':
        {'ID': 'SB_P', 'Number': 1, 'Type': 'Float', 'Description': "\"p-value or effect size by Fisher's exact test of strand bias for ref and alt alleles\""},
        'sb_read12_p':
        {'ID': 'SB_READ12_P', 'Number': 1, 'Type': 'Float', 'Description': "\"p-value or effect size by Fisher's exact test of read1/read2 bias for ref and alt alleles\""},
        'mosaic_likelihood':
        {'ID': 'MOSAIC_LIKELIHOOD', 'Number': 1, 'Type': 'Float', 'Description': '"mosaic genotype likelihood calculated (assuming uniform distribution of mosaics allele fraction from 0-1)"'},
        'het_likelihood':
        {'ID': 'HET_LIKELIHOOD', 'Number': 1, 'Type': 'Float', 'Description': '"Genotype likelihood of the variant being germline heterozygous"'},
        'refhom_likelihood':
        {'ID': 'REFHOM_LIKELIHOOD', 'Number': 1, 'Type': 'Float', 'Description': '"reference-homozygous genotype likelihood"'}
        }

__filter__ = {\
        'het':
        {'ID': 'het', 'Description': '"To be specified"'},
        'refhom':
        {'ID': 'refhom', 'Description': '"To be specified"'},
        'repeat':
        {'ID': 'repeat', 'Description': '"To be specified"'}}

__format__ = {\
        'ref_softclip':
        {'ID': 'REF_SOFTCLIP', 'Number': 1, 'Type': 'Float', 'Description': '"proportion of soft-clipped reads for ref reads"'},
        'alt_softclip':
        {'ID': 'ALT_SOFTCLIP', 'Number': 1, 'Type': 'Float', 'Description': '"proportion of soft-clipped reads for alt reads"'},
        'leftpos_p':
        {'ID': 'LEFTPOS_P', 'Number': 1, 'Type': 'Float', 'Description': "\"p-value or effect size by wilcoxon's rank sum test of left-most positions of ref and alt reads\""},
        'ref_baseq1b_p':
        {'ID': 'REF_BASEQ1B_P', 'Number': 1, 'Type': 'Float', 'Description': "\"p-value or effect size by Wilcoxon's rank sum test of base qualities from ref reads at mutant position, compared with base qualities from ref reads at 1bp downtream of the mutant position\""},
        'ref_baseq1b_t':
        {'ID': 'REF_BASEQ1B_T', 'Number': 1, 'Type': 'Float', 'Description': "\"The test statistic under the large-sample approximation that the rank sum statistic is normally distributed (wilcox rank sum test of base qualities from ref reads at mutant position, compared with base qualities from ref reads at 1bp downtream of the mutant position)\""},
        'alt_baseq1b_p':
        {'ID': 'ALT_BASEQ1B_P', 'Number': 1, 'Type': 'Float', 'Description': "\"p-value or effect size by Wilcoxon's rank sum test of base qualities from alt reads at mutant position, compared with base qualities from alt reads at 1bp downtream of the mutant position\""},
        'alt_baseq1b_t':
        {'ID': 'ALT_BASEQ1B_T', 'Number': 1, 'Type': 'Float', 'Description': "\"The test statistic under the large-sample approximation that the rank sum statistic is normally distributed (wilcox rank sum test of base qualities from alt reads at mutant position, compared with base qualities from alt reads at 1bp downtream of the mutant position)\""},
        'major_mismatches_mean':
        {'ID': 'MAJOR_MISMATCHES_MEAN', 'Number': 1, 'Type': 'Float', 'Description': '"average mismatches per ref reads"'},
        'minor_mismatches_mean':
        {'ID': 'MINOR_MISMATCHES_MEAN', 'Number': 1, 'Type': 'Float', 'Description': '"average mismatches per alt reads"'},
        'mismatches_p':
        {'ID': 'MISMATCHES_P', 'Number': 1, 'Type': 'Float', 'Description': "\"p-value or effect size by Wilcoxon's rank sum test of mismatches per ref reads vs. mismatches per alt reads\""},
        'mapq_p':
        {'ID': 'MAPQ_P', 'Number': 1, 'Type': 'Float', 'Description': "\"p-value or effect size by Wilcoxon's rank sum test of mapping qualities of ref and alt reads\""},
        'mapq_difference':
        {'ID': 'MAPQ_DIFFERENCE', 'Number': 1, 'Type': 'Float', 'Description': "\"difference of average map quality per alt reads vs. average map quality per ref reads\""},
        'AF':
        {'ID': 'AF', 'Number': 1, 'Type': 'Float', 'Description': '"variant allele fraction"'},
        'dp':
        {'ID': 'DP', 'Number': 1, 'Type': 'Integer', 'Description': '"read depth at mutant position"'},
        'dp_diff':
        {'ID': 'DP_DIFF', 'Number': 1, 'Type': 'Float', 'Description': '"difference of average read depths of local (<200bp) and distant (>2kb) regions"'},
        #'dp_p':
        #{'ID': 'DP_P', 'Number': 1, 'Type': 'Float', 'Description': "\"p-value or effect size by Wilcoxon's rank sum test of read depths sampled within 200bp window surrounding the mutant position vs. read depths sampled in distant regions from the mutant position (>2kb)\""},
        'conflict_num':
        {'ID': 'CONFLICT_NUM', 'Number': 1, 'Type': 'Integer', 'Description': "\"number of read pairs supporting both ref and alt alleles\""}
        }

import pandas as pd
import os.path
import subprocess
import tempfile
import os
import sys


def read_mf_predictions(mfpred_path):
    mfpred = pd.read_csv(mfpred_path, sep='\t')
    return(mfpred)


def file_format_metainfo():
    fileformat='VCFv4.3'
    val = '##fileformat=' + fileformat
    return(val)


def reference_metainfo(refpath='/home/attila/data/refgenome/GRCh37/dna/hs37d5.fa'):
    refpath = os.path.realpath(refpath)
    val = '##reference=file://' + refpath
    return(val)


def source_metainfo():
    val = '##source=MosaicForecast'
    return(val)


def contig_metainfo(reffaipath='/home/attila/data/refgenome/GRCh37/dna/hs37d5.fa.fai'):
    fai = pd.read_csv(reffaipath, sep='\t', names=['ID', 'length'], usecols=[0, 1])
    def helper(ix):
        ID = str(fai.iloc[ix, 0])
        length = str(fai.iloc[ix, 1])
        res = '##contig=<ID=' + ID + ',length=' + length + '>'
        return(res)
    val = [helper(i) for i in fai.index]
    val = '\n'.join(val)
    return(val)


def info_format_metainfo(kind='info'):
    if kind == 'info':
        dic = __info__
    elif kind == 'filter':
        dic = __filter__
    elif kind == 'format':
        dic = __format__
    def one_info_metainfo(mfid):
        def one_part(part):
            return(part + '=' + str(dic[mfid][part]))
        res = [one_part(y) for y in dic[mfid].keys()]
        res = ','.join(res)
        res = '##' + kind.upper() + '=<' + res + '>'
        return(res)
    res = [one_info_metainfo(k) for k in dic.keys()]
    res = '\n'.join(res)
    return(res)


def metainfo(refpath='/home/attila/data/refgenome/GRCh37/dna/hs37d5.fa'):
    '''
    Create meta-information lines (i.e. the header) of VCF

    Parameter(s):
    refpath: filepath to reference sequence

    Returns:
    all meta-information lines, separated by \n, in a single str
    '''
    if not os.path.isfile(refpath):
        error = 'Cannot find reference sequence at path: "' + refpath + '"'
        raise ValueError(error)
    reffaipath = refpath + '.fai'
    if not os.path.isfile(reffaipath):
        error = 'Cannot find reference sequence index at expected path: "' + reffaipath + '"'
        raise ValueError(error)
    reference = lambda : reference_metainfo(refpath)
    contig = lambda: contig_metainfo(reffaipath=reffaipath)
    infofun = lambda : info_format_metainfo(kind='info')
    filterfun = lambda : info_format_metainfo(kind='filter')
    formatfun = lambda : info_format_metainfo(kind='format')
    funs = [file_format_metainfo, source_metainfo, reference, contig, infofun,
            filterfun, formatfun]
    val = [f() for f in funs]
    val = '\n'.join(val) + '\n'
    return(val)


def chrom_alt_fields(mfpred):
    '''
    Create #CHROM, POS, REF, and ALT fields from "id"

    Parameter:
    mfpred: the pandas DataFrame from read_mf_predictions

    Returns:
    a pandas DataFrame containing the fields mentioned and, in addition, a
    sample field

    Details:
    There must be one and only one sample name in the sample field otherwise
    ValueError is raised.
    '''
    mat = [y.split(sep='~') for y in mfpred['id']]
    ncol = len(mat[0])
    fields = ['sample', '#CHROM', 'POS', 'REF', 'ALT']
    dtypes = ['category', 'category', 'int64', 'category', 'category']
    if ncol != len(fields):
        raise ValueError('"id" has unexpected number of subfields.')
    # transpose mat to make dictionary from columns
    dic = {fields[i]: [row[i] for row in mat] for i in range(ncol)}
    df = pd.DataFrame(dic)
    d = {y[0]: y[1] for y in zip(fields, dtypes)}
    df = df.astype(d)
    if len(df['sample'].cat.categories) != 1:
        raise ValueError('There must be one and only one sample.')
    def id_field():
        df1 = df[['sample', '#CHROM', 'POS']]
        df2 = pd.DataFrame({'ID': ['MISSING'] * len(df)})
        df3 = df[['REF', 'ALT']]
        res = pd.concat([df1, df2, df3], axis=1)
        return(res)
    df = id_field()
    return(df)


def qual_filter_fields(mfpred):
    '''
    Create QUAL, and FILTER fields from mfpred

    Parameter:
    mfpred: the pandas DataFrame from read_mf_predictions

    Returns:
    a pandas DataFrame containing the fields mentioned
    '''
    FILTER = [y.replace('mosaic', 'PASS') for y in mfpred['prediction']]
    df = pd.DataFrame({'QUAL': 'MISSING', 'FILTER': FILTER})
    return(df.astype('category'))


def info_fields(mfpred):
    '''
    Create the INFO fields from mfpred

    Parameter:
    mfpred: the pandas DataFrame from read_mf_predictions

    Returns:
    a pandas DataFrame containing the INFO fields separated by a ";"

    Details:
    The __info__  variable (a dict) is used to define the ID, Number, Type and
    Description of a particular INFO field.  The __format_spec__ variable is
    used to format all INFO fields of Type "float".
    '''
    def one_info_field(mfid='mappability'):
        ID = __info__[mfid]['ID']
        Type = __info__[mfid]['Type']
        def info_float(x):
            val = ID + '=' + str(format(x, __format_spec__))
            return(val)
        def info_integer(x):
            val = ID + '=' + str(x)
            return(val)
        def info_string(x):
            val = ID + '=' + x
            return(val)
        if Type == 'Float':
            helper = info_float
        elif Type == 'Integer':
            helper = info_integer
        elif Type == 'String':
            helper = info_string
        elif Type == 'Flag':
            raise ValueError('Flag INFO field type is not yet implemented')
        data = [helper(y) for y in mfpred[mfid]]
        return(data)
    fields = __info__.keys()
    val = [one_info_field(mfid=f) for f in fields]
    val = [';'.join(y) for y in zip(*val)]
    df = pd.DataFrame({'INFO': val})
    return(df)


def format_sample_fields(mfpred):
    '''
    Create the FORMAT and sample fields from mfpred

    Parameter:
    mfpred: the pandas DataFrame from read_mf_predictions

    Returns:
    a pandas DataFrame containing the FORMAT fields separated by a ";" and the
    corresponding sample field

    Details:
    The __format__  variable (a dict) is used to define the ID, Number, Type and
    Description of a particular FORMAT field.  The __format_spec__ variable is
    used to format all FORMAT fields of Type "float".
    '''
    def one_sample_field(mfpred, mfid='AF'):
        ID = __format__[mfid]['ID']
        Type = __format__[mfid]['Type']
        def sample_float(x):
            return(str(format(x, __format_spec__)))
        def sample_general(x):
            return(str(x))
        if Type == 'Float':
            helper = sample_float
        else:
            helper = sample_general
        val = [helper(y) for y in mfpred[mfid]]
        return(val)
    fields = __format__.keys()
    ids = [__format__[f]['ID'] for f in fields]
    FORMAT = ':'.join(ids)
    vals = [one_sample_field(mfpred, mfid=f) for f in fields]
    vals = [':'.join(y) for y in zip(*vals)]
    df = pd.DataFrame({'FORMAT': FORMAT, 'SAMPLE': vals})
    return(df)



def data_lines(mfpred):
    '''
    Creates all data lines (i.e. the body) of the VCF

    Parameter:
    mfpred: the pandas DataFrame from read_mf_predictions

    Returns:
    a pandas DataFrame containing all fields
    '''
    funs = [chrom_alt_fields, qual_filter_fields, info_fields, format_sample_fields]
    l = [f(mfpred) for f in funs]
    df = pd.concat(l, axis=1)
    sample = df['sample'][0]
    df1 = df[df.columns[1:-1]]
    df2 = pd.DataFrame({sample: df['SAMPLE']}, index=range(len(df)))
    val = pd.concat([df1, df2], axis=1)
    return(val)


def main(mfpred_path, refpath, vcfgz_path):
    mfpred = read_mf_predictions(mfpred_path=mfpred_path)
    head = metainfo(refpath=refpath)
    body_df = data_lines(mfpred)
    body = body_df.to_csv(sep='\t', header=True, index=False)
    VCF = head + body
    args0 = ['bcftools', 'sort', '-Oz', '-o', vcfgz_path]
    subprocess.run(args0, input=VCF, encoding="utf-8")
    args1 = ['bcftools', 'index', '--tbi', vcfgz_path]
    subprocess.run(args1)
    return(None)


if __name__ == '__main__':
    mfpred_path = sys.stdin
    if len(sys.argv) != 3:
        print(
    '''
    Convert MosaicForecast output MF.predictions into a sorted, gzipped, indexed VCF

    Usage:

    mf2vcf.py <input.MF.predictions refseq.fa output.vcf.gz
    or
    cat input.MF.predictions | mf2vcf.py refseq.fa output.vcf.gz


    Details:

    Currently gzip compressed VCF is the only supported output type.

    For each INFO, FILTER, and FORMAT field the appropriate ID,
    Number, Type, and 'Description' must be specified by editing the __info__,
    __filter__, and __format__ dictionaries in mf2vcf.py.

    The pandas Python package must be installed.
    '''
        )
    else:
        refpath = str(sys.argv[1])
        vcfgz_path = str(sys.argv[2])
        main(mfpred_path, refpath, vcfgz_path)
