# -*- coding: utf-8 -*-
import logging
import sys
import argparse

from . import phase
from . import extract_feature
from . import predict
from . import mosaicpc_util as mutil

VERSION = "1.0"
VERSION_DATE = "2018/10/09"
opt_silence = False
PROG = 'mosaicpc'


def pl(cont):
    mutil.print_log({'silence':opt_silence},cont)


def print_option(opt):
    pl("#========== option ==========")
    for k in ['subcommand']:
        pl("# " + k + " : " + str(opt[k]))
    for k in opt.keys():
        if k != 'subcommand':
            pl("# " + k + " : " + str(opt[k]))
    pl("#============================")


def dispatch_job(opt):
    if "phase" == opt['subcommand']:
        phase.run(opt)
    if "extract_feature" == opt['subcommand']:
        extract_feature.run(opt)
    if "predict" == opt['subcommand']:
        predict.run(opt)
    pass

def cli():
    global VERSION, VERSION_DATE, opt_silence, PROG

    parser = argparse.ArgumentParser(usage='%(prog)s <sub-command> [options]', description='%(prog)s '+VERSION+" ("+VERSION_DATE+")"+': A mosaic SNV detecting software based on phasing and random forest')
    parser.add_argument('-v','--version', action='version', version='%(prog)s '+VERSION+" ("+VERSION_DATE+")")

    ################################################################
    ##### PHASING #####
    ################################################################
    subparsers = parser.add_subparsers(title="sub-commands", dest="subcommand",metavar='', prog=PROG)
    p1 = subparsers.add_parser('phase', help='read-based phasing', description='read-based phasing')
    p1.add_argument('-b','--bam_dir', dest='bam_dir', default=None, help="directory path which contains bam files")
    p1.add_argument('-p','--input_pos', dest='input_pos', default="", help="file path which include target positions. (file format:chr pos-1 pos ref alt sample, sep=\\t)")
    p1.add_argument('-o','--output_dir', dest='output_dir', default="mosaicpc_output", help="directory path for outputs (default: mosaicpc_output)")
    p1.add_argument('-r','--ref_fasta', dest='ref_fasta', default="./", help="reference fasta file")
    p1.add_argument('-n','--n_jobs', dest='n_jobs', type=int, default=1, help='remove bsub command for running on Orchestra server (default: 1)')
    p1.add_argument('-d','--min_dp_inforSNPs', dest='min_dp_inforSNPs',type=int, default=20, help="minimum depth for informative SNPs (default:20)")
    p1.add_argument('-m','--mapq', dest='mapq_threshold',type=int, default=20, help="mapping quality score threshold (default:20)")
    p1.add_argument('-q','--bqse_quality', dest='bqse_quality_threshold',type=int, default=20, help="base quality score threshold (default:20)")
    p1.add_argument('-l','--log', dest='log_file', default="multiple_inforSNPs.log", help="log file name (default: multiple_inforSNPs.log)")
    p1.add_argument('-s','--silence', dest='silence', action='store_true', default=False, help='silence')

    ################################################################
    ##### FEATURE EXTRACT #####
    ################################################################
    p1 = subparsers.add_parser('extract_feature', help='extract features')
    p1.add_argument('-b','--bam_dir', dest='bam_dir', default=None, help="directory path which contains bam files")
    p1.add_argument('-p','--input_pos', dest='input_pos', default="", help="file path which include target positions. (file format:chr pos-1 pos ref alt sample, sep=\\t)")
    p1.add_argument('-o','--output_feature', dest='output_feature', default="mosaicpc_output.feature", help="directory path for outputs (default: mosaicpc_output)")
    p1.add_argument('-r','--ref_fasta', dest='ref_fasta', default="./", help="reference fasta file")
    p1.add_argument('-s','--silence', dest='silence', action='store_true', default=False, help='silence')
    
    ################################################################
    ##### PREDICT #####
    ################################################################
    MODEL_LIST = ["50x_rf_PCAandPhase_30mtry.rds","100x_rf_PCAandPhase_30mtry.rds","150x_rf_PCAandPhase_30mtry.rds","200x_rf_PCAandPhase_30mtry.rds","250x_rf_PCAandPhase_30mtry.rds"]

    p1 = subparsers.add_parser('predict', help='predict phasing of variants', description='predict phasing of variants')
    p1.add_argument('-i','--input', dest='input_feature_list', default=None, help='input feature list')
    p1.add_argument('-m','--model', dest='model', default="", help='model')
    p1.add_argument('-o','--output', dest='output', default="", help="output")
    p1.add_argument('-s','--silence', dest='silence', action='store_true', default=False, help='silence')
    
    ################################################################
    ##### MODELING #####
    ################################################################
    p1 = subparsers.add_parser('model', help='training', description='training')
    p1.add_argument('-conf', dest='conf', default=None, help='configuration file')
    p1.add_argument('-aligner', dest='aligner', choices=["bwa"], default="bwa", help='aligner')
    p1.add_argument('-wd', dest='wd', default="./", help="working space (./)")
    p1.add_argument('-out', dest='out', default=None, help="out script file (default:print screen)")
    p1.add_argument('-s','--silence', dest='silence', action='store_true', default=False, help='silence')
    
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1][0] != '-'):
        sys.argv.append('-h')
    opt = vars(parser.parse_args())
    opt_silence = opt['silence']
    print_option(opt)
    dispatch_job(opt)





