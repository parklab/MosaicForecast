import sys
from . import mosaicpc_util as mutil

def run(opt):
    print ("predict>run")
    if not option_check(opt):
        sys.exit(1)
    step1_predict(opt['input_feature_list'], opt['output'], opt['model'])

def option_check(opt):
    flag = True
    return flag

def step1_predict(input_feature_list, output, model):
    cmd = "Rscript mosaicpc/predict.R " + input_feature_list + " " + model + " " + output
    print (cmd)
    mutil.run_cmd(cmd)
