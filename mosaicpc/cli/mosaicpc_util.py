import os
import sys
from subprocess import *

def mkdir(dir_path):
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)

def strip_endslash(path1):
    path2 = ""
    if path1.endswith('/'):
        path2=path1[:-1]
    else:
        path2=path1
    return path2

def print_error(msg):
    print ("ERROR: " + msg)
    sys.exit(1)

def file_check(filepath):
    if not os.path.exists(filepath):
        print_error(filepath+" is not exist.")

def run_cmd(cmd):
    Popen(cmd, shell=True, stdout=PIPE).communicate()

def get_chr_sizes(fai_file):
    chr_sizes=dict()
    for line in open(fai_file):
        line=line.rstrip()
        fields=line.split('\t')
        chr_sizes[fields[0]]=chr_sizes.get(fields[0],fields[1])
    return chr_sizes

def print_log(opt,cont):
    if not opt['silence']:
        print (cont)

def arr_join(arr, delimiter=" "):
    cont = ""
    for a1 in arr:
        if cont != "":
            cont += delimiter
        cont += str(a1)
    return cont

def fileSave (path, cont, opt, gzip_flag = "n"):
    if gzip_flag == "gz":
        import gzip
        f = gzip.open(path, opt)
        f.write(cont)
        f.close()
    else:
        f = open (path, opt)
        f.write(cont)
        f.close