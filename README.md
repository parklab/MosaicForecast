# MosaicPC-Calico
A mosaicSNV detecting software based on phasing and random forest

# Versions:
Python version 3.6.1
R version 3.4.1

# Packages need to be installed:
## Python packages need to import or install:
os\
sys\
collections\
itertools\
subprocess\
multiprocessing\
regex (2017.9.23)\
numpy (1.12.1)\
pyfaidx (0.5.3.1)\
pysam (0.11.2.2)\
pysamstats (1.0.1)\
scipy (0.19.0)\
pandas (0.20.1)
## R packages need to import or install:
caret\
e1071\
glmnet\
nnet

# Usage:
## Phasing:
python phase\_WGS.python3.6.1\_multiplesamples.py\
Usage: python phase\_WGS.py bam\_dir data\_dir ref\_fasta n\_jobs input\_positions(file format:sample chr pos ref alt, sep=\t) min\_dp\_inforSNPs(int, can be set to 20)
## Feature extraction:
### 1st step:
python feature\_extraction.python3.6.1.py\
Usage: python(v3.6.1) feature\_extraction.py input\_bed(file\_format: chr pos-1 pos ref alt sample, sep="\t") output\_features bam\_dir reference\_fasta\
Note: 1. name of bam files are "sample.bam" by default. 2.there should be a fai file under the same dir of the fasta file (samtools faidx input.fa) 3. we did not use dbSNP AF as an feature, but you can use it to train your model if you have interest in common variants.
### 2nd step:
Rscript feature\_extraction.R\
Usage: Rscript feature\_extraction.R input\_file(feature\_list\_frompython) output\_file(feature\_list\_R)

## Training:
### Method1. train based on phasable sites:\
hap=2: het\
hap=3: mosaic\
hap>3: cnv/repeat\

### Method2. train based on validated and phasable sites:
1st step:\

2nd step:\

## Prediction:
Rscript prediction.R\
Error: Rscript prediction.R input\_file(feature\_list) model\_trained output\_file(predictions)

