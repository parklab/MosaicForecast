# MosaicPC-Calico
`Mosaicpc` is a mosaic SNV detecting software based on phasing and random forest

## Installation

Requirements:

* Python version 2.7/3.4+
* Python packages : `collections`, `itertools`, `subprocess`, `multiprocessing`, `regex (2017.9.23)`, `numpy (1.12.1)`, `pyfaidx (0.5.3.1)`, `pysam (0.11.2.2)`, `pysamstats (1.0.1)`, `scipy (0.19.0)`, `pandas (0.20.1)`
 
* R version 3.4.1
* R packages : `caret (6.0-78)`, `e1071 (1.6-8)`, `glmnet (2.0-13)`, `nnet (7.3-12)`, `randomForest (4.6-14)`

MosaicPC is available through pypi. To install, type:

```sh
$ pip install mosaicpc
```


## Usage (Command line interface)

The `mosaicpc` package provides command line tools for phasing, modeling, and predicting variants.

```sh
$ mosaicpc [subcommand] [options]
```

```
usage: mosaicpc <sub-command> [options]

mosaicpc 0.1 (2018/10/09): A mosaic SNV detecting software based on phasing and
random forest

optional arguments:
  -h, --help       show this help message and exit
  -v, --version    show program's version number and exit

sub-commands:

    phase          read-based phasing
    extract_feature
                   extract features
    predict        predict phasing of variants
    model          training
```


### 1. Phasing:

**Usage:** 

```sh
$ mosaicpc phase 

usage: mosaic phase [-h] [-b BAM_DIR] [-p INPUT_POS] [-o OUTPUT_DIR]
                   [-r REF_FASTA] [-n N_JOBS] [-d MIN_DP_INFORSNPS]
                   [-m MAPQ_THRESHOLD] [-q BQSE_QUALITY_THRESHOLD]
                   [-l LOG_FILE] [-s]

read-based phasing

optional arguments:
  -h, --help            show this help message and exit
  -b BAM_DIR, --bam_dir BAM_DIR
                        directory path which contains bam files
  -p INPUT_POS, --input_pos INPUT_POS
                        file path which include target positions. (file
                        format:chr pos-1 pos ref alt sample, sep=\t)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        directory path for outputs (default: mosaicpc_output)
  -r REF_FASTA, --ref_fasta REF_FASTA
                        reference fasta file
  -n N_JOBS, --n_jobs N_JOBS
                        remove bsub command for running on Orchestra server
                        (default: 1)
  -d MIN_DP_INFORSNPS, --min_dp_inforSNPs MIN_DP_INFORSNPS
                        minimum depth for informative SNPs (default:20)
  -m MAPQ_THRESHOLD, --mapq MAPQ_THRESHOLD
                        mapping quality score threshold (default:20)
  -q BQSE_QUALITY_THRESHOLD, --bqse_quality BQSE_QUALITY_THRESHOLD
                        base quality score threshold (default:20)
  -l LOG_FILE, --log LOG_FILE
                        log file name (default: multiple_inforSNPs.log)
  -s, --silence         silence
```


**Note:** 

1. Name of bam files should be "sample.bam" under the `BAM_DIR`. 
2. There should be a fai file under the same dir of the fasta file (`samtools faidx input.fa`). 
3. The `MIN_DP_INFORSNPS` is the minimum depth of coverage of trustworthy neaby het SNPs.

```
hap=2: hete
hap=3: mosaic
hap>=4: cnv/repeat
```

**Input files:**

**Output files:**

**Demo:**

```
mosaicpc phase -b demo -o demo/phasing -r ${reference_dir}/human_g1k_v37_decoy.fasta -n 2 -p demo/test.input -d 20
```



### 2. Feature extraction:

**Usage:** 

```
$ mosaicpc extract_feature

usage: mosaicpc extract_feature [-h] [-b BAM_DIR] [-p INPUT_POS]
                             [-o OUTPUT_FEATURE] [-r REF_FASTA] [-s]

optional arguments:
  -h, --help            show this help message and exit
  -b BAM_DIR, --bam_dir BAM_DIR
                        directory path which contains bam files
  -p INPUT_POS, --input_pos INPUT_POS
                        file path which include target positions. (file
                        format:chr pos-1 pos ref alt sample, sep=\t)
  -o OUTPUT_FEATURE, --output_feature OUTPUT_FEATURE
                        directory path for outputs (default: mosaicpc_output)
  -r REF_FASTA, --ref_fasta REF_FASTA
                        reference fasta file
  -s, --silence         silence
```

**Note:** 

1. Name of bam files should be "sample.bam" under the `BAM_DIR`. 
2. There should be a fai file under the same dir of the fasta file (`samtools faidx input.fa`). 
3. We did not use dbSNP AF as an feature because we only focus on ultra-rare mosaic mutations, but you can use it to train your model if you have interest in common variants.

**Input files:**

**Ouput files:**

**Demo:**

```
mosaicpc extract_feature -p demo/test.input -o demo/test.features -b demo -r ${reference_dir}/human_g1k_v37_decoy.fasta
```


### 3. Prediction:

**Usage:** 

```sh
$ mosaicpc predict
usage: mosaicpc predict [-h] [-i INPUT_FEATURE_LIST] [-m MODEL] [-o OUTPUT] [-s]

predict phasing of variants

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FEATURE_LIST, --input INPUT_FEATURE_LIST
                        input feature list
  -m MODEL, --model MODEL
                        model
  -o OUTPUT, --output OUTPUT
                        output
  -s, --silence         silence
```

**Demo:**

```
moasicpc predict -i demo/test_feature_list_R -m models_trained/250x_rf_PCAandPhase_30mtry.rds -o demo/test_predictions
```


> You may use our models trained with WGS data (150bp paired-end data, given different read depths):
>
> * 50x\_rf\_PCAandPhase\_30mtry.rds
> * 100x\_rf\_PCAandPhase\_30mtry.rds
> * 150x\_rf\_PCAandPhase\_30mtry.rds
> * 200x\_rf\_PCAandPhase\_30mtry.rds
> * 250x\_rf\_PCAandPhase\_30mtry.rds

### 4. Make model:
You could also train models using your own data.

#### Method1. train models only based on phasable sites:

**Usage:** 

```
$ mosaicpc model
```

**Demo:**

```
mosaicpc model -f demo/all_putative_mosaics_feature.list.features_addphasing_addvalidation -o demo/test.rds 
```



#### Method2. correct genotypes from phasing and train models based on corrected phasable sites:

Recommended when you have >=150 phasable sites validated orthogonally.

##### 1st step: correct genotype labels of phasable sites using experimentally validated sites:

**Usage:** 

```
Rscript phasing_correction_train.R trainset prediction_model_phasingcorrection output_file_phasingcorrected
```

**Demo:**

```
Rscript phasing_correction_train.R demo/all_putative_mosaics_feature.list.features_addphasing_addvalidation demo/prediction_model_phasingcorrection.rds demo/test_phasingcorrected
```

##### 2nd step: train based on corrected phasable sites:
**Usage:** 

```
Rscript train_on_correctedphasing.R trainset prediction_model
```

**Demo:**

```
Rscript train_on_correctedphasing.R demo/test_phasingcorrected demo/test2.rds
```

