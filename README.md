# MosaicForecast
A machine learning method that leverages read-based phasing and read-level features to accurately detect mosaic single-nucleotide variants (SNVs) from NGS data.

# Required Interpreter Versions:
* Python version 3.4+
* R version 3.4+

# Packages need to be installed:
## Python packages:
* collections
* itertools
* subprocess
* multiprocessing
* regex (2017.9.23)
* numpy (1.12.1)
* pyfaidx (0.5.3.1)
* pysam (0.11.2.2)
* pysamstats (1.0.1)
* scipy (0.19.0)
* pandas (0.20.1)
## R packages:
* caret (6.0-78)
* e1071 (1.6-8)
* glmnet (2.0-13)
* nnet (7.3-12)
* ggbiplot (0.55)

# Usage:
## Phasing:
**Usage:** 

```
python Phase.py bam_dir output_dir ref_fasta n_jobs_parallel input_positions min_dp_inforSNPs
```

**Note:** 

1. Name of bam files should be "sample.bam" under the bam\_dir. 
2. There should be a fai file under the same dir of the fasta file (samtools faidx input.fa).
3. File format of the input\_positions: chr pos-1 pos ref alt sample, sep=\t 
4. The "min\_dp\_inforSNPs" is the minimum depth of coverage of trustworthy neaby het SNPs, can be set to 20.

**Output:**
```
output_dir/all.phasing
hap=2: het
hap=3: mosaic
hap>3: cnv/repeat
```

**Demo:**

```
python Phase.py demo demo/phasing ${reference_dir}/human_g1k_v37_decoy.fasta 1 demo/test.input 20
```

## Extraction of read-level features:
**Usage:**

### 1st step:
```
python ReadLevel_Features_extraction.py input.bed output.features bam_dir reference_fasta
```
**Note:**
1. Name of bam files are "sample.bam" by default.
2. There should be a fai file under the same dir of the fasta file (samtools faidx input.fa).
3. File format of the input.bed: chr pos-1 pos ref alt sample, sep=\t 

### 2nd step:
```
Rscript Rscript ReadLevel_Features_extraction.R input_file output_file read_length(integer) type(pvalue||effectsize)
```
**Note:**
1. Use "pvalue" when your data has relatively even read coverage (i.e. WGS data) or the training sample size is big (i.e., >10000 sites);
2. Use "effectsize" when your data has extrmely un-even read coverage and small training sample size. The "effectsize" mode is relatively slow.
3. The input\_file is the file from the 1st step.

**Output:**
```
A list of read-level features for each input site.
```
**Demo:**
```
python ReadLevel_Features_extraction.py demo/test.input demo/test.features demo ${reference_fasta}
Rscript ReadLevel_Features_extraction.R demo/test.features demo/test.features_R 150 pvalue 
```


## Prediction:

**Usage:**

Rscript Prediction.R input\_file(feature\_list) model\_trained output\_file(predictions)


**Demo:**
Rscript Prediction.R demo/test.features\_R models\_trained/brain\_MT2-PON.250x\_MosaicForecast-Refine\_pvalue.rds  demo/test\_predictions

> You may use our models trained with brain WGS data (paired-end read at 50-250X read depths, we train our models based on Mutect2-PON callings):
>
> * models\_trained/50-200X\_brainWGS/brain\_MT2-PON.50x.rds
> * models\_trained/50-200X\_brainWGS/brain\_MT2-PON.100x.rds
> * models\_trained/50-200X\_brainWGS/brain\_MT2-PON.150x.rds
> * models\_trained/50-200X\_brainWGS/brain\_MT2-PON.200x.rds
> * models\_trained/brain\_MT2-PON.250x\_MosaicForecast-Refine\_pvalue.rds

## You could also train models using your own data:
**Usage:**
```
Rscript Train_RFmodel.R input(trainset) output(prediction_model) type(Phase|Refine)
```
**Note:** 

1. You could choose to train your model based on Phasing (hap=2, hap=3, hap>3, type in "Phase") or Refined genotypes ("mosaic","het","refhom","repeat", type in "Refine").
2. The input file should be a list of pre-generated read-level features, adding a column termed "phase" (Phase model) or "phase\_model\_corrected" (Refined genotypes model). 


**Demo:**
```
Rscript Train_RFmodel.R demo/phasable_trainset demo/Phase_model.rds Phase
Rscript Train_RFmodel.R demo/phasable_trainset demo/Refine_model.rds Refine
```

## Convert phasing to four-category genotypes based on experimental data:

#### (Recommended when you have >=100 hap>=3 sites validated orthogonally)
#### 1st step: correct genotype labels of phasable sites using experimentally validated sites:
```
Rscript phasing_correction_train.R
```
**Usage:** Rscript training\_phasing\_correction.R trainset prediction\_model\_phasingcorrection output\_file\_phasingcorrected

**Demo:**\
Rscript phasing\_correction\_train.R demo/
all\_putative\_mosaics\_feature.list.features\_addphasing\_addvalidation demo/prediction\_model\_phasingcorrection.rds demo/test\_phasingcorrected


#### 2nd step: train based on corrected phasable sites:
```
Rscript train_on_correctedphasing.R
```
**Usage:** Rscript training\_on\_phasing.R trainset prediction\_model

**Demo:**
Rscript train\_on\_correctedphasing.R demo/test\_phasingcorrected demo/test2.rds


