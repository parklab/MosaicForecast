# MosaicForecast
A machine learning method that leverages read-based phasing and read-level features to accurately detect mosaic single-nucleotide variants (SNVs) from NGS data.

## Required Interpreter Versions:
* Python version 3.4+
* R version 3.4+

## Packages need to be installed:
### Python packages:
* collections
* itertools
* subprocess
* multiprocessing
* regex 
* numpy (1.16.1)
* pandas (0.20.1)
* pyfaidx (0.5.3.1)
* pysam (0.11.2.2)
* pysamstats (1.0.1)
* scipy (1.2.1)
### R packages:
* caret (6.0-78)
* e1071 (1.6-8)
* glmnet (2.0-13)
* nnet (7.3-12)
* ggbiplot (0.55)
### Others:
* samtools (1.9)
* wigToBigWig (v4)
* bigWigAverageOverBed (v2)

## Resources:
#### Mappability score: 
* Umap: https://bismap.hoffmanlab.org/
#### Regions to filter out:
* Segmental Duplication regions: http://humanparalogy.gs.washington.edu/ 
* Regions enriched for SNPs with >=3 haplotypes: https://github.com/parklab/MosaicForecast/tree/master/resources
#### Population allale frequency
* Gnomad datasets: https://gnomad.broadinstitute.org/downloads


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

**Demo:**

```
python Phase.py demo demo/phasing ${reference_dir}/human_g1k_v37_decoy.fasta 2 demo/test.input 20
```

**Output:**
```
output_dir/all.phasing
hap=2: likely het variants
hap=3: likely mosaic variants
hap>3: likely cnv/repeat

Intermediate files:
1. output_dir/all.merged.inforSNPs.pos: all nearby inforSNPs of candidate mosaics.
2. output_dir/all_2x2table: 2x2 tables by all nearby inforSNPs.
3. output_dir/all.phasing_2by2: Phasing results of mosaics and all nearby inforSNPs (2x2 table).
4. output_dir/multiple_inforSNPs.log: Phasing results of different pairs of inforSNPs.

```

## Extraction of read-level features:
**Usage:**

### 1st step:
```
python ReadLevel_Features_extraction.py input.bed output_features bam_dir ref.fa
Umap_mappability(bigWig file,k=24) n_jobs_parallel

```
**Note:**
1. Names of bam files should be "sample.bam" under the bam_dir.
2. There should be a fai file under the same dir of the fasta file (samtools faidx input.fa)
3. File format of the input.bed: chr pos-1 pos ref alt sample, sep=\t 
4. We did not use gnomad population AF as an feature (instead we use it to filter), but you can use it to train your model if you have interest in common variants
5. The program to extract mappability score: "bigWigAverageOverBed" could be downloaded here at http://hgdownload.soe.ucsc.edu/admin/exe/, the program to convert wiggle file to BigWig file "wigToBigWig", and the "fetchChromSizes" script to create the chrom.sizes file for the UCSC database with which you are working (e.g., hg19) could be downloaded from the same directory. The wiggle file containing mappability score (Umap,k=24) could be downloaded here: https://bismap.hoffmanlab.org/

### 2nd step:
```
Rscript Rscript ReadLevel_Features_extraction.R input_file output_file read_length(integer) 
type(pvalue||effectsize)
```
**Note:**
1. Use "pvalue" when your data has relatively even read coverage (i.e. WGS data) or the training sample size is big (i.e., >10000 sites);
2. Use "effectsize" when your data has extrmely un-even read coverage and small training sample size. The "effectsize" mode is relatively slow.
3. The input\_file is the file from the 1st step.

**Demo:**
```
python ReadLevel_Features_extraction.py demo/test.input demo/test.features demo ${ref.fa} ${Umap.bw} 2
Rscript ReadLevel_Features_extraction.R demo/test.features demo/test.features_R 150 pvalue
 
```

**Output:**
```
A list of read-level features for each input site.
```

## Genotype Prediction:

**Usage:**

Rscript Prediction.R input\_file(feature\_list) model\_trained output\_file(predictions)

**Note:**
1. The input\_file is a list of read-level features.
2. The model\_trained is the pre-trained RF model to predict genotypes.

**Demo:**

Rscript Prediction.R demo/test.features\_R models\_trained/brain\_MT2-PON.250x\_MosaicForecast-Refine\_pvalue.rds  demo/test\_predictions

**Output:**
```
Genotype predictions for all input sites.
```


> You may use our models trained with brain WGS data (paired-end read at 50-250X read depths, we train our models based on Mutect2-PON callings):
>
> * models\_trained/50-200X\_brainWGS/brain\_MT2-PON.50x.rds
> * models\_trained/50-200X\_brainWGS/brain\_MT2-PON.100x.rds
> * models\_trained/50-200X\_brainWGS/brain\_MT2-PON.150x.rds
> * models\_trained/50-200X\_brainWGS/brain\_MT2-PON.200x.rds
> * models\_trained/brain\_MT2-PON.250x\_MosaicForecast-Refine\_pvalue.rds
> * models\_trained/brain\_MT2-PON.250x\_MosaicForecast-Phase\_pvalue.rds


## You could also train RF models using your own data:
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

**Output:**
```
Random Forest prediction model
```

## Convert phasing to four-category genotypes based on experimental data:

#### (Recommended when you have >=100 hap>=3 sites validated orthogonally)
**Usage:**
```
Rscript Phasing_Refine_Multinomial_Logistic_Regression.R input(trainset) output1(model) output2(converted genotypes) read_length(int) type(pvalue|effectsize)
```

**Note:**

1. Use "pvalue" when your data has relatively even read coverage (i.e. WGS data) or the training sample size is big (i.e., >10000 sites);
2. Use "effectsize" when your data has extrmely un-even read coverage and small training sample size.
3. The input file should be a list of pre-generated read-level features for all phasable sites, adding a column termed "phase", containing the pre-generated haplotype number for each site (hap=2, hap=3, hap>3), and a column termed "validation", containing the orthogonally validation results. The un-evalulated sites shoule be "NA" in the "validation" column.
4. The output1 is the multinomial regression model, the output2 is the extraplolated four-category genotypes for all phasable sites.

**Demo:**
```
Rscript PhasingRefine_MultinomialLogisticRegression.R demo/phasable_trainset demo/model_phasingcorrection.rds demo/phasable_sites_convertedgenotypes 150 pvalue
```

**Output:**
```
1. Four-category genotypes extrapolated based on phasing and read-level features
2. "phasablesites_PCA.pdf", showing positions of different phasable sites in the PCA space constructed with read-level features. 
```


