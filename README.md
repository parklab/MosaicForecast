# MosaicForecast
A machine learning method that leverages read-based phasing and read-level features to accurately detect mosaic SNVs (SNPs, small indels) from NGS data.

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
* Umap score (k=24): https://bismap.hoffmanlab.org/
#### Regions to filter out:
* Segmental Duplication regions: http://humanparalogy.gs.washington.edu/ 
* Regions enriched for SNPs with >=3 haplotypes: https://github.com/parklab/MosaicForecast/tree/master/resources
* Simple repeats (shoule be removed if call mosaic INDELS): https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=716260999_jBJij5JXiQiuykIobdBExCLj0XEf
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
python Phase.py demo demo/phasing ${human_g1k_v37_decoy.fasta} 2 demo/test.input 20
```

**Output:**
```
output_dir/all.phasing
```

| sample | chr | pos | ref | alt | phasing | conflicting_reads |
| --- | --- | --- | --- | --- | --- | --- |
| test | 12 | 52644508 | C | T | hap=3 | 0 |
| test | 15 | 75918044 | G | A | hap=3 | 0 |
| test | 1 | 1004865 | G | C | hap=3 | 0 |
| test | 1 | 2591769 | AG | A | hap>3 | 2 |
| test | 1 | 33801576 | TTTGTTG | T | hap=3 | 0 |

```
hap=2: likely het variants
hap=3: likely mosaic variants
hap>3: likely cnv/repeat
conflicting_reads: number of read pairs supporting both ref and alt alleles.

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
python ReadLevel_Features_extraction.py demo/test.input demo/test.features_forR demo ${ref.fa}
${k24.umap.wg.bw} 2
Rscript ReadLevel_Features_extraction.R demo/test.features demo/test.features 150 pvalue
 
```

**Output:**
```
A list of read-level features for each input site.
```

| id | dp_p | conflict_num | mappability | type | length | GCcontent | ref_softclip | alt_softclip | querypos_p | leftpos_p | seqpos_p | mapq_p | baseq_p | baseq_t | ref_baseq1b_p | ref_baseq1b_t | alt_baseq1b_p | alt_baseq1b_t | sb_p | context | major_mismatches_mean | minor_mismatches_mean | mismatches_p | AF | dp | mosaic_likelihood | het_likelihood | refhom_likelihood | althom_likelihood | mapq_difference | sb_read12_p | dp_diff |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| test\~11\~40316580\~C\~T | 0.3183008162818 | 0 | 1 | SNP | 0 | 0.476190476190476 | 0.0240384615384615 | 0 | 0.1582 | 0.16521 | 0.68821 | NA | 0.91657 | -0.57364 | 0.98911 | 0.21893 | 0.67576 | -0.8528 | 0.69934 | GGA | 0.00878 | 0.0144466666666667 | 0.29396 | 0.028 | 214 | 0.999414559235067 | 3.90999117967593e-49 | 0.000585440764932926 | 0 | 0 | 0.69142 | -12.8571 |
| test\~12\~52644508\~C\~T | 0.197545792452075 | 0 | 1 | SNP | 0 | 0.571428571428571 | 0.0208333333333333 | 0 | 0.19325 | 0.20057 | 0.88251 | NA | 0.11764 | -0.95448 | 0.31536 | 0.6827 | 0.31601 | 0.58756 | 0.13401 | CGC | 0.01236 | 0.0127266666666667 | 0.17424 | 0.054 | 203 | 0.999999999985368 | 5.11687178601205e-39 | 1.46319954019795e-11 | 0 | 0 | 0.36124 | -12.8571 |


## Genotype Prediction:

**Usage:**

Rscript Prediction.R input\_file(feature\_list) model\_trained output\_file(predictions)

**Note:**
1. The input\_file is a list of read-level features.
2. The model\_trained is the pre-trained RF model to predict genotypes.

**Demo:**

Rscript Prediction.R demo/test.SNP.features models\_trained/250xRFmodel\_addRMSK\_Refine.rds  demo/test.SNP.predictions
Rscript Prediction.R demo/test.DEL.features models\_trained/250xRFmodel\_addRMSK\_Refine.rds  demo/test.DEL.predictions

**Output:**
```
Genotype predictions for all input sites.
```
> You may use our models trained with brain WGS data (paired-end read at 50-250X read depths, we train our models based on Mutect2-PON callings):
>
> * models\_trained/50xRFmodel\_addRMSK\_Refine.rds
> * models\_trained/100xRFmodel\_addRMSK\_Refine.rds
> * models\_trained/150xRFmodel\_addRMSK\_Refine.rds
> * models\_trained/200xRFmodel\_addRMSK\_Refine.rds
> * models\_trained/250xRFmodel\_addRMSK\_Refine.rds


| id | prediction | het | mosaic | refhom | repeat |
| --- | --- | --- | --- | --- | --- |
| test\~11\~40316580\~C\~T | mosaic | 0.002 | 0.958 | 0 | 0.04 |
| test\~12\~52644508\~C\~T | mosaic | 0.002 | 0.982 | 0 | 0.016 |
| test\~15\~75918044\~G\~A | mosaic | 0.006 | 0.812 | 0 | 0.182 |
| test\~1\~1004865\~G\~C | mosaic | 0.006 | 0.988 | 0 | 0.006 |


## You could also train RF models using your own data:
**Usage:**
```
Rscript Train_RFmodel.R input(trainset) output(prediction_model) type_model(Phase|Refine) type_variant(SNP|INS|DEL)
```
**Note:** 

1. You could choose to train your model based on Phasing (hap=2, hap=3, hap>3, type in "Phase") or Refined genotypes ("mosaic","het","refhom","repeat", type in "Refine").
2. The input file should be a list of pre-generated read-level features, adding a column termed "phase" (Phase model) or "phase\_model\_corrected" (Refined genotypes model). 


**Demo:**
```
Rscript Train_RFmodel.R demo/phasable_trainset demo/Phase_model.rds Phase SNP
Rscript Train_RFmodel.R demo/phasable_trainset demo/Refine_model.rds Refine DEL
Rscript Train_RFmodel.R demo/deletions_phasable_trainset demo/Deletions_Refine_model.rds Refine DEL
```

**Output:**
```
Random Forest prediction model
```

## Convert phasing to four-category genotypes based on experimental data:

#### (Recommended when you have >=100 hap>=3 sites validated orthogonally)
**Usage:**
```
Rscript Phasing_Refine_Multinomial_Logistic_Regression.R input(trainset) output1(model)
output2(converted genotypes) read_length(int) type(pvalue|effectsize)
```

**Note:**

1. Use "pvalue" when your data has relatively even read coverage (i.e. WGS data) or the training sample size is big (i.e., >5000 sites);
2. Use "effectsize" when your data has extrmely un-even read coverage and small training sample size.
3. The input file should be a list of pre-generated read-level features for all phasable sites, adding a column termed "phase", containing the pre-generated haplotype number for each site (hap=2, hap=3, hap>3), and a column termed "validation", containing the orthogonally validation results. The un-evalulated sites shoule be "NA" in the "validation" column.
4. The output1 is the multinomial regression model, the output2 is the extraplolated four-category genotypes for all phasable sites.

**Demo:**
```
Rscript PhasingRefine_MultinomialLogisticRegression.R demo/phasable_trainset
demo/model_phasingcorrection.rds demo/phasable_sites_convertedgenotypes 150 pvalue
```

**Output:**
```
1. Four-category genotypes extrapolated based on phasing and read-level features
2. "phasablesites_PCA.pdf", showing positions of different phasable sites in the PCA space constructed with read-level features. 
```


