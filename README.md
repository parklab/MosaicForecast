# MosaicForecast
A machine learning method that leverages read-based phasing and read-level features to accurately detect mosaic SNVs (SNPs, small indels) from NGS data.

![MF_pipeline](https://user-images.githubusercontent.com/8002850/55032948-61016880-4fe8-11e9-8cf8-343fd3cdd26e.png)


## Dependency:
### Required Interpreter Versions:
* Python version 3.4+
* R version 3.4+
### Python packages:
* collections
* itertools
* subprocess
* multiprocessing
* regex
* uuid
* math 
* numpy (1.16.1)
* pandas (0.20.1)
* pyfaidx (0.5.3)
* pysam (0.15.2)
* pysamstats (1.1.2)
* scipy (1.2.1)
### R packages:
* caret (6.0-80)
* e1071 (1.7-1)
* glmnet (2.0-16)
* nnet (7.3-12)
* mlr (2.13)
* RColorBrewer (1.1.2)
* ggplot2 ()
### Other softwares:
* samtools (1.9): https://sourceforge.net/projects/samtools/files/samtools/1.9/
* wigToBigWig (v4): http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
* bigWigAverageOverBed (v2): http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
### Install Dependencies:
1. We have created a docker image with all dependencies installed:
	https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh  
2. You could also install conda first, and then install the dependencies as described in the Dockerfile.
	https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh  

## Resources:
#### Mappability score: 
* Umap score (k=24): https://bismap.hoffmanlab.org/
#### Regions to filter out:
* Segmental Duplication regions: http://humanparalogy.gs.washington.edu/ 
* Regions enriched for SNPs with >=3 haplotypes: https://github.com/parklab/MosaicForecast/tree/master/resources
* Simple repeats (should be removed before calling mosaic INDELS): https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=716260999_jBJij5JXiQiuykIobdBExCLj0XEf
#### Population allale frequency
* Gnomad datasets: https://gnomad.broadinstitute.org/downloads

# Usage:
## Phasing:
**Usage:** 

python Phase.py bam_dir output_dir ref_fasta n_jobs_parallel input_positions min_dp_inforSNPs

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

python ReadLevel_Features_extraction.py input.bed output_features bam_dir ref.fa Umap_mappability(bigWig file,k=24) read_length n_jobs_parallel

**Note:**
1. Names of bam files should be "sample.bam" under the bam_dir.
2. There should be a fai file under the same dir of the fasta file (samtools faidx input.fa)
3. File format of the input.bed: chr pos-1 pos ref alt sample, sep=\t 
4. We did not use gnomad population AF as an feature (instead we use it to filter), but you can use it to train your model if you have interest in common variants
5. The program to extract mappability score: "bigWigAverageOverBed" could be downloaded here at http://hgdownload.soe.ucsc.edu/admin/exe/, the program to convert wiggle file to BigWig file "wigToBigWig", and the "fetchChromSizes" script to create the chrom.sizes file for the UCSC database with which you are working (e.g., hg19) could be downloaded from the same directory. The wiggle file containing mappability score (Umap,k=24) could be downloaded here: https://bismap.hoffmanlab.org/

**Demo:**
```
python ReadLevel_Features_extraction.py demo/test.input demo/test.features_forR demo ${ref.fa} ${k24.umap.wg.bw} 150 2  
```
**Output:**
```
A list of read-level features for each input site.
```

| id | dp_p | conflict_num | mappability | type | length | GCcontent | ref_softclip | alt_softclip | querypos_p | leftpos_p | seqpos_p | mapq_p | baseq_p | baseq_t | ref_baseq1b_p | ref_baseq1b_t | alt_baseq1b_p | alt_baseq1b_t | sb_p | context | major_mismatches_mean | minor_mismatches_mean | mismatches_p | AF | dp | mosaic_likelihood | het_likelihood | refhom_likelihood | althom_likelihood | mapq_difference | sb_read12_p | dp_diff |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| test\~11\~40316580\~C\~T | 0.3183008162818 | 0 | 1 | SNP | 0 | 0.476190476190476 | 0.0240384615384615 | 0 | 0.1582 | 0.16521 | 0.68821 | NA | 0.91657 | -0.57364 | 0.98911 | 0.21893 | 0.67576 | -0.8528 | 0.69934 | GGA | 0.00878 | 0.0144466666666667 | 0.29396 | 0.028 | 214 | 0.999414559235067 | 3.90999117967593e-49 | 0.000585440764932926 | 0 | 0 | 0.69142 | -12.8571 |
| test\~12\~52644508\~C\~T | 0.197545792452075 | 0 | 1 | SNP | 0 | 0.571428571428571 | 0.0208333333333333 | 0 | 0.19325 | 0.20057 | 0.88251 | NA | 0.11764 | -0.95448 | 0.31536 | 0.6827 | 0.31601 | 0.58756 | 0.13401 | CGC | 0.01236 | 0.0127266666666667 | 0.17424 | 0.054 | 203 | 0.999999999985368 | 5.11687178601205e-39 | 1.46319954019795e-11 | 0 | 0 | 0.36124 | -12.8571 |

```
1. id: uniq ID of the input candidate sites.
2. mappability: UMAP mappability score at the candidate site (k=24).
3. type: type of the candidate mutation (SNP, MNP, INS or DEL).
4. length: difference of base pair lengh of ref and alt allele for candidate sites.
5. GCcontent: 20-bp local GCcontent.
6. ref_softclip: proportion of soft-clipped reads for ref reads.
7. alt_softclip: proportion of soft-clipped reads for alt reads.
8. querypos_p: p-value or effect size by wilcoxon's rank sum test of base query positions of ref and alt alleles.
9. leftpos_p: p-value or effect size by wilcoxon's rank sum test of left-most positions of ref and alt reads.
10. seqpos_p: p-value or effect size by wilcoxon's rank sum test of base sequencing cycles of ref and alt alleles.
11. baseq_p: p-value or effect size by Wilcoxon's rank sum test of base qualities of ref and alt alleles.
12. baseq_t: The test statistic under the large-sample approximation that the rank sum statistic is normally distributed (wilcox rank sum test of base qualites of alt alleles vs. ref alleles).
13. ref_baseq1b_p: p-value or effect size by Wilcoxon's rank sum test of base qualities from ref reads at mutant position, compared with base qualities from ref reads at 1bp downtream of the mutant position.
14. ref_baseq1b_t: The test statistic under the large-sample approximation that the rank sum statistic is normally distributed (wilcox rank sum test of base qualities from ref reads at mutant position, compared with base qualities from ref reads at 1bp downtream of the mutant position).
15. alt_baseq1b_p: p-value or effect size by Wilcoxon's rank sum test of base qualities from alt reads at mutant position, compared with base qualities from alt reads at 1bp downtream of the mutant position.
16. alt_baseq1b_t: The test statistic under the large-sample approximation that the rank sum statistic is normally distributed (wilcox rank sum test of base qualities from alt reads at mutant position, compared with base qualities from alt reads at 1bp downtream of the mutant position).
17. context: three-nucleotide base context on the reads surrounding the mutant position.
18. major_mismatches_mean: average mismatches per ref reads.
19. minor_mismatches_mean: average mismatches per alt reads.
20. mismatches_p: p-value or effect size by Wilcoxon's rank sum test of mismatches per ref reads vs. mismatches per alt reads.
21. sb_p: p-value or effect size by Fisher's exact test of strand bias for ref and alt alleles.
22. sb_read12_p: p-value or effect size by Fisher's exact test of read1/read2 bias for ref and alt alleles.
23. mosaic_likelihood: mosaic genotype likelihood calculated (assuming uniform distribution of mosaics allele fraction from 0-1).
24. het_likelihood: Genotype likelihood of the variant being germline heterozygous.
25. refhom_likelihood: reference-homozygous genotype likelihood.
26. mapq_p: p-value or effect size by Wilcoxon's rank sum test of mapping qualities of ref and alt reads.
27. mapq_difference: difference of average map quality per alt reads vs. average map quality per ref reads.
28. AF: variant allele fraction.
29. dp: read depth at mutant position.
30. dp_diff: difference of average read depths of local (<200bp) and distant (>2kb) regions.
31. dp_p: p-value or effect size by Wilcoxon's rank sum test of read depths sampled within 200bp window surrounding the mutant position vs. read depths sampled in distant regions from the mutant position (>2kb).
32. conflict_num: number of read pairs supporting both ref and alt alleles.
```

## Genotype Prediction:

**Usage:**

Rscript Prediction.R input\_file(feature\_list) model\_trained output\_file(predictions)

**Note:**
1. The "input\_file" is a list of read-level features obtained in the last step.
2. The "model\_trained" is the pre-trained RF model to predict genotypes.

> You may use our models trained with brain WGS data for SNPs (paired-end read at 50-250X read depths, we train our models based on Mutect2-PON callings. To our experience, the models were pretty robust across different depths, but the best strategy would be using a model with similar depth with your data):
>
> * models\_trained/50xRFmodel\_addRMSK\_Refine.rds
> * models\_trained/100xRFmodel\_addRMSK\_Refine.rds
> * models\_trained/150xRFmodel\_addRMSK\_Refine.rds
> * models\_trained/200xRFmodel\_addRMSK\_Refine.rds
> * models\_trained/250xRFmodel\_addRMSK\_Refine.rds
>
> We also pre-trained a model for mosaic deletions (using paired-end read at 250X, with phasing information):
> * models\_trained/deletions\_250x.RF.rds

**Demo:**
```
Rscript Prediction.R demo/test.SNP.features models\_trained/250xRFmodel\_addRMSK\_Refine.rds  demo/test.SNP.predictions   
Rscript Prediction.R demo/test.DEL.features models\_trained/250xRFmodel\_addRMSK\_Refine.rds  demo/test.DEL.predictions
```
**Output:**
```
Genotype predictions for all input sites.
```

| id | AF | dp | prediction | het | mosaic | refhom | repeat |
| --- | --- | --- | --- | --- | --- | --- | --- |
| test\~11\~40316580\~C\~T | 0.028 | 214 | mosaic | 0.002 | 0.958 | 0 | 0.04 |
| test\~12\~52644508\~C\~T | 0.054 | 203 | mosaic | 0.002 | 0.982 | 0 | 0.016 |
| test\~15\~75918044\~G\~A | 0.036 | 193 | mosaic | 0.006 | 0.812 | 0 | 0.182 |
| test\~1\~1004865\~G\~C | 0.085 | 212 | mosaic | 0.006 | 0.988 | 0 | 0.006 |

```
1. prediction: genotype predictions including refhom, het, mosaic and repeat.
2. het/mosaic/refhom/repeat: genotyping probabilities for each genotype.
```

## You could also train RF models using your own data:
**Usage:**

Rscript Train_RFmodel.R input(trainset) output(prediction_model) type_model(Phase|Refine) type_variant(SNP|INS|DEL)

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

#### (Recommended when you have >=100 hap>=3 sites validated orthogonally or checked manually)
**Usage:**

Rscript Phasing_Refine_Multinomial_Logistic_Regression.R input(trainset) output1(model) output2(converted genotypes) read_length(int) type(pvalue|effectsize)

**Note:**

1. Use "pvalue" when your data has relatively even read coverage (i.e. WGS data) or the training sample size is big (i.e., >5000 sites);
2. Use "effectsize" when your data has extrmely un-even read coverage and small training sample size.
3. The input file should be a list of pre-generated read-level features for all phasable sites, adding a column termed "phase", containing the pre-generated haplotype number for each site (hap=2, hap=3, hap>3), and a column termed "validation", containing the orthogonally validation results. The un-evalulated sites shoule be "NA" in the "validation" column.
4. The output1 is the multinomial regression model, the output2 is the extraplolated four-category genotypes for all phasable sites.

**Demo:**
```
Rscript PhasingRefine.R demo/phasable_trainset demo/model_phasingcorrection.rds demo/phasable_sites_convertedgenotypes 150 demo/phasable_sites_Refine.pdf 
```
**Output:**
```
A list of extrapolated genotypes based on Phasing, Readlevel features and orthogonal validations.
```

| id | phase | validation | phase_model_corrected | pc1 | pc2 | pc3 | pc4 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1465\~7\~64358306\~C\~T | hap=3 | repeat | repeat | 1.76125132511822 | -0.0974980761360579 | 0.040830632773886 | 1.76651681595286 |
| 1465\~2\~10704065\~T\~C | hap=3 | repeat | repeat | -0.0124184653486693 | -0.289637541460141 | -2.01395435019693 | 1.5135587692184 |
| 1465\~10\~42529522\~G\~C | hap=3 | NA | repeat | 0.707159292549739 | -3.91325643368487 | 1.54152147288395 | -1.10846627458624 |
| 1465\~7\~122592074\~G\~T | hap=3 | mosaic | mosaic | -0.104370300644773 | 2.33641117703566 | -0.0299470543274239 | 1.40486521150486 |
| 1465\~X\~61712742\~A\~G | hap=3 | repeat | repeat | 0.574318366975694 | 1.16511088082416 | 1.24290479319458 | 1.24880945486403 |
| 1465\~10\~42544320\~C\~T | hap=3 | NA | repeat | 0.544602308768669 | 1.64954352441225 | 0.28095361817584 | 0.744802179936821 |

```
1. phase_model_corrected: Four-category genotypes extrapolated based on phasing and read-level features.
2. pc1/pc2/pc3/pc4/pc5: the first five PCA components constructed with read-level features.
3. demo/phasable_sites_Refine.pdf: A plot showing the genotype extrapolation from phasing to 4-category genotypes.
```
![genotype_extrapolation_phase](https://user-images.githubusercontent.com/8002850/55031935-239bdb80-4fe6-11e9-9a2f-76aa869d33ca.png)


## Contact:
If you have any questions please contact us:

Yanmei Dou: yanmei_dou@hms.harvard.edu, douyanmei@gmail.com  
Peter J Park: peter_park@hms.harvard.edu


