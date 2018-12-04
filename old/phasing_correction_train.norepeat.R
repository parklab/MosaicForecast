#!/usr/bin/env Rscript
.libPaths( c( .libPaths(), "/n/data1/hms/dbmi/park/yanmei/tools/R_packages/") )

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
	stop("Rscript phasing_correction_train.norepeat.R trainset prediction_model_phasingcorrection output_file_phasingcorrected read_length(int)", call.=FALSE)
} else if (length(args)==4) {
	train_file <- args[1]
	prediction_model <- args[2]
	output_file <- args[3]
	read_length <- as.numeric(args[4])
}

library(stats)
library(caret)
library(nnet)
library(glmnet)
library(ggbiplot)

#head train_phasable_sites
#chr     pos     ref     alt     MAF     id      dp_p    querypos_p      leftpos_p       seqpos_p        mapq_p  baseq_p baseq_t ref_baseq1b_p   ref_baseq1b_t   alt_baseq1b_p   alt_baseq1b_t   sb_p    context     major_mismatches_mean   minor_mismatches_mean   mismatches_p    AF      dp      mosaic_likelihood       het_likelihood  refhom_likelihood       althom_likelihood       mapq_difference sb_read12_pdp_diff  repeats validation      phase conflicting_reads       phase_corrected
#10      10009041        G       T       0       fb0c6353-a90c-45e2-9355-7cd16cf756ff_10_10009041_G_T    0.243532076735  0.52723 0.51954 0.04323 0.31514 0.15322 -1.60007        0.93261 -0.40356        0.49078     0.36508 1       TCA     0.855   1.923   0.04837 0.339   115     0.954322027344281       0.0456779726557187      2.52244597973152e-116   2.97816045946957e-244   -1.57895        1       -4.89285700000001   rmsk    TP      hap=3   0       mosaic

input <- read.delim(train_file,header=TRUE)
all_phasable <- subset(input, phase != "notphased")
all_phasable <-all_phasable[!is.na(all_phasable$mosaic_likelihood),]
all_phasable$mapq_p[is.na(all_phasable$mapq_p)]<-1
all_phasable <- all_phasable[complete.cases(all_phasable[,seq(1,28)]),]

all_phasable.2 <- subset(all_phasable, select=c(querypos_p,leftpos_p, seqpos_p, mapq_p, baseq_p, baseq_t, ref_baseq1b_p, ref_baseq1b_t, alt_baseq1b_p, alt_baseq1b_t, sb_p, major_mismatches_mean, minor_mismatches_mean, mismatches_p, AF, dp, mosaic_likelihood, het_likelihood, refhom_likelihood, mapq_difference, sb_read12_p, dp_diff, conflict_num))
all_phasable.3<- all_phasable.2
all_phasable.3$querypos_p=log(all_phasable.3$querypos_p+1e-7)
all_phasable.3$leftpos_p=log(all_phasable.3$leftpos_p+1e-7)
all_phasable.3$seqpos_p=log(all_phasable.3$seqpos_p+1e-7)
all_phasable.3$mapq_p=log(all_phasable.3$mapq_p+1e-7)
all_phasable.3$baseq_p=log(all_phasable.3$baseq_p+1e-7)
all_phasable.3$ref_baseq1b_p=log(all_phasable.3$ref_baseq1b_p+1e-7)
all_phasable.3$alt_baseq1b_p=log(all_phasable.3$alt_baseq1b_p+1e-7)
all_phasable.3$sb_p=log(all_phasable.3$sb_p+1e-7)
all_phasable.3$mismatches_p=log(all_phasable.3$mismatches_p+1e-7)
all_phasable.3$sb_read12_p=log(all_phasable.3$sb_read12_p+1e-7)
all_phasable.3$major_mismatches_mean=all_phasable.3$major_mismatches_mean*read_length
all_phasable.3$minor_mismatches_mean=all_phasable.3$minor_mismatches_mean*read_length
##all_phasable.3$dp=log(all_phasable.3$dp)

pc<-prcomp(all_phasable.3,
       center = TRUE,
       scale. = TRUE) 

pdf("phasablesites_PCA.pdf")
g <- ggbiplot(pc, obs.scale = 1, var.scale = 1, 
	groups = all_phasable$phase, ellipse = TRUE, 
	circle = TRUE)+
	theme(legend.direction ='horizontal', 
	legend.position = 'top')
g
dev.off()

all_phasable$pc1 <- pc$x[,1]
all_phasable$pc2 <- pc$x[,2]
all_phasable$pc3 <- pc$x[,3]
all_phasable$pc4 <- pc$x[,4]
#all_phasable$pc5 <- pc$x[,5]

#dp_p    querypos_p      leftpos_p       seqpos_p        mapq_p  baseq_p  ref_baseq1b_p   alt_baseq1b_p  sb_p   mismatches_p    sb_read12_p     

set.seed(123)
all_train <- all_phasable[!is.na(all_phasable$validation),]
#all_train.2 <- subset(all_train, select=c(phase, validation, pc1, pc2, pc3, pc4, pc5))
all_train.2 <- subset(all_train, select=c(phase, validation, pc1, pc2, pc3, pc4))
all_train.2 <- subset(all_train.2, phase!="hap=2")
all_train.2$phase <- as.factor(all_train.2$phase)

model <- train(validation ~ ., all_train.2, method="glmnet",tuneGrid=expand.grid(.alpha=0:1, .lambda=0:30/10))
saveRDS(model,prediction_model)


all_phasable.4 <- subset(all_phasable, select=c(phase, validation, pc1, pc2, pc3, pc4))
colnames(all_phasable.4) <- c("phase","validation","pc1","pc2","pc3","pc4")
all_phasable_nonhet <- subset(all_phasable.4, phase!="hap=2")
all_phasable_nonhet$phase <- as.factor(all_phasable_nonhet$phase)

nonhet_phasable <- subset(all_phasable, phase!="hap=2")
het_phasable <- subset(all_phasable, phase=="hap=2")
nonhet_phasable$phase_model_corrected <- predict(model, all_phasable_nonhet)
het_phasable$phase_model_corrected <- "het"
all_phasable <- rbind(het_phasable, nonhet_phasable)

write.table(all_phasable, output_file,sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)


