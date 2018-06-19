#!/usr/bin/env Rscript
.libPaths( c( .libPaths(), "/n/data1/hms/dbmi/park/yanmei/tools/R_packages/") )

args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) {
	#stop("Rscript training_phasing_correction.R validated_phasble_sites_input output_model_phasingcorrection output_file_phasingcorrected", call.=FALSE)
	stop("Rscript training_phasing_correction.R trainset output_model_phasingcorrection output_file_phasingcorrected", call.=FALSE)
} else if (length(args)==3) {
	train_file <- args[1]
	output_model <- args[2]
	output_file <- args[3]
}

library(stats)
library(caret)
library(nnet)
library(glmnet)
library(ggbiplot)

input <- read.delim(train_file,header=TRUE)
all_phasable <- subset(input, phase != "notphased")
all_phasable <-all_phasable[!is.na(all_phasable$mosaic_likelihood),]
all_phasable$mapq_p[is.na(all_phasable$mapq_p)]<-1
#all_phasable$phase <- as.factor(gsub('hap>3', 'repeat', gsub('hap=3', 'mosaic', gsub('hap=2','het',all_phasable$phase))))
all_phasable$phase <- gsub('hap>3', 'repeat', gsub('hap=3', 'mosaic', gsub('hap=2','het',all_phasable$phase)))

all_phasable.2 <- subset(all_phasable, select=-c(althom_likelihood, id, context, phase, validation, dp_p))
#all_phasable.2 <- subset(all_phasable, select=-c(althom_likelihood, id, context, phase, validation))
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
##all_phasable.3$dp=log(all_phasable.3$dp)

pc<-prcomp(all_phasable.3,
       center = TRUE,
       scale. = TRUE) 

pdf("demo/phasablesites_PCA.pdf")
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

#dp_p    querypos_p      leftpos_p       seqpos_p        mapq_p  baseq_p  ref_baseq1b_p   alt_baseq1b_p  sb_p   mismatches_p    sb_read12_p     
all_train <- all_phasable[!is.na(all_phasable$validation),]
all_train.2 <- subset(all_train, select=c(phase, validation, pc1, pc2, pc3, pc4))
#all_train.2 <- subset(all_train, select=-c(althom_likelihood, id, context, phase, validation))

model <- train(validation ~ ., all_train.2, method="glmnet",tuneGrid=expand.grid(.alpha=0:1, .lambda=0:30/10))
saveRDS(model,output_model)


all_phasable.4 <- subset(all_phasable, select=c(phase, validation, pc1, pc2, pc3, pc4))
colnames(all_phasable.4) <- c("phase","validation","pc1","pc2","pc3","pc4")

all_phasable$phase_corrected <- predict(model, all_phasable.4)
write.table(all_phasable, output_file,sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)


