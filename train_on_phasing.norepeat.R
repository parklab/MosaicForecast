#!/usr/bin/env Rscript
.libPaths( c( .libPaths(), "/n/data1/hms/dbmi/park/yanmei/tools/R_packages/") )

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
	stop("Rscript training_on_phasing.R trainset prediction_model", call.=FALSE)
} else if (length(args)==2) {
	input_file <- args[1]
	prediction_model <- args[2]
}

library(caret)
library(e1071)

#head demo/test_phasingcorrected
#id      dp_p    querypos_p      leftpos_p       seqpos_p        mapq_p  baseq_p baseq_t ref_baseq1b_p   ref_baseq1b_t   alt_baseq1b_p   alt_baseq1b_t   sb_p    context major_mismatches_mean   minor_mismatches_mean       mismatches_p    AF      dp      mapq_difference sb_read12_p     dp_diff mosaic_likelihood       het_likelihood  refhom_likelihood       althom_likelihood       phase   validation      pc1pc2      pc3     pc4     phase_corrected

input <- read.delim(input_file, header=TRUE)
input$mapq_p[is.na(input$mapq_p)]<-1
all_train <- input
all_train <- subset(input, phase != "notphased")
all_train$phase <- as.factor(as.character(all_train$phase))
all_train <-all_train[!is.na(all_train$mosaic_likelihood),]
#all_train.2 <- subset(all_train, select=-c(althom_likelihood, id, validation, dp_p, pc1, pc2, pc3, pc4, phase))
#all_train.2 <- subset(all_train, select=c(querypos_p,leftpos_p,seqpos_p,mapq_p,baseq_p,baseq_t,ref_baseq1b_p,ref_baseq1b_t,alt_baseq1b_p,alt_baseq1b_t,sb_p,context,major_mismatches_mean,minor_mismatches_mean,mismatches_p,AF,dp,mapq_difference,sb_read12_p,dp_diff,mosaic_likelihood,het_likelihood,refhom_likelihood,phasing))
all_train.2 <- subset(all_train, select=c(querypos_p,leftpos_p,seqpos_p,mapq_p,baseq_p,baseq_t,ref_baseq1b_p,ref_baseq1b_t,alt_baseq1b_p,alt_baseq1b_t,sb_p,context,major_mismatches_mean,minor_mismatches_mean,mismatches_p,AF,dp,mapq_difference,sb_read12_p,dp_diff,mosaic_likelihood,het_likelihood,refhom_likelihood,phase))

control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
tunegrid <- expand.grid(.mtry=30)
metric <- "Accuracy"
rf_gridsearch <- train(phase ~., data=all_train.2, method="rf", metric=metric,tuneGrid=tunegrid, trControl=control,na.action=na.exclude)
saveRDS(rf_gridsearch,file=prediction_model)

#input$prediction_phasing <- predict(rf_gridsearch, input)
#write.table(input, "test.prediction",sep="\t",quote=FALSE,row.names=FALSE, col.names=TRUE)




