#!/usr/bin/env Rscript
.libPaths( c( .libPaths(), "/n/data1/hms/dbmi/park/yanmei/tools/R_packages/") )

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
	stop("Rscript Train_RFmodel.R trainset prediction_model type(Phase|Refine)
	Note:
	The \"Phase\" model indicates the RF model trained on phasing (hap=2, hap=3, hap>3); 
	The \"Refine\" model indicates the RF model trained on Refined-genotypes from the multinomial logistic regression model (het, mosaic, repeat, refhom)
", call.=FALSE)
} else if (length(args)==3) {
	input_file <- args[1]
	prediction_model <- args[2]
	type <- as.character(args[3])
}

library(caret)
library(e1071)

if (type=="Phase") {
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
} else if (type=="Refine"){
	input <- read.delim(input_file, header=TRUE)
	input$mapq_p[is.na(input$mapq_p)]<-1
	all_train <- input
	all_train <- subset(input, phase != "notphased")
	all_train <-all_train[!is.na(all_train$mosaic_likelihood),]
	#if(sum(all_train$MAF==".")>0){
	#        all_train$MAF<-0
	#}
	#all_train$MAF[is.na(all_train$MAF)]<-0
	all_train.2 <- subset(all_train, select=c(querypos_p,leftpos_p,seqpos_p,mapq_p,baseq_p,baseq_t,ref_baseq1b_p,ref_baseq1b_t,alt_baseq1b_p,alt_baseq1b_t,sb_p,context,major_mismatches_mean,minor_mismatches_mean,mismatches_p,AF,dp,mapq_difference,sb_read12_p,dp_diff,mosaic_likelihood,het_likelihood,refhom_likelihood,phase_model_corrected,conflict_num))
	#all_train.2 <- subset(all_train, select=c(querypos_p,leftpos_p,seqpos_p,mapq_p,baseq_p,baseq_t,ref_baseq1b_p,ref_baseq1b_t,alt_baseq1b_p,alt_baseq1b_t,sb_p,context,major_mismatches_mean,minor_mismatches_mean,mismatches_p,AF,dp,mapq_difference,sb_read12_p,dp_diff,mosaic_likelihood,het_likelihood,refhom_likelihood,phase_corrected,MAF,repeats,ECNT,HCNT))
	
	all_train.2$sb_p[all_train.2$sb_p=="Inf"]<- 100
	all_train.2$sb_read12_p[all_train.2$sb_read12_p=="Inf"]<- 100
	
	control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
	tunegrid <- expand.grid(.mtry=30)
	metric <- "Accuracy"
	rf_gridsearch <- train(phase_model_corrected ~., data=all_train.2, method="rf", metric=metric,tuneGrid=tunegrid, trControl=control,na.action=na.exclude)
	saveRDS(rf_gridsearch,file=prediction_model)

#write.table(input, "test.prediction",sep="\t",quote=FALSE,row.names=FALSE, col.names=TRUE)
}



