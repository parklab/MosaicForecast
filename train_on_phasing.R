#!/usr/bin/env Rscript
.libPaths( c( .libPaths(), "/n/data1/hms/dbmi/park/yanmei/tools/R_packages/") )

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
	stop("Rscript training_on_phasing.R trainset(demo/all_putative_mosaics_feature.list.features_addphasing_addvalidation) prediction_model", call.=FALSE)
} else if (length(args)==2) {
	input_file <- args[1]
	prediction_model <- args[2]
}

library(caret)
library(e1071)

input <- read.delim(input_file, header=TRUE)
input$mapq_p[is.na(input$mapq_p)]<-1
input$phase <- gsub('hap>3', 'repeat', gsub('hap=3', 'mosaic', gsub('hap=2','het',input$phase)))

all_train <- subset(input, phase != "notphased")
all_train <-all_train[!is.na(all_train$mosaic_likelihood),]
all_train.2 <- subset(all_train, select=-c(althom_likelihood, id, context, validation,dp_p))

control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
tunegrid <- expand.grid(.mtry=30)
metric <- "Accuracy"
rf_gridsearch <- train(phase~., data=all_train.2, method="rf", metric=metric,tuneGrid=tunegrid, trControl=control)
saveRDS(rf_gridsearch,file=prediction_model)

#input$prediction_phasing <- predict(rf_gridsearch, input)
#write.table(input, "test.prediction",sep="\t",quote=FALSE,row.names=FALSE, col.names=TRUE)




