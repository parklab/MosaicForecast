#!/usr/bin/env Rscript

.libPaths( c( .libPaths(), "/n/data1/hms/dbmi/park/yanmei/tools/R_packages/") )


args<-commandArgs(TRUE)
if (length(args)!=3) {
  stop("Rscript Prediction.R input_file(feature_list) model_trained output_file(predictions)", call.=FALSE)
} else if (length(args)==3) {
	input_file=args[1]
	model=args[2]
	output_file=args[3]
}
verbose=TRUE

library(caret)
library(e1071)

M <- readRDS(model)
input<-read.delim(input_file,header=TRUE,sep="\t")
input <- input[!is.na(input$mosaic_likelihood),]
input <- input[!is.na(input$baseq_t),]
input$mapq_p[is.na(input$mapq_p)]<-1
input[is.na(input)]<-0
#if(sum(input$MAF==".")>0){
#	input[input$MAF==".",]$MAF <-0
#}
#input$MAF <- as.character(input$MAF)
#input$MAF <- as.numeric(input$MAF)


output<- input
output$prediction <- predict(M,input)
prediction_probs <- predict(M,input,type="prob")
output <- cbind(output, prediction_probs)

write.table(output,row.names=FALSE, col.names=TRUE,sep="\t",file=output_file,quote=FALSE)



