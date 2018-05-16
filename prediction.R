#!/usr/bin/env Rscript

.libPaths( c( .libPaths(), "/n/data1/hms/dbmi/park/yanmei/tools/R_packages/") )


args<-commandArgs(TRUE)
if (length(args)!=3) {
  stop("Rscript prediction.R input_file(feature_list) model_trained output_file(predictions)", call.=FALSE)
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
output<- input
output$prediction <- predict(M,input)

write.table(output,row.names=FALSE, col.names=TRUE,sep="\t",file=output_file,quote=FALSE)



