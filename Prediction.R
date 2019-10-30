#!/usr/bin/env Rscript

.libPaths( c( .libPaths(), "/n/data1/hms/dbmi/park/yanmei/tools/R_packages/") )


args<-commandArgs(TRUE)
if (length(args)!=4) {
  stop("Rscript Prediction.R input_file(feature_list) model_trained model_type(Phase|Refine) output_file(predictions)", call.=FALSE)
} else if (length(args)==4) {
	input_file=args[1]
	model=args[2]
	model_type=args[3]
	output_file=args[4]
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
input$sb_p[input$sb_p=="Inf"]<- 100
input$sb_read12_p[input$sb_read12_p=="Inf"]<- 100


output<- input
output$prediction <- predict(M,input)
prediction_probs <- predict(M,input,type="prob")
output <- cbind(output, prediction_probs)
output <- subset(output,mappability>0)

df1 <- subset(output,(type=="SNP" | type=="MNP") & indel_proportion_SNPonly<0.3)
df2 <- subset(output, (type!="SNP" & type!="MNP") )
output <- rbind(df1,df2)

autosomal <- c(seq(1,23),paste0("chr",seq(1,22)))
chrX <- c("chrX","X")
chrY <- c("chrY","Y")

my_classify <- function(x){
	chr=unlist(strsplit(as.character(x[1]),"~"))[2]
	##if(grepl("X",chr)){
	if(chr %in% chrX){
		return("X")
	}else if(chr %in% chrY){
		return("Y")
	}else if(chr %in% autosomal){
		return("autosomal")
	}else{
		return("others")
	}
}
output$chromosome <- apply(output,1,my_classify)

if(model_type=="Refine"){
	predictedothers <- subset(output, prediction !="mosaic")
	predictedmosaics <- subset(output, prediction=="mosaic")
	predictedmosaics <- subset(predictedmosaics, ref_softclip<0.1)
	predictedmosaics$prediction <- as.character(predictedmosaics$prediction)
	##optional filters (i.e., ultra-high depth):
	if (nrow(predictedmosaics[predictedmosaics$AF<0.01,])>=1){
		predictedmosaics[predictedmosaics$AF<0.01,]$prediction <- paste0(subset(predictedmosaics,AF<0.01)$prediction,";cautious:AF<0.01")
	}
	if (nrow(predictedmosaics[predictedmosaics$AF*predictedmosaics$dp<2,])>=1){
		predictedmosaics[predictedmosaics$AF*predictedmosaics$dp<2,]$prediction <- paste0(subset(predictedmosaics,AF*dp<2)$prediction,";cautious:only-1-altallele")
	}
	if (nrow(predictedmosaics[predictedmosaics$dp>=mean(output$dp)*2,])>=1){
		predictedmosaics[predictedmosaics$dp>=mean(output$dp)*2,]$prediction <- paste0(subset(predictedmosaics,dp>=mean(output$dp)*2)$prediction,";low-confidence:extra-high-coverage")
	}
	df1 <- subset(predictedmosaics, !(dp >= mean(output$dp)*1.5 & AF>=0.2))
	df2 <- subset(predictedmosaics, (dp >= mean(output$dp)*1.5 & AF>=0.2))
	if(nrow(df2)>=1){
		df2$prediction <- paste0(df2$prediction,";low-confidence:likelyCNV")
	}
	output <- rbind(df1,df2, predictedothers)
}else if(model_type=="Phase"){
	predictedothers <- subset(output, prediction !="hap=3")
	predictedmosaics <- subset(output, prediction=="hap=3")
	predictedmosaics <- subset(predictedmosaics, ref_softclip<0.1)
	predictedmosaics$prediction <- as.character(predictedmosaics$prediction)
	##optional filters (i.e., ultra-high depth):
	if (nrow(predictedmosaics[predictedmosaics$AF<0.01,])>=1){
		predictedmosaics[predictedmosaics$AF<0.01,]$prediction <- paste0(subset(predictedmosaics,AF<0.01)$prediction,";cautious:AF<0.01")
	}
	if (nrow(predictedmosaics[predictedmosaics$AF*predictedmosaics$dp<2,])>=1){
		predictedmosaics[predictedmosaics$AF*predictedmosaics$dp<2,]$prediction <- paste0(subset(predictedmosaics,AF*dp<2)$prediction,";cautious:only-1-altallele")
	}
	if (nrow(predictedmosaics[predictedmosaics$dp>=mean(output$dp)*2,])>=1){
		predictedmosaics[predictedmosaics$dp>=mean(output$dp)*2,]$prediction <- paste0(subset(predictedmosaics,dp>=mean(output$dp)*2)$prediction,";low-confidence:extra-high-coverage")
	}
	df1 <- subset(predictedmosaics, !(dp >= mean(output$dp)*1.5 & AF>=0.2))
	df2 <- subset(predictedmosaics, (dp >= mean(output$dp)*1.5 & AF>=0.2))
	if(nrow(df2)>=1){
		df2$prediction <- paste0(df2$prediction,";low-confidence:likelyCNV")
	}
	output <- rbind(df1,df2, predictedothers)
}

write.table(output,row.names=FALSE, col.names=TRUE,sep="\t",file=output_file,quote=FALSE)



