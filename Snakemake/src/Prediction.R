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

autosomal <- c(seq(1,22),paste0("chr",seq(1,22)))
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
output_autosomal<- subset(output, chromosome=="autosomal")
output_X<- subset(output, chromosome=="X")
output_Y<- subset(output, chromosome=="Y")

if(model_type=="Refine"){
df1 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)

df2 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)

df3 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)

df4 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)

df5 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)

df6 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)

	predictedothers <- subset(output, prediction !="mosaic")

	if(nrow(output_autosomal)>=1){
		predictedmosaics_autosomal <- subset(output, prediction=="mosaic" & chromosome=="autosomal")
		predictedmosaics_autosomal <- subset(predictedmosaics_autosomal, ref_softclip<0.1)
		predictedmosaics_autosomal$prediction <- as.character(predictedmosaics_autosomal$prediction)
		##optional filters (i.e., ultra-high depth):
		if (nrow(predictedmosaics_autosomal[predictedmosaics_autosomal$AF<0.01,])>=1){
			predictedmosaics_autosomal[predictedmosaics_autosomal$AF<0.01,]$prediction <- paste0(subset(predictedmosaics_autosomal,AF<0.01)$prediction,";cautious:AF<0.01")
		}
		if (nrow(predictedmosaics_autosomal[predictedmosaics_autosomal$AF*predictedmosaics_autosomal$dp<2,])>=1){
			predictedmosaics_autosomal[predictedmosaics_autosomal$AF*predictedmosaics_autosomal$dp<2,]$prediction <- paste0(subset(predictedmosaics_autosomal,AF*dp<2)$prediction,";cautious:only-1-altallele")
		}
		if (nrow(predictedmosaics_autosomal[predictedmosaics_autosomal$dp>=mean(output_autosomal$dp)*2,])>=1){
			predictedmosaics_autosomal[predictedmosaics_autosomal$dp>=mean(output_autosomal$dp)*2,]$prediction <- paste0(subset(predictedmosaics_autosomal,dp>=mean(output_autosomal$dp)*2)$prediction,";low-confidence:extra-high-coverage")
		}
		df1 <- subset(predictedmosaics_autosomal, !(dp >= mean(output_autosomal$dp)*1.5 & AF>=0.2))
		df2 <- subset(predictedmosaics_autosomal, (dp >= mean(output_autosomal$dp)*1.5 & AF>=0.2))
		if(nrow(df2)>=1){
			df2$prediction <- paste0(df2$prediction,";low-confidence:likelyCNV")
		}
	}

	if(nrow(output_X)>=1){
		predictedmosaics_X <- subset(output, prediction=="mosaic" & chromosome=="X")
		predictedmosaics_X <- subset(predictedmosaics_X, ref_softclip<0.1)
		predictedmosaics_X$prediction <- as.character(predictedmosaics_X$prediction)
		##optional filters (i.e., ultra-high depth):
		if (nrow(predictedmosaics_X[predictedmosaics_X$AF<0.01,])>=1){
			predictedmosaics_X[predictedmosaics_X$AF<0.01,]$prediction <- paste0(subset(predictedmosaics_X,AF<0.01)$prediction,";cautious:AF<0.01")
		}
		if (nrow(predictedmosaics_X[predictedmosaics_X$AF*predictedmosaics_X$dp<2,])>=1){
			predictedmosaics_X[predictedmosaics_X$AF*predictedmosaics_X$dp<2,]$prediction <- paste0(subset(predictedmosaics_X,AF*dp<2)$prediction,";cautious:only-1-altallele")
		}
		if (nrow(predictedmosaics_X[predictedmosaics_X$dp>=mean(output_X$dp)*2,])>=1){
			predictedmosaics_X[predictedmosaics_X$dp>=mean(output_X$dp)*2,]$prediction <- paste0(subset(predictedmosaics_X,dp>=mean(output_X$dp)*2)$prediction,";low-confidence:extra-high-coverage")
		}
		df3 <- subset(predictedmosaics_X, !(dp >= mean(output_X$dp)*1.5 & AF>=0.2))
		df4 <- subset(predictedmosaics_X, (dp >= mean(output_X$dp)*1.5 & AF>=0.2))
		if(nrow(df4)>=1){
			df4$prediction <- paste0(df4$prediction,";low-confidence:likelyCNV")
		}
	}


	if(nrow(output_Y)>=1){
		predictedmosaics_Y <- subset(output, prediction=="mosaic" & chromosome=="Y")
		predictedmosaics_Y <- subset(predictedmosaics_Y, ref_softclip<0.1)
		predictedmosaics_Y$prediction <- as.character(predictedmosaics_Y$prediction)
		##optional filters (i.e., ultra-high depth):
		if (nrow(predictedmosaics_Y[predictedmosaics_Y$AF<0.01,])>=1){
			predictedmosaics_Y[predictedmosaics_Y$AF<0.01,]$prediction <- paste0(subset(predictedmosaics_Y,AF<0.01)$prediction,";cautious:AF<0.01")
		}
		if (nrow(predictedmosaics_Y[predictedmosaics_Y$AF*predictedmosaics_Y$dp<2,])>=1){
			predictedmosaics_Y[predictedmosaics_Y$AF*predictedmosaics_Y$dp<2,]$prediction <- paste0(subset(predictedmosaics_Y,AF*dp<2)$prediction,";cautious:only-1-altallele")
		}
		if (nrow(predictedmosaics_Y[predictedmosaics_Y$dp>=mean(output_Y$dp)*2,])>=1){
			predictedmosaics_Y[predictedmosaics_Y$dp>=mean(output_Y$dp)*2,]$prediction <- paste0(subset(predictedmosaics_Y,dp>=mean(output_Y$dp)*2)$prediction,";low-confidence:extra-high-coverage")
		}
		df5 <- subset(predictedmosaics_Y, !(dp >= mean(output_Y$dp)*1.5 & AF>=0.2))
		df6 <- subset(predictedmosaics_Y, (dp >= mean(output_Y$dp)*1.5 & AF>=0.2))
		if(nrow(df6)>=1){
			df6$prediction <- paste0(df6$prediction,";low-confidence:likelyCNV")
		}
	}
	output <- rbind(df1,df2,df3,df4,df5,df6, predictedothers)
}else if(model_type=="Phase"){

df1 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)

df2 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)

df3 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)

df4 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)

df5 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)

df6 <- data.frame(Date=as.Date(character()),
                 File=character(),
                 User=character(),
                 stringsAsFactors=FALSE)


	predictedothers <- subset(output, prediction !="hap=3")


	if(nrow(output_X)>=1){
		predictedmosaics_X <- subset(output, prediction=="hap=3" & chromosome=="X")
		predictedmosaics_X <- subset(predictedmosaics_X, ref_softclip<0.1)
		predictedmosaics_X$prediction <- as.character(predictedmosaics_X$prediction)
		##optional filters (i.e., ultra-high depth):
		if (nrow(predictedmosaics_X[predictedmosaics_X$AF<0.01,])>=1){
			predictedmosaics_X[predictedmosaics_X$AF<0.01,]$prediction <- paste0(subset(predictedmosaics_X,AF<0.01)$prediction,";cautious:AF<0.01")
		}
		if (nrow(predictedmosaics_X[predictedmosaics_X$AF*predictedmosaics_X$dp<2,])>=1){
			predictedmosaics_X[predictedmosaics_X$AF*predictedmosaics_X$dp<2,]$prediction <- paste0(subset(predictedmosaics_X,AF*dp<2)$prediction,";cautious:only-1-altallele")
		}
		if (nrow(predictedmosaics_X[predictedmosaics_X$dp>=mean(output_X$dp)*2,])>=1){
			predictedmosaics_X[predictedmosaics_X$dp>=mean(output_X$dp)*2,]$prediction <- paste0(subset(predictedmosaics_X,dp>=mean(output_X$dp)*2)$prediction,";low-confidence:extra-high-coverage")
		}
		df1 <- subset(predictedmosaics_X, !(dp >= mean(output_X$dp)*1.5 & AF>=0.2))
		df2 <- subset(predictedmosaics_X, (dp >= mean(output_X$dp)*1.5 & AF>=0.2))
		if(nrow(df2)>=1){
			df2$prediction <- paste0(df2$prediction,";low-confidence:likelyCNV")
		}
	}

	if(nrow(output_Y)>=1){
		predictedmosaics_Y <- subset(output, prediction=="hap=3" & chromosome=="Y")
		predictedmosaics_Y <- subset(predictedmosaics_Y, ref_softclip<0.1)
		predictedmosaics_Y$prediction <- as.character(predictedmosaics_Y$prediction)
		##optional filters (i.e., ultra-high depth):
		if (nrow(predictedmosaics_Y[predictedmosaics_Y$AF<0.01,])>=1){
			predictedmosaics_Y[predictedmosaics_Y$AF<0.01,]$prediction <- paste0(subset(predictedmosaics_Y,AF<0.01)$prediction,";cautious:AF<0.01")
		}
		if (nrow(predictedmosaics_Y[predictedmosaics_Y$AF*predictedmosaics_Y$dp<2,])>=1){
			predictedmosaics_Y[predictedmosaics_Y$AF*predictedmosaics_Y$dp<2,]$prediction <- paste0(subset(predictedmosaics_Y,AF*dp<2)$prediction,";cautious:only-1-altallele")
		}
		if (nrow(predictedmosaics_Y[predictedmosaics_Y$dp>=mean(output_Y$dp)*2,])>=1){
			predictedmosaics_Y[predictedmosaics_Y$dp>=mean(output_Y$dp)*2,]$prediction <- paste0(subset(predictedmosaics_Y,dp>=mean(output_Y$dp)*2)$prediction,";low-confidence:extra-high-coverage")
		}
		df3 <- subset(predictedmosaics_Y, !(dp >= mean(output_Y$dp)*1.5 & AF>=0.2))
		df4 <- subset(predictedmosaics_Y, (dp >= mean(output_Y$dp)*1.5 & AF>=0.2))
		if(nrow(df4)>=1){
			df4$prediction <- paste0(df4$prediction,";low-confidence:likelyCNV")
		}
	}

	if(nrow(output_autosomal)>=1){
		predictedmosaics_autosomal <- subset(output, prediction=="hap=3" & chromosome=="autosomal")
		predictedmosaics_autosomal <- subset(predictedmosaics_autosomal, ref_softclip<0.1)
		predictedmosaics_autosomal$prediction <- as.character(predictedmosaics_autosomal$prediction)
		##optional filters (i.e., ultra-high depth):
		if (nrow(predictedmosaics_autosomal[predictedmosaics_autosomal$AF<0.01,])>=1){
			predictedmosaics_autosomal[predictedmosaics_autosomal$AF<0.01,]$prediction <- paste0(subset(predictedmosaics_autosomal,AF<0.01)$prediction,";cautious:AF<0.01")
		}
		if (nrow(predictedmosaics_autosomal[predictedmosaics_autosomal$AF*predictedmosaics_autosomal$dp<2,])>=1){
			predictedmosaics_autosomal[predictedmosaics_autosomal$AF*predictedmosaics_autosomal$dp<2,]$prediction <- paste0(subset(predictedmosaics_autosomal,AF*dp<2)$prediction,";cautious:only-1-altallele")
		}
		if (nrow(predictedmosaics_autosomal[predictedmosaics_autosomal$dp>=mean(output_autosomal$dp)*2,])>=1){
			predictedmosaics_autosomal[predictedmosaics_autosomal$dp>=mean(output_autosomal$dp)*2,]$prediction <- paste0(subset(predictedmosaics_autosomal,dp>=mean(output_autosomal$dp)*2)$prediction,";low-confidence:extra-high-coverage")
		}
		df5 <- subset(predictedmosaics_autosomal, !(dp >= mean(output_autosomal$dp)*1.5 & AF>=0.2))
		df6 <- subset(predictedmosaics_autosomal, (dp >= mean(output_autosomal$dp)*1.5 & AF>=0.2))
		if(nrow(df6)>=1){
			df6$prediction <- paste0(df6$prediction,";low-confidence:likelyCNV")
		}
	}
	output <- rbind(df1,df2,df3,df4,df5,df6, predictedothers)
}

write.table(output,row.names=FALSE, col.names=TRUE,sep="\t",file=output_file,quote=FALSE)



