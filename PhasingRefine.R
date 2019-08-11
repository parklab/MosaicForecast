#!/usr/bin/env Rscript
.libPaths( c( .libPaths(), "/n/data1/hms/dbmi/park/yanmei/tools/R_packages/") )

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
	stop("Rscript PhasingRefine.R trainset prediction_model_phasingcorrection output_file_phasingcorrected read_length(int) pdf(PhasingRefine plot) variant_type(SNP|INDEL)
", call.=FALSE)
} else if (length(args)==6) {
	train_file <- args[1]
	prediction_model <- args[2]
	output_file <- args[3]
	read_length <- as.numeric(args[4])
	output_pdf <- args[5]
	variant_type <- args[6]
	#type <- as.character(args[5])
	#model_type <- as.character(args[6])
}

#library(ggbiplot)
library(stats)
library(caret)
library(nnet)
library(glmnet)
library(e1071)
library(ggplot2)
library(RColorBrewer)

model_type="glmnet"
type="pvalue"

my_chrXY <- function(x){
  !(strsplit(x,"~")[[1]][2]=="X"||strsplit(x,"~")[[1]][2]=="Y")
}

if (variant_type=="SNP") {
	#head train_phasable_sites
	#chr     pos     ref     alt     MAF     id      dp_p    querypos_p      leftpos_p       seqpos_p        mapq_p  baseq_p baseq_t ref_baseq1b_p   ref_baseq1b_t   alt_baseq1b_p   alt_baseq1b_t   sb_p    context     major_mismatches_mean   minor_mismatches_mean   mismatches_p    AF      dp      mosaic_likelihood       het_likelihood  refhom_likelihood       althom_likelihood       mapq_difference sb_read12_pdp_diff  repeats validation      phase conflicting_reads       phase_corrected
	#10      10009041        G       T       0       fb0c6353-a90c-45e2-9355-7cd16cf756ff_10_10009041_G_T    0.243532076735  0.52723 0.51954 0.04323 0.31514 0.15322 -1.60007        0.93261 -0.40356        0.49078     0.36508 1       TCA     0.855   1.923   0.04837 0.339   115     0.954322027344281       0.0456779726557187      2.52244597973152e-116   2.97816045946957e-244   -1.57895        1       -4.89285700000001   rmsk    TP      hap=3   0       mosaic
	
	all_input <- read.delim(train_file,header=TRUE)
	all_input <- all_input[apply(all_input,1,my_chrXY),]
	all_input <-all_input[!is.na(all_input$mosaic_likelihood),]
	all_input$mapq_p[is.na(all_input$mapq_p)]<-1
	all_input <- all_input[complete.cases(all_input[,seq(1,28)]),]
	
	all_input.2 <- subset(all_input, select=c(querypos_p,leftpos_p, seqpos_p, mapq_p, baseq_p, baseq_t, ref_baseq1b_p, ref_baseq1b_t, alt_baseq1b_p, alt_baseq1b_t, sb_p, major_mismatches_mean, minor_mismatches_mean, mismatches_p, AF, dp, mosaic_likelihood, het_likelihood, refhom_likelihood, mapq_difference, sb_read12_p, dp_diff, conflict_num, mappability,ref_softclip, alt_softclip, indel_proportion_SNPonly, alt2_proportion_SNPonly))
	all_input.3<- all_input.2
	all_input.3$querypos_p=log(all_input.3$querypos_p+1e-7)
	all_input.3$leftpos_p=log(all_input.3$leftpos_p+1e-7)
	all_input.3$seqpos_p=log(all_input.3$seqpos_p+1e-7)
	all_input.3$mapq_p=log(all_input.3$mapq_p+1e-7)
	all_input.3$baseq_p=log(all_input.3$baseq_p+1e-7)
	all_input.3$ref_baseq1b_p=log(all_input.3$ref_baseq1b_p+1e-7)
	all_input.3$alt_baseq1b_p=log(all_input.3$alt_baseq1b_p+1e-7)
	all_input.3$sb_p=log(all_input.3$sb_p+1e-7)
	all_input.3$mismatches_p=log(all_input.3$mismatches_p+1e-7)
	all_input.3$sb_read12_p=log(all_input.3$sb_read12_p+1e-7)
	all_input.3$indel_proportion_SNPonly=log(all_input.3$indel_proportion_SNPonly+1e-7)
	all_input.3$alt2_proportion_SNPonly=log(all_input.3$alt2_proportion_SNPonly+1e-7)
	all_input.3$major_mismatches_mean=all_input.3$major_mismatches_mean*read_length
	all_input.3$minor_mismatches_mean=all_input.3$minor_mismatches_mean*read_length
	##all_phasable.3$dp=log(all_phasable.3$dp)

	drops <- c("conflict_num")
	all_input.3<-all_input.3[ , !(names(all_input.3) %in% drops)]
	
	pc<-prcomp(all_input.3,
	       center = TRUE,
	       scale. = TRUE) 
	all_input$pc1 <- pc$x[,1]
	all_input$pc2 <- pc$x[,2]
	all_input$pc3 <- pc$x[,3]
	all_input$pc4 <- pc$x[,4]
	all_input$pc5 <- pc$x[,5]
	
	all_phasable <- subset(all_input, phase != "notphased")
	
	#dp_p    querypos_p      leftpos_p       seqpos_p        mapq_p  baseq_p  ref_baseq1b_p   alt_baseq1b_p  sb_p   mismatches_p    sb_read12_p     
	set.seed(123)
	all_train <- all_phasable[!is.na(all_phasable$validation),]
	all_train.2 <- subset(all_train, select=c(phase, validation, pc1, pc2, pc3, pc4, pc5))
	all_train.2 <- subset(all_train.2, phase!="hap=2")


##add a line here to balance the number of differnet validated sites:
	if (model_type=="glmnet"){
		num_het<- sum(all_train.2$validation=="het")
		num_mosaic<- sum(all_train.2$validation=="mosaic")
		num_refhom<- sum(all_train.2$validation=="refhom")
		num_repeat<- sum(all_train.2$validation=="repeat")
		if (num_repeat>(num_het+num_mosaic+num_refhom)*1){
			num_sample=round((num_het+num_mosaic+num_refhom)*1)
			all_train.2 <- rbind(subset(all_train.2,validation!="repeat"),
					subset(all_train.2,validation=="repeat")[sample(num_repeat,num_sample),])	
		}
		

		all_train.2$phase <- as.factor(all_train.2$phase)		
		model <- train(validation ~ ., all_train.2, method="glmnet",tuneGrid=expand.grid(.alpha=0:1, .lambda=0:30/10))
		saveRDS(model,prediction_model)

	}else if (model_type=="naivebayes"){	
		model <- naiveBayes(validation ~ ., all_train.2)
		saveRDS(model,prediction_model)
	}
	
	
	all_phasable.4 <- subset(all_phasable, select=c(phase, validation, pc1, pc2, pc3, pc4, pc5))
	colnames(all_phasable.4) <- c("phase","validation","pc1","pc2","pc3","pc4","pc5")
	all_phasable_nonhet <- subset(all_phasable.4, phase!="hap=2")
	all_phasable_nonhet$phase <- as.factor(all_phasable_nonhet$phase)
	
	nonhet_phasable <- subset(all_phasable, phase!="hap=2")
	het_phasable <- subset(all_phasable, phase=="hap=2")
	nonhet_phasable$phase_model_corrected <- predict(model, all_phasable_nonhet)
	het_phasable$phase_model_corrected <- "het"
	all_phasable <- rbind(het_phasable, nonhet_phasable)
	
	write.table(all_phasable, output_file,sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)
} else if (variant_type=="INDEL") {
	#head train_phasable_sites
	#chr     pos     ref     alt     MAF     id      dp_p    querypos_p      leftpos_p       seqpos_p        mapq_p  baseq_p baseq_t ref_baseq1b_p   ref_baseq1b_t   alt_baseq1b_p   alt_baseq1b_t   sb_p    context     major_mismatches_mean   minor_mismatches_mean   mismatches_p    AF      dp      mosaic_likelihood       het_likelihood  refhom_likelihood       althom_likelihood       mapq_difference sb_read12_pdp_diff  repeats validation      phase conflicting_reads       phase_corrected
	#10      10009041        G       T       0       fb0c6353-a90c-45e2-9355-7cd16cf756ff_10_10009041_G_T    0.243532076735  0.52723 0.51954 0.04323 0.31514 0.15322 -1.60007        0.93261 -0.40356        0.49078     0.36508 1       TCA     0.855   1.923   0.04837 0.339   115     0.954322027344281       0.0456779726557187      2.52244597973152e-116   2.97816045946957e-244   -1.57895        1       -4.89285700000001   rmsk    TP      hap=3   0       mosaic
	
	all_input <- read.delim(train_file,header=TRUE)
	all_input <- all_input[apply(all_input,1,my_chrXY),]
	all_input <-all_input[!is.na(all_input$mosaic_likelihood),]
	all_input$mapq_p[is.na(all_input$mapq_p)]<-1
	all_input <- all_input[complete.cases(all_input[,seq(1,28)]),]
	
	all_input.2 <- subset(all_input, select=c(querypos_p,leftpos_p, seqpos_p, mapq_p, baseq_p, baseq_t, ref_baseq1b_p, ref_baseq1b_t, alt_baseq1b_p, alt_baseq1b_t, sb_p, major_mismatches_mean, minor_mismatches_mean, mismatches_p, AF, dp, mosaic_likelihood, het_likelihood, refhom_likelihood, mapq_difference, sb_read12_p, dp_diff, conflict_num, mappability,ref_softclip, alt_softclip,length))
	all_input.3<- all_input.2
	all_input.3$querypos_p=log(all_input.3$querypos_p+1e-7)
	all_input.3$leftpos_p=log(all_input.3$leftpos_p+1e-7)
	all_input.3$seqpos_p=log(all_input.3$seqpos_p+1e-7)
	all_input.3$mapq_p=log(all_input.3$mapq_p+1e-7)
	all_input.3$baseq_p=log(all_input.3$baseq_p+1e-7)
	all_input.3$ref_baseq1b_p=log(all_input.3$ref_baseq1b_p+1e-7)
	all_input.3$alt_baseq1b_p=log(all_input.3$alt_baseq1b_p+1e-7)
	all_input.3$sb_p=log(all_input.3$sb_p+1e-7)
	all_input.3$mismatches_p=log(all_input.3$mismatches_p+1e-7)
	all_input.3$sb_read12_p=log(all_input.3$sb_read12_p+1e-7)
	all_input.3$major_mismatches_mean=all_input.3$major_mismatches_mean*read_length
	all_input.3$minor_mismatches_mean=all_input.3$minor_mismatches_mean*read_length
	##all_phasable.3$dp=log(all_phasable.3$dp)
	
	pc<-prcomp(all_input.3,
	       center = TRUE,
	       scale. = TRUE) 
	all_input$pc1 <- pc$x[,1]
	all_input$pc2 <- pc$x[,2]
	all_input$pc3 <- pc$x[,3]
	all_input$pc4 <- pc$x[,4]
	all_input$pc5 <- pc$x[,5]
	
	all_phasable <- subset(all_input, phase != "notphased")
	
	#dp_p    querypos_p      leftpos_p       seqpos_p        mapq_p  baseq_p  ref_baseq1b_p   alt_baseq1b_p  sb_p   mismatches_p    sb_read12_p     
	set.seed(123)
	all_train <- all_phasable[!is.na(all_phasable$validation),]
	all_train.2 <- subset(all_train, select=c(phase, validation, pc1, pc2, pc3, pc4, pc5))
	all_train.2 <- subset(all_train.2, phase!="hap=2")


##add a line here to balance the number of differnet validated sites:
	if (model_type=="glmnet"){
		num_het<- sum(all_train.2$validation=="het")
		num_mosaic<- sum(all_train.2$validation=="mosaic")
		num_refhom<- sum(all_train.2$validation=="refhom")
		num_repeat<- sum(all_train.2$validation=="repeat")
		if (num_repeat>(num_het+num_mosaic+num_refhom)*1){
			num_sample=round((num_het+num_mosaic+num_refhom)*1)
			all_train.2 <- rbind(subset(all_train.2,validation!="repeat"),
					subset(all_train.2,validation=="repeat")[sample(num_repeat,num_sample),])	
		}
		

		all_train.2$phase <- as.factor(all_train.2$phase)		
		model <- train(validation ~ ., all_train.2, method="glmnet",tuneGrid=expand.grid(.alpha=0:1, .lambda=0:30/10))
		saveRDS(model,prediction_model)

	}else if (model_type=="naivebayes"){	
		model <- naiveBayes(validation ~ ., all_train.2)
		saveRDS(model,prediction_model)
	}
	
	
	all_phasable.4 <- subset(all_phasable, select=c(phase, validation, pc1, pc2, pc3, pc4, pc5))
	colnames(all_phasable.4) <- c("phase","validation","pc1","pc2","pc3","pc4","pc5")
	all_phasable_nonhet <- subset(all_phasable.4, phase!="hap=2")
	all_phasable_nonhet$phase <- as.factor(all_phasable_nonhet$phase)
	
	nonhet_phasable <- subset(all_phasable, phase!="hap=2")
	het_phasable <- subset(all_phasable, phase=="hap=2")
	nonhet_phasable$phase_model_corrected <- predict(model, all_phasable_nonhet)
	het_phasable$phase_model_corrected <- "het"
	all_phasable <- rbind(het_phasable, nonhet_phasable)
	
	write.table(all_phasable, output_file,sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)
}



 



library(mlr)
df <- all_train.2
df$validation <- gsub("repeat","repeats",df$validation)
#df$phase <- as.character(df$phase)
df_mosaic <- subset(df,phase=="hap=3")
df_repeat <- subset(df,phase=="hap>3")
df_het<- subset(all_train, select=c(phase, validation, pc1, pc2, pc3, pc4, pc5))
df_het <- subset(all_train,phase=="hap=2")
df_het$validation <- gsub("repeat","repeats",df_het$validation)
df_mosaic$phase <- as.factor(as.character(df_mosaic$phase))
df_repeat$phase <- as.factor(as.character(df_repeat$phase))
df_het$phase <- as.factor(as.character(df_het$phase))


learnerGLMN=makeLearner(id="Elasticnet","classif.glmnet", predict.type = "prob")
taskmosaic=makeClassifTask(data=df_mosaic,target="validation")
taskhet=makeClassifTask(data=df_het,target="validation")
taskrepeat=makeClassifTask(data=df_repeat,target="validation")

pdf(output_pdf,width=6, height=5)

plotLearnerPrediction(learnerGLMN,taskmosaic,features=c("pc1","pc2"),cv=100L,gridsize=100)+
scale_fill_manual(values=c(mosaic="#ffae00", het=brewer.pal(9,"Set3")[9], refhom=brewer.pal(9,"Set3")[4], repeats=brewer.pal(8,"Set3")[5]))+
	  	theme_bw()+
		ggtitle("Hap=3")
plotLearnerPrediction(learnerGLMN,taskhet,features=c("pc1","pc2"),cv=100L,gridsize=100)+
scale_fill_manual(values=c(mosaic="#ffae00", het=brewer.pal(9,"Set3")[9], refhom=brewer.pal(9,"Set3")[4], repeats=brewer.pal(8,"Set3")[5]))+
	  	theme_bw()+
		ggtitle("Hap=2")
plotLearnerPrediction(learnerGLMN,taskrepeat,features=c("pc1","pc2"),cv=100L,gridsize=100)+
scale_fill_manual(values=c(mosaic="#ffae00", het=brewer.pal(9,"Set3")[9], refhom=brewer.pal(9,"Set3")[4], repeats=brewer.pal(8,"Set3")[5]))+
	  	theme_bw()+
		ggtitle("Hap>3")

dev.off()
