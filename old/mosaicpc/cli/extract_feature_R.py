
RSCRIPT='''#!/usr/bin/env Rscript
# .libPaths( c( .libPaths(), "/n/data1/hms/dbmi/park/yanmei/tools/R_packages/") )

args<-commandArgs(TRUE)
if (length(args)!=1) {
    stop("Rscript extraction_feature.R input_file(feature_list_from python)", call.=FALSE)
} else if (length(args)==1) {
    verbose=TRUE
    input_file=args[1]
}
output_file=paste(input_file,".list_R",sep="")

#id      querypos_major  querypos_minor  leftpos_major   leftpos_minor   seqpos_major    seqpos_minor    mapq_major      mapq_minor      baseq_major     baseq_minor     baseq_major_near1b      baseq_minor+near1b  major_plus      major_minus     minor_plus      minor_minus     context1        context2        context1_count  context2_count  mismatches_major        mismatches_minor        major_read1     major_read2 minor_read1     minor_read2     dp_near dp_far  dp_p

input=read.delim(input_file,header=TRUE,sep=" ")
input <- subset(input, (((querypos_minor!="," & seqpos_minor!=",") & seqpos_major!="," )  & baseq_minor_near1b!=",") & leftpos_minor!=",")

my.wilcox.p <- function(x)
{
    ref_pos=as.numeric(unlist(strsplit(as.character(x[2]),",")))
    alt_pos=as.numeric(unlist(strsplit(as.character(x[3]),",")))
    round(wilcox.test(ref_pos,alt_pos)$p.value,5)
}
input$querypos_p=apply(input,1,my.wilcox.p)

my.wilcox.p <- function(x)
{
    ref_pos=as.numeric(unlist(strsplit(as.character(x[4]),",")))
    alt_pos=as.numeric(unlist(strsplit(as.character(x[5]),",")))
    round(wilcox.test(ref_pos,alt_pos)$p.value,5)
}
input$leftpos_p=apply(input,1,my.wilcox.p)

my.wilcox.p <- function(x)
{
    ref_pos=as.numeric(unlist(strsplit(as.character(x[6]),",")))
    alt_pos=as.numeric(unlist(strsplit(as.character(x[7]),",")))
    round(wilcox.test(ref_pos,alt_pos)$p.value,5)
}
input$seqpos_p=apply(input,1,my.wilcox.p)

my.wilcox.p <- function(x)
{
    ref_mapq=as.numeric(unlist(strsplit(as.character(x[8]),",")))
    alt_mapq=as.numeric(unlist(strsplit(as.character(x[9]),",")))
    round(wilcox.test(ref_mapq,alt_mapq)$p.value,5)
}
input$mapq_p=apply(input,1,my.wilcox.p)


my.t.test.t.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
    #if (is(obj, "try-error")) return(NA) else return(obj$statistic)
    if (attr(obj,"class")=="try-error") return(NA) else return(obj$statistic) 
}

#my.t.test <- function(x)
#{
#	ref_mapq=as.numeric(unlist(strsplit(as.character(x[8]),",")))
#	alt_mapq=as.numeric(unlist(strsplit(as.character(x[9]),",")))
#	round(my.t.test.t.value(ref_mapq,alt_mapq),5)
#}
#input$mapq_t=apply(input,1,my.t.test)

my.wilcox.p <- function(x)
{
    ref_baseq=as.numeric(unlist(strsplit(as.character(x[10]),",")))
    alt_baseq=as.numeric(unlist(strsplit(as.character(x[11]),",")))
    round(wilcox.test(ref_baseq,alt_baseq)$p.value,5)
}
input$baseq_p=apply(input,1,my.wilcox.p)

my.t.test <- function(x)
{
    ref_baseq=as.numeric(unlist(strsplit(as.character(x[10]),",")))
    alt_baseq=as.numeric(unlist(strsplit(as.character(x[11]),",")))
    round(my.t.test.t.value(ref_baseq,alt_baseq),5)	
}
input$baseq_t=apply(input,1,my.t.test)

my.wilcox.p <- function(x)
{
    ref_baseq=as.numeric(unlist(strsplit(as.character(x[10]),",")))
    ref_baseq1b=as.numeric(unlist(strsplit(as.character(x[12]),",")))
    round(wilcox.test(ref_baseq,ref_baseq1b)$p.value,5)
}
input$ref_baseq1b_p=apply(input,1,my.wilcox.p)

my.t.test <- function(x)
{
    ref_baseq=as.numeric(unlist(strsplit(as.character(x[10]),",")))
    ref_baseq1b=as.numeric(unlist(strsplit(as.character(x[12]),",")))
    round(my.t.test.t.value(ref_baseq,ref_baseq1b),5)
}
input$ref_baseq1b_t=apply(input,1,my.t.test)

my.wilcox.p <- function(x)
{
    alt_baseq=as.numeric(unlist(strsplit(as.character(x[11]),",")))
    alt_baseq1b=as.numeric(unlist(strsplit(as.character(x[13]),",")))
    round(wilcox.test(alt_baseq,alt_baseq1b)$p.value,5)
}
input$alt_baseq1b_p=apply(input,1,my.wilcox.p)

my.t.test <- function(x)
{
    alt_baseq=as.numeric(unlist(strsplit(as.character(x[11]),",")))
    alt_baseq1b=as.numeric(unlist(strsplit(as.character(x[13]),",")))
    round(my.t.test.t.value(alt_baseq,alt_baseq1b),5)
}
input$alt_baseq1b_t=apply(input,1,my.t.test)


my.fisher.p<- function(x)
{
    ref_plus <- as.numeric(unlist(as.character(x[14])))
    ref_minus <- as.numeric(unlist(as.character(x[15])))
    alt_plus <- as.numeric(unlist(as.character(x[16])))
    alt_minus <- as.numeric(unlist(as.character(x[17])))
    round(fisher.test(cbind(c(ref_plus,ref_minus),c(alt_plus,alt_minus)))$p.value,5)
}
input$sb_p=apply(input,1,my.fisher.p)

my.context_selection <- function(x)
{
    if (x[20]>=x[21]) {
    x[18]
    } else{
    x[19]
    }
}
input$context=apply(input,1,my.context_selection)

my.mean <-function(x)
{
    mismatch_major=as.numeric(unlist(strsplit(as.character(x[22]),",")))
    round(mean(mismatch_major),3)
}
input$major_mismatches_mean=apply(input,1,my.mean)

my.mean <-function(x)
{
    mismatch_minor=as.numeric(unlist(strsplit(as.character(x[23]),",")))
    round(mean(mismatch_minor),3)
}
input$minor_mismatches_mean=apply(input,1,my.mean)

my.wilcox.p <- function(x)
{
    mismatch_major=as.numeric(unlist(strsplit(as.character(x[22]),",")))
    mismatch_minor=as.numeric(unlist(strsplit(as.character(x[23]),",")))-1
    round(wilcox.test(mismatch_major,mismatch_minor)$p.value,5)
}
input$mismatches_p=apply(input,1,my.wilcox.p)

my.AF <- function(x)
{
    depth <- sum(as.numeric(unlist(as.character(x[14]))), as.numeric(unlist(as.character(x[15]))), as.numeric(unlist(as.character(x[16]))), as.numeric(unlist(as.character(x[17]))))
    alt <- sum(as.numeric(unlist(as.character(x[16]))), as.numeric(unlist(as.character(x[17]))))
    round(alt/depth,3)
}
input$AF=apply(input,1,my.AF)

my.dp <- function(x)
{
    depth <- sum(as.numeric(unlist(as.character(x[14]))), as.numeric(unlist(as.character(x[15]))), as.numeric(unlist(as.character(x[16]))), as.numeric(unlist(as.character(x[17]))))
    depth
}
input$dp=apply(input,1,my.dp)

my.mosaic_likelihood <- function(x)
{
    depth <- sum(as.numeric(unlist(as.character(x[14]))), as.numeric(unlist(as.character(x[15]))), as.numeric(unlist(as.character(x[16]))), as.numeric(unlist(as.character(x[17]))))
    alt <- sum(as.numeric(unlist(as.character(x[16]))), as.numeric(unlist(as.character(x[17]))))
    baseq_major <- as.numeric(unlist(strsplit(as.character(x[10]),",")))
    baseq_minor <- as.numeric(unlist(strsplit(as.character(x[11]),",")))
    r=0
    for (i in 1: length(baseq_major)){r=r+0.1^(baseq_major[i]/10)}
    for (i in 1: length(baseq_minor)){r=r+1-0.1^(baseq_minor[i]/10)}
    log10(beta(r+1,depth-r+1))
}
input$mosaic_likelihood=apply(input,1,my.mosaic_likelihood)

my.het_likelihood <- function(x)
{
    depth <- sum(as.numeric(unlist(as.character(x[14]))), as.numeric(unlist(as.character(x[15]))), as.numeric(unlist(as.character(x[16]))), as.numeric(unlist(as.character(x[17]))))
    alt <- sum(as.numeric(unlist(as.character(x[16]))), as.numeric(unlist(as.character(x[17]))))
    log10(0.5^depth)
}
input$het_likelihood=apply(input,1,my.het_likelihood)

my.refhom_likelihood <- function(x)
{
    baseq_major <- as.numeric(unlist(strsplit(as.character(x[10]),",")))
    baseq_minor <- as.numeric(unlist(strsplit(as.character(x[11]),",")))
    q=log10(1)
    for (i in 1: length(baseq_major)){q=q+log10(1-0.1^(baseq_major[i]/10))}
    for (i in 1: length(baseq_minor)){q=q+log10(0.1^(baseq_minor[i]/10))}
    q
}
input$refhom_likelihood=apply(input,1,my.refhom_likelihood)

my.althom_likelihood <- function(x)
{
    baseq_major <- as.numeric(unlist(strsplit(as.character(x[10]),",")))
    baseq_minor <- as.numeric(unlist(strsplit(as.character(x[11]),",")))
    q=log10(1)
    for (i in 1: length(baseq_minor)){q=q+log10(1-0.1^(baseq_minor[i]/10))}
    for (i in 1: length(baseq_major)){q=q+log10(0.1^(baseq_major[i]/10))}
    q
}
input$althom_likelihood=apply(input,1,my.althom_likelihood)

my.mean <- function(x)
{
    ref_mapq=as.numeric(unlist(strsplit(as.character(x[8]),",")))
    alt_mapq=as.numeric(unlist(strsplit(as.character(x[9]),",")))
    round(mean(ref_mapq)-mean(alt_mapq),5)
}
input$mapq_difference=apply(input,1,my.mean)

my.fisher.p<- function(x)
{
    ref_read1 <- as.numeric(unlist(as.character(x[24])))
    ref_read2 <- as.numeric(unlist(as.character(x[25])))
    alt_read1 <- as.numeric(unlist(as.character(x[26])))
    alt_read2 <- as.numeric(unlist(as.character(x[27])))
    round(fisher.test(cbind(c(ref_read1,ref_read2),c(alt_read1,alt_read2)))$p.value,5)
}
input$sb_read12_p=apply(input,1,my.fisher.p)

my.diff <- function(x)
{
    dp_near=as.numeric(unlist(as.character(x[28])))
    dp_far=as.numeric(unlist(as.character(x[29])))
    dp_near-dp_far
}
input$dp_diff=apply(input,1,my.diff)


output=input[,c(1,seq(30,54))]

mosaic_p <- 10^(as.numeric(as.vector(output$mosaic_likelihood)))
het_p <- 10^(as.numeric(as.vector(output$het_likelihood)))
refhom_p <- 10^(as.numeric(as.vector(output$refhom_likelihood)))
althom_p <- 10^(as.numeric(as.vector(output$althom_likelihood)))
input$normalize <- mosaic_p+het_p+refhom_p+althom_p
output$mosaic_likelihood <- mosaic_p/input$normalize
output$het_likelihood <- het_p/input$normalize
output$refhom_likelihood <- refhom_p/input$normalize
output$althom_likelihood <- althom_p/input$normalize

write.table(output,output_file,quote=F,sep="\t",row.names=F,col.names=T)
#head -1 all_putative_mosaics_feature.list.forR
#id      querypos_major  querypos_minor  leftpos_major   leftpos_minor   seqpos_major    seqpos_minor    mapq_major      mapq_minor      baseq_major     baseq_minor     baseq_major_near1b      baseq_minor+near1b  major_plus      major_minus     minor_plus      minor_minus     context_reference_forward       context_reference_reverse       context_antireference_forward   context_antireference_reverse   context_reference_forward_count     context_reference_reverse_count context_antireference_forward_count     context_antireference_reverse_count     mismatches_major        mismatches_minor        major_read1major_read2      minor_read1     minor_read2     dp_near dp_far  dp_p


'''



