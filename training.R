
setwd("/Users/yanmei/Documents/RIPs/forR")


trimend_mismatches_1465<- read.delim("1465_trimend_mismatches",header=TRUE)



install.packages("rARPACK")
library(rARPACK)

setwd("/Users/yanmei/Documents/RIPs/forR/")
all_train <- read.delim("all_phasable_addcluster.2",header=TRUE)
all_train.2 <- subset(all_train,select=-c(althom_likelihood,phase,id,prediction,SV,phyloP,phastCons,
                                          diagnosis,cluster,Func.refGene,Gene.refGene,
                                          ExonicFunc.refGene,gene,pLI_score,Cytoband,
                                          start,end,obs,expected,var,context,kmeans_cluster,
                                          validation))

all_train.3<- all_train.2
all_train.3$querypos_p=log(all_train.3$querypos_p+1e-7)
all_train.3$leftpos_p=log(all_train.3$leftpos_p+1e-7)
all_train.3$seqpos_p=log(all_train.3$seqpos_p+1e-7)
all_train.3$mapq_p=log(all_train.3$mapq_p+1e-7)
all_train.3$baseq_p=log(all_train.3$baseq_p+1e-7)
all_train.3$ref_baseq1b_p=log(all_train.3$ref_baseq1b_p+1e-7)
all_train.3$alt_baseq1b_p=log(all_train.3$alt_baseq1b_p+1e-7)
all_train.3$sb_p=log(all_train.3$sb_p+1e-7)
all_train.3$mismatches_p=log(all_train.3$mismatches_p+1e-7)
all_train.3$sb_read12_p=log(all_train.3$sb_read12_p+1e-7)


all.non_train<- read.delim("all.non_train",header=TRUE)
all.non_train.2 <- subset(all.non_train,select=-c(althom_likelihood,phase,id,phase,SV,phyloP,
                                                  phastCons,
                                          diagnosis,cluster,Func.refGene,Gene.refGene,
                                          ExonicFunc.refGene,gene,pLI_score,Cytoband,
                                          start,end,obs,expected,var,context))

all.non_train.3<- subset(all.non_train.2,select=-c(prediction))
all.non_train.3$querypos_p=log(all.non_train.3$querypos_p+1e-7)
all.non_train.3$leftpos_p=log(all.non_train.3$leftpos_p+1e-7)
all.non_train.3$seqpos_p=log(all.non_train.3$seqpos_p+1e-7)
all.non_train.3$mapq_p=log(all.non_train.3$mapq_p+1e-7)
all.non_train.3$baseq_p=log(all.non_train.3$baseq_p+1e-7)
all.non_train.3$ref_baseq1b_p=log(all.non_train.3$ref_baseq1b_p+1e-7)
all.non_train.3$alt_baseq1b_p=log(all.non_train.3$alt_baseq1b_p+1e-7)
all.non_train.3$sb_p=log(all.non_train.3$sb_p+1e-7)
all.non_train.3$mismatches_p=log(all.non_train.3$mismatches_p+1e-7)
all.non_train.3$sb_read12_p=log(all.non_train.3$sb_read12_p+1e-7)
all.non_train.3$dp=log(all.non_train.3$dp)






pc<-prcomp(all_train.3,
       center = TRUE,
       scale. = TRUE) 

library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(pc, obs.scale = 1, var.scale = 1, 
              groups = all_train$phase, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

##cluster:
g <- ggbiplot(pc, obs.scale = 1, var.scale = 1, 
              groups = as.factor(all_train$kmeans_cluster), ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

##cluster:
g <- ggbiplot(pc, obs.scale = 1, var.scale = 1, 
              groups = as.factor(all_train$validation), ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)



install.packages("tsne")
install.packages("rgl")
install.packages("FactoMineR")
install.packages("vegan")

library(tsne)
library(rgl)
library(FactoMineR)
library(vegan)

pc3d<-cbind(pc$x[,1], pc$x[,2], pc$x[,3])
plot3d(pc3d, col = as.numeric(as.factor(all_train$phase)),type="s",size=1,scale=0.2)

pc3d<-cbind(pc$x[,1], pc$x[,2], pc$x[,3])
plot3d(pc3d, col = as.numeric(all_train$kmeans_cluster,type="s",size=1,scale=0.2))

##validated sites:
pc3d<-cbind(pc$x[,1], pc$x[,2], pc$x[,3])
plot3d(pc3d, col = as.numeric(all_train$validation,type="s",size=1,scale=0.2))



all_train.4 <- all_train.3
all_train.4$phase <- as.numeric(as.factor(all_train$phase))
modFit <-train(phase ~., data=all_train.4, method="lasso")
plot.enet(modFit$finalModel)








pc2<-prcomp(all.non_train.3,
           center = TRUE,
           scale. = TRUE) 

library(devtools)
library(ggbiplot)
g <- ggbiplot(pc2, obs.scale = 1, var.scale = 1, 
              groups = all.non_train$prediction, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

library(tsne)
library(rgl)
library(FactoMineR)
library(vegan)

pc3d<-cbind(pc2$x[,1], pc2$x[,2], pc2$x[,3])
plot3d(pc3d, col = as.numeric(as.factor(all.non_train$prediction)),type="s",size=1,scale=0.2)


all_train <- read.delim("all.train",header=TRUE)
all_train.2 <- subset(all_train,select=-c(althom_likelihood,prediction,SV,phyloP,phastCons,
                                          diagnosis,cluster,Func.refGene,Gene.refGene,
                                          ExonicFunc.refGene,gene,pLI_score,Cytoband,
                                          start,end,obs,expected,var,context))
all_train_kmeans<- all_train.2
all_train_kmeans$querypos_p<-kmeans(all_train.2$querypos_p,2)[[1]]
all_train_kmeans$leftpos_p<- kmeans(all_train_kmeans$leftpos_p,2)[[1]]
all_train_kmeans$seqpos_p<-kmeans(all_train_kmeans$seqpos_p,2)[[1]]
all_train_kmeans$mapq_p<- kmeans(all_train_kmeans$mapq_p,2)[[1]]
all_train_kmeans$baseq_p<-kmeans(all_train_kmeans$baseq_p,2)[[1]]
all_train_kmeans$baseq_t<-kmeans(all_train_kmeans$baseq_t,2)[[1]]
all_train_kmeans$ref_baseq1b_p<-kmeans(all_train_kmeans$ref_baseq1b_p,2)[[1]]
all_train_kmeans$ref_baseq1b_t<-kmeans(all_train_kmeans$ref_baseq1b_t,3)[[1]]
all_train_kmeans$alt_baseq1b_p<-kmeans(all_train_kmeans$alt_baseq1b_p,2)[[1]]
all_train_kmeans$alt_baseq1b_t<-kmeans(all_train_kmeans$alt_baseq1b_t,3)[[1]]
all_train_kmeans$sb_p<-kmeans(all_train_kmeans$sb_p,2)[[1]]
all_train_kmeans$major_mismatches_mean<-kmeans(all_train_kmeans$major_mismatches_mean,2)[[1]]
all_train_kmeans$minor_mismatches_mean<-kmeans(all_train_kmeans$minor_mismatches_mean,2)[[1]]
all_train_kmeans$mismatches_p<- kmeans(all_train_kmeans$mismatches_p,2)[[1]]
all_train_kmeans$AF<-kmeans(all_train_kmeans$AF,3)[[1]]
all_train_kmeans$dp<-kmeans(all_train_kmeans$dp,3)[[1]]
all_train_kmeans$mapq_difference<- kmeans(all_train_kmeans$mapq_difference,3)[[1]]
all_train_kmeans$sb_read12_p<-kmeans(all_train_kmeans$sb_read12_p,2)[[1]]
all_train_kmeans$dp_diff<-kmeans(all_train_kmeans$dp_diff,3)[[1]]
all_train_kmeans$mosaic_likelihood<-kmeans(all_train_kmeans$mosaic_likelihood,3)[[1]]
all_train_kmeans$het_likelihood<-kmeans(all_train_kmeans$het_likelihood,3)[[1]]
all_train_kmeans$refhom_likelihood<-kmeans(all_train_kmeans$refhom_likelihood,2)[[1]]

all_train_kmeans<-subset(all_train_kmeans,select=-c(leftpos_p))
install.packages("lsa")
library(lsa)
cosine_similariy <- cosine(t(as.matrix(all_train_kmeans[,-c(1,2)])))
text2vec v0.5.0
by Dmitriy Selivanov
install.packages("proxy")
library(proxy)
dat2<- as.matrix(all_train_kmeans[,-c(1,2)])
jaccard_similarity <- dist(dat2,by_rows=TRUE,method="Jaccard")


library(data.table)
data <- matrix(round(runif(25*25,0,25)),ncol=25)
tdata <- t(data)

same=0
n=0
for (i in 1:nrow(all_train_kmeans)){
  for (j in i:nrow(all_train_kmeans)){
    n=n+1
    same[n] <- sum(all_train_kmeans[i,-c(1,2)]==all_train_kmeans[j,-c(1,2)])
  }
}


install.packages("rARPACK")
##use this package above to calculate svds quickly.
sim.matrix.2<- read.delim("sim.matrix.2",header=FALSE)


#https://stackoverflow.com/questions/25978069/creating-a-matrix-from-a-list-of-pairwise-comparisons



##kmeans of each feature && cluster:
sim.matrix.3 <- read.delim("sim.matrix.3",header=FALSE)
similarity <- matrix(sim.matrix.3[,3],3849,3849)

rownames(similarity)<- sim.matrix.3[1:3849,2]
colnames(similarity)<- sim.matrix.3[1:3849,2]

write.table(similarity,col.names = TRUE,row.names = TRUE,quote=FALSE,sep="\t","similarity.text")
tmp <- svds(similarity,4)
similarity_svds <- tmp$u
write.table(similarity_svds,col.names=TRUE,row.names=TRUE,"similarity_4svds.txt",sep="\t",quote=FALSE)
rownames(similarity_svds)<- sim.matrix.3[1:3849,2]
colnames(similarity_svds)<- paste0("v",seq(1:4))
cluster_svds <- kmeans(similarity_svds,4)
write.table(cluster_svds$cluster,col.names=TRUE,row.names=TRUE,"cluster_4svds_kmeans.txt",quote=FALSE,sep="\t")

my_distance <- function(x){
  d1<- sqrt((x[1]-cluster_svds$centers[1,1])^2+(x[2]-cluster_svds$centers[1,2])^2+
              (x[3]-cluster_svds$centers[1,3])^2+(x[4]-cluster_svds$centers[1,4])^2)
  d2 <- sqrt( (x[1]-cluster_svds$centers[2,1])^2 +(x[2]-cluster_svds$centers[2,2])^2+ (x[3]-cluster_svds$centers[2,3])^2
              + (x[4]-cluster_svds$centers[2,4])^2)
  d3 <- sqrt((x[1]-cluster_svds$centers[3,1])^2+(x[2]-cluster_svds$centers[3,2])^2+(x[3]-cluster_svds$centers[3,3])^2+
               (x[4]-cluster_svds$centers[3,4])^2 )
  d4 <- sqrt((x[1]-cluster_svds$centers[4,1])^2+(x[2]-cluster_svds$centers[4,2])^2+
               (x[3]-cluster_svds$centers[4,3])^2+(x[4]-cluster_svds$centers[4,4])^2)
  list(d1,d2,d3,d4,min(abs(d2-d3),abs(d1-d2),abs(d1-d3)) )
}

mat <- do.call("rbind",apply(tmp$u,1,my_distance))
cluster_distance <- cbind(cluster_svds$cluster,mat)
colnames(cluster_distance)<- c("kmeans_cluster","d1","d2","d3","d4","dd")
write.table(cluster_distance,col.names=TRUE,row.names=TRUE,"cluster_kmeans_adddistance.txt",quote=FALSE,sep="\t")





##graph laplician
D<-diag(1/sqrt(apply(similarity,1,sum)),nrow(similarity),ncol(similarity))
A<- similarity
similarity_normalized <- D%*%A%*%D
tmp <- svds(similarity_normalized,4)
similarity_svds <- tmp$u

rownames(similarity_svds)<- sim.matrix.3[1:3849,2]
colnames(similarity_svds)<- paste0("v",seq(1:4))
write.table(similarity_svds,col.names=TRUE,row.names=TRUE,"similarity_normalized.txt",sep="\t",quote=FALSE)
cluster_svds <- kmeans(similarity_svds,4)

my_distance <- function(x){
  d1<- sqrt((x[1]-cluster_svds$centers[1,1])^2+(x[2]-cluster_svds$centers[1,2])^2+
                                          (x[3]-cluster_svds$centers[1,3])^2+(x[4]-cluster_svds$centers[1,4])^2)
  d2 <- sqrt( (x[1]-cluster_svds$centers[2,1])^2 +(x[2]-cluster_svds$centers[2,2])^2+ (x[3]-cluster_svds$centers[2,3])^2
              + (x[4]-cluster_svds$centers[2,4])^2)
  d3 <- sqrt((x[1]-cluster_svds$centers[3,1])^2+(x[2]-cluster_svds$centers[3,2])^2+(x[3]-cluster_svds$centers[3,3])^2+
               (x[4]-cluster_svds$centers[3,4])^2 )
  d4 <- sqrt((x[1]-cluster_svds$centers[4,1])^2+(x[2]-cluster_svds$centers[4,2])^2+
               (x[3]-cluster_svds$centers[4,3])^2+(x[4]-cluster_svds$centers[4,4])^2)
  list(d1,d2,d3,d4,min(abs(d2-d3),abs(d1-d2),abs(d1-d3)) )
}

mat <- do.call("rbind",apply(tmp$u,1,my_distance))


cluster_distance <- cbind(cluster_svds$cluster,mat)
colnames(cluster_distance)<- c("kmeans_cluster","d1","d2","d3","d4","dd")
write.table(cluster_distance,col.names=TRUE,row.names=TRUE,"cluster_kmeans_normalized_adddistance.txt",quote=FALSE,sep="\t")




##multinomial logistic regression:
#https://www.reddit.com/r/statistics/comments/14fv87/help_needed_multinomial_logistic_regression_in_r/
library(nnet)
library(glmnet)
library(caret)

##only use kmeans_cluster to train the model:
all_phasable_validation <- read.delim("all_phasable_validation",header=TRUE)
validation <-subset(all_phasable_validation,select=c(kmeans_cluster,phase,validation))
set.seed(123)
model <- train(validation ~ ., validation, method="glmnet",tuneGrid=expand.grid(.alpha=0:1, .lambda=0:30/10))
confusionMatrix(predict(model,validation),validation$validation)

##use both kmeans_cluser and distance to train the model:
input_log <- read.delim("input_log",header=TRUE)
validation <-subset(input_log,select=c(kmeans_cluster,phase,validation,d1,d2,d3,d4))
set.seed(123)
model <- train(validation ~ ., validation, method="glmnet",tuneGrid=expand.grid(.alpha=0:1, .lambda=0:30/10))
confusionMatrix(predict(model,validation),validation$validation)


##use normalized kmeans_cluster and distance to train the model:
input_log_normalized <- read.delim("input_log_normalized",header=TRUE)
validation <-subset(input_log_normalized,select=c(kmeans_cluster,phase,validation,d1,d2,d3,d4,dd))
set.seed(123)
model <- train(validation ~ ., validation, method="glmnet",tuneGrid=expand.grid(.alpha=0:1, .lambda=0:30/10))
confusionMatrix(predict(model,validation),validation$validation)

##used PCA results to train the model:
write.table(pc$x,sep="\t",)

all_train <- read.delim("all_phasable_addcluster.2",header=TRUE)
all_train.2 <- subset(all_train,select=-c(althom_likelihood,phase,id,prediction,SV,phyloP,phastCons,
                                          diagnosis,cluster,Func.refGene,Gene.refGene,
                                          ExonicFunc.refGene,gene,pLI_score,Cytoband,
                                          start,end,obs,expected,var,context,kmeans_cluster,
                                          validation))

all_train.3<- all_train.2
all_train.3$querypos_p=log(all_train.3$querypos_p+1e-7)
all_train.3$leftpos_p=log(all_train.3$leftpos_p+1e-7)
all_train.3$seqpos_p=log(all_train.3$seqpos_p+1e-7)
all_train.3$mapq_p=log(all_train.3$mapq_p+1e-7)
all_train.3$baseq_p=log(all_train.3$baseq_p+1e-7)
all_train.3$ref_baseq1b_p=log(all_train.3$ref_baseq1b_p+1e-7)
all_train.3$alt_baseq1b_p=log(all_train.3$alt_baseq1b_p+1e-7)
all_train.3$sb_p=log(all_train.3$sb_p+1e-7)
all_train.3$mismatches_p=log(all_train.3$mismatches_p+1e-7)
all_train.3$sb_read12_p=log(all_train.3$sb_read12_p+1e-7)

pc<-prcomp(all_train.3,
           center = TRUE,
           scale. = TRUE)
pc_all <- as.data.frame(pc$x)
pc_all$id<- all_train$id  
write.table(pc_all,sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE,"all_phasable_pc")


pc3d<-cbind(pc$x[,1], pc$x[,2], pc$x[,3])
plot3d(pc3d, col = as.numeric(all_train$kmeans_cluster,type="s",size=1,scale=0.2))

plot(model)
coef(model$finalModel, s=model$bestTune$.lambda)

all_phasable_validation_pc4<- read.delim("all_phasable_validation_pc4",header=TRUE)

set.seed(123)
flds<-createFolds(validation$validation, k = 5, list = TRUE, returnTrain = FALSE)
#validation <-subset(all_phasable_validation_pc4,select=c(phase,validation,kmeans_cluster,PC1,PC2,PC3,PC4))
validation <-subset(all_phasable_validation_pc4,select=c(phase,validation,PC1,PC2,PC3,PC4))

training <- validation[-flds[[4]],]
testing <- validation[flds[[4]],]
model <- train(validation ~ ., training, method="glmnet",tuneGrid=expand.grid(.alpha=0:1, .lambda=0:30/10))
confusionMatrix(predict(model,testing),testing$validation)

model <- train(validation ~ ., validation, method="glmnet",tuneGrid=expand.grid(.alpha=0:1, .lambda=0:30/10))
confusionMatrix(predict(model,validation),validation$validation)





all_phasable_validation_pc4_addpredict<- all_phasable_validation_pc4
all_phasable_validation_pc4_addpredict$predict<- predict(model,validation)

pc3d<-cbind(all_phasable_validation_pc4_addpredict$PC1,
            all_phasable_validation_pc4_addpredict$PC2,
            all_phasable_validation_pc4_addpredict$PC3)
plot3d(pc3d, col = as.numeric(all_phasable_validation_pc4_addpredict$predict,type="s",size=1,scale=0.2))
plot3d(pc3d, col = as.numeric(all_phasable_validation_pc4_addpredict$phase,type="s",size=1,scale=0.2))

plot3d(pc3d, col = as.numeric(all_phasable_validation_pc4_addpredict$validation,type="s",size=1,scale=0.2))


##re-predict all phasable sites:
all_phasable_pc4<- read.delim("all_phasable_pc4.2",header=TRUE)
predict(model, all_phasable_pc4)
all_phasable_pc4$predict_PC<-predict(model, all_phasable_pc4)

pc3d<-cbind(all_phasable_pc4$PC1,
            all_phasable_pc4$PC2,
            all_phasable_pc4$PC3)
plot3d(pc3d, col = as.numeric(all_phasable_pc4$predict_PC,type="s",size=1,scale=0.2))
plot3d(pc3d, col = as.numeric(all_phasable_pc4$validation,type="s",size=1,scale=0.2))

plot3d(pc3d, col = as.numeric(all_phasable_pc4$phase,type="s",size=1,scale=0.2))


pc<-prcomp(all_train.3,
           center = TRUE,
           scale. = TRUE)
pc_all <- as.data.frame(pc$x)
pc_all$id<- all_train$id  
write.table(all_phasable_pc4,sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE,"all_phasable_pc4_addprediction.txt")



x<-c(0.45,81/84,14/31)
x<- t(as.matrix(x))
colnames(x)<- c("mosaic","het","mult_hap")
rownames(x)<- "precision"
par(mfrow=c(1,2))
barplot(x,col = c("lightblue", "mistyrose", "lightcyan"),beside=TRUE,ylim=c(0,1),
        ylab="precision",
        main="Phasable sites")
x<-c(0.54,0.957,0.6621)
x<- t(as.matrix(x))
colnames(x)<- c("mosaic","het","mult_hap")
rownames(x)<- "precision"
barplot(x,col = c("lightblue", "mistyrose", "lightcyan"),beside=TRUE,ylim=c(0,1),
        ylab="precision",
        main="Predicted sites")




x<-c(0.614,96/99,0.6875)
x<- t(as.matrix(x))
colnames(x)<- c("mosaic","het","refhom/mult")
rownames(x)<- "precision"
par(mfrow=c(1,2))
barplot(x,col = c("lightblue", "mistyrose", "lightcyan"),beside=TRUE,ylim=c(0,1),
        ylab="precision",
        main="Phasable sites")
x<-c(0.6471,0.9545,0.6944)
x<- t(as.matrix(x))
colnames(x)<- c("mosaic","het","refhom/mult")
rownames(x)<- "precision"
barplot(x,col = c("lightblue", "mistyrose", "lightcyan"),beside=TRUE,ylim=c(0,1),
        ylab="precision",
        main="Predicted sites")

par(mfrow=c(1,1))
x<-c(0.7317,0.9620,0.6363)
x<- t(as.matrix(x))
colnames(x)<- c("mosaic","het","refhom/mult")
rownames(x)<- "precision"
barplot(x,col = c("lightblue", "mistyrose", "lightcyan"),beside=TRUE,ylim=c(0,1),
        ylab="precision",
        main="Predicted sites without clustered sites")




