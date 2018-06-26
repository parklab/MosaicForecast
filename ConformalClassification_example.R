source("ConformalClassification.R")

# Define:
X_train=
Y_train=
X_train=

algorithm <- "rf"

# trControl <- trainControl(method="repeatedcv", 
#                           number=10, repeats=3, 
#                           search="grid",
#                           savePredictions=TRUE)
# 
# trControl <- trainControl(method = "cv",  number=5,savePredictions=TRUE)
# tunegrid <- expand.grid(.mtry=30)
# metric <- "Accuracy"
# 
# nb_trees=500
# rf_gridsearch <- train(LogSDescsTrain, LogSTrain, method="rf", #metric=metric,
#                        #tuneGrid=tunegrid, 
#                        trControl=trControl,
#                        type="classification", predict.all=TRUE,
#                        keep.forest=TRUE,norm.votes=TRUE,ntree=nb_trees)
# 

kfolds=10
trControl <- trainControl(method = "cv",  number=10,savePredictions=TRUE)
set.seed(3)
nb_trees <- 100
model <- train(X_train, Y_train,
               algorithm,type="classification",
               trControl=trControl,predict.all=TRUE,
               keep.forest=TRUE,norm.votes=TRUE,
               ntree=nb_trees)


# Instantiate the class and get the p.values
example <- ConformalClassification$new()
example$CalculateCVScores(model=model)
example$CalculatePValues(new.data=X_test)
# we get the p.values:
example$p.values$P.values
# we get the significance of these p.values.
example$p.values$Significance_p.values
example$ClassPredictions$aggregate
example$ClassPredictions$individual

model$trainingData$.outcome

