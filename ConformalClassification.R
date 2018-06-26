#############################################
### Conformal Prediction for Classification
#############################################

ConformalClassification <- setRefClass(
  "ConformalClassification",
  fields = list(
    ClassificationModel = "ANY",
    confidence = "numeric",
    data.new = "ANY",
    NonconformityScoresMatrix ="ANY",
    ClassPredictions = "ANY",
    p.values = "ANY"
  ),
  methods = list(
    initialize = function(confi = 0.8)
    {
      "This method is called when you create an instance of the class."
      if (confi > 1 || confi < 0)
        stop("Confidence must be between 0 and 1")
      confidence <<- confi
      cat("Conformal Prediction Class for Classification Instantiated")
      cat("\n")
    },
    CalculateCVScores = function(model=NULL)
    {
      if(is.null(model))
        stop("To calculate the alphas, a point prediction model and an error model 
             need to be suppplied")
      if(model$modelType != "Classification" )
        stop("The model type needs to be 'Classification'")     
      ClassificationModel <<- model
      print("Calculating the vector of nonconformity measures for the CV predictions (label wise Mondrian ICP)..")
      cat('\n')
      MondrianICP <- data.frame( model$finalModel$votes )
      names(MondrianICP) = levels(model$trainingData$.outcome)
      # add true labels
      #MondrianICP <- as.data.frame(MondrianICP)
      MondrianICP$label <- as.vector( model$trainingData$.outcome )
      # consider only the values for one class
      for (lab in levels(model$trainingData$.outcome) ){
        MondrianICP[which(MondrianICP$label != lab),lab] = NA
        MondrianICP[,lab] <- sort(MondrianICP[,lab], decreasing=F,na.last = T)
      }
      ###MondrianICP[,1:2] <- apply(MondrianICP[,1:2], 2, sort, decreasing=FALSE)
      MondrianICP$label=NULL
      NonconformityScoresMatrix <<- MondrianICP
    },
    CalculatePValues = function(new.data=NULL)
    {
      if (is.null(new.data)){
        stop("\nArgument 'data.new' cannot be empty.\nNew datapoints are required as input\n")
      }
      else{
        data.new <<- new.data
      }
      #require("caret") || stop("Pacakge 'caret' is required to make new predictions")
      
      print("Classifying the input data..")
      cat('\n')
      pred <- predict(ClassificationModel$finalModel, newdata = new.data,predict.all=TRUE) # individual or aggregate
      ClassPredictions <<- pred
      ntrees <- ClassificationModel$finalModel$ntree
      
      library(plyr)
      votes <- alply(pred$individual, 1, function(x) { table(x) })
      ##votes <- apply(pred$individual,1,function(x){table(x)})
      
      out<-c()
      for (i in colnames(NonconformityScoresMatrix)){
        out<-cbind(out,sapply(votes,function(x) x[i]))
      }
      out[is.na(out)] <- 0
      out <- out/ntrees
      colnames(out) <- colnames(NonconformityScoresMatrix)
      
      pval <- t(apply(out,1,function(x){ apply(do.call(rbind, lapply(as.data.frame(t(NonconformityScoresMatrix)), "<", x)),2,sum,na.rm=T)    }))
      pval <- pval / nrow(NonconformityScoresMatrix)
      # this also works but is slower
      # library(plyr)
      # now <- t(apply(out,1,function(x){ apply(aaply(NonconformityScoresMatrix, 1, "<", x),2,sum)    }))
      # http://stackoverflow.com/questions/20596433/how-to-divide-each-row-of-a-matrix-by-elements-of-a-vector-in-r
      pval_signif <- (pval > (1-confidence))*1
      p.values <<- list(P.values = pval,Significance_p.values = pval_signif)
    }
    
  )
)
