CV_TDR <- function(data,outcome,loops=50,scale=TRUE,Decor=TRUE,IDSample=FALSE,filterMethod=univariate_Wilcoxon,...)
{
  SelectedTrainFeatures <- list()
  SelectedTestFeatures <- list()
  SelectedLASSOFeatures <- list()
  TDR <- numeric(loops)
  ROCAUC <- numeric(loops)
  topROCAUC <- numeric(loops)
  
  DeSelectedTrainFeatures <- list()
  DeSelectedTestFeatures <- list()
  DeSelectedLASSOFeatures <- list()
  DeTDR <- numeric(loops)
  DeROCAUC <- numeric(loops)
  DetopROCAUC <- numeric(loops)
  DeFraction <- numeric(loops)
  DeFractionLasso <- numeric(loops)
  Avlen <- numeric(loops)
  NumLatent <- numeric(loops)
  totCorrelated <- numeric(loops)
  rocTest <- list()
  datadim <- NULL
  testPredict <- NULL
  DeTestPredict <- NULL
  outvalues <- unique(data[,outcome])
  
  ### Some global cleaning
  sdiszero <- apply(data,2,var) > 1.0e-16
  data <- data[,sdiszero]

  varlist <- colnames(data)[colnames(data) != outcome]
  tokeep <- c(as.character(correlated_Remove(data,varlist,thr=0.9999)),outcome)
  data <- data[,tokeep]
  
  
  for (lp in c(1:loops))
  {
    cat("<")
    if (IDSample)
    {
      IDs <- unique(data$ID)
      sampleID <- sample(length(IDs),length(IDs)/2)
      trainIDs <- IDs[sampleID]
      testIDs <- IDs[-sampleID]
      trainSet <- data[data$ID %in% trainIDs,]
      testSet <- data[data$ID %in% testIDs,]
      trainSet$ID <- NULL
      testSet$ID <- NULL
    }
    else
    {
      firstSet <- data[data[,outcome]==outvalues[1],]
      secondSet <- data[data[,outcome]==outvalues[2],]
      sampleFirst <- sample(nrow(firstSet),nrow(firstSet)/2)
      sampleSecond <- sample(nrow(secondSet),nrow(secondSet)/2)
      trainSet <- rbind(firstSet[sampleFirst,],secondSet[sampleSecond,])
      testSet <- rbind(firstSet[-sampleFirst,],secondSet[-sampleSecond,])
      firstSet <- NULL
      secondSet <- NULL
    }
    
    datadim <- rbind(datadim,c(trainRows=nrow(trainSet),
                               testRows=nrow(testSet),
                               dataColums=ncol(trainSet),
                               table(trainSet[,outcome]),
                               table(testSet[,outcome])))
    
    
    if (scale)
    {
      trainScale <- FRESAScale(trainSet,method="OrderLogit")
      trainSet <- trainScale$scaledData
      testSet <- FRESAScale(testSet,method="OrderLogit",refMean=trainScale$refMean,refDisp=trainScale$refDisp)$scaledData
      trainScale <- NULL
    }
    ml <- LASSO_MIN(formula(paste(outcome,"~ .")),trainSet,family="binomial")
    SelectedLASSOFeatures[[lp]] <- ml$selectedfeatures
    SelectedTrainFeatures[[lp]] <- unique(c(names(filterMethod(trainSet,outcome,pvalue=0.2,limit=-1)),SelectedLASSOFeatures[[lp]]))

    
    mlt <- LASSO_MIN(formula(paste(outcome,"~ .")),testSet,family="binomial")
    SelectedTestFeatures[[lp]] <- unique(c(names(filterMethod(testSet,outcome,pvalue=0.2,limit=-1)),mlt$selectedfeatures))
    TDR[lp] <- sum(SelectedTrainFeatures[[lp]] %in% SelectedTestFeatures[[lp]])/length(SelectedTrainFeatures[[lp]])
    
    topROCAUC[lp] <- pROC::roc(testSet[,outcome],testSet[,SelectedTrainFeatures[[lp]][1]],auc=TRUE,quiet=TRUE)$auc
    pred <- as.vector(predict(ml,testSet))
    rocRaw <- pROC::roc(testSet[,outcome],pred,auc=TRUE,quiet=TRUE)
    testPredict <- rbind(testPredict,cbind(rownames(testSet),testSet[,outcome],pred))
    ROCAUC[lp] <- rocRaw$auc
    if (Decor)
    {
      cat("(")
      trainSet <- IDeA(trainSet,...)
      totCorrelated[lp] <- sqrt(attr(trainSet,"totCorrelated"))

      testSet <- predictDecorrelate(trainSet,testSet)
      cat(")")
      
      ml <- LASSO_MIN(formula(paste(outcome,"~ .")),trainSet,family="binomial")
      DeSelectedLASSOFeatures[[lp]] <- ml$selectedfeatures
      DeSelectedTrainFeatures[[lp]] <- unique(c(names(filterMethod(trainSet,outcome,pvalue=0.2,limit=-1)),DeSelectedLASSOFeatures[[lp]]))
      
      mlt <- LASSO_MIN(formula(paste(outcome,"~ .")),testSet,family="binomial")
      DeSelectedTestFeatures[[lp]] <- unique(c(names(filterMethod(testSet,outcome,pvalue=0.2,limit=-1)),mlt$selectedfeatures))
      DeTDR[lp] <- sum(DeSelectedTrainFeatures[[lp]] %in% DeSelectedTestFeatures[[lp]])/length(DeSelectedTrainFeatures[[lp]])
      
      DetopROCAUC[lp] <- pROC::roc(testSet[,outcome],testSet[,DeSelectedTrainFeatures[[lp]][1]],auc=TRUE,quiet=TRUE)$auc
      
      pred <- as.vector(predict(ml,testSet))
      rocDe <- pROC::roc(testSet[,outcome],pred,auc=TRUE,quiet=TRUE)
      DeTestPredict <- rbind(DeTestPredict,cbind(rownames(testSet),testSet[,outcome],pred))
      
      DeROCAUC[lp] <- rocDe$auc
      rocTest[[lp]] <- roc.test(rocRaw,rocDe)$p.value
      DeFraction[lp] <- sum(str_detect(DeSelectedTrainFeatures[[lp]],"La_"))/length(DeSelectedTrainFeatures[[lp]])
      if (length(DeSelectedLASSOFeatures[[lp]])>0)
        DeFractionLasso[lp] <- sum(str_detect(DeSelectedLASSOFeatures[[lp]],"La_"))/length(DeSelectedLASSOFeatures[[lp]])
      else
        DeFractionLasso[lp] <- 0;
      ltvar <- getLatentCoefficients(trainSet);
      Avlen[lp] <- mean(sapply(ltvar,length))
      NumLatent[lp] <- length(ltvar)
      
    }
    cat(">")
    if ((lp %% 10) == 0)  cat(lp,"\n")
  }
  cat("\n")
  medTestRaw <- boxplot(as.numeric(testPredict[,3])~testPredict[,1],plot=FALSE)
  medTestclass <- boxplot(as.numeric(testPredict[,2])~testPredict[,1],plot=FALSE)
  medTestRaw <- cbind(medTestclass$stats[3,],medTestRaw$stats[3,]) 
  DeMedTestRaw <- boxplot(as.numeric(DeTestPredict[,3])~DeTestPredict[,1],plot=FALSE)
  medTestclass <- boxplot(as.numeric(DeTestPredict[,2])~DeTestPredict[,1],plot=FALSE)
  DeMedTestRaw <- cbind(medTestclass$stats[3,],DeMedTestRaw$stats[3,]) 
  FDRAnalysis <- list(SelectedTrainFeatures=SelectedTrainFeatures,
                      SelectedTestFeatures=SelectedTestFeatures,
                      SelectedLASSOFeatures=SelectedLASSOFeatures,
                      TDR=TDR,
                      ROCAUC=ROCAUC,
                      topROCAUC=topROCAUC,
                      DeSelectedTrainFeatures=DeSelectedTrainFeatures,
                      DeSelectedTestFeatures=DeSelectedTestFeatures,
                      DeSelectedLASSOFeatures=DeSelectedLASSOFeatures,
                      DeTDR=DeTDR,
                      DeROCAUC=DeROCAUC,
                      DetopROCAUC=DetopROCAUC,
                      DeFraction=DeFraction,
                      DeFractionLasso=DeFractionLasso,
                      datadim=datadim,
                      rocTest=rocTest,
                      medTestRaw=medTestRaw,
                      DeMedTestRaw=DeMedTestRaw,
                      testPredict=testPredict,
                      DeTestPredict=DeTestPredict,
                      Avlen=Avlen,
                      NumLatent=NumLatent,
                      totCorrelated=totCorrelated
                      )
  return (FDRAnalysis)
}