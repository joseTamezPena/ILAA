



CV_TDR <- function(data,outcome,loops=50,scale=TRUE,Decor=TRUE,IDSample=FALSE,
                   MLMethod=LASSO_MIN,MLMethodParams=list(family="binomial"),
                   filterMethod=univariate_Wilcoxon,DecorPCAEFA=TRUE,...)
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
  
  PCASelectedTrainFeatures <- list()
  PCASelectedTestFeatures <- list()
  PCATDR <- numeric(loops)
  PCAROCAUC <-numeric(loops)
  PCAtopROCAUC <- numeric(loops)
    
  EFASelectedTrainFeatures <- list()
  EFASelectedTestFeatures <- list()
  EFATDR <- numeric(loops)
  EFAROCAUC <-numeric(loops)
  EFAtopROCAUC <- numeric(loops)
  
  DeROCAUC <- numeric(loops)
  DetopROCAUC <- numeric(loops)
  DeFraction <- numeric(loops)
  DeFractionLasso <- numeric(loops)
  Avlen <- numeric(loops)
  NumLatent <- numeric(loops)
#  totCorrelated <- numeric(loops)
  rocTest <- list()
  datadim <- NULL
  testPredict <- NULL
  DeTestPredict <- NULL
  PCATestPredict <- NULL
  EFATestPredict <- NULL
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
    ml <- MLMethod(formula(paste(outcome,"~ .")),trainSet,MLMethodParams)
    SelectedLASSOFeatures[[lp]] <- ml$selectedfeatures
    SelectedTrainFeatures[[lp]] <- unique(c(names(filterMethod(trainSet,outcome,pvalue=0.2,limit=-1)),SelectedLASSOFeatures[[lp]]))

    
    mlt <- MLMethod(formula(paste(outcome,"~ .")),testSet,MLMethodParams)
    SelectedTestFeatures[[lp]] <- unique(c(names(filterMethod(testSet,outcome,pvalue=0.2,limit=-1)),mlt$selectedfeatures))
    TDR[lp] <- sum(SelectedTrainFeatures[[lp]] %in% SelectedTestFeatures[[lp]])/length(SelectedTrainFeatures[[lp]])
    
    topROCAUC[lp] <- pROC::roc(testSet[,outcome],testSet[,SelectedTrainFeatures[[lp]][1]],auc=TRUE,quiet=TRUE)$auc
    pred <- as.vector(predict(ml,testSet))
    rocRaw <- pROC::roc(testSet[,outcome],pred,auc=TRUE,quiet=TRUE)
    testPredict <- rbind(testPredict,cbind(rownames(testSet),testSet[,outcome],pred))
    ROCAUC[lp] <- rocRaw$auc
    
    if (Decor)
    {
      ## ILAA decorrelation
      cat("(")
      DEtrainSet <- ILAA(trainSet,...)
#      totCorrelated[lp] <- sqrt(attr(DEtrainSet,"totCorrelated"))

      DEtestSet <- predictDecorrelate(DEtrainSet,testSet)
      cat(")")
      
      ml <- MLMethod(formula(paste(outcome,"~ .")),DEtrainSet,MLMethodParams)
      DeSelectedLASSOFeatures[[lp]] <- ml$selectedfeatures
      DeSelectedTrainFeatures[[lp]] <- unique(c(names(filterMethod(DEtrainSet,outcome,pvalue=0.2,limit=-1)),DeSelectedLASSOFeatures[[lp]]))
      
      mlt <- MLMethod(formula(paste(outcome,"~ .")),DEtestSet,MLMethodParams)
      DeSelectedTestFeatures[[lp]] <- unique(c(names(filterMethod(DEtestSet,outcome,pvalue=0.2,limit=-1)),mlt$selectedfeatures))
      DeTDR[lp] <- sum(DeSelectedTrainFeatures[[lp]] %in% DeSelectedTestFeatures[[lp]])/length(DeSelectedTrainFeatures[[lp]])
      
      DetopROCAUC[lp] <- pROC::roc(DEtestSet[,outcome],DEtestSet[,DeSelectedTrainFeatures[[lp]][1]],auc=TRUE,quiet=TRUE)$auc
      
      pred <- as.vector(predict(ml,DEtestSet))
      rocDe <- pROC::roc(DEtestSet[,outcome],pred,auc=TRUE,quiet=TRUE)
      DeTestPredict <- rbind(DeTestPredict,cbind(rownames(DEtestSet),DEtestSet[,outcome],pred))
      
      DeROCAUC[lp] <- rocDe$auc
      rocTest[[lp]] <- roc.test(rocRaw,rocDe)$p.value
      DeFraction[lp] <- sum(str_detect(DeSelectedTrainFeatures[[lp]],"La_"))/length(DeSelectedTrainFeatures[[lp]])
      if (length(DeSelectedLASSOFeatures[[lp]])>0)
        DeFractionLasso[lp] <- sum(str_detect(DeSelectedLASSOFeatures[[lp]],"La_"))/length(DeSelectedLASSOFeatures[[lp]])
      else
        DeFractionLasso[lp] <- 0;
      ltvar <- getLatentCoefficients(DEtrainSet);
      Avlen[lp] <- mean(sapply(ltvar,length))
      NumLatent[lp] <- length(ltvar)

      if (DecorPCAEFA)
      {
          
   ### Preparation for PCA and FCA
        cat("|")
        iscontinous <- sapply(apply(trainSet,2,unique),length) >= 5 ## Only variables with enough samples
        ndf <- nrow(trainSet)-2
        tvalue <- qt(1.0 - 0.001,ndf)
        rcrit <- tvalue/sqrt(ndf + tvalue^2)
        cmat <- abs(cor(trainSet[,iscontinous]));
        diag(cmat) <- 0;
        maxcor <- apply((cmat > rcrit),2,sum);
        isCorrContinous <- names(maxcor[maxcor > 1])
        idNotInCorCon <- colnames(trainSet)[!(colnames(trainSet) %in% isCorrContinous)]
        cat(round(rcrit,2),",",length(isCorrContinous))
        
   ### PCA Decorrelation
        pc <- prcomp(trainSet[,isCorrContinous],center = TRUE,scale. = TRUE,tol=0.01)   #principal components
        PCApred <- predict(pc,trainSet[,isCorrContinous])
        PCA_Train <- as.data.frame(cbind(PCApred,trainSet[,idNotInCorCon]))
        colnames(PCA_Train) <- c(colnames(PCApred),idNotInCorCon) 
        PCApred <- predict(pc,testSet[,isCorrContinous])
        PCA_Test <- as.data.frame(cbind(PCApred,testSet[,idNotInCorCon]))
        colnames(PCA_Test) <- c(colnames(PCApred),idNotInCorCon) 
        ml <- MLMethod(formula(paste(outcome,"~ .")),PCA_Train,MLMethodParams)
        pred <- as.vector(predict(ml,PCA_Test))
        rocPCA <- pROC::roc(PCA_Test[,outcome],pred,auc=TRUE,quiet=TRUE)
        
        PCAROCAUC[lp] <- rocPCA$auc
        PCATestPredict <- rbind(PCATestPredict,cbind(rownames(PCA_Test),PCA_Test[,outcome],pred))
        PCASelectedTrainFeatures[[lp]] <- unique(c(names(filterMethod(PCA_Train,outcome,pvalue=0.2,limit=-1)),ml$selectedfeatures))
        mlt <- MLMethod(formula(paste(outcome,"~ .")),PCA_Test,MLMethodParams)
        PCASelectedTestFeatures[[lp]] <- unique(c(names(filterMethod(PCA_Test,outcome,pvalue=0.2,limit=-1)),mlt$selectedfeatures))
        PCATDR[lp] <- sum(PCASelectedTrainFeatures[[lp]] %in% PCASelectedTestFeatures[[lp]])/length(PCASelectedTrainFeatures[[lp]])
        
        PCAtopROCAUC[lp] <- pROC::roc(PCA_Test[,outcome],PCA_Test[,PCASelectedTrainFeatures[[lp]][1]],auc=TRUE,quiet=TRUE)$auc
        
        cat("|")
        
   ### EFA Decorrelation
        cat("[")
        if (length(isCorrContinous)<1000)
        {
          topred <- min(length(isCorrContinous),nrow(trainSet),(ncol(PCA_Train)/2+1))
          uls <- try(fa(trainSet[,isCorrContinous],topred,rotate="varimax",warnings=FALSE))  # EFA analysis
          if (!inherits(uls,"try-error"))
          {
            EFApred <- predict(uls,trainSet[,isCorrContinous])
            EFA_Train <- as.data.frame(cbind(EFApred,trainSet[,idNotInCorCon]))
            colnames(EFA_Train) <- c(colnames(EFApred),idNotInCorCon) 
            EFApred <- predict(uls,testSet[,isCorrContinous])
            EFA_Test <- as.data.frame(cbind(EFApred,testSet[,idNotInCorCon]))
            colnames(EFA_Test) <- c(colnames(EFApred),idNotInCorCon) 
            ml <- MLMethod(formula(paste(outcome,"~ .")),EFA_Train,MLMethodParams)
            pred <- as.vector(predict(ml,EFA_Test))
            pred[is.na(pred)] <- 0
            pred[is.infinite(pred)] <- 0
#            print(table(EFA_Test[,outcome]))
            rocEFA <- pROC::roc(EFA_Test[,outcome],pred,auc=TRUE,quiet=TRUE)
            
            EFAROCAUC[lp] <- rocEFA$auc
            EFATestPredict <- rbind(EFATestPredict,cbind(rownames(EFA_Test),EFA_Test[,outcome],pred))
            EFASelectedTrainFeatures[[lp]] <- unique(c(names(filterMethod(EFA_Train,outcome,pvalue=0.2,limit=-1)),ml$selectedfeatures))
            mlt <- MLMethod(formula(paste(outcome,"~ .")),EFA_Test,MLMethodParams)
            EFASelectedTestFeatures[[lp]] <- unique(c(names(filterMethod(EFA_Test,outcome,pvalue=0.2,limit=-1)),mlt$selectedfeatures))
            EFATDR[lp] <- sum(EFASelectedTrainFeatures[[lp]] %in% EFASelectedTestFeatures[[lp]])/length(EFASelectedTrainFeatures[[lp]])
      
            EFAtopROCAUC[lp] <- pROC::roc(EFA_Test[,outcome],EFA_Test[,EFASelectedTrainFeatures[[lp]][1]],auc=TRUE,quiet=TRUE)$auc
          }
          else
          {
            EFAROCAUC[lp] <- 0.5
            EFAtopROCAUC[lp] <- 0.5
            EFATDR[lp] <- 0;
            EFATestPredict <- rbind(EFATestPredict,cbind(rownames(EFA_Test),EFA_Test[,outcome],rep(-1,nrow(EFA_Test))))
            EFASelectedTrainFeatures[[lp]] <- ""
            EFASelectedTestFeatures[[lp]] <- ""
          }
        }
      }
      
      cat("]")
    }
    cat(">")
    if ((lp %% 10) == 0)  cat(lp,"\r\n")
  }
  cat("\n")
  medTestRaw <- boxplot(as.numeric(testPredict[,3])~testPredict[,1],plot=FALSE)
  medTestclass <- boxplot(as.numeric(testPredict[,2])~testPredict[,1],plot=FALSE)
  medTestRaw <- cbind(medTestclass$stats[3,],medTestRaw$stats[3,]) 
  DeMedTestRaw <- boxplot(as.numeric(DeTestPredict[,3])~DeTestPredict[,1],plot=FALSE)
  medTestclass <- boxplot(as.numeric(DeTestPredict[,2])~DeTestPredict[,1],plot=FALSE)
  DeMedTestRaw <- cbind(medTestclass$stats[3,],DeMedTestRaw$stats[3,]) 
  
  PCAMedTestRaw <- NULL
  PCAmedTestclass <- NULL
  EFAMedTestRaw <- NULL
  EFAmedTestclass <- NULL
  if (DecorPCAEFA)
  {
    PCAMedTestRaw <- boxplot(as.numeric(PCATestPredict[,3])~PCATestPredict[,1],plot=FALSE)
    PCAmedTestclass <- boxplot(as.numeric(PCATestPredict[,2])~PCATestPredict[,1],plot=FALSE)
    PCAMedTestRaw <- cbind(PCAmedTestclass$stats[3,],PCAMedTestRaw$stats[3,]) 

    if (length(isCorrContinous)<1000)
    {
      EFAMedTestRaw <- boxplot(as.numeric(EFATestPredict[,3])~EFATestPredict[,1],plot=FALSE)
      EFAmedTestclass <- boxplot(as.numeric(EFATestPredict[,2])~EFATestPredict[,1],plot=FALSE)
      EFAMedTestRaw <- cbind(EFAmedTestclass$stats[3,],EFAMedTestRaw$stats[3,]) 
    }
    else ## Fake numbers
    {
      EFATestPredict <- PCATestPredict
      EFAMedTestRaw <- PCAMedTestRaw
      EFAROCAUC <- PCAROCAUC
      EFASelectedTrainFeatures <- PCASelectedTrainFeatures
      EFASelectedTestFeatures <- PCASelectedTestFeatures
      EFATDR <- PCATDR
      EFAtopROCAUC <- PCAtopROCAUC
    }
  }
  
  FDRAnalysis <- list(SelectedTrainFeatures=SelectedTrainFeatures,
                      SelectedTestFeatures=SelectedTestFeatures,
                      SelectedLASSOFeatures=SelectedLASSOFeatures,
                      TDR=TDR,
                      ROCAUC=ROCAUC,
                      topROCAUC=topROCAUC,
                      testPredict=testPredict,
                      medTestRaw=medTestRaw,
                      
                      DeSelectedTrainFeatures=DeSelectedTrainFeatures,
                      DeSelectedTestFeatures=DeSelectedTestFeatures,
                      DeSelectedLASSOFeatures=DeSelectedLASSOFeatures,
                      DeTestPredict=DeTestPredict,
                      DeMedTestRaw=DeMedTestRaw,
                      
                      DeTDR=DeTDR,
                      DeROCAUC=DeROCAUC,
                      DetopROCAUC=DetopROCAUC,
                      DeFraction=DeFraction,
                      DeFractionLasso=DeFractionLasso,
                      rocTest=rocTest,
                      
                      PCASelectedTrainFeatures=PCASelectedTrainFeatures,
                      PCASelectedTestFeatures=PCASelectedTestFeatures,
                      PCATestPredict=PCATestPredict,
                      PCAMedTestRaw=PCAMedTestRaw,
                      PCATDR=PCATDR,
                      PCAROCAUC=PCAROCAUC,
                      PCAtopROCAUC=PCAtopROCAUC,
                      
                      EFASelectedTrainFeatures=EFASelectedTrainFeatures,
                      EFASelectedTestFeatures=EFASelectedTestFeatures,
                      EFATestPredict=EFATestPredict,
                      EFAMedTestRaw=EFAMedTestRaw,
                      EFATDR=EFATDR,
                      EFAROCAUC=EFAROCAUC,
                      EFAtopROCAUC=EFAtopROCAUC,
                      
                      Avlen=Avlen,
                      NumLatent=NumLatent,
#                      totCorrelated=totCorrelated,
                      datadim=datadim
  )
  return (FDRAnalysis)
}