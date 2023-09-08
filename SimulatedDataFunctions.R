## Creates a dataframe with a subset of correlated features
## At return:
## data=the synthetic dataframe,
## rotMatrix= the rotation matrix,
## onlyRotatedData= the dataframe with only rotated features

SyntheticData <- function(NumOfObs=1000,NumOfCorFeatures=30,randomToCorrelRatio=2,pOfnoZero=0.35,noiseLevel=0.2,maxLatSize=5)
{
  NRandom <- NumOfCorFeatures
  NObser <- NumOfObs
  FacRand <- randomToCorrelRatio
  
  ## Create a set of random data frame with normal and uniform distributions
  randdata <- matrix(0,nrow=NObser,ncol=NRandom)
  
  for (i in 1:(NRandom/2))
  {
    randdata[,i*2-1] <- rnorm(NObser)
    randdata[,i*2] <- 3*(runif(NObser,0,1) - 0.5)
  }
  
  ## Create the rotation matrix
  rotmat <- matrix(0,nrow=NRandom,ncol=NRandom);
  pupdate <- pOfnoZero
  isLatent <- runif(NRandom) <= sqrt(pupdate)
  isBasis <- (runif(NRandom) <= pupdate) & !isLatent
  latsize <- numeric(NRandom)
#  print(isLatent)
#  print(isBasis)
  for (i in 1:NRandom)
  {
    pin <- pupdate
    if (isBasis[i])  pin <- sqrt(pupdate)       
    if (i <= NRandom/2)
    {
      for (j in 1:(NRandom/2))
      {
        if (isLatent[j] && (latsize[j] < maxLatSize))
        {
          if (runif(1) <= pin)
          {
            rotmat[i,j] <- rnorm(1)
            if (abs(rotmat[i,j]) < 0.2)
            {
              rotmat[i,j] <- 0.2*sign(rotmat[i,j])
            }
            latsize[j] <- latsize[j] + 1;
          }
        }
      }
    }
    else
    {
      for (j in (NRandom/2 + 1):NRandom)
      {
        if (isLatent[j] && (latsize[j] < maxLatSize))
        {
          if (runif(1) <= pin)
          {
            rotmat[i,j] <- rnorm(1)
            if (abs(rotmat[i,j]) < 0.2)
            {
              rotmat[i,j] <- 0.2*sign(rotmat[i,j])
            }
            latsize[j] <- latsize[j] + 1;
          }
        }
      }
    }
  }
  diag(rotmat) <- 1.0
  
  ## Rotate the random matrix
  rotmat <- as.data.frame(rotmat)
  rownames(rotmat) <- paste("V",rownames(rotmat),sep="")
  colnames(rotmat) <- paste("La_",colnames(rotmat),sep="")
  rotmat <- as.matrix(rotmat)
  
  irotmat <- solve(rotmat)

  MCorrelObs <- as.data.frame(randdata %*% irotmat);

  ## Add noise features
  MCorrelObs <- MCorrelObs + matrix(rnorm(NRandom*NObser,0,noiseLevel),nrow=NObser,ncol=NRandom)
  ## Add random features
  CorrelObs <- cbind(MCorrelObs,matrix(rnorm(FacRand*NRandom*NObser),nrow=NObser,ncol=FacRand*NRandom))
  colnames(CorrelObs) <- paste("ER",colnames(CorrelObs),sep="")
  colnames(CorrelObs) <- str_replace_all(colnames(CorrelObs),"ERV","V")

  
  
  result <- list(data=CorrelObs,rotMatrix=rotmat,onlyRotatedData=MCorrelObs)
  return (result)
}

#############################
## Estimates the IDeA rotation and compares with ground truth
## At return:
## FalseDiscovery
## Accuracy
## Correlation on test set
## detectedFeatures
## UPSTM
##rcrit
####################
IDeAEvaluation <- function(traindata,testdata,trueRotation,method=method,corRank=corRank,thr=thr,type=type)
{
  
#  IdeT <- IDeA(traindata,method=method,corRank=corRank,thr=thr,type=type,verbose=TRUE)
  IdeT <- IDeA(traindata,method=method,corRank=corRank,thr=thr,type=type)
  UPSTM <- attr(IdeT,"UPSTM")
  rmat <- UPSTM
  
  if (method=="fast") method="pearson"
  cormat <- abs(cor(IdeT,method=method));
  cormat[is.na(cormat)] <- 0;
  diag(cormat) <- 0;
  trainCorrelation <- max(cormat)
  ideac <- min(attr(IdeT,"IDeAEvolution")$Corr)
#  cat(trainCorrelation,"|")
  if (trainCorrelation > thr)
  {
    if (ideac < (0.95*trainCorrelation))
    {
      mcor <- apply(cormat,2,max)
      cat("\n(",trainCorrelation,":",ideac,")\n");
      print(mcor[mcor==trainCorrelation])
      print(colnames(rmat))
    }
  }
  cormat <- cor(predictDecorrelate(IdeT,testdata),method=method);
  cormat[is.na(cormat)] <- 0;
  diag(cormat) <- 0;
  testCorrelation <- max(abs(cormat))
  falsedisc <- str_detect(rownames(rmat),"R")
  rmat <- rmat[!falsedisc,]
  falsedisc <- str_detect(colnames(rmat),"R")
  rmat <- rmat[,!falsedisc]
  falseDiscovery <- sum(falsedisc)

  rnrotmat <- trueRotation
  colnames(rnrotmat) <- str_remove_all(colnames(rnrotmat),"La_")
  trueUsedFeatures <- rnrotmat != 0
  trueUsedFeatures <- 1*(trueUsedFeatures | t(trueUsedFeatures))
  
  rnrmat <- rmat
  colnames(rnrmat) <- str_remove_all(colnames(rnrmat),"La_")
  detectedFeatures <- rnrmat != 0
  detectedFeatures <- 1*(detectedFeatures | t(detectedFeatures))
  
  
  aux <- detectedFeatures
  detectedFeatures <- 0*trueUsedFeatures
  fullRot <- detectedFeatures
  for (rn in rownames(aux))
  {
    for (cn in colnames(aux))
    {
      detectedFeatures[rn,cn] <- aux[rn,cn];
    }
  }
  diag(detectedFeatures) <- 1
  for (rn in rownames(rnrmat))
  {
    for (cn in colnames(rnrmat))
    {
      fullRot[rn,cn] <- rnrmat[rn,cn];
    }
  }
  diag(fullRot) <- 1
  

  pairAsosTable <- NULL
  for (rn in c(1:(nrow(trueUsedFeatures)-1)))
  {
    for (cn in c((rn+1):ncol(trueUsedFeatures)))
    {
      pairAsosTable <- rbind(pairAsosTable,c(detectedFeatures[rn,cn],trueUsedFeatures[rn,cn]))
    }
  }
  colnames(pairAsosTable) <- c("Estimated","True")
  pairAsosTable <- as.data.frame(pairAsosTable)
  pairTable <- table(pairAsosTable)
#  npairAsosTable <- as.data.frame(pairAsosTable==0)
#  pair_Analysis <- epiR::epi.tests(table(npairAsosTable))
#  sen <- pair_Analysis$detail[3,2]
#  spe <- pair_Analysis$detail[4,2]
#  acc <- pair_Analysis$detail[5,2]
#  pander::pander(c(Accuracy=acc,Sensitivity=sen,Specificity=spe))
  sen <- sum((pairAsosTable$Estimated==pairAsosTable$True)&pairAsosTable$True)/sum(pairAsosTable$True)
  spe <- sum((pairAsosTable$Estimated==pairAsosTable$True)&(pairAsosTable$True==0))/sum(pairAsosTable$True==0)
  acc <- sum(pairAsosTable$Estimated==pairAsosTable$True)/nrow(pairAsosTable)
#  pander::pander(c(Accuracy=acc,Sensitivity=sen,Specificity=spe))
  
  
  
  
  result <- list(trainCorrelation=trainCorrelation,
                 testCorrelation=testCorrelation,
                 falseDiscovery=falseDiscovery,
                 pairTable=pairTable,
                 DiscoveryAccuracy=acc,
                 DiscoverySen=sen,
                 DiscoverySpe=spe,
                 detectedFeatures=detectedFeatures,
                 UPSTM=fullRot,
                 UsedBetas=1*(fullRot != 0),
                 rcrit=attr(IdeT,"R.critical")
                 )
  return(result)
}
  
  