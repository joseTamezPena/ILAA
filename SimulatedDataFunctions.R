## Creates a dataframe with a subset of correlated features
## At return:
## data=the synthetic dataframe,
## rotMatrix= the rotation matrix,
## onlyRotatedData= the dataframe with only rotated features

SyntheticData <- function(NumOfObs=1000,NumOfCorFeatures=30,randomToCorrelRatio=2,pOfnoZero=0.1,noiseLevel=0.1)
{
  NRandom <- NumOfCorFeatures
  NObser <- NumOfObs
  FacRand <- randomToCorrelRatio
  
  ## Create a set of random data frame with normla and uniform distributions
  randdata <- matrix(0,nrow=NObser,ncol=NRandom)
  
  for (i in 1:(NRandom/2))
  {
    randdata[,i*2-1] <- rnorm(NObser)
    randdata[,i*2] <- 3*(runif(NObser,0,1) - 0.5)
  }
  
  ## Create the rotation matrix
  rotmat <- matrix(0,nrow=NRandom,ncol=NRandom);
  pupdate <- pOfnoZero
  for (i in 1:NRandom)
  {
    if (i <= NRandom/2)
    {
      for (j in 1:(NRandom/2))
      {
        if (runif(1) <= pupdate)
        {
          rotmat[i,j] <- 0.5*rnorm(1)
        }
      }
    }
    else
    {
#      if (i <= 2*NRandom/3)
      {
        for (j in (NRandom/2 + 1):NRandom)
        {
          if (runif(1) <= pupdate)
          {
            rotmat[i,j] <- rnorm(1)
          }
        }
      }
#      else
#      {
#        for (j in (2*NRandom/3 + 1):NRandom)
#        {
#          if (runif(1) < pupdate)
#          {
#            rotmat[i,j] <- rnorm(1)
#          }
#        }
#      }
    }
  }
  rotmat[abs(rotmat) < 0.01] <- 0
  diag(rotmat) <- 0.5 + runif(NRandom)
  
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

## Estimates the IDeA rotation and compares with ground truth
## At return:
## FalseDiscovery
## Accuracy
## Correlation on test set

IDeAEvaluation <- function(traindata,testdata,trueRotation,method=method,corRank=corRank,thr=thr,type=type)
{
  
#  IdeT <- IDeA(traindata,method=method,corRank=corRank,thr=thr,type=type,verbose=TRUE)
  IdeT <- IDeA(traindata,method=method,corRank=corRank,thr=thr,type=type)
  rmat <- attr(IdeT,"UPSTM")
  
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
  GTassVar <- rnrotmat != 0
  GTassVar <- 1*(GTassVar | t(GTassVar))
  
  rnrmat <- rmat
  colnames(rnrmat) <- str_remove_all(colnames(rnrmat),"La_")
  IDassVar <- rnrmat != 0
  IDassVar <- 1*(IDassVar | t(IDassVar))
  
  
  aux <- IDassVar
  IDassVar <- 0*GTassVar
  for (rn in rownames(aux))
  {
    for (cn in colnames(aux))
    {
      IDassVar[rn,cn] <- aux[rn,cn];
    }
  }
  diag(IDassVar) <- 1
  DiscoveryAccuracy <- sum(GTassVar==IDassVar)/(ncol(rnrotmat)^2)
  result <- list(trainCorrelation=trainCorrelation,
                 testCorrelation=testCorrelation,
                 falseDiscovery=falseDiscovery,
                 DiscoveryAccuracy=DiscoveryAccuracy,
                 detectedFeatures=IDassVar
                 )
  return(result)
}
  
  