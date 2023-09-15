tvals <- function(x,target)
{
  lmt <- lm(target ~ x)
  sm <- summary(lmt)
  tval <- 0
  if (nrow(sm$coef)>1)
  {
    tval <- sm$coefficients[2,3]
  }
  return(tval)
}

data <- dataframe
CV_IDeA <- function(data,outcome,loops=50,...)
{
  varlist <- colnames(data)
  varlist <- varlist[varlist != outcome]
  
  rawCorrel <- numeric(loops)
  deCorrel <- numeric(loops)
  rawtvalues <- NULL
  detvalues <- NULL
  l <- 1
  for (l in c(1:loops))
  {
    cat(".")
    if (l %% 10 ==0) cat("\n")

    trainsamples <- sample(nrow(data),nrow(data)/2)
    
    trainingset <- data[trainsamples,]
    testingset <- data[-trainsamples,]
    
    
    decormat <- IDeA(trainingset[,varlist],thr=thro)
    
    trainage_DE <- as.data.frame(cbind(decormat,trainingset[,outcome]))
    colnames(trainage_DE) <- c(colnames(decormat),outcome)
    testage_DE <- predictDecorrelate(decormat,testingset)
    
    outcomeModel <- LASSO_1SE(formula(paste(outcome,"~.")),trainingset,...);
#    outcomeModel <- LASSO_1SE(formula(paste(outcome,"~.")),trainingset);
    predOutcome <- predict(outcomeModel,testingset)
    
    outcomeModel_DE <- LASSO_1SE(formula(paste(outcome,"~.")),trainage_DE,...);
#    outcomeModel_DE <- LASSO_1SE(formula(paste(outcome,"~.")),trainage_DE);
    predOutcome_DE <- predict(outcomeModel_DE,testage_DE)

    testOutcome <- testingset[,outcome]
    names(testOutcome) <- rownames(testingset)
    
    ## Testing correlation
    rawCorrel[l] <- cor(predOutcome,testOutcome)
    deCorrel[l] <- cor(predOutcome_DE,testOutcome)
    
    ## Model feature significance
    rawunittvalues <- apply(testingset[,varlist],2,tvals,testOutcome)
    varlistDe <- colnames(testage_DE)[colnames(testage_DE) != outcome]
    tstde <- testage_DE[,varlistDe]
    deunittvalues <- apply(tstde,2,tvals,testOutcome)
    
    lmod <- lm(paste(outcome,"~."),testingset[,c(outcome,names(outcomeModel$coef)[-1])],...)
#    lmod <- lm(paste(outcome,"~."),testingset[,c(outcome,names(outcomeModel$coef)[-1])])
    sm <- summary(lmod)$coef
    coefnames <- rownames(sm)[rownames(sm) %in% names(rawunittvalues)]
    rawtvalues <- rbind(rawtvalues,cbind(rawunittvalues[coefnames],sm[coefnames,3]))
    
    lmod_DE <- lm(paste(outcome,"~."),testage_DE[,c(outcome,names(outcomeModel_DE$coef)[-1])],...)
#    lmod_DE <- lm(paste(outcome,"~."),testage_DE[,c(outcome,names(outcomeModel_DE$coef)[-1])])
    sm <- summary(lmod_DE)$coef
    coefnames <- rownames(sm)[rownames(sm) %in% names(deunittvalues)]
    detvalues <- rbind(detvalues,cbind(deunittvalues[coefnames],sm[coefnames,3]))
  }
  rawtvalues <- as.data.frame(rawtvalues)
  colnames(rawtvalues) <- c("Uni","Model")
  detvalues <- as.data.frame(detvalues)
  colnames(detvalues) <- c("Uni","Model")
  result <- list(testRawCorrelations = rawCorrel,
                 testDeCorrelations = deCorrel,
                 rawtValues = rawtvalues,
                 detValues= detvalues)
  return(result)
}
  