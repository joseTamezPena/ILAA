# Testing IDeA
# ILAA is a wrapper of IDeA
library("FRESA.CAD")
data('iris')

colors <- c("red","green","blue")
names(colors) <- names(table(iris$Species))
classcolor <- colors[iris$Species]

#Decorrelating with usupervised basis and correlation goal set to 0.25
system.time(irisDecor <- IDeA(iris,thr=0.25))

## The transformation matrix is stored at "UPLTM" attribute
UPLTM <- attr(irisDecor,"UPLTM")
print(UPLTM)

#Decorrelating with supervised basis and correlation goal set to 0.25
system.time(irisDecorOutcome <- IDeA(iris,Outcome="Species",thr=0.25))
## The transformation matrix is stored at "UPLTM" attribute
UPLTM <- attr(irisDecorOutcome,"UPLTM")
print(UPLTM)

## Compute PCA 
features <- colnames(iris[,sapply(iris,is,"numeric")])
irisPCA <- prcomp(iris[,features]);
## The PCA transformation
print(irisPCA$rotation)

## Plot the transformed sets
plot(iris[,features],col=classcolor,main="Raw IRIS")

plot(as.data.frame(irisPCA$x),col=classcolor,main="PCA IRIS")

featuresDecor <- colnames(irisDecor[,sapply(irisDecor,is,"numeric")])
plot(irisDecor[,featuresDecor],col=classcolor,main="Outcome-Blind IDeA IRIS")


featuresDecor <- colnames(irisDecorOutcome[,sapply(irisDecorOutcome,is,"numeric")])
plot(irisDecorOutcome[,featuresDecor],col=classcolor,main="Outcome-Driven IDeA IRIS")

data('iris')
# Lets model setosa using logistic regression
setosaData <- iris
setosaData$setosa <- 1.0*(as.character(iris$Species)=="setosa")
setosaData$Species <- NULL
setosaILAA <- ILAA(setosaData,thr=0.25)
modelSetosa <- glm(setosa~.,setosaILAA,family="binomial")
#The model cofficients in the ERT space
print(modelSetosa$coefficients)
#Get the observed Coefficients
observedCoef <- getObservedCoef(setosaILAA,modelSetosa)
#The model cofficients in the Observed space
print(observedCoef)



