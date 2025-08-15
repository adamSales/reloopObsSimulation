source ('simulationFunctions.r')

packs <- c("MASS","optmatch","SuperLearner","snow","mvtnorm","remotes","parallel","randomForest","glmnet","tmle","doMC")
install.packages(setdiff(packs,installed.packages()),repos="https://cran.rstudio.com/")

if(!require("dRCT")) remotes::install_github("manncz/dRCT")

library(MASS)
library(optmatch)
library(SuperLearner)
#library(Matching)
#library(survey)
#library(cem)
library(mvtnorm)
library(parallel)
library(randomForest)
library(glmnet)
library(dRCT)
library(tmle)
library(utils)

sock <- TRUE
clustLoad <- FALSE
if(!is.element('cl',ls())) {clustLoad <- TRUE} else if(!'cluster'%in%class(cl)) clustLoad <- TRUE

if(clustLoad){
    if(!is.element('numClus',ls())) numClus <- 10
    if(sock|Sys.info()['sysname']=='Windows'){
        library(snow)
        library(doSNOW)
        cl <- makeCluster(numClus)
        registerDoSNOW(cl)
        clusterEvalQ(cl,source('simulationFunctions.r'))

    } # else{
    #     library(doMC)
    #     registerDoMC(numClus)
    # }
}

cf1 <- function(x,y){
    a <- quantile(y,0.025)
    b <- coef(lm(y~x))[2]
    x*(y-2*a)/abs(2*a*b)
}

cf2 <- function(x,y){
    qq <- quantile(x,0.025)
    x*(x-2*qq)/abs(2*qq)
}

cf3 <- function(x,y){
    q5 <- median(x)
    q7 <- quantile(x,0.7)
    (x-q7)^2/(2*(q5-q7))
}

resLinear <- justPSMsim(500,curved=FALSE,curveFun=\(x,y) x,parr=TRUE)
save(list=ls(),file=paste0("output/linearResults",Sys.Date(),".RData"))

resCurved <- justPSMsim(500,curved = TRUE,parr=TRUE)
save(list=setdiff(ls(),"resLinear"),
   file=paste0("output/curvedResults",Sys.Date(),".RData"))

resCF1 <- justPSMsim(500, curved=FALSE,curveFun=cf1,parr=TRUE)
save(resCF1,file=paste0("output/cf1results",Sys.Date(),".RData"))

resCF2 <- justPSMsim(500, curved=FALSE,curveFun=cf2,parr=TRUE)
save(resCF2,file=paste0("output/cf2results",Sys.Date(),".RData"))

resCF3 <- justPSMsim(500,curved=FALSE,curveFun=cf3,parr=TRUE)
save(resCF3,file=paste0("output/cf3results",Sys.Date(),".RData"))

stopCluster(cl)
