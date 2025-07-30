source ('simulationFunctions.r')

#packs <- c("MASS","optmatch","SuperLearner","mvtnorm","parallel","randomForest","glmnet","tmle","doMC")
#install.packages(setdiff(packs,installed.packages()),repos="https://mirror.its.umich.edu/cran" )

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
        clusterEvalQ(cl,source ('simulationFunctions.r'))

    } # else{
    #     library(doMC)
    #     registerDoMC(numClus)
    # }
}


res <- justPSMjustBad(500, parr=TRUE)
save(res,file=paste0("output/fullBadResults",Sys.Date(),".RData"))

if(is.element('cl',ls())&inherits(cl,"cluster")) stopCluster(cl)