source('sourceCode/simulationFunctions.r')

library(MASS)
library(optmatch)
library(SuperLearner)
library(Matching)
library(survey)
#library(cem)
library(mvtnorm)
library(parallel)
library(randomForest)
library(glmnet)
library(loop.estimator)
library(tmle)

clustLoad <- FALSE
if(!is.element('cl',ls())) {clustLoad <- TRUE} else if(!'cluster'%in%class(cl)) clustLoad <- TRUE

if(clustLoad){
    if(!is.element('numClus',ls())) numClus <- 10
    if(Sys.info()['sysname']=='Windows'){
        library(snow)
        library(doSNOW)
        cl <- makeCluster(numClus)
        registerDoSNOW(cl)
    } else{
        library(doMC)
        registerDoMC(numClus)
    }
}

smallsim <- function(B=200,X,bg,nt,curved=FALSE, SL.library,SL.libraryZ){
    if(Sys.info()['sysname']=='Windows'){
        print('windows')
        clusterExport(cl,list=c('X','bg','nt','curved','SL.library','SL.libraryZ'),
                      envir=environment())
        clusterExport(cl,list=as.vector(lsf.str(envir=.GlobalEnv)),envir=.GlobalEnv)

        return(
            foreach(i = 1:B,.combine=rbind,
                    .packages=c('optmatch','SuperLearner','tmle',
                                'Matching','survey',#'cem',
                                'glmnet','randomForest',
                                'MASS'))%dopar%{
                                        #cat(i,' ')
                                    simOne(X=X,bg=bg,nt=nt,curved=curved,
                                           SL.library=SL.library,
                                           SL.libraryZ=SL.libraryZ)
                }
        )
    }
    simFun <- function(i) simOne(X=X,bg=bg,nt=nt,curved=curved,SL.library=SL.library,SL.libraryZ=SL.libraryZ)
    res <- mclapply(1:B,simFun,mc.cores=min(parallel::detectCores(),10))
    do.call('rbind',res)
}


smallsimSer <- function(B,X,bg,nt,curved,SL.library,SL.libraryZ){
    simFun <- function(i) simOne(X=X,bg=bg,nt=nt,curved=curved,SL.library=SL.library)
    res <- lapply(1:B,simFun)
    do.call('rbind',res)
}


simTotRR <- function(B,n=400,p=600,nt=50,gm=c(0,0.1,0.5),DECAY=c(0,0.004,0.05),
    SL.library=c('SL.glmnet','SL.randomForest','SL.ridge'),
    SL.libraryZ=setdiff(SL.library,'SL.ridge'),parr=TRUE){

    startTime <- Sys.time()
    runSimCurrent <- readLines('sourceCode/runSim.r')
    simulationFunctionsCurrent <- readLines('sourceCode/simulationFunctions.r')
    CALL <- match.call()

    if(!parr) smallsim <- smallsimSer
    resultsGood <- list()
    X <- lapply(DECAY,function(decay) makeX(ev=exp(-decay*c(1:p))))
    BG <- coefs(gm,n,p)
    tau <- 0 # for now, at least; I don't think this matters

    for(d in 1:length(DECAY))
        for(g in 1:length(gm)){
            cat('gamma=',gm[g],' tau=',tau,' decay=',DECAY[d],'\n')
            s <- smallsim(B=B,X=X[[d]],bg=BG[[g]],nt=nt,curved=FALSE,SL.library=SL.library,SL.libraryZ=SL.libraryZ)
            resultsGood[[paste(gm[g],'_',DECAY[d],sep='')]] <- s
            #save(list=c(as.vector(lsf.str(envir=.GlobalEnv)),'resultsGood','startTime','runSimCurrent','simulationFunctionsCurrent','X','BG','CALL',as.character(CALL)[-c(1,2)]),
             #    file=paste('output/simGood',Sys.Date(),'.RData',sep=''))
        }
    list(resultsGood=resultsGood,BG=BG,X=X)
}

justPSMjustBad <- function(B,n=400,p=600,nt=50,gm=c(0,0.5),DECAY=c(0,0.05),parr=FALSE){

    ## for reproducibility
    startTime <- Sys.time()
    runSimCurrent <- readLines('sourceCode/runSim.r')
    simulationFunctionsCurrent <- readLines('sourceCode/simulationFunctions.r')
    CALL <- match.call()

    if(!parr) smallsim <- smallsimSer
    resultsBadLasso <- list()
    resultsBadRF <- list()

    X <- lapply(DECAY,function(decay) makeX(ev=exp(-decay*c(1:p))))
    BG <- coefs(gm,n,p)


    for(d in 1:length(DECAY))
        for(g in 1:length(gm)){
            cat('lambda=',gm[g],' decay=',DECAY[d],'\n')
            cat('lasso\n')
            print(Sys.time())

            lasso <- psmSim(B=B,X=X[[d]],bg=BG[[g]],nt=nt,curved=TRUE,SL.library=c('SL.glmnet'),mv=FALSE)
            resultsBadLasso[[paste(gm[g],'_',DECAY[d],'_','lasso',sep='')]] <- lasso
            save(list=c(as.vector(lsf.str(envir=.GlobalEnv)),'resultsBadLasso','startTime','runSimCurrent','simulationFunctionsCurrent','X','BG','CALL'),
                 file=paste0('output/simBadLasso',Sys.Date(),'.RData'))

            rf <- psmSim(B=B,X=X[[d]],bg=BG[[g]],nt=nt,curved=TRUE,SL.library=c('SL.randomForest'),mv=FALSE)
            resultsBadRF[[paste(gm[g],'_',DECAY[d],'_','rf',sep='')]] <- rf
            save(list=c(as.vector(lsf.str(envir=.GlobalEnv)),'X','BG','resultsBadRF','startTime','runSimCurrent','simulationFunctionsCurrent','CALL'),
                 file=paste0('output/simBadRF',Sys.Date(),'.RData'))
        }
    list(resultsBadLasso,resultsBadRF)
}



psmSim <- function(B,X,bg,nt,curved,SL.library,mv=FALSE,parr=TRUE){

    if(!parr){
        simFun <- function(i) justPSM(X=X,bg=bg,nt=nt,curved=curved,SL.library=SL.library,mv=mv)
        res <- lapply(1:B,simFun)
        return(do.call('rbind',res))
    }

    if(Sys.info()['sysname']=='Windows'){
        print('windows')
        clusterExport(cl,list=c('X','bg','nt','curved','SL.library','mv'),
                      envir=environment())
        clusterExport(cl,list=as.vector(lsf.str(envir=.GlobalEnv)),envir=.GlobalEnv)

        return(
            foreach(i = 1:B,.combine=rbind,
                    .packages=c('optmatch','SuperLearner',
                                'glmnet','randomForest',
                                'MASS'))%dopar%{
                                    justPSM(X=X,bg=bg,nt=nt,curved=curved,SL.library=SL.library,mv=mv)
                                }
        )
    }

    simFun <- function(i) simOne(X=X,bg=bg,nt=nt,curved=curved,SL.library=SL.library,SL.libraryZ=SL.libraryZ)
    res <- mclapply(1:B,simFun,mc.cores=min(parallel::detectCores(),10))
    do.call('rbind',res)
}


runSim <- function(){
    st <- simTotRR(B=1000)
    save(list=c(as.vector(lsf.str(envir=.GlobalEnv)),'st'),
         file=paste0('output/simResultsGood',Sys.Date()))
    bad <- justPSMjustBad(B=1000)
    save(list=c(as.vector(lsf.str(envir=.GlobalEnv)),'bad'),
         file=paste0('output/simResultsBad',Sys.Date()))
    list(st,bad)
}

