### estimation functions
rfRemMod <- function(X,dat,match){
  remnantMod <- randomForest(x=X[dat$z==0&is.na(match),],y=dat$y[dat$z==0&is.na(match)],
        xtest=X)
  remnantMod$test$predicted
}

psmEst <- function(dat,X,pscores,progMod=rfRemMod,...){ #,SL.library,runmv=FALSE){

    match <- PSmatch(dat,pscores)
    est <- matchEst(dat$y,dat$z,match)

    prog <- progMod(X,dat,match)

    e <- dat$y-prog

    raw <- mean(dat$y[dat$z==1])-mean(dat$y[dat$z==0])

    rebar <- matchEst(e,dat$z,match)

    reloopOLS <- p_loop(Y=dat$y[!is.na(match)],Tr=dat$z[!is.na(match)],
                        Z=cbind(prog[!is.na(match)]),
                        P=as.numeric(match[!is.na(match)]),
                        pred=p_ols_interp)
    reloopOLSv12 <- p_loop(Y=dat$y[!is.na(match)],Tr=dat$z[!is.na(match)],
                           Z=cbind(prog[!is.na(match)]),
                           P=as.numeric(match[!is.na(match)]),
                           pred=p_ols_v12)

    reloopOLSpo <- p_loop(Y=dat$y[!is.na(match)],Tr=dat$z[!is.na(match)],
                           Z=cbind(prog[!is.na(match)]),
                           P=as.numeric(match[!is.na(match)]),
                           pred=p_ols_po)

    # reloopOLSplus <- p_loop(Y=dat$y[!is.na(match)],Tr=dat$z[!is.na(match)],
    #                     Z=cbind(dat$x[!is.na(match),1:5],prog[!is.na(match)]),
    #                     P=as.numeric(match[!is.na(match)]),
    #                     pred=p_ols_interp)
    reloopRF <- p_loop(Y=dat$y[!is.na(match)],Tr=dat$z[!is.na(match)],
                        Z=cbind(dat$x[!is.na(match),1:5],prog[!is.na(match)]),
                        P=as.numeric(match[!is.na(match)]),
                        pred=p_rf_interp)


    mv <- NA
    #if(runmv) mv <- MV(dat,X,match,pscores,SL.library)
    R2 <- 0 #rf$rsq[length(rf$rsq)] #1-min(progMod$cvRisk)/var(dat$y[is.na(match)])

    setNames(c(raw,est,rebar,
              reloopOLS[1],
              reloopOLSplus[1],
              reloopRF[1],sd(dat$y),
              mv,R2), c("raw",'psm','psm.rebar','reloopOLS','reloopOLSplus','reloopRF',"sdY",'MV','R2'))
}

PSmatch <- function(dat,pscores){
  if(all(pscores>0 & pscores<1)) pscores <- qlogis(pscores)
  m <- pairmatch(dat$z~pscores,remove=TRUE,data=dat)
  m
}

matchEst <- function(outcome,z,m)
  lm(outcome~z+m)$coef[2]


propmodel <- function(dat,X){
  glm(z~X[,1:5],data=dat,family=binomial(logit))
}


#### data generating functions

# from https://stat.ethz.ch/pipermail/r-help/2008-February/153708.html
Posdef <- function (n, ev = runif(n, 0, 10))
{
  Z <- matrix(ncol=n,rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

makeX <- function(n=400,p=600,decay=0.1,ev=exp(-decay*c(1:p))){
  sig <- Posdef(p,ev)
  #diag(sig) <- 1
  X <- rmvnorm(n,rep(0,p),sig)
  X <- scale(X)
  X
}


coefs <- function(gm,n,p){
  beta <- c(rep(1,5),rexp(p-5,5))
  BG <- lapply(gm, function(gammaRate){
    gamma <- beta
    gamma[6:p] <- gammaRate*beta[6:p]
    cbind(beta,gamma)
  })
  BG
}


makeDataCurved <- function(X,bg,nt,justTrt=FALSE){
  n <- nrow(X)

  linPred <- X%*%bg[,"gamma"]

  cc <- sort(linPred)
  thresh <- sort(cc,decreasing=TRUE)[2*nt]
  matchedZ <- sample(rep(c(0,1),nt))
  Z <- rep(0,n)
  matched <- linPred>=thresh
  Z[matched] <- matchedZ

  Yclin <- crossprod(t(X),bg[,'beta'])
  Yclin[matched] <- mean(Yclin)-crossprod(t(X),bg[,'beta'])[matched] ## orig version
      #           2*mean(Yclin[matched])-Yclin[matched]#crossprod(t(X),bg[,'beta'])[matched]
  Y <- Yclin + rnorm(n)
  out <- data.frame(y=Y,z=Z)
  attr(out,"matched") <- matched
  #cat("e")
  out
}


makeDataCurved2 <- function(X,bg,nt,justTrt=FALSE){
  n <- nrow(X)
   #cat("c")
  linPred <- rowSums(X[,1:5])
      #X%*%bg[,"gamma"] #crossprod(t(X),bg[,'gamma'])
   #cat("d")
  cc <- sort(linPred)
  thresh <- sort(cc,decreasing=TRUE)[2*nt]
  matchedZ <- sample(rep(c(0,1),nt))
  Z <- rep(0,n)
  matched <- linPred>=thresh
  Z[matched] <- matchedZ

  Yclin <- crossprod(t(X),bg[,'beta'])
  Yclin[matched] <- 2*mean(Yclin[matched])-Yclin[matched]#crossprod(t(X),bg[,'beta'])[matched]
  Y <- Yclin + rnorm(n)
  out <- data.frame(y=Y,z=Z)
  attr(out,"matched") <- matched
  #cat("e")
  out
}


makeData <- function(X,bg,nt,trtCurve=FALSE,curveFun=function(x,y) x){
  n <- nrow(X)

  linPred <- crossprod(t(X),bg[,'gamma'])

  fakeOut <- c(rep(1,nt),rep(0,n-nt))
  mmm <- glm(fakeOut~1,family=binomial,offset=linPred)
  ps <- predict(mmm,type='response')

  Z <- rbinom(n,1,ps)
  Yclin <- crossprod(t(X),bg[,'beta'])
  if(trtCurve) Yclin[Z==1] <- 2*mean(Yclin[Z==1])-Yclin[Z==1]
  Y <- curveFun(Yclin,linPred) + rnorm(n)
  data.frame(y=Y,z=Z)

}


justPSM <- function(X,bg,nt,curved,
                    #SL.library,mv=FALSE,
                    trtCurve=FALSE,
                    curveFun=\(x,y) x){
    stopifnot(!(curved&trtCurve))
    #cat("a")
    if(curved) dat <- makeDataCurved(X,bg,nt)
    else dat <- makeData(X,bg,nt,trtCurve = trtCurve,curveFun = curveFun)
    #cat("b")
    pmod <- propmodel(dat,X)
    #cat(".")
    c(psmEst(dat=dat,X=X,pscores=pmod$linear),SD=sd(dat$y))
}


###########################################
## functions to run simulation
###########################################


psmSim <- function(B,X,bg,nt,curved,trtCurve=FALSE,curveFun=\(x,y) x, parr=TRUE,...){

    simFun <- function(i){
            #cat(i," ")
            justPSM(X=X,bg=bg,nt=nt,curved=curved,trtCurve=trtCurve,curveFun=curveFun)
        }

    if(!parr){
        res <- lapply(1:B,simFun)
        return(do.call('rbind',res))
    }

    if(sock|Sys.info()['sysname']=='Windows'){
        print('windows')
        pb <- txtProgressBar(max = B, style = 3)
        progress_fun <- function(nn) setTxtProgressBar(pb, nn)
        opts <- list(progress = progress_fun)

        #clusterExport(cl,list=as.list(match.call())[-1],
        #              envir=environment())
        #clusterExport(cl,list=as.vector(lsf.str(envir=.GlobalEnv)),envir=.GlobalEnv)
        clusterExport(cl,list=c("X","bg","nt","curved","trtCurve"),envir=environment())
        return(
            foreach(i = 1:B,.combine=rbind,
                    .packages=c('optmatch','SuperLearner',
                                'glmnet','randomForest',"dRCT",
                                'MASS'),
                                .options.snow =opts)%dopar%{
                                simFun(i)
                                }
        )
    }

    res <- mclapply(1:B,simFun,mc.cores=5)
    do.call('rbind',res)
}



justPSMsim <- function(B,n=400,p=600,nt=50,gm=c(0,0.5),DECAY=c(0,0.05),curved=TRUE,trtCurve=FALSE,curveFun=\(x,y) x,parr=FALSE){

    ## for reproducibility
    startTime <- Sys.time()
    runSimCurrent <- readLines('runSim.r')
    simulationFunctionsCurrent <- readLines('simulationFunctions.r')
    CALL <- match.call()

    #if(!parr) smallsim <- smallsimSer
    #resultsBadLasso <- list()
    resultsBadRF <- list()

    X <- lapply(DECAY,function(decay) makeX(ev=exp(-decay*c(1:p))))
    BG <- coefs(gm,n,p)

    for(d in 1:length(DECAY))
        for(g in 1:length(gm)){
            cat('lambda=',gm[g],' decay=',DECAY[d],'\n')
            #cat('lasso\n')
            print(Sys.time())

            rf <- psmSim(B=B,X=X[[d]],bg=BG[[g]],nt=nt,curved=curved,trtCurve = trtCurve, curveFun = curveFun,parr=parr)
            resultsBadRF[[paste(gm[g],'_',DECAY[d],'_','rf',sep='')]] <- rf
            save(list=c(as.vector(lsf.str(envir=.GlobalEnv)),'X','BG','resultsBadRF','startTime','runSimCurrent','simulationFunctionsCurrent','CALL'),
                 file=paste0('output/simBadRF',Sys.Date(),'.RData'))
        }
    resultsBadRF
}

justPSMtrtBad <- function(B,n=400,p=600,nt=50,gm=c(0,0.5),DECAY=c(0,0.05),parr=FALSE){

    ## for reproducibility
    startTime <- Sys.time()
    runSimCurrent <- readLines('runSim.r')
    simulationFunctionsCurrent <- readLines('simulationFunctions.r')
    CALL <- match.call()

    #if(!parr) smallsim <- smallsimSer
    resultsBadLasso <- list()
    resultsBadRF <- list()

    X <- lapply(DECAY,function(decay) makeX(ev=exp(-decay*c(1:p))))
    BG <- coefs(gm,n,p)


    for(d in 1:length(DECAY))
        for(g in 1:length(gm)){
            cat('lambda=',gm[g],' decay=',DECAY[d],'\n')
            #cat('lasso\n')
            print(Sys.time())

            rf <- psmSim(B=B,X=X[[d]],bg=BG[[g]],nt=nt,curved=FALSE,trtCurve=TRUE,parr=parr)
            resultsBadRF[[paste(gm[g],'_',DECAY[d],'_','rf',sep='')]] <- rf
            save(list=c(as.vector(lsf.str(envir=.GlobalEnv)),'X','BG','resultsBadRF','startTime','runSimCurrent','simulationFunctionsCurrent','CALL'),
                 file=paste0('output/simTrtCurve',Sys.Date(),'.RData'))
        }
    list(resultsBadLasso,resultsBadRF)
}

