### estimation functions

MLprop <- function(dat,X,SL.libraryZ){
    mod <- SuperLearner(dat$z,as.data.frame(X),family=binomial,SL.library=SL.libraryZ)
    mod$SL.predict

}
MLest <- function(dat,X,SL.library,SL.libraryZ){
    pscores <- MLprop(dat,X,SL.libraryZ)
    psmEst(dat,X,pscores,SL.library,runmv=FALSE)
}

TMLE <- function(dat,X,SL.library,pscores){
  mod=try(tmle(
    Y=dat$y,
    A=dat$z,
    W=as.data.frame(X),
    Q.SL.library=SL.library,
    g1W=pscores))
 if(inherits(mod,"try-error")){
  print("NO TMLE THIS TIME")
  return(NA)
 }
 mod$estimates$ATE$psi
}

### nearest neighbor matching with or without bias adjustment, rebar
NNmatch <- function(propscore,dat,X,SL.library){
  match <- Match(dat$y,dat$z,propscore,Z=X[,1:5],M=1,ties=TRUE,BiasAdjust = TRUE,estimand = 'ATT',replace = TRUE)
  est.adj <- match$est
  match <- Match(dat$y,dat$z,propscore,Z=X[,1:5],M=1,ties=TRUE,BiasAdjust = FALSE,estimand = 'ATT',replace = TRUE)
  est.unadj <- match$est

  ### prog score adjustment to both

  matched <- c(match$index.treated,match$index.control)
  m <- rep(NA,nrow(X))
  m[matched] <- 1
  progs <- progEst(dat,X,m,returnMod=TRUE,SL.library=SL.library)$SL.predict

  NNrebar <- Match(dat$y-progs,dat$z,propscore,Z=X[,1:5],M=1,ties=TRUE,
                   BiasAdjust = FALSE,estimand = 'ATT',replace = TRUE)$est
  NNairebar <- Match(dat$y-progs,dat$z,propscore,Z=X[,1:5],M=1,ties=TRUE,
                     BiasAdjust = TRUE,estimand = 'ATT',replace = TRUE)$est
  setNames(c(NNunadj=est.unadj,NNadj=est.adj,NNrebar=NNrebar,NNairebar=NNairebar),
           c('NN','NN.adj','NN.rebar','NN.adj.rebar'))
}

cemEst <- function(dat,X,numcut=5,SL.library){
  cemdat <- data.frame(z=dat$z,X[,1:5])
  match <- try(cem('z',cemdat,cutpoints=list(X1=numcut,X2=numcut,X3=numcut,X4=numcut,X5=numcut))$mstrata)

  progs <- progEst(dat,X,match,returnMod=TRUE,SL.library=SL.library)$SL.predict
  e <- dat$y-progs

  setNames(c(matchEst(outcome=dat$y,z=dat$z,m=match),
           coef(lm(y~dat$z+X[,1:5],data=dat,subset=!is.na(match)))[2],
           coef(lm(e~dat$z+X[,1:5],data=dat,subset=!is.na(match)))[2],
           matchEst(outcome=e,z=dat$z,m=match)),
           c('CEM','CEM.adj','CEM.adj.rebar','CEM.rebar'))
}


progEst <- function(dat,X,m,ps,returnMod=FALSE,SL.library){

  if(!missing(ps)) X <- model.matrix(~X*ps)
  progDat <- dat[dat$z==0 & is.na(m),]
  progX <- X[dat$z==0 & is.na(m),]

  progMod <- SuperLearner(progDat$y,as.data.frame(progX),as.data.frame(X),
                             SL.library=SL.library)

  dat$progScores <- progMod$SL.predict

  est <-  lm(y-progScores~z+m,data=dat)$coef[2]
  if(returnMod){
      return(progMod)
  }
  c(est,R2=progMod$cvRisk)
}

psmEst <- function(dat,X,pscores,SL.library,runmv=FALSE){

    match <- PSmatch(dat,pscores)
    est <- matchEst(dat$y,dat$z,match)
    #progMod <- progEst(dat=dat,X=X,m=match,returnMod=TRUE,SL.library=SL.library)
    #prog <- progMod$SL.predict
    
    #rf=randomForest(X[dat$z==0&is.na(match),],dat$y[dat$z==0&is.na(match)])
    #prog=predict(rf,X)
    remnantMod <- SuperLearner(
      Y=dat$y[dat$z==0&is.na(match)],
      X=as.data.frame(X[dat$z==0&is.na(match),]),
      newX=as.data.frame(X),
      SL.library=SL.library)
    
    prog <- remnantMod$SL.predict

    e <- dat$y-prog
    rebar <- matchEst(e,dat$z,match)

    reloopOLS <- p_loop(Y=dat$y[!is.na(match)],Tr=dat$z[!is.na(match)],
                        Z=cbind(prog[!is.na(match)]),
                        P=as.numeric(match[!is.na(match)]),
                        pred=p_ols_interp)
    reloopOLSplus <- p_loop(Y=dat$y[!is.na(match)],Tr=dat$z[!is.na(match)],
                        Z=cbind(dat$x[!is.na(match),1:5],prog[!is.na(match)]),
                        P=as.numeric(match[!is.na(match)]),
                        pred=p_ols_interp)
    reloopRF <- p_loop(Y=dat$y[!is.na(match)],Tr=dat$z[!is.na(match)],
                        Z=cbind(dat$x[!is.na(match),1:5],prog[!is.na(match)]),
                        P=as.numeric(match[!is.na(match)]),
                        pred=p_rf_interp)


    mv <- NA
    if(runmv) mv <- MV(dat,X,match,pscores,SL.library)
    R2 <- rf$rsq[length(rf$rsq)] #1-min(progMod$cvRisk)/var(dat$y[is.na(match)])

    setNames(c(est,rebar,
              reloopOLS[1],
              reloopOLSplus[1],
              reloopRF[1],
              mv,R2), c('psm','psm.rebar','reloopOLS','reloopOLSplus','reloopRF','MV','R2'))
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

invLogit <- function(lp){
  odds <- exp(lp)
  odds/(1+odds)
}


makeDataCurved <- function(X,bg,nt){
  n <- nrow(X)

  linPred <- crossprod(t(X),bg[,'gamma'])

  cc <- sort(linPred)
  thresh <- sort(cc,decreasing=TRUE)[2*nt]
  matchedZ <- sample(rep(c(0,1),nt))
  Z <- rep(0,n)
  matched <- linPred>=thresh
  Z[matched] <- matchedZ

  Yclin <- crossprod(t(X),bg[,'beta'])
  Yclin[matched] <- mean(Yclin)-crossprod(t(X),bg[,'beta'])[matched]
  Y <- Yclin + rnorm(n)
  data.frame(y=Y,z=Z)

}

makeData <- function(X,bg,nt){
  n <- nrow(X)

  linPred <- crossprod(t(X),bg[,'gamma'])

  fakeOut <- c(rep(1,nt),rep(0,n-nt))
  mmm <- glm(fakeOut~1,family=binomial,offset=linPred)
  ps <- predict(mmm,type='response')

  Z <- rbinom(n,1,ps)
   Yclin <- crossprod(t(X),bg[,'beta'])
  Y <- Yclin + rnorm(n)
  data.frame(y=Y,z=Z)

}


simOne <- function(X,bg,nt,curved,SL.library,SL.libraryZ){
    if(curved) dat <- makeDataCurved(X,bg,nt)
    else dat <- makeData(X,bg,nt)

    pmod <- propmodel(dat,X)

    pscoreML <- MLprop(dat,X,SL.libraryZ)

    if(curved){
        tauhats <- c(
            mean(dat$y[dat$z==1])-mean(dat$y[dat$z==0]), #unmatched
            psmEst(dat,X,pmod$linear,SL.library), #PSM
            )
        } else{
            tauhats <- c(
                unadj=mean(dat$y[dat$z==1])-mean(dat$y[dat$z==0]), #unmatched
                psm=psmEst(dat,X,pmod$linear,SL.library),
                psmML=psmEst(dat,X,pscoreML,SL.library),
                tmleEst=TMLE(dat,X,SL.library,pscoreML)
                )#, #PSM
                #NNmatch(pmod$linear.predictors,dat,X,SL.library=SL.library))#, #NN
                #cemEst(dat,X,SL.library=SL.library), #CEM
                #MLest(dat,X,SL.library,SL.libraryZ)) ## ML + PSM
        }

    c(tauhats,sd(dat$y))
}

simOnePsmML <- function(X,bg,nt,curved,SL.library,SL.libraryZ){
    if(curved) dat <- makeDataCurved(X,bg,nt)
    else dat <- makeData(X,bg,nt)

    pmod <- propmodel(dat,X)

    tauhats <- c(#psmEst(dat,X,pmod$linear,SL.library), #PSM
                 #NNmatch(pmod$linear.predictors,dat,X,SL.library=SL.library), #NN
                 #cemEst(dat,X,SL.library=SL.library), #CEM
                 MLest(dat,X,SL.library,SL.libraryZ)) ## ML + PSM


    c(tauhats,sd(dat$y))
}

simOneUnmatched <- function(X,bg,nt,curved){
    if(curved) dat <- makeDataCurved(X,bg,nt)
    else dat <- makeData(X,bg,nt)
    c(mean(dat$y[dat$z==1])-mean(dat$y[dat$z==0]),sd(dat$y))
}

plotCurved <- function(gm=0,decay=0.05,curved=TRUE){
    X <- makeX(decay=decay)
    bg <- coefs(gm,400,600)[[1]]

    if(curved) dat <- makeDataCurved(X,bg)
    else dat <- makeData(X,bg)

    pmod <- propmodel(dat,X)
    m <- PSmatch(dat,X,pmod)

    yhat <- progEst(dat,X,m,returnMod=TRUE)$SL.predict

    xyplot(yhat~dat$y|matched(m),groups=dat$z,auto.key=TRUE)
}


MV <- function(dat,X,m,pscores,SL.library){
    mBig <- fullmatch(dat$z~pscores,max=5,data=dat)
    train <- !matched(m)
    mBig <- mBig[train]
    yt <- dat$y[train]
    Xt <- as.data.frame(X[train,])
    mvMod <- SuperLearner(yt[is.na(mBig)],Xt[is.na(mBig),],Xt,SL.library=SL.library)
    MV <- mean((yt[matched(mBig)]-mvMod$SL.predict[matched(mBig)])^2)/min(mvMod$cvRisk)
    MV
}



recoverSim <- function(){
    lll <- list.files()[grep('sim',list.files())]

    st <- list()


    for(l in lll){
        print(l)
        params <- strsplit(l,c('sim|_|FALSE|TRUE|.RData'))[[1]]
        print(params)
        load(l)
        st[[paste(params,sep='_')]] <- s
    }
    st
}


ntFun <- function(i){
    dat <- makeData(X,bg)
    sum(dat$z)
}


SL.ranger <- function(Y,X,newX,family,
                      mtry = ifelse(family$family == "gaussian",
                          floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)),
                      ntree = 1000, nodesize = ifelse(family$family == "gaussian", 5,1),...){
    #.SL.require('ranger')
    dat <- as.data.frame(X)
    dat$Y <- Y
    newDat <- as.data.frame(newX)

    fit.rf <- ranger::ranger(Y~.,data=dat,num.trees=ntree,mtry=mtry,min.node.size=nodesize)
    pred <- predict(fit.rf,newDat)$predictions
    fit <- list(object=fit.rf)
    out <- list(pred=pred,fit=fit)
    class(out$fit) <- c("SL.ranger")
    return(out)
}


justPSM <- function(X,bg,nt,curved,SL.library,mv=TRUE){
    if(curved) dat <- makeDataCurved(X,bg,nt)
    else dat <- makeData(X,bg,nt)

    pmod <- propmodel(dat,X)

    c(psmEst(dat=dat,X=X,pscores=pmod$linear,SL.library=SL.library,runmv=mv),SD=sd(dat$y))
}

psmEstNoMV <- function(dat,X,pmod,SL.library){

    match <- PSmatch(dat,X,pmod)
    est <- matchEst(dat$y,dat$z,match)
    progMod <- progEst(dat=dat,X=X,m=match,returnMod=TRUE,SL.library=SL.library)
    e <- dat$y-progMod$SL.predict
    rebar <- matchEst(e,dat$z,match)

    R2 <- 1-min(progMod$cvRisk)/var(dat$y[is.na(match)])

    setNames(c(est,rebar,R2), c('psm','psm.rebar','R2'))
}
