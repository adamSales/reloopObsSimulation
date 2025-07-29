library(ggplot2)
library(Hmisc)
#library(grDevices)
#################################
#### box-plots for "royal rumble" simulation
################################


oneRun <- function(s,...){
    tot <- reformRun(s)

    boxplot(est~method,data=tot,las=2,mar=c(6.1,4.1,4.1,2.1),...)

}


reformRun <- function(run){

                                        # colnames(run) <- fixNames(colnames(run))
    SD <- if('SD'%in%colnames(run)) run[,'SD'] else run[,ncol(run)]
    colnames(run) <- tolower(colnames(run))

    colnames(run)[which(colnames(run)=='psm')[2]] <- 'psmml'
    colnames(run)[which(colnames(run)=='psm.rebar')[2]] <- 'psmml.rebar'
    if(colnames(run)[1]=='') colnames(run)[1] <- 'diffinmeans'
    colnames(run)[which(colnames(run)=='reloopols')] <- 'psm.reloop'
    colnames(run)[which(colnames(run)=='reloopolsplus')] <- 'psm.reloopPlusOLS'
    colnames(run)[which(colnames(run)=='relooprf')] <- 'psm.reloopPlusRF'



    methods <- c('diffinmeans','psm','psm.rebar',
                 'nn','nn.adj','nn.adj.rebar',
                 'cem','cem.adj','cem.adj.rebar','psmml','psmml.rebar','psm.reloop',
                 'psm.reloopPlusOLS','psm.reloopPlusRF')
    run <- run[,colnames(run)%in%methods]
    for(cc in 1:ncol(run)) run[,cc] <- run[,cc]/SD
    tot <- data.frame(est=do.call('c',as.data.frame(run)),
                      method=rep(colnames(run),each=nrow(run)))

    tot
}

makeNum <- function(run){
    as.data.frame(
        lapply(as.data.frame(run,stringsAsFactors=FALSE),as.numeric))
}

fixNames <- function(nm){
    nm <- tolower(nm)
    nm <- sub('unadj','',nm)
    facs <- c('psm','nn','cem','adj','rebar','mv','r2')
    wfacs <- lapply(facs,function(x) grep(x,nm))
    names(wfacs) <- facs

    new <- character(length(nm))
    for(f in facs)
        new[wfacs[[f]]] <- paste(new[wfacs[[f]]],f,sep='.')
    gsub('^.','',new)
}


plotBadGG <- function(stLasso,stRF){#,sdYc){

    nsim <- nrow(stLasso[[1]])
    ncond <- length(stLasso)

    stTot <- as.data.frame(rbind(do.call('rbind',stLasso),
                   do.call('rbind',stRF)),row.names=seq(2*nsim*ncond))#[,-which(colnames(stRF[[1]])=='MV')]))

    stTot$adjMeth <- c(rep('Lasso',nsim*ncond),rep('RF',nsim*ncond))

    rfCond <- strsplit(names(stRF),'_')
    lasCond <- strsplit(names(stLasso),'_')
    stTot$rho <- c(rep(vapply(lasCond,function(x) x[2],'a'),each=nsim),
                   rep(vapply(rfCond,function(x) x[2],'a'),each=nsim))
    stTot$kappa <- c(rep(vapply(lasCond,function(x) x[1],'a'),each=nsim),
                   rep(vapply(rfCond,function(x) x[1],'a'),each=nsim))

    stTot <- rbind(
        cbind(subset(stTot,select=-psm),Estimator='Rebar'),
        setNames(cbind(subset(stTot,adjMeth!='RF',select=-psm.rebar),'Matching'),
                 c('psm.rebar',names(stTot)[-c(1:2)],'Estimator')))

    stTot$psm.rebar <- -stTot$psm.rebar/stTot$SD

    stTot$adjMeth[stTot$Estimator=='Matching'] <- 'PSM'
    stTot$adjMeth <- factor(stTot$adjMeth,levels=c('PSM','Lasso','RF'))

    stTot$Estimator <- relevel(stTot$Estimator,ref='Matching')

    p <- ggplot(stTot,aes(adjMeth,psm.rebar,fill=Estimator))+geom_boxplot()+
        geom_hline(yintercept=0,lty=2)+ labs(y='Estimate',x='')+
            facet_grid(kappa~rho,labeller=labeller(
                                     kappa=setNames(paste('Unmatched\n Confounding:\n',
                                         c('None','High')),c('0','0.5')),
                                     rho=with(stTot,
                                         setNames(c(paste('$R^2_{\\text{remnant}}(Lasso)\\approx$',
                                                        round(mean(
                                                            R2[adjMeth=='Lasso' & rho=='0'],
                                                            na.rm=TRUE),2),'\n',
                                                        '$R^2_{\\text{remnant}}(RF)\\approx$',
                                                        round( mean(R2[adjMeth=='RF' & rho=='0'],
                                                        na.rm=TRUE),2)),
                                         paste('$R^2_{\\text{remnant}}(Lasso)\\approx$',
                                                        round(mean(
                                                            R2[adjMeth=='Lasso' & rho=='0.05'],
                                                            na.rm=TRUE),2),'\n',
                                                        '$R^2_{\\text{remnant}}(RF)\\approx$',
                                                        round(mean(R2[adjMeth=='RF'&rho=='0.05'],
                                                                   na.rm=TRUE),2))),
                                                  c('0','0.05')))))+
                theme(legend.position='top')
    p+scale_fill_manual(values=c("white","#4e2c89"))+scale_color_manual(values=c("black","#4e2c89"))


}




plotSimGG <- function(st,tikz=FALSE,simple=FALSE){

    if(!is.numeric(st[[1]])){
         st <- lapply(st,makeNum)
    }

    levs <- t(as.data.frame(strsplit(names(st),'_')))


    R2 <- vapply(unique(levs[,2]),function(rh)
        mean(unlist(lapply(st[which(levs[,2]==rh)],function(xx) xx[,'R2'])),na.rm=TRUE),1)

    st <- lapply(st,reformRun)
    dat <- do.call("rbind",st)

    dat$kappa <- rep(levs[,1],each=nrow(st[[1]]))
    dat$rho <- rep(levs[,2],each=nrow(st[[1]]))

    dat$method <- as.character(dat$method)

    if(simple) dat$matchMeth=dat$method else{
        dat$matchMeth <- factor(gsub('\\.[A-Za-z\\.]+','',dat$method))
        levels(dat$matchMeth) <- list(`Optimal Pairs\n Logistic PS`='psm',
                                  `Nearest Neighbors \n Logistic PS`='nn',
                                  `Coursened Exact \n Matches`='cem',
                                  `Optimal pairs \n SuperLearner PS`='psmml')
    }
    dat$Adjustment <- factor(ifelse(grepl('rebar',dat$method),'Rebar',
                             ifelse(grepl('.adj',dat$method),'Within-Sample',
                             ifelse(grepl('reloop',dat$method),'ReLOOP','None'))),
                             levels=c('None','Within-Sample','Rebar','ReLOOP'))

    if(tikz){
     rhoHead <- setNames(paste('$R^2_{\\text{remnant}}\\approx$',round(R2,2)),
                         names(R2))
     } else rhoHead <- setNames(paste("R2=",round(R2*100),"%"),names(R2)) #expression(paste(R[remnant]^2%~~%0,'.',round(R2*100),sep='')),
#                                         names(R2))

    p <- ggplot(dat,aes(matchMeth,est,fill=Adjustment,color=Adjustment))

    p <- p+geom_boxplot(notch=FALSE,outlier.size=0.5) + geom_boxplot(color='black',notch=FALSE,show.legend=FALSE,coef=0,outlier.shape=NA)#outlier.colour=pdat$fill)

    p <- p + geom_hline(yintercept=0,lty=2)+
        labs(y=ifelse(tikz, 'Standardized Bias $\\hat{\\tau}/SD(y_C)$',
                     expression(paste('Standardized Bias',hat(tau)/SD(y[C])))),x='')+
            facet_grid(kappa~rho,labeller=labeller(
                                 kappa=setNames(paste('Unmatched\n Confounding:\n',c('None','Low','High')),
                                                sort(unique(levs[,1]))),
                                     rho=rhoHead))+
        theme(legend.position='top',axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))

    pdat <- ggplot_build(p)$data[[1]]
    aaa <- expand.grid(sort(names(R2)),sort(unique(levs[,1])))
    aaa <- aaa[as.numeric(pdat$PANEL),]
    pdat$rho <- aaa[[1]]
    pdat$kappa <- aaa[[2]]
    p <- p + geom_segment(data=pdat, aes(x=xmin, xend=xmax,
                         y=middle, yend=middle), colour="black", size=1,inherit.aes=FALSE)

    p+scale_fill_manual(values=c("white","#f6b111","#4e2c89"))+scale_color_manual(values=c("black","#f6b111","#4e2c89"))

}



### check for weirdness
ppp <- function(weird){
    par(mfrow=c(3,4))
    numer <- which(!is.na(vapply(weird[1,],as.numeric,1)))
    for(i in numer) plot(as.numeric(weird[,i]),main=colnames(weird)[i])
    }
