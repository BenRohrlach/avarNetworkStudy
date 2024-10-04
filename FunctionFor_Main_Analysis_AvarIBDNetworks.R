# Functions for network analysis of Avar Cemetries 
# Created by: AB Rohrlach
# Created on: 11-05-2023


#### packages and custom functions ----
if(!require("pacman")) install.packages("pacman")
pacman::p_load(ggmap,
               tidyverse,
               ggplot2,
               spaa,
               network,
               ggnetwork,
               ergm,
               GERGM,
               philentropy,
               stringi)

df2dist <- function(idi,idj,dij,default=Inf){
  if(length(idi)==length(idj)&length(idi)==length(dij)){
    ids <- unique(c(idi,idj))
    n <- length(ids)
  }else{
    stop('Input vectors are not the same length.')
  }
  
  Dij <- matrix(rep(default,n^2),nrow=n,dimnames=list(ids,ids))
  for(k in 1:length(idi)){
    IDi <- idi[k]
    rowi <- which(ids==IDi)
    IDj <- idj[k]
    colj <- which(ids==IDj)
    Dij[rowi,colj] <- Dij[colj,rowi] <- dij[k]
  }
  for(k in 1:n){
    Dij[k,k] <- 0
  }
  return(Dij)
}

# makeZero <- function(x,v){
#   return(x*v)
# }

# maxNorm <- function(x){
#   y <- max(x)-x
#   return(y)
# }

bestModelErgm <- function(ergm.list){
  aic_val <- ergm.list %>% 
    lapply(AIC) %>%
    unlist()
  aic_inf <- ergm.list %>%
    lapply(coef) %>%
    lapply(is.infinite) %>%
    lapply(any) %>%
    unlist()
  tibble(aic_inf=aic_inf,
         aic_val=aic_val) %>%
    dplyr::mutate(n=1:n()) %>%
    dplyr::filter(!aic_inf) %>%
    dplyr::filter(aic_val==min(aic_val)) %>%
    dplyr::pull(n) %>%
    return()
}

distRules <- function(D,R){
  n <- nrow(D)
  Y <- matrix(rep(0,n^2),nrow=n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      Y[i,j] <- Y[j,i] <- ifelse(R[i,j],D[i,j],0)
    }
  }
  Y
  return(Y)
}

ERGMS.aicPlot <- function(ergm.list,form_v,ang=35,plotmax=20,tsize=12,printRes=F){
  aic_inf <- lapply(ergm.list,coef) %>%
    lapply(is.infinite) %>%
    lapply(any) %>%
    unlist() %>%
    ifelse(Inf,1)
  aic_v <- lapply(ergm.list, AIC) %>% unlist()
  coefs_v <- lapply(ergm.list, function(x){ return(length(x$coefficients))}) %>% unlist()
  aic_tib <- tibble(AIC=aic_v*aic_inf,n=coefs_v,call=form_v) %>%
    dplyr::filter(!is.infinite(AIC)) %>%
    dplyr::arrange(AIC) %>%
    dplyr::slice(1:plotmax) %>%
    dplyr::arrange(n) %>%
    dplyr::mutate(xrel=1:n())
  
  rect_tib <- aic_tib %>%
    dplyr::group_by(n) %>%
    dplyr::summarise(xmin=min(xrel)-0.5,xmax=max(xrel)+0.5)
  
  outplot <- aic_tib %>%
    dplyr::mutate(best=AIC==min(AIC)) %>%
    ggplot(aes(x=xrel,y=AIC,group=1))+
    theme_bw()+
    geom_rect(inherit.aes=F,
              data=rect_tib,
              aes(ymin=-Inf,ymax=Inf,xmin=xmin,xmax=xmax,alpha=(n)))+
    geom_line(linetype='dashed')+
    geom_point(pch=19,size=6)+
    geom_text(aes(label=n),col='white')+
    geom_point(data=aic_tib%>%dplyr::filter(AIC==min(AIC)),
               pch=19,size=5,col='red')+
    geom_text(data=aic_tib%>%dplyr::filter(AIC==min(AIC)),
              aes(label=n),col='black')+
    scale_x_continuous(breaks=aic_tib$xrel,
                       labels=aic_tib$call,
                       expand=c(0,0))+
    scale_alpha_continuous(range=c(0.1,0.5),
                           name='# coeffcients')+
    theme(axis.text.x=element_text(angle=ang,hjust=1),
          legend.position='none',
          axis.text=element_text(size=tsize))+
    xlab(NULL)
  if(printRes){
    aic_tib %>%
      dplyr::arrange(AIC) %>%
      print(n=1e4)
  }
  return(outplot)
}

# makeFactorsLast <- function(x,lastx){
#   lastx <- c('Female')
#   return(factor(x,levels=c(sort(setdiff(x,lastx)),lastx)))
# }

coef.gergm <- function(g.obj,plot_coef=TRUE){
  c.obj <- t(g.obj@lambda.coef) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("var") %>%
    as_tibble() %>%
    dplyr::slice(c(grep('intercept$',var),
                   grep('dispersion$',var),
                   grep('dispersion$|intercept$',var,invert=T))) %>%
    dplyr::mutate(fold=exp(est),z=est/se,p=2*pnorm(-abs(z))) %>%
    dplyr::mutate(sig=case_when(p<=0.001 ~ '***',
                                p<=0.01  ~ '**',
                                p<=0.05 ~ '*',
                                p<=0.1   ~ '.',
                                T ~ ''))
  if(plot_coef){
    v1 <- c('intercept','dispersion')
    v2 <- setdiff(c.obj$var,v1)
    c.gg <- c.obj %>%
      dplyr::mutate(var=factor(var,levels=c(v1,v2))) %>%
      dplyr::filter(!grepl('dispersion',var)) %>%
      dplyr::mutate(sig=p<0.05,
                    lwr=z-2*1,
                    upr=z+2*1) %>%
      ggplot(aes(y=var,x=z))+
      theme_bw()+
      geom_vline(xintercept=0,linetype='dashed')+
      geom_errorbarh(aes(xmin=lwr,xmax=upr),height=0.1)+
      geom_point(aes(col=sig),size=2)+
      scale_y_discrete(limits=rev)+
      scale_colour_manual(values=c('black','red'))+
      theme(legend.position='none')+
      ylab(NULL)+
      xlab('Z-score')
    plot(c.gg)
  }
  print(c.obj)
}

fixF <- function(x){
  X <- strsplit(x,'\\+') %>% 
    unlist() %>%
    str_trim() %>% 
    stri_remove_empty_na() %>%
    paste(collapse='+')
  return(X)
}

makeModelList <- function(responseString,varStr,cnst=NULL,plot_coef=TRUE){
  levs <- c('nodefactor','nodematch','nodemix')
  n <- length(varStr)
  terms <- expand.grid(levs,varStr) %>%
    dplyr::mutate(terms=paste0(Var1,"('",Var2,"')")) %>%
    dplyr::pull(terms)
  
  # Formulae
  fs <- split(terms,ceiling(seq_along(terms)/3)) %>%
    lapply(function(x){return(c('',x))}) %>%
    expand.grid() %>%
    dplyr::slice(-1) %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    tidyr::unite(col=formula,sep=' + ') %>%
    rowwise() %>% 
    dplyr::mutate(formula=fixF(formula)) %>% 
    dplyr::mutate(f=paste0(responseString,' ~ edges + ',formula)) %>%
    dplyr::pull(f) %>%
    c(paste0(responseString,' ~ edges'),.)
  if(!is.null(cnst)){
    fs <- fs %>%
      paste0('+',cnst)
  }
  K <- length(fs)
  model.out <- vector(mode='list',length=K)
  for(k in 1:K){
    model.out[[k]] <- ergm(formula=formula(fs[[k]]),verbose=0)
  }
  
  return(list(models=model.out,formulae=fs))
}