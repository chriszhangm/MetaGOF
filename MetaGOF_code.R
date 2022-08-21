#install.packages("rstan")
#install.packages("mvtnorm")
library(rstan)
library(mvtnorm)


#xc a (non-empty) integer vector of events in the control group
#xt a (non-empty) integer vector of events in the treatment group
#nc a (non-empty) integer vector of the total number of patients in the control group
#nt a (non-empty) integer vector of the total number of patients in the treatment group
#method a character string specifying the method of Goodness-of-fit test, must be one of "IPQ","Chen","Wang" and "Naive". See details about these options.
#cov.prior Given method = "IPQ", a character string specifying the choice of covariance priors, must be one of "IW","HIW","SIW" and "IND".See details about these options.
#
# @details
# method = "IPQ": our proposed method.
# method = "Chen": A bootstrap-based approach from Chen et al.(2015).
# method = "Wang": A method based on the standardrization framework from Wang et al.(2019a).
# method = "Naive": Use Shaprio-Wilks test to estimated log odds ratio directly.
#
# cov.prior = "IW": Inverse-Wishart prior.
# cov.prior = "HIW": Huangs Inverse-Wishart prior proposed by Huang and Wand (2013).
# cov.prior = "SIW": Scaled Inverse-Whshart prior developed by O’Malley and Zaslavsky (2008).
# cov.prior = "IND": Independent prior.
#

# examples
# load(t2d.data.rda)
# xc = t2d.data$xc
# xt = t2d.data$xt
# nc = t2d.data$nc
# nt = t2d.data$nt
# MetaGOF(xc,xt,nc,nt)

MetaGOF = function(xc,xt,nc,nt,
                   method='IPQ',
                   cov.prior = 'IND'){
  if (length(xt)==0|length(xc)==0|length(nt)==0|length(nc)==0) {
    stop("Please input the data whose length is at least 1.")
  }
  #define functions
  #Simple Average Method (Initial Guess of log odds ratio)
  SA<-function(Xc,Xt,nc,nt)
  {
    thetaa<-log((Xt+1/2)/(nt-Xt+1/2))-log((Xc+1/2)/(nc-Xc+1/2))
    return(mean(thetaa))
  }
  #Initial guess of log odds of the control group
  mu0<-function(Xc,Xt,nc,nt)
  {
    return(mean(log((Xc+1/2)/(nc-Xc+1/2))))
  }
  #DSL approach for estimating between-study variance tauS
  tSDSL<-function(Xc,Xt,nc,nt) 
  {
    pt<-(Xt+1/2)/(nt+2*1/2)
    pc<-(Xc+1/2)/(nc+2*1/2)
    stS0<-1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
    wtS0<-1/stS0
    thetaa<-log((Xt+1/2)/(nt-Xt+1/2))-log((Xc+1/2)/(nc-Xc+1/2))
    thetaw0a<-sum(wtS0*thetaa)/sum(wtS0)
    Q<-sum(wtS0*(thetaa-thetaw0a)^2)
    tSDSL<-max(0,(Q-(length(Xc)-1))/(sum(wtS0)-sum(wtS0^2)/sum(wtS0)))
    return(tSDSL)
  }
  ###For Chen, Wang, and Naive methods only##
  #correction procedure is needed if we use chen or wang's approach. Otherwise, we just need skip the following code
  xc_copy = xc
  xt_copy = xt
  nc_copy = nc
  nt_copy = nt
  
  id_correct = which(xc_copy==0)
  xc_copy[id_correct] = 0.5
  nc_copy[id_correct] = nc_copy[id_correct] + 1
  
  id_correct = which(xt_copy==0)
  xt_copy[id_correct] = 0.5
  nt_copy[id_correct] = nt_copy[id_correct] + 1
  yi = log((xt_copy)/(nt_copy-xt_copy))-log((xc_copy)/(nc_copy-xc_copy))
  varOR = 1/(xt_copy) + 1/(xc_copy) + 1/(nt_copy-xt_copy) + 1/(nc_copy-xc_copy)
  ###

  
  if (method=="IPQ") {
    if(cov.prior=="IND"){
      stan_script = 'IND.stan'
    }
    else if(cov.prior=="IW"){
      stan_script = 'IW.stan'
    }
    else if(cov.prior=="HIW"){
      stan_script = 'HIW.stan'  
    }
    else if(cov.prior=="SIW"){
      stan_script = 'SIW.stan'
    }
    #fitting model
    k= length(xt)
    thetaraw = SA(xc,xt,nc,nt)
    xcraw = mu0(xc,xt,nc,nt)
    mydata= list(K=k,
                 J=2,
                 XC=xc,
                 XT=xt,
                 NXC=nc-xc,
                 NXT=nt-xt,
                 xcraw = xcraw,
                 thetaraw = thetaraw)
    fit = stan_model(stan_script)
    fit1 = sampling(fit,
                    data = mydata,
                    warmup =5000,
                    iter = 10000,
                    chains = 1,cores = 1,
                    refresh=0,
                    pars=c('theta','mu','overeff','tauS','rho','Sigma'))
    rrfit=summary(fit1)
    rr=rrfit$summary
    #posterior draws of all parameters
    phidif = as.matrix(extract(fit1,'theta')$theta[,,2] - extract(fit1,'theta')$theta[,,1])
    mu1draw = as.numeric(extract(fit1,'mu')$mu[,1])
    mudif = as.numeric(extract(fit1,'overeff')$overeff)
    tausamp = as.numeric(extract(fit1,'tauS')$tauS)
    rhodraw = as.numeric(extract(fit1,'rho')$rho)
    sigma1draw = as.numeric(extract(fit1,'Sigma')$Sigma[,1,1])
    sigma2draw = as.numeric(extract(fit1,'Sigma')$Sigma[,2,2])
    #output
    ind = 10000-5000
    zsamp = (phidif - mudif)/tausamp
    ptot = apply(X = zsamp, MARGIN = 1,FUN=function(x) shapiro.test(x)$p.value) #check normality
    #cauchy combination test
    cauchyt=sum((1/ind)*tan((0.5-ptot)*pi))
    IPQ=pcauchy(cauchyt,lower.tail = F)
    return(IPQ)
  }
  else if(method=="Chen"){
    
    s0 = as.numeric(shapiro.test(yi)$stat)
    tauSh = tSDSL(xc,xt,nc,nt)
    covhat = diag(tauSh+varOR)
    #based on 1000 replicates
    yi_rep = mvtnorm::rmvnorm(1000,sigma = covhat)
    sB = sapply(1:1000, FUN = function(x) return(as.numeric(shapiro.test(yi_rep[x,])$stat)))
    psw = sum(sB<s0)/1000
    return(psw)
  }
  else if(method=="Wang"){
    #wang's approach
    metawl <- function(new_mean, sd) {
      num <- length(new_mean)
      theta <- sum(new_mean*sd^(-2))/sum(sd^(-2))
      thetav <- rep(theta, num)
      Q <- sum((new_mean-thetav)^2*sd^(-2))
      tau2 <- (Q-(num-1))/(sum(sd^(-2))-sum(sd^(-4))/sum(sd^(-2)))
      thetar <- sum(new_mean/(sd^2+tau2))/sum((sd^2+tau2)^(-1))
      seTE2 <-
        sum((new_mean-thetar)^2*((sd^2+tau2)^(-1)/sum((sd^2+tau2)^(-1))))/(num-1)
      return(list(thetar=thetar, seTE2=seTE2, tau2=tau2))
    }
    sswtest <- function(new_mean, sd) {
      num <- length(new_mean)
      variance <- sd^2
      tau2 <- metawl(new_mean,sd)$tau2
      if (any(sd<0)) stop ("Standard error(s) < 0")
      if (tau2<=0) return(list(pvalue=999))
      stdmean <- rep(0, num)
      for (i in 1:num) {
        thetarJ <- metawl(new_mean[-i], sd[-i])$thetar
        seTE2J <- metawl(new_mean[-i], sd[-i])$seTE2
        stdmean[i] <- (new_mean[i]-thetarJ)*(1/(tau2+seTE2J+variance[i]))^0.5
      }
      pvalue <- shapiro.test(stdmean)$p.value
      result <- list(stdmean=stdmean, pvalue=pvalue)
      return(result)
    }
    wlt = sswtest(yi,sqrt(varOR))$pvalue
    return(wlt)
  }
  else if(method=="Naive"){
    naive = shapiro.test(yi)$p.value
    return(naive)
  }
}

