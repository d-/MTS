require(fGarch)
"ncAR" <- function(xt,order=c(1,1),include.mean=TRUE,method="Nelder-Mead"){
### Perform quick estimation of non-causal AR models
nT=length(xt)
p1=order[1]; p2=order[2]; p=p1+p2
## The model is factored as ncAR * cAR = at with at being Student-t innovation.
##
para=NULL; loB=NULL; upB=NULL
meanX=mean(xt); xtm=xt
if(include.mean){para=c(para,meanX); xtm=xt-meanX}
if(p1 > 0){
 m1=arima(xtm,order=c(p1,0,0),include.mean=F)
 coef=m1$coef[1:p1]; para=c(para,car=coef)
 }
if(p2 > 0){
 coef=rep(1.6,p2); coef.se=rep(0.2,p2)
 para=c(para,ncar=coef)
 }
para=c(para,sig=log(sqrt(var(xtm))),v=5)
##
ncARlike <- function(para,xt,order,include.mean){
 p1=order[1]; p2=order[2]; nT=length(xt); yt=xt; icnt=0
 if(include.mean){icnt=1; yt=xt-para[1]}
 if(p1 > 0){cphi=para[(icnt+1):(icnt+p1)]; icnt=icnt+p1}
 if(p2 > 0){ncphi=para[(icnt+1):(icnt+p2)]; icnt=icnt+p2}
 Ut = yt
 if(p2 > 0){
  ist = p2+1
  for (i in 1:p2){
   Ut[(p2+1):nT] = Ut[ist:nT]-yt[(ist-i):(nT-i)]*ncphi[i]
   }
  }
  ist=p1+p2+1
  At=Ut
  if(p1 > 0){
   for (i in 1:p1){
     At[ist:nT]=At[ist:nT]-Ut[(ist-i):(nT-i)]*cphi[i]
    }
   }    
 sig=exp(para[icnt+1]); v=para[icnt+2]
 lastcoef=1
 if(p2 > 0)lastcoef=abs(ncphi[p2])
 like = sum(dstd(At[ist:nT]/sig,mean=0,sd=1,nu=v,log=T)) - (nT-ist+1)*log(sig)
 f <- -like-(nT-ist+1)*log(lastcoef)
 }
##
 fit <- optim(para, ncARlike,method=method,hessian=TRUE, xt=xt, order=order, include.mean=include.mean)
     if (fit$convergence) 
            warning("optimization may not have succeeded")
        par.ests <- fit$par
 cat("est:", par.ests, "\n")
        converged <- fit$convergence
        nllh.final <- fit$value
     varcov=solve(fit$hessian)
 cat("est.se:", sqrt(diag(varcov)), "\n")
 cat("log likelihood:", -nllh.final, "\n")
 
### calculate the residuals
 yt=xt; icnt=0
 if(include.mean){yt=xt-par.ests[1]; icnt=1}
 if(p1 > 0){cphi=par.ests[(icnt+1):(icnt+p1)]; icnt=icnt+p1}
 if(p2 > 0){ncphi=par.ests[(icnt+1):(icnt+p2)]; icnt=icnt+p2}
 Ut=yt; ist=1
 if(p2 > 0){ ist=p2+1
  for (i in 1:p2){
   Ut[ist:nT]=Ut[ist:nT]-ncphi[i]*yt[(ist-i):(nT-i)]
   }
  }
  at = Ut
  if( p1 > 0){ ist=p1+p2+1
   for (i in 1:p1){
     at[ist:nT]=at[ist:nT]-cphi[i]*Ut[(ist-i):(nT-i)]
     }
   }
 sig=exp(par.ests[icnt+1])
 cat("Standard error of the innovation: ", sig,"\n")

 ncAR <- list(par=par.ests, varcov=varcov, residuals=at[ist:nT])

}

#### Generating ncAR models
"gncAR" <- function(nobs,poly,cpoly=NULL,v=5,sig=1){
## Generate non-causual AR models
## poly denoets the coefficients of the AR polynomial, starting with 1.
## For instance, poly=c(1,-2) denotes (1-2B)
## cpoly: causal AR polynomial, starting with 1.
## For instance, cpoly=c(1,-0.5) denotes (1-0.5B)
nT=nobs+100
at=rstd(nT,mean=0,sd=1,nu=v)*sig
## check the causal part
yt=at ## initialization
p1=length(cpoly)
if(p1 > 1){
coef=-cpoly[2:p1]
yt=filter(at,coef,method="r",init=rep(0,p1-1))
}
##
xt=yt
p=length(poly)
if(p > 0){
phi=poly/poly[p]
coefat=phi[1]; coef=-rev(phi[1:(p-1)])
p=p-1  # the actual order
xt[(nT-p+1):nT] = 0  # start with zero initial values
for (i in 1:(nT-p)){
 it = nT-p -i + 1
 tmp = coefat * yt[it+p]
 for (j in 1:p){
  tmp=tmp+ coef[j]*xt[it+j]
  }
 xt[it]=tmp
 }
}
xt=xt[1:(nT-100)]

gncAR <- list(at = at[1:(nT-100)], series=xt)
}

#### Use Forward operators, i.e. F = B^{-1}.
####
"ncARf" <- function(xt,order=c(1,1),include.mean=TRUE,method="Nelder-Mead"){
### Perform quick estimation of non-causal AR models
nT=length(xt)
p1=order[1]; p2=order[2]; p=p1+p2
## The model is factored as phi(B)*psi(F)xt = at with at being Student-t innovation.
##
para=NULL 
meanX=mean(xt); xtm=xt
if(include.mean){para=c(para,meanX); xtm=xt-meanX}
if(p1 > 0){
 m1=arima(xtm,order=c(p1,0,0),include.mean=F)
 coef=m1$coef[1:p1]; para=c(para,car=coef)
 }
if(p2 > 0){
 coef=rep(0.2,p2); para=c(para,ncar=coef)
 }
para=c(para,sig=log(sqrt(var(xtm))),v=5)
##
ncARt <- function(para,xt,order,include.mean){
 p1=order[1]; p2=order[2]; nT=length(xt); yt=xt; icnt=0; iend=nT
 if(include.mean){icnt=1; yt=xt-para[1]}
 if(p1 > 0){cphi=para[(icnt+1):(icnt+p1)]; icnt=icnt+p1}
 if(p2 > 0){ncphi=para[(icnt+1):(icnt+p2)]; icnt=icnt+p2}
 Ut = yt
 if(p2 > 0){
  ist = 1; iend=nT-p2
  for (i in 1:p2){
   Ut[ist:iend] = Ut[ist:iend]-yt[(ist+i):(iend+i)]*ncphi[i]
   }
  }
  ist=p1+1
  At=Ut
  if(p1 > 0){
   for (i in 1:p1){
     At[ist:iend]=At[ist:iend]-Ut[(ist-i):(iend-i)]*cphi[i]
    }
   }    
 sig=exp(para[icnt+1]); v=para[icnt+2]
 like = sum(dstd(At[ist:iend]/sig,mean=0,sd=1,nu=v,log=T)) - (iend-ist+1)*log(sig)
 f <- -like
 }
##
 fit <- optim(para, ncARt,method=method,hessian=TRUE, xt=xt, order=order, include.mean=include.mean)
     if (fit$convergence) 
            warning("optimization may not have succeeded")
        par.ests <- fit$par
 cat("est:", par.ests, "\n")
        converged <- fit$convergence
        nllh.final <- fit$value
     varcov=solve(fit$hessian)
 cat("est.se:", sqrt(diag(varcov)), "\n")
 cat("log likelihood:", -nllh.final, "\n")
 
### calculate the residuals
 yt=xt; icnt=0
 if(include.mean){yt=xt-par.ests[1]; icnt=1}
 if(p1 > 0){cphi=par.ests[(icnt+1):(icnt+p1)]; icnt=icnt+p1}
 if(p2 > 0){ncphi=par.ests[(icnt+1):(icnt+p2)]; icnt=icnt+p2}
 Ut=yt; iend=nT-p2; ist=1
 if(p2 > 0){
  for (i in 1:p2){
   Ut[ist:iend]=Ut[ist:iend]-ncphi[i]*yt[(ist+i):(iend+i)]
   }
  }
  at = Ut
  if( p1 > 0){ist=p1+1
   for (i in 1:p1){
     at[ist:iend]=at[ist:iend]-cphi[i]*Ut[(ist-i):(iend-i)]
     }
   }
 sig=exp(par.ests[icnt+1])
 cat("Standard error of the innovations: ",sig,'\n')
 ncARf <- list(par=par.ests, varcov=varcov, residuals=at[ist:iend])

}

### Genreate ncARa models
### 
"gncARf" <- function(nT,cpoly=NULL,poly=NULL, sig=1, v=5,skip=100){
##
## The model is poly(F)*cpoly(B)xt = at
## where at = sig*et, with et being a standardized Student-t innovation with v degrees of freedom.
##
## poly = (1, -0.5) means (1-0.5F); cpoly=c(1,-0.5) means (1-0.5B)
### skip: the number of data points dropped to remove the impact of initial values
##
p2 = length(poly)-1; p1=length(cpoly)-1
n = nT+2*skip
at = sig*rstd(n,mean=0,sd=1,nu=v)
 yt=at
if(p1 > 0){
  f1=-cpoly[2:(p1+1)]
  yt=filter(at,f1,method="r",init=rep(0,p1))
 }
 xt=yt
if(p2 > 0){
 nphi = -poly[2:(p2+1)]
 iend=n-p2
 for (it in 1:iend){
  itt = iend-it+1
  tmp = yt[itt]
   for (i in 1:p2){
    tmp = tmp + nphi[i]*xt[itt+i]
    }
   xt[itt]=tmp
   }
  }
 xt = xt[(skip+1):(skip+nT)]
 at = at[(skip+1):(skip+nT)]

gncARf <- list(series=xt,at=at)
}
  