
##### Seasonal Models
"sVARMACpp" <- function(da,order,sorder,s,include.mean=T,fixed=NULL,prelim=F,details=F,thres=0.8,switch=F){
   # Estimation of a multiplicative vector ARMA model using conditional MLE (Gaussian dist)
   #  This program is modified from VARMA.R. It was created first on April 21, 2011.
   #### Modified on September 19, 2012 by adding "switch=F" to constrol Theta(B)*theta(B)
   #
   # The following comments are for VARMA.R program.
   # April 16: remove some insigificant parameters at the preliminary estimation so that
   #           the computational speed can be improved.
   # April 18: add subcommand "prelim" to see simplification after the AR approximation.
   # When prelim=TRUE, fixed is assigned based on the results of AR approximation.
   # Here "thres" is used only when prelim = TRUE.
   #
   # The seasonal order is (AR,MA) = (P,Q). The seasonality is s.
   #
   if(!is.matrix(da))da=as.matrix(da)
   p=order[1];d=order[2];q=order[3];P=sorder[1];D=sorder[2];Q=sorder[3]
   nT=dim(da)[1]; k=dim(da)[2]
   # basic setup.
   if(p < 0)p=0; if(q < 0)q=0; if(P < 0) P = 0; if(Q < 0) Q = 0; if(s < 0) s=-s
   if(d > 1){
      cat("Regular difference is adjusted to d=1","\n")
      d=1
   }
   if(D > 1){
      cat("Seasonal difference is adjusted to D=1","\n")
      D=1
   }
   ### number of coefficient parameters
   kp=k*p
   kq=k*q
   kP=k*P
   kQ=k*Q
   # Take care of the difference
   if(d==1){
      X=NULL
      MEAN=rep(0,k)
      for (j in 1:k){
         X=cbind(X,diff(da[,j]))
         t1=t.test(X[,j])
         if(t1$p.value < 0.05)MEAN[j]=1
      }
      if(sum(MEAN) < 1)include.mean=F
   }
   else{
      X=da
   }
   #
   if(D==1){
      DX=NULL
      Smean=rep(0,k)
      for (j in 1:k){
         DX=cbind(DX,diff(X[,j],s))
         t1=t.test(DX[,j])
         if(t1$p.value < 0.05)Smean[j]=1
      }
      if(sum(Smean) < 1)include.mean=F
   }
   else{
      DX=X
   }
   #### sample size after differencing.
   nT=dim(DX)[1]
   arlags=NULL
   if(p > 0){
      arlags=c(1:p)
      if(P > 0)arlags=c(arlags,c(1:P)*s,c(1:P)*s+c(1:p))
   }
   else{
      if(P > 0)arlags=c(1:P)*s
   }
   malags=NULL
   if(q > 0){
      malags=c(1:q)
      if(Q > 0)malags=c(malags,c(1:Q)*s,c(1:Q)*s+c(1:q))
   }
   else{
      if(Q > 0)malags=c(1:Q)*s
   }
   # number of AR and MA lags of the model
   nar=length(arlags)
   nma=length(malags)
   idim=k*(nar+nma)
   if(include.mean)idim=idim+1
   if(length(fixed)==0)fixed=matrix(1,idim,k)
   ### Assign the data globally for estimation purpose.
   Mtdata <<- DX
   ARlags <<- arlags
   MAlags <<- malags
   Order <<- c(order,sorder)
   inc.mean <<- include.mean
   fix1 <<- fixed
   swi <<- switch
   ####
   phi=NULL; sphi=NULL; sephi=NULL; sesphi=NULL
   if(p > 0)phi=matrix(0,k,k*p); sephi=phi
   if(P > 0)sphi=matrix(0,k,k*P);sesphi=sphi
   theta=NULL; stheta=NULL;setheta=NULL; sestheta=NULL
   if(q > 0)theta=matrix(0,k,k*q);setheta=theta
   if(Q > 0)stheta=matrix(0,k,k*Q);sestheta=stheta
   ## Obtain initial estimates of the component parameters using univariate models.
   ### For cross-series initial estimates, we use linear models with univariate at-series
   resi=NULL
   for (j in 1:k){
      m1=arima(DX[,j],order=c(p,0,q),seasona=list(order=c(P,0,Q),period=s))
      resi=cbind(resi,m1$residuals)
      seest=sqrt(diag(m1$var.coef))
      ##print(m1$coef)
      icnt=0
      # locate regular AR coef
      if(p > 0){
         for (i in 1:p){
            icnt=icnt+1
            ii=(i-1)*k
            phi[j,(ii+j)]=m1$coef[icnt]
            sephi[j,(ii+j)]=seest[icnt]
         }
      }
      # locate regular MA coef.
      if(q > 0){
         for (i in 1:q){
            ii=(i-1)*k
            icnt=icnt+1
            theta[j,(ii+j)]=-m1$coef[icnt]
            setheta[j,(ii+j)]=seest[icnt]
         }
      }
      # locate seasonal AR coef.
      if(P > 0){
         for (i in 1:P){
            icnt=icnt+1
            ii=(i-1)*k
            sphi[j,(ii+j)]=m1$coef[icnt]
            sesphi[j,(ii+j)]=seest[icnt]
         }
      }
      # locate seasonal MA coef
      if(Q > 0){
         for (i in 1:Q){
            ii=(i-1)*k
            icnt=icnt+1
            stheta[j,(ii+j)]=-m1$coef[icnt]
            sestheta[j,(ii+j)]=seest[icnt]
         }
      }
      # end of the j-loop
   }
   #### Obtain estimates of cross-series parameters, using Least-Squares approximation.
   m2=siniEST(DX,resi,arlags,malags,include.mean)
   #### Fill in the coefficient matrices
   beta=t(m2$estimates)
   sebeta=t(m2$se)
   ##
   icnst=0
   if(include.mean)icnst=1
   if(nar > 0){
      if(p > 0){
         for (i in 1:p){
            idx=(i-1)*k
            for (ii in 1:k){
               jdx=c(1:k)[-ii]
               phi[jdx,(idx+ii)]=beta[jdx,(icnst+idx+ii)]
               sephi[jdx,(idx+ii)]=sebeta[jdx,(icnst+idx+ii)]
            }
            #end of i-loop
         }
         #end of if(p > 0) statement
      }
      if(P > 0){
         for (i in 1:P){
            kdx=(i-1)*k
            idx=k*p+kdx
            for (ii in 1:k){
               jdx=c(1:k)[-ii]
               sphi[jdx,(kdx+ii)]=beta[jdx,(icnst+idx+ii)]
               sesphi[jdx,(kdx+ii)]=sebeta[jdx,(icnst+idx+ii)]
            }
         }
      }
      #end of if(nar > 0) statement
   }
   if(nma > 0){
      if(q > 0){
         for (i in 1:q){
            kdx=(i-1)*k
            idx=nar*k+kdx
            for (ii in 1:k){
               jdx=c(1:k)[-ii]
               theta[jdx,(kdx+ii)]=-beta[jdx,(icnst+idx+ii)]
               setheta[jdx,(kdx+ii)]=sebeta[jdx,(icnst+idx+ii)]
            }
         }
      }
      if(Q > 0){
         for (i in 1:Q){
            kdx=(i-1)*k
            idx=(nar+q)*k+kdx
            for (ii in 1:k){
               jdx=c(1:k)[-ii]
               stheta[jdx,(kdx+ii)]=-beta[jdx,(icnst+idx+ii)]
               sestheta[jdx,(kdx+ii)]=sebeta[jdx,(icnst+idx+ii)]
            }
         }
      }
      #end of the statement if(nma > 0)
   }
   # Identify parameters to be estimated.
   par=NULL
   separ=NULL
   ist=0
   ## We took the transpose of beta and sebeta after siniEST program.
   if(include.mean){
      jdx=c(1:k)[fix1[1,]==1]
      if(length(jdx) > 0){
         par=beta[jdx,1]
         separ=sebeta[jdx,1]
      }
      ist=1
   }
   if(nar > 0){
      if(p > 0){
         for (j in 1:k){
            idx=c(1:kp)[fix1[(ist+1):(ist+kp),j]==1]
            if(length(idx) > 0){
               par=c(par,phi[j,idx])
               separ=c(separ,sephi[j,idx])
            }
            #end of j-loop
         }
         ist=ist+kp
      }
      ###print(ist)
      if(P > 0){
         for (j in 1:k){
            idx=c(1:kP)[fix1[(ist+1):(ist+kP),j]==1]
            if(length(idx) > 0){
               par=c(par,sphi[j,idx])
               separ=c(separ,sesphi[j,idx])
            }
         }
         #end of if(P > 0) statement
         ist=ist+kP
      }
      # end of if(nar > 0) statement
   }
   #
   if(nma > 0){
      if(q > 0){
         for (j in 1:k){
            idx=c(1:kq)[fix1[(ist+1):(ist+kq),j]==1]
            if(length(idx) > 0){
               par=c(par,theta[j,idx])
               separ=c(separ,setheta[j,idx])
            }
            #end of j-loop
         }
         ist=ist+kq
      }
      if(Q > 0){
         for (j in 1:k){
            idx=c(1:kQ)[fix1[(ist+1):(ist+kQ),j]==1]
            if(length(idx) > 0){
               par=c(par,stheta[j,idx])
               separ=c(separ,sestheta[j,idx])
            }
         }
         #end of if(Q > 0) statement
      }
      #end of (if nma > 0) statement
   }
   #### keep the first few residuals to be used in likelihood evaluation to compute "at".
   jst=max(arlags[nar],malags[nma])
   Sresi <<- resi[1:jst,]
   #########
   cat("Number of parameters: ",length(par),"\n")
   cat("initial estimates: ",par,"\n")
   ### Set up lower and upper bounds
   lowerBounds=par; upperBounds=par
   for (j in 1:length(par)){
      lowerBounds[j] = par[j]-2*separ[j]
      upperBounds[j] = par[j]+2*separ[j]
   }
   ###mm=optim(par,LLKvma,method=c("L-BFGS-B"),lower=lowerBounds,upper=upperBounds,hessian=TRUE)
   ###mm=optim(par,LLKvma,method=c("BFGS"),hessian=TRUE)
   ##est=mm$par
   ##H=mm$hessian
   # Step 5: Estimate Parameters and Compute Numerically Hessian:
   if(details){
      fit = nlminb(start = par, objective = LLKsvarmaCpp,
      lower = lowerBounds, upper = upperBounds, control = list(trace=3))
   }
   else {
      fit = nlminb(start = par, objective = LLKsvarmaCpp, lower = lowerBounds, upper = upperBounds)
   }
   epsilon = 0.0001 * fit$par
   npar=length(par)
   Hessian = matrix(0, ncol = npar, nrow = npar)
   for (i in 1:npar) {
      for (j in 1:npar) {
         x1 = x2 = x3 = x4  = fit$par
         x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
         x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
         x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
         x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
         Hessian[i, j] = (LLKsvarmaCpp(x1)-LLKsvarmaCpp(x2)-LLKsvarmaCpp(x3)+LLKsvarmaCpp(x4))/
         (4*epsilon[i]*epsilon[j])
      }
   }
   # Step 6: Create and Print Summary Report:
   se.coef = sqrt(diag(solve(Hessian)))
   tval = fit$par/se.coef
   matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
   dimnames(matcoef) = list(names(tval), c(" Estimate",
   " Std. Error", " t value", "Pr(>|t|)"))
   cat("\nCoefficient(s):\n")
   printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
   est=fit$par
   #
   ### restore estimates to the format of unconstrained case for printing purpose.
   #### icnt: parameter count
   #### ist: location count
   ist=0
   icnt = 0
   Ph0=rep(0,k)
   sePh0=rep(0,k)
   beta=NULL
   sebeta=NULL
   if(include.mean){
      idx=c(1:k)[fix1[1,]==1]
      icnt=length(idx)
      if(icnt > 0){
         Ph0[idx]=est[1:icnt]
         sePh0[idx]=se.coef[1:icnt]
      }
      ist=1
      beta=rbind(beta,Ph0)
      sebeta=rbind(sebeta,sePh0)
   }
   PH=NULL; sePH=NULL; sPH=NULL; sesPH=NULL
   if(p > 0){
      PH=matrix(0,kp,k)
      sePH=matrix(0,kp,k)
      for (j in 1:k){
         idx=c(1:kp)[fix1[(ist+1):(ist+kp),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            PH[idx,j]=est[(icnt+1):(icnt+jdx)]
            sePH[idx,j]=se.coef[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
         # end of j-loop
      }
      #end of if (p > 0)
      ist=ist+kp
      beta=rbind(beta,PH)
      sebeta=rbind(sebeta,sePH)
   }
   if(P > 0){
      sPH=matrix(0,kP,k)
      sesPH=matrix(0,kP,k)
      for (j in 1:k){
         idx=c(1:kP)[fix1[(ist+1):(ist+kP),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            sPH[idx,j]=est[(icnt+1):(icnt+jdx)]
            sesPH[idx,j]=se.coef[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
      }
      #end of if (P > 0) statement
      ist=ist+kP
      beta=rbind(beta,sPH)
      sebeta=rbind(sebeta,sesPH)
   }
   #
   TH=NULL;seTH=NULL; sTH=NULL; sesTH=NULL
   if(q > 0){
      TH=matrix(0,kq,k)
      seTH=matrix(0,kq,k)
      for (j in 1:k){
         idx=c(1:kq)[fix1[(ist+1):(ist+kq),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            TH[idx,j]=est[(icnt+1):(icnt+jdx)]
            seTH[idx,j]=se.coef[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
         # end of j-loop
      }
      # end of if(q > 0).
      ist=ist+kq
      beta=rbind(beta,-TH)
      sebeta=rbind(sebeta,seTH)
   }
   if(Q > 0){
      sTH=matrix(0,kQ,k)
      sesTH=matrix(0,kQ,k)
      for (j in 1:k){
         idx=c(1:kQ)[fix1[(ist+1):(ist+kQ),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            sTH[idx,j]=est[(icnt+1):(icnt+jdx)]
            sesTH[idx,j]=se.coef[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
      }
      beta=rbind(beta,-sTH)
      sebeta=rbind(sebeta,sesTH)
   }
   ###
   cat("---","\n")
   cat("Estimates in matrix form:","\n")
   if(include.mean){
      cat("Constant term: ","\n")
      cat("Estimates: ",Ph0,"\n")
   }
   if(p > 0){
      cat("Regular AR coefficient matrix","\n")
      jcnt=0
      for (i in 1:p){
         cat("AR(",i,")-matrix","\n")
         ph=t(PH[(jcnt+1):(jcnt+k),])
         print(ph,digits=3)
         jcnt=jcnt+k
      }
      # end of if (p > 0)
   }
   ### Seasonal AR part
   if(P > 0){
      cat("Seasonal AR coefficient matrix","\n")
      jcnt=0
      for (i in 1:P){
         cat("AR(",i*s,")-matrix","\n")
         ph=t(sPH[(jcnt+1):(jcnt+k),])
         print(ph,digits=3)
         jcnt=jcnt+k
      }
   }
   ## Regular MA part
   if(q > 0){
      cat("Regular MA coefficient matrix","\n")
      icnt=0
      for (i in 1:q){
         cat("MA(",i,")-matrix","\n")
         the=t(TH[(icnt+1):(icnt+k),])
         print(the,digits=3)
         icnt=icnt+k
      }
      # end of the statement if(q > 0)
   }
   # Seasonal MA part
   if(Q > 0){
      cat("Seasonal MA coefficient matrix","\n")
      icnt=0
      for (i in 1:Q){
         cat("MA(",i*s,")-matrix","\n")
         the=t(sTH[(icnt+1):(icnt+k),])
         print(the,digists=3)
         icnt=icnt+k
      }
   }
   ######### Obtain product coefficient matrices
   if((p > 0)&&(P > 0)){
      if(switch){
         Phi=t(Mtxprod1(t(PH),t(sPH),p,P))
      }
      else{
         Phi=t(Mtxprod(t(PH),t(sPH),p,P))
      }
   }
   #
   if((p > 0)&&(P==0))Phi=PH
   if((p==0)&&(P > 0))Phi=sPH
   if((q > 0)&&(Q > 0)){
      if(switch){
         Theta=t(Mtxprod1(t(TH),t(sTH),q,Q))
      }
      else{
         Theta=t(Mtxprod(t(TH),t(sTH),q,Q))
      }
   }
   #
   if((q > 0)&&(Q==0))Theta=TH
   if((q==0)&&(Q > 0))Theta=sTH
   ##### Compute the residuals
   zt=Mtdata
   pqmax=max(ARlags[nar],MAlags[nma])
   ist=pqmax+1
   #### consider the case t from ist to T
   at=Sresi[1:pqmax,]
   for (t in ist:nT){
      tmp=zt[t,]-Ph0
      if(nar > 0){
         for (j in 1:nar){
            jj=ARlags[j]
            jdx=(j-1)*k
            ph=Phi[(jdx+1):(jdx+k),]
            tmp=tmp-matrix(zt[(t-jj),],1,k)%*%ph
         }
      }
      if(nma > 0){
         for (j in 1:nma){
            jj=MAlags[j]
            jdx=(j-1)*k
            th=Theta[(jdx+1):(jdx+k),]
            tmp=tmp+matrix(at[(t-jj),],1,k)%*%th
         }
      }
      at=rbind(at,tmp)
   }
   #
   at=at[(ist:nT),]
   c1 = rep("resi",k)
   colnames(at) <- c1
   sig=t(at)%*%at/(nT-pqmax)
   ##
   cat(" ","\n")
   cat("Residuals cov-matrix:","\n")
   print(sig)
   dd=det(sig)
   d1=log(dd)
   aic=d1+2*npar/nT
   bic=d1+log(nT)*npar/nT
   cat("----","\n")
   cat("aic= ",round(aic,4),"\n")
   cat("bic= ",round(bic,4),"\n")
   ### prepare for output storage.
   if(length(PH) > 0)PH=t(PH)
   if(length(sPH) > 0)sPH=t(sPH)
   if(length(TH) > 0)TH=t(TH)
   if(length(sTH) > 0)sTH=t(sTH)
   
   
   sVARMACpp <- list(data=da,order=order,sorder=sorder,period=s,cnst=include.mean,coef=beta,secoef=sebeta,residuals=at,Sigma=sig,aic=aic,bic=bic,regPhi=PH,seaPhi=sPH, regTheta=TH, seaTheta=sTH, Ph0=Ph0,switch=switch)
}
