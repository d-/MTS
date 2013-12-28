LLKvmaCpp <- function(par){
   zt=VMAdata
   k=ncol(zt)
   nT=nrow(zt)
   q <- MAq
   kq=k*q
   fix <- fix1
   include.mean <- inc.mean
   
   ListObs = .Call("GetVMAObs", VMAdata, fix1, par, MAq, inc.mean)
   zt  = do.call(rbind, ListObs)
   
   #### Slow!
   ListTH = .Call("GetVMATH", zt, fix1, par, MAq, inc.mean)
   TH  = do.call(rbind, ListTH)
       
   mm=eigen(TH)
   V1=mm$values
   P1=mm$vectors
   v1=Mod(V1)
   ich=0
   for (i in 1:kq){
      if(v1[i] > 1)V1[i]=1/V1[i]
      ich=1
   }
   if(ich > 0){
      ###cat("Invertibility checked and adjusted: ","\n")
      P1i=solve(P1)
      GG=diag(V1)
      TH=Re(P1%*%GG%*%P1i)
      Theta=t(TH[1:k,])
      ##cat("adjusted Theta","\n")
      ##print(TH[1:k,])
      ### re-adjust the MA parameter
      icnt <- VMAcnt
      ist=0
      if(icnt > 0)ist=1
      for (j in 1:k){
         idx=c(1:kq)[fix[(ist+1):(ist+kq),j]==1]
         jcnt=length(idx)
         if(jcnt > 0){
            par[(icnt+1):(icnt+jcnt)]=TH[j,idx]
            icnt=icnt+jcnt
         }
      }
      ##
      ParMA <- par
   }
   
   ##
   at=mFilter(zt,t(Theta))
   #
   sig=t(at)%*%at/nT
   ##ll=dmnorm(at,rep(0,k),sig)
   ll=dmvnorm(at,rep(0,k),sig)
   LLKvmaCpp =-sum(log(ll))
   LLKvmaCpp
}


