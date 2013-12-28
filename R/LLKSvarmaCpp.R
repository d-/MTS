LLKsvarmaCpp <- function(par){
   ## recall the relevant information.
   zt=Mtdata
   k=dim(zt)[2]
   nT=dim(zt)[1]
   nar=length(ARlags); nma=length(MAlags)
   istart=max(ARlags[nar],MAlags[nma])+1
   
   ListResiduals = .Call("GetSVarmaResiduals", Vtsdata, fix1, par, Order, ARlags, MAlags, Sresi, swi, inc.mean)		
   at  = do.call(rbind,ListResiduals)

   at=at[(istart:nT),]
   sig=t(at)%*%at/(nT-istart+1)
   ###ll=dmnorm(at,rep(0,k),sig)
   ll=dmvnorm(at,rep(0,k),sig)
   LLKsvarmaCpp=-sum(log(ll))
   
   LLKsvarmaCpp
}