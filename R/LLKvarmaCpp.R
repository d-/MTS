LLKvarmaCpp <- function(par){
zt    = Vtsdata
 k    = dim(Vtsdata)[2]
nT    = dim(Vtsdata)[1]
pqmax = max(ARp, MAq)

ListResiduals = .Call("GetVarmaResiduals", Vtsdata, fix1, par, ARp, MAq, inc.mean)		

at  = do.call(rbind,ListResiduals)
sig = t(at) %*% at/(nT-pqmax)
ll  = dmvnorm(at,rep(0,k),sig)

LLKvarmaCpp =-sum(log(ll))   
LLKvarmaCpp
}


