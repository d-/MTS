"BVAR" <- function(z,p=1,C,V0,n0=5,Phi0=NULL,include.mean=T){
## Perform Bayesian estimation of a VAR(p) model
##
## z: time series (T-by-k)
## p: AR order
## phi0: prior mean for coefficient matrix [k-by-(kp+1)]
## C: precision matrix of coefficient matrix. [(kp+1)-by-(kp+1)]
## (V0,n0): prior input for Sigma_a (inverted Wishart parameters)
##
if(!is.matrix(z))z=as.matrix(z)
if(p < 1) p=1
if(!is.matrix(C))C=as.matrix(C)
if(!is.matrix(V0))V0=as.matrix(V0)
if(n0 < 1)n0=1
##
T=dim(z)[1]
k=dim(z)[2]
idim=k*p+1
if(length(Phi0) <= 0)Phi0=matrix(0,idim,k)
X=NULL
ne=T-p
if(include.mean)X=rep(1,ne)
for (i in 1:p){
X=cbind(X,z[(p+1-i):(T-i),])
}
Z=as.matrix(z[(p+1):T,])
X=as.matrix(X)
### 
XpX=t(X)%*%X
XpY=t(X)%*%Z
## Bayesian Estimate
WpW=XpX+C
WpWinv=solve(WpW)
WpY=XpY+C%*%Phi0
Bbhat=WpWinv%*%WpY
bAhat=Z-X%*%Bbhat
bB=Bbhat-Phi0
S=t(bAhat)%*%bAhat +t(bB)%*%C%*%bB
BSig=(V0+S)/(n0+ne-k-1)
SD=kronecker(BSig,WpWinv)
phi=c(Bbhat)
se=sqrt(diag(SD))
Est=cbind(phi,se,phi/se)
colnames(Est) <- c("Est","s.e.","t-ratio")
cat("Bayesian estimate:","\n")
print(Est)
cat("Covariance matrix: ","\n")
print(BSig)
}

