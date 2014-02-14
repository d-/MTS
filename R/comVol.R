"comVol" <- function(rtn,m=10,p=1,stand=FALSE){
# checking for common volatility components
#
if(!is.matrix(rtn))rtn=as.matrix(rtn)
# Fit a VAR(p) model to remove any serial correlations in the data.
x=VARfitC(rtn,p)$residuals
#
nT=dim(x)[1]
k=dim(x)[2]
#
if(m < 1)m=1
# standardize the returns
# mean of x is zero because VARfit employs a constant term. 
V1=cov(x)
##print(V1,digits=3)
m1=eigen(V1)
D1=diag(1/sqrt(m1$values))
P1=m1$vectors
Shalf=P1%*%D1%*%t(P1)
x1=x%*%Shalf
#
A=matrix(0,k,k)
for (h in 1:m){
ist=h+1
for (i in 1:k){
for (j in i:k){
Cmtx=matrix(0,k,k)
y2=x1[(ist-h):(nT-h),i]*x1[(ist-h):(nT-h),j]
for (ii in 1:k){
for (jj in ii:k){
y1=x1[ist:nT,ii]*x1[ist:nT,jj]
Cmtx[ii,jj]=cov(y1,y2)*(nT-h)/nT
Cmtx[jj,ii]=Cmtx[ii,jj]
}
}
Cmtx=Cmtx*((nT-h)/nT)
A= A+Cmtx%*%Cmtx
#end of j
}
#end of i
}
#end of h
}
#print(Cmtx)
if(stand){
dd=diag(A)
D=diag(1/sqrt(dd))
A=D%*%A%*%D
}
else{
A=A/(k*(k+1)/2)
}
m2=eigen(A)
Valu=m2$values
Prop=Valu/sum(Valu)
cat("eigen-values: ",Valu,"\n")
cat("proportion:   ",Prop,"\n")
Vec=m2$vectors
Mmtx=Shalf%*%Vec
# normalize each column of Mmtx
for (j in 1:k){
Mmtx[,j]=Mmtx[,j]/sqrt(sum(Mmtx[,j]^2))
}
#
# Perform ARCH tests for each transformed series
Tst=NULL
Tx = x%*%Mmtx
for (i in 1:k){
TT=NULL
mtst10=archTstC(Tx[,i],10)
TT=c(TT,mtst10)
mtst20=archTstC(Tx[,i],20)
TT=c(TT,mtst20)
mtst30=archTstC(Tx[,i],30)
TT=c(TT,mtst30)
Tst=rbind(Tst,c(i,TT))
}
cat("Checking: ","\n")
cat("Results of individual F-test for ARCH effect","\n")
cat("Numbers of lags used: 10, 20, 30","\n")
cat("Component,(F-ratio P-val) (F-ratio P-val) (F-ratio P-Val)","\n")
print(Tst,digits=3)

comVol <- list(residuals=x,values=m2$values,vectors=m2$vectors,M=Mmtx)
}

"VARfitC" <- function(xt,p){
#fit a AVR(p) model to the data
if(!is.matrix(xt))xt=as.matrix(xt)
nT=dim(xt)[1]
k=dim(xt)[2]
Resi=xt
beta=NULL
#
if(p > 0){
ist=p+1
Y=xt[ist:nT,]
effN=nT-p
X=matrix(1,effN,1)
for (h in 1:p){
X=cbind(X,xt[(ist-h):(nT-h),])
}
XtX=t(X)%*%X
XtY=t(X)%*%Y
XtXinv=solve(XtX)
beta=XtXinv%*%XtY
Resi=Y-X%*%beta
}

VARfitC <- list(coefs=beta,residuals=Resi)

}

#
"archTstC" <- function(x,m){
# perform F-test for ARCH effect using x^2 series 
# m*Fratio is asymptotically chi-square with m degrees of freedom.
#
if(m < 1)m=1
nT=length(x)
ist=m+1
EffN = nT-m
Xmtx=matrix(1,EffN,1)
Y=matrix(x[ist:nT]^2,EffN,1)
for (j in 1:m){
Xmtx=cbind(Xmtx,x[(ist-j):(nT-j)]^2)
}
XtX=t(Xmtx)%*%Xmtx
XtY=t(Xmtx)%*%Y
XtXinv=solve(XtX)
beta=XtXinv%*%XtY
Resi=Y-Xmtx%*%beta
Ywm=Y-mean(Y)
SSR=sum(Resi^2)
deg=EffN-m-1
Fratio=((sum(Ywm^2)-SSR)/m)/(SSR/deg)
pv=1-pf(Fratio,m,deg)
result=c(Fratio,pv)

result
}