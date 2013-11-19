\name{MBcluster}
\alias{MBcluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Model-Based Clustering Analysis
%%  ~~function to do ... ~~
}
\description{Perform clustering analysis via the MCMC methods
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MBcluster(data, p, K, mcmc, prior = NULL, differ = FALSE, start.values = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A T-by-k data matrix ofk-dimensional  time series
%%     ~~Describe \code{data} here~~
}
  \item{p}{Univariate AR order
%%     ~~Describe \code{p} here~~
}
  \item{K}{The number of clusters
%%     ~~Describe \code{K} here~~
}
  \item{mcmc}{A vector consisting of the number of burn-in and the total 
number of itertations of MCMC methods.
%%     ~~Describe \code{mcmc} here~~
}
  \item{prior}{prior distribution 
%%     ~~Describe \code{prior} here~~
}
  \item{differ}{Differencing the data
%%     ~~Describe \code{differ} here~~
}
  \item{start.values}{Initial values
%%     ~~Describe \code{start.values} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Parameter estimates for each clusters and members of each cluster
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 6)
%% ~put references to the literature/web site here ~
}
\author{Yongning Wang
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-unempstatesAdj.txt",header=T)
zt=diffM(da)
iter=c(2000,5000)
m1=MBcluster(zt,4,4,mcmc=iter)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, p, K, mcmc, prior = NULL, differ = FALSE, start.values = NULL) 
{
    data = as.matrix(data)
    N = length(data[1, ])
    colname = rep(0, N)
    if (length(colnames(data)) == 0) {
        for (n in 1:N) {
            colname[n] = paste("Series ", n)
        }
    }
    else {
        colname = colnames(data)
    }
    if (differ == TRUE) {
        data = apply(data, 2, diff)
    }
    T = length(data[, 1])
    if (is.null(prior)) {
        cat("Use default priors", "\n")
        phi0 = rep(0, p + 1)
        Sigma0_phi = 10^4 * diag(p + 1)
        n0 = 2
        sig02 = 0.5
        e0 = rep(3, K)
        prior = list(phi0 = phi0, Sigma0_phi = Sigma0_phi, n0 = n0, 
            sig02 = sig02, e0 = e0)
    }
    s = rep(0, N)
    S = matrix(0, (mcmc$burnin + mcmc$rep), N)
    phi = matrix(0, K * (mcmc$burnin + mcmc$rep), (p + 1))
    sigma2 = rep(0, K * (mcmc$burnin + mcmc$rep))
    pr = matrix(0, (mcmc$burnin + mcmc$rep), K)
    Sigma_inv = NULL
    Sigma1_inv_phi = NULL
    eta = rdirichlet(1, prior$e0)
    p0 = log(eta)
    n0 = prior$n0
    sig02 = prior$sig02
    phi0 = prior$phi0
    Sigma0_inv = kron(diag(K), solve(prior$Sigma0_phi))
    n0sig02 = n0 * sig02
    if (is.null(start.values)) {
        init = kmeans(t(data), K)
        center = init$center
        for (k in 1:K) {
            mfit = arima(center[k, ], order = c(p, 0, 0), method = "CSS")
            phi[k, 1] = mfit$coef[p + 1]
            phi[k, 2:(p + 1)] = mfit$coef[1:p]
        }
        sigma2[1:K] = rigamma(1, n0/2, n0sig02/2)
    }
    else {
        phi[1:K, ] = start.values$phi
        sigma2[1:K] = start.values$sigma2
    }
    for (run in 2:(mcmc$burnin + mcmc$rep)) {
        lik = matrix(0, K, N)
        post = matrix(0, K, N)
        for (n in 1:N) {
            y = data[T:(p + 1), n]
            X = rep(1, T - p)
            for (q in 1:p) {
                X = cbind(X, data[(T - q):(p + 1 - q), n])
            }
            for (k in 1:K) {
                lik[k, n] = sum(dnorm(y, X \%*\% phi[(run - 2) * 
                  K + k, ], sqrt(sigma2[(run - 2) * K + k]), 
                  log = TRUE))
                post[k, n] = p0[k] + lik[k, n]
            }
            post[, n] = post[, n] - max(post[, n])
            s[n] = sample(1:K, size = 1, prob = exp(post[, n]), 
                replace = TRUE)
        }
        ss = unique(s)
        a = rep(0, N)
        k = 1
        while ((length(ss) < K) && (k <= K)) {
            if (sum(ss == k) == 0) {
                s1 = sort.list(post[k, a == 0], decreasing = TRUE)[1:floor(N/K)]
                s[s1] = k
                a[s1] = 1
            }
            ss = unique(s)
            k = sample(1:K, 1)
        }
        phi1 = rep(0, (p + 1) * K)
        Sigma_inv[[run]] = matrix(0, (p + 1) * K, (p + 1) * K)
        Sigma1_inv = matrix(0, (p + 1) * K, (p + 1) * K)
        for (k in 1:K) {
            data_pool = data[, s == k]
            n1 = n0 + (T - p)
            n1sig12 = n0sig02
            pooldim = sum(s == k)
            data_pool = matrix(data_pool, T, pooldim)
            y = c(data_pool[T:(p + 1), ])
            X = rep(1, T - p)
            for (q in 1:p) {
                X = cbind(X, data_pool[(T - q):(p + 1 - q), 1])
            }
            if (pooldim > 1) {
                for (m in 2:pooldim) {
                  Xm = rep(1, T - p)
                  for (q in 1:p) {
                    Xm = cbind(Xm, data_pool[(T - q):(p + 1 - 
                      q), m])
                  }
                  X = rbind(X, Xm)
                  n1 = n1 + (T - p)
                }
            }
            n1sig12 = n1sig12 + t(y - X \%*\% phi[((run - 2) * 
                K + k), ]) \%*\% (y - X \%*\% phi[((run - 2) * K + 
                k), ])
            sigma2[(run - 1) * K + k] = rigamma(1, n1/2, n1sig12/2)
            Sigma1_inv[((k - 1) * (p + 1) + 1):(k * (p + 1)), 
                ((k - 1) * (p + 1) + 1):(k * (p + 1))] = Sigma0_inv[((k - 
                1) * (p + 1) + 1):(k * (p + 1)), ((k - 1) * (p + 
                1) + 1):(k * (p + 1))] + t(X) \%*\% X/sigma2[(run - 
                1) * K + k]
            Sigma1_inv_phi[[k]] = Sigma0_inv[((k - 1) * (p + 
                1) + 1):(k * (p + 1)), ((k - 1) * (p + 1) + 1):(k * 
                (p + 1))] \%*\% phi0 + t(X) \%*\% y/sigma2[(run - 
                1) * K + k]
            phi1[((k - 1) * (p + 1) + 1):(k * (p + 1))] = solve(Sigma1_inv[((k - 
                1) * (p + 1) + 1):(k * (p + 1)), ((k - 1) * (p + 
                1) + 1):(k * (p + 1))]) \%*\% Sigma1_inv_phi[[k]]
            phi[((run - 1) * K + k), ] = rmvnorm(1, phi1[((k - 
                1) * (p + 1) + 1):(k * (p + 1))], solve(Sigma1_inv[((k - 
                1) * (p + 1) + 1):(k * (p + 1)), ((k - 1) * (p + 
                1) + 1):(k * (p + 1))]))
            Sigma_inv[[run]][((k - 1) * (p + 1) + 1):(k * (p + 
                1)), ((k - 1) * (p + 1) + 1):(k * (p + 1))] = Sigma1_inv[((k - 
                1) * (p + 1) + 1):(k * (p + 1)), ((k - 1) * (p + 
                1) + 1):(k * (p + 1))]
        }
        S[run, ] = s
    }
    phiGibbs = phi[-c(1:(K * mcmc$burnin)), ]
    sigma2Gibbs = sigma2[-c(1:(K * mcmc$burnin))]
    SGibbs = S[-c(1:mcmc$burnin), ]
    km = kmeans(cbind(phiGibbs, sqrt(sigma2Gibbs)), K, nstart = 25)
    estim = km$centers
    cluster = matrix(0, mcmc$rep, N)
    pr = matrix(0, mcmc$rep, K)
    SigmaGibbs = NULL
    for (i in 1:mcmc$rep) {
        switch = sort.list(km$cluster[((i - 1) * K + 1):(i * 
            K)], decreasing = FALSE)
        phitemp = phiGibbs[((i - 1) * K + 1):(i * K), ]
        phiGibbs[((i - 1) * K + 1):(i * K), ] = phitemp[switch, 
            ]
        sigma2temp = sigma2Gibbs[((i - 1) * K + 1):(i * K)]
        sigma2Gibbs[((i - 1) * K + 1):(i * K)] = sigma2temp[switch]
        SigmaGibbs[[i]] = solve(Sigma_inv[[i + mcmc$burnin]])
        for (j in 1:N) {
            cluster[i, j] = km$cluster[(i - 1) * K + SGibbs[i, 
                j]]
        }
    }
    Sigmahat = matrix(0, (p + 1) * K, (p + 1) * K)
    for (i in 1:mcmc$rep) {
        Sigmatemp = SigmaGibbs[[i]]
        for (k in 1:K) {
            pr[i, k] = sum(cluster[i, ] == k)/N
            SigmaGibbs[[i]][((k - 1) * (p + 1) + 1):(k * (p + 
                1)), ((k - 1) * (p + 1) + 1):(k * (p + 1))] = Sigmatemp[(switch[k] * 
                (p + 1) - p):(switch[k] * (p + 1)), (switch[k] * 
                (p + 1) - p):(switch[k] * (p + 1))]
        }
        Sigmahat = Sigmahat + SigmaGibbs[[i]]
    }
    prhat = apply(pr, 2, mean)
    Sigmahat = Sigmahat/mcmc$rep
    prob = matrix(0, K, N)
    for (k in 1:K) {
        prob[k, ] = apply(cluster == k, 2, mean)
    }
    logp = rep(0, mcmc$rep)
    for (i in 1:mcmc$rep) {
        if (min(pr[i, ]) > 0) {
            L = matrix(0, N, K)
            for (j in 1:N) {
                y = data[T:(p + 1), j]
                X = rep(1, T - p)
                for (q in 1:p) {
                  X = cbind(X, data[(T - q):(p + 1 - q), j])
                }
                for (k in 1:K) {
                  L[j, k] = log(pr[i, k]) + sum(dnorm(y, X \%*\% 
                    phiGibbs[(i - 1) * K + k, ], sqrt(sigma2Gibbs[(i - 
                    1) * K + k]), log = TRUE))
                }
            }
            Lmax = max(L)
            Ld = L - Lmax
            logp[i] = logp[i] + N * Lmax + sum(log(apply(exp(Ld), 
                1, sum)))
        }
    }
    logp = logp[logp > 0]
    A = max(logp)
    logphat = log(1/mean(1/exp(logp - A))) + A
    rn = rep(0, K)
    for (i in 1:K) {
        rn[i] = paste("Cluster", i)
    }
    rownames(estim) = rn
    cn = rep(0, p + 2)
    cn[p + 2] = paste("sigma")
    for (i in 1:(p + 1)) {
        cn[i] = paste("phi", i - 1)
    }
    colnames(estim) = cn
    rn = rep(0, K)
    for (i in 1:K) {
        rn[i] = paste("Cluster", i)
    }
    rownames(prob) = rn
    cn = rep(0, N)
    for (i in 1:N) {
        cn[i] = colname[i]
    }
    colnames(prob) = cn
    cls = vector("list", K)
    for (k in 1:K) {
        if (length(which(prob[k, ] == apply(prob, 2, max))) > 
            0) {
            cls[[k]] = colname[which(prob[k, ] == apply(prob, 
                2, max))]
        }
    }
    clsn = rep(0, N)
    for (i in 1:N) {
        clsn[i] = which(prob[, i] == apply(prob, 2, max)[i])[1]
    }
    cat("\n Estimation for Cluster Parameters:\n")
    cat(" Number of Clusters: K=", K, "\n")
    cat(" Number of Lags in AR model: p=", p, "\n")
    print(estim, digits = 5)
    cat("\n Classification Probabilities:\n")
    print(t(prob), digits = 3)
    cat("\n Classification:\n")
    for (k in 1:K) {
        nk1 = length(cls[[k]])
        if (nk1 > 0) {
            cat("Cluster ", k, ":\n")
            cat("Number of members: ", nk1, "\n")
            cat(cls[[k]], "\n")
        }
    }
    cat("\n Marginal LogLikelihood:", logphat, "\n")
    result = list(estim = estim, prob = prob, pr = pr, cluster = cls, 
        clsindex = clsn, phiGibbs = phiGibbs, sigma2Gibbs = sigma2Gibbs, 
        SigmaInv_Gibbs = Sigma_inv, S_draw = SGibbs, MarginalLik = logphat)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
