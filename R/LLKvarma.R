#' Likelihood function for VARMA
LLKvarma <- function(par, zt, p, q, include.mean, fixed) {
  nT <- dim(zt)[1]
  k <- dim(zt)[2]
  pqmax <- max(p, q)
  kp <- k * p
  kq <- k * q
  beta <- NULL
  ist <- 0
  icnt <- 0
  Ph0 <- rep(0, k)
  if (include.mean) {
    idx <- seq_len(k)[fixed[1, ] == 1]
    icnt <- length(idx)
    if (icnt > 0) {
      Ph0[idx] <- par[1:icnt]
    }
    ist <- 1
    beta <- rbind(beta, Ph0)
  }
  PH <- matrix(0, kp, k)
  if (p > 0) {
    for (j in 1:k) {
      idx <- seq_len(kp)[fixed[(ist + 1):(ist + kp), j] == 1]
      jdx <- length(idx)
      if (jdx > 0) {
        PH[idx, j] <- par[(icnt + 1):(icnt + jdx)]
        icnt <- icnt + jdx
      }
    }
    ist <- ist + kp
    beta <- rbind(beta, PH)
  }
  TH <- matrix(0, kq, k)
  if (q > 0) {
    for (j in 1:k) {
      idx <- c(1:kq)[fixed[(ist + 1):(ist + kq), j] == 1]
      jdx <- length(idx)
      if (jdx > 0) {
        TH[idx, j] <- par[(icnt + 1):(icnt + jdx)]
        icnt <- icnt + jdx
      }
    }
    beta <- rbind(beta, TH)
  }

  # Calculate the residuals
  at <- varmaResiduals(zt, Ph0, PH, TH, p, q, include.mean, beta)

  sig <- t(at) %*% at / (nT - pqmax)
  ll <- dmvnorm(at, rep(0, k), sig)
  LLKvarma <- -sum(log(ll))
  LLKvarma
}
