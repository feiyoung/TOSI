# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

cor.mat <- function (p, rho, type = "toeplitz")
{
  mat <- diag(p)
  if(p == 1) return(mat)
  if (type == "toeplitz") {
    for (i in 2:p) {
      for (j in 1:i) {
        mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
      }
    }
  }
  if (type == "identity") {
    mat[mat == 0] <- rho
  }
  return(mat)
}


# Factor model ------------------------------------------------------------


gendata_Fac <- function(n, p, seed=1, q=6, pzero= floor(p/4), sigma2=0.1, gamma=1, heter=F, rho=1){
  # sigma2 <- 0.1; heter=F
  set.seed(1)

  psub <- p - pzero
  p1_bar <- floor(psub/q)
  sk <- seq(1.5, 0.3, length=q)
  B0 <- matrix(0, p, q)
  for(k in 1:q){
    B0[(p1_bar *(k-1)+1):(k*p1_bar), k] <- runif(p1_bar) + sk[k]
  }
  nzind <- 1: (p1_bar* q);
  sB = numeric(q)
  for(k in 1:q){
    b0k <- B0[,k]
    sB[k] <- sign(b0k[b0k>0][1])
  }
  B0 <- rho*sapply(1:q, function(k) B0[,k]*sB[k])

  set.seed(seed)
  H <- MASS::mvrnorm(n, rep(0,q), cor.mat(q, 0.5))
  covH <- cov(H)
  eigcH <- eigen(covH)
  if(q==1){
    covH12 <- eigcH$values^(-1/2)
  }else{
    covH12 <- eigcH$vectors %*% diag(eigcH$values^(-1/2))%*%eigcH$vectors
  }

  H0 <- (H-matrix(colMeans(H), n, q, byrow=T)) %*% covH12
  if(heter){
    Sigma <- diag(gamma + runif(p))
  }else{
    Sigma <- sigma2*diag(rep(1,p))
  }

  X <- H0 %*% t(B0) +  MASS::mvrnorm(n, rep(0,p), Sigma)
  return(list(X=X, H0=H0, B0=B0, ind_nz = nzind))
}

signrevise <- function(A1, A2){
  nzid1 <- which(rowSums(A1^2)> 1e-5)[1]
  q <- ncol(A1)
  A <- sapply(1:q, function(k){
    if(sign(A1[nzid1,k]) != sign(A2[nzid1,k]))
      return(-A1[,k])
    return(A1[,k])
  })
  return(A)
}

gsspFactorm <- function(X, q=NULL, lambda1=nrow(X)^(1/4), lambda2=nrow(X)^(1/4)){

  # lambda1 <- 0.05; lambda2 <- 0.5
  n <- nrow(X)
  p <- ncol(X)
  if(p >n){
    svdX <- eigen(X%*%t(X))
    evalues <- svdX$values
    eigrt <- evalues[1:(21-1)]/evalues[2:21]
    if(is.null(q)){
      q <- which.max(eigrt)
    }

    hatF <- as.matrix(svdX$vector[, 1:q] * sqrt(n))
    B2  <- n^(-1)*t(X) %*% hatF
    hB2 <- B2 * matrix(sign(B2[1,]), nrow=p, ncol=q, byrow=T)

    Lvec <- sqrt(rowSums(B2^2))
    w1 <- 1 / Lvec
    w2 <- 1 / abs(B2)
    A <- sign(B2) * pmax(abs(B2) - lambda2 * w2 /(2*n), 0)
    rm(Lvec); Lvec <- sqrt(rowSums(A^2))
    B3 <- matrix(pmax(1 - 1/2*lambda1 * w1/ Lvec, 0), p, q) * A
    rm(B2)
    B2 <- B3
    #  cbind(B3[1:10, 1:3], B0[1:10,1:3])
    Bqr <- qr(B3)
    if(q == 1){
      B4 <- qr.Q(Bqr) %*% sqrt(eigen(t(B3)%*%B3)$values)
    }else{
      B4 <- qr.Q(Bqr) %*% sqrt(diag(eigen(t(B3)%*%B3)$values))
    }
    # B4 <- B3
    B4[abs(B4)< 1e-5] <- 0
    ind <- which(rowSums(B4^2)> 1e-5) # if
    hB <- B4
    hH <- hatF
    if(length(ind)>0){
      jnz <- ind[1]
      for(k in 1:q){
        if(sum(abs(B3[,k]))> 0){
          knz <- which(abs(B3[,k]) > 1e-5)[1]
          hB[,k] <- B4[,k] * sign(B3[knz, k])
          hH[,k] <- hatF[,k] * sign(B3[knz, k])
        }
      }
    }else{
      hB <-B4
      hH <- hatF
    }

  }else{
    svdX <- eigen(t(X)%*%X)
    evalues <- svdX$values
    eigrt <- evalues[1:(21-1)]/evalues[2:21]
    if(is.null(q)){
      q <- which.max(eigrt)
    }
    hB1 <- as.matrix(svdX$vector[, 1:q])

    hH1  <- n^(-1)* X %*% hB1
    svdH <- svd(hH1)
    hH2 <- signrevise(svdH$u *sqrt(n), hH1)
    if(q == 1){
      hB1 <- hB1 %*% svdH$d[1:q] *sqrt(n)
    }else{
      hB1 <- hB1 %*% diag(svdH$d[1:q]) *sqrt(n)
    }


    hB2 <- hB1 * matrix(sign(hB1[1,]), nrow=p, ncol=q, byrow=T)

    Lvec <- sqrt(rowSums(hB2^2)) # lambda2 <- 0
    w1 <- 1 / Lvec
    w2 <- 1 / abs(hB2)
    A <- sign(hB2) * pmax(abs(hB2) - lambda2 * w2 /(2* n), 0)
    rm(Lvec); Lvec <- sqrt(rowSums(A^2))
    B3 <- matrix(pmax(1 - 1/2*lambda1 * w1/ Lvec, 0), p, q) * A
    B2 <- B3
    Bqr <- qr(B3)
    if(q ==1){
      B4 <- qr.Q(Bqr) %*% sqrt(eigen(t(B3)%*%B3)$values)
    }else{
      B4 <- qr.Q(Bqr) %*% sqrt(diag(eigen(t(B3)%*%B3)$values))
    }

    B4 <- signrevise(B4,hB1)
    B4[abs(B4)< 1e-5] <- 0
    ind <- which(rowSums(B4^2)> 1e-5)
    hB <- B4
    hH <- hH2
    if(length(ind)>0){
      jnz <- ind[1]
      for(k in 1:q){
        if(sum(abs(B3[,k]))> 0){
          knz <- which(abs(B3[,k]) > 1e-5)[1]
          hB[,k] <- B4[,k] * sign(B3[knz, k])
          hH[,k] <- hH[,k] * sign(B3[knz, k])
        }
      }
    }else{
      hB <-B4
      hH <- hH2
    }
  }

  res <- list()
  res$hH <- hH
  res$sphB <- hB
  # res$sphB2 <- B2
  res$hB <- hB2
  res$q <- q
  res$propvar <- sum(evalues[1:q]) / sum(evalues)
  res$egvalues <- evalues
  attr(res, 'class') <- 'fac'
  return(res)
}


# Assess function
ccorFun <- function(hH, H){
  q <- ncol(H)
  cancor(hH,H)$cor[q]
}

assessBsFun <- function(hB, B0){
  # row sparsity
  pred <- (rowSums(hB^2)>1e-5)
  true <- (rowSums(B0^2)>1e-5)

  precision <- sum(pred & true) / sum(pred)
  scr <- recall <- sum(pred & true) / sum(true)

  Fmeasure <- 2 * precision * recall / (precision + recall)
  # entry sparsity
  pred2 <- (hB^2>1e-5)
  true2 <- (B0^2 > 1e-5)
  pre2 <- sum(pred2 & true2) / sum(pred2)
  scr2 <- rec2 <- sum(pred2 & true2) / sum(true2)
  fm2 <- 2 * pre2 * rec2 / (pre2 + rec2)
  return(c(rwo_scr=scr, row_fmea=Fmeasure, entry_scr=scr2,
           entry_fmea=fm2,ccorB=ccorFun(hB,B0)))
}
bic.fun1 <- function(X, c1_set, C0=4){
  nlam <- length(c1_set)
  n <- nrow(X); p <- ncol(X)
  lambda1_set <- c1_set * n^(1/4)
  BICv <- numeric(nlam)
  l2pen <- matrix(0, nlam,2)
  for(j in 1:nlam){
    spfac <- gsspFactorm(X, lambda1=lambda1_set[j], lambda2 = 0)
    hnz <- sum(rowSums(spfac$sphB^2)> 1e-5)
    l2 <- min(norm(X- spfac$hH %*% t(spfac$sphB), 'F'), norm(X + spfac$hH %*% t(spfac$sphB), 'F'))
    l2pen[j,] <- c(l2, C0*(n+p)/(n*p)*log(n*p/(n+p))* hnz)
    BICv[j] <- l2pen[j,1] + l2pen[j,2]
  }
  return(list(BIC=BICv, lambda1_set = lambda1_set, l2pen=l2pen) )
}


bic.fun2 <- function(X, c2_set, lambda1.min=0.2*nrow(X)^(1/4),C0=4){
  nlam <- length(c2_set)
  n <- nrow(X); p <- ncol(X)
  lambda2_set <- c2_set * n^(1/4)
  BICv <- numeric(nlam)
  l2pen <- matrix(0, nlam,2)
  for(j in 1:nlam){
    # j<- 1
    spfac <- gsspFactorm(X, lambda2=lambda2_set[j], lambda1 = lambda1.min)
    hnz <- sum(spfac$sphB^2> 1e-5)
    l2 <- min(norm(X- spfac$hH %*% t(spfac$sphB), 'F'), norm(X + spfac$hH %*% t(spfac$sphB), 'F'))
    l2pen[j,] <- c(l2, C0*(n+p)/(n*p)*log(n*p/(n+p))* hnz)
    BICv[j] <- l2pen[j,1] + l2pen[j,2]
  }
  return(list(BIC=BICv, lambda2_set = lambda2_set, l2pen=l2pen) )
}


bic.spfac <- function(X, c1.max= 10, nlamb1=10, C10=4, c2.max=10, nlamb2=10, C20=4){
  # nlambda <- 10; c.max <- 10
  c1_set <- exp(seq(log(c1.max), log(0.001 * c1.max),len = nlamb1 - 1))
  c2_set <- exp(seq(log(c2.max), log(0.001 * c2.max),len = nlamb2 - 1))
  bic1list <- bic.fun1(X, c1_set, C0=C10)
  lambda1.min <- bic1list$lambda1_set[which.min(bic1list$BIC)]
  bic2list <- bic.fun2(X, c2_set,lambda1.min, C0=C20)
  lambda2.min <- bic2list$lambda2_set[which.min(bic2list$BIC)]
  biclist <- list()
  biclist$lambda1.min <- lambda1.min
  biclist$lambda2.min <- lambda2.min
  biclist$bic1 <- cbind(c1=c1_set, lambda1=bic1list$lambda1_set, bic1=bic1list$BIC)
  biclist$bic2 <- cbind(c2=c2_set, lambda2=bic2list$lambda2_set, bic2=bic2list$BIC)
  class(biclist) <- c('pena_info','BIC')
  return(biclist)
}
cv.fun <- function(Xtr, hHtr, Xts, hHts, lambda1_set,lambda2_set){
  n <- nrow(Xtr); q <- ncol(hHtr)
  nts <- nrow(Xts); p <- ncol(Xts)
  B2  <- t(qr.solve(t(hHtr)%*%hHtr)%*% t(hHtr)%*% Xtr)

  Lvec <- sqrt(rowSums(B2^2))
  w1 <- 1 / Lvec
  w2 <- 1 / abs(B2)
  nlam1 <- length(lambda1_set)
  nlam2 <- length(lambda2_set)
  nprod12 <- nlam1*nlam2
  lamMat <- cbind(rep(lambda1_set, each=nlam2),
                  rep(lambda2_set, length=nprod12))
  cVal <- sapply(1: nprod12,  function(j){
    A <- sign(B2) * pmax(abs(B2) - lamMat[j,2] * w2 /(2*n), 0)
    Lvec <- sqrt(rowSums(A^2))
    B3 <- matrix(pmax(1 - 1/2*lamMat[j,1] * w1/ Lvec, 0), p, q) * A
    norm(Xts - hHts %*% t(B3), 'F')^2 / (nts*p)
  } )
  return(cVal)
}

cv.spfac <- function(X, lambda1_set, lambda2_set, nfolds=5){
  spfac <- gsspFactorm(X) # choose q
  n <- nrow(X)
  hH <- spfac$hH
  fold <- ceiling(sample(1:n)/(n + sqrt(.Machine$double.eps)) *
                    nfolds)
  n <- nrow(X); p <- ncol(X)
  # lambda1_set <- c1_set * n^(1/4)
  # lambda2_set <- c2_set * n^(1/4)
  CVs <- sapply(1:nfolds, function(j){
    cv.fun(X[fold!=j,], hH[fold!=j,], X[fold==j,], hH[fold==j,],
           lambda1_set=lambda1_set, lambda2_set=lambda2_set)
  })
  mCVs = apply(CVs, 1, mean)
  nlam1 <- length(lambda1_set)
  nlam2 <- length(lambda2_set)
  nprod12 <- nlam1*nlam2
  lamMat <- cbind(rep(lambda1_set, each=nlam2),
                  rep(lambda2_set, length=nprod12))
  lamcvMat <- cbind(lamMat, mCVs)
  lamcv.min <- lamcvMat[which.min(mCVs),]
  names(lamcv.min) <- c('lambda1.min', 'lambda2.min', 'cv.min')
  return(list(lamcv.min=lamcv.min,  lamcvMat = lamcvMat,
              lambda1_set = lambda1_set, lambda2_set = lambda2_set
  ) )
}


Factorm <- function(X, q=NULL){
  n <- nrow(X)
  p <- ncol(X)
  if(p >n){
    svdX <- eigen(X%*%t(X))
    evalues <- svdX$values
    eigrt <- evalues[1:(21-1)]/evalues[2:21]
    if(is.null(q)){
      q <- which.max(eigrt)
    }

    hatF <- as.matrix(svdX$vector[, 1:q] * sqrt(n))
    B2  <- n^(-1)*t(X) %*% hatF

    sB <- sign(B2[1,])
    hB <- B2 * matrix(sB, nrow=p, ncol=q, byrow=T)
    hH <- sapply(1:q, function(k) hatF[,k]*sign(B2[1,])[k])
  }else{
    svdX <- eigen(t(X)%*%X)
    evalues <- svdX$values
    eigrt <- evalues[1:(21-1)]/evalues[2:21]
    if(is.null(q)){
      q <- which.max(eigrt)
    }
    hB1 <- as.matrix(svdX$vector[, 1:q])

    hH1  <- n^(-1)* X %*% hB1
    svdH <- svd(hH1)
    hH2 <- signrevise(svdH$u *sqrt(n), hH1)
    if(q == 1){
      hB1 <- hB1 %*% svdH$d[1:q] *sqrt(n)
    }else{
      hB1 <- hB1 %*% diag(svdH$d[1:q]) *sqrt(n)
    }

    sB <- sign(hB1[1,])
    hB <- hB1 * matrix(sB, nrow=p, ncol = q, byrow = T)
    hH <- sapply(1:q, function(j) hH2[,j]*sB[j])
  }
  sigma2vec <- colMeans( (X-hH %*% t(hB))^2)

  res <- list()
  res$hH <- hH
  res$hB <- hB
  res$q <- q
  res$sigma2vec <- sigma2vec
  res$propvar <- sum(evalues[1:q]) / sum(evalues)
  res$egvalues <- evalues
  attr(res, 'class') <- 'fac'
  return(res)
}


#  Two-Stage Maximum Row Test method for rows of loading matrix in factor model
TSrowMaxST <- function(X, G1,  alpha=0.05, seed=1, sub.frac=0.5){

  fac <- Factorm(X);  q <- fac$q
  n <- nrow(X)

  ns <- round(n* sub.frac)
  set.seed(seed)
  ids <- sample(n, ns)
  hB <- Factorm(X[ids, ], q=q)$hB
  hBG1Mat <- matrix(hB[G1, ], nrow=length(G1), ncol=q)
  norm1bG1 <- apply(hBG1Mat,1, function(x) sum(abs(x)))
  K1 <- min(1, length(G1) )
  id1 <- order(norm1bG1, decreasing = T)[1:K1]
  G1 <- G1[id1]

  idt <- setdiff(1:n, ids)
  nt <- length(idt)

  fac <- Factorm(X[idt, ], q = q)
  dLam1 <- sqrt(fac$sigma2vec[G1])
  hBG1 <- fac$hB[G1, ]

  maxC1 <- qchisq(1-alpha, q)
  T1 <- nt * sum(hBG1*hBG1)/dLam1^2
  PV <-  1- pchisq(T1, q)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
  row.names(pMat) <- c('chiq_test')


  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Max-test'
  return(pMat)
}

#  Two-Stage Minimum Row Test method for rows of loading matrix in factor model
TSrowMinST <- function(X, G2, alpha=0.05, seed=1, sub.frac=0.5){

  fac <- Factorm(X);   q <- fac$q
  n <- nrow(X)

  ns <- round(n* sub.frac)
  set.seed(seed)
  ids <- sample(n, ns)
  hB <- Factorm(X[ids, ], q=q)$hB
  hBG2Mat <- matrix(hB[G2, ], nrow=length(G2), ncol=q)
  norm1bG2 <- apply(hBG2Mat, 1, function(x) sum(abs(x)))
  K2 <- min(1, length(G2))
  id2 <- order(norm1bG2)[1:K2]
  G2 <- G2[id2]

  idt <- setdiff(1:n, ids)
  nt <- length(idt)



  fac <- Factorm(X[idt,], q=q)
  dLam2 <- sqrt(fac$sigma2vec[G2])
  hBG2 <- fac$hB[G2,]

  minC2 <- qchisq(1-alpha, q)
  R2 <- nt * sum(hBG2*hBG2)/dLam2^2
  PV <-  1- pchisq(R2, q)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(minC2, R2, R2 > minC2, PV)
  row.names(pMat) <- c('chiq_test')


  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Min-test'
  return(pMat)
}

# Transform bi-index  to single index for  matrix.
indMat2vecFun <- function(S, nrow){
  ns <- nrow(S)
  sapply(1:ns, function(j) S[j,1]-1 + (S[j,2]-1)*nrow + 1)
}
## Transform  single index  to bi-index for  matrix.
indvec2matFun <- function(vec, nrow){
  nvec <- length(vec)
  S <- sapply(1:nvec, function(j) {
    j1 <- (vec[j]-1) %% nrow +1
    k1 <- floor((vec[j]-1)/ nrow) + 1
    return(c(j1,k1))
  })
  return(t(S))
}

#  Two-Stage Maximum Entry Test method for rows of loading matrix in factor model
TSentryMaxST <- function(X, S1,  alpha=0.05, seed=1, sub.frac=0.5, q=NULL){

  if(!is.matrix(S1)) S1 <- matrix(S1, 1,2)
  fac <- Factorm(X);
  if(is.null(q)) q <- fac$q
  n <- nrow(X); p <- ncol(X)

  ns <- round(n* sub.frac)
  set.seed(seed)
  ids <- sample(n, ns)
  hB <- Factorm(X[ids, ], q=q)$hB
  S1vec <- indMat2vecFun(S1, p)
  hBG1vec <- hB[S1vec]
  K1 <- min(1, length(S1vec))
  id1 <- order(abs(hBG1vec), decreasing = T)[1:K1]
  G1 <- S1vec[id1]
  #indvec2matFun(G1, p); datlist1$B0[G1]

  idt <- setdiff(1:n, ids)
  nt <- n - ns

  fac <- Factorm(X[idt, ], q = q)
  dLam1 <- sqrt(fac$sigma2vec[S1[id1,1]])
  hBG1 <- fac$hB[G1]
  maxC1 <- qchisq(1-alpha, 1)
  T1 <- nt * sum(hBG1*hBG1)/dLam1^2
  PV <-  1- pchisq(T1, 1)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
  row.names(pMat) <- c('chiq_test')

  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Max-test'
  return(pMat)
}

# Two-Stage Minimum Entry Test method for rows of loading matrix in factor model
TSentryMinST <- function(X, S2,  alpha=0.05, seed=1, sub.frac=0.5, q= NULL){

  if(!is.matrix(S2)) S2 <- matrix(S2, 1,2)
  fac <- Factorm(X);
  if(is.null(q)) q <- fac$q
  n <- nrow(X); p <- ncol(X)

  ns <- round(n* sub.frac)
  set.seed(seed)
  ids <- sample(n, ns)
  hB <- Factorm(X[ids, ], q=q)$hB
  S2vec <- indMat2vecFun(S2, p)
  hBG1vec <- hB[S2vec]
  K1 <- min(1, length(S2vec))
  id1 <- order(abs(hBG1vec), decreasing = F)[1:K1]
  G1 <- S2vec[id1]

  idt <- setdiff(1:n, ids)
  nt <- n - ns

  fac <- Factorm(X[idt, ], q = q)
  dLam1 <- sqrt(fac$sigma2vec[S2[id1,1]])
  hBG1 <- fac$hB[G1]

  maxC1 <- qchisq(1-alpha, 1)
  T1 <- nt * sum(hBG1*hBG1)/dLam1^2
  PV <-  1- pchisq(T1, 1)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
  row.names(pMat) <- c('chiq_test')


  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Max-test'
  return(pMat)
}

# Multi-split Two-Stage Maximum Row Test method for rows of loading matrix in factor model
MultiTSrowMaxST <- function(X, G1,  alpha=0.05, Nsplit= 5, sub.frac=0.5){

  fac <- Factorm(X);  q <- fac$q
  n <- nrow(X)
  ns <- round(n* sub.frac)
  T1vec <- numeric(Nsplit)
  for(im in 1:Nsplit){
    # im <- 3
    set.seed(im)
    ids <- sample(n, ns)
    hB <- Factorm(X[ids, ], q=q)$hB
    hBG1Mat <- matrix(hB[G1, ], nrow=length(G1), ncol=q)
    norm1bG1 <- apply(hBG1Mat,1, function(x) sum(abs(x)))
    K1 <- min(1, length(G1) )
    id1 <- order(norm1bG1, decreasing = T)[1:K1]
    G11 <- G1[id1]

    idt <- setdiff(1:n, ids)
    nt <- length(idt)
    fac <- Factorm(X[idt, ], q = q)
    dLam1 <- sqrt(fac$sigma2vec[G11])
    hBG1 <- fac$hB[G11, ]
    T1vec[im] <- nt*sum(hBG1*hBG1)/dLam1^2
  }

  T1 <- median(T1vec)
  maxC1 <- qchisq(1-alpha, q)

  PV <-  1- pchisq(T1, q)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
  row.names(pMat) <- c('chiq_test')


  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Max-test'
  return(pMat)
}

# Multi-split Two-Stage Minimum Row Test method for rows of loading matrix in factor model
MultiTSrowMinST <- function(X, G2,  alpha=0.05, Nsplit= 5, sub.frac=0.5){

  fac <- Factorm(X);  q <- fac$q
  n <- nrow(X)
  ns <- round(n* sub.frac)
  T1vec <- numeric(Nsplit)
  for(im in 1:Nsplit){
    # im <- 1
    set.seed(im)
    ids <- sample(n, ns)
    hB <- Factorm(X[ids, ], q=q)$hB
    hBG1Mat <- matrix(hB[G2, ], nrow=length(G2), ncol=q)
    norm1bG1 <- apply(hBG1Mat,1, function(x) sum(abs(x)))
    K1 <- min(1, length(G2) )
    id1 <- order(norm1bG1, decreasing = F)[1:K1]
    G21 <- G2[id1]

    idt <- setdiff(1:n, ids)
    nt <- length(idt)



    fac <- Factorm(X[idt, ], q = q)
    dLam1 <- sqrt(fac$sigma2vec[G21])
    hBG1 <- fac$hB[G21, ]
    T1vec[im] <- nt*sum(hBG1*hBG1)/dLam1^2
  }

  T1 <- mean(T1vec)
  minC1 <- qchisq(1-alpha, q)

  PV <-  1- pchisq(T1, q)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(minC1, T1, T1 > minC1, PV)
  row.names(pMat) <- c('chiq_test')

  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'min-test'
  return(pMat)
}

# Mean testing ------------------------------------------------------------
gendata_Mean <- function(n, p, s0= floor(p/2), seed=1, rho= 1, tau=1){
  mu <- rep(0, p)
  set.seed(1)
  mu[1:s0] <- runif(s0)* rho
  set.seed(seed)
  X <- mvrnorm(n=n, mu=mu, Sigma = tau*cor.mat(p, rho=0.5))
  return(list(X=X, mu=mu, p0=s0))
}

MeanMax <- function(X, G1, alpha=0.05,frac.size=0.5, seed=1, standardized=F){
  n <- nrow(X)

  if(length(G1) <= 10){
    frac.size <- 0.1
  }
  ns <- round(n*frac.size)
  set.seed(seed)
  ids <- sample(n, ns)
  hmu <- colMeans(X[ids,])

  abs_muG <- abs(hmu[G1])
  if(standardized){
    std <- apply(as.matrix(X[ids, G1], ncol=length(G1)), 2, sd)
    abs_muG <- abs_muG / std
  }
  K <- min(1, length(G1) )
  id1 <- order(abs_muG, decreasing = T)[1:K]
  test.set <- G1[id1]
  tsXid <- setdiff(1:n, ids)
  nt <- length(tsXid)
  X1 <- X[tsXid, test.set]

  hmu <- mean(X1)
  hsigma2 <- var(X1)


  maxC1 <- qchisq(1-alpha, 1)
  T1 <- nt * hmu^2/ hsigma2
  PV <-  1- pchisq(T1, 1)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
  row.names(pMat) <- c('chiq_test')


  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Max-test'
  return(pMat)
}

MeanMin <- function(X, G2,  alpha=0.05,frac.size=0.5, seed=1, standardized=F){
  n <- nrow(X)
  if(length(G2) <= 10){
    frac.size = 0.1
  }
  ns <- round(n*frac.size)
  set.seed(seed)
  ids <- sample(n, ns)
  hmu <- colMeans(X[ids,])
  abs_muG <- abs(hmu[G2])
  if(standardized){
    std <- apply(as.matrix(X[ids, G2], ncol=length(G2)), 2, sd)
    abs_muG <- abs_muG / std
  }
  K <- min(1, length(G2) )
  id1 <- order(abs_muG)[1:K]
  test.set <- G2[id1]
  tsXid <- setdiff(1:n, ids)
  nt <- length(tsXid)
  X1 <- X[tsXid, test.set]

  hmu <- mean(X1)
  hsigma2 <- var(X1)


  maxC1 <- qchisq(1-alpha, 1)
  T1 <- nt * hmu^2/ hsigma2
  PV <-  1- pchisq(T1, 1)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
  row.names(pMat) <- c('chiq_test')


  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Min-test'
  return(pMat)
}
MultiMeanMax <- function(X, G1, Nsplit = 5, alpha=0.05,frac.size=0.5, seed=1, standardized=F){
  n <- nrow(X)
  test.set <- G1
  if(length(test.set)<=10){
    frac.size <- 0.1
  }
  ns <- round(n*frac.size)
  Tvec <- numeric(Nsplit)
  for(im in 1: Nsplit){

    set.seed(im+seed)
    ids <- sample(n, ns)
    hmu <- colMeans(X[ids,])
    abs_muG <- abs(hmu[test.set])
    if(standardized){
      std <- apply(as.matrix(X[ids, test.set], ncol=length(test.set)), 2, sd)
      abs_muG <- abs_muG / std
    }
    K <- min(1, length(test.set))
    id1 <- order(abs_muG, decreasing = T)[1:K]
    test.set1 <- test.set[id1]
    its <- setdiff(1:n, ids)
    n0 <- length(its)

    hmu <- mean(X[its, test.set1])
    hsd <- sd(X[its, test.set])
    Tvec[im] <- sqrt(n0) * abs(hmu / hsd)
  }
  T1 <- median(Tvec)
  minC1 <- qnorm(1-alpha/2)

  PV <-  2*(1- pnorm(T1))
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(minC1, T1, T1 > minC1, PV)
  row.names(pMat) <- c('Z_test')

  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'max-test'
  return(pMat)
}
MultiMeanMin <- function(X, G2, Nsplit = 5, alpha=0.05,frac.size=0.5, seed=1,standardized=F){
  n <- nrow(X)
  test.set <- G2
  if(length(test.set) <= 10){
    frac.size = 0.1
  }
  ns <- round(n*frac.size)
  Tvec <- numeric(Nsplit)
  for(im in 1: Nsplit){

    set.seed(im+seed)
    ids <- sample(n, ns)
    hmu <- colMeans(X[ids,])
    abs_muG <- abs(hmu[test.set])
    if(standardized){
      std <- apply(as.matrix(X[ids, test.set], ncol=length(test.set)), 2, sd)
      abs_muG <- abs_muG / std
    }
    K <- min(1, length(test.set))
    id1 <- order(abs_muG, decreasing = F)[1:K]
    test.set1 <- test.set[id1]
    its <- setdiff(1:n, ids)
    n0 <- length(its)

    hmu <- mean(X[its, test.set1])
    hsd <- sd(X[its, test.set])
    Tvec[im] <- sqrt(n0) * abs(hmu / hsd)
  }
  T1 <- median(Tvec)
  minC1 <- qnorm(1-alpha/2)

  PV <-  2*(1- pnorm(T1))
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(minC1, T1, T1 > minC1, PV)
  row.names(pMat) <- c('Z_test')

  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'min-test'
  return(pMat)
}

# Regssion testing --------------------------------------------------------
gendata_Reg <- function(n=100, p = 20, s0=5, rho=1, seed=1){
  set.seed(1)
  beta <- rep(0,p)
  beta[1:s0] <- runif(s0,0,2) *rho
  set.seed(seed)
  Sigma <- matrix(NA, p, p)
  for (i in 1:p) Sigma[i,] <- 0.8^(abs(i-(1:p)))
  X <- matrix(rnorm(n*p), n, p)
  X <- t(t(chol(Sigma))%*%t(X))

  # Only first three coefficients are not zeros
  Y <- X %*%beta+rt(n,4)/sqrt(2)
  return(list(Y=Y, X=X, beta0=beta, index_nz=1:s0))
}

RegMax <-  function(X, Y,  G1, alpha=0.05, seed=1, sub.frac=0.5, standardized=F) {
  require(glmnet)
  require(hdi)
  require(SIS)
  require(scalreg)
  n <- dim(X)[1]
  p <- dim(X)[2]
  # Devide sample into two parts
  if(length(G1) <10){
    sub.frac <- 0.2
  }
  n1 <- floor(sub.frac*n)
  n0 <- n-floor(n1)
  set.seed(seed)
  S1 <- sample(1:n, n1, replace=FALSE)
  X.sub <- X[S1,]
  if(standardized){
    X.sub <- scale(X.sub)
  }

  Y.sub <- Y[S1]
  cvfit <- cv.glmnet(X.sub, Y.sub, intercept=FALSE)
  cf <- as.numeric(coef(cvfit, s="lambda.min"))[-1]
  # cf <- as.numeric(glmnet(X.sub, Y.sub, lambda=1e-6)$beta)
  cf_testset <- cf[G1]
  # sort the coefficients
  K <- min(1, length(G1) ) # only caputure the maximum coefs.
  id_test.set <- order(abs(cf_testset), decreasing = T)[1:K]
  test.set <- G1[id_test.set]
  tsXid <- setdiff(1:n, S1)
  Xts <- X[tsXid, ]; n0 <- nrow(Xts)
  Yts <- Y[tsXid]


  score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")
  node <- score.nodewiselasso(Xts, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
                              parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
  Theta <- node$out
  Gram<- t(Xts)%*%Xts/n0

  sreg <- scalreg::scalreg(Xts,Yts)
  beta.hat <- sreg$coefficients
  sigma.sq <- sum((Yts-Xts%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))

  index <- test.set

  Omega <- (t(Theta[,index])%*%Gram%*%Theta[,index])*sigma.sq
  beta.db <- beta.hat[index]+Theta[index,]%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
  T1 <- n*beta.db^2 / Omega;

  maxC1 <- qchisq(1-alpha, 1)
  PV <-  1- pchisq(T1, 1)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
  row.names(pMat) <- c('chiq_test')


  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Max-test'
  return(pMat)
}


RegMin <-  function(X, Y,  G2, alpha=0.05, seed=1, sub.frac=0.5, standardized=F) {
  require(glmnet)
  require(hdi)
  require(SIS)
  require(scalreg)
  n <- dim(X)[1]
  p <- dim(X)[2]
  # Devide sample into two parts
  if(length(G2) <10){
    sub.frac <- 0.2
  }
  n1 <- floor(sub.frac*n)
  n0 <- n-floor(n1)
  set.seed(seed)
  S1 <- sample(1:n, n1, replace=FALSE)
  X.sub <- X[S1,]
  if(standardized){
    X.sub <- scale(X.sub)
  }
  Y.sub <- Y[S1]
  lambda <- 1e-7
  glmnet1 <- glmnet(X.sub, Y.sub, lambda=lambda)
  cf <- coef(glmnet1)[-1]
  cf_testset <- cf[G2]
  # sort the coefficients
  K <- min(1, length(G2) ) # only caputure the minimum coef.
  id_test.set <- order(abs(cf_testset))[1:K]
  test.set <- G2[id_test.set]
  tsXid <- setdiff(1:n, S1)
  Xts <- X[tsXid, ]; n0 <- nrow(Xts)
  Yts <- Y[tsXid]


  score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")
  node <- score.nodewiselasso(Xts, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
                              parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
  Theta <- node$out
  Gram<- t(Xts)%*%Xts/n0

  sreg <- scalreg::scalreg(Xts,Yts)
  beta.hat <- sreg$coefficients
  sigma.sq <- sum((Yts-Xts%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))

  index <- test.set

  Omega <- (t(Theta[,index])%*%Gram%*%Theta[,index])*sigma.sq
  beta.db <- beta.hat[index]+Theta[index,]%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
  T1 <- n * beta.db^2 / Omega;

  maxC1 <- qchisq(1-alpha, 1)
  PV <-  1- pchisq(T1, 1)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
  row.names(pMat) <- c('chiq_test')


  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Min-test'
  return(pMat)
}

linReg_TOSI <- function(X, Y, alpha=0.05, seed=1, sub.frac=0.3, standardized=F){
  library(glmnet)
  p <- ncol(X)
  cvlist <- cv.glmnet(X,Y)

  id_lambdamin <- which.min(cvlist$cvm)
  glmlist <- glmnet(X, Y,lambda=cvlist$lambda[id_lambdamin])
  index_nz <- which(glmlist$beta !=0)
  stmin1 <- RegMin(X, Y,index_nz, alpha, seed, sub.frac, standardized);
  Tmin <- stmin1[1,4] > alpha
  index_zero <- setdiff(1:p, index_nz)
  stmax1 <- RegMax(X, Y,  index_zero, alpha, seed, sub.frac, standardized);
  Tmax <-  stmax1[1,4] > alpha
  # Seek the optimal lambda by Dichotomy optimization
  leftside <- 0; rightside <- id_lambdamin
  k <- 1
  while(Tmin ==T || Tmax==F){
    if(k > 5) break
    id_lambda <- round((leftside + rightside)/2)

    glmlist <- glmnet(X, Y,lambda=cvlist$lambda[id_lambda])
    index_nz <- which(glmlist$beta !=0)
    stmin1 <- RegMin(X, Y,index_nz, alpha, seed, sub.frac, standardized);
    Tmin <- stmin1[1,4] > alpha
    if(Tmin){
      rightside <- id_lambda
    }
    if(!Tmin){
      leftside <- id_lambda
      index_zero <- setdiff(1:p, index_nz)
      stmax1 <- RegMax(X, Y,  index_zero, alpha, seed, sub.frac, standardized);
      Tmax <-  stmax1[1,4] > alpha
    }
    cat('Searching penalty parameter for', k, '-th round\n')
    k <- k +1
  }
  Xind <- X[,index_nz]
  datXY <- data.frame(Y, Xind)
  colnames(datXY) <- c('y', paste0('x',index_nz))
  lmreg <- lm(y~., data=datXY)
  slmreg <- list()
  slmreg$lmreg <- lmreg
  slmreg$glmnet <- glmlist
  slmreg$Infer <- list(min_test=stmin1, max_test = stmax1)

  return(slmreg)
}


MultiRegMax <- function(X, Y,  G1, Nsplit = 5, alpha=0.05, seed=1, sub.frac=0.5, standardized=F){
  require(glmnet)
  require(hdi)
  require(SIS)
  require(scalreg)
  n <- dim(X)[1]
  p <- dim(X)[2]
  if(length(G1)<10){
    sub.frac <- 0.2
  }
  n1 <- floor(sub.frac*n)
  n0 <- n-floor(n1)

  Tvec <- numeric(Nsplit)
  for(i in 1:Nsplit){
    # Devide sample into two parts
    set.seed(seed+i)
    S1 <- sample(1:n, n1, replace=FALSE)
    X.sub <- X[S1,]
    if(standardized){
      X.sub <- scale(X.sub)
    }
    Y.sub <- Y[S1]
    cvfit <- cv.glmnet(X.sub, Y.sub, intercept=FALSE)
    cf <- as.numeric(coef(cvfit, s="lambda.min"))[-1]
    cf_testset <- cf[G1]
    # sort the coefficients
    K <- min(1, length(G1) ) # only caputure the maximum coefs.
    id_test.set <- order(abs(cf_testset), decreasing = T)[1:K]
    test.set <- G1[id_test.set]
    tsXid <- setdiff(1:n, S1)
    Xts <- X[tsXid, ]; n0 <- nrow(Xts)
    Yts <- Y[tsXid]


    score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")
    node <- score.nodewiselasso(Xts, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
                                parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
    Theta <- node$out
    Gram<- t(Xts)%*%Xts/n0

    sreg <- scalreg::scalreg(Xts,Yts)
    beta.hat <- sreg$coefficients
    sigma.sq <- sum((Yts-Xts%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))

    index <- test.set

    Omega <- (t(Theta[,index])%*%Gram%*%Theta[,index])*sigma.sq
    beta.db <- beta.hat[index]+Theta[index,]%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
    Tvec[i] <- n*beta.db^2 / Omega;
  }
  T1 <- median(Tvec)
  maxC1 <- qchisq(1-alpha, 1)
  PV <-  1- pchisq(T1, 1)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
  row.names(pMat) <- c('chiq_test')


  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Max-test'
  return(pMat)
}


MultiRegMin <- function(X, Y,  G2, Nsplit = 5, alpha=0.05, seed=1, sub.frac=0.5, standardized=F){
  require(glmnet)
  require(hdi)
  require(SIS)
  require(scalreg)
  n <- dim(X)[1]
  p <- dim(X)[2]
  if(length(G2)<10){
    sub.frac <- 0.2
  }
  n1 <- floor(sub.frac*n)
  n0 <- n-floor(n1)
  Tvec <- numeric(Nsplit)
  for(i in 1:Nsplit){
    # Devide sample into two parts
    set.seed(i+seed)
    S1 <- sample(1:n, n1, replace=FALSE)
    X.sub <- X[S1,]
    if(standardized){
      X.sub <- scale(X.sub)
    }
    Y.sub <- Y[S1]
    lambda <- 1e-7
    glmnet1 <- glmnet(X.sub, Y.sub, lambda=lambda)
    cf <- coef(glmnet1)[-1]
    cf_testset <- cf[G2]
    # sort the coefficients
    K <- min(1, length(G2) ) # only caputure the minimum coef.
    id_test.set <- order(abs(cf_testset))[1:K]
    test.set <- G2[id_test.set]
    tsXid <- setdiff(1:n, S1)
    Xts <- X[tsXid, ]; n0 <- nrow(Xts)
    Yts <- Y[tsXid]


    score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")
    node <- score.nodewiselasso(Xts, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
                                parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
    Theta <- node$out
    Gram<- t(Xts)%*%Xts/n0

    sreg <- scalreg::scalreg(Xts,Yts)
    beta.hat <- sreg$coefficients
    sigma.sq <- sum((Yts-Xts%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))

    index <- test.set

    Omega <- (t(Theta[,index])%*%Gram%*%Theta[,index])*sigma.sq
    beta.db <- beta.hat[index]+Theta[index,]%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
    Tvec[i] <- n*beta.db^2 / Omega;
  }
  T1 <- median(Tvec)
  maxC1 <- qchisq(1-alpha, 1)
  PV <-  1- pchisq(T1, 1)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
  row.names(pMat) <- c('chiq_test')


  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Max-test'
  return(pMat)
}

