###Update Date: 2020-04-03
### Written by Wei Liu in 2020-01-16

library(MASS)
# Generate correlation matrix
cor.mat <- function (p, rho, type = "toeplitz") {
  mat <- diag(p)
  if(p == 1) return(mat)
  if (type == "toeplitz"){
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
# Factor model estimation for high dimensional data
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
    B2  <- matrix(n^(-1)*t(X) %*% hatF, ncol=q)
    
    sbMat <- if(q==1){
      sign(B2[1,])
    }else{
      diag(sign(B2[1,])) 
    }
    hB <- B2 %*% sbMat
    hH <- hatF %*% sbMat
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
    if(q == 1){
      hB1 <- hB1 %*% svdH$d[1:q] *sqrt(n)
    }else{
      hB1 <- hB1 %*% diag(svdH$d[1:q]) *sqrt(n)
    }
    
    sbMat <- if(q==1){
      sign(hB1[1,])
    }else{
      diag(sign(hB1[1,])) 
    }
    hB <- hB1 %*% sbMat
    hH <- hH1 %*% sbMat
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
# Generate data with H0 coupling with identifiability condition
gendatass1 <- function(n, p, seed=1, q=6, pzero= floor(p/4), sigma2=0.1, gamma=1, heter=F, rho=1){
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
# Generate data with H0 not normalized form.
gendatass2 <- function(n, p, seed=1, q=6, pzero= floor(p/4), sigma2=0.1, gamma=1, heter=F, rho=1){
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
  H0 <- MASS::mvrnorm(n, rep(0,q), cor.mat(q, 0))
  if(heter){
    Sigma <- diag(gamma + runif(p))
  }else{
    Sigma <- sigma2*diag(rep(1,p))
  }
  
  X <- H0 %*% t(B0) +  MASS::mvrnorm(n, rep(0,p), Sigma)
  return(list(X=X, H0=H0, B0=B0, ind_nz = nzind))
}

#  Two-Stage Maximum Row Test method for rows of loading matrix in factor model
TSrowMaxST <- function(X, G1, q=NULL, alpha=0.05, seed=1, sub.frac=0.5, standardized=F){
  
  if(is.null(q)){
    fac <- Factorm(X);  q <- fac$q
  }
  
  
  n <- nrow(X)
  ns <- round(n* sub.frac)
  set.seed(seed)
  ids <- sample(n, ns)
  fac_test <- Factorm(X[ids, ], q=q)
  hB <- fac_test$hB
  
  hBG1Mat <- matrix(hB[G1, ], nrow=length(G1), ncol=q)
  norm1bG1 <- apply(hBG1Mat,1, function(x) sum(x^2))
  if(standardized){
    norm1bG1 <-  norm1bG1 / fac_test$sigma2vec[G1]
  }
  K1 <- min(1, length(G1))
  id1 <- order(norm1bG1, decreasing = T)[1:K1]
  G1 <- G1[id1]
  
  idt <- setdiff(1:n, ids)
  nt <- length(idt)
  
  fac <- Factorm(X[idt, ], q = q)
  dLam1 <- fac$sigma2vec[G1]
  hBG1 <- fac$hB[G1, ]
  
  maxC1 <- qchisq(1-alpha, q)
  T1 <- nt * sum(hBG1*hBG1)/dLam1
  PV <-  1- pchisq(T1, q)
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
  row.names(pMat) <- c('chiq_test')
  
  
  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'Max-test'
  return(pMat)
}

#  Two-Stage Minimum Row Test method for rows of loading matrix in factor model
TSrowMinST <- function(X, G2, alpha=0.05, seed=1, sub.frac=0.5, standardized=F){
  
  fac <- Factorm(X);   q <- fac$q
  n <- nrow(X)
  
  ns <- round(n* sub.frac)
  set.seed(seed)
  ids <- sample(n, ns)
  fac_test <- Factorm(X[ids, ], q=q)
  hB <- fac_test$hB
  hBG2Mat <- matrix(hB[G2, ], nrow=length(G2), ncol=q)
  norm1bG2 <- apply(hBG2Mat, 1, function(x) sum(abs(x)))
  if(standardized){
    norm1bG2 <-  norm1bG2 / sqrt(fac_test$sigma2vec[G2])
  }
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

# Transform two-index  to single index for  matrix.
indMat2vecFun <- function(S, nrow){
  ns <- nrow(S)
  sapply(1:ns, function(j) S[j,1]-1 + (S[j,2]-1)*nrow + 1)
} 


# Multi-split Two-Stage Maximum Row Test method for rows of loading matrix in factor model
MultiTSrowMaxST <- function(X, G1, q=NULL,  alpha=0.05, Nsplit= 5, standardized=F,sub.frac=0.5, r=0.5, BH=F, seed=1){
  
  if(is.null(q)){
    fac <- Factorm(X);  q <- fac$q
  }
  n <- nrow(X)
  ns <- round(n* sub.frac)
  Pvec <- numeric(Nsplit)
  for(im in 1:Nsplit){
    # im <- 3
    set.seed(im+seed)
    ids <- sample(n, ns)
    fac_test <- Factorm(X[ids, ], q=q)
    hB <- fac_test$hB
    hBG1Mat <- matrix(hB[G1, ], nrow=length(G1), ncol=q)
    norm1bG1 <- apply(hBG1Mat,1, function(x) sum(x^2))
    if(standardized){
      norm1bG1 <-  norm1bG1 / fac_test$sigma2vec[G1]
    }
    K1 <- min(1, length(G1) )
    id1 <- order(norm1bG1, decreasing = T)[1:K1]
    G11 <- G1[id1]
    
    idt <- setdiff(1:n, ids)
    nt <- length(idt)
    fac <- Factorm(X[idt, ], q = q)
    dLam1 <- fac$sigma2vec[G11]
    hBG1 <- fac$hB[G11, ]
    T1 <- nt*sum(hBG1*hBG1)/dLam1
    Pvec[im] <-  1- pchisq(T1, q)
  }
  
  
  nr <- length(r)
  reject <- numeric(nr)
  for(i in 1:nr){
    gamma <- alpha*r[i]
    if(mean(Pvec <= gamma) >= r[i]){
      reject[i] <- 1
    }
  }
  nme <- as.character(r)
  if(BH){
    Pvec <- p.adjust(Pvec, method="BH")
    gamma <- alpha
    tmp <- 0
    if(mean(Pvec <= gamma) >= 1/Nsplit){
      tmp <- 1
    }
    reject <- c(reject, tmp)
    nme <- c(nme, "BH")
    adj.pval <- min(Pvec)
    attr(reject, "adj.pval") <- adj.pval
  }
  names(reject) <- nme
  return(reject)
}

# Multi-split Two-Stage Minimum Row Test method for rows of loading matrix in factor model
MultiTSrowMinST <- function(X, G2,  q=NULL, alpha=0.05, Nsplit= 5, standardized=T,sub.frac=0.5, r=0.5, BH=F){
  
  if(is.null(q)){
    fac <- Factorm(X);  q <- fac$q
  }
  n <- nrow(X)
  ns <- round(n* sub.frac)
  Pvec <- numeric(Nsplit)
  
  
  for(im in 1:Nsplit){
    # im <- 1
    set.seed(im)
    ids <- sample(n, ns)
    fac_test <- Factorm(X[ids, ], q=q)
    hB <- fac_test$hB
    hBG1Mat <- matrix(hB[G2, ], nrow=length(G2), ncol=q)
    norm1bG1 <- apply(hBG1Mat,1, function(x) sum(x^2))
    if(standardized){
      norm1bG1 <-  norm1bG1 / sqrt(fac_test$sigma2vec[G2])
    }
    K1 <- min(1, length(G2) )
    id1 <- order(norm1bG1, decreasing = F)[1:K1]
    G21 <- G2[id1]
    
    idt <- setdiff(1:n, ids)
    nt <- length(idt)
    
    
    
    fac <- Factorm(X[idt, ], q = q)
    dLam1 <- sqrt(fac$sigma2vec[G21])
    hBG1 <- fac$hB[G21, ]
    
    T1 <- nt*sum(hBG1*hBG1)/dLam1^2
    Pvec[im] <-  1- pchisq(T1, q)
}


  nr <- length(r)
  reject <- numeric(nr)
  for(i in 1:nr){
    gamma <- alpha*r[i]
    if(mean(Pvec <= gamma) >= r[i]){
      reject[i] <- 1
    }
  }
  nme <- as.character(r)
  if(BH){
    Pvec <- p.adjust(Pvec, method="BH")
    gamma <- alpha
    tmp <- 0
    if(mean(Pvec <= gamma) >= 1/Nsplit){
      tmp <- 1
    }
    reject <- c(reject, tmp)
    nme <- c(nme, "BH")
    adj.pval <- min(Pvec)
    attr(reject, "adj.pval") <- adj.pval
  }
  names(reject) <- nme
  return(reject)
}

## Some Summary functions
indvec2matFun <- function(vec, nrow){
  nvec <- length(vec)
  S <- sapply(1:nvec, function(j) {
    j1 <- (vec[j]-1) %% nrow +1
    k1 <- floor((vec[j]-1)/ nrow) + 1
    return(c(j1,k1))
  })
  return(t(S))
}
getTSrowMinST <- function(X, G2list, ...){
  
  nG2 <- length(G2list)
  tmp <- sapply(1:nG2, function(l) TSrowMinST(X, G2=G2list[[l]],...), simplify = "array")
  t(tmp[1,,])
}
getMTSrowMinST <- function(X, G2list, ...){
  
  nG1 <- length(G2list)
  tmp <- sapply(1:nG1, function(l) MultiTSrowMinST(X, G2=G2list[[l]],...), simplify = "array")
  t(tmp[1,,])
}
getTSrowMaxST <- function(X, G1list, ...){
  
  nG1 <- length(G1list)
  tmp <- sapply(1:nG1, function(l) TSrowMaxST(X, G1=G1list[[l]],...), simplify = "array")
  t(tmp[1,,])
}
getMTSrowMaxST <- function(X, G1list, ...){
  
  nG1 <- length(G1list)
  tmp <- sapply(1:nG1, function(l) MultiTSrowMaxST(X, G1=G1list[[l]],...), simplify = "array")
  t(tmp[1,,])
}


getmultiRowMax <- function(X,G1list,...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- MultiTSrowMaxST(X, G1 =G1list[[l]], ...)
    return(lreg)
  })
  if(is.vector(tmp)) names(tmp) <- paste0('G1', 1:nS1)
  if(is.matrix(tmp)) colnames(tmp) <- paste0('G1', 1:nS1)
  return(tmp)
}
getmultiRowMin <- function(X,G2list,...){
  nS1 <- length(G2list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- MultiTSrowMinST(X, G2 =G2list[[l]], ...)
    return(lreg)
  })
  if(is.vector(tmp)) names(tmp) <- paste0('G2', 1:nS1)
  if(is.matrix(tmp)) colnames(tmp) <- paste0('G2', 1:nS1)
  return(tmp)
}



# Mean testing ------------------------------------------------------------



getMeanMax <- function(X,G1list,...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- MeanMax(X, G1 =G1list[[l]], ...)
    return(lreg)
  })
  colnames(tmp) <- paste0('G1', 1:nS1)
  row.names(tmp) <- c('CriticalValue', 'TestStatistic',
                      'reject_status',    ' p-value')
  return(tmp)
}
getmaxMeantest <- function(X,G1list,...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- maxMeanTest(X, test.set  =G1list[[l]], ...)
    return(lreg[[1]])
  })
  colnames(tmp) <- paste0('G1', 1:nS1)
  row.names(tmp) <- c('CriticalValue', 'TestStatistic',
                      'reject_status',    ' p-value')
  return(tmp)
}


getmultiMeanMax2 <- function(X,G1list,...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- MultiMeanMax2(X, test.set =G1list[[l]], ...)
    return(lreg)
  })
  if(is.vector(tmp)) names(tmp) <- paste0('G1', 1:nS1)
  if(is.matrix(tmp)) colnames(tmp) <- paste0('G1', 1:nS1)
  return(tmp)
}


MeanMax <- function(X, G1, alpha=0.05,frac.size=0.5, seed=1, standardized=F){
  n <- nrow(X)
  
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



MultiMeanMax2 <- function(X, test.set, Nsplit = 5, alpha=0.05,frac.size=0.5, standardized=F, r=0.5, BH=F){
  
  ## Conservative inference by combing p-values based on quantile.
  
  n <- nrow(X)
  
  # if(length(test.set)<=10){
  #   frac.size <- 0.1
  # }
  ns <- round(n*frac.size)
  Pvec <- numeric(Nsplit)
  for(im in 1: Nsplit){
    
    set.seed(im)
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
    T1 <- sqrt(n0) * abs(hmu / hsd)
    Pvec[im] <- 2*(1- pnorm(T1))
  }
  nr <- length(r)
  reject <- numeric(nr)
  for(i in 1:nr){
    gamma <- alpha*r[i]
    if(mean(Pvec <= gamma) >= r[i]){
      reject[i] <- 1
    }
  }
  nme <- as.character(r)
  if(BH){
    Pvec <- p.adjust(Pvec, method="BH")
    gamma <- alpha
    tmp <- 0
    if(mean(Pvec <= gamma) >= 1/Nsplit){
      tmp <- 1
    }
    reject <- c(reject, tmp)
    nme <- c(nme, "BH")
    adj.pval <- min(Pvec)
    attr(reject, "adj.pval") <- adj.pval
  }
  names(reject) <- nme
  return(reject)
}


getMeanMin <- function(X,G1list,...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- MeanMin(X, G2 =G1list[[l]], ...)
    return(lreg)
  })
  colnames(tmp) <- paste0('G2', 1:nS1)
  row.names(tmp) <- c('CriticalValue', 'TestStatistic',
                      'reject_status',    ' p-value')
  return(tmp)
}



getmultiMeanMin2 <- function(X,G1list,...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- MultiMeanMin2(X, test.set =G1list[[l]], ...)
    return(lreg)
  })
  if(is.vector(tmp)) names(tmp) <- paste0('G1', 1:nS1)
  if(is.matrix(tmp)) colnames(tmp) <- paste0('G1', 1:nS1)
  return(tmp)
}

MeanMin <- function(X, G2,  alpha=0.05,frac.size=0.5, seed=1, standardized=F){
  n <- nrow(X)
  
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

MultiMeanMin2 <- function(X, test.set, Nsplit = 5, alpha=0.05,frac.size=0.5, standardized=F, r=0.5, BH=F){
  
  ## Conservative inference by combing p-values based on quantile.
  n <- nrow(X)
  
  ns <- round(n*frac.size)
  Pvec <- numeric(Nsplit)
  for(im in 1: Nsplit){
    
    set.seed(im)
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
    T1 <- sqrt(n0) * abs(hmu / hsd)
    Pvec[im] <-  2*(1- pnorm(T1))
  }
  
  nr <- length(r)
  reject <- numeric(nr)
  for(i in 1:nr){
    gamma <- alpha*r[i]
    if(mean(Pvec <= gamma) >= r[i]){
      reject[i] <- 1
    }
  }
  nme <- as.character(r)
  if(BH){
    Pvec <- p.adjust(Pvec, method="BH")
    gamma <- alpha
    tmp <- 0
    if(mean(Pvec <= gamma) >= 1/Nsplit){
      tmp <- 1
    }
    reject <- c(reject, tmp)
    nme <- c(nme, "BH")
    adj.pval <- min(Pvec)
    attr(reject, "adj.pval") <- adj.pval
  }
  names(reject) <- nme
  return(reject)
}




## Chernozhukov et al. (2013) 
maxMeanTest <- function(X, test.set, M= 500, alpha=0.05,frac.size=0.5, seed=1, screening=F){
  n <- nrow(X)
  if(screening){
    ns <- round(n*frac.size)
    set.seed(1)
    ids <- sample(n, ns)
    hmu <- colMeans(X[ids,])
    abs_muG <- abs(hmu[ test.set])
    K <- min(1, length(test.set) )
    id1 <- order(abs_muG, decreasing = T)[1:K]
    test.set <- test.set[id1]
    tsXid <- setdiff(1:n, ids)
    X <- X[tsXid, ]
  }
  hmu <- colMeans(X)
  n0 <- nrow(X)
  ntest <- length(test.set)
  XtsMat <- matrix(X[,test.set], nrow=n0, ncol=ntest)
  omega <- apply(XtsMat, 2, sd)
  set.seed(seed)
  
  E <- matrix(rnorm(n0*M),n0,M)
  cvMat <- t(XtsMat)%*%E /sqrt(n0)
  stat.boot.nst <- apply(cvMat, 2, function(x) max(abs(x)))
  # hist(stat.boot.nst)
  stat.boot.st <- apply(cvMat, 2, function(x) max(abs(x/ omega)))
  # hist(stat.boot.st)
  
  stat.nst <- max(abs(hmu[test.set])) * sqrt(n0)
  stat.st <- max(abs(hmu[test.set]/omega)) * sqrt(n0)
  
  # calculate the p-values
  pv.nst <- 1-ecdf(stat.boot.nst)(stat.nst)
  pv.st <-  1- ecdf(stat.boot.st)(stat.st)
  cri_val.nst <- quantile(stat.boot.nst,1-alpha)
  cri_val.st <- quantile(stat.boot.st,1-alpha)
  if (stat.nst> cri_val.nst) rej.nst <- 1 else rej.nst <-0
  if (stat.st> cri_val.st) rej.st <- 1 else rej.st <- 0
  stvec <- c(cri_val.st,stat.st,rej.st, pv.st)
  names(stvec) <- c("critical value","test statistic","reject", "p-value")
  nstvec <- c(stat.nst,cri_val.nst,rej.nst ,pv.nst)
  names(nstvec) <- c("critical value","test statistic","reject", "p-value")
  result <- list(stvec, nstvec)
  attr(result, 'test') <- 'max-type'
  return(result)
}


# Regssion testing --------------------------------------------------------
score.partialnodewiselasso <- function (x, index_set,wantTheta = FALSE, verbose = FALSE, lambdaseq = "quantile",
          parallel = FALSE, ncores = 8, oldschool = FALSE, lambdatuningfactor = 1,
          cv.verbose = FALSE, do.ZnZ = TRUE){
  
  
  nodewise.getlambdasequence = getFromNamespace("nodewise.getlambdasequence", "hdi")
  nodewise.getlambdasequence.old = getFromNamespace("nodewise.getlambdasequence.old", "hdi")
  cv.nodewise.bestlambda = getFromNamespace("cv.nodewise.bestlambda", "hdi")
  improve.lambda.pick <- getFromNamespace("improve.lambda.pick", "hdi")
  score.getZforlambda <- getFromNamespace("score.getZforlambda", "hdi")
  
  lambdas <- switch(lambdaseq, quantile = nodewise.getlambdasequence(x),
                    linear = nodewise.getlambdasequence.old(x, verbose),
                    stop("invalid 'lambdaseq': ", lambdaseq))
  if (verbose) {
    cat("Using the following lambda values:", lambdas,
        "\n")
  }
  cvlambdas <- cv.nodewise.bestlambda(lambdas = lambdas, x = x,
                                      parallel = parallel, ncores = ncores, oldschool = oldschool,
                                      verbose = cv.verbose)
  if (verbose) {
    cat(paste("lambda.min is", cvlambdas$lambda.min),
        "\n")
    cat(paste("lambda.1se is", cvlambdas$lambda.1se),
        "\n")
  }
  if (do.ZnZ) {
    bestlambda <- improve.lambda.pick(x = x, parallel = parallel,
                                      ncores = ncores, lambdas = lambdas, bestlambda = cvlambdas$lambda.min,
                                      verbose = verbose)
    if (verbose) {
      cat("Doing Z&Z technique for picking lambda\n")
      cat("The new lambda is", bestlambda, "\n")
      cat("In comparison to the cross validation lambda, lambda = c * lambda_cv\n")
      cat("c=", bestlambda/cvlambdas$lambda.min,
          "\n")
    }
  }
  else {
    if (lambdatuningfactor == "lambda.1se") {
      if (verbose)
        cat("lambda.1se used for nodewise tuning\n")
      bestlambda <- cvlambdas$lambda.1se
    }
    else {
      if (verbose)
        cat("lambdatuningfactor used is", lambdatuningfactor,
            "\n")
      bestlambda <- cvlambdas$lambda.min * lambdatuningfactor
    }
  }
  if (verbose) {
    cat("Picked the best lambda:", bestlambda, "\n")
  }
  if (wantTheta) {
    out <- score.getpartialThetaforlambda(x = x, index_set=index_set,lambda = bestlambda,
                                          parallel = parallel, ncores = ncores, oldschool = TRUE,
                                          verbose = verbose)
  }
  else {
    Z <- score.getZforlambda(x = x, lambda = bestlambda,
                             parallel = parallel, ncores = ncores, oldschool = oldschool)
    out <- Z
  }
  return.out <- list(out = out, bestlambda = bestlambda)
  return(return.out)
}
score.getpartialThetaforlambda <- function (x,index_set, lambda, parallel = FALSE, ncores = 8, oldschool = FALSE,
          verbose = FALSE, oldtausq = TRUE){
  message("Calculating Thetahat by doing nodewise regressions and dropping the unpenalized intercept")
  n <- nrow(x)
  p <- ncol(x)
  p1 <- length(index_set)
  C <- matrix(0, p, p1)
  
  T2 <- numeric(p1)
  if (oldschool) {
    
    message("doing getThetaforlambda oldschool")
    for (ii in 1:p1) {
      i <- index_set[ii]
      C[i,ii] <- 1
      glmnetfit <- glmnet(x[, -i], x[, i])
      coeffs <- as.vector(predict(glmnetfit, x[, -i], type = "coefficients",
                                  s = lambda))[-1]
      C[-i, ii] <- -as.vector(coeffs)
      if (oldtausq) {
        T2[ii] <- as.numeric(crossprod(x[, i])/n - x[,
                                                     i] %*% (x[, -i] %*% coeffs)/n)
      }
      else {
        T2[ii] <- as.numeric((x[, i] %*% (x[, i] - predict(glmnetfit,
                                                           x[, -i], s = lambda)))/n)
      }
    }
  }
  else {
    stop("not implemented yet!")
  }
  if(p1==1){
    T2 <- matrix(1/T2, 1,1)
  }else{
    T2 <- diag(1/T2)
  }
  thetahat <- C %*% T2
  return(thetahat)
}
getTheta <- function(Xts, index_set){
  node2 <- score.partialnodewiselasso(Xts, index_set=index_set, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
                                      parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
  Theta <- node2$out
  return(Theta)
}
RegMax <-  function(X, Y,  G1, alpha=0.05, seed=1, sub.frac=0.5, standardized=F) {
    # alpha=0.05; seed=1; sub.frac=0.5; standardized=F
    require(glmnet)
    require(hdi)
    require(SIS)
    require(scalreg)
    n <- dim(X)[1]
    p <- dim(X)[2]
    # Devide sample into two parts
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
    
    
    # score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")
    # node <- score.nodewiselasso(Xts, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
    #                             parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
    # Theta <- node$out
    Theta_index <- getTheta(Xts, index_set = test.set)
    Gram<- t(Xts)%*%Xts/n0
    
    sreg <- scalreg::scalreg(Xts,Yts)
    beta.hat <- sreg$coefficients
    sigma.sq <- sum((Yts-Xts%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))
    
    index <- test.set
    
    # Omega <- (t(Theta[,index])%*%Gram%*%Theta[,index])*sigma.sq
    # beta.db <- beta.hat[index]+Theta[index,]%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
    Omega <- (t(Theta_index)%*%Gram%*%Theta_index)*sigma.sq
    beta.db <- beta.hat[index]+t(Theta_index)%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
    T1 <- n0 * beta.db^2 / Omega;
    
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
    
    
    # score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")
    # node <- score.nodewiselasso(Xts, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
    #                             parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
    # Theta <- node$out
    Theta_index <- getTheta(Xts, index_set = test.set)
    Gram<- t(Xts)%*%Xts/n0
    
    sreg <- scalreg::scalreg(Xts,Yts)
    beta.hat <- sreg$coefficients
    sigma.sq <- sum((Yts-Xts%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))
    
    index <- test.set
    
    # Omega <- (t(Theta[,index])%*%Gram%*%Theta[,index])*sigma.sq
    # beta.db <- beta.hat[index]+Theta[index,]%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
    Omega <- (t(Theta_index)%*%Gram%*%Theta_index)*sigma.sq
    beta.db <- beta.hat[index]+t(Theta_index)%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
    T1 <- n0 * beta.db^2 / Omega;
    
    maxC1 <- qchisq(1-alpha, 1)
    PV <-  1- pchisq(T1, 1)
    pMat <- matrix(0,1,4)
    pMat[1,] <- c(maxC1, T1, T1 > maxC1, PV)
    row.names(pMat) <- c('chiq_test')
    
    
    colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
    class(pMat) <- 'Min-test'
    return(pMat)
  }



multiRegMax <- function(X, Y,  G1, Nsplit = 5, alpha=0.05, seed=1, sub.frac=0.5, standardized=F, r= 0.5, BH=F){
  require(glmnet)
  require(hdi)
  require(SIS)
  require(scalreg)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  n1 <- floor(sub.frac*n) 
  n0 <- n-floor(n1)
  
  Pvec <- numeric(Nsplit)
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
    
    
    # score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")
    # node <- score.nodewiselasso(Xts, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
    #                             parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
    # Theta <- node$out
    Theta_index <- getTheta(Xts, index_set = test.set)
    Gram<- t(Xts)%*%Xts/n0
    
    sreg <- scalreg::scalreg(Xts,Yts)
    beta.hat <- sreg$coefficients
    sigma.sq <- sum((Yts-Xts%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))
    
    index <- test.set
    
    # Omega <- (t(Theta[,index])%*%Gram%*%Theta[,index])*sigma.sq
    # beta.db <- beta.hat[index]+Theta[index,]%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
    Omega <- (t(Theta_index)%*%Gram%*%Theta_index)*sigma.sq
    beta.db <- beta.hat[index]+t(Theta_index)%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
    T1 <- n0 * beta.db^2 / Omega;
    Pvec[i] <-  1- pchisq(T1, 1)
  }
  nr <- length(r)
  reject <- numeric(nr)
  for(i in 1:nr){
    gamma <- alpha*r[i]
    if(mean(Pvec <= gamma) >= r[i]){
      reject[i] <- 1
    }
  }
  nme <- as.character(r)
  if(BH){
    Pvec <- p.adjust(Pvec, method="BH")
    gamma <- alpha
    tmp <- 0
    if(mean(Pvec <= gamma) >= 1/Nsplit){
      tmp <- 1
    }
    
    reject <- c(reject, tmp)
    nme <- c(nme, "BH")
    adj.pval <- min(Pvec)
    attr(reject, "adj.pval") <- adj.pval
  }
  names(reject) <- nme
  return(reject)
}


multiRegMin <- function(X, Y,  G2, Nsplit = 5, alpha=0.05, seed=1, sub.frac=0.5, standardized=F, r=0.5, BH=F){
  require(glmnet)
  require(hdi)
  require(SIS)
  require(scalreg)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  n1 <- floor(sub.frac*n) 
  n0 <- n-floor(n1)
  Pvec <- numeric(Nsplit)
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
    
    
    # score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")
    # node <- score.nodewiselasso(Xts, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
    #                             parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
    # Theta <- node$out
    Theta_index <- getTheta(Xts, index_set = test.set)
    Gram<- t(Xts)%*%Xts/n0
    
    sreg <- scalreg::scalreg(Xts,Yts)
    beta.hat <- sreg$coefficients
    sigma.sq <- sum((Yts-Xts%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))
    
    index <- test.set
    
    # Omega <- (t(Theta[,index])%*%Gram%*%Theta[,index])*sigma.sq
    # beta.db <- beta.hat[index]+Theta[index,]%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
    Omega <- (t(Theta_index)%*%Gram%*%Theta_index)*sigma.sq
    beta.db <- beta.hat[index]+t(Theta_index)%*%t(Xts)%*%(Yts-Xts%*%beta.hat)/n0
    T1 <- n0*beta.db^2 / Omega;
    Pvec[i] <-  1- pchisq(T1, 1)
  }
  nr <- length(r)
  reject <- numeric(nr)
  for(i in 1:nr){
    gamma <- alpha*r[i]
    if(mean(Pvec <= gamma) >= r[i]){
      reject[i] <- 1
    }
  }
  nme <- as.character(r)
  if(BH){
    Pvec <- p.adjust(Pvec, method="BH")
    gamma <- alpha
    tmp <- 0
    if(mean(Pvec <= gamma) >= 1/Nsplit){
      tmp <- 1
    }
    reject <- c(reject, tmp)
    nme <- c(nme, "BH")
    adj.pval <- min(Pvec)
    attr(reject, "adj.pval") <- adj.pval
  }
  names(reject) <- nme
  return(reject)
}

regST_max <- function(X, Y,  test.set, M=500, alpha=0.05, seed=1) {
  
  ## Zhang, Cheng, 2017, JASA 
  
  n0 <- dim(X)[1]
  p <- dim(X)[2]
  
  score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")
  node <- score.nodewiselasso(X, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
                              parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
  Theta <- node$out
  Gram<- t(X)%*%X/n0
  
  sreg <- scalreg::scalreg(X,Y)
  beta.hat <- sreg$coefficients
  sigma.sq <- sum((Y-X%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))
  
  index <- test.set
  
  
  
  Omega <- diag(Theta%*%Gram%*%t(Theta))*sigma.sq
  beta.db <- beta.hat+Theta%*%t(X)%*%(Y-X%*%beta.hat)/n0
  margin.st <- sqrt(n0)*abs(beta.db[index])/(sqrt(Omega[index]) + 1e-7)
  margin.nst <- sqrt(n0)*abs(beta.db[index])
  stat.st <- max(margin.st)
  stat.nst <- max(margin.nst)
  
  set.seed(seed)
  E <- matrix(rnorm(n0*M),n0,M)
  cvMat <- Theta[index,]%*%t(X)%*%E /sqrt(n0)
  stat.boot.nst <- apply(cvMat, 2, function(x) max(abs(x)))
  stat.boot.st <- apply(cvMat, 2, function(x) max(abs(x/ sqrt(Omega[index]))))
  
  # calculate the p-values
  pv.nst <- 1-ecdf(stat.boot.nst)(stat.nst)
  pv.st <-  1- ecdf(stat.boot.st)(stat.st)
  cri_val.nst <- quantile(stat.boot.nst,1-alpha)
  cri_val.st <- quantile(stat.boot.st,1-alpha)
  if (stat.nst> cri_val.nst) rej.nst <- 1 else rej.nst <-0
  if (stat.st> cri_val.st) rej.st <- 1 else rej.st <- 0
  stvec <- c(cri_val.st,stat.st,rej.st, pv.st)
  names(stvec) <- c("critical value","test statistic","reject", "p-value")
  nstvec <- c(cri_val.nst, stat.nst,rej.nst ,pv.nst)
  names(nstvec) <- c("critical value","test statistic","reject", "p-value")
  result <- list(stvec, nstvec)
  attr(result, 'test') <- 'max-type'
  return(result)
}


lasso_IC <- function(x, y, criteria="BIC"){
  require(glmnet)
  fit1 = cv.glmnet(x,y)
  # calculate REM_d for each lambda
  hbetaMat <- fit1$glmnet.fit$beta
  dfVec <- fit1$glmnet.fit$df
  
  getIC <- function(j,cri='BIC',  hbetaMat, dfVec){
    ntmp <- length(y)
    REM_d_tmp <- sum((y - x %*% hbetaMat[,j])^2)/ ntmp
    if(tolower(cri) == 'bic'){
      ic <- REM_d_tmp*(1 + log(ntmp)*dfVec[j]/(ntmp-dfVec[j]))
    }else if(tolower(cri) == 'aic'){
      ic <- REM_d_tmp*(1 + 2*dfVec[j]/(ntmp-dfVec[j]))
    }else{
      stop("Unsuppoted criteria! ")
    }
    return(ic)
  }
  
  criVec <- sapply(1:length(dfVec), getIC, cri=criteria,  hbetaMat=hbetaMat, dfVec =dfVec)
  
  hbetaMat[,which.min(criVec)]
}


vars_comp2 <- function(X, Y, X.ind, Y.ind, true_index, alpha=0.05, seed=1, sub.frac=0.3, standardized=F){ 
  
  # X <- dat1$X; Y <-dat1$Y; X.ind<-dat$X; Y.ind <- dat$Y;  true_index<- dat1$index_nz
  # alpha=0.05; seed=1; sub.frac=0.5; standardized=F
  library(glmnet)
  p <- ncol(X)
  Xall <- rbind(X,X.ind); Yall <- c(Y, Y.ind)
  timeVec <- rep(NA, 5)
  names(timeVec) <- c("TOSI", "CV Lasso", "BIC", "AIC", "scaled Lasso")
  
  
  ## CV Lasso
  tic <- proc.time()
  cvlist <- cv.glmnet(Xall,Yall)
  toc <- proc.time()
  timeVec[2] <- toc[3] - tic[3]
  
  
  ## TOSI
  tic <- proc.time()
  id_lambdamin <- which.min(cvlist$cvm)
  glmlist <- glmnet(Xall, Yall,lambda=cvlist$lambda[id_lambdamin])
  index_nz <- which(glmlist$beta !=0)
  cv_lasso_index <- index_nz
  glmlist <- glmnet(X.ind, Y.ind,lambda=cvlist$lambda[id_lambdamin])
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
    
    glmlist <- glmnet(X.ind, Y.ind,lambda=cvlist$lambda[id_lambda])
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
  toc <- proc.time()
  timeVec[1] <- toc[3] - tic[3]
  tosi_index <- index_nz
  
  ## BIC
  tic <- proc.time()
  hbeta_bic <- lasso_IC(Xall,Yall, criteria = 'BIC')
  bic_index <- which(hbeta_bic != 0)
  toc <- proc.time()
  timeVec[3] <- toc[3] - tic[3]
  
  ## AIC
  tic <- proc.time()
  hbeta_aic <- lasso_IC(Xall,Yall, criteria = 'AIC')
  aic_index <- which(hbeta_aic != 0)
  toc <- proc.time()
  timeVec[4] <- toc[3] - tic[3]
  
  # scaled Lasso
  library(scalreg)
  tic <- proc.time()
  reslist1 <- scalreg(Xall,Yall)
  hbeta_scaleLasso <- reslist1$coefficients
  scaleLasso_index <- which(hbeta_scaleLasso!=0)
  toc <- proc.time()
  timeVec[5] <- toc[3] - tic[3]
  
  cat(timeVec, '\n')
  
  
  
  
  indexList <- list(tosi_index, cv_lasso_index, bic_index, aic_index, scaleLasso_index)
  
  measureMat <- matrix(0, length(timeVec), 3)
  colnames(measureMat) <- c('NF', 'IN', 'CS')
  row.names(measureMat) <-  names(timeVec)
  measureMat[,1] <- sapply(indexList, length)
  for(jj in 1:length(indexList)){
    tmp1 <- intersect(indexList[[jj]], true_index)
    if(length(tmp1) == length(true_index)){
      measureMat[jj,2] <- 1
      if(length(indexList[[jj]]) == length(true_index))
        measureMat[jj,3] <- 1
    } 
  }
  measureMat <- cbind(measureMat, time=timeVec)
  
  return(measureMat)
}

TOSI_selectLambda <- function(X, Y, X.ind, Y.ind, true_index, alpha=0.05, seed=1, sub.frac=0.3, standardized=F){ 
  
  # X <- dat1$X; Y <-dat1$Y; X.ind<-dat$X; Y.ind <- dat$Y;  true_index<- dat1$index_nz
  # alpha=0.05; seed=1; sub.frac=0.5; standardized=F
  library(glmnet)
  p <- ncol(X)
  Xall <- rbind(X,X.ind); Yall <- c(Y, Y.ind)
  
  
  ## CV Lasso
  tic <- proc.time()
  cvlist <- cv.glmnet(Xall,Yall)
  toc <- proc.time()
  
  
  
  ## TOSI
  tic <- proc.time()
  id_lambdamin <- which.min(cvlist$cvm)
  glmlist <- glmnet(Xall, Yall,lambda=cvlist$lambda[id_lambdamin])
  index_nz <- which(glmlist$beta !=0)
  cv_lasso_index <- index_nz
  glmlist <- glmnet(X.ind, Y.ind,lambda=cvlist$lambda[id_lambdamin])
  index_nz <- which(glmlist$beta !=0)
  # stmin1 <- RegMin(X, Y,index_nz, alpha, seed, sub.frac, standardized);
  # Tmin <- stmin1[1,4] > alpha
  Tmin <- T
  index_zero <- setdiff(1:p, index_nz)
  # stmax1 <- RegMax(X, Y,  index_zero, alpha, seed, sub.frac, standardized); 
  # Tmax <-  stmax1[1,4] > alpha
  Tmax <- T
  # Seek the optimal lambda by Dichotomy optimization
  leftside <- 0; rightside <- id_lambdamin
  index_nz_old <- index_nz
  k <- 1
  while(Tmin ==T || Tmax==F){
    if(k > 3) break
    id_lambda <- round((leftside + rightside)/2)
    
    glmlist <- glmnet(X.ind, Y.ind,lambda=cvlist$lambda[id_lambda])
    index_nz <- which(glmlist$beta !=0)
    if(length(index_nz_old)==length(index_nz)) break;
    index_nz_old <- index_nz
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
    message('Searching penalty parameter for', k, '-th round')
    k <- k +1
  }
  message("Fishi bisection searching!")
  toc <- proc.time()
  timeVec<- toc[3] - tic[3]
  tosi_index <- index_nz
  
  
  indexList <- list(tosi_index)
  
  measureMat <- matrix(0,1, 3)
  colnames(measureMat) <- c('NF', 'IN', 'CS')
  row.names(measureMat) <-  "TOSI"
  measureMat[,1] <- sapply(indexList, length)
  for(jj in 1:length(indexList)){
    tmp1 <- intersect(indexList[[jj]], true_index)
    if(length(tmp1) == length(true_index)){
      measureMat[jj,2] <- 1
      if(length(indexList[[jj]]) == length(true_index))
        measureMat[jj,3] <- 1
    } 
  }
  measureMat <- cbind(measureMat, time=timeVec)
  
  return(measureMat)
  
}

# TOSI_selectLambda2 <- function(X, Y, X.ind, Y.ind, true_index, alpha=0.05, seed=1, sub.frac=0.3, standardized=F){ 
#   
#   # X <- dat1$X; Y <-dat1$Y; X.ind<-dat$X; Y.ind <- dat$Y;  true_index<- dat1$index_nz
#   # alpha=0.05; seed=1; sub.frac=0.5; standardized=F
#   library(glmnet)
#   p <- ncol(X)
#   Xall <- rbind(X,X.ind); Yall <- c(Y, Y.ind)
#   
#   require(glmnet)
#   cvlist = cv.glmnet(Xall, Yall)
#   # calculate REM_d for each lambda
#   hbetaMat <- cvlist$glmnet.fit$beta
#   dfVec <- cvlist$glmnet.fit$df
#   
#   getIC <- function(j,  hbetaMat, dfVec){
#     ntmp <- length(Yall)
#     REM_d_tmp <- sum((Yall - Xall %*% hbetaMat[,j])^2)/ ntmp
#    
#     ic <- REM_d_tmp*(1 + log(ntmp)*dfVec[j]/(ntmp-dfVec[j]))
#     return(ic)
#   }
#   
#   criVec <- sapply(1:length(dfVec), getIC,  hbetaMat=hbetaMat, dfVec =dfVec)
#   
#   
#   ## TOSI
#   tic <- proc.time()
#   id_lambdamin <- which.min(criVec)
#   ## set the searching interval
#   if(id_lambdamin >4 && id_lambdamin <= length(criVec)-4){
#     
#     idxSet <- seq(id_lambdamin-4,  id_lambdamin+4)
#     id_lambdamin <- 5
#   }else if(id_lambdamin <=4){
#     idxSet <- seq(1:(id_lambdamin+4))
#   }else if(id_lambdamin > length(criVec)-4){
#     idxSet <- seq((id_lambdamin-4): length(criVec))
#     id_lambdamin <- 5
#   }
#   lambdaset <- cvlist$lambda[idxSet]
#   
#   glmlist <- glmnet(Xall, Yall,lambda=lambdaset[id_lambdamin])
#   index_nz <- which(glmlist$beta !=0)
#   cv_lasso_index <- index_nz
#   glmlist <- glmnet(X.ind, Y.ind,lambda=lambdaset[id_lambdamin])
#   index_nz <- which(glmlist$beta !=0)
#   stmin1 <- RegMin(X, Y,index_nz, alpha, seed, sub.frac, standardized);
#   Tmin <- stmin1[1,4] > alpha
#   index_zero <- setdiff(1:p, index_nz)
#   # stmax1 <- RegMax(X, Y,  index_zero, alpha, seed, sub.frac, standardized); 
#   # Tmax <-  stmax1[1,4] > alpha
#   Tmax <- T
#   # Seek the optimal lambda by bisection searching
#   leftside <- 1; rightside <- id_lambdamin
#   k <- 1
#   index_nz_old <- index_nz
#   while(Tmin ==T || Tmax==F){
#     if(k > 3) break
#     id_lambda <- round((leftside + rightside)/2)
#     
#     glmlist <- glmnet(X.ind, Y.ind,lambda=lambdaset[id_lambda])
#     # if(sum(glmlist$beta!=0) == 0){
#     #   rightside <- id_lambda
#     #   k <- k + 1
#     #   next
#     # }
#     index_nz <- which(glmlist$beta !=0)
#     if(length(index_nz_old)==length(index_nz)) break;
#     index_nz_old <- index_nz
#     stmin1 <- RegMin(X, Y,index_nz, alpha, seed, sub.frac, standardized);
#     Tmin <- stmin1[1,4] > alpha
#     if(Tmin){
#       rightside <- id_lambda
#     }
#     if(!Tmin){
#       leftside <- id_lambda
#       index_zero <- setdiff(1:p, index_nz)
#       stmax1 <- RegMax(X, Y,  index_zero, alpha, seed, sub.frac, standardized); 
#       Tmax <-  stmax1[1,4] > alpha
#     }
#     cat('Searching penalty parameter for', k, '-th round\n')
#     k <- k +1
#   }
#   toc <- proc.time()
#   timeVec <- toc[3] - tic[3]
#   tosi_index <- index_nz
#   
#   
#   
#   indexList <- list(tosi_index)
#   
#   measureMat <- matrix(0,1, 3)
#   colnames(measureMat) <- c('NF', 'IN', 'CS')
#   row.names(measureMat) <-  "TOSI"
#   measureMat[,1] <- sapply(indexList, length)
#   for(jj in 1:length(indexList)){
#     tmp1 <- intersect(indexList[[jj]], true_index)
#     if(length(tmp1) == length(true_index)){
#       measureMat[jj,2] <- 1
#       if(length(indexList[[jj]]) == length(true_index))
#         measureMat[jj,3] <- 1
#     } 
#   }
#   measureMat <- cbind(measureMat, time=timeVec)
#   
#   return(measureMat)
# }


vars_comp <- function(X, Y, X.ind, Y.ind, true_index, alpha=0.05, seed=1, sub.frac=0.3, standardized=F){ 
  library(glmnet)
  p <- ncol(X)
  Xall <- rbind(X,X.ind); Yall <- c(Y, Y.ind)
  cvlist <- cv.glmnet(Xall,Yall)
  
  id_lambdamin <- which.min(cvlist$cvm)
  glmlist <- glmnet(Xall, Yall,lambda=cvlist$lambda[id_lambdamin])
  index_nz <- which(glmlist$beta !=0)
  cv_lasso_index <- index_nz
  
  glmlist <- glmnet(X.ind, Y.ind,lambda=cvlist$lambda[id_lambdamin])
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
    
    glmlist <- glmnet(X.ind, Y.ind,lambda=cvlist$lambda[id_lambda])
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
  tosi_index <- index_nz
  measureMat <- matrix(0, 2, 3)
  colnames(measureMat) <- c('NF', 'IN', 'CS')
  row.names(measureMat) <- c('cv_lasso', 'TOSI')
  measureMat[,1] <- c(length(cv_lasso_index), length(tosi_index))
  tmp1 <- intersect(cv_lasso_index, true_index)
  if(length(tmp1) == length(true_index)){
    measureMat[1,2] <- 1
    if(length(cv_lasso_index) == length(true_index))
      measureMat[1,3] <- 1
  } 
  tmp2 <- intersect(tosi_index, true_index)
  if(length(tmp2) == length(true_index)){
    measureMat[2,2] <- 1
    if(length(tosi_index) == length(true_index))
      measureMat[2,3] <- 1
  } 
  
  return(measureMat)
}

# test CV scad
vars_scad <- function(X, Y, true_index){ 
  library(ncvreg)
  p <- ncol(X)
  cvlist <- cv.ncvreg(X,Y)
  
  glmlist <- ncvreg(X, Y,lambda=cvlist$lambda.min)
  index_nz <- which(glmlist$beta !=0)
  cv_scad_index <- index_nz
  
  
  measureMat <- numeric(3)
  names(measureMat) <- c('NF', 'IN', 'CS')
  measureMat[1] <- length(cv_scad_index)
  tmp1 <- intersect(cv_scad_index, true_index)
  if(length(tmp1) == length(true_index)){
    measureMat[2] <- 1
    if(length(cv_scad_index) == length(true_index))
      measureMat[3] <- 1
  } 
  
  return(measureMat)
}


#### This function is from SILM R package.
#### Zhang, X., and Cheng, G. (2017) Simultaneous Inference for High-dimensional Linear Models, Journal of the American Statistical Association, 112, 757-768.
ST <- function (X.f, Y.f, sub.size, test.set, M = 500, alpha = 0.05) {
  require(glmnet)
  require(hdi)
  require(SIS)
  require(scalreg)
  n <- dim(X.f)[1]
  p <- dim(X.f)[2]
  n1 <- sub.size
  n0 <- n - floor(n1)
  S1 <- sample(1:n, floor(n1), replace = FALSE)
  X.sub <- X.f[S1, ]
  Y.sub <- Y.f[S1]
  cvfit <- cv.glmnet(X.sub, Y.sub, intercept = FALSE)
  cf <- as.numeric(coef(cvfit, s = "lambda.min"))[-1]
  set1 <- (1:p)[abs(cf) > 0]
  resi <- Y.sub - X.sub %*% cf
  beta.m <- t(standardize(X.sub[, -set1])) %*% resi
  screen.set <- sort(order(abs(beta.m), decreasing = TRUE)[1:(n0 - 
                                                                1 - length(set1))])
  a <- (1:p)[-set1]
  screen.set <- union(a[screen.set], set1)
  X <- X.f[-S1, screen.set]
  Y <- Y.f[-S1]
  score.nodewiselasso = getFromNamespace("score.nodewiselasso", 
                                         "hdi")
  node <- score.nodewiselasso(X, wantTheta = TRUE, verbose = FALSE, 
                              lambdaseq = "quantile", parallel = FALSE, ncores = 2, 
                              oldschool = FALSE, lambdatuningfactor = 1)
  Theta <- node$out
  Gram <- t(X) %*% X/n0
  sreg <- scalreg(X, Y)
  beta.hat <- sreg$coefficients
  sigma.sq <- sum((Y - X %*% beta.hat)^2)/(n0 - sum(abs(beta.hat) > 
                                                      0))
  test.set.i <- intersect(screen.set, test.set)
  index <- screen.set %in% test.set.i
  Omega <- diag(Theta %*% Gram %*% t(Theta)) * sigma.sq
  beta.db <- beta.hat + Theta %*% t(X) %*% (Y - X %*% beta.hat)/n0
  margin.st <- sqrt(n0) * abs(beta.db[index])/sqrt(Omega[index])
  margin.nst <- sqrt(n0) * abs(beta.db[index])
  stat.st <- max(margin.st)
  stat.nst <- max(margin.nst)
  stat.boot.st <- stat.boot.nst <- rep(NA, M)
  for (i in 1:M) {
    e <- rnorm(n0)
    xi.boot <- Theta[index, ] %*% t(X) %*% e * sqrt(sigma.sq)/sqrt(n0)
    stat.boot.nst[i] <- max(abs(xi.boot))
    stat.boot.st[i] <- max(abs(xi.boot/sqrt(Omega[index])))
  }
  if (stat.nst > quantile(stat.boot.nst, 1 - alpha)) 
    rej.nst <- 1
  else rej.nst <- 0
  if (stat.st > quantile(stat.boot.st, 1 - alpha)) 
    rej.st <- 1
  else rej.st <- 0
  result <- c(rej.st, rej.nst)
  names(result) <- c("rej.st", "rej.nst")
  return(result)
}

Reg_BC <- function(X,Y,  G1, alpha= 0.05, method="bonferroni"){
  require(hdi)
  lassop <- lasso.proj(X, Y)
  
  Pvec <- lassop$pval[G1]
  adjPval <- p.adjust(Pvec, method = method)
  reject <- 0
  if(any(adjPval < alpha)) reject <- 1
  return(reject)
}
getReg_BC <- function(X,Y, G1list, ...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- Reg_BC(X,Y, G1 =G1list[[l]], ...)
    return(lreg)
  })
  return(tmp)
}

RegMin_BC <- function(X,Y,  G1, alpha= 0.05, method="bonferroni"){
  require(hdi)
  lassop <- lasso.proj(X, Y)
  
  Pvec <- lassop$pval[G1]
  adjPval <- p.adjust(Pvec, method = method)
  reject <- 0
  if(all(adjPval < alpha)) reject <- 1
  return(reject)
}
getRegMin_BC <- function(X,Y, G1list, ...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- RegMin_BC(X,Y, G1 =G1list[[l]], ...)
    return(lreg)
  })
  return(tmp)
}

getRegSTMax <- function(X,Y, G1list,...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- regST_max(X, Y, test.set =G1list[[l]], ...)
    return(lreg[[1]])
  })
  colnames(tmp) <- paste0('G1', 1:nS1)
  row.names(tmp) <- c('CriticalValue', 'TestStatistic',
                      'reject_status',    ' p-value')
  return(tmp)
}

getRegSTMax2 <- function(X, Y, G1list,...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- ST(X, Y,sub.size=nrow(X)*0.3, test.set  =G1list[[l]], ...)
    return(lreg)
  })
  colnames(tmp) <- paste0('G1', 1:nS1)
  row.names(tmp) <- c("rej.st", "rej.nst")
  return(tmp)
}

getRegMax <- function(X, Y, G1list,...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- RegMax(X, Y, G1 =G1list[[l]], ...)
    return(lreg)
  })
  colnames(tmp) <- paste0('G1', 1:nS1)
  row.names(tmp) <- c('CriticalValue', 'TestStatistic',
                      'reject_status',    ' p-value')
  return(tmp)
}
getmultiRegMax <- function(X,Y, G1list,...){
  nS1 <- length(G1list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- multiRegMax(X, Y, G1 =G1list[[l]], ...)
    return(lreg)
  })
  if(is.vector(tmp)) names(tmp) <- paste0('G1', 1:nS1)
  if(is.matrix(tmp)) colnames(tmp) <- paste0('G1', 1:nS1)
  return(tmp)
}

getRegMin <- function(X,Y, G2list,...){
  nS1 <- length(G2list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- RegMin(X, Y,G2 =G2list[[l]], ...)
    return(lreg)
  })
  colnames(tmp) <- paste0('G2', 1:nS1)
  row.names(tmp) <- c('CriticalValue', 'TestStatistic',
                      'reject_status',    ' p-value')
  return(tmp)
}
getmultiRegMin <- function(X, Y, G2list,...){
  nS1 <- length(G2list)
  tmp <- sapply(1:nS1, function(l){
    lreg <- multiRegMin(X, Y, G2 =G2list[[l]], ...)
    return(lreg)
  })
  if(is.vector(tmp)) names(tmp) <- paste0('G2', 1:nS1)
  if(is.matrix(tmp)) colnames(tmp) <- paste0('G2', 1:nS1)
  return(tmp)
}
