## Date: 2022-06-20

# Functions ---------------------------------------------------------------


gendata_logit <- function(n=100, p = 20, s0=3, rho=1, seed=1){
  set.seed(1)
  beta <- rep(0,p)
  beta[1:s0] <- runif(s0,0,2) *rho
  set.seed(seed)
  Sigma <- matrix(NA, p, p)
  for (i in 1:p) Sigma[i,] <- 0.2^(abs(i-(1:p)))
  X <- matrix(rnorm(n*p), n, p)
  X <- t(t(chol(Sigma))%*%t(X))
  
  eta <- X %*%beta
  Pi <- 1/(1+ exp(-eta))
  # Only first three coefficients are not zeros
  y <- rbinom(n, 1, Pi)
  while(sum(y)/n<0.02 | sum(y)/n>0.98 ){
    y=rbinom(n,1,Pi)
  }
  
  return(list(Y=y, X=X, beta0=beta, index_nz=1:s0))
}
RegMin_binary <- function(X, Y,  G1, alpha=0.05, seed=1, sub.frac=0.5) {
  # alpha=0.05; seed=1; sub.frac=0.5; standardized=F
  require(glmnet)
  # require(hdi)
  # require(SIS)
  # require(scalreg)
  n <- dim(X)[1]
  p <- dim(X)[2]
  # Devide sample into two parts
  n1 <- floor(sub.frac*n) 
  n0 <- n-floor(n1)
  set.seed(seed)
  S1 <- sample(1:n, n1, replace=FALSE)
  X.sub <- X[S1,]
  
  Y.sub <- Y[S1]
  cf <- as.numeric(glmnet(X.sub, Y.sub, lambda=1e-5,  intercept=TRUE, family = binomial)$beta)
  cf_testset <- cf[G1]
  # sort the coefficients
  K <- min(1, length(G1) ) # only caputure the maximum coefs.
  id_test.set <- order(abs(cf_testset), decreasing = FALSE)[1:K]
  test.set <- G1[id_test.set]
  tsXid <- setdiff(1:n, S1)
  Xts <- X[tsXid, ]; n0 <- nrow(Xts)
  Yts <- Y[tsXid]
  
  
  
  
  index <- test.set
  
  # Est <- get_debias_infer(Xts, Yts, index)
  Est = GLM_binary(X =Xts, y = Yts, index = index)
  est_bias_correct <- Est$prop.est
  se <- Est$se
  T1 <- abs(est_bias_correct/se)
  PV <- 2*(1-pnorm(T1))
  maxC1 <- 1.96
  
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1,PV< alpha, PV)
  row.names(pMat) <- c('Ztest')
  
  
  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'ToMax'
  return(pMat)
}
MultiRegMin_binary <- function(X, Y,  G2, alpha=0.05, seed=1, sub.frac=0.5, Nsplit = 5, r= 0.5, BH=T){
  require(glmnet)
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
    Y.sub <- Y[S1]
    cf <- as.numeric(glmnet(X.sub, Y.sub, lambda=1e-5,  intercept=TRUE, family = binomial)$beta)
    cf_testset <- cf[G2]
    # sort the coefficients
    K <- min(1, length(G2) ) # only caputure the minimum coef.
    id_test.set <- order(abs(cf_testset))[1:K]
    test.set <- G2[id_test.set]
    tsXid <- setdiff(1:n, S1)
    Xts <- X[tsXid, ]; n0 <- nrow(Xts)
    Yts <- Y[tsXid]
    index <- test.set
    
    # Est <- get_debias_infer(Xts, Yts, index)
    Est = GLM_binary(X =Xts, y = Yts, index = index)
    est_bias_correct <- Est$prop.est
    se <- Est$se
    T1 <- abs(est_bias_correct/se)
    PV <- 2*(1-pnorm(T1))
    
    Pvec[i] <-  PV
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


get_debias_infer <- function(X,y, index, intercept=TRUE){
  data <- na.omit(data.frame(y, X))
  X <- as.matrix(data[, -1])
  y <- as.vector(data[, 1])
  p <- ncol(X)
  n <- nrow(X)
  mean = colMeans(X)
  M = matrix(rep(mean, nrow(X)), byrow = T, nrow = nrow(X), 
             ncol = ncol(X))
  X = X - M
  
  ej = rep(0, p)
  ej[index] = 1
  Est <- SIHR::LF_logistic(X = X, y = y, loading = ej, 
                           weight = NULL, trans = FALSE, intercept = intercept, 
                           intercept.loading = FALSE)
  
  
  return(Est)
}
RegMax_binary <- function(X, Y,  G1, alpha=0.05, seed=1, sub.frac=0.5) {
  # alpha=0.05; seed=1; sub.frac=0.5; standardized=F
  require(glmnet)
  require(SIHR)
  # require(SIS)
  # require(scalreg)
  n <- dim(X)[1]
  p <- dim(X)[2]
  # Devide sample into two parts
  n1 <- floor(sub.frac*n) 
  n0 <- n-floor(n1)
  set.seed(seed)
  S1 <- sample(1:n, n1, replace=FALSE)
  X.sub <- X[S1,]
  
  Y.sub <- Y[S1]
  tic <- proc.time()
  cvfit <- cv.glmnet(X.sub, Y.sub, intercept=TRUE, lambda=c(1, 0.1, 1e-2,1e-5),family = binomial)
  toc <- proc.time()
  cf <- as.numeric(coef(cvfit, s="lambda.min"))[-1]
  # cf <- as.numeric(glmnet(X.sub, Y.sub, lambda=1e-4,  intercept=TRUE, family = binomial)$beta)
  cf_testset <- cf[G1]
  # sort the coefficients
  K <- min(1, length(G1) ) # only caputure the maximum coefs.
  id_test.set <- order(abs(cf_testset), decreasing = T)[1:K]
  test.set <- G1[id_test.set]
  tsXid <- setdiff(1:n, S1)
  Xts <- X[tsXid, ]; n0 <- nrow(Xts)
  Yts <- Y[tsXid]
  
  
  
  
  index <- test.set
  
  # Est <- get_debias_infer(Xts, Yts, index)
  Est = GLM_binary(X =Xts, y = Yts, index = index)
  est_bias_correct <- Est$prop.est
  se <- Est$se
  T1 <- abs(est_bias_correct/se)
  PV <- 2*(1-pnorm(T1))
  maxC1 <- 1.96
  
  pMat <- matrix(0,1,4)
  pMat[1,] <- c(maxC1, T1,PV< alpha, PV)
  row.names(pMat) <- c('Ztest')
  
  
  colnames(pMat) <- c('CriticalValue', 'TestStatistic', 'reject_status', 'p-value')
  class(pMat) <- 'ToMax'
  return(pMat)
}
MultiRegMax_binary <- function(X, Y,  G1, alpha=0.05, seed=1, sub.frac=0.5, Nsplit = 5, r= 0.5, BH=T){
  require(glmnet)
  require(SIHR)
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
    Y.sub <- Y[S1]
    tic <- proc.time()
    cvfit <- cv.glmnet(X.sub, Y.sub, intercept=TRUE, lambda=c(1, 0.1, 1e-2,1e-5),family = binomial)
    toc <- proc.time()
    cf <- as.numeric(coef(cvfit, s="lambda.min"))[-1]
    
    cf_testset <- cf[G1]
    # sort the coefficients
    K <- min(1, length(G1) ) # only caputure the maximum coefs.
    id_test.set <- order(abs(cf_testset), decreasing = T)[1:K]
    test.set <- G1[id_test.set]
    tsXid <- setdiff(1:n, S1)
    Xts <- X[tsXid, ]; n0 <- nrow(Xts)
    Yts <- Y[tsXid]
    
    index <- test.set
    
    # Est <- get_debias_infer(Xts, Yts, index)
    Est = GLM_binary(X =Xts, y = Yts, index = index)
    est_bias_correct <- Est$prop.est
    se <- Est$se
    T1 <- abs(est_bias_correct/se)
    PV <- 2*(1-pnorm(T1))
    
    Pvec[i] <-  PV 
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


Binary_infer <- function(X, y, Index, alpha = 0.05){
  
  nIndex <- length(Index)
  P.vals <- rep(NA, nIndex)
  for(i in 1:nIndex){
    Est = GLM_binary(X =X, y = y, index = Index[i])
    est_bias_correct <- Est$prop.est
    se <- Est$se
    T1 <- abs(est_bias_correct/se)
    PV <- 2*(1-pnorm(T1))
    P.vals[i] <- PV
  }
  
  return(P.vals)
  
}




# Simulation --------------------------------------------------------------
library(SIHR)

#--Check the functions---
dat1 <- gendata_logit(rho=1, n=100, p=30)
dat1$beta0
str(dat1)

Est = GLM_binary(X = dat1$X, y = dat1$Y, index = 1:3)
Est[1:4]

X <- dat1$X; Y <- dat1$Y
G1 <- 1:3

RegMax_binary(X, Y, c(1,3))
(res <- MultiRegMax_binary(X, Y, c(1,3),  Nsplit = 3, r=1/3, BH=T))

RegMin_binary(X, Y, 5:10)
MultiRegMin_binary(X, Y, 5:10, Nsplit = 3, r=1/3, BH=T)

## Method in Cai (2021)

Binary_infer(X, Y, 10:20)

## Method in  van de Geer, S., Bühlmann, P., Ritov, Y. and Dezeure, R. (2014)
fit.lasso <- hdi:::lasso.proj(x = X, y = Y, standardize = F, family = "binomial")
fit.lasso$pval.corr

## Method in Bulman
fit.ridge <- hdi:::ridge.proj(x = X, y = Y, standardize = F, family = "binomial")
fit.ridge$pval.corr

#--Finish  the check  of functions---


########### VBRD-14 and B-09

n <- 100
p <- 30; s0 <- 5
G1list <- list(c(p-1, p), (p/2):p, (s0+1):p) # True H0
G2list <-  list(c(2, s0+1), c(3, (s0+1):p),c(3,4, (s0+1):p) ) # False H0

rho <- 0.5
Glist <- c(G1list, G2list)
alpha <- 0.05
N <- 100


powersizeArray <- array(dim=c(2, 6, N))
for(i in 1:N){
  # i <- 1
  message("i = ", i)
  dat1 <- gendata_logit(rho=rho, n=n, p=p, s0=s0, seed = i)
  X <- dat1$X; Y <- dat1$Y
  ## Method in  van de Geer, S., Bühlmann, P., Ritov, Y. and Dezeure, R. (2014)
  fit.lasso <- hdi:::lasso.proj(x = X, y = Y, standardize = F, family = "binomial")
  pc <- fit.lasso$pval.corr
  reject_vec <- rep(NA, 6)
  
  for(ii in 1:6){
    reject_vec[ii] <- any(pc[Glist[[ii]]]<alpha)
  }
  
  powersizeArray[1,,i] <- reject_vec
  ## Method in Bulman
  fit.ridge <- hdi:::ridge.proj(x = X, y = Y, standardize = F, family = "binomial")
  
  pc <- fit.ridge$pval.corr
  reject_vec <- rep(NA, 6)
  for(ii in 1:6){
    reject_vec[ii] <- any(pc[Glist[[ii]]]<alpha)
  }
  
  powersizeArray[2,,i] <- reject_vec
  
 
  
}
apply(powersizeArray, c(1,2), mean, na.rm=T)
save(powersizeArray, file= 'Proj_rho05_n100_p30_N100.rds')


### RegMax/ MultiRegMax
N <- 500
powersizeArray <- array(dim=c(2, 6, N))
for(i in 1:N){
  # i <- 1
  tic <- proc.time()
  message("i = ", i)
  dat1 <- gendata_logit(rho=rho, n=n, p=p, s0=s0, seed = i)
  X <- dat1$X; Y <- dat1$Y
  
  reject_vec <- rep(NA, 6)
  for(ii in 1:6){
    # ii <- 5
    message("RegMax_binary: ii=", ii)
    re <- RegMax_binary(X, Y, Glist[[ii]])
    reject_vec[ii] <- re[3]
  }
  powersizeArray[1,,i] <- reject_vec
  
  toc <- proc.time()
  toc-tic
  
  Ns <- 2
  reject_vec <- rep(NA, 6)
  for(ii in 1:6){
    message("RegMax_binary: ii=", ii)
    re <- MultiRegMax_binary(X, Y, Glist[[ii]],   Nsplit = Ns, r=1/Ns, BH=T)
    reject_vec[ii] <- re[2]
  }
  powersizeArray[2,,i] <- reject_vec
  
  
  
}
apply(powersizeArray, c(1,2), mean, na.rm=T)

save(powersizeArray, file= 'ToMax1_2_rho05_n100_p30_N200.rds')

###MultiRegMax(5)
N <- 500
powersizeArray <- array(dim=c(1, 6, N))
for(i in 1:N){
  # i <- 1
  tic <- proc.time()
  message("i = ", i)
  dat1 <- gendata_logit(rho=rho, n=n, p=p, s0=s0, seed = i)
  X <- dat1$X; Y <- dat1$Y
  
  Ns <- 5
  reject_vec <- rep(NA, 6)
  for(ii in 1:6){
    message("RegMax_binary: ii=", ii)
    re <- MultiRegMax_binary(X, Y, Glist[[ii]],   Nsplit = Ns, r=1/Ns, BH=T)
    reject_vec[ii] <- re[2]
  }
  powersizeArray[1,,i] <- reject_vec
  
  tmp <- apply(powersizeArray, c(1,2), mean, na.rm=T)
  cat(tmp, '\n')
  
}
apply(powersizeArray, c(1,2), mean, na.rm=T)
save(powersizeArray, file= 'ToMax5_rho05_n100_p30_N200.rds')


##
###MultiRegMax(8)
N <- 500
powersizeArray <- array(dim=c(1, 6, N))
for(i in 1:N){
  # i <- 1
  tic <- proc.time()
  message("i = ", i)
  dat1 <- gendata_logit(rho=rho, n=n, p=p, s0=s0, seed = i)
  X <- dat1$X; Y <- dat1$Y
  
  Ns <- 8
  reject_vec <- rep(NA, 6)
  for(ii in 1:6){
    message("RegMax_binary: ii=", ii)
    re <- MultiRegMax_binary(X, Y, Glist[[ii]],   Nsplit = Ns, r=1/Ns, BH=T)
    reject_vec[ii] <- re[2]
  }
  powersizeArray[1,,i] <- reject_vec
  
  tmp <- apply(powersizeArray, c(1,2), mean, na.rm=T)
  cat(tmp, '\n')
  
}
apply(powersizeArray, c(1,2), mean, na.rm=T)

save(powersizeArray, file= 'ToMax8_rho05_n100_p30_N200.rds')


### Cai et al. 2021 (CZM-21)
powersizeArray <- array(dim=c(1, 6, N))
for(i in 1:N){
  # i<- 1
  message("i = ", i)
  tic <- proc.time()
  dat1 <- gendata_logit(rho=rho, n=n, p=p, s0=s0, seed = i)
  X <- dat1$X; Y <- dat1$Y
  
  pvec_tmp <- Binary_infer(X, Y, 1:p)
  pc <- p.adjust(pvec_tmp, method='BY')
  reject_vec <- rep(NA, 6)
  for(ii in 1:6){
    reject_vec[ii] <- any(pc[Glist[[ii]]]<alpha)
  }
  powersizeArray[1,,i] <- reject_vec
  
  toc <- proc.time()
  cat(toc[3]- tic[3], '\n')
  
  tmp <- apply(powersizeArray, c(1,2), mean, na.rm=T)
  cat(tmp, '\n')
}
apply(powersizeArray, c(1,2), mean, na.rm=T)

save(powersizeArray, file= 'Cai2021_rho05_n100_p30_N100.rds')

# ## ToMin test for logistic regression ----------------------------------------------------------


library(SIHR)
n <- 100 
p <- 30; s0 <- 5

G21list <- list(c(p-1, p), (p/2):p, (s0+1):p)
G22list <- list(c(1,2 ), 1:4, 1:s0)

G2list <- c(G21list, G22list)

rho <- 0.5

alpha <- 0.05
N <- 500
powersizeArray <- array(dim=c(2,6, N))
for(i in 1:N){
  # i <- 1
  message("i = ", i)
  dat1 <- gendata_logit(rho=rho, n=n, p=p, s0=s0, seed = i)
  X <- dat1$X; Y <- dat1$Y
  dat1$beta0[1:s0]
  reject_vec <- rep(NA, 6)
  for(ii in 1:6){
    message("RegMin_binary: ii=", ii)
    re <- RegMin_binary(X, Y, G2list[[ii]], sub.frac = 0.3)
    reject_vec[ii] <- re[3]
  }
  powersizeArray[1,,i] <- reject_vec
  
  
  Ns <- 2
  reject_vec <- rep(NA, 6)
  for(ii in 1:6){
    message("RegMin_binary: ii=", ii)
    re <- MultiRegMin_binary(X, Y, G2list[[ii]],seed=0, sub.frac = 0.3,  Nsplit = Ns, r=1/Ns, BH=T)
    reject_vec[ii] <- re[2]
  }
  powersizeArray[2,,i] <- reject_vec
  
  tmp <- apply(powersizeArray, c(1,2), mean, na.rm=T)
  cat(tmp, '\n')
}

apply(powersizeArray, c(1,2), mean, na.rm=T)

save(powersizeArray, file= 'ToMin1_2_rho05_n100_p30_N200.rds')

## MultiMin repeat 5:
alpha <- 0.05
N <- 500
powersizeArray <- array(dim=c(1,6, N))
for(i in 1:N){
  message("i = ", i)
  dat1 <- gendata_logit(rho=rho, n=n, p=p, s0=s0, seed = i)
  X <- dat1$X; Y <- dat1$Y
  dat1$beta0[1:s0]
  
  Ns <- 5
  reject_vec <- rep(NA, 6)
  for(ii in 1:6){
    message("RegMin_binary: ii=", ii)
    re <- MultiRegMin_binary(X, Y, G2list[[ii]],seed=0, sub.frac = 0.3,  Nsplit = Ns, r=1/Ns, BH=T)
    reject_vec[ii] <- re[2]
  }
  powersizeArray[1,,i] <- reject_vec
  
  tmp <- apply(powersizeArray, c(1,2), mean, na.rm=T)
  cat(tmp, '\n')
}

apply(powersizeArray, c(1,2), mean, na.rm=T)

save(powersizeArray, file= 'ToMin5_rho05_n100_p30_N200.rds')


## MultiMin repeat 8:
alpha <- 0.05
N <- 500
powersizeArray <- array(dim=c(1,6, N))
for(i in 1:N){
  
  message("i = ", i)
  dat1 <- gendata_logit(rho=rho, n=n, p=p, s0=s0, seed = i)
  X <- dat1$X; Y <- dat1$Y
  dat1$beta0[1:s0]
  
  Ns <- 8
  reject_vec <- rep(NA, 6)
  for(ii in 1:6){
    message("RegMin_binary: ii=", ii)
    re <- MultiRegMin_binary(X, Y, G2list[[ii]],seed=0, sub.frac = 0.3,  Nsplit = Ns, r=1/Ns, BH=T)
    reject_vec[ii] <- re[2]
  }
  powersizeArray[1,,i] <- reject_vec
  
  tmp <- apply(powersizeArray, c(1,2), mean, na.rm=T)
  cat(tmp, '\n')
}

apply(powersizeArray, c(1,2), mean, na.rm=T)
save(powersizeArray, file= 'ToMin8_rho05_n100_p30_N200.rds')



