
setwd("D:\\LearnFiles\\Research paper\\Sparse factor inference\\Rcode\\Rcode_submits\\")
source("definedFunctions.R")
# Data generating function ------------------------------------------------


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
#----------------




# Start comprison in simulation -------------------------------------------

# ToMax, ZC1-17, and ZZ-14 -------------------------------------------------------------------

####rho=0.3/0.5/ 0.8, True H0-----------------------------------------------------

parafun1 <- function(i, rho=0.5){
  ## i <- 1
  cat('i=', i, '\n')
  p <- 50
  s0 <-  5 
  
  G1list <- list(c(p-1, p), (p/2):p, (s0+1):p)
  
  
  Ns_set <- c(1, 2, 5, 8)# , 15, 20) to speed up
  n1 <- 50; n2 <- 100
  cat('i=', i, '\n')
  tmp1Mat <- matrix(NA, length(Ns_set)+2, 3)
  dat1 <- gendata_Reg(n=n1, p, s0, rho, seed=i)
  # 1. ToMax method: Choose different split times
  for(js in 1:length(Ns_set)){
    # js <- 1
    tmp1Mat[js, ] <- getmultiRegMax(dat1$X, dat1$Y,G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  # 2. ZC1-17
  tmp1Mat[js+1, ] <-  getRegSTMax(dat1$X, dat1$Y,G1list)[3, ]
  # 3. ZZ-14
  tmp1Mat[js+2, ] <- getReg_BC(dat1$X, dat1$Y, G1list, method = 'BY')
  
  tmp2Mat <- matrix(NA, length(Ns_set)+2, 3)
  
  ## Run n=n2
  dat1 <- gendata_Reg(n=n2, p, s0, rho, seed=i)
  # Choose different split times
  for(js in 1:length(Ns_set)){
    tmp2Mat[js, ] <- getmultiRegMax(dat1$X, dat1$Y,G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  tmp2Mat[js+1, ] <- getRegSTMax(dat1$X, dat1$Y,G1list)[3, ]
  tmp2Mat[js+2, ] <- getReg_BC(dat1$X, dat1$Y, G1list, method = 'BY')
  cbind(tmp1Mat, tmp2Mat) # 6*6 matrix
}

parafun1(1)
# Parallel
library(parallel)
ncore <- detectCores()
cl <- makeCluster(20) # ncore/2
clusterExport(cl, list("gendata_Reg", "getRegMax","RegMax", "multiRegMax",
                       "regST_max","getReg_BC", "getTheta" ,"cv.nodewise.bestlambda", 
                       "score.partialnodewiselasso", "nodewise.getlambdasequence", 
                       "improve.lambda.pick", "score.getpartialThetaforlambda", 
                       "getRegSTMax",  "getmultiRegMax", "Reg_BC" ))
# clusterCall(cl, function() library(hdi))
N <- 500
tic <- proc.time()
iArray <- parSapply(cl, 1:N, parafun1, rho=0.8,  # 0.5
                    simplify = 'array')
toc <- proc.time() - tic; toc
time_used <- toc[3]

# save(iArray, time_used, file='rho05_regMax_trueNULL.rds')
save(iArray, time_used, file='rho08_regMax_trueNULL.rds')



####rho=0.3/0.5/0.8, False H0-----------------------------------------------------
parafun1 <- function(i, rho=0.5){
  cat('i=', i, '\n')
  p <- 50
  s0 <-  5 
  
  G1list <-  list(c(2, s0+1), c(3, (s0+1):p),c(3,4, (s0+1):p) )
  
  
  Ns_set <- c(1, 2, 5, 8)# , 15, 20) to speed up
  n1 <- 50; n2 <- 100
  cat('i=', i, '\n')
  tmp1Mat <- matrix(NA, length(Ns_set)+2, 3)
  dat1 <- gendata_Reg(n=n1, p, s0, rho, seed=i)
  # Choose different split times
  for(js in 1:length(Ns_set)){
    # js <- 1
    tmp1Mat[js, ] <- getmultiRegMax(dat1$X, dat1$Y,G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  tmp1Mat[js+1, ] <- getRegSTMax(dat1$X, dat1$Y,G1list)[3, ]
  tmp1Mat[js+2, ] <- getReg_BC(dat1$X, dat1$Y, G1list, method = 'BY')
  
  tmp2Mat <- matrix(NA, length(Ns_set)+2, 3)
  dat1 <- gendata_Reg(n=n2, p, s0, rho, seed=i)
  # Choose different split times
  for(js in 1:length(Ns_set)){
    tmp2Mat[js, ] <- getmultiRegMax(dat1$X, dat1$Y,G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  tmp2Mat[js+1, ] <- getRegSTMax(dat1$X, dat1$Y,G1list)[3, ]
  tmp2Mat[js+2, ] <- getReg_BC(dat1$X, dat1$Y, G1list, method = 'BY')
  cbind(tmp1Mat, tmp2Mat) # 6*6 matrix
}

# Parallel
library(parallel)
ncore <- detectCores()
cl <- makeCluster(15) # ncore/2
clusterExport(cl, list("gendata_Reg", "getRegMax","RegMax", "multiRegMax",
                       "regST_max","getReg_BC", "getTheta" ,"cv.nodewise.bestlambda", 
                       "score.partialnodewiselasso", "nodewise.getlambdasequence", 
                       "improve.lambda.pick", "score.getpartialThetaforlambda", 
                       "getRegSTMax",  "getmultiRegMax", "Reg_BC" ))

N <- 500
tic <- proc.time()
iArray <- parSapply(cl, 1:N, parafun1, rho=0.8, # 0.5
                    simplify = 'array')
toc <- proc.time() - tic; toc
time_used <- toc[3]
save(iArray, time_used, file='rho08_regMax_falseNULL.rds')



#### rho=0.3/0.5/0.8, ToMin(L), True  H0-----------------------------
parafun1 <- function(i, rho=0.5){
  
  p <- 50
  s0 <-  5 
  
  G1list <-  list(c(p-1, p), (p/2):p, (s0+1):p)
  
  
  Ns_set <- c(1, 2, 5, 8)# , 15, 20)
  n1 <- 50; n2 <- 100
  cat('i=', i, '\n')
  tmp1Mat <- matrix(NA, length(Ns_set), 3)
  dat1 <- gendata_Reg(n=n1, p, s0, rho, seed=i)
  # Choose different split times
  for(js in 1:length(Ns_set)){
    #js <- 1
    tmp1Mat[js, ] <- getmultiRegMin(dat1$X, dat1$Y,G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  
  
  tmp2Mat <- matrix(NA, length(Ns_set), 3)
  dat1 <- gendata_Reg(n=n2, p, s0, rho, seed=i)
  # Choose different split times
  for(js in 1:length(Ns_set)){
    tmp2Mat[js, ] <-  getmultiRegMin(dat1$X, dat1$Y,G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  
  cbind(tmp1Mat, tmp2Mat) # 6*6 matrix
}

# Parallel
library(parallel)
ncore <- detectCores()
cl <- makeCluster(20) # ncore/2
clusterExport(cl, list("gendata_Reg", "getRegMax","RegMax", "multiRegMax",
                       "regST_max","getReg_BC", "getTheta" ,"cv.nodewise.bestlambda", 
                       "score.partialnodewiselasso", "nodewise.getlambdasequence", 
                       "improve.lambda.pick", "score.getpartialThetaforlambda", 
                       "getRegSTMax",  "getmultiRegMax", "Reg_BC", "getmultiRegMin",
                       "multiRegMin"))

N <- 500
tic <- proc.time()
iArray <- parSapply(cl, 1:N, parafun1, rho=0.8,  # 
                    simplify = 'array')
toc <- proc.time() - tic; toc
time_used <- toc[3]
save(iArray, time_used, file='rho08_regMin_trueNULL.rds')
# save(iArray, time_used, file='rho05_regMin_trueNULL.rds')

#### rho=0.5, ToMin(L), False H0-----------------------------

parafun1 <- function(i, rho=0.5){
  
  p <- 50
  s0 <-  5 
  G1list <-  list(c(1,2 ), 1:4, 1:s0)
  
  
  Ns_set <- c(1, 2, 5, 8) # c(2, 5, 8, 15, 20)
  n1 <- 50; n2 <- 100
  cat('i=', i, '\n')
  tmp1Mat <- matrix(NA, length(Ns_set), 3)
  dat1 <- gendata_Reg(n=n1, p, s0, rho, seed=i)
  # Choose different split times
  for(js in 1:length(Ns_set)){
    # js <- 1
    tmp1Mat[js, ] <- getmultiRegMin(dat1$X, dat1$Y,G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  
  
  tmp2Mat <- matrix(NA, length(Ns_set), 3)
  dat1 <- gendata_Reg(n=n2, p, s0, rho, seed=i)
  # Choose different split times
  for(js in 1:length(Ns_set)){
    tmp2Mat[js, ] <-  getmultiRegMin(dat1$X, dat1$Y,G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  cbind(tmp1Mat, tmp2Mat) # 6*6 matrix
}

# Parallel
library(parallel)
ncore <- detectCores()
cl <- makeCluster(20) # ncore/2
clusterExport(cl, list("gendata_Reg", "getRegMax","RegMax", "multiRegMax",
                       "regST_max","getReg_BC", "getTheta" ,"cv.nodewise.bestlambda", 
                       "score.partialnodewiselasso", "nodewise.getlambdasequence", 
                       "improve.lambda.pick", "score.getpartialThetaforlambda", 
                       "getRegSTMax",  "getmultiRegMax", "Reg_BC", "getmultiRegMin",
                       "multiRegMin"))

N <- 500
tic <- proc.time()
iArray <- parSapply(cl, 1:N, parafun1, rho=0.8,  # 0.5, 
                    simplify = 'array')
toc <- proc.time() - tic; toc
time_used <- toc[3]

save(iArray, time_used, file='rho08_regMin_falseNULL.rds')






## rho control the signal strength
rho_set <- c(0.3, 0.5, 0.8)
rho <- 0.8 
library(hdi)
p <- 50
s0 <-  5 
G1list <- list(c(p-1, p), (p/2):p, (s0+1):p)
G2list <-  list(c(2, s0+1), c(3, (s0+1):p),c(3,4, (s0+1):p) )

alpha <- 0.05
N <- 500




# ZC3-17(1/3) and ZC3-17(1/5):  three-step procecedure in Zhang and Cheng (2017) -------------------------------
library(SILM)
rho <-  0.8 #0.3 #0.5 #  
n1 <-    50 #100 #
library(hdi)
p <- 50
s0 <-  5 
G1list <- list(c(p-1, p), (p/2):p, (s0+1):p)
G2list <-  list(c(2, s0+1), c(3, (s0+1):p),c(3,4, (s0+1):p) )

transferST <- function(res1){
  ## tranfer the rejection to 1, non-rejection to 0.
  res <- c(0, 0)
  names(res) <- names(res1)[c(2,4)]
  if(res1[[2]] == "fail to reject"){
    res[1] <- 0
  }else{
    res[1] <- 1
  }
  
  if(res1[[4]] == "fail to reject"){
    res[2] <- 0
  }else{
    res[2] <- 1
  }
  return(res)
}

alpha <- 0.05
N <- 500
subfracVec <- c(1/5, 1/3)

tic_total <- proc.time()
powersizeArray <- array(dim=c(6,4, N)) # 6 rows denotes the true hypotheises with first three, and false hypotheises with last three. 
for(i in 1:N){
  # i <- 1
  message("i = ", i)
  
  set.seed(i)
  dat1 <- gendata_Reg(n=n1, p, s0, rho, seed=i)
  
  tic <- proc.time()
  jj <- 1 # sub.frac = 1/5
  for(r in 1:length(G1list)){
    # r <- 2
    res1 <- c(NA, NA)
    try({res1 <- SILM::ST(dat1$X, dat1$Y, sub.size=n1*subfracVec[jj], G1list[[r]])
    res1 <- transferST(res1)
    }, silent=T)
    powersizeArray[r,1,i] <- res1[1] # standardized
    powersizeArray[r,2,i] <- res1[2] # non-standardized
  }
  for(r in 1:length(G2list)){
    res1 <- c(NA, NA)
    try({res2 <- SILM::ST(dat1$X, dat1$Y, sub.size=n1*subfracVec[jj], G2list[[r]])
    res21 <- transferST(res2)
    }, silent=T)
    powersizeArray[r+3,1,i] <- res21[1]
    powersizeArray[r+3,2,i] <- res21[2]
  }
  
  jj <- 2 # sub.frac = 1/3
  for(r in 1:length(G1list)){
    # r <- 1
    res1 <- c(NA, NA)
    try({res1 <- SILM::ST(dat1$X, dat1$Y, sub.size=n1*subfracVec[jj], G1list[[r]])
    res1 <- transferST(res1)
    }, silent=T)
    powersizeArray[r,3,i] <- res1[1] # standardized
    powersizeArray[r,4,i] <- res1[2] # non-standardized
  }
  for(r in 1:length(G2list)){
    res1 <- c(NA, NA)
    try({res2 <- SILM::ST(dat1$X, dat1$Y, sub.size=n1*subfracVec[jj], G2list[[r]])
    res21 <- transferST(res2)
    }, silent=T)
    powersizeArray[r+3,3,i] <- res21[1]
    powersizeArray[r+3,4,i] <- res21[2]
  }
  
  
  
  toc <- proc.time()
  toc - tic
  
  save(powersizeArray, file=paste0("n", n1, "rho", rho, "_STsubfrac2", ".rds"))
}

toc_total <- proc.time()
time_total <- toc_total[3] - tic_total[3]
save(powersizeArray, time_total,file=paste0("n", n1, "rho", rho, "_STsubfrac2", ".rds"))



#B-13 p-values corrected method  based on the ridge projection method proposed by Bühlmann (2013) --------


powersizeArray <- array(dim=c(6,2, N)) # 6 rows denotes the true hypotheises with first three, and false hypotheises with last three. 
rho <- 0.8; # 0.3, 0.5
for(i in 1:N){
  # i <- 1
  message("i = ", i)
  n1 <- 50
  set.seed(i)
  dat1 <- gendata_Reg(n=n1, p, s0, rho, seed=i)
  res1 <- ridge.proj(dat1$X, dat1$Y)
  pcor_tmp <- res1$pval.corr
  
  for(r in 1:length(G1list)){
    powersizeArray[r,1,i] <- any(pcor_tmp[G1list[[r]]]< alpha)
    
  }
  for(r in 1:length(G2list)){
    powersizeArray[r,2,i] <- any(pcor_tmp[G2list[[r]]]< alpha)
  }
  
  n1 <- 100
  set.seed(i)
  dat1 <- gendata_Reg(n=n1, p, s0, rho, seed=i)
  res1 <- ridge.proj(dat1$X, dat1$Y)
  pcor_tmp <- res1$pval.corr
  
  for(r in 1:length(G1list)){
    powersizeArray[r+3,1,i] <- any(pcor_tmp[G1list[[r]]]< alpha)
    
  }
  for(r in 1:length(G2list)){
    powersizeArray[r+3,2,i] <- any(pcor_tmp[G2list[[r]]]< alpha)
  }
  
}
save(powersizeArray, file=paste0("n2", "rho", rho, "_ridge.rds"))


# MMB-09: $p$-values corrected method based on the sample multi-splitting approach in Meinshausen,  Meier, and Bühlmann (2009) --------------------------------------------



rho <-  0.8 # 0.3 #0.5# 

library(hdi)
p <- 50
s0 <-  5 
G1list <- list(c(p-1, p), (p/2):p, (s0+1):p)
G2list <-  list(c(2, s0+1), c(3, (s0+1):p),c(3,4, (s0+1):p) )


alpha <- 0.05
N <- 500
powersizeArray <- array(dim=c(6,2, N)) 
for(i in 1:N){
  # i <- 1
  message("i = ", i)
  n1 <- 50
  set.seed(i)
  dat1 <- gendata_Reg(n=n1, p, s0, rho, seed=i)
  res_mulsplit <- multi.split(dat1$X, dat1$Y, fraction = 0.7)
  pcor_mulsplit <- res_mulsplit$pval.corr
  
  for(r in 1:length(G1list)){
    powersizeArray[r,1,i] <- any(pcor_mulsplit[G1list[[r]]]< alpha)
  }
  for(r in 1:length(G2list)){
    powersizeArray[r+3,1,i] <- any(pcor_mulsplit[G2list[[r]]]< alpha)
  }
  
  n1 <- 100
  set.seed(i)
  dat1 <- gendata_Reg(n=n1, p, s0, rho, seed=i)
  res_mulsplit <- multi.split(dat1$X, dat1$Y, fraction = 0.7)
  pcor_mulsplit <- res_mulsplit$pval.corr
  
  for(r in 1:length(G1list)){
    powersizeArray[r,2,i] <- any(pcor_mulsplit[G1list[[r]]]< alpha)
  }
  for(r in 1:length(G2list)){
    powersizeArray[r+3,2,i] <- any(pcor_mulsplit[G2list[[r]]]< alpha)
  }
  save(powersizeArray, file=paste0("n2","rho", rho, "_Multi.rds"))
}
save(powersizeArray, file=paste0("n2","rho", rho, "_Multi.rds"))



