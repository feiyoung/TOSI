

rm(list=ls())
source('definedFunctions.R')



# Max-test for TRUE H0--------------------------------------------------------------
###  mFacMax
set.seed(1)
n <- 200; p = 150
q <- 1
datlist1 <- gendatass2(n= n, p = p, q=1,sigma2=2, rho=0.8, seed= 1)
s0 <- tail(datlist1$ind_nz,1)
ind_z <- setdiff(1:p, datlist1$ind_nz)

G1list <- list(c(p-1,p), floor(0.9*p):p,ind_z)
N <- 500
test_stat_val_mat <- matrix(0, N, length(G1list))
Ns_set <- c(2, 5, 8, 15, 20)
iArray <- array(dim=c(N, length(Ns_set), length(G1list)))
for(i in 1:N){
  # i <- 1
  cat('i=', i, '\n')
  dat1 <- gendatass2(n= n, p = p, q=q,sigma2=1, rho=0.3, seed= i)
  
  # Choose different split times
  for(js in 1:length(Ns_set)){
    iArray[i,js, ] <- getmultiRowMax(dat1$X, G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
}
resultTest1 <- apply(iArray, c(2,3), mean)
row.names(resultTest1) <- c(paste0("split=",Ns_set))
colnames(resultTest1) <- paste0('G1', 1:3)
resultTest1




set.seed(1)
n <- 400; p = 150
datlist1 <- gendatass2(n= n, p = p, q=1,sigma2=2, rho=0.8, seed= 1)
s0 <- tail(datlist1$ind_nz,1)
ind_z <- setdiff(1:p, datlist1$ind_nz)

G1list <- list(c(p-1,p), floor(0.9*p):p,ind_z)
N <- 500
test_stat_val_mat <- matrix(0, N, length(G1list))
Ns_set <- c(2, 5, 8, 15, 20)
iArray <- array(dim=c(N, length(Ns_set), length(G1list)))
for(i in 1:N){
  # i <- 1
  cat('i=', i, '\n')
  dat1 <- gendatass2(n= n, p = p, q=q,sigma2=1, rho=0.3, seed= i)
  
  # Choose different split times
  for(js in 1:length(Ns_set)){
    iArray[i,js, ] <-  getmultiRowMax(dat1$X, G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  iArray[i,,]
}
resultTest2 <- apply(iArray, c(2,3), mean)
row.names(resultTest2) <- c(paste0("split=",Ns_set))
colnames(resultTest2) <- paste0('G1', 1:3)
resultTest2


# Max-test for False H0--------------------------------------------------------------
### False H0;
## setting 2：n = 100
set.seed(1)
n <- 100; p = 150
datlist1 <- gendatass2(n= n, p = p, q=1,sigma2=2, rho=0.8, seed= 1)
s0 <- tail(datlist1$ind_nz,1)
ind_z <- setdiff(1:p, datlist1$ind_nz)
G1list <- list(c(2,s0+1), c(3, ind_z),
               c(3,4, ind_z))
N <- 500
test_stat_val_mat <- matrix(0, N, length(G1list))
Ns_set <- c(2, 5, 8, 15, 20)
iArray <- array(dim=c(N, length(Ns_set), length(G1list)))
for(i in 1:N){
  # i <- 1
  cat('i=', i, '\n')
  dat1 <- gendatass2(n= n, p = p, q=q,sigma2=3, rho=0.3, seed= i)
  
  # Choose different split times
  for(js in 1:length(Ns_set)){
    iArray[i,js, ] <- 0 # getmultiRowMax(dat1$X, G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  iArray[i,js+1, ] <- getMaxrowST_BC(dat1$X, G1list, p.adjust.methods="bonferroni")
  iArray[i,js+2, ] <- getMaxrowST_BC(dat1$X, G1list, p.adjust.methods="BY")
}
resultTest1 <- apply(iArray, c(2,3), mean)
row.names(resultTest1) <- c(paste0("split=",Ns_set))
colnames(resultTest1) <- paste0('G1', 4:6)
resultTest1


## setting 2： n=200
n <- 200
G1list <- list(c(2,s0+1), c(3, ind_z),
               c(3,4, ind_z))
N <- 500
test_stat_val_mat <- matrix(0, N, length(G1list))
Ns_set <- c(2, 5, 8, 15, 20)
iArray <- array(dim=c(N, length(Ns_set), length(G1list)))
for(i in 1:N){
  # i <- 1
  cat('i=', i, '\n')
  dat1 <- gendatass2(n= n, p = p, q=q,sigma2=3, rho=0.3, seed= i)
  
  # Choose different split times
  for(js in 1:length(Ns_set)){
    iArray[i,js, ] <-  getmultiRowMax(dat1$X, G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  iArray[i,,]
}
resultTest2 <- apply(iArray, c(2,3), mean)
row.names(resultTest2) <- c(paste0("split=",Ns_set))
colnames(resultTest2) <- paste0('G1', 4:6)
resultTest2


# Min-test for TRUE H0--------------------------------------------------------------
### 
set.seed(1)
n <- 200; p = 150
datlist1 <- gendatass2(n= n, p = p, q=1,sigma2=2, rho=0.8, seed= 1)
s0 <- tail(datlist1$ind_nz,1)
ind_z <- setdiff(1:p, datlist1$ind_nz)
G2list <- list(c(p-2, p-1), (s0+1):p, 1:p)

N <- 500
test_stat_val_mat <- matrix(0, N, length(G1list))
Ns_set <- c(2, 5, 8, 15, 20)
iArray <- array(dim=c(N, length(Ns_set), length(G1list)))
for(i in 1:N){
  cat('i=', i, '\n')
  dat1 <- gendatass2(n= n, p = p, q=1,sigma2=1, rho=1, seed= i)
  
  # Choose different split times
  for(js in 1:length(Ns_set)){
    iArray[i,js, ] <- 0 # getmultiRowMin(dat1$X, G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
}
resultTest1 <- apply(iArray, c(2,3), mean)
row.names(resultTest1) <- c(paste0("split=",Ns_set))
colnames(resultTest1) <- paste0('G2', 1:3)
resultTest1




n <- 400
N <- 500
test_stat_val_mat <- matrix(0, N, length(G1list))
Ns_set <- c(2, 5, 8, 15, 20)
iArray <- array(dim=c(N, length(Ns_set), length(G1list)))
for(i in 1:N){
  cat('i=', i, '\n')
  dat1 <- gendatass2(n= n, p = p, q=q,sigma2=1, rho=1, seed= i)
  
  # Choose different split times
  for(js in 1:length(Ns_set)){
    iArray[i,js, ] <- getmultiRowMin(dat1$X, G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
}
resultTest2 <- apply(iArray, c(2,3), mean)
row.names(resultTest2) <- c(paste0("split=",Ns_set))
colnames(resultTest2) <- paste0('G2', 1:3)
resultTest2


# Min-test for False H0--------------------------------------------------------------
### False H0;
## setting 2：n = 100
set.seed(1)
n <- 100; p = 150
datlist1 <- gendatass2(n= n, p = p, q=1,sigma2=2, rho=0.8, seed= 1)
s0 <- tail(datlist1$ind_nz,1)
ind_z <- setdiff(1:p, datlist1$ind_nz)
G2list <- list(1:2, 1:4, datlist1$ind_nz)
N <- 500
q <- 1
test_stat_val_mat <- matrix(0, N, length(G2list))
Ns_set <- c(2, 5, 8, 15, 20)
iArray <- array(dim=c(N, length(Ns_set), length(G2list)))
for(i in 1:N){
  # i <- 1
  cat('i=', i, '\n')
  dat1 <- gendatass2(n= n, p = p, q=q,sigma2=4, rho=0.3, seed= i)
  
  # Choose different split times
  for(js in 1:length(Ns_set)){
    iArray[i,js, ] <-  getmultiRowMin(dat1$X, G2list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  iArray[i,,]
}
resultTest1 <- apply(iArray, c(2,3), mean)
row.names(resultTest1) <- c(paste0("split=",Ns_set))
colnames(resultTest1) <- paste0('G2', 4:6)
resultTest1


## setting 2： n=200
n <- 200
N <- 500
test_stat_val_mat <- matrix(0, N, length(G2list))
Ns_set <- c(2, 5, 8, 15, 20)
iArray <- array(dim=c(N, length(Ns_set), length(G2list)))
for(i in 1:N){
  # i<- 1
  cat('i=', i, '\n')
  dat1 <- gendatass2(n= n, p = p, q=q,sigma2=4, rho=0.3, seed= i)
  
  # Choose different split times
  for(js in 1:length(Ns_set)){
    iArray[i,js, ] <- getmultiRowMin(dat1$X, G2list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  iArray[i,,]
}
resultTest2 <- apply(iArray, c(2,3), mean)
row.names(resultTest2) <- c(paste0("split=",Ns_set))
colnames(resultTest2) <- paste0('G2', 4:6)
resultTest2

