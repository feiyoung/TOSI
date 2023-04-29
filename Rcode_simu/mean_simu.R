


rm(list=ls())
source('definedFunctions.R')

gendata_Mean <- function(n, p, s0= floor(p/2), seed=1, rho= 1, tau=1, cor=0.5){
  mu <- rep(0, p)
  set.seed(1)
  mu[1:s0] <- runif(s0)* rho
  set.seed(seed)
  X <- mvrnorm(n=n, mu=mu, Sigma = tau*cor.mat(p, rho=cor))
  X[,1:(floor(p)/2)] <- - X[,1:(floor(p)/2)]
  mu[1:(floor(p)/2)] <- -mu[1:(floor(p)/2)] 
  return(list(X=X, mu=mu, p0=s0))
}

###  mMeanMax
n <- 100; p <- 100
s0 <- 5  # 5 #  fixed s0
#s0 <- p/2
rho <- 0.5;  # 0.5 or 1
tau <- 1
G1list <- list((p-1):p, (p/2):p,(s0+1):p, c(2, s0+1), c(3, (s0+1):p),c(3,4, (s0+1):p))
N <- 10 ## take 10 for exmaple. set N=500 for formal run.
test_stat_val_mat <- matrix(0, N, length(G1list))
Ns_set <- c(2, 5, 8)# , 15, 20)
iArray <- array(dim=c(N, length(Ns_set)+1, length(G1list)))
for(i in 1:N){
  # i <- 1
  cat('i=', i, '\n')
  dat1 <- gendata_Mean(n, p, s0, seed=i, rho, tau,cor=0.9)
  
  # Choose different split times
  for(js in 1:length(Ns_set)){
    iArray[i,js,] <-  getmultiMeanMax2(dat1$X, G1list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }
  
  # max-type test
  iArray[i,js+1,] <- getmaxMeantest(dat1$X, G1list)[3,]
  iArray[i,,]
}
resultTest <- apply(iArray, c(2,3), mean)
row.names(resultTest) <- c(paste0("split=",Ns_set), 'max-type')
colnames(resultTest) <- paste0('G1', 1:6)
resultTest
simutool::latexTable(resultTest, paraname = NULL,digits = 3)


### mMeanMin
n <- 200; ## 100 or 200
p <- 100
s0 <- 5 # fixed s0
#s0 <- 10
rho <- 1;  # 0.5 0r 1
tau <- 1;
# G2list <-list((p-1):p, 1:p, 1:(s0+3), 1:(s0+1))
G2list <-list((p-4):p, (s0+1):p, 1:p, 1:2, 1:4, 1:s0)
N <- 10
Ns_set <- c(2, 5, 8, 15, 20)
iArray <- array(dim=c(N, length(Ns_set), length(G2list)))
for(i in 1:N){
  
  cat('i=', i, '\n')
  dat1 <- gendata_Mean(n, p, s0, seed=i, rho, tau)
  for(js in 1:length(Ns_set)){
    iArray[i,js,] <-   getmultiMeanMin2(dat1$X, G2list, standardize=T, Nsplit = Ns_set[js], r=1/Ns_set[js], BH=T)[2,]
  }

}
resultTest <- apply(iArray, c(2,3), mean)
row.names(resultTest) <- paste0("split=",Ns_set)
colnames(resultTest) <- paste0('G2', 1:6)
resultTest


simutool::latexTable(resultTest, paraname = NULL,digits = 3)
