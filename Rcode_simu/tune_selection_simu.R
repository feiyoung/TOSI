
### Compare the tuning parameter selection------------------------------

source("definedFunctions.R")
n <- 50;  # 50 or 100
p <- 50
s0 <-  3#
rho <-  3 # 2 or 3
i <- 1
dat <- gendata_Reg(50, p, s0, rho, seed=i)
dat1 <- gendata_Reg(n, p, s0, rho, seed=i)



# Repeat
N <- 500
varsArray <- array(dim=c(5, 4, N))
for(i in 1:N){
  # i <- 1
  cat('i=', i, '\n')
  
  dat1 <- gendata_Reg(n, p, s0, seed=i, rho)
  try({
    varsArray[,,i] <- vars_comp2(dat1$X, dat1$Y,dat$X,dat$Y, dat1$index_nz, alpha=0.05, seed=1, sub.frac=0.5, standardized=T)
    
  }, silent = T)
  
}

save(varsArray, file = paste0("n",n, "rho", rho, "CompareMethod5.rds"))





# Test the sensitivity to the significance levels of TOSI -------------------------


n_set <- c(50, 100)
n <- 50; # 100
p <- 50
s0 <-  3#
rho <-  3 #3
i <- 1
dat <- gendata_Reg(50, p, s0, rho, seed=i)
dat1 <- gendata_Reg(n, p, s0, rho, seed=i)

## 
alpha_set <- c(0.1, 0.05, 0.01)
n_alpha <- length(alpha_set)


N <- 500
varsArray_tosi <- array(dim=c(n_alpha, 4, N))
for(i in 1:N){
  # i <- 1
  cat('i=', i, '\n')
  
  dat1 <- gendata_Reg(n, p, s0, seed=i, rho)
  
  for(j in 1:n_alpha){
    # j <- 1
    message("alpha = ", alpha_set[j])
    try({
      varsArray_tosi[j,,i] <- TOSI_selectLambda(dat1$X, dat1$Y,dat$X,dat$Y, dat1$index_nz, 
                                                alpha=alpha_set[j], seed=1, sub.frac=0.5, standardized=T)
      
    }, silent = T)
  }
  
  
}


save(varsArray_tosi, file = paste0("./Rdata_R2/n",n, "rho", rho, "_selectLambda_TOSI_Alpha3.rds"))
