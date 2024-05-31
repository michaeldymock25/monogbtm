## this script contains examples of fitting group based trajectory models with gbtm_monotone

## this function simulates a data set

sim_data <- function(n, TT, x, group_mem, thetas, sd){  
  dat <- matrix(NA, nrow = n, ncol = TT)                           ## matrix to contain simulated observations
  p <- ncol(thetas)                                                ## maximum polynomial degree
  x_scl <- (x - mean(x))/sd(x)                                     ## scale the time varying covariate for computation
  X_des <- matrix(t(matrix(x_scl, nrow = TT, ncol = p))^(0:(p-1)), ## create design matrix X_des
                  nrow = n*TT, ncol = p, byrow = T)
  order <- sample(1:n, n)            
  csums <- cumsum(group_mem)         
  r_samps <- rnorm(n*TT, mean = 0, sd = sd)
  count <- 0
  j <- 1
  for(i in 1:n){
    dat[i,] <- X_des[(1+TT*(order[i]-1)):(TT*order[i]),]%*%thetas[j,] + r_samps[(1+TT*(i-1)):(TT*i)]
    if(count >= csums[j]) j <- j + 1
    count <- count + 1
  }
  return(dat)
}

## example without monotonicity constraints

n_ex1 <- 800                                               ## number of individuals n
TT_ex1 <- 5                                                ## number of time points TT
J_ex1 <- 4                                                 ## number of groups J
x_ex1 <- 1:TT_ex1                                          ## setup simple design matrix
group_mem_ex1 <- c(200, 300, 150, 150)                     ## number of individuals belonging to each group
thetas_ex1 <- matrix(c(10,   0,    0, 0,                   ## parameters governing group trajectories
                       13,  -3,   -2, 0,
                        5,   2, -0.5, 1,
                        1, 0.4,  1.3, 0),
                     nrow = J_ex1, ncol = 4, byrow = TRUE)    
sd_ex1 <- 2                                                ## standard deviation
poly_degs_ex1 <- c(0, 2, 2, 3)                             ## degrees for polynomials
covs_ex1 <- cbind(rep(1, n_ex1), rpois(n_ex1, lambda = 5)) ## random covariates
                  
dat_ex1 <- sim_data(n_ex1, TT_ex1, x_ex1, group_mem_ex1, thetas_ex1, sd_ex1)

## fit models under each setup

gbtm_monotone(dat_ex1, J_ex1, x_ex1, poly_degs_ex1, monotone = FALSE)                        ## no covariates
gbtm_monotone(dat_ex1, J_ex1, x_ex1, poly_degs_ex1, monotone = FALSE, covariates = covs_ex1) ## covariates

## example with monotonicity constraints (increasing)

n_ex2 <- 600                                                ## number of individuals n
TT_ex2 <- 6                                                 ## number of time points TT
J_ex2 <- 3                                                  ## number of groups J
x_ex2 <- 1:TT_ex2                                           ## setup simple design matrix
group_mem_ex2 <- c(200, 200, 200)                           ## number of individuals belonging to each group
thetas_ex2 <- matrix(c(10, 1,   0, 1,                       ## parameters governing group trajectories
                       11, 2, 0.1, 2,
                       13, 3, 0.2, 3),      
                     nrow = J_ex2, ncol = 4, byrow = TRUE)    
sd_ex2 <- 2.5                                               ## standard deviation
poly_degs_ex2 <- c(3, 3, 3)                                 ## degrees for polynomials
covs_ex2 <- cbind(rep(1, n_ex2), rpois(n_ex2, lambda = 10)) ## random covariates

dat_ex2 <- sim_data(n_ex2, TT_ex2, x_ex2, group_mem_ex2, thetas_ex2, sd_ex2)   

## fit models under each setup

gbtm_monotone(dat_ex2, J_ex2, x_ex2, poly_degs_ex2)                        ## no covariates
gbtm_monotone(dat_ex2, J_ex2, x_ex2, poly_degs_ex2, covariates = covs_ex2) ## covariates
gbtm_monotone(dat_ex2, J_ex2, x_ex2, poly_degs_ex2, conf_level = 95)       ## no covariates, 95% confidence bands

## example with monotonicity constraints (decreasing)

n_ex3 <- 500                                                 ## number of individuals n
TT_ex3 <- 6                                                  ## number of time points T
J_ex3 <- 4                                                   ## number of groups J
x_ex3 <- 1:TT_ex3                                            ## setup simple design matrix
group_mem_ex3 <- c(100, 200, 100, 100)                       ## number of individuals belonging to each group
thetas_ex3 <- matrix(c(10, -1, -0.5, -2, 0,  0,              ## parameters governing group trajectories
                       12, -2,   -1, -3, 0,  0,
                       14, -3, -1.5, -4, 0,  0,
                       16, -4,   -2, -4, 3, -3), 
                     nrow = J_ex3, ncol = 6, byrow = TRUE)    
sd_ex3 <- 3                                                  ## standard deviation
poly_degs_ex3 <- c(3,3,3,5)                                  ## degrees for polynomials

dat_ex3 <- sim_data(n_ex3, TT_ex3, x_ex3, group_mem_ex3, thetas_ex3, sd_ex3) 
gbtm_monotone(dat_ex3, J_ex3, x_ex3, poly_degs_ex3)       

