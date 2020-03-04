## this script contains examples of fitting group based trajectory models with gbtm_monotone

sim_data <- function(n, dim, x, group_mem, thetas, sd, response){  
  
  ## this function simulates a data set
  
  dat <- matrix(NA, nrow = n, ncol = dim)         ## matrix to contain simulated observations
  if(is.vector(x)){                                                ## if time-varying covariate is the same for all individuals
    x_scl <- (x - mean(x))/sd(x)                                   ## scale the time varying covariate for computation
    X_des <- matrix(t(matrix(x_scl, nrow = dim,                    ## create design matrix X_des
                             ncol = (ncol(thetas))))^(0:(ncol(thetas)-1)), nrow = n*dim, ncol = ncol(thetas), byrow = T)  
  } else if(is.matrix(x)){                                         ## if time-varying covariate is not the same for all individuals
    X_des_unscl <- as.vector(t(x))                                 ## start creating design matrix X_des (unscaled)
    x_unique <- unique(X_des_unscl)
    x_scl <- (x_unique - mean(x_unique))/sd(x_unique)              ## scale the unique times
    X_des <- vector(length = length(X_des_unscl))                  ## set up scaled design matrix
    for(i in 1:length(x_scl)){
      X_des[X_des_unscl == sort(x_unique)[i]] <- x_scl[i]          ## set each element to appropriate scaled time
    }
    X_des <- t(t(matrix(X_des, nrow = length(X_des), ncol = ncol(thetas)))^(0:(ncol(thetas)-1)))  
  }
  order <- sample(1:n, n)            ## randomly order individuals to ensure proper mixing between simulated data sets 
  csums <- cumsum(group_mem)         ## (only required if time-varying covariates are different between individuals
  if(response == "gaussian"){        ## simulate gaussian responses
     r_samps <- rnorm(n*dim, mean = 0, sd = sd)
     count <- 0
     j <- 1
     for(i in 1:n){
       dat[i,] <- X_des[(1+dim*(order[i]-1)):(dim*order[i]),]%*%thetas[j,] + r_samps[(1+dim*(i-1)):(dim*i)]
       if(count >= csums[j]){
         j <- j + 1
       }
       count <- count + 1
     }
  }
  else if(response == "poisson"){    ## simulate count responses
    count <- 0
    j <- 1
    for(i in 1:n){
      dat[i,] <- rpois(dim, lambda = exp(X_des[(1+dim*(order[i]-1)):(dim*order[i]),]%*%thetas[j,]))
      if(count >= csums[j]){
        j <- j + 1
      }
      count <- count + 1
    }
  }
  else if(response == "bernoulli"){  ## simulate binary responses
    count <- 0
    j <- 1
    for(i in 1:n){
      p <- exp(X_des[(1+dim*(order[i]-1)):(dim*order[i]),]%*%thetas[j,])/(1+exp(X_des[(1+dim*(order[i]-1)):(dim*order[i]),]%*%thetas[j,]))
      dat[i,] <- rbinom(dim, size = 1, prob = p)
      if(count >= csums[j]){
        j <- j + 1
      }
      count <- count + 1
    }
  }
  return(dat)
}

## gaussian response example without monotonicity constraints (variations for design matrix, covariates and polynomial degrees)

n_ex1 <- 800                                                 ## number of individuals n
dim_ex1 <- 5                                                 ## number of time points T
x_ex1.1 <- 1:5                                               ## setup simple design matrix
x_ex1.2 <- matrix(1:5, nrow = 800, ncol = 5, byrow = T)      ## setup complex design matrix
x_ex1.2[sample(800, 100),] <- c(3:7)
group_mem_ex1 <- c(200, 300, 150, 150)                       ## number of individuals belonging to each group
thetas_ex1 <- matrix(c(10,0,0,0,13,-3,-2,0,5,2,-0.5,1,1,0.4, ## parameters governing group trajectories
                       1.3,0), nrow = 4, ncol = 4, byrow = T)    
sd_ex1 <- 2                                                  ## standard deviation
poly_degs_ex1.1 <- c(3,3,3,3)                                ## degrees for polynomials
poly_degs_ex1.2 <- c(0,2,2,3)
covs_ex1 <- cbind(1, matrix(rpois(800, lambda = 5), nrow = 800, ncol = 1))   ## random covariates

dat_ex1.1 <- sim_data(n_ex1, dim_ex1, x_ex1.1, group_mem_ex1, thetas_ex1, sd_ex1, 'gaussian')   ## simulate data sets under both settings
dat_ex1.2 <- sim_data(n_ex1, dim_ex1, x_ex1.2, group_mem_ex1, thetas_ex1, sd_ex1, 'gaussian')

## fit models under each setup

gbtm_monotone(dat_ex1.1, 4, x_ex1.1, poly_degs_ex1.1, monotone = FALSE)       ## simple design, no covariates, all cubics
gbtm_monotone(dat_ex1.2, 4, x_ex1.2, poly_degs_ex1.1, monotone = FALSE)       ## complex design, no covariates, all cubics
gbtm_monotone(dat_ex1.1, 4, x_ex1.1, poly_degs_ex1.2, monotone = FALSE)       ## simple design, no covariates, differing degrees
gbtm_monotone(dat_ex1.1, 4, x_ex1.1, poly_degs_ex1.1, monotone = FALSE, covariates = covs_ex1) ## simple design, covariates, all cubics

## gaussian response example with monotonicity constraints (increasing) (variations for covariates, polynomial degrees and confidence bands)

n_ex2 <- 600                                                 ## number of individuals n
dim_ex2 <- 6                                                 ## number of time points T
x_ex2 <- 1:6                                                 ## setup simple design matrix
group_mem_ex2 <- c(200, 200, 200)                            ## number of individuals belonging to each group
thetas_ex2 <- matrix(c(10,1,0,1,11,2,0.1,2,13,3,0.2,3),      ## parameters governing group trajectories
                       nrow = 3, ncol = 4, byrow = T)    
sd_ex2 <- 2.5                                                ## standard deviation
poly_degs_ex2.1 <- c(3,3,3)                                  ## degrees for polynomials
poly_degs_ex2.2 <- c(1,3,5)
covs_ex2 <- cbind(1, matrix(rpois(600, lambda = 10), nrow = 600, ncol = 1))   ## random covariates

dat_ex2 <- sim_data(n_ex2, dim_ex2, x_ex2, group_mem_ex2, thetas_ex2, sd_ex2, 'gaussian')   ## simulate data set

## fit models under each setup

gbtm_monotone(dat_ex2, 3, x_ex2, poly_degs_ex2.1)       ## no covariates, all cubics
gbtm_monotone(dat_ex2, 3, x_ex2, poly_degs_ex2.2)       ## no covariates, differing degrees
gbtm_monotone(dat_ex2, 3, x_ex2, poly_degs_ex2.1, covariates = covs_ex2) ## covariates, all cubics
gbtm_monotone(dat_ex2, 3, x_ex2, poly_degs_ex2.1, conf_level = 95)       ## no covariates, all cubics, confidence bands 95%

## gaussian response example with monotonicity constraints (decreasing) (variations for polynomial degrees)

n_ex3 <- 500                                                 ## number of individuals n
dim_ex3 <- 6                                                 ## number of time points T
x_ex3 <- 1:6                                                 ## setup simple design matrix
group_mem_ex3 <- c(100, 200, 100, 100)                       ## number of individuals belonging to each group
thetas_ex3 <- matrix(c(10,-1,-0.5,-2,0,0,12,-2,-1,-3,0,0,14, ## parameters governing group trajectories
                       -3,-1.5,-4,0,0,16,-4,-2,-4,3,-3), nrow = 4, ncol = 6, byrow = T)    
sd_ex3 <- 3                                                  ## standard deviation
poly_degs_ex3.1 <- c(3,3,3,5)                                ## degrees for polynomials
poly_degs_ex3.2 <- c(1,3,5,5)

dat_ex3 <- sim_data(n_ex3, dim_ex3, x_ex3, group_mem_ex3, thetas_ex3, sd_ex3, 'gaussian')   ## simulate data set

## fit models under each setup

gbtm_monotone(dat_ex3, 4, x_ex3, poly_degs_ex3.1)       ## no covariates, four cubics one quintic
gbtm_monotone(dat_ex3, 4, x_ex3, poly_degs_ex3.2)       ## no covariates, differing degrees

