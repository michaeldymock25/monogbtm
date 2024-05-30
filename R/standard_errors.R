## functions used to compute standard errors for EM algorithm
## see standard errors chapter for reference

std_errors <- function(data, weights, TT, n_gps, estimates, X_des, covariates, response){  
  
  ## computes standard errors by taking the inverse of the observed information matrix
  
  Jc_mat <- com_inf(weights, TT, n_gps, estimates, X_des, covariates, response)      ## complete data information matrix
  Jm_mat <- cov_score(data, TT, weights, estimates, X_des, covariates, response)     ## covariance of the score vector
  inf_mat <- Jc_mat - Jm_mat                                                         ## observed information matrix
  return(solve(inf_mat))   
}

com_inf <- function(weights, TT, n_gps, estimates, X_des, Z, response){   
  
  ## computes complete data information matrix (block matrix)
  ## note Z is a matrix containing the covariates
  
  if(is.null(Z)){
    com_inf_gpprops_mat <- com_inf_pi(weights, n_gps, estimates[["pis"]])   ## if there are no covariates compute matrix for pis (group proportions)
  } else {                                                                  ## otherwise compute matrix for gammas
    com_inf_gpprops_mat <- matrix(NA, nrow = ((n_gps-1)*(ncol(Z))), ncol = ((n_gps-1)*(ncol(Z))))
    for(j in 2:n_gps){
      for(k in 2:n_gps){
        com_inf_gpprops_mat[((1+ncol(Z)*(j-2)):(ncol(Z)*(j-1))), ((1+ncol(Z)*(k-2)):(ncol(Z)*(k-1)))] <- com_inf_gamma(weights, estimates[["gammas"]], j, k, Z)
      }
    }
  }
  theta_mats <- lapply(1:n_gps, function(gp) com_inf_theta(weights[,gp], TT, X_des, estimates[["thetas"]][[gp]], estimates[["sigma2"]], response)) ## compute matrix for each group (j in 1 to J)
  com_inf_mat <- Matrix::bdiag(list(com_inf_gpprops_mat, Matrix::bdiag(theta_mats)))                                                               ## combine all above matrices as a block diagonal matrix
  return(com_inf_mat)
}

com_inf_pi <- function(weights, n_gps, pis){
  
  ## computes complete data information matrix for pi
  
  pi_sums <- sum(weights[,n_gps]/pis[n_gps]^2)                                ## store outside loop to speed up computations
  mat <- matrix(pi_sums, nrow = n_gps-1, ncol = n_gps-1)                      ## fix off-diagonal elements
  for(j in 1:ncol(mat)) mat[j,j] <- sum(weights[,j]/pis[j]^2) + pi_sums       ## fix diagonal elements
  return(mat)
}

com_inf_gamma <- function(weights, gammas, g1, g2, Z){  
  
  ## computes complete data information matrix for gamma
  
  mat <- matrix(NA, nrow = ncol(Z), ncol = ncol(Z))
  gam_sums <- colSums(exp(gammas%*%t(Z)))      ## store outside loop to speed up computations
  if(g1 == g2){                                ## if g1 = g2 compute I_c(gamma_g1)
    for(j in 1:nrow(mat)){
      for(k in j:ncol(mat)){
        mat[j,k] <- sum(Z[,j]*Z[,k]*exp(gammas[g1,]%*%t(Z))*(gam_sums-exp(gammas[g1,]%*%t(Z)))/gam_sums^2)
        mat[k,j] <- mat[j,k]
      }
    } 
  } else {                                     ## otherwise compute I_c(gamma_g1,gamma_g2)
    for(j in 1:nrow(mat)){
      for(k in j:ncol(mat)){
        mat[j,k] <- -sum(Z[,j]*Z[,k]*exp(gammas[g1,]%*%t(Z))*exp(gammas[g2,]%*%t(Z))/gam_sums^2)
        mat[k,j] <- mat[j,k]
      }
    }
  }
  return(mat)
}

com_inf_theta <- function(weights, TT, X_des, theta, sigma2, response){  
  
  ## computes complete data information matrix for theta
  
  X_des <- X_des[,1:length(theta)] ## adjust design matrix for each group theta (as polynomial degree may differ)
  if(length(theta)!= 1){               
    X_des_temp <- X_des[,2]   
  } else {  
    X_des_temp <- X_des       ## if theta polynomial is only a constant
  }
  mat <- matrix(NA, nrow = length(theta), ncol = length(theta))
  weights_ext <- rep(weights, each = TT)    ## extend weights to match length of design matrix (n x T)
  if(response == "gaussian"){
    for(j in 1:nrow(mat)){
      for(k in j:ncol(mat)){
        mat[j,k] <- 1/sigma2*sum(weights_ext*X_des_temp^(j+k-2))
        mat[k,j] <- mat[j,k]
      }
    }
  } else if(response == "poisson"){
    theta_const<- exp(X_des%*%theta)  ## store outside loop to speed up computations
    for(j in 1:nrow(mat)){
      for(k in j:ncol(mat)){
        mat[j,k] <- sum(weights_ext*X_des_temp^(j+k-2)*theta_const) 
        mat[k,j] <- mat[j,k]
      }
    }
  } else if(response == "bernoulli"){
    theta_const<- exp(X_des%*%theta)  ## store outside loop to speed up computations
    for(j in 1:nrow(mat)){
      for(k in j:ncol(mat)){
        mat[j,k] <- sum(weights_ext*X_des_temp^(j+k-2)*theta_const/(1+theta_const)^2) 
        mat[k,j] <- mat[j,k]
      }
    }
  } 
  return(mat)
}

cov_score <- function(data, TT, n_gps, weights, estimates, X_des, Z, response){  
  
  ## computes covariance of the score vector
  ## note Z is a matrix containing the covariates
  ## note theta contains the polynomial for a group whilst Theta is an object for standard error computation
  
  thetas <- estimates[["thetas"]]                                       ## extract all thetas
  sigma2 <- estimates[["sigma2"]]                                       ## extract estimated variance
  Thetas <- Thetas(data, TT, n_gps, X_des, thetas, sigma2, response)    ## compute matrix containing the Thetas
  if(is.null(Z)){   
    Z <- matrix(NA, ncol = 1)                                                                   ## set ncol(Z) to one for later indexing
    cov_score_mat <- matrix(NA, nrow = (n_gps+length(unlist(thetas))-1),                        ## setup dimensions of block matrix and compute the
                                ncol = (n_gps+length(unlist(thetas))-1))                        ## covariance of the score vector with respect to pi
    cov_score_mat[1:(n_gps-1),1:(n_gps-1)] <- cov_score_pi(estimates[["pis"]], weights, n_gps)  
    count <- 0  
    for(j in 1:n_gps){
      temp_mat <- cov_score_pi_theta(Thetas[[j]], estimates[["pis"]], j, weights, n_gps)        ## compute covariance of the score vector 
      cov_score_mat[1:(n_gps-1), (n_gps+count):(n_gps+count+length(thetas[[j]])-1)] <- temp_mat ## with respect to pi and each group theta
      cov_score_mat[(n_gps+count):(n_gps+count+length(thetas[[j]])-1), 1:(n_gps-1)] <- t(temp_mat)
      count <- count + length(thetas[[j]])
    }
  } else {
    cov_score_mat <- matrix(NA, nrow = ((n_gps-1)*ncol(Z)+length(unlist(thetas))),     ## setup dimensions of block matrix and compute the
                                ncol = ((n_gps-1)*ncol(Z)+length(unlist(thetas))))     ## covariance of the score vector with respect to gamma
    cov_score_mat[1:((n_gps-1)*ncol(Z)),1:((n_gps-1)*ncol(Z))] <- cov_score_gamma(Z, weights, n_gps)
    count <- 0
    for(j in 1:n_gps){
      temp_mat <- cov_score_gamma_theta(Thetas[[j]], weights, Z, j, n_gps)             ## compute covariance of the score vector
      cov_score_mat[1:((n_gps-1)*ncol(Z)),                                             ## with respect to gamma and each group theta
                    (1+(n_gps-1)*ncol(Z)+count):(1+(n_gps-1)*ncol(Z)+count+length(thetas[[j]])-1)] <- temp_mat                                                           
      cov_score_mat[(1+(n_gps-1)*ncol(Z)+count):(1+(n_gps-1)*ncol(Z)+count+length(thetas[[j]])-1),
                    1:((n_gps-1)*ncol(Z))] <- t(temp_mat)
      
      count <- count + length(thetas[[j]])
    }
  }
  count1 <- 0
  for(j in 1:n_gps){
    count2 <- 0
    for(k in 1:n_gps){
      if(j == k){
        cov_score_mat[(1+(n_gps-1)*ncol(Z)+count1):((n_gps-1)*ncol(Z)+count1+length(thetas[[k]])),  ## compute covariance of the score vector with respect to theta
                      (1+(n_gps-1)*ncol(Z)+count1):((n_gps-1)*ncol(Z)+count1+length(thetas[[k]]))] <- cov_score_theta(Thetas[[j]], weights, j)
      } else if(k > j){
        temp_mat <- cov_score_theta_theta(Thetas[[j]], Thetas[[k]], weights, j, k)   ## compute covariance of the score vector with respect to two different groups                                    
        cov_score_mat[(1+(n_gps-1)*ncol(Z)+count1):((n_gps-1)*ncol(Z)+count1+length(thetas[[j]])),
                      (1+(n_gps-1)*ncol(Z)+count2):((n_gps-1)*ncol(Z)+count2+length(thetas[[k]]))] <- temp_mat
        cov_score_mat[(1+(n_gps-1)*ncol(Z)+count2):((n_gps-1)*ncol(Z)+count2+length(thetas[[k]])),
                      (1+(n_gps-1)*ncol(Z)+count1):((n_gps-1)*ncol(Z)+count1+length(thetas[[j]]))] <- t(temp_mat)
      }
      count2 <- count2 + length(thetas[[k]])
    }
    count1 <- count1 + length(thetas[[j]])
  }
  return(cov_score_mat)
}

Thetas <- function(data, TT, n_gps, X_des, thetas, sigma2, response){
  
  ## computes Thetas outside loops (in other functions) to speed up computations
  ## note that Thetas is a list containing an (n x p_j) matrix for each group
  
  Thetas <- list()                         ## intialise list
  obs <- as.vector(t(data))                ## extract observations and store in an (n x T) vector
  for(j in 1:n_gps){
    Theta_mat <- matrix(NA, nrow = nrow(data), ncol = length(thetas[[j]]))
    X_gp <- X_des[,1:length(thetas[[j]])]  ## adjust design matrix for each group theta (as polynomial degree may differ)
    if(length(thetas[[j]])==1){
      X_gp <- as.matrix(X_gp)
    }
    if(response == "gaussian"){
      theta_const <- X_gp%*%thetas[[j]]        ## store outside loop to speed up computations
      Theta_mat <- sapply(1:length(thetas[[j]]), function(p) colSums(matrix(1/sigma2*X_gp[,p]*(obs-theta_const), nrow = TT, ncol = nrow(data))))
    } else if(response == "poisson"){
      theta_const <- exp(X_gp%*%thetas[[j]])   ## store outside loop to speed up computations
      Theta_mat <- sapply(1:length(thetas[[j]]), function(p) colSums(matrix(X_gp[,p]*(obs-theta_const), nrow = TT, ncol = nrow(data))))
    } else if(response == "bernoulli"){
      theta_const <- exp(X_gp%*%thetas[[j]])/(1+exp(X_gp%*%thetas[[j]]))   ## store outside loop to speed up computations
      Theta_mat <- sapply(1:length(thetas[[j]]), function(p) colSums(matrix(X_gp[,p]*(obs-theta_const), nrow = TT, ncol = nrow(data))))
    }
    Thetas[[j]] <- Theta_mat
  }
  return(Thetas)
}

cov_score_pi <- function(pis, weights, n_gps){ 
  
  ## computes covariance of the score vector with respect to pi
  
  mat <- matrix(NA, nrow = n_gps-1, ncol = n_gps-1)
  for(j in 1:ncol(mat)){
    for(k in j:ncol(mat)){
      if(j == k){
        mat[j,j] <- sum(weights[,j]*(1-weights[,j])/pis[j]^2+weights[,n_gps]*(1-weights[,n_gps])/pis[n_gps]^2+2*weights[,j]*weights[,n_gps]/(pis[j]*pis[n_gps]))
      } else {
        mat[j,k] <- sum(-weights[,j]*weights[,k]/(pis[j]*pis[k])+weights[,n_gps]*(1-weights[,n_gps])/pis[n_gps]^2+weights[,j]*weights[,n_gps]/(pis[j]*pis[n_gps])+weights[,k]*weights[,n_gps]/(pis[k]*pis[n_gps]))
        mat[k,j] <- mat[j,k]
      }
    }
  }
  return(mat)
}

cov_score_pi_theta <- function(Theta, pis, gp, weights, n_gps){ 
  
  ## computes covariance of the score vector with respect to pi and each group theta
  
  mat <- matrix(NA, nrow = n_gps-1, ncol = ncol(Theta))
  if(gp == n_gps){
    for(k in 1:nrow(mat)){
      mat[k,] <- t(Theta)%*%(-weights[,gp]*(weights[,k]/pis[k]+(1-weights[,n_gps])/pis[n_gps]))
    }
  } else {
    for(k in 1:nrow(mat)){
      if(k == gp){
        mat[k,] <- t(Theta)%*%(weights[,gp]*((1-weights[,k])/pis[k]+weights[,n_gps]/pis[n_gps]))
      } else {
        mat[k,] <- t(Theta)%*%(weights[,gp]*(-weights[,k]/pis[k]+weights[,n_gps]/pis[n_gps]))
      }
    }
  }
  return(mat)
}

cov_score_gamma <- function(Z, weights, n_gps){ 
  
  ## computes covariance of the score vector with respect to gamma
  ## computes each smaller matrix temp_mat and stores in block matrix mat (note that temp_mat is always symmetric)
  
  mat <- matrix(NA, nrow = ((n_gps-1)*ncol(Z)), ncol = ((n_gps-1)*ncol(Z)))
  vars <- weights*(1-weights)                                
  for(j in 2:n_gps){         
    for(k in j:n_gps){       
      if(j == k){
        temp_mat <- matrix(NA, nrow = ncol(Z), ncol = ncol(Z))
        for(m in 1:nrow(temp_mat)){
          temp_mat[m,] <- colSums(Z[,m]*Z*vars[,j])
        }
        mat[(1+(j-2)*ncol(Z)):((j-1)*ncol(Z)),(1+(j-2)*ncol(Z)):((j-1)*ncol(Z))] <- temp_mat
      } else {
        covs <- -weights[,j]*weights[,k]
        temp_mat <- matrix(NA, nrow = ncol(Z), ncol = ncol(Z))
        for(m in 1:nrow(temp_mat)){
          temp_mat[m,] <- colSums(Z[,m]*Z*covs)
        }
        mat[(1+(j-2)*ncol(Z)):((j-1)*ncol(Z)),(1+(k-2)*ncol(Z)):((k-1)*ncol(Z))] <- temp_mat
        mat[(1+(k-2)*ncol(Z)):((k-1)*ncol(Z)),(1+(j-2)*ncol(Z)):((j-1)*ncol(Z))] <- temp_mat
      }
    }
  }
  return(mat)
}  

cov_score_gamma_theta <- function(Theta, weights, Z, gp, n_gps){ 
  
  ## computes covariance of the score vector with respect to gamma and each group theta
  ## computes each smaller matrix temp_mat and stores in block matrix mat
  
  mat <- matrix(NA, nrow = ((n_gps-1)*ncol(Z)), ncol = ncol(Theta))
  for(j in 2:n_gps){
    temp_mat <- matrix(NA, nrow = ncol(Z), ncol = ncol(Theta))
    if(j == gp){
      vars <- weights[,gp]*(1-weights[,gp])       
      for(m in 1:ncol(temp_mat)){
        temp_mat[,m] <- colSums(Theta[,m]*Z*vars)
      }
    } else {
      covs <- -weights[,j]*weights[,gp]
      for(m in 1:ncol(temp_mat)){
        temp_mat[,m] <- colSums(Theta[,m]*Z*covs)
      }
    }
    mat[(1+(j-2)*ncol(Z)):((j-1)*ncol(Z)),] <- temp_mat
  }
  return(mat)
}  

cov_score_theta <- function(Theta, weights, gp){ 
  
  ## computes covariance of the score vector with respect to theta

  mat <- t(apply(Theta, 2, function(x) colSums(x*Theta*weights[,gp]*(1-weights[,gp]))))
  return(mat)
}  

cov_score_theta_theta <- function(Theta_g1, Theta_g2, weights, gp1, gp2){ 
  
  ## computes covariance of the score vector with respect two different groups
  
  temp <- Theta_g2*weights[,gp1]*weights[,gp2]
  mat <- t(apply(Theta_g1, 2, function(x) -colSums(x*temp)))
  return(mat)
}
