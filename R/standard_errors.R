## functions used to compute standard errors for EM algorithm
## see standard errors chapter for reference

std_errors <- function(data, dim, estimates, X_des, covariates, response){  
  
  ## computes standard errors by taking the inverse of the observed information matrix
  
  weights <- data[,(dim+1):ncol(data)]                                                ## extract weights
  Jc_mat <- com_inf(weights, dim, estimates, X_des, covariates, response)             ## complete data information matrix
  Jm_mat <- cov_score(data, dim, weights, estimates, X_des, covariates, response)     ## covariance of the score vector
  inf_mat <- Jc_mat - Jm_mat                                                          ## observed information matrix
  return(solve(inf_mat))   
}

com_inf <- function(weights, dim, estimates, X_des, Z, response){   

  ## computes complete data information matrix (block matrix)
  ## note Z is a matrix containing the covariates
  
  if(is.null(Z)){
   com_inf_gpprops_mat <- com_inf_pi(weights, estimates[[1]])   ## if there are no covariates compute matrix for pis (group proportions)
  } else {                                                      ## otherwise compute matrix for gammas
   com_inf_gpprops_mat <- matrix(NA, nrow = ((ncol(weights)-1)*(ncol(Z))), ncol = ((ncol(weights)-1)*(ncol(Z))))
   for(j in 2:ncol(weights)){
     for(k in 2:ncol(weights)){
       com_inf_gpprops_mat[((1+ncol(Z)*(j-2)):(ncol(Z)*(j-1))),
                           ((1+ncol(Z)*(k-2)):(ncol(Z)*(k-1)))] <- com_inf_gamma(weights, estimates[[4]], j, k, Z)
     }
   }
  }
  theta_mats <- list()
  for(j in 1:ncol(weights)){                                          ## compute matrix for each group (j in 1 to J)
   theta_mats[[j]] <- com_inf_theta(weights[,j], dim, X_des, estimates[[2]][[j]], estimates[[3]], response)
  }
  com_inf_mat <- Matrix::bdiag(list(com_inf_gpprops_mat, Matrix::bdiag(theta_mats)))  ## combine all above matrices as a block diagonal matrix
  return(com_inf_mat)
}

com_inf_pi <- function(weights, pis){
  
  ## computes complete data information matrix for pi
  
  pi_sums <- sum(weights[,ncol(weights)]/pis[length(pis)]^2)                  ## store outside loop to speed up computations
  mat <- matrix(pi_sums, nrow = (ncol(weights)-1), ncol = (ncol(weights)-1))  ## fix off-diagonal elements
  for(j in 1:ncol(mat)){
    mat[j,j] <- sum(weights[,j]/pis[j]^2) + pi_sums                           ## fix diagonal elements
  }
  return(mat)
}

com_inf_gamma <- function(weights, gammas, g1, g2, Z){  
  
  ## computes complete data information matrix for gamma
  
  mat <- matrix(NA, nrow = ncol(Z), ncol = ncol(Z))
  gam_sums <- colSums(exp(gammas%*%t(Z)))      ## store outside loop to speed up computations
  if(g1 == g2){                                ## if g1 = g2 compute I_c(gamma_g1)
    temp_g1 <- exp(gammas[g1,]%*%t(Z))
    for(j in 1:nrow(mat)){
      for(k in j:ncol(mat)){
        mat[j,k] <- sum(Z[,j]*Z[,k]*temp_g1*(gam_sums-temp_g1)/gam_sums^2)
        mat[k,j] <- mat[j,k]
      }
    } 
  } else {                                     ## otherwise compute I_c(gamma_g1,gamma_g2)
    temp_g1 <- exp(gammas[g1,]%*%t(Z))
    temp_g2 <- exp(gammas[g2,]%*%t(Z))
    for(j in 1:nrow(mat)){
      for(k in j:ncol(mat)){
        mat[j,k] <- -sum(Z[,j]*Z[,k]*tempg1*tempg2/gam_sums^2)
        mat[k,j] <- mat[j,k]
      }
    }
  }
  return(mat)
}

com_inf_theta <- function(weights, dim, X_des, theta, sigma2, response){  
  
  ## computes complete data information matrix for theta
  
  X_des <- X_des[,1:length(theta)]    ## adjust design matrix for each group theta (as polynomial degree may differ)
  ind <- 2                            ## set an indicator variable
  if(length(theta)==1){               ## change indicator if theta polynomial is only a constant
    X_des <- as.matrix(X_des)
    ind <- 1
  }
  mat <- matrix(NA, nrow = length(theta), ncol = length(theta))
  weights_ext <- rep(weights, each = dim)    ## extend weights to match length of design matrix (n x T)
  if(response == "gaussian"){
    for(j in 1:nrow(mat)){
      for(k in j:ncol(mat)){
        mat[j,k] <- 1/sigma2*sum(weights_ext*X_des[,ind]^(j+k-2))
        mat[k,j] <- mat[j,k]
      }
    }
  } else if(response == "poisson"){
    theta_const<- exp(X_des%*%theta)  ## store outside loop to speed up computations
    for(j in 1:nrow(mat)){
      for(k in j:ncol(mat)){
        mat[j,k] <- sum(weights_ext*X_des[,ind]^(j+k-2)*theta_const) 
        mat[k,j] <- mat[j,k]
      }
    }
  } else if(response == "bernoulli"){
    theta_const<- exp(X_des%*%theta)  ## store outside loop to speed up computations
    for(j in 1:nrow(mat)){
      for(k in j:ncol(mat)){
        mat[j,k] <- sum(weights_ext*X_des[,ind]^(j+k-2)*theta_const/(1+theta_const)^2) 
        mat[k,j] <- mat[j,k]
      }
    }
  } 
  return(mat)
}

cov_score <- function(data, dim, weights, estimates, X_des, Z, response){  
  
  ## computes covariance of the score vector
  ## note Z is a matrix containing the covariates
  ## note theta contains the polynomial for a group whilst Theta is an object for standard error computation
  
  thetas <- estimates[[2]]                                       ## extract all thetas
  sigma2 <- estimates[[3]]                                       ## extract estimated variance
  Thetas <- Thetas(data, dim, X_des, thetas, sigma2, response)   ## compute matrix containing the Thetas
  if(is.null(Z)){
    Z <- matrix(NA, ncol = 1)                                                       ## set ncol(Z) to one for later indexing
    cov_score_mat <- matrix(NA, nrow = (ncol(weights)+length(unlist(thetas))-1),    ## setup dimensions of block matrix and compute the
                            ncol = (ncol(weights)+length(unlist(thetas))-1))        ## covariance of the score vector with respect to pi
    cov_score_mat[1:(ncol(weights)-1),1:(ncol(weights)-1)] <- cov_score_pi(estimates[[1]], weights)  
    count <- 0  
    for(j in 1:ncol(weights)){
      temp_mat <- cov_score_pi_theta(Thetas[[j]], estimates[[1]], j, weights)   ## compute covariance of the score vector 
      cov_score_mat[1:(ncol(weights)-1),                                        ## with respect to pi and each group theta
                    (ncol(weights)+count):(ncol(weights)+count+length(thetas[[j]])-1)] <- temp_mat
      cov_score_mat[(ncol(weights)+count):(ncol(weights)+count+length(thetas[[j]])-1),
                   1:(ncol(weights)-1)] <- t(temp_mat)
      count <- count + length(thetas[[j]])
    }
  } else {
    cov_score_mat <- matrix(NA, nrow = ((ncol(weights)-1)*ncol(Z)+length(unlist(thetas))), ## setup dimensions of block matrix and compute the
                            ncol = ((ncol(weights)-1)*ncol(Z)+length(unlist(thetas))))     ## covariance of the score vector with respect to gamma
    cov_score_mat[1:((ncol(weights)-1)*ncol(Z)),1:((ncol(weights)-1)*ncol(Z))] <- cov_score_gamma(Z, weights)
    count <- 0
    for(j in 1:ncol(weights)){
      temp_mat <- cov_score_gamma_theta(Thetas[[j]], weights, Z, j)   ## compute covariance of the score vector
      cov_score_mat[1:((ncol(weights)-1)*ncol(Z)),                    ## with respect to gamma and each group theta
                    (1+(ncol(weights)-1)*ncol(Z)+count):(1+(ncol(weights)-1)*ncol(Z)+count+length(thetas[[j]])-1)] <- temp_mat                                                           
      cov_score_mat[(1+(ncol(weights)-1)*ncol(Z)+count):(1+(ncol(weights)-1)*ncol(Z)+count+length(thetas[[j]])-1),
                    1:((ncol(weights)-1)*ncol(Z))] <- t(temp_mat)
      
      count <- count + length(thetas[[j]])
    }
  }
  count1 <- 0
  for(j in 1:ncol(weights)){
    count2 <- 0
    for(k in 1:ncol(weights)){
      if(j == k){
        cov_score_mat[(1+(ncol(weights)-1)*ncol(Z)+count1):((ncol(weights)-1)*ncol(Z)+count1+length(thetas[[k]])),  ## compute covariance of the score vector with respect to theta
                      (1+(ncol(weights)-1)*ncol(Z)+count1):((ncol(weights)-1)*ncol(Z)+count1+length(thetas[[k]]))] <- cov_score_theta(Thetas[[j]], weights, j)
      } else if(k > j){
        temp_mat <- cov_score_theta_theta(Thetas[[j]], Thetas[[k]], weights, j, k)   ## compute covariance of the score vector with respect to two different groups                                    
        cov_score_mat[(1+(ncol(weights)-1)*ncol(Z)+count1):((ncol(weights)-1)*ncol(Z)+count1+length(thetas[[j]])),
                      (1+(ncol(weights)-1)*ncol(Z)+count2):((ncol(weights)-1)*ncol(Z)+count2+length(thetas[[k]]))] <- temp_mat
        cov_score_mat[(1+(ncol(weights)-1)*ncol(Z)+count2):((ncol(weights)-1)*ncol(Z)+count2+length(thetas[[k]])),
                      (1+(ncol(weights)-1)*ncol(Z)+count1):((ncol(weights)-1)*ncol(Z)+count1+length(thetas[[j]]))] <- t(temp_mat)
      }
      count2 <- count2 + length(thetas[[k]])
    }
    count1 <- count1 + length(thetas[[j]])
  }
  return(cov_score_mat)
}

Thetas <- function(data, dim, X_des, thetas, sigma2, response){
  
  ## computes Thetas outside loops (in other functions) to speed up computations
  ## note that Thetas is a list containing an (n x p_j) matrix for each group
  
  Thetas <- list()                         ## intialise list
  obs <- as.vector(t(data[,1:dim]))        ## extract observations and store in an (n x T) vector
  for(j in 1:(ncol(data)-dim)){
    Theta_mat <- matrix(NA, nrow = nrow(data), ncol = length(thetas[[j]]))
    X_gp <- X_des[,1:length(thetas[[j]])]  ## adjust design matrix for each group theta (as polynomial degree may differ)
    if(length(thetas[[j]])==1){
      X_gp <- as.matrix(X_gp)
    }
    if(response == "gaussian"){
      theta_const <- X_gp%*%thetas[[j]]        ## store outside loop to speed up computations
      for(p in 1:length(thetas[[j]])){
        Theta_mat[,p] <- colSums(matrix(1/sigma2*X_gp[,p]*(obs-theta_const), nrow = dim, ncol = nrow(data)))
      }
    } else if(response == "poisson"){
      theta_const <- exp(X_gp%*%thetas[[j]])   ## store outside loop to speed up computations
      for(p in 1:length(thetas[[j]])){
        Theta_mat[,p] <- colSums(matrix(X_gp[,p]*(obs-theta_const), nrow = dim, ncol = nrow(data)))
      }
    } else if(response == "bernoulli"){
      theta_const <- exp(X_gp%*%thetas[[j]])/(1+exp(X_gp%*%thetas[[j]]))   ## store outside loop to speed up computations
      for(p in 1:length(thetas[[j]])){
        Theta_mat[,p] <- colSums(matrix(X_gp[,p]*(obs-theta_const), nrow = dim, ncol = nrow(data)))
      }
    }
    Thetas[[j]] <- Theta_mat
  }
  return(Thetas)
}

cov_score_pi <- function(pis, weights){ 
  
  ## computes covariance of the score vector with respect to pi
  
  mat <- matrix(NA, nrow = (ncol(weights)-1), ncol = (ncol(weights)-1))
  for(j in 1:ncol(mat)){
    for(k in j:ncol(mat)){
      if(j == k){
        mat[j,j] <- sum(weights[,j]*(1-weights[,j])/pis[j]^2+weights[,ncol(weights)]*(1-weights[,ncol(weights)])
                        /pis[length(pis)]^2+2*weights[,j]*weights[,ncol(weights)]/(pis[j]*pis[length(pis)]))
      } else {
        mat[j,k] <- sum(-weights[,j]*weights[,k]/(pis[j]*pis[k])+weights[,ncol(weights)]*(1-weights[,ncol(weights)])
                        /pis[length(pis)]^2+weights[,j]*weights[,ncol(weights)]/(pis[j]*pis[length(pis)])
                        +weights[,k]*weights[,ncol(weights)]/(pis[k]*pis[length(pis)]))
        mat[k,j] <- mat[j,k]
      }
    }
  }
  return(mat)
}

cov_score_pi_theta <- function(Theta, pis, gp, weights){ 
  
  ## computes covariance of the score vector with respect to pi and each group theta
  
  mat <- matrix(NA, nrow = (ncol(weights)-1), ncol = ncol(Theta))
  if(gp == ncol(weights)){
    for(k in 1:nrow(mat)){
      mat[k,] <- t(Theta)%*%(-weights[,gp]*(weights[,k]/pis[k]+(1-weights[,ncol(weights)])/pis[length(pis)]))
    }
  } else {
    for(k in 1:nrow(mat)){
      if(k == gp){
        mat[k,] <- t(Theta)%*%(weights[,gp]*((1-weights[,k])/pis[k]+weights[,ncol(weights)]/pis[length(pis)]))
      } else {
        mat[k,] <- t(Theta)%*%(weights[,gp]*(-weights[,k]/pis[k]+weights[,ncol(weights)]/pis[length(pis)]))
      }
    }
  }
  return(mat)
}

cov_score_gamma <- function(Z, weights){ 
  
  ## computes covariance of the score vector with respect to gamma
  ## computes each smaller matrix temp_mat and stores in block matrix mat (note that temp_mat is always symmetric)
  
  mat <- matrix(NA, nrow = ((ncol(weights)-1)*ncol(Z)), ncol = ((ncol(weights)-1)*ncol(Z)))
  vars <- weights*(1-weights)                                
  for(j in 2:ncol(weights)){         
    for(k in 2:ncol(weights)){       
      if(j == k){
        temp_mat <- matrix(NA, nrow = ncol(Z), ncol = ncol(Z))
        for(m in 1:nrow(temp_mat)){
          temp_mat[m,] <- colSums(Z[,m]*Z*vars[,j])
        }
        mat[(1+(j-2)*ncol(Z)):((j-1)*ncol(Z)),(1+(j-2)*ncol(Z)):((j-1)*ncol(Z))] <- temp_mat
      } else if(k > j){
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

cov_score_gamma_theta <- function(Theta, weights, Z, gp){ 
  
  ## computes covariance of the score vector with respect to gamma and each group theta
  ## computes each smaller matrix temp_mat and stores in block matrix mat
  
  mat <- matrix(NA, nrow = ((ncol(weights)-1)*ncol(Z)), ncol = ncol(Theta))
  for(j in 2:ncol(weights)){
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
  
  mat <- matrix(NA, nrow = ncol(Theta), ncol = ncol(Theta))
  for(m in 1:nrow(mat)){
    mat[m,] <- colSums(Theta[,m]*Theta*weights[,gp]*(1-weights[,gp]))
  }
  return(mat)
}  

cov_score_theta_theta <- function(Theta_g1, Theta_g2, weights, gp1, gp2){ 
  
  ## computes covariance of the score vector with respect two different groups
  
  mat <- matrix(NA, nrow = ncol(Theta_g1), ncol = ncol(Theta_g2))
  for(m in 1:nrow(mat)){
    mat[m,] <- -colSums(Theta_g1[,m]*Theta_g2*weights[,gp1]*weights[,gp2])
  }
  return(mat)
} 

