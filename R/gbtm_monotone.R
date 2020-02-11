requireNamespace("gtools", quietly = TRUE)
requireNamespace("Matrix", quietly = TRUE)
requireNamespace("matrixStats", quietly = TRUE)
requireNamespace("MonoPoly", quietly = TRUE)

## functions required to estimate group based trajectory models with monotonicity constraints (with standard errors script)
## to estimate a model, run gbtm_monotone noting the following specifications:

## data - a (n x T) matrix containing the observations for n individuals over T periods of time
## n_gps - the number of groups to fit
## x - if the time-varying covariate, i.e. time or age, is the same for all individuals then the user can input a vector
## containing these times, e.g. x = c(10,11,12) if all individuals have measurements at ages 10, 11 and 12
## otherwise the user is required to input a (n x T) design matrix containing all times for all individuals
## poly_degs - the degrees of the polynomials to fit to each group, e.g. c(2,2,3) will fit two quadratic groups
## and one cubic group (note the number of polynomial degrees must match the number of groups)
## n_starts - number of sets of initial values to generate (would recommend 5-10 sets unless the data is extremely complex)
## monotone - is TRUE if fitting monotone trajectories, FALSE otherwise
## covariates - if using time-stable covariates input a (n x (q+1)) matrix of q covariates
## plot_it - if TRUE will plot the estimated trajectories
## conf_level - computes confidence bands (at level conf_level) and plots with trajectories (note plot_it must also be TRUE)
## boot - LEAVE as 'off', it will turn to 'on' if computing bootstrapped iterations for confidence bands
## response - data type of observations (can be gaussian, poisson or bernoulli for gaussian, count and binary data respectively)

## for examples see gbtm_monotone_examples script (after running all functions)

gbtm_monotone <- function(data, n_gps, x, poly_degs = rep(3, n_gps), n_starts = 10, monotone = TRUE, covariates = NULL, plot_it = TRUE, conf_level = NULL, boot = 'off', response = 'gaussian'){

  ## main function to estimate group based trajectory models with monotonicity constraints
  ## first the format of the function's inputs are checked and corrected if required then
  ## several sets of initial values are generated (see paper for reference)
  ## the EM algorithm is run followed by computation of standard errors and confidence bands if required

  start.time <- Sys.time()
  if(boot == 'off'){                                                 ## if we are not in a bootstrapped run (see conf_boot)
    data <- as.matrix(data)                                          ## data required to be an n x T matrix
    dim <- ncol(data)                                                ## each individual has a T dimensional data sequence
    data <- cbind(data, matrix(0, nrow = nrow(data), ncol = n_gps))  ## columns added to store weights for each group
    if(is.vector(x)){                                                ## if time-varying covariate is the same for all individuals
      x_scl <- (x - mean(x))/sd(x)                                   ## scale the time varying covariate for computation
      X_des <- matrix(t(matrix(x_scl, nrow = dim,                    ## create design matrix X_des
                               ncol = (max(poly_degs)+1)))^(0:max(poly_degs)), nrow = dim*nrow(data), ncol = (max(poly_degs)+1), byrow = T)
    } else if(is.matrix(x)){                                         ## if time-varying covariate is not the same for all individuals
      X_des_unscl <- as.vector(t(x))                                 ## start creating design matrix X_des (unscaled)
      x_unique <- unique(X_des_unscl)
      x_scl <- (x_unique - mean(x_unique))/sd(x_unique)              ## scale the unique times
      X_des <- vector(length = length(X_des_unscl))                  ## set up scaled design matrix
      for(i in 1:length(x_scl)){
        X_des[X_des_unscl == sort(x_unique)[i]] <- x_scl[i]          ## set each element to appropriate scaled time
      }
      X_des <- t(t(matrix(X_des, nrow = length(X_des), ncol = (max(poly_degs)+1)))^(0:max(poly_degs)))
    }
    if(!is.null(covariates)){
      covariates <- as.matrix(covariates)                            ## covariates required to be an n x q+1 matrix
      if(sum(covariates[,1] == 1) != nrow(covariates)){
        covariates <- cbind(1, covariates)                           ## ensure the first column of covariates is a vector
      }                                                              ## of ones (for the intercept term)
    }
  } else if(boot == 'on'){
    X_des <- x
    dim <- ncol(data) - n_gps
  }
  if(response == "gaussian"){
    sigma2 <-  mean(diag(var(data[,1:dim])))                       ## if required estimate the initial value of sigma2
  }
  else {sigma2 <- NULL}
  mus <- list()                                                    ## create lists to store initial values
  pis <- list()
  poly_degs_temp <- matrix(NA, nrow = n_starts, ncol = n_gps)
  log_liks <- vector(length = n_starts)                                                                    ## list all permutations of the polynomial
  poly_degs_orders <- unique(gtools::permutations(n_gps, n_gps, poly_degs, set = FALSE, repeats.allowed = FALSE))  ## degrees to test during initial value search
  for(i in 1:n_starts){                                           ## generate initial values using kmeans
    init_values <- kmeans(data[,1:dim], n_gps)                    ## and store estimates do this n_starts times and proceed
    mus[[i]] <- init_values$centers
    pis[[i]] <- matrix(rep(init_values$size, nrow(data)), nrow = nrow(data), byrow = TRUE)/nrow(data)
    init_ests <- list(pis[[i]], mus[[i]], sigma2)
    log_liks_orders <- vector(length = nrow(poly_degs_orders))    ## for each permutation of the polynomial degrees run the
    data_temp <- e_step_init(data, dim, init_ests, response)      ## EM algorithm one iteration and record the log-likelihood
    for(j in 1:nrow(poly_degs_orders)){
      estimates_temp <- m_step(data_temp, dim, X_des, poly_degs_orders[j,], monotone, covariates, response)
      log_liks_orders[j] <- log_lik(data_temp, dim, estimates_temp, X_des, poly_degs_orders[j,], response)
    }
    poly_degs_temp[i,] <- poly_degs_orders[which.max(log_liks_orders),]   ## for each set of initial values store the permutation
    log_liks[i] <- max(log_liks_orders)                                   ## of the polynomial degrees that maximising the log-likelihood
    cat("Initial Values Generated",i,"\n")
  }
  poly_degs <- poly_degs_temp[which.max(log_liks),]                                  ## choose the set of initial values with corresponding
  init_ests <- list(pis[[which.max(log_liks)]], mus[[which.max(log_liks)]], sigma2)  ## degree permutation that maximises the log-likelihood
  data <- e_step_init(data, dim, init_ests, response)                                ## proceed once more to the first EM iteration
  estimates <- m_step(data, dim, X_des, poly_degs, monotone, covariates, response)
  log_lik <- log_lik(data, dim, estimates, X_des, poly_degs, response)
  cat("Iteration: ","0   ","Log-Likelihood: ",log_lik, "\n")              ## log-likelihood at optimal initial values
  log_lik_old <- -Inf
  count <- 0
  while(log_lik > log_lik_old + 0.0001 & log_lik > 1.0001*log_lik_old & count < 10000){    ## criteria for convergence
    count <- count + 1
    log_lik_old <- log_lik
    data <- e_step(data, dim, estimates, X_des, poly_degs, response)                  ## expectation step
    estimates <- m_step(data, dim, X_des, poly_degs, monotone, covariates, response)  ## maximisation step
    log_lik <- log_lik(data, dim, estimates, X_des, poly_degs, response)              ## compute log-likelihood
    cat("Iteration: ",as.character(count),"  ","Log-Likelihood: ",log_lik, "\n")      ## log-likelihood at current iteration
  }
  thetas <- estimates[[2]]                                                            ## store polynomial coefficients
  if(boot == 'on') return(thetas)                                                               ## if in a bootstrapped run end here and return thetas
  cat("Converged","\n")                                                                         ## declare convergence
  standard_errors <- sqrt(diag(std_errors(data, dim, estimates, X_des, covariates, response)))  ## compute standard errors
  if(plot_it == TRUE){
    cols <- c("blue","red","green","black","purple","yellow","orange","pink","darkgreen","brown","gray")        ## line colours for groups
    x_pred <- seq(from = min(X_des[,2]), to  = max(X_des[,2]), by = 0.01)                                       ## range to predict over
    X_pred <- t(matrix(x_pred, nrow = (max(poly_degs)+1), ncol = length(x_pred), byrow = T)^(0:max(poly_degs))) ## design matrix for prediction
    y_pred <- matrix(NA, nrow = length(x_pred), ncol = n_gps)
    for(j in 1:n_gps){                                                                                          ## compute fitted values
      if(response == 'gaussian'){
        y_pred[,j] <- X_pred[,1:length(thetas[[j]])]%*%thetas[[j]]
      } else if(response == 'poisson'){
        y_pred[,j] <- exp(X_pred[,1:length(thetas[[j]])]%*%thetas[[j]])
      } else if(response == 'bernoulli'){
        y_pred[,j] <- exp(X_pred[,1:length(thetas[[j]])]%*%thetas[[j]])/(1+exp(X_pred[,1:length(thetas[[j]])]%*%thetas[[j]]))
      }
    }
    if(!is.null(conf_level)){
      cat("Computing confidence intervals ...", "\n")                      ## compute confidence bands at level conf_level
      conf_bands <- conf_boot(data, estimates, 250, conf_level, n_gps, dim, X_des, n_starts, X_pred, poly_degs, covariates, monotone, response)
      plot(y_pred[,1]~x_pred, type = 'l', ylim = c((min(y_pred)-0.1), (max(y_pred)+0.1)), col = cols[1], xlab = "Scaled Time", ylab = "Response")     ## plot fitted values
      conf_cols <- c("gray50","gray55","gray60","gray65","gray70","gray75")                                  ## confidence band colours
      for(j in 1:n_gps){                                                                                     ## add confidence bands
        polygon(c(x_pred, rev(x_pred)), c(conf_bands[[2]][j,], rev(conf_bands[[1]][j,])), col = conf_cols[j], border = NA)
      }
    } else {
      plot(y_pred[,1]~x_pred, type = 'l', ylim = c((min(y_pred)-0.1), (max(y_pred)+0.1)), col = cols[1], xlab = "Scaled Time", ylab = "Response")     ## plot fitted values
    }
    for(j in 1:n_gps){
      lines(y_pred[,j]~x_pred, col = cols[j])
    }
  }
  out <- data.frame(Group = factor(rep('',length(unlist(thetas))), levels = c('',1:n_gps)))
  param_names <- c("Intercept","Linear","Quadratic","Cubic","Quartic","Quintic")
  out$Parameter <- 0
  count <- 1
  for(j in 1:n_gps){
    out$Group[count] <- j
    out$Parameter[count:(count+length(thetas[[j]])-1)] <- param_names[1:length(thetas[[j]])]
    count <- count + length(thetas[[j]])
  }
  out$Coefficient <- round(unlist(thetas),4)
  if(is.null(covariates)){
    out$Standard_Error <- round(standard_errors[-c(1:(n_gps-1))],4)
  } else {
    out$Standard_Error <- round(standard_errors[-c(1:(ncol(covariates)*(n_gps-1)))],4)
  }
  out$t_value <- round(unlist(thetas)/out$Standard_Error,4)
  out$p_value <- round(2*pt(abs(out$t_value), df = length(unlist(thetas)), lower.tail = FALSE),4)
  if(is.null(covariates)){
    out_covs <- data.frame(Group = c(1:n_gps), Membership = estimates[[1]])
  } else {
    out_covs <- data.frame(Group = factor(rep('',(1+ncol(covariates)*(n_gps-1))), levels = c('',1:n_gps)))
    covs_names <- vector(length = (ncol(covariates)-1))
    for(i in 1:length(covs_names)){
      covs_names[i] <- paste("Covariate",i)
    }
    count <- 1
    for(j in 1:n_gps){
      out_covs$Group[count] <- j
      if(j == 1){
        out_covs$Parameter[count] <- 'Constant'
        count <- count + 1
      } else {
        out_covs$Parameter[count:(count + ncol(covariates) - 1)] <- c('Constant',covs_names)
        count <- count + ncol(covariates)
      }
    }
    out_covs$Coefficient <- round(c(0,t(estimates[[4]][-1,])),4)
    out_covs$Standard_Error <- c('',round(standard_errors[1:(ncol(covariates)*(n_gps-1))],4))
    out_covs$t_value <- c('',round(out_covs$Coefficient[-1]/as.numeric(out_covs$Standard_Error[-1]),4))
    out_covs$p_value <- c('',round(2*pt(abs(as.numeric(out_covs$t_value[-1])), df = length(unlist(thetas)), lower.tail = FALSE),4))
  }
  memberships <- data.frame(matrix(NA, nrow = nrow(data), ncol = 0))
  pis <- data.frame(matrix(NA, nrow = ifelse(is.matrix(estimates[[1]]), nrow(estimates[[1]]), 1) , ncol = 0))
  for(i in 1:n_gps){
    memberships[,paste("Group",i)] <- round(data[, dim+i],10)
    pis[,paste("Group",i)] <- ifelse(is.matrix(estimates[[1]]), estimates[[1]][,i], estimates[[1]][i])
    thetas[[i]] <- as.vector(thetas[[i]])
  }
  if(!is.null(covariates)){
    gammas <- data.frame(Intercept = estimates[[4]][,1])
    for(i in 2:ncol(covariates)) gammas[, paste("Covariate", i)] <- estimates[[4]][,i]
  } else { gammas <- NULL }
  end.time <- Sys.time()
  time.diff <- end.time - start.time
  return(list(memberships = memberships, pis = pis, thetas = thetas, sigma2 = estimates[[3]], gammas = gammas, log_likelihood = log_lik, time = time.diff, summary = out,  covariates_summary = out_covs))
}

e_step_init <- function(data, dim, estimates, response){

  ## expectation step (used strictly for initial values)
  ## computes the group membership probabilities (weights)
  ## for individual i in group j given estimated group means from kmeans

  weights <- matrix(NA, nrow = nrow(data), ncol = (ncol(data) - dim))  ## to store group weights
  obs <- data[,1:dim]                                                  ## extract observations and store in an (n x T) vector
  for(j in 1:ncol(weights)){
    if(response == "gaussian"){
      weights[,j] <- matrixStats::rowProds(dnorm(obs, mean = matrix(estimates[[2]][j,], nrow = nrow(data), ncol = dim, byrow = T), sd = sqrt(estimates[[3]])))
    }
    else if(response == "poisson"){
      weights[,j] <- matrixStats::rowProds(dpois(obs, lambda = matrix(estimates[[2]][j,], nrow = nrow(data), ncol = dim, byrow = T)))
    }
    else if(response == "bernoulli"){
      weights[,j] <- matrixStats::rowProds(dbinom(obs, size = 1, prob = matrix(estimates[[2]][j,], nrow = nrow(data), ncol = dim, byrow = T)))
    }
  }
  weights <- estimates[[1]]*weights      ## multiply elementwise by group proportions
  weights <- weights/rowSums(weights)    ## standardise weights (individual rows sum to one)
  return(cbind(data[,1:dim], weights))   ## return data containing observations and weights
}

e_step <- function(data, dim, estimates, X_des, poly_degs, response){

  ## expectation step (used during standard EM steps) computes the
  ## group membership probabilities (weights) for individual i in group j
  ## given parameter estimates from the maximisation step

  weights <- matrix(NA, nrow = nrow(data), ncol = (ncol(data) - dim))     ## to store group weights
  obs <- data[,1:dim]                                                     ## extract observations and store in an (n x T) vector
  for(j in 1:ncol(weights)){
    if(response == "gaussian"){
      weights[,j] <- matrixStats::rowProds(dnorm(obs, mean = matrix(X_des[,1:(poly_degs[j]+1)]%*%estimates[[2]][[j]], nrow = nrow(data), byrow = T), sd = sqrt(estimates[[3]])))
    }
    else if(response == "poisson"){
      weights[,j] <- matrixStats::rowProds(dpois(obs, lambda = exp(matrix(X_des[,1:(poly_degs[j]+1)]%*%estimates[[2]][[j]], nrow = nrow(data), byrow = T))))
    }
    else if(response == "bernoulli"){
      probs_trans <- matrix(X_des[,1:(poly_degs[j]+1)]%*%estimates[[2]][[j]], nrow = nrow(data), byrow = T)
      weights[,j] <- matrixStats::rowProds(dbinom(obs, size = 1, prob = exp(probs_trans)/(1+exp(probs_trans))))
    }
  }
  if(is.matrix(estimates[[1]])){            ## multiply elementwise by individual group weightings
    weights <- estimates[[1]]*weights       ## method depends on if covariates are provided or not
  } else {                                  ## see computation of pis in m_step
    weights <- t(estimates[[1]]*t(weights))
  }
  weights <- weights/rowSums(weights)       ## standardise weights (individual rows sum to one)
  return(cbind(data[,1:dim], weights))      ## return data containing observations and weights
}

m_step <- function(data, dim, X_des, poly_degs, monotone, covariates, response){

  ## maximisation step to maximise the log-likelihood given group memberships
  ## from the expectation step. Computes parameter estimates for pis (and/or gammas)
  ## theta and sigma2 (for gaussian response only)

  n_gps <- ncol(data) - dim
  if(is.null(covariates)){
    n_effs <- colSums(data[,-c(1:dim)])     ## effective number of individuals per group
    pis <- n_effs/nrow(data)                ## estimated group membership proportions
    gammas <- NULL
  } else {                                  ## estimate gammas with group proportions function group_props
    opt <- optim(par = rep(0,n_gps*ncol(covariates)), fn = group_props, weights = data[,-c(1:dim)], n_gps = n_gps, covariates = covariates,
                 control = list(maxit = 1000, REPORT = T, fnscale = -1), lower = rep(c(-1e-16,rep(-Inf,(n_gps-1))),ncol(covariates)),
                 upper = rep(c(1e-16,rep(Inf,(n_gps-1))),ncol(covariates)), method = 'L-BFGS-B')
    gammas <- matrix(opt$par, nrow = n_gps)
    pis <- t(exp(gammas%*%t(covariates)))/rowSums(t(exp(gammas%*%t(covariates))))
  }
  thetas <- list()                       ## initialise list for ML estimates
  obs <- as.vector(t(data[,1:dim]))      ## extract observations and store in an (n x T) vector
  if(response == "gaussian"){
    cumss <- 0                           ## start cumulative sums of squares to estimate variance
    for(j in 1:n_gps){
      weights_ext <- rep(data[, (dim+j)], each = dim)   ## extend weights to match length of design matrix (n x T)
      if(monotone == TRUE){
        if(poly_degs[j]%%2 == 0){      ## if fitting quadratics or quartics set appropriate lower bound
          deg.is.odd <- FALSE
          type <- 1
          K <- trunc((poly_degs[j]-1)/2) + 1
          lower <- min(X_des[,2])      ## then use monpol to estimate theta under monotonicity constraints
        } else {deg.is.odd <- TRUE; type <- 0; K <- (poly_degs[j]-1)/2;
          lower <- -Inf}
        model <- MonoPoly::SOSpol.fit(x = X_des[,2], y = obs, w = weights_ext, deg.is.odd = deg.is.odd,
                            a = lower, b = Inf, K = K, type = type)
        #model <- MonoPoly::monpol(obs ~ X_des[,2], weights = weights_ext, degree = poly_degs[j], a = lower, b = Inf)
        #thetas[[j]] <- as.matrix(coef(model))
        thetas[[j]] <- as.matrix(model$beta.raw)
        #cumss <- cumss + sum(residuals(model)^2*weights_ext)
        cumss <- cumss + sum(weights_ext*(obs-X_des[,1:(poly_degs[j]+1)]%*%thetas[[j]])^2)   ## compute sums of squares
      } else {
        model <- lm.wfit(x = as.matrix(X_des[,1:(poly_degs[j]+1)]), y = obs, w = weights_ext) ## if optimisation is unconstrained
        thetas[[j]] <- as.matrix(coef(model))                                                 ## use lm instead of monpol
        cumss <- cumss + sum(resid(model)^2*weights_ext)
      }
    }
    df <- length(obs) - sum(poly_degs+1)  ## compute degrees of freedom
    sigma2 <- cumss/df                    ## estimate sigma2
  }
  ## if response is poisson or bernoulli use glm instead of lm and disregard sigma2
  else if(response == "poisson"){
    for(j in 1:n_gps){
      weights_ext <- rep(data[, (dim+j)], each = dim)
      model <- glm.fit(x = as.matrix(X_des[,1:(poly_degs[j]+1)]), y = obs, weights = weights_ext, family = poisson(link = "log"))
      thetas[[j]] <- as.matrix(coef(model))
    }
    sigma2 <- NULL
  }
  else if(response == "bernoulli"){
    for(j in 1:n_gps){
      weights_ext <- rep(data[, (dim+j)], each = dim)
      model <- glm.fit(x = as.matrix(X_des[,1:(poly_degs[j]+1)]), y = obs, weights = weights_ext, family = quasibinomial(link = "logit"))
      thetas[[j]] <- as.matrix(coef(model))
    }
    sigma2 <- NULL
  }
  return(list(pis, thetas, sigma2, gammas))
}

group_props <- function(gammas, weights, n_gps, covariates){

  ## function to optimise for group proportion estimation
  ## note that covariates is an n x q+1 matrix and
  ## gammas is a J x q+1 matrix of initial values for gammas

  gammas <- matrix(gammas, nrow = n_gps, ncol = ncol(covariates))
  opt <- sum(t(weights)*(gammas%*%t(covariates) - matrix(rep(log(colSums(exp(gammas%*%t(covariates)))), each = n_gps), nrow = n_gps)))
  return(opt)
}

log_lik <- function(data, dim, estimates, X_des, poly_degs, response){

  ## computes the log-likelihood given the parameter estimates from maximisation step

  n_gps <- ncol(data) - dim
  obs <- data[,1:dim]                                   ## extract observations
  mat <- matrix(NA, nrow = nrow(data), ncol = n_gps)    ## to store components of log-likelihood
  if(response == "gaussian"){
    for(j in 1:n_gps){
      mat[,j] <- matrixStats::rowProds(dnorm(obs, mean = matrix(X_des[,1:(poly_degs[j]+1)]%*%estimates[[2]][[j]], nrow = nrow(data), byrow = T), sd = sqrt(estimates[[3]])))
    }
  } else if(response == "poisson"){
    for(j in 1:n_gps){
      mat[,j] <- matrixStats::rowProds(dpois(obs, lambda = exp(matrix(X_des[,1:(poly_degs[j]+1)]%*%estimates[[2]][[j]], nrow = nrow(data), byrow = T))))
    }
  } else if(response == "bernoulli"){
    for(j in 1:n_gps){
      probs_trans <- matrix(X_des[,1:(poly_degs[j]+1)]%*%estimates[[2]][[j]], nrow = nrow(data), byrow = T)
      mat[,j] <- matrixStats::rowProds(dbinom(obs, size = 1, prob = exp(probs_trans)/(1+exp(probs_trans))))
    }
  }
  if(is.matrix(estimates[[1]])){       ## multiply elementwise by individual group weightings
    mat <- estimates[[1]]*mat          ## method depends on if covariates are provided or not
  } else {                             ## see computation of pis in m_step
    mat <- t(estimates[[1]]*t(mat))
  }
  return(sum(log(rowSums(mat))))
}

conf_boot <- function(data, estimates, nboot, conf_level, n_gps, dim, X_des, n_starts, X_pred, poly_degs, covariates, monotone, response){

  ## computes (conf_level) confidence bands using a cases bootstrap
  ## the entire algorithm is run nboot times on resampled data sets

  true_coefs <- list()
  for(j in 1:n_gps){
    true_coefs[[j]] <- matrix(0, nrow = (max(poly_degs)+1), ncol = 1)     ## list of current estimated coefficients each with
    true_coefs[[j]][1:(poly_degs[j]+1)] <- estimates[[2]][[j]]            ## the same length (higher order terms set to zero)
  }
  coefs <- list()
  for(j in 1:n_gps){
    coefs[[j]] <- matrix(NA, nrow = nboot, ncol = (max(poly_degs)+1))     ## matrices to store bootstrapped coefficients
  }
  if(is.null(covariates)){
    covariates_temp <- NULL
  }
  print_k <- 0
  for(i in 1:nboot){
    ## sample new data set to bootstrap, run algorithm and extract estimated coefficients
    boot_inds <- sample(1:nrow(data), size = nrow(data), replace = T)
    boot_dat <- cbind(data[boot_inds, 1:dim], matrix(0, nrow = nrow(data), ncol = n_gps))
    if(length(unique(X_des[,2])) == dim){     ## if all individuals have same time_varying covariates
      X_des_temp <- X_des                     ## then the design matrix can stay unchanged
    } else {                                  ## otherwise reformat design matrix to match bootstrapped individuals
      X_des_temp <- matrix(X_des[,2], nrow = nrow(data), ncol = dim, byrow = T)[boot_inds,]
      X_des_temp <- t(t(matrix(as.vector(t(X_des_temp)), nrow = nrow(data)*dim, ncol = (max(poly_degs)+1)))^(0:max(poly_degs)))
    }
    if(!is.null(covariates)){
      covariates_temp <- covariates[boot_inds,]
    }
    invisible(capture.output(est_coefs <- gbtm_monotone(boot_dat, n_gps, x = X_des_temp, poly_degs = poly_degs, n_starts = n_starts, monotone = monotone, covariates = covariates_temp, plot_it = FALSE, conf_level = NULL, boot = 'on', response = response)))
    ## to avoid label switching problem we first set estimated coefficients to have the same length
    ## by setting the higher order terms to zero
    ## we then reorder the groups to match the order of true_coefs using minimum euclidean distance between estimates
    ## note that each set of estimates may only be chosen once

    coefs_temp <- list()
    for(j in 1:n_gps){
      coefs_temp[[j]] <- matrix(0, nrow = (max(poly_degs)+1), ncol = 1)
      coefs_temp[[j]][1:length(est_coefs[[j]])] <- est_coefs[[j]]
    }
    dists_old <- NULL
    for(j in 1:n_gps){
      dists <- sapply(coefs_temp, function(x) sqrt(sum((x-true_coefs[[j]])^2)))
      dists[dists_old] <- Inf
      coefs[[j]][i,] <- coefs_temp[[which.min(dists)]]
      dists_old <- cbind(dists_old, which.min(dists))
    }
    if(i >= nboot*1/5 & print_k == 0){
      cat("20%", "\n")
      print_k <- 1
    }
    if(i >= nboot*2/5 & print_k == 1){
      cat("40%", "\n")
      print_k <- 2
    }
    if(i >= nboot*3/5 & print_k == 2){
      cat("60%", "\n")
      print_k <- 3
    }
    if(i >= nboot*4/5 & print_k == 3){
      cat("80%", "\n")
      print_k <- 4
    }
    if(i == nboot) cat("100%")
  }
  lower <- matrix(NA, nrow = n_gps, ncol = nrow(X_pred))
  upper <- matrix(NA, nrow = n_gps, ncol = nrow(X_pred))
  for(j in 1:n_gps){
    if(response == 'gaussian'){                    ## compute fitted values using each group's
     fit_vals <- X_pred%*%t(coefs[[j]])            ## bootstrapped coefficients and take appropriate quantiles
    } else if(response == 'poisson'){
     fit_vals <- exp(X_pred%*%t(coefs[[j]]))
    } else if(response == 'bernoulli'){
     fit_vals <- exp(X_pred%*%t(coefs[[j]]))/(1+exp(X_pred%*%t(coefs[[j]])))
    }
    lower[j,] <- apply(fit_vals, 1, function(x) quantile(x, probs = (1-conf_level/100)/2))
    upper[j,] <- apply(fit_vals, 1, function(x) quantile(x, probs = (1-(1-conf_level/100)/2)))
  }
  return(list(lower, upper))
}

