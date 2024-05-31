requireNamespace("gtools", quietly = TRUE)
requireNamespace("Matrix", quietly = TRUE)
requireNamespace("matrixStats", quietly = TRUE)
requireNamespace("MonoPoly", quietly = TRUE)
requireNamespace("ggplot2", quietly = TRUE)

## functions required to estimate group based trajectory models with monotonicity constraints (with standard errors script)
## to estimate a model, run gbtm_monotone noting the following specifications:

## data - a (n x T) matrix containing the observations for n individuals over T periods of time
## n_gps - the number of groups to fit
## x - if the time-varying covariate, i.e. time or age, is the same for all individuals then the user can input a vector
## containing these times, e.g. x = c(10,11,12) if all individuals have measurements at ages 10, 11 and 12
## otherwise the user is required to input a (n x T) design matrix containing all times for all individuals
## poly_degs - the degrees of the polynomials to fit to each group, e.g. c(2,2,3) will fit two quadratic groups
## and one cubic group (note the number of polynomial degrees must match the number of groups)
## n_starts - number of sets of initial values to generate (5 sets is recommended unless the dataset is extremely complex)
## monotone - is TRUE if fitting monotone trajectories, FALSE otherwise
## covariates - if using time-stable covariates input a (n x (q+1)) matrix of q covariates
## plot_it - if TRUE will plot the estimated trajectories
## conf_level - computes confidence bands (at level conf_level) and plots with trajectories (note plot_it must also be TRUE)
## boot - LEAVE as 'off', it will turn to 'on' if computing bootstrapped iterations for confidence bands
## sims - LEAVE as 'off' (required for testing)
## response - data type of observations (can be gaussian, poisson or bernoulli for gaussian, count and binary data respectively)

## for examples see gbtm_monotone_examples script (after running all functions)

gbtm_monotone <- function(data, n_gps, x, poly_degs = rep(3, n_gps), n_starts = 5, monotone = TRUE, covariates = NULL, plot_it = TRUE, conf_level = NULL, boot = 'off', sims = 'off', response = 'gaussian'){
  
  ## main function to estimate group based trajectory models with monotonicity constraints
  ## first the format of the function's inputs are checked and corrected if required then
  ## several sets of initial values are generated 
  ## the EM algorithm is run followed by computation of standard errors and confidence bands if required
  
  start.time <- Sys.time()
  if(boot == 'off'){                                                   ## if we are not in a bootstrapped run (see conf_boot)
    data <- as.matrix(data)                                            ## data required to be an n x TT matrix
    if(any(is.na(data))) stop("Failed: Missing Data Detected")         ## stop if there is any missing data
    n <- nrow(data)                                                    ## number of individuals
    TT <- ncol(data)                                                   ## number of time points
    weights <- matrix(0, nrow = n, ncol = n_gps)                       ## store weights for each group
    p <- max(poly_degs)                                                ## maximum polynomial degree
    if(is.vector(x)){                                                  ## if time-varying covariate is the same for all individuals
      if(length(x) != TT) stop("Failed: Dimension of data not equal to dimension of time varying covariate")
      x_scl <- (x - mean(x))/sd(x)                                     ## scale the time varying covariate for computation
      X_des <- matrix(t(matrix(x_scl, nrow = TT, ncol = p+1))^(0:p),   ## create design matrix X_des        
                      nrow = TT*n, 
                      ncol = p+1,
                      byrow = TRUE)
    } else if(is.matrix(x)){                                           ## if time-varying covariate is not the same for all individuals
      if(ncol(x) != TT) stop("Failed: Dimension of data not equal to dimension of time varying covariate")
      X_des_unscl <- as.vector(t(x))                                   ## start creating design matrix X_des (unscaled)
      x_unique <- unique(X_des_unscl)
      x_scl <- (x_unique - mean(x_unique))/sd(x_unique)                ## scale the unique times
      X_des <- vector(length = length(X_des_unscl))                    ## set up scaled design matrix
      for(i in 1:length(x_scl)){
        X_des[X_des_unscl == sort(x_unique)[i]] <- x_scl[i]            ## set each element to appropriate scaled time
      }
      X_des <- t(t(matrix(X_des, nrow = length(X_des), ncol = p+1))^(0:p))
    }
    if(length(poly_degs) != n_gps) stop("Failed: Number of polynomials not equal to number of groups")
    if(!is.null(covariates)){
      covariates <- as.matrix(covariates)                              ## covariates required to be an n x q+1 matrix
      if(!all(covariates[,1] == 1)) covariates <- cbind(1, covariates) ## ensure the first column of covariates is a vector of ones (for the intercept term)
    } 
    cat("Starting Generation of Initial Values", "\n")
  } else if(boot == 'on'){                                             ## store key values from bootstrapped data
    X_des <- x
    n <- nrow(data)
    TT <- ncol(data) - n_gps
  }
  if(response == "gaussian"){
    sigma2 <-  mean(diag(var(data)))                                  ## if required estimate the initial value of sigma2
  } else {
    sigma2 <- NULL
  }
  init_values <- lapply(1:n_starts, function(i) kmeans(data, n_gps))                                                   ## run kmeans n_starts times
  mus <- lapply(init_values, function(x) x$centers)                                                                    ## store centers
  pis <- lapply(init_values, function(x) matrix(rep(x$size, n), nrow = n, byrow = TRUE)/n)                             ## store sizes
  init_ests <- lapply(1:n_starts, function(i) list(pis = pis[[i]], mus = mus[[i]], sigma2 = sigma2))                   ## store initial values
  weights <- lapply(1:n_starts, function(i) e_step_init(data, n, TT, init_ests[[i]], response))                        ## run E step once for each set of initial values
  poly_degs_orders <- unique(gtools::permutations(n_gps, n_gps, poly_degs, set = FALSE, repeats.allowed = FALSE))      ## degrees to test during initial value search
  init_log_liks <- sapply(weights, function(w){                                                                        ## run M step once for each combination then store log-likelihood
                      apply(poly_degs_orders, 1, function(order){
                        ests_tmp <- m_step(data, w, n, TT, n_gps, X_des, order, monotone, covariates, response)
                        log_lik(data, n, TT, n_gps, ests_tmp, X_des, order, response)})}) 
  if(is.matrix(init_log_liks)){                                                                                        
    max_log_lik_pos_arr <- which(init_log_liks == max(init_log_liks), arr.ind = TRUE)                                      
    max_log_lik_pos <- max_log_lik_pos_arr[,"col"]                                                                     ## position of optimal values
    poly_degs <- poly_degs_orders[max_log_lik_pos_arr[,"row"],]                                                        ## degree permutation that maximises the log-likelihood 
  } else {
    max_log_lik_pos <- which.max(init_log_liks)                                                                        ## position of optimal values
    poly_degs <- as.vector(poly_degs_orders)                                                                           ## degree permutation that maximises the log-likelihood 
  }                
  init_ests <- list(pis[[max_log_lik_pos]], mus[[max_log_lik_pos]], sigma2)                                            ## optimal initial estimates
  weights <- weights[[max_log_lik_pos]]                                                                                ## weights at optimal initial estimates
  if(boot == 'off') cat("Finished Generation of Initial Values", "\n")
  estimates <- m_step(data, weights, n, TT, n_gps, X_des, poly_degs, monotone, covariates, response)                   ## proceed with initial maximisation step
  log_lik_c <- log_lik(data, n, TT, n_gps, estimates, X_des, poly_degs, response)                                      ## log-likelihood at optimal initial estimates
  if(boot == 'off') cat("Iteration: ", "0   ", "Log-Likelihood: ", log_lik_c, "\n")                      
  log_lik_old <- -Inf
  count <- 0
  while(log_lik_c > log_lik_old + 0.01 & log_lik_c > 1.001*log_lik_old & count < 10000){                               ## criteria for convergence
    count <- count + 1
    log_lik_old <- log_lik_c
    weights <- e_step(data, n, TT, n_gps, estimates, X_des, poly_degs, response)                                       ## expectation step
    estimates <- m_step(data, weights, n, TT, n_gps, X_des, poly_degs, monotone, covariates, response)                 ## maximisation step
    log_lik_c <- log_lik(data, n, TT, n_gps, estimates, X_des, poly_degs, response)                                    ## compute log-likelihood
    if(boot == 'off') cat("Iteration: ", as.character(count), "  ", "Log-Likelihood: ", log_lik_c, "\n")               ## log-likelihood at current iteration
  }
  thetas <- estimates[["thetas"]]                                                                                      ## store polynomial coefficients
  if(boot == 'on') return(thetas)                                                                                      ## if in a bootstrapped run end here and return thetas
  if(sims == "on") return(list(thetas = thetas, log_likelihood = log_lik_c, time = Sys.time() - start.time))           ## if testing return output here                   
  cat("Converged","\n")                                                                                                ## declare convergence
  standard_errors <- sqrt(diag(std_errors(data, weights, TT, n_gps, estimates, X_des, covariates, response)))          ## compute standard errors
  if(plot_it == TRUE){
    x_pred <- seq(from = min(X_des[,2]), to  = max(X_des[,2]), by = 0.01)                                              ## range to predict over
    X_pred <- t(matrix(x_pred, nrow = p+1, ncol = length(x_pred), byrow = TRUE)^(0:p))                                 ## design matrix for prediction
    if(response == 'gaussian')  y_pred <- sapply(thetas, function(gp)     X_pred[,1:length(gp)]%*%gp)                  ## compute fitted values
    if(response == 'poisson')   y_pred <- sapply(thetas, function(gp) exp(X_pred[,1:length(gp)]%*%gp))
    if(response == 'bernoulli') y_pred <- sapply(thetas, function(gp) exp(X_pred[,1:length(gp)]%*%gp)/(1+exp(X_pred[,1:length(gp)]%*%gp)))
    dat_vis <- data.frame(x = rep(x_pred, n_gps), y = as.vector(y_pred), Group = as.factor(rep(as.character(1:n_gps), each = length(x_pred))))
    g_plot <- ggplot2::ggplot(data = dat_vis, ggplot2::aes(x = x, y = y, group = Group)) + 
              ggplot2::geom_line(ggplot2::aes(colour = Group), size = 1) + 
              ggplot2::xlab("Scaled Time") + 
              ggplot2::ylab("Response") + 
              ggplot2::scale_color_manual(values = c("blue","red","green","black","purple","yellow","orange","pink","darkgreen","brown","gray")) +
              ggplot2::theme_bw()
    if(!is.null(conf_level)){
      cat("Computing confidence bands ...", "\n")                                                                                                 
      conf_bands <- conf_boot(data, n, TT, n_gps, p, thetas, conf_level, X_des, n_starts, X_pred, poly_degs, covariates, monotone, response)
      dat_vis <- cbind(dat_vis, data.frame(lower = as.vector(conf_bands[[1]]), upper = as.vector(conf_bands[[2]])))
      g_plot <- g_plot + 
                ggplot2::geom_ribbon(data = dat_vis, ggplot2::aes(ymin = lower, ymax = upper, fill = Group)) + 
                ggplot2::geom_line(ggplot2::aes(colour = Group), size = 1) + 
                ggplot2::scale_fill_manual(values = c("gray50","gray60","gray70","gray80","gray90"), drop = TRUE) 
    } 
    print(g_plot)
  }
  param_names <- c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic")
  out <- data.frame(Group = rep('',length(unlist(thetas))), 
                    Parameter = as.vector(unlist(sapply(thetas, function(x) param_names[1:length(x)]))), 
                    Coefficient = round(unlist(thetas),4))
  out[out$Parameter == "Intercept",]$Group <- 1:n_gps
  if(is.null(covariates)){
    out$Standard_Error <- round(standard_errors[-c(1:(n_gps-1))],4)
    out_covs <- data.frame(Group = c(1:n_gps), Membership = estimates[["pis"]])
  } else{
    out$Standard_Error <- round(standard_errors[-c(1:(ncol(covariates)*(n_gps-1)))],4)
    covs_names <- paste("Covariate",1:(ncol(covariates)-1))
    out_covs <- data.frame(Group = rep('',(1+ncol(covariates)*(n_gps-1))),
                           Parameter = c("Constant", rep(c("Constant", covs_names), n_gps - 1)), 
                           Coefficient = round(c(0,t(estimates[["gammas"]][-1,])),4), Standard_Error = c('',round(standard_errors[1:(ncol(covariates)*(n_gps-1))],4)))
    out_covs[out_covs$Parameter == "Constant",]$Group <- 1:n_gps
    out_covs$t_value <- c('',round(out_covs$Coefficient[-1]/as.numeric(out_covs$Standard_Error[-1]),4))
    out_covs$p_value <- c('',round(2*pt(abs(as.numeric(out_covs$t_value[-1])), df = length(unlist(thetas)), lower.tail = FALSE),4))
  }
  out$t_value <- round(unlist(thetas)/out$Standard_Error,4)
  out$p_value <- round(2*pt(abs(out$t_value), df = length(unlist(thetas)), lower.tail = FALSE),4)
  memberships <- data.frame(round(weights, 10))
  colnames(memberships) <- paste("Group", 1:n_gps)
  pis <- as.data.frame(matrix(estimates[["pis"]], ncol = n_gps))
  colnames(pis) <- paste("Group", 1:n_gps)
  if(!is.null(covariates)){
    gammas <- data.frame(estimates[["gammas"]])
    colnames(gammas) <- c("Intercept", paste("Covariate", 2:ncol(covariates)))
  } else {
    gammas <- NULL
  }
  end.time <- Sys.time()
  time.diff <- end.time - start.time
  return(list(memberships = memberships, pis = pis, thetas = thetas, sigma2 = estimates[["sigma2"]], gammas = gammas, log_likelihood = log_lik_c, time = time.diff, summary = out,  covariates_summary = out_covs))
}

e_step_init <- function(data, n, TT, estimates, response){
  
  ## expectation step (used strictly for initial values)
  ## computes the group membership probabilities (weights)
  ## for individual i in group j given estimated group means from kmeans

  if(response == "gaussian")  weights <- apply(estimates[["mus"]], 1, function(mus) matrixStats::rowProds(dnorm(data,            mean = matrix(mus, nrow = n, ncol = TT, byrow = TRUE), sd = sqrt(estimates[["sigma2"]]))))
  if(response == "poisson")   weights <- apply(estimates[["mus"]], 1, function(mus) matrixStats::rowProds(dpois(data,          lambda = matrix(mus, nrow = n, ncol = TT, byrow = TRUE))))
  if(response == "bernoulli") weights <- apply(estimates[["mus"]], 1, function(mus) matrixStats::rowProds(dbinom(data, size = 1, prob = matrix(mus, nrow = n, ncol = TT, byrow = TRUE))))
  weights <- estimates[["pis"]]*weights  ## multiply elementwise by group proportions
  weights <- weights/rowSums(weights)    ## standardise weights (individual rows sum to one)
  return(weights)                        ## return weights
}

e_step <- function(data, n, TT, n_gps, estimates, X_des, poly_degs, response){
  
  ## expectation step (used during standard EM steps) computes the
  ## group membership probabilities (weights) for individual i in group j
  ## given parameter estimates from the maximisation step
  
  if(response == "gaussian")  weights <- sapply(1:n_gps, function(gp)  matrixStats::rowProds(dnorm(data, mean =       matrix(X_des[,1:(poly_degs[gp]+1)]%*%estimates[["thetas"]][[gp]], nrow = n, byrow = TRUE), sd = sqrt(estimates[["sigma2"]]))))
  if(response == "poisson")   weights <- sapply(1:n_gps, function(gp)  matrixStats::rowProds(dpois(data, lambda = exp(matrix(X_des[,1:(poly_degs[gp]+1)]%*%estimates[["thetas"]][[gp]], nrow = n, byrow = TRUE)))))
  if(response == "bernoulli") weights <- sapply(1:n_gps, function(gp){ 
                                            probs_trans <- matrix(X_des[,1:(poly_degs[gp]+1)]%*%estimates[["thetas"]][[gp]], nrow = n, byrow = TRUE)
                                            matrixStats::rowProds(dbinom(data, size = 1, prob = exp(probs_trans)/(1+exp(probs_trans))))})
                                                
  if(is.matrix(estimates[["pis"]])){            ## multiply elementwise by individual group weightings
    weights <- estimates[["pis"]]*weights       ## method depends on if covariates are provided or not
  } else {                                      ## see computation of pis in m_step
    weights <- t(estimates[["pis"]]*t(weights))
  }
  weights <- weights/rowSums(weights)           ## standardise weights (individual rows sum to one)
  return(weights)                               ## return weights
}

m_step <- function(data, weights, n, TT, n_gps, X_des, poly_degs, monotone, covariates, response){
  
  ## maximisation step to maximise the log-likelihood given group memberships from the expectation step
  ## computes parameter estimates for pis (and/or gammas), theta and sigma2 (for gaussian response only)
  
  if(is.null(covariates)){
    pis <- colSums(weights)/n               ## estimated group membership proportions               
    gammas <- NULL
  } else {                                  ## estimate gammas with group proportions function group_props
    opt <- optim(par = rep(0,(n_gps-1)*ncol(covariates)), 
                 fn = group_props, weights = weights, n_gps = n_gps, covariates = covariates,
                 control = list(maxit = 1000, REPORT = T, fnscale = -1), method = 'L-BFGS-B')
    gammas <- rbind(0, matrix(opt$par, nrow = n_gps - 1))
    pis <- t(exp(gammas%*%t(covariates)))/rowSums(t(exp(gammas%*%t(covariates))))
  }
  if(response == "gaussian"){
    weights_ext <- lapply(1:n_gps, function(gp) rep(weights[,gp], each = TT))   ## extend weights to match length of design matrix (n x T)
    if(monotone == TRUE){
      models <- lapply(1:n_gps, function(gp){
                   if(poly_degs[gp]%%2 == 0){      
                     deg.is.odd <- FALSE
                     type <- 1
                     K <- trunc((poly_degs[gp]-1)/2) + 1
                    lower <- min(X_des[,2])      
                   } else {
                     deg.is.odd <- TRUE
                     type <- 0
                     K <- (poly_degs[gp]-1)/2
                    lower <- -Inf
                    }
                    MonoPoly::SOSpol.fit(x = X_des[,2], y = as.vector(t(data)), w = weights_ext[[gp]], deg.is.odd = deg.is.odd, a = lower, b = Inf, K = K, type = type)})
      thetas <- lapply(models, function(x) as.matrix(x$beta.raw))
      cumss <- sum(sapply(1:n_gps, function(gp) sum(weights_ext[[gp]]*(as.vector(data)-X_des[,1:(poly_degs[gp]+1)]%*%thetas[[gp]])^2))) 
    } else {
      models <- lapply(1:n_gps, function(gp) lm.wfit(x = as.matrix(X_des[,1:(poly_degs[gp]+1)]), y = as.vector(t(data)), w = weights_ext[[gp]])) 
      thetas <- lapply(models, function(x) as.matrix(coef(x)))                                            
      cumss <- sum(sapply(1:n_gps, function(gp) sum(resid(models[[gp]])^2*weights_ext[[gp]])))
    }
    sigma2 <- cumss/(n*TT - sum(poly_degs+1))     
  }
  if(response == "poisson"){
    thetas <- lapply(1:n_gps, function(gp){
                 weights_ext <- rep(weights[,gp], each = TT)
                 as.matrix(coef(glm.fit(x = as.matrix(X_des[,1:(poly_degs[gp]+1)]), y = as.vector(t(data)), weights = weights_ext, family = poisson(link = "log"))))})
    sigma2 <- NULL
  }
  if(response == "bernoulli"){
    thetas <- lapply(1:n_gps, function(gp){
                 weights_ext <- rep(weights[,gp], each = TT)
                 as.matrix(coef(glm.fit(x = as.matrix(X_des[,1:(poly_degs[gp]+1)]), y = as.vector(t(data)), weights = weights_ext, family = quasibinomial(link = "logit"))))})
    sigma2 <- NULL
  }
  return(list(pis = pis, thetas = thetas, sigma2 = sigma2, gammas = gammas))
}

group_props <- function(gammas, weights, n_gps, covariates){
  
  ## function to optimise for group proportion estimation
  ## note that covariates is an n x q+1 matrix and
  ## gammas is a J x q+1 matrix of initial values for gammas
  
  gammas <- rbind(0, matrix(gammas, nrow = n_gps - 1, ncol = ncol(covariates)))
  opt <- sum(t(weights)*(gammas%*%t(covariates) - matrix(rep(log(colSums(exp(gammas%*%t(covariates)))), each = n_gps), nrow = n_gps)))
  return(opt)
}

log_lik <- function(data, n, TT, n_gps, estimates, X_des, poly_degs, response){
  
  ## computes the log-likelihood given the parameter estimates from maximisation step
  
  if(response == "gaussian") mat <- sapply(1:n_gps, function(gp) matrixStats::rowProds(dnorm(data, mean = matrix(X_des[,1:(poly_degs[gp]+1)]%*%estimates[["thetas"]][[gp]], nrow = n, byrow = TRUE), sd = sqrt(estimates[["sigma2"]]))))
  if(response == "poisson")  mat <- sapply(1:n_gps, function(gp) matrixStats::rowProds(dpois(data, lambda = exp(matrix(X_des[,1:(poly_degs[gp]+1)]%*%estimates[["thetas"]][[gp]], nrow = n, byrow = TRUE)))))
  if(response == "bernoulli"){
    mat <- sapply(1:n_gps, function(gp){
              probs_trans <- matrix(X_des[,1:(poly_degs[gp]+1)]%*%estimates[["thetas"]][[gp]], nrow = n, byrow = TRUE)
              matrixStats::rowProds(dbinom(data, size = 1, prob = exp(probs_trans)/(1+exp(probs_trans))))})
  }
  if(is.matrix(estimates[["pis"]])){       ## multiply elementwise by individual group weightings
    mat <- estimates[["pis"]]*mat          ## method depends on if covariates are provided or not
  } else {                                 ## see computation of pis in m_step
    mat <- t(estimates[["pis"]]*t(mat))
  }
  return(sum(log(rowSums(mat))))
}

conf_boot <- function(data, n, TT, n_gps, p, thetas, conf_level, X_des, n_starts, X_pred, poly_degs, covariates, monotone, response){
  
  ## computes (conf_level) confidence bands using a cases bootstrap
  ## the entire algorithm is run nboot times on resampled data sets

  n_boot <- 500
  true_coefs <- sapply(1:n_gps, function(gp){
                  mat <- matrix(0, nrow = p+1, ncol = 1)                                                      ## list of current estimated coefficients each with
                  mat[1:(poly_degs[gp]+1)] <- estimates[["thetas"]][[gp]]                                     ## the same length (higher order terms set to zero)
                  mat}) 
  boot_inds <- lapply(1:nboot, function(i) sample(1:n, size = n, replace = TRUE))                             ## resample data
  boot_dat <- lapply(boot_inds, function(inds) list(data = data[inds,], weights = matrix(0, nrow = n, ncol = n_gps)))
  if(length(unique(X_des[,2])) == TT){
    X_des_temp <- lapply(1:nboot, function(i) X_des)                                                          ## if all individuals have same time_varying covariates then the design matrix can stay unchanged
  } else {                                                                                                    ## otherwise reformat design matrix to match bootstrapped individuals
    X_des_temp <- lapply(boot_inds, function(inds){
                     mat <- matrix(X_des[,2], nrow = n, ncol = TT, byrow = TRUE)[inds,]  
                     t(t(matrix(as.vector(t(mat)), nrow = n*TT, ncol = p+1))^(0:p))})
  }
  if(is.null(covariates)){                                                                                    ## do the same for the covariates (if they exist)
    covariates_temp <- rep(NULL, nboot)
  } else {
    covariates_temp <- lapply(boot_inds, function(inds) covariates[inds,])
  }
  coefs <- lapply(1:nboot, function(sim){
                            est_coefs <- gbtm_monotone(boot_dat[[sim]], n_gps, x = X_des_temp[[sim]], poly_degs = poly_degs, n_starts = n_starts, monotone = monotone, covariates = covariates_temp[[sim]], plot_it = FALSE, conf_level = NULL, boot = 'on', response = response)
                            ## to avoid label switching problem we first set estimated coefficients to have the same length
                            ## by setting the higher order terms to zero
                            ## we then reorder the groups to match the order of true_coefs using minimum euclidean distance between estimates
                            ## note that each set of estimates may only be chosen once
                            coefs_temp <- sapply(1:n_gps, function(gp){
                                            mat <- matrix(0, nrow = p+1, ncol = 1)
                                            mat[1:length(est_coefs[[gp]])] <- est_coefs[[gp]]
                                            mat})
                            dists_old <- NULL
                            coefs_out <- coefs_temp
                            for(j in 1:n_gps){
                              dists <- apply(coefs_temp, 2, function(x) sqrt(sum((x-true_coefs[,j])^2)))
                              dists[dists_old] <- Inf
                              dists_old <- cbind(dists_old, which.min(dists))
                              coefs_out[,j] <- coefs_temp[,which.min(dists)]
                            }
                            coefs_out})
  coefs <- lapply(1:n_gps, function(gp) do.call(cbind, lapply(coefs, function(x) x[,gp])))
  if(response == 'gaussian')  fit_vals <- lapply(coefs, function(ests)     X_pred%*%ests)
  if(response == 'poisson')   fit_vals <- lapply(coefs, function(ests) exp(X_pred%*%ests))
  if(response == 'bernoulli') fit_vals <- lapply(coefs, function(ests) exp(X_pred%*%ests)/(1+exp(X_pred%*%ests)))
  CB <- lapply(fit_vals, function(x) matrixStats::rowQuantiles(x, probs = c((1-conf_level/100)/2, (1-(1-conf_level/100)/2))))
  lower <- as.vector(do.call(cbind, lapply(CB, function(x) x[,1])))
  upper <- as.vector(do.call(cbind, lapply(CB, function(x) x[,2])))
  return(list(lower = lower, upper = upper))
}
