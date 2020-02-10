\name{gbtm_monotone}
\alias{gbtm_monotone}
\title{Fits Monotone Polynomials to Group Based Trajectory Models}
\usage{
gbtm_monotone()
}
\arguments{
  \item{data}{a (n x T) matrix containing the observations for n individuals over T periods of time}
  \item{n_gps}{the number of groups to fit}
  \item{x}{if the time-varying covariate, i.e. time or age, is the same for all individuals then the user can input a vector containing these times, e.g. x = c(10,11,12) if all individuals have measurements at ages 10, 11 and 12 otherwise the user is required to input a (n x T) design matrix containing all times for all individuals}
  \item{poly_degs}{the degrees of the polynomials to fit to each group, e.g. c(2,2,3) will fit two quadratic groups and one cubic group (note the number of polynomial degrees must match the number of groups)}
  \item{n_starts}{number of sets of initial values to generate (would recommend 5-10 sets unless the data is extremely complex)}
  \item{monotone}{is TRUE if fitting monotone trajectories, FALSE otherwise}
  \item{covariates}{if using time-stable covariates input a (n x (q+1)) matrix of q covariates}
  \item{plot_it}{if TRUE will plot the estimated trajectories}
  \item{conf_level}{computes confidence bands (at level conf_level) and plots with trajectories (note plot_it must also be TRUE)}
  \item{boot}{LEAVE as 'off', it will turn to 'on' if computing bootstrapped iterations for confidence bands}
  \item{response}{data type of observations (can be gaussian, poisson or bernoulli for gaussian, count and binary data respectively}
}
\value{
  \item{memberships}{matrix containing data with the estimated group membership probabilities}
  \item{pis}{matrix containing the estimated group proportions}
  \item{thetas}{list containing the parameter estimates for the polynomials}
  \item{sigma2}{estimates for the variance parameter}
  \item{gammas}{matrix containing estimates for the parameters governing the covariate structure}
  \item{log_likelihood}{log-likelihood value of the final model}
  \item{time}{length of time taken to run model}
  \item{summary}{tabulated summary of the output}
  \item{covariates_summary}{tabulated summary of the covariate structure output}
}
\description{
This function allows the user to fit group based trajectory models using
    monotonically constrained polynomials. This is useful when fitting models
    where the response curve is known to be monotone.
}
\examples{
gbtm_monotone(data, 4, x, c(1,2,3,3), monotone = FALSE)
gbtm_monotone(data, 3, x, c(3,3,3), conf_level = 95)
}
