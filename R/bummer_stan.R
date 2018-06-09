#' Bayesian BUMMER with Stan
#'
#' @export
#' @param y Numeric matrix of compositional count values.
#' @param X Numeric vector of climate values.
#' @param output_samples Number of samples passed to Stan
#' @param algorithm Variational algorithm of meanfield or fullrank for Stan
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#'
bummer_stan <- function(y, X, y_pred, n_grid=1000, 
                        n_samples=1000, algorithm="fullrank", 
                        pooled=TRUE, vb=TRUE, mcore=TRUE, ...) {
  if (mcore) {
    options(mc.cores = parallel::detectCores())
  }
  s_X <- sd(X)
  ## setup grid for importance sampling
  X_grid <- seq(from=min(X) - 1.25*s_X, 
                to=max(X) + 1.25*s_X, length=n_grid)
  
  dat_fit <- list(y=y, X=X, N=nrow(y), d=ncol(y), y_pred=y_pred, 
                  N_pred=nrow(y_pred), X_grid=X_grid, N_grid=length(X_grid))
  # str(dat_fit)
  if (vb) {
    ## variational model fit
    if (pooled) {
      inits <- list(censored_a=runif(dat_fit$d, 0.5, 1.5), 
                    mu_a = 0, sigma_a = 1, 
                    censored_sigma=runif(dat_fit$d, 0.5, 1.5), 
                    mu_sigma = 0, sigma_sigma = 1,
                    mu=runif(dat_fit$d, -0.5, 0.5))
      
      ## fit model
      out <- rstan::vb(stanmodels$bummer_pooled, data=dat_fit,
                       pars=c("alpha", "alpha_pred"), include=FALSE,
                       output_samples=n_samples, algorithm=algorithm,
                       init=inits, ...)
    } else {
      out <- rstan::vb(stanmodels$bummer, data=dat_fit,
                       pars=c("alpha", "alpha_pred"), include=FALSE,
                       output_samples=n_samples, algorithm=algorithm, ...)
    }
    # out <- sampling(stanmodels$bummer, data=dat_fit, ...)
  } else {
    ## MCMC model fit
    if (pooled) {
      ## Hierarchically pooled model
      inits <- list(censored_a=runif(dat_fit$d, 0.5, 1.5),
                    mu_a = 0, sigma_a = 1,
                    censored_sigma=runif(dat_fit$d, 0.5, 1.5),
                    mu_sigma = 0, sigma_sigma = 1,
                    mu=runif(dat_fit$d, -0.5, 0.5))
      
      out <- rstan::sampling(stanmodels$bummer_pooled, data=dat_fit,
                             pars=c("alpha", "alpha_pred"), include=FALSE,
                             chains=4, iter=n_samples, ...)
    } else {
      out <- rstan::sampling(stanmodels$bummer, data=dat_fit,
                             pars=c("alpha", "alpha_pred"), include=FALSE,
                             chains=4, iter=n_samples,  ...)
    }
  }
    return(out)
}
