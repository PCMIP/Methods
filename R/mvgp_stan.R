#' Bayesian MVGP with Stan
#'
#' @export
#' @param y Numeric matrix of compositional count values.
#' @param X Numeric vector of climate values.
#' @param y_pred Numeric matrix of compositional count values for prediction.
#' @param df Integer for the degrees of freedom in the basis expansion
#' @param n_grid Number of grid locations for importance sampling
#' @param n_samples Number of samples passed to Stan
#' @param algorithm Variational algorithm of meanfield or fullrank for Stan
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#'
mvgp_stan <- function(y, X, y_pred, df=6, n_grid=1000, 
                      n_samples=1000, algorithm="fullrank", 
                      pooled=TRUE, vb=TRUE, mcore=TRUE, ...) {
  s_X <- sd(X)
  knots <- seq(min(X) - 0.125 * s_X, max(X) + 0.125 * s_X, length=df-2)
  interior_knots <- knots[-c(1, df-2)]
  Boundary_knots <- knots[c(1, df-2)]
  
  ## setup basis expansion
  Xbs  <- bs(X, intercept = TRUE,
             knots = interior_knots,
             Boundary.knots = Boundary_knots)
  
  ## setup grid for importance sampling
  X_grid <- seq(from=min(X) - 0.125 * s_X, 
                to=max(X) + 0.125 * s_X, length=n_grid)
  Xbs_grid  <- bs(X_grid, intercept = TRUE,
                  knots = interior_knots,
                  Boundary.knots = Boundary_knots)
  
  dat_fit <- list(y=y, Xbs=Xbs, N=nrow(y), d=ncol(y), y_pred=y_pred, 
                  N_pred=nrow(y_pred), X_grid=X_grid, N_grid=length(X_grid), 
                  Xbs_grid=Xbs_grid, df=df)
  # str(dat_fit)
  # out <- rstan::vb(stanmodels$mvgp, data=dat_fit,
  #                  output_samples=output_samples, algorithm=algorithm, ...)
  # out <- sampling(stanmodels$mvgp, data=dat_fit, ...)
  
  if (vb) {
    ## variational model fit
    if (pooled) {
      ## fit model
      out <- rstan::vb(stanmodels$mvgp_pooled, data=dat_fit,
                       pars=c("alpha", "alpha_pred"), include=FALSE,
                       output_samples=n_samples, algorithm=algorithm,
                       init=inits, ...)
    } else {
      out <- rstan::vb(stanmodels$mvgp, data=dat_fit,
                       pars=c("alpha", "alpha_pred"), include=FALSE,
                       output_samples=n_samples, algorithm=algorithm, ...)
    }
    # out <- sampling(stanmodels$bummer, data=dat_fit, ...)
  } else {
    ## MCMC model fit
    if (pooled) {
      ## Hierarchically pooled model
      out <- rstan::sampling(stanmodels$mvgp_pooled, data=dat_fit,
                             pars=c("alpha", "alpha_pred"), include=FALSE,
                             chains=4, iter=n_samples, ...)
    } else {
      out <- rstan::sampling(stanmodels$mvgp, data=dat_fit,
                             pars=c("alpha", "alpha_pred"), include=FALSE,
                             chains=4, iter=n_samples,  ...)
    }
  }
  return(out)
}
