#' Bayesian BUMMER with Stan
#'
#' @export
#' @param y Numeric matrix of compositional count values.
#' @param X Numeric vector of climate values.
#' @param output_samples Number of samples passed to Stan
#' @param algorithm Variational algorithm of meanfield or fullrank for Stan
#' @param n_samples Number of samples
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
# fits BUMMER model on cross-validation data
# y_train A data frame of compositional count abundances of size N_site by N_taxa 
# X_train A data frame of the climate covariate of size N_site by 1
#
# returns a list of cross-validation statistics for each held out sample
# MSPE Squared Prediction Error
# MAE Absolute Error
# CRPS Continuous Ranked Probability Score
# coverage Empirical 95% coverage rate

fit_BUMMER_CV <- function(y_train, y_test, X_train, X_test, 
                          n_samples=1000, n_grid=1000,
                          algorithm="fullrank",
                          pooled=TRUE, vb=TRUE, mcore=TRUE) {#, ...) {
  library(rstan)
  library(scoringRules)
  if (mcore) {
    options(mc.cores = parallel::detectCores())
  }
  ## STAN BUMMER
  ## default to variational estimation for large datasets
  fit <- bummer_stan(y=y_train, X=X_train, y_pred=y_test,
                     n_samples=n_samples, n_grid=n_grid,
                     algorithm=algorithm,
                     pooled=pooled, vb=vb, mcore=mcore)
  e <- rstan::extract(fit)
  X_pred <- e$X_pred
  
  # CRPS <- makeCRPSGauss(t(matrix(predict(rf, test, predict.all=TRUE)$individual, 
  #                           length(idx_test), 500)), X_test)
  # CRPS <- abs(predict(rf, test) - X_test)
  # CRPS <- makeCRPS(X_pred, X_test)
  ## CRPS using scoringRules package
  CRPS <- crps_sample(X_test, t(X_pred))
  MSPE <- (apply(X_pred, 2, mean) - X_test)^2
  MAE <- abs(apply(X_pred, 2, median) - X_test)
  CI <- t( apply( X_pred, 2,
                  function(x) {
                    quantile(x, c(0.025,0.975))
                  }))
  coverage <- ( (X_test >= CI[, 1]) & (X_test <= CI[, 2]) )
  out <- list(MSPE=MSPE, MAE=MAE, CRPS=CRPS, coverage=coverage, 
              observations=X_test, mu=apply(X_pred, 2, mean), 
              sd = apply(X_pred, 2, sd), X_pred=t(X_pred))
  
  ## weird trick to allow a matrix in a data.frame
  ## https://stackoverflow.com/questions/6143697/data-frame-with-a-column-containing-a-matrix-in-r?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
  class(out) <- "data.frame"
  names(out) <- c("MSPE", "MAE", "CRPS", "coverage", "observations", 
                  "mu", "sd", "X_pred")
  row.names(out) <- 1:length(X_test)
  return(out)
}
