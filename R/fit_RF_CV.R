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
# fits random forest model on cross-validation data
# y_train_prop A data frame of relative abundances of size N_site by N_taxa 
# X_train A data frame of the climate covariate of size N_site by 1
# sse Sum of Squared Error ...
# nboot Number of bootstrapped samples for prediction
#
# returns a list of cross-validation statistics for each held out sample
# MSPE Squared Prediction Error
# MAE Absolute Error
# CRPS Continuous Ranked Probability Score
# coverage Empirical 95% coverage rate

fit_RF_CV <- function(y_train_prop, y_test_prop, X_train, X_test, 
                      sse=TRUE, nboot=1000) {#, ...) {
  library(randomForest)
  library(scoringRules)
  ## Random Forest
  train <- data.frame(covar=X_train, y_train_prop)
  test <- data.frame(y_test_prop)
  rf <- randomForest(covar ~ ., data = train)
  # CRPS <- makeCRPSGauss(t(matrix(predict(rf, test, predict.all=TRUE)$individual, 
  #                           length(idx_test), 500)), X_test)
  # CRPS <- abs(predict(rf, test) - X_test)
  CRPS <- crps_sample(X_test, predict(rf, test, predict.all=TRUE)$individual)
  MSPE <- (predict(rf, test) - X_test)^2
  MAE <- abs(predict(rf, test) - X_test)
  rf_CI <- t( apply( predict(rf, test, predict.all=TRUE)$individual, 1,
                     function(x) {
                       quantile(x, c(0.025,0.975))
                     }))
  rf_sd <- apply( predict(rf, test, predict.all=TRUE)$individual, 1,
                     function(x) {
                       sd(x)
                     })
  coverage <- ( (X_test >= rf_CI[, 1]) & (X_test <= rf_CI[, 2]) )
  return(data.frame(MSPE=MSPE, MAE=MAE, CRPS=CRPS, coverage=coverage, observations=X_test, 
                    mu=predict(rf, test), sd=rf_sd))
}