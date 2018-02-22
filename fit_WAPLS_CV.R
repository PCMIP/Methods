# fits weighted-averaging partial least squares model on cross-validation data
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

fit_WAPLS_CV <- function(y_train_prop, X_train, sse=TRUE, nboot=1000, ...) {
  ## WAPLS reconstruction - subset to deal with all zero occurrence species
  zeros_idx <- which(colSums(y_train_prop) == 0)
  if (length(zeros_idx) > 0) {
    modWAPLS <- rioja::WAPLS(y_train_prop[, - zeros_idx], X_train)     
    predWAPLS <- predict(modWAPLS, y_test_prop, sse=TRUE, nboot=1000)
  } else {
    modWAPLS <- rioja::WAPLS(y_train_prop, X_train)     
    predWAPLS <- predict(modWAPLS, y_test_prop, sse=TRUE, nboot=1000)
  }
  CRPS <- makeCRPSGauss(predWAPLS$fit[, 1], sqrt(predWAPLS$v1.boot[, 1]),
                        X_test)
  MSPE <- (predWAPLS$fit[, 1] - X_test)^2
  MAE <- abs(predWAPLS$fit[, 1] - X_test)
  coverage <- (
    X_test >=
      (predWAPLS$fit[, 1] - 2*sqrt(predWAPLS$v1.boot[, 1]))) & 
    (X_test <= 
       (predWAPLS$fit[, 1] + 2 * sqrt(predWAPLS$v1.boot[, 1])))
  return(list(MSPE=MSPE, MAE=MAE, CRPS=CRPS, coverage=coverage))
}