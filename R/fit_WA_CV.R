# fits weighted-averaging model on cross-validation data
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

fit_WA_CV <- function(y_train_prop, y_test_prop, X_train, X_test, 
                      sse=TRUE, nboot=1000) {#, ...) {
  ## WA reconstruction - subset to deal with all zero occurrence species
  zeros_idx <- which(colSums(y_train_prop) == 0)
  if (length(zeros_idx) > 0) {
    modWA <- rioja::WA(y_train_prop[, - zeros_idx], X_train)#, ...)
    predWA <- predict(modWA, y_test_prop[, - zeros_idx], sse=sse, nboot=nboot)#, ...)
  } else {
    ## no data to subset
    modWA <- rioja::WA(y_train_prop, X_train)#, ...)
    predWA <- predict(modWA, y_test_prop, sse=sse, nboot=nboot)#, ...)
  }
  CRPS <- makeCRPSGauss(predWA$fit.boot[, 1],
                        sqrt(predWA$v1.boot[, 1]^2 + predWA$v2.boot[1]^2), X_test)
  MAE <- abs(predWA$fit.boot[, 1] - X_test)
  MSPE <- (predWA$fit.boot[, 1] - X_test)^2
  coverage <- (X_test >=
                 (predWA$fit.boot[, 1] - 2*sqrt(predWA$v1.boot[, 1]^2 + predWA$v2.boot[1]^2)) &
                 (X_test <= (predWA$fit.boot[, 1] + 2*sqrt(predWA$v1.boot[, 1]^2 + predWA$v2.boot[1]^2))))
  return(data.frame(MSPE=MSPE, MAE=MAE, CRPS=CRPS, coverage=coverage, observations=X_test, 
                    mu=predWA$fit.boot[,1], sd=sqrt(predWA$v1.boot[, 1]^2 + predWA$v2.boot[1]^2)))
}

