# fits MLRC on cross-validation data
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

fit_MLRC_CV <- function(y_train_prop, y_test_prop, X_train, X_test, sse=TRUE, nboot=nboot, ...) {
  
  ## MLRC reconstruction - subset to deal with all zero occurrence species
  zeros_idx <- which(colSums(y_train_prop) == 0)
  if (length(zeros_idx) > 0) {
    modMLRC <- rioja::MLRC(y_train_prop[, - zeros_idx], as.numeric(X_train))
    predMLRC <- predict(modMLRC, y_test_prop[, - zeros_idx],
                        sse=sse, nboot=nboot)
  } else {
    modMLRC <- rioja::MLRC(y_train_prop, as.numeric(X_train))
    predMLRC <- predict(modMLRC, y_test_prop, sse=sse, nboot=nboot)
  }
  CRPS <- makeCRPSGauss(predMLRC$fit[, 1],
                        sqrt(predMLRC$v1.boot[, 1]^2 + predMLRC$v2.boot[1]^2),
                        X_test)
  MSPE <- (predMLRC$fit[, 1] - X_test)^2
  MAE <- abs(predMLRC$fit[, 1] - X_test)
  coverage <- ( X_test >= (predMLRC$fit[, 1] - 
                             2*sqrt(predMLRC$v1.boot[, 1]^2 + predMLRC$v2.boot[1]^2))) & 
    (X_test <= (predMLRC$fit[, 1] + 
                  2 * sqrt(predMLRC$v1.boot[, 1]^2 + predMLRC$v2.boot[1]^2)))
  return(data.frame(MSPE=MSPE, MAE=MAE, CRPS=CRPS, coverage=coverage, observations=X_test, mu=predMLRC$fit[, 1], 
              sd=sqrt(predMLRC$v1.boot[, 1]^2 + predMLRC$v2.boot[1]^2)))
}
