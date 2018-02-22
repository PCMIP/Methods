# fits modern-analog model on cross-validation data
# y_train_prop A data frame of relative abundances of size N_site by N_taxa 
# X_train A data frame of the climate covariate of size N_site by 1
# k is the number of analogs
# sse Sum of Squared Error ...
# n.boot Number of bootstrapped samples for prediction
# lean verbose model output
#
# returns a list of cross-validation statistics for each held out sample
# MSPE Squared Prediction Error
# MAE Absolute Error
# CRPS Continuous Ranked Probability Score
# coverage Empirical 95% coverage rate

fit_MAT_CV <- function(y_train_prop, X_train, k=k, sse=TRUE, n.boot=1000, lean=FALSE, ...) {
  ## Modern analogue technique
  modMAT <- rioja::MAT(as.data.frame(y_train_prop), X_train, k=k, lean=lean, ...)
  predMAT <- predict(modMAT, as.data.frame(y_test_prop), k=k, sse=sse, n.boot=n.boot, ...)
  CRPS <- makeCRPSGauss(
    predMAT$fit.boot[, 2],
    sqrt(predMAT$v1.boot[, 2]^2+ predMAT$v2.boot[2]),X_test)
  MSPE <- (predMAT$fit.boot[, 2] - X_test)^2
  MAE <- abs(predMAT$fit.boot[, 2] - X_test)
  coverage <-
    (X_test >= (predMAT$fit.boot[, 2] -
                    2 * sqrt(predMAT$v1.boot[, 2]^2+ predMAT$v2.boot[2])) &
        (X_test <= (predMAT$fit.boot[, 2] +
                      2*  sqrt(predMAT$v1.boot[, 2]^2+ predMAT$v2.boot[2]))))
  return(list(MSPE=MSPE, MAE=MAE, CRPS=CRPS, coverage=coverage))
}
