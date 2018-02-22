fit_MAT_CV <- function(y_train_prop, X_train, k=k, sse=TRUE, nboot=1000, lean=FALSE, ...) {
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
