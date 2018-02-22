##
##
##
##
##
##
make_CV <- function(y, X, model_name, k_fold) {

  N <- nrow(y)
  ## randomly permute the sample indices into folds
  folds <- cut(sample(1:N, N), breaks=k_fold, labels=FALSE)

  ## run the cross-validation subroutine, can parallelize this later
  out <- sapply(1:k_fold, make_CV_fold, model_name=model_name, y=y, X=X,
                folds=folds)

  return(list(CRPS=out$CRPS, MSPE=out$MSPE, MAE=out$MAE, coverage=out$coverage))
}



make_CV_fold <- function (i, model_name=model_name, y=y, X=X, folds) {

  library(rioja)
  # library(analogue)

  ## setup cross-validation training and test data
  idx_test <- which(folds == i, arr.ind=TRUE)
  y_train <- as.matrix(y[-idx_test, ])
  y_train_prop <- counts2proportions(y_train)
  y_test <- as.matrix(y[idx_test, ])
  y_test_prop <- counts2proportions(y_test)
  X_train <- c(X_cv[-idx_test])
  X_test <- c(X_cv[idx_test])

  if (model_name=="WA") {
    out <- fit_WA_CV(y_train_prop, X_train, sse=TRUE, nboot=1000, ...)
  } else  if (model_name=="MAT") {
    out <- fit_MAT_CV(y_train_prop, X_train, sse=TRUE, nboot=1000, ...)
  }
  return(list(CRPS=out$CRPS, MSPE=out$MSPE, MAE=out$MAE, coverage=out$coverage))
}
