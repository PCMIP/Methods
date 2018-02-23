##
##
##
##
##
##
make_CV <- function(y, X, model_name, kfold) {

  # y=y[1:10,]
  # X=X[1:10]
  
  N <- nrow(y)
  ## randomly permute the sample indices into folds
  folds <- cut(sample(1:N, N), breaks=kfold, labels=FALSE)

  ## run the cross-validation subroutine, can parallelize this later
  # foo <- seq(1, kfold)
  library(dplyr)
  library(purrr)
  out <- (1:kfold)  %>% 
    purrr::map(function(i){ make_CV_fold(i = i, model_name=model_name, y=y, X=X,
                folds=folds)}) %>% bind_rows()
  
  # out <- make_CV_fold(1, model_name = model_name, y=y, X=X, folds=folds)
  ## 

  return(out)
}



make_CV_fold <- function (i, model_name, y, X, folds) {#, ...) {

  library(rioja)
  # library(analogue)

  ## setup cross-validation training and test data
  idx_test <- which(folds == i, arr.ind=TRUE)
  y_train <- as.matrix(y[-idx_test, ])
  y_train_prop <- counts2proportions(y_train)
  y_test <- as.matrix(y[idx_test, ])
  y_test_prop <- counts2proportions(y_test)
  X_train <- c(X[-idx_test])
  X_test <- c(X[idx_test])

  if (model_name=="WA") {
    out <- fit_WA_CV(y_train_prop, y_test_prop, X_train, X_test, sse=TRUE, 
                     nboot=100)#, ...)
  } else  if (model_name=="MAT") {
    out <- fit_MAT_CV(y_train_prop, y_test_prop, X_train, X_test, k=3,
                      sse=TRUE, nboot=100, lean=FALSE)#, ...)
  } else if (model_name=='MLRC'){
    out <- fit_MLRC_CV(y_train_prop, y_test_prop, X_train, X_test, sse=TRUE, nboot=100)
  } else if (model_name=='RF'){
    out <- fit_RF_CV(y_train_prop, y_test_prop, X_train, X_test, sse=TRUE, nboot=100)
  }
  return(data.frame(CRPS=out$CRPS, MSPE=out$MSPE, MAE=out$MAE, coverage=out$coverage, observations=out$observations, mu=out$mu, sd=out$sd))
}
