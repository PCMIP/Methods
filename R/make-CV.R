#' Cross-validation of Compositional datasets
#'
#' @export
#' @param y Numeric matrix of compositional count values.
#' @param X Numeric vector of climate values.
#' @param y_pred Numeric matrix of compositional count values for prediction.
#' @param df Integer for the degrees of freedom in the basis expansion
#' @param n_samples Number of samples passed to Stan
#' @param n_grid Number of grid points for importance sampling prediction
#' @param algorithm Variational algorithm of meanfield or fullrank for Stan
#' @param pooled Hierarchically pooled version of algorithm 
#' @param vb Run variational Bayes or NUTS HMC sampler
#' @param mcore Enable multicore on a local maching
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of 
#'
#'
#'
make_CV <- function(y, X, model_name, kfold,
                    ## Bayesian model options
                    n_samples=1000, n_grid=1000,
                    algorithm="fullrank",
                    pooled=TRUE, vb=TRUE, mcore=TRUE) {
  library(magrittr)
  library(dplyr)
  # y=y[1:10,]
  # X=X[1:10]
  mu_X <- mean(X)
  sd_X <- sd(X)
  X_center <- (X - mu_X) / sd_X
  
  N <- nrow(y)
  ## randomly permute the sample indices into folds
  folds <- cut(sample(1:N, N), breaks=kfold, labels=FALSE)
  
  ## run the cross-validation subroutine, can parallelize this later
  # foo <- seq(1, kfold)
  # library(dplyr)
  # library(purrr)
  
  # join_lists <- function(mylist) {
  #   do.call(mapply, c(cbind, mylist))
  # }
  # join_lists <- function(mylist) {
  #   Reduce(function(x,y) Map(cbind, x, y), mylist)
  # }

  
  out <- (1:kfold)  %>%
    purrr::map(function(i){
      make_CV_fold(i = i, model_name=model_name, y=y, X=X_center,
                   folds=folds, n_samples=n_samples, n_grid=n_grid,
                   algorithm=algorithm,
                   pooled=pooled, vb=vb, mcore=mcore)
    # }) %>% join_lists()#bind_rows()
    })
  #%>% bind_rows()
  
  ## deal with different size of Bayesian prediction object
  if (model_name == "BUMMER" | model_name=="MVGP") {
    X_pred <- do.call("rbind", lapply(out, "[[", 6))
    out <- lapply(out, function(x){ 
      x[[6]] <- NULL
      return(x)
    })
  }
  
  ## bind rows of out to produce a final data.frame
  out <- bind_rows(out)
  
  ## add in posterior predictive distribution
  if (model_name == "BUMMER" | model_name=="MVGP") {
    out$X_pred <- X_pred
  }
  

  
  # out <- make_CV_fold(1, model_name = model_name, y=y, X=X, folds=folds)
  ## 
  
  return(out)
}


#' Cross-validation of Compositional datasets
#'
#' @export
#' @param i Cross-validation fold
#' @param model_name Model being fit
#' @param y Numeric matrix of compositional count values.
#' @param X Numeric vector of climate values.
#' @param y_pred Numeric matrix of compositional count values for prediction.
#' @param df Integer for the degrees of freedom in the basis expansion
#' @param n_samples Number of samples passed to Stan
#' @param n_grid Number of grid points for importance sampling prediction
#' @param algorithm Variational algorithm of meanfield or fullrank for Stan
#' @param pooled Hierarchically pooled version of algorithm 
#' @param vb Run variational Bayes or NUTS HMC sampler
#' @param mcore Enable multicore on a local maching
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of 
#'
#'
#'
make_CV_fold <- function (i, model_name, y, X, folds, 
                          n_samples=1000, n_grid=1000,
                          algorithm="fullrank",
                          pooled=TRUE, vb=TRUE, mcore=TRUE) {#, ...) {
  
  # library(rioja)
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
                     nboot=1000)#, ...)
  } else if (model_name=="WAPLS") {
    out <- fit_WAPLS_CV(y_train_prop, y_test_prop, X_train, X_test, sse=TRUE, 
                        nboot=1000, npls=5)#, ...)
  } else if (model_name=="MAT") {
    out <- fit_MAT_CV(y_train_prop, y_test_prop, X_train, X_test, k=3,
                      sse=TRUE, nboot=1000, lean=FALSE)#, ...)
  } else if (model_name=='MLRC'){
    out <- fit_MLRC_CV(y_train_prop, y_test_prop, X_train, X_test, sse=TRUE, nboot=0100)
  } else if (model_name=='RF'){
    out <- fit_RF_CV(y_train_prop, y_test_prop, X_train, X_test, sse=TRUE, nboot=1000)
  } else if (model_name=='BUMMER') {
    out <- fit_BUMMER_CV(y_train, y_test, X_train, X_test, 
                         n_samples=n_samples, n_grid=n_grid,
                         algorithm=algorithm, pooled=pooled, 
                         vb=vb, mcore=mcore)
  } else if (model_name =='MVGP') {
    
  } else {
    stop("Only available models are WA, MAT, MLRC, RF, BUMMER, or MVGP")
  }
  if (model_name=='BUMMER' | model_name=='MVGP') {
    return(list(CRPS=out$CRPS, MSPE=out$MSPE, MAE=out$MAE, 
                      coverage=out$coverage, observations=out$observations, 
                      X_pred=out$X_pred))    
  } else {
    return(list(CRPS=out$CRPS, MSPE=out$MSPE, MAE=out$MAE, 
                      coverage=out$coverage, observations=out$observations, 
                      mu=out$mu, sd=out$sd))
  }
}
