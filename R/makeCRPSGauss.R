#' Function to evalutate the Continous Ranked Probability Score
#' for a prediction that has a Gaussian distribution
#' This function is useful for WA, WAPLS, and MAT 
#' (and other potential methods that use bootstrapped uncertainty)
#'
#' @export
#' @param mu Numeric vVctor of predictive means for the hold out data
#' @param sdev Numeric Vector of predictive standard deviations for the hold out data
#' @param truth Numeric Vector of held out climate variables and is the target of the prediction
## of samples passed to Stan
#' @return A Numeric Vector of each held-out sample's CRPS score using a Gaussian predictive distribution 
#'
makeCRPSGauss <- function(mu, sdev, truth) {
  - sdev * (1 / sqrt(pi) - 2 * dnorm((truth- mu)/sdev) -
            (truth- mu)/sdev * (2 * pnorm((truth- mu)/sdev) - 1))
}
