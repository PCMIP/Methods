## Function to evalutate the Continous Ranked Probability Score
## for a prediction that has a Gaussian distribution
## This function is useful for WA, WAPLS, and MAT 
## (and other potential methods that use bootstrapped uncertainty)
##
## mu is a vector of predictive means for the hold out data
## sdev is a vector of predictive standard deviations for the hold out data
## truth is a vector of held out climate variables and is the target of the prediction
##
## The function outputs the sample specifice CRPS score
makeCRPSGauss <- function(mu, sdev, truth) {
  - sdev * (1 / sqrt(pi) - 2 * dnorm((truth- mu)/sdev) -
            (truth- mu)/sdev * (2 * pnorm((truth- mu)/sdev) - 1))
}
