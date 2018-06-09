set.seed(401)
library(Methods)
library(MCMCpack)
library(splines)
N <- 800
d <- 8
df <- 4
N_pred <- 50
X <- runif(N, -1, 1)
X_pred <- runif(N_pred, -1, 1)
s_X <- sd(X)
knots <- seq(min(X) - 0.125 * s_X, max(X) + 0.125 * s_X, length=df-2)
interior_knots <- knots[-c(1, df-2)]
Boundary_knots <- knots[c(1, df-2)]
Xbs <- bs(X, intercept = TRUE,
          knots = interior_knots,
          Boundary.knots = Boundary_knots)
Xbs_pred <- bs(X_pred, intercept = TRUE,
               knots = interior_knots,
               Boundary.knots = Boundary_knots)

# beta <- matrix(rnorm(df*d, 0, 1), df, d)
beta <- cbind(matrix(rnorm(df*(d-1), 0, 2), df, d-1), rep(0, d))

alpha <- matrix(0, N, d)
p <- matrix(0, N, d)
y <- matrix(0, N, d)
alpha_pred <- matrix(0, N_pred, d)
p_pred <- matrix(0, N_pred, d)
y_pred <- matrix(0, N_pred, d)

for (i in 1:N) {
  alpha[i, ] <- exp(Xbs[i, ] %*% beta)
  # alpha[i, ] <- pmax(mu + Xbs[i, ] %*% beta, 1e-10)
  # Dirichlet-Multinomial
  p[i, ] <- rdirichlet(1, alpha[i, ])
  y[i, ] <- rmultinom(1, 100, p[i, ])
}

for (i in 1:N_pred) {
  alpha_pred[i, ] <- exp(Xbs_pred[i, ] %*% beta)
  # alpha_pred[i, ] <- pmax(mu + Xbs_pred[i, ] %*% beta, 1e-10)
  p_pred[i, ] <- rdirichlet(1, alpha_pred[i, ])
  y_pred[i, ] <- rmultinom(1, 100, p_pred[i, ])
}


##
## Plot data
##

y_prop <- y
p_alpha <- alpha
for (i in 1:N) {
  y_prop[i, ] <- y_prop[i, ] / sum(y_prop[i, ])
  p_alpha[i, ] <- p_alpha[i, ] / sum(p_alpha[i, ])
}


plotData <- data.frame(
  species = as.factor(rep(1:d, each=N)),
  Count   = c(as.matrix(y_prop)), 
  Wetness = rep(X, times=d), 
  alpha        = c(p_alpha))
library(ggplot2)
ggplot(plotData, aes(x=Wetness, y=Count, color=species, group=species)) +
  geom_point(alpha=0.25) +
  geom_line(aes(x=Wetness, y=alpha, col = species), plotData, lwd=1.25) + 
  theme(legend.position="none") + ggtitle("Training Data") + 
  labs(x="Water Table Depth", y="Composition")  

## Fit the model

fit <- mvgp_stan(y, X, y_pred, df=df, n_grid=100, output_samples=500)
print(fit)
e <- rstan::extract(fit)

plot(apply(e$X_pred, 2, mean) ~ X_pred)
abline(a=0, b=1, col="red")



##
## Plot Model Fit
##

p_alpha_post <- array(0, dim(e$alpha))
for (k in 1:dim(e$alpha)[1]) {
  for (i in 1:N) {
    p_alpha_post[k, i, ] <- e$alpha[k, i, ] / sum(e$alpha[k, i, ])
  }  
}

plotFit <- data.frame(
  species = as.factor(rep(1:d, each=N)),
  Count   = c(as.matrix(y_prop)), 
  Wetness = rep(X, times=d), 
  alpha   = c(apply(p_alpha_post, c(2, 3), mean)),
  alpha_lower   = c(apply(p_alpha_post, c(2, 3), quantile, prob=0.025)),
  alpha_upper   = c(apply(p_alpha_post, c(2, 3), quantile, prob=0.975)))

library(ggplot2)
ggplot(plotFit, aes(x=Wetness, y=Count, color=species, group=species)) +
  geom_point(alpha=0.25) +
  # geom_line(aes(x=Wetness, y=alpha, col=species), plotFit, lwd=1.25) +
  geom_line(aes(x=Wetness, y=alpha, col=species), plotData, lwd=2.25) +
  geom_ribbon(aes(ymin=alpha_lower, ymax=alpha_upper, fill=species, group=species),
              linetype=0, alpha=0.5) + 
  theme(legend.position="none") + ggtitle("Training Data") + 
  labs(x="Water Table Depth", y="Composition")  

# timesTwo(100)
# 
# N <- 100
# y=matrix(rnorm(N^2), N, N) 
# X=rnorm(N)
# makeCRPS(y, X, N)
# 
# 
