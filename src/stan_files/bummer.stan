data {
  int<lower=1> N;
  int<lower=1> d;
  vector[N] X;
  matrix[N, d] y;
  int<lower=1> N_pred;
  matrix[N_pred, d] y_pred;
  int<lower=1> N_grid;
  vector[N_grid] X_grid;
} 
parameters {
  vector<lower=0>[d] a;
  vector[d] mu;
  vector<lower=0>[d] sigma2;  
}
transformed parameters {
  matrix<lower=0>[N, d] alpha;
  for (i in 1:N) {
    for (j in 1:d) {
      alpha[i, j] =  a[j] * exp( - pow(mu[j] - X[i], 2.0) / sigma2[j]);
    }
  }  

}
model {
  a ~ normal(0, 1);
  mu ~ normal(0, 1);
  sigma2 ~ cauchy(0, 1);
  for (i in 1:N) {
    real alpha_sum;
    vector[d] alpha_plus_y;
    vector[d] log_gamma_alpha_plus_y;
    vector[d] log_gamma_alpha;

    alpha_sum = sum(alpha[i, ]);
    alpha_plus_y = alpha[i, ]' + y[i, ]';
    for (j in 1:d) {
      log_gamma_alpha_plus_y[j] = lgamma(alpha_plus_y[j]);
      log_gamma_alpha[j] = lgamma(alpha[i, j]);
    }
    target += lgamma(alpha_sum) + sum(log_gamma_alpha_plus_y) - 
              lgamma(alpha_sum + sum(y[i, ])) - sum(log_gamma_alpha);
  }
}
generated quantities {
  matrix<lower=0>[N_grid, d] alpha_pred;
  vector[N_pred] X_pred;
  
  for (i in 1:N_grid) {
    for (j in 1:d) {
      alpha_pred[i, j] =  a[j] * exp( - pow(mu[j] - X_grid[i], 2.0) / sigma2[j]);
    }
  }  

  for (i in 1:N_pred) {
    int idx_pred;
    vector[N_grid] ps;
    for (j in 1:N_grid) {
      real alpha_pred_sum;
      vector[d] alpha_pred_plus_y_pred;
      vector[d] log_gamma_alpha_pred_plus_y_pred;
      vector[d] log_gamma_alpha_pred;
      alpha_pred_sum = sum(alpha_pred[j, ]);
      alpha_pred_plus_y_pred = alpha_pred[j, ]' + y_pred[i, ]'; 
      for (k in 1:d) {
        log_gamma_alpha_pred_plus_y_pred[k] = lgamma(alpha_pred_plus_y_pred[k]);
        log_gamma_alpha_pred[k] = lgamma(alpha_pred[j, k]);
      }
      ps[j] = lgamma(alpha_pred_sum) + sum(log_gamma_alpha_pred_plus_y_pred) - 
                  lgamma(alpha_pred_sum + sum(y_pred[i, ])) - 
                  sum(log_gamma_alpha_pred);
    }
    idx_pred = categorical_rng(softmax(ps));
    X_pred[i] = X_grid[idx_pred];
  }
}
