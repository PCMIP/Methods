data {
  int<lower=1> N;
  int<lower=1> d;
  int<lower=1> df;
  matrix[N, df] Xbs;
  matrix[N, d] y;
  int<lower=1> N_pred;
  matrix[N_pred, d] y_pred;
  int<lower=1> N_grid;
  vector[N_grid] X_grid;
  matrix[N_grid, df] Xbs_grid;
} 
parameters {
  matrix[df, d-1] beta;
  vector[df] mu_beta;
  cholesky_factor_corr[df] L_Omega_beta;
  vector<lower=0>[df] sigma_beta;
  // matrix[df, d] beta;
}
transformed parameters {
  matrix<lower=0>[N, d] alpha;
  for (i in 1:N) {
    for (j in 1:(d-1)) {
      alpha[i, j] = exp(Xbs[i, ] * beta[, j]);
    } 
    alpha[i, d] = 1.0;
    // for (j in 1:d) {
    //   alpha[i, j] = exp(Xbs[i, ] * beta[, j]);
      
      // // tobit link, can add exponential later using if statement
      // real link = Xbs[i, ] * beta[, j];
      // if (link > 1e-10) {
      //   alpha[i, j] = link;
      // } else {
      //   alpha[i, j] = 1e-10;
      // }
    // }
  }  
}
model {
  matrix[df, df] L_Sigma_beta;
  for (j in 1:(d-1)) {
    beta[, j] ~ multi_normal_cholesky(mu_beta, L_Sigma_beta);
  }
  mu_beta ~ normal(0, 1);
  L_Omega_beta ~ lkj_corr_cholesky(2);
  sigma_beta ~ cauchy(0, 1);
  L_Sigma_beta = diag_pre_multiply(sigma_beta, L_Omega_beta);
  
  
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
  // vector[N_pred] idx_pred;
  vector[N_pred] X_pred;
  
  for (i in 1:N_grid) {
    for (j in 1:(d-1)) {
      alpha_pred[i, j] = exp(Xbs_grid[i, ] * beta[, j]);
    } 
    alpha_pred[i, d] = 1.0;
    // for (j in 1:d) {
    //   alpha_pred[i, j] = exp(Xbs_grid[i, ] * beta[, j]);
      
      // real link_pred = Xbs_grid[i, ] * beta[, j];
      // if (link_pred > 1e-10) {
      //   alpha_pred[i, j] = link_pred;
      // } else {
      //   alpha_pred[i, j] = 1e-10;
      // }
    // }
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
    // idx_pred[i] = categorical(softmax(pis[i, ]));
    // X_pred[i] = X_grid[idx_pred[i]]
    idx_pred = categorical_rng(softmax(ps));
    X_pred[i] = X_grid[idx_pred];
  }
}
