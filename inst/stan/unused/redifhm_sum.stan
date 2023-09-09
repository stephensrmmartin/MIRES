functions {
#include /functions/common.stan
#include /functions/sum_to_zero.stan
#include /functions/ud.stan
  
}

data {
  int N; // Total number of observations
  int J; // Total number of indicators
  int K; // Number of groups

  array[N] int group; // Group indicator array

  matrix[N, J] x; // Indicator values

  int<lower=0, upper = 1> prior_only; // Whether to sample from prior only.
  int eta_cor_nonmi; // Not used for unidimensional.
  real hmre_mu;
  real<lower=0> hmre_scale;
}

transformed data {
  int total = 3*J;
  array[total] int hm_item_index = gen_item_indices(J);
  array[total] int hm_param_index = gen_param_indices(J);
  array[3, J] int lamResNu_indices = gen_lamResNu_indices(J);
  vector[N*J] x_vector = to_vector(x);
}

parameters {
  // Fixed Effects
  row_vector[J] lambda_log;
  row_vector[J] nu;
  row_vector[J] resid_log;

  // Random Effects (Lambda, resid, nu, eta mean, eta logsd)
  matrix[K, total] random_z; // 3J Vectors of uncor, std. REs.
  cholesky_factor_corr[total] random_L; // Chol. factor of RE correlations
  vector<lower=0>[total] random_sigma; // SD of REs

  // Latent
  vector[N] eta_z; // N-length vector of latent factor scores; std.
  vector[K-1] eta_mean_s;  // K-1 latent means; K-th is determined so that mean(eta_mean) = 0
  vector<lower=0>[K-1] eta_sd_s; // K-1 latent SDs; K-th is determined so that prod(eta_sd) = 1

  // Hierarchical Inclusion Model
  real hm_tau;
  vector[3] hm_param;
  vector[J] hm_item;
  vector[total] hm_lambda;
  
}

transformed parameters {
  matrix[K, total] random = z_to_random(random_z, random_sigma, random_L);
  matrix[K, J] lambda_random = random[, lamResNu_indices[1]];
  matrix[K, J] resid_random = random[, lamResNu_indices[2]];
  matrix[K, J] nu_random = random[, lamResNu_indices[3]];
  row_vector[J] lambda_lowerbound = compute_lambda_lowerbounds(lambda_random);
  row_vector[J] lambda = exp(lambda_log) + lambda_lowerbound;
  vector[K] eta_mean = eta_means_stz(eta_mean_s);
  vector[K] eta_sd = eta_sds_pto(eta_sd_s);
  vector[N] eta = eta_mean[group] + eta_z .* eta_sd[group];
}

model {
  // Declarations
  vector[total] hm_hat = exp(hm_tau + hm_param[hm_param_index] + hm_item[hm_item_index] + hm_lambda);
  matrix[N, J] xhat = rep_matrix(nu, N) + eta*lambda; // Fixed effects
  matrix[N, J] s_loghat = rep_matrix(resid_log, N); // Fixed effects
  xhat += nu_random[group] + rep_matrix(eta, J) .* lambda_random[group];
  s_loghat += resid_random[group];

  // Priors
  /* Lambda:
     p(Lambda) = p(exp(lam) + lam_lb)|dl/dL|
     = p(exp(lam) + lam_lb)exp(lam)
     We want p(Lambda) to be N^+(0, 1).
     The following is the same as lambda ~ N^+(0, 1)T[lowerbound,] with jacobian
   */
  target += normal_lpdf(lambda | 0, 1) - normal_lccdf(lambda_lowerbound | 0, 1) + sum(lambda_log);
  resid_log ~ normal(0, 1);
  nu ~ normal(0, 1);

  to_vector(random_z) ~ std_normal();
  random_L ~ lkj_corr_cholesky(1);

  eta_z ~ std_normal();
  eta_mean_s ~ std_normal();
  eta_sd_s ~ std_normal();

  hm_tau ~ normal(hmre_mu, hmre_scale);
  hm_param ~ normal(hmre_mu, hmre_scale);
  hm_item ~ normal(hmre_mu, hmre_scale);
  hm_lambda ~ normal(hmre_mu, hmre_scale);

  // Hierarchical inclusion
  random_sigma ~ normal(0, hm_hat);


  // Likelihood
  if(!prior_only) {
    x_vector ~ normal(to_vector(xhat), to_vector(exp(s_loghat)));
  }
  
}


generated quantities {
  /* corr_matrix[total] RE_cor = L_to_cor(random_L); */
  matrix[total,total] RE_cor = L_to_cor(random_L);
}
