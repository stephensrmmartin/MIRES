functions {
#include /functions/common.stan
#include /functions/md.stan
  
}

data {
  int N; // Total number of observations
  int J; // Total number of indicators
  int K; // Number of groups
  int F; // Number of latent factors

  int group[N]; // Group indicator array

  matrix[N, J] x; // Indicator values

  int J_f[F]; // Number of indicators per factor
  int F_ind[F, J]; // Which indicators for each factor.

  // Options
  int<lower=0, upper=1> prior_only; // Whether to sample from prior only.
  int<lower=0, upper=1> eta_cor_nonmi; // Whether correlations between factors 1:F is not invariant across 1:K groups.
  
}

transformed data {
  int total_lambda = sum(J_f);
  int total_param = total_lambda + 2*J;

  int hm_item_index[total_param] = gen_item_indices_md(J, F, J_f, F_ind);
  int hm_param_index[total_param] = gen_param_indices_md(J, J_f);
  int lamResNu_bounds[3, 2] = gen_lamResNu_bounds(J, J_f);

  vector[N*J] x_vector = to_vector(x);
}


parameters {
  // Fixed Effects
  row_vector[total_lambda] lambda_log_est; // TODO: Separate cross-loadings if possible; these should be able to go negative. Or, free all and allow post-process in GQ to ensure a reference loading is positive (note: Not *first* loading, because that could be a cross-loading).
  row_vector[J] nu;
  row_vector[J] resid_log;

  // Random Effects
  matrix[K, total_param + 2*F] lambda_resid_nu_mean_logsd_z;
  cholesky_factor_corr[total_param + 2*F] lambda_resid_nu_mean_logsd_L;
  vector<lower = 0>[total_param + 2*F] lambda_resid_nu_mean_logsd_sigma; // Assess whether it's ID'd to estimate mean_logsd sigma

  // Latent
  matrix[N, F] eta_z; // NxF matrix of latent factor scores.
  cholesky_factor_corr[F] eta_L_fixed; // Latent factor correlations
  // If not invariant, define:
  cholesky_factor_corr[F] eta_L_random[K * eta_cor_nonmi];
  vector<lower=0, upper = 1>[eta_cor_nonmi] eta_L_random_weight;


  // Hierarchical Inclusion Model
  real hm_tau;
  vector[3] hm_param;
  vector[J] hm_item;
  vector[total_param] hm_lambda;
  // vector[F] hm_factor; 
  /* hm_factor: Two ideas. 1) Whether the RESD in question is relevant to a factor. Loadings dictate this. 2) Would encode which factor the param belongs, but only loadings differ between factors, so without lots of cross-loadings, it'd be unidentifiable.
     Best plan is just to use it to encode which factor(s) are relevant for a given item. A 'factor' effect. Really like two hm_tau's.
   */
  
}

transformed parameters {
  matrix[K, total_param + 2*F] lambda_resid_nu_mean_logsd_random = z_to_random(lambda_resid_nu_mean_logsd_z, lambda_resid_nu_mean_logsd_sigma, lambda_resid_nu_mean_logsd_L);
  matrix[K, total_lambda] lambda_est_random = lambda_resid_nu_mean_logsd_random[, lamResNu_bounds[1,1]:lamResNu_bounds[1,2]];
  matrix[K, J] resid_random = lambda_resid_nu_mean_logsd_random[, lamResNu_bounds[2,1]:lamResNu_bounds[2,2]];
  matrix[K, J] nu_random = lambda_resid_nu_mean_logsd_random[, lamResNu_bounds[3,1]:lamResNu_bounds[3,2]];
  matrix[K, F] eta_mean = lambda_resid_nu_mean_logsd_random[, (lamResNu_bounds[3,2] + 1):(lamResNu_bounds[3,2] + F)];
  matrix[K, F] eta_sd = exp(lambda_resid_nu_mean_logsd_random[, (total_param + F + 1):(total_param + 2*F)]);
  // Compute etas; pre-compute the eta covariances, and allow the eta_cor_nonmi option.
  matrix[2*F, 2*F] eta_cov_U[K];
  matrix[N, F] eta = eta_mean[group]; // Initialize to means; eta = 0 + mean + stoch. error
  row_vector[total_lambda] lambda_lowerbound = compute_lambda_lowerbounds(lambda_est_random);
  row_vector[total_lambda] lambda_est = exp(lambda_log_est) + lambda_lowerbound;
  matrix[F, J] lambda = lambda_mat(J_f, F_ind, lambda_est);
  matrix[F, J] lambda_random[K];
  // TODO: Compute lambda_lowerbounds using vector; then construct lambda_mat and lambda_mat_random

  // Compute stochastic latent errors
  for(k in 1:K) { // Compute upper-triangular cholesky-covariance for each group.
    if(eta_cor_nonmi) { // Non-invariant correlations
      eta_cov_U[k] = (diag_pre_multiply(eta_sd[k], convex_combine_Ls(eta_L_fixed, eta_L_random[k], eta_L_random_weight)))';
    } else {
      eta_cov_U[k] = (diag_pre_multiply(eta_sd[k], eta_L_fixed))';
    }
  }
  for(n in 1:N) {
    eta[n] += eta_z[n] * eta_cov_U[group[n]];
  }

  // Restructure lambda_est_random to matrices.
  for(k in 1:K) {
    lambda_random[k] = lambda_mat(J_f, F_ind, lambda_est_random[k]);
  }
}

model {
  // Declarations
  vector[total_param] hm_hat = exp(hm_tau + hm_param[hm_param_index] + hm_item[hm_item_index] + hm_lambda);
  matrix[N, J] xhat = rep_matrix(nu, N) + eta * lambda; // Fixed effects
  matrix[N, J] s_loghat = rep_matrix(resid_log, N);
  // Random effects
  xhat += nu_random[group];
  s_loghat += resid_random[group];
  for(n in 1:N) {
    xhat[n] += eta[n] * lambda_random[group[n]];
  }

  // Priors
  target += normal_lpdf(lambda_est | 0, 1) - normal_lccdf(lambda_lowerbound | 0, 1) + sum(lambda_log_est);
  resid_log ~ std_normal();
  nu ~ normal(0, 1);

  to_vector(lambda_resid_nu_mean_logsd_z) ~ std_normal();
  lambda_resid_nu_mean_logsd_L ~ lkj_corr_cholesky(1);
  lambda_resid_nu_mean_logsd_sigma ~ std_normal();

  to_vector(eta_z) ~ std_normal();
  eta_L_fixed ~ lkj_corr_cholesky(1);
  if(eta_cor_nonmi) {
    // eta_L_random_weight ~ uniform(0,1)
    // Maybe use ~ beta() to prefer pooling.
    for(k in 1:K) {
      /*
	1 / alpha: When alpha |-> 1, lkj is uniform, and no pooling occurs.
	When alpha |-> 0, lkj is very peaked over identity, and full pooling occurs.
       */
      eta_L_random[k] ~ lkj_corr_cholesky(1.0 / eta_L_random_weight[1]); // 
    }
  }

  hm_tau ~ std_normal();
  hm_param ~ std_normal();
  hm_item ~ std_normal();
  hm_lambda ~ std_normal();

  // Hierarchical inclusion
  lambda_resid_nu_mean_logsd_sigma[1:total_param] ~ normal(0, hm_hat);

  // Priors on the SD of eta-mean and logsd.
  lambda_resid_nu_mean_logsd_sigma[(total_param + 1) : rows(lambda_resid_nu_mean_logsd_sigma)] ~ normal(0, 1);

  if(!prior_only) {
    x_vector ~ normal(to_vector(xhat), to_vector(exp(s_loghat)));
  }
  
}

generated quantities {
  
}
