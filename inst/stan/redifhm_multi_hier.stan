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
  vector[total_lambda] lambda_log_est;
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
  matrix[K, total_lambda] lambda_random = lambda_resid_nu_mean_logsd_random[, lamResNu_bounds[1,1]:lamResNu_bounds[1,2]];
  matrix[K, J] resid_random = lambda_resid_nu_mean_logsd_random[, lamResNu_bounds[2,1]:lamResNu_bounds[2,2]];
  matrix[K, J] nu_random = lambda_resid_nu_mean_logsd_random[, lamResNu_bounds[3,1]:lamResNu_bounds[3,2]];
  matrix[K, F] eta_mean = lambda_resid_nu_mean_logsd_random[, (lamResNu_bounds[3,2] + 1):(lamResNu_bounds[3,2] + F)];
  matrix[K, F] eta_sd = exp(lambda_resid_nu_mean_logsd_random[, (total_param + F + 1):(total_param + 2*F)]);
  // Compute etas; pre-compute the eta covariances, and allow the eta_cor_nonmi option.
  matrix[2*F, 2*F] eta_cov_U[K];
  matrix[N, F] eta = eta_mean[group]; // Initialize to means; eta = 0 + mean + stoch. error
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
}

model {
  
}

generated quantities {
  
}
