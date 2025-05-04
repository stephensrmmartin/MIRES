functions {
#include /functions/common.stan
#include /functions/ud.stan
#include /functions/marg.stan
#include /functions/sum_to_zero.stan
  
}

data {
  int N; // Total number of observations (Rows)
  int J; // Total number of indicators (Cols)
  int K; // Number of groups

  array[N] int group; // Group indicator array

  matrix[N, J] x; // Indicator values

  // Model configuration
  int<lower=0, upper=1> prior_only; // Whether to sample from prior only

  int eta_cor_nonmi; // Not used for unidimensional.

  int<lower=0, upper=1> use_hmre; // Whether hier. inclusion model is used.
  real hmre_mu; // hmre location
  real<lower=0> hmre_scale; // hmre scale

  int<lower=0, upper=1> marginalize; // Whether to marginalize over eta (1), or estimate them (0).

  int<lower=0, upper=1> sum_coding; // Whether to use sum-to-one/prod-to-1 constraint (1) or hierarchicalize (0) on the latent means and SDs.
}

transformed data {
  int total = 3 * J; // Number of random measurement model effects (Not including latent mu/sd)
  array[total] int hm_item_index = gen_item_indices(J);
  array[total] int hm_param_index = gen_param_indices(J);
  array[3, J] int lamResNu_indices = gen_lamResNu_indices(J);
  int save_scores = 1 - marginalize; // Inverse boolean for easy reading
  int hier_coding = 1 - sum_coding; // Inverse boolean for easy reading
  vector[N * J] x_vector;
  array[N] row_vector[J] x_sorted;
  array[K,2] int x_sorted_indices;

  // Optimizations
  if(save_scores) { // Vectorize input data in col-major order. x_vector ~ N(to_vector(xhat), to_vector(shat)).
    x_vector = to_vector(x);
  } else { // Sort data into K chunks for a loop over K vectorized multi_normal.
    x_sorted = sort_data_by_group(x, group);
    x_sorted_indices = sort_data_by_group_indices(group);
  }
  
}

parameters {
  // Fixed Effects
  row_vector[J] lambda_log; // Loadings
  row_vector[J] nu; // Intercepts
  row_vector[J] resid_log; // (Log) Residual SD

  // Random Effects (Lambda, resid, nu, [eta-mean, eta-logsd])
  matrix[K, total + 2*hier_coding] random_z; // 3J (+2) vectors of uncorrelated, Std. REs.
  cholesky_factor_corr[total + 2*hier_coding] random_L; // Chol. factor of RE cors.
  vector<lower=0>[total] random_sigma; // Not including RE-SDs of latent mean/sd


  // Latent [Optional]
  vector[N*save_scores] eta_z; // N-length vector of untransformed latent factor scores.
  //// Hierarchical [Optional]
  vector<lower=0>[2*hier_coding] eta_random_sigma; // RE-SDs of random latent means and (log) SDs.

  //// Sum-coded [Optional]
  vector[(K-1)*sum_coding] eta_mean_s; // K-1 eta_means
  vector<lower=0>[(K-1)*sum_coding] eta_sd_s; // K-1 eta_SDs

  // Hierarchical Inclusion model
  real hm_tau; // Global regularization. No way to disable this *currently*; will ignore for non-re models. Only one param; shouldn't slow.
  vector[3*use_hmre] hm_param; // Parameter-wise regularization. [Optional]
  vector[J*use_hmre] hm_item; // Item-wise regularization. [Optional]
  vector[total] hm_lambda; // Parameter X Item regularization.
}

transformed parameters {
  vector<lower=0>[total + 2*hier_coding] random_sigma_all = append_row(random_sigma, eta_random_sigma);
  matrix[K, total + 2*hier_coding] random = z_to_random(random_z, random_sigma_all, random_L);
  matrix[K, J] lambda_random = random[, lamResNu_indices[1]];
  matrix[K, J] resid_random = random[, lamResNu_indices[2]];
  matrix[K, J] nu_random = random[, lamResNu_indices[3]];
  vector[K] eta_mean = sum_coding ? eta_means_stz(eta_mean_s) : random[, total + 1];
  vector[K] eta_sd = sum_coding ? eta_sds_pto(eta_sd_s) : exp(random[, total + 2]);
  row_vector[J] lambda_lowerbound = compute_lambda_lowerbounds(lambda_random);
  row_vector[J] lambda = exp(lambda_log) + lambda_lowerbound;
  vector[N * save_scores] eta; // Declare
  matrix[K * marginalize, J * marginalize] multi_normal_mu;
  array[K * marginalize] matrix[J * marginalize, J * marginalize] multi_normal_sigma;

  if(save_scores) { // Compute latent score.
    eta = eta_mean[group] + eta_z .* eta_sd[group];
  } else {
    // TODO: Optimization needed: Make marg_expect_uni return an array of vectors for mu.
    multi_normal_mu = marg_expect_uni(lambda, nu, lambda_random, nu_random, eta_mean);
    multi_normal_sigma = marg_cov_uni(lambda, resid_log, lambda_random, resid_random, eta_sd);
  }
}

model {
  // Declarations
  matrix[N * save_scores, J * save_scores] xhat;
  matrix[N * save_scores, J * save_scores] s_loghat;
  vector[total] hm_hat;
  if(save_scores) { // Construct xhat and s_loghat.
    xhat = rep_matrix(nu, N) + eta*lambda; // Fixed effects
    s_loghat = rep_matrix(resid_log, N); // Fixed effects
    xhat += nu_random[group] + rep_matrix(eta, J) .* lambda_random[group];
    s_loghat += resid_random[group];
  }

  if(use_hmre) { // Construct dependent hierarchical inclusion model.
    hm_hat = exp(hm_tau + hm_param[hm_param_index] + hm_item[hm_item_index] + hm_lambda);
  } else { // Construct independent hierarchical inclusion model.
    // In the future, may want to change this to use, say, a constant, or integrated, prior.
    hm_hat = exp(hm_lambda); // For random-effect models; they use a non-dependent regularization.
  }

  // Priors
  /* Lambda:
     p(Lambda) = p(exp(lam) + lam_lb)|dl/dL|
     = p(exp(lam) + lam_lb)exp(lam)
     We want p(Lambda) to be N^+(0, 1).
     The following is the same as lambda ~ N^+(0, 1)T[lowerbound,] with jacobian
   */
  target += normal_lpdf(lambda | 0, 1) - normal_lccdf(lambda_lowerbound | 0, 1) + sum(lambda_log);
  //// (Log) Residual SD
  resid_log ~ std_normal();
  //// Intercepts
  nu ~ std_normal();

  //// Random effects
  to_vector(random_z) ~ std_normal();
  random_L ~ lkj_corr_cholesky(1);

  //// Latent scores
  if(save_scores) { // Std. latent scores
    eta_z ~ std_normal();
  }

  //// Latent Means/SDs
  if(sum_coding) { // Priors on K-1 eta mean/sds
    eta_mean_s ~ std_normal();
    eta_sd_s ~ std_normal();
  } else { // Else: In Random_z; need priors on RE-SDs
    eta_random_sigma ~ std_normal();
  }

  //// Hierarchical Inclusion Params.
  hm_tau ~ normal(hmre_mu, hmre_scale); // Always included for now; not always used.
  if(use_hmre) { // HMRE priors
    hm_param ~ normal(hmre_mu, hmre_scale);
    hm_item ~ normal(hmre_mu, hmre_scale);
    hm_lambda ~ normal(hmre_mu, hmre_scale); // Always inculded for now.
  } else {
    hm_lambda ~ normal(4 * hmre_mu, 2 * hmre_scale); // Implied: LN(exp(hm_lambda) | 4mu, sqrt(4 scale))
  }

  //// Random Effect SDs
  random_sigma ~ normal(0, hm_hat);

  // Likelihood
  if(!prior_only) {
    if(save_scores) { // Use vectorized data optimization.
      x_vector ~ normal(to_vector(xhat), to_vector(exp(s_loghat)));
    } else { // Marginalized: Use sorted-by-group vectorized multi normal
      for(k in 1:K) {
	x_sorted[x_sorted_indices[k,1]:x_sorted_indices[k,2]] ~ multi_normal(multi_normal_mu[k], multi_normal_sigma[k]);
      }
    }
  }

}

generated quantities {
  matrix[total, total] RE_cor = L_to_cor(random_L)[1:total, 1:total]; // Excluding eta_{mean,sd} cors.
}

