

    /*
      Get the marginal covariance matrices for groups from model parameters.
      This lets us forgo computing the \sum_k^K n_k eta values, and instead compute K covariance matrices.
      Y_ik | k = nu_k + eta_ik * Lambda_k + epsilon_ik, epsilon_ik |k ~ N(0, Sigma_k)
      E(Y_ik | k) = nu_k + E(eta_ik | k) * Lambda_k
      Cov(Y_ik | k) = Lambda'_k V(eta_ik | k) Lambda_k + diag(exp(2 * resid_log_k))
    */
array[] matrix marg_cov_uni(
		      row_vector lambda,
                      row_vector resid_log,
                      matrix lambda_random,
                      matrix resid_random,
		      vector eta_sd
		 ) {
  int J = cols(lambda);
  int K = rows(lambda_random);
  array[K] matrix[J,J] cov_out;
  matrix[K, J] lambda_k = (rep_matrix(lambda, K) + lambda_random) .* rep_matrix(eta_sd, J);
  matrix[K, J] resid_k = exp(2 * (rep_matrix(resid_log, K) + resid_random));
  for(k in 1:K) {
    /* row_vector[J] lambda_k = (lambda + lambda_random[k]); */
    /* matrix[1, J] lambda_k = to_matrix((lambda + lambda_random[k]) * eta_sd[k]); */
    // Optimize this; really sure there's an optimize fn for these ops.
    // Might be able to analytically get cholesky factor from this...
    // No (row)vector crossproducts in Stan.
    // Would be faster in likelihood if instead of n in 1:N, it was k in 1:K, with sorted data.
    /* cov_out[k] = lambda_k' * eta_sd[k]^2 * lambda_k + diag_matrix(exp(2 * (resid_log + resid_random[k]))'); */
    /* cov_out[k] = crossprod(lambda_k) + diag_matrix(exp(2 * (resid_log + resid_random[k]))'); */
    cov_out[k] = crossprod(to_matrix(lambda_k[k,])) + diag_matrix(resid_k[k,]');
  }
  return(cov_out);
}

matrix marg_expect_uni(
			 row_vector lambda,
			 row_vector nu,
			 matrix lambda_random,
			 matrix nu_random,
			 vector eta_mean
			 ) {
  int J = cols(lambda);
  int K = rows(lambda_random);
  matrix[K, J] exp_out;
  exp_out = rep_matrix(nu, K) + nu_random + rep_matrix(eta_mean, J) .* (rep_matrix(lambda, K) + lambda_random);
  /* for(k in 1:K) { */
  /*   exp_out[k] = (nu + nu_random[k]) + eta_mean[k] * (lambda + lambda_random[k]); */
  /* } */
  return(exp_out);
}

array[,] int sort_data_by_group_indices(array[] int group) {
  int K = max(group);
  int N = num_elements(group);
  array[N] int group_sorted = sort_asc(group);
  array[K,2] int out;
  array[K] int n_k = rep_array(0, K);
  int index = 1;
  for(n in 1:N) { // Count observations of each group.
    n_k[group[n]] += 1;
  }
  for(k in 1:K) {
    out[k,1] = index;
    out[k,2] = index + n_k[k] - 1;
    index += n_k[k];
  }
  return(out);
}

array[] row_vector sort_data_by_group(matrix x, array[] int group) {
  int K = max(group);
  int N = rows(x);
  int J = cols(x);
  array[N] row_vector[J] out;
  /* matrix[N, J] out; */
  array[N] int group_ordered = sort_indices_asc(group);

  /* int index = 1; */
  /* for(k in 1:K) { */
  /*   for(n in 1:N) { */
  /*     if(group[n] == k) { */
  /* 	group_ordered[index] = n; */
  /* 	index += 1; */
  /*     } */
  /*   } // N */
  /* } // K */
  /* out = x[group_ordered]; */
  for(n in 1:N) {
    out[n] = x[group_ordered[n]];
  }
  return(out);
}

