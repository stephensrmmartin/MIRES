functions {
  vector break_that_stick(vector stick_slices) {
    int K = num_elements(stick_slices) + 1;
    vector[K] pi;
    pi[1] = stick_slices[1];
    for(k in 2:(K - 1)) {
      pi[k] = stick_slices[k] * prod(1 - stick_slices[1:(k - 1)]);
    }
    // pi[k] = 1 - sum(pi[1:(K - 1)]);
    pi[K] = prod(1 - stick_slices[1:(K - 1)]);

    return(pi);
  }
}

data {
  int N; // Number of values for DP density estimation
  int K; // Maximum number of DP clusters
  vector[N] y; // Data
}

transformed data {
  
}

parameters {
  // DP params
  real<lower = 0> alpha; // Dirichlet-beta prior value.
  vector<lower = 0, upper = 1>[K - 1] stick_slices;

  // Cluster params
  vector[K] mu;
  vector<lower = 0>[K] sigma;
}

transformed parameters {
  vector<lower = 0, upper = 1>[K] pi = sort_desc(break_that_stick(stick_slices));
  
}

model {
  /* real sigma = 1.0; // Just use a N(mu, 1) for DP components for now. */
  vector[K] log_pi = log(pi);

  // Priors
  mu ~ normal(0, 3);
  sigma ~ normal(0, 2);
  alpha ~ gamma(2, 2);
  stick_slices ~ beta(1, alpha);

  // Model
  /*
    p(y|...) = sum_k pi_k p(y | theta_k)
    log(p(y|...)) = log(sum_k pi_k p(y | theta_k))
                  = log(sum_k exp(log(pi_k) + log_p(y | theta_k)))
   */ 
  for(n in 1:N) {
    vector[K] lp_y = log_pi;
    for(k in 1:K) {
      lp_y[k] += normal_lpdf(y[n] | mu[k], sigma[k]);
    }
    target += log_sum_exp(lp_y);
  }
  
}

generated quantities {
  
}
