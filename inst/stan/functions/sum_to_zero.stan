
  /*
    Compute K eta_means from K-1 eta_mean_s
    Assumes sum-to-zero (STZ) constraint.
   */
  vector eta_means_stz(vector eta_mean_s) {
    int K = num_elements(eta_mean_s) + 1;
    vector[K] eta_mean;
    eta_mean[2:K] = eta_mean_s;
    eta_mean[1] = -sum(eta_mean_s);

    return(eta_mean);
    
  }

  vector eta_sds_pto(vector eta_sd_s) {
    int K = num_elements(eta_sd_s) + 1;
    vector[K] eta_sd;
    eta_sd[2:K] = eta_sd_s;
    eta_sd[1] = exp(-sum(log(eta_sd_s)));

    return(eta_sd);
    
  }
