
  /*
    Convert z-matrix, scales, and (cholesky) corr matrix to REs
    @param matrix[K, num] of MVN(0, I) zs; where K is number of groups; num is number of REs.
    @param vector[num] sds; scales of REs
    @param matrix[num, num]; (cholesky) corr matrix of REs.
    @return matrix[K, num]; matrix of REs.
   */ 
  matrix z_to_random(matrix z, vector sds, matrix L) {
    int K = rows(z);
    int num = cols(z);
    matrix[K, num] re = z * (diag_pre_multiply(sds, L))';

    return(re);
    
  }

  /*
    Goal: lambda[jk] > 0, for all k.
    Solution: lambda[jk] = lambda[j] + u_lambda[jk] >= 0; met when lambda[j] + min(u_lambda[jk])_k  >= 0
    Therefore: lambda[j] >= -min(u_lambda[jk])_k
    @param matrix[K, J] lambda_random;
    @return row_vector[J] lambda lower bounds.
   */
  row_vector compute_lambda_lowerbounds(matrix lambda_random){
    int J = cols(lambda_random);
    row_vector[J] lambda_lowerbound;

    for(j in 1:J) {
      lambda_lowerbound[j] = -min(lambda_random[, j]);
    }

    return(lambda_lowerbound);
    
  }

  matrix L_to_cor(matrix L) {
    return(multiply_lower_tri_self_transpose(L));
  }
