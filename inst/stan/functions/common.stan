
  /*
    Generates [1,...,J, 1, ..., J, 1, ..., J]
    1 = indicator 1
    J = indicator J
   */ 
  int[] gen_item_indices(int J) {
    int hm_item_index[3*J];

    int base = 0;
    for(j in 1:(3*J)) {
      hm_item_index[j] = j - base;
      if(j - base == J) {
	base += J;
      }
    }

    return(hm_item_index);
    
  }

  /*
    Generates [1, 1, ..., 2, 2, ..., 3, 3, ...]
    1 = lambda
    2 = resid
    3 = nu
   */ 
  int[] gen_param_indices(int J) {
    int hm_param_index[3*J];

    int base = 1;
    for(j in 1:(3*J)) {
      hm_param_index[j] = base;
      if(j % J == 0) {
	base += 1;
      }
    }

    return(hm_param_index);
    
  }

  /*
    Generates 3xJ matrix; e.g., for J = 4 items:
    |1  2  3  4|
    |5  6  7  8|
    |9 10 11 12|
    Intended usecase is for unrolling random effects from 1:(3J) vector into 1:J vectors.
    Instead of using (1:J); ((J+1) : (J*2)); ((J*2 + 1):(J*3)), one can just do:
    lamResNu_indices[1], lamResNu_indices[2], lamResNu_indices[3]
   */
  int[,] gen_lamResNu_indices(int J) {
    int lamResNu_indices[3, J];

    for(p in 1:3) {
      for (j in 1:J) {
	lamResNu_indices[p, j] = j + J*(p -1);
      }
    }

    return(lamResNu_indices);
  }

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
    int J = rows(lambda_random);
    row_vector[J] lambda_lowerbound;

    for(j in 1:J) {
      lambda_lowerbound[j] = -min(lambda_random[, j]);
    }

    return(lambda_lowerbound);
    
  }

  matrix L_to_cor(matrix L) {
    return(multiply_lower_tri_self_transpose(L));
  }
