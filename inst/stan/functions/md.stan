  /*
    Generates item indices for each param, but for multidimensional models.
    MD models require more care, because of the possibility of cross-loadings.
    When you have J indicators, there aren't J loadings; it depends on the spec.
    Nus and sigmas constitute 2*J indices, but not loadings.
    Currently, loadings are organized in (row-major) order of lambdaEst.
   */
  int[] gen_item_indices_md(int J, int F, int[] J_f, int[,] F_ind) {
    int total_lambda = sum(J_f);
    int resid_nu[2*J];
    int lambda[total_lambda];
    int lambda_resid_nu[total_lambda + 2*J];

    // Resid and nu
    int base = 0;
    for(j in 1:(2*J)) {
      resid_nu[j] = j - base;
      if(j - base == J) {
	base += J;
      }
    }

    // Loadings
    base = 1;
    for(f in 1:F) {
      lambda[base:(base - 1 + J_f[f])] = F_ind[f, 1:J_f[f]];
      base += J_f[f];
    }

    // Combine
    lambda_resid_nu = append_array(lambda, resid_nu);

    return(lambda_resid_nu);
  }

  int[] gen_param_indices_md(int J, int[] J_f) {
    int resid[J] = rep_array(2, J);
    int nu[J] = rep_array(3, J);
    int lambda_total = sum(J_f);
    int lambda[lambda_total] = rep_array(1, lambda_total);
    int lambda_resid_nu[lambda_total + 2*J] = append_array(lambda, append_array(resid, nu));

    return(lambda_resid_nu);
  }

  /*
    Generates 3x2 array, to ease indexing later.
    E.g.,
    [1, 7] // Loading locations are vec[1:7]
    [8, 12] // Resid locations are vec[8:12]
    [13, 17] // Nu locations are vec[13:27]
   */
  int[,] gen_lamResNu_bounds(int J, int[] J_f){
    int out[3, 2];
    int lambda_total = sum(J_f);
    out[1] = {1, lambda_total};
    out[2] = {lambda_total + 1, lambda_total + J};
    out[3] = {lambda_total + J + 1, lambda_total + 2*J};

    return(out);
  }

  matrix lambda_mat(int[] J_f, int[,] F_ind, row_vector lambdaEst) {
    int F_J[2] = dims(F_ind);  
    int F = F_J[1];
    int J = F_J[2];
    int tot = sum(J_f);
    matrix[F, J] lambda = rep_matrix(rep_vector(0.0, F), J);

    int count = 1;
    for(f in 1:F) {
      for(jj in F_ind[f, 1:J_f[f]]) {
	lambda[f, jj] = lambdaEst[count];
	count += 1;
      }
    }

    return(lambda);
  }

/* 
   Convex combination of two L cor matrices.
   Converts the Ls to cor matrices, convexly combines them, returns a cholesky-decomposed cor matrix.
   TODO: Check whether we need to convert to Cors first, or whether convex combinations of cholesky-decomposed R's yields a valid cholesky-decomposed R.
   @param matrix L1: A lower-diagonal matrix.
   @param matrix L2: A lower-diagonal matrix.
   @param vector[1] weight_L2: The weight for LL2.
   @return matrix (1 - w)LL1 + w(LL2)
 */
matrix convex_combine_Ls(matrix L1, matrix L2, vector weight_L2) {
  int R = rows(L1);
  int C = cols(L1);
  matrix[R, C] LL1 = multiply_lower_tri_self_transpose(L1);
  matrix[R, C] LL2 = multiply_lower_tri_self_transpose(L2);
  matrix[R, C] outLL = (1 - weight_L2[1]) * LL1 + (weight_L2[1]) * LL2;
  
  matrix[R, C] outL = cholesky_decompose(outLL);

  return(outL);
  
}
