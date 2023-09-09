
  /*
    Generates [1,...,J, 1, ..., J, 1, ..., J]
    1 = indicator 1
    J = indicator J
   */ 
  array[] int gen_item_indices(int J) {
    array[3*J] int hm_item_index;

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
  array[] int gen_param_indices(int J) {
    array[3*J] int hm_param_index;

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
  array[,] int gen_lamResNu_indices(int J) {
    array[3, J] int lamResNu_indices;

    for(p in 1:3) {
      for (j in 1:J) {
	lamResNu_indices[p, j] = j + J*(p -1);
      }
    }

    return(lamResNu_indices);
  }
