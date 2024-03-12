## function to run glm of epigenomic model, l_s:learning_size, t:tolerance
do_ep_glm <- function(gb, l_s, t){
  
  # covariate
  Yji <- gb %>% 
    dplyr::select(5:last_col())
  
  # calculation of once computed variables: lambda & SBj & gene_order & TBj
  once_compute = calculate_onceCompute(gb)
  lambda = once_compute$lambda
  SBj = once_compute$SBj
  TBj = once_compute$TBj
  gene_order = once_compute$gene_order
  
  ##### initialize all values
  k = rep(0.0, ncol(Yji) - 2)
  
  expNdot <- calculate_expNdot(k, Yji)
  
  UBj = calculate_UBj(expNdot, gene_order)
  
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  
  VBj = calculate_VBj(expNdot, Yji, gene_order)
  
  L0 = calculate_likelihood(SBj, k, TBj, UBj)
  
  g = calculate_gradient(lambda, alphaj, VBj, TBj)
  
  
  print("finish initializing")

  
  ###################### original GA ####################################
  learning_size = l_s #1e-6
  
  increase_cut <- t #1e-2
  
  go_next <- T
  
  total_l = c(L0)
  total_g <- c(g)
  total_k <- c(k)
  
  
  while(go_next == T){
    
    # Propose next kappa
    k1 = g*learning_size + k
    #initialize change_step for each iteration
    change_step <- F
    
    ## calculation for new log likelihood
    expNdot <- calculate_expNdot(k1, Yji)
    UBj = calculate_UBj(expNdot, gene_order)
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj(expNdot, Yji, gene_order)
    
    L = calculate_likelihood(SBj, k1, TBj, UBj)
    print("Proposal Likelihood:")
    print(L)
    
    ## compare old likelihood and new likelihood
    while(L < L0){
      print("Decrease learning_size")
      change_step <- T
      learning_size = learning_size/2
      print("learning_size:")
      print(learning_size)
      
      # Propose next kappa
      k1 = g*learning_size + k
      
      expNdot <- calculate_expNdot(k1, Yji)
      UBj = calculate_UBj(expNdot, gene_order)
      alphaj = calculate_alphaj(lambda, SBj, UBj)
      VBj = calculate_VBj(expNdot, Yji, gene_order)
      
      L = calculate_likelihood(SBj, k1, TBj, UBj)
      
      print("Likelihood increment:")
      print(L-L0)
      
    }
    
    if((L - L0) < increase_cut){ # to accelerate the ga
      print("Stop!")
      go_next <- F
    }
    
    
    k = k1
    L0 = L
    
    g = calculate_gradient(lambda, alphaj, VBj, TBj)
    
    #record log likelihood
    total_l = c(total_l, L)
    total_k = c(total_k, k)
    total_g = c(total_g, g)
    
  }
  
  cur_k = as.data.frame(t(k)) %>% 
    tibble::as_tibble()
  
  return(cur_k)
}


## function to run glm of kmer model only, specifically for simulation
do_sim_kmer_glm <- function(lambda1, Yji, gb, l_s, t){
  
  # calculation of once computed variables: lambda & SBj & gene_order & TBj
  once_compute = calculate_onceCompute(gb, Yji)
  lambda = once_compute$lambda
  SBj = once_compute$SBj
  gene_order = once_compute$gene_order
  TBj = once_compute$TBj
  
  # initialize k, and other items
  k = rep(0.0, ncol(Yji))
  n = nrow(Yji)
  
  expNdot <- calculate_expNdot(k, Yji)
  UBj = calculate_UBj(expNdot, gene_order)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  VBj = calculate_VBj(expNdot, Yji, gene_order)
  L0  = calculate_lasso_likelihood(SBj, k, TBj, UBj, lambda1, n)
  g = calculate_lasso_gradient(lambda, alphaj, VBj, TBj, lambda1, k, n)
  
  
  ##### GA #####
  learning_size = l_s
  
  # This value sets a bound of parameter precision
  tolerance <- t
  
  go_next <- T
  
  total_l = c(L0)
  total_g <- c(g)
  total_k <- c(k)
  
  
  lastL_decrease <- F
  
  while(go_next == T){
    
    # Propose next kappa
    k1 = g*learning_size + k
    #initialize change_step for each iteration
    change_step <- F
    
    
    ## calculation for new log likelihood
    expNdot <- calculate_expNdot(k1, Yji)
    UBj = calculate_UBj(expNdot, gene_order)
    
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj(expNdot, Yji, gene_order)
    
    L = calculate_lasso_likelihood(SBj, k1, TBj, UBj, lambda1, n)
    print("Proposal Likelihood:")
    print(L)
    
    
    ## compare old likelihood and new likelihood
    if (lastL_decrease){
      while(L < L0){
        print("Decrease learning_size")
        change_step <- T
        learning_size = learning_size/2
        print("learning_size:")
        print(learning_size)
        
        # Propose next kappa
        k1 = g*learning_size + k
        
        expNdot <- calculate_expNdot(k1, Yji)
        UBj = calculate_UBj(expNdot, gene_order)
        alphaj = calculate_alphaj(lambda, SBj, UBj)
        VBj = calculate_VBj(expNdot, Yji, gene_order)
        
        L = calculate_lasso_likelihood(SBj, k1, TBj, UBj, lambda1, n)
        
        print("Likelihood increment:")
        print(L-L0)
        
      }
      lastL_decrease = F
    }
    if (!lastL_decrease & L < L0) {
      lastL_decrease = T
    }
    
    if (L > L0){
      lastL_decrease = F
    }
    
    if((L-L0) < t & (L-L0) > 0 ){
      print("Stop!")
      go_next <- F
    }
    
    
    k = k1
    L0 = L
    
    g = calculate_lasso_gradient(lambda, alphaj, VBj, TBj, lambda1, k, n)
    
    #record log likelihood
    total_l = c(total_l, L)
  }
  
  cur_k = tibble::tibble(kappa = k)
  
  return(cur_k)
}

## function to run glm of kmer model only, adding linear transformation
do_allmer_glm <- function(lambda1, Yji, gb, l_s, t){
  
  # calculate c1 and c2 for linear transformation 
  norm_item = calculate_norm_item(Yji)
  c1 = norm_item$c1
  c2 = norm_item$c2
  
  # calculation of once computeted variables: lambda & SBj & gene_order & TBj
  once_compute_norm = calculate_onceCompute_norm(gb, Yji, c1, c2)
  lambda = once_compute_norm$lambda
  SBj = once_compute_norm$SBj
  gene_order = once_compute_norm$gene_order
  TBj = once_compute_norm$TBj
  
  # initialize k, and other items
  k = rep(0.01, ncol(Yji))
  n = nrow(Yji)                             
  
  # initiate iterating variables
  expNdot = calculate_expNdot_norm(k, Yji, c1, c2)
  UBj = calculate_UBj(expNdot, gene_order)
  alphaj = calculate_alphaj(lambda, SBj, UBj)
  VBj = calculate_VBj_norm(expNdot, Yji, gene_order, c1, c2, UBj)
  L0  = calculate_lasso_likelihood(SBj, k, TBj, UBj, lambda1, n)
  g = calculate_lasso_gradient(lambda, alphaj, VBj, TBj, lambda1, k, n)
  
  
  ##### GA #####
  learning_size = l_s
  
  # This value sets a bound of parameter precision
  tolerance <- t
  
  go_next <- T
  
  total_l = c(L0)
  total_g <- c(g)
  total_k <- c(k)
  
  lastL_decrease <- F
  
  while(go_next == T){
    
    # Propose next kappa
    k1 = g*learning_size + k
    #initialize change_step for each iteration
    change_step <- F
    
    
    ## calculation for new log likelihood
    expNdot <- calculate_expNdot_norm(k1, Yji, c1,c2)
    UBj = calculate_UBj(expNdot, gene_order)
    alphaj = calculate_alphaj(lambda, SBj, UBj)
    VBj = calculate_VBj_norm(expNdot, Yji, gene_order, c1, c2, UBj)
    
    L = calculate_lasso_likelihood(SBj, k1, TBj, UBj, lambda1, n)
    print("Proposal Likelihood:")
    print(L)
    
    
    ## compare old likelihood and new likelihood
    if (lastL_decrease){
      while(L < L0){
        print("Decrease learning_size")
        change_step <- T
        learning_size = learning_size/2
        print("learning_size:")
        print(learning_size)
        
        # Propose next kappa
        k1 = g*learning_size + k
        
        expNdot <- calculate_expNdot_norm(k1, Yji, c1,c2)
        UBj = calculate_UBj(expNdot, gene_order)
        alphaj = calculate_alphaj(lambda, SBj, UBj)
        VBj = calculate_VBj_norm(expNdot, Yji, gene_order, c1, c2, UBj)
        
        L = calculate_lasso_likelihood(SBj, k1, TBj, UBj, lambda1, n)
        
        print("Likelihood increment:")
        print(L-L0)

      }
      lastL_decrease = F
    }
    if (!lastL_decrease & L < L0) {
      lastL_decrease = T
    }
    
    if (L > L0){
      lastL_decrease = F
    }
    
    if((L - L0) < tolerance & (L-L0) > 0 ){
      print("Stop!")
      go_next <- F
    }
    
    
    k = k1
    L0 = L
    
    g = calculate_lasso_gradient(lambda, alphaj, VBj, TBj, lambda1, k, n)
    
    #record log likelihood
    total_l = c(total_l, L)
  }
  
  cur_k = tibble::tibble(kappa = k)
  
  return(cur_k)
}




## function to run glm of combined model (epigenomic + kmer model)
do_allmer_ep_glm = function(grid, y1, y2, gb, l_s, t){
  
  # grid = all_grid
  # l_s = 1e-6
  # t = 1
  
  # current lasso hyperparameter
  lambda1 = grid[1] %>% as.numeric() # convert to numeric, to avoid producing name attached number
  
  # total data point
  n = nrow(gb)
  
  # initialize k 
  k = rep(0.01, ncol(y1) + ncol(y2))
  
  # y1 linear transformation
  norm_item = calculate_norm_item(y1)
  c1 = norm_item$c1
  c2 = norm_item$c2
  
  # y1 based once compute quantities
  once_compute_1 = calculate_onceCompute_norm(gb, y1, c1, c2)
  lambda_1 = once_compute_1$lambda
  SBj_1 = once_compute_1$SBj
  gene_order_1 = once_compute_1$gene_order
  TBj_1 = once_compute_1$TBj
  
  # print("done: TBj_1")
  
  # y2 based once compute quantities
  once_compute_2 = calculate_onceCompute(gb, y2)
  TBj_2 = once_compute_2$TBj
  
  # print("done: TBj_2")
  
  # y1-y2 combined once compute 
  TBj_1_2 = Matrix::cbind2(TBj_1, TBj_2)
  
  # print("done: TBj_1_2")
  
  ## remove useless object and release memory
  rm(TBj_1)
  rm(TBj_2)
  gc()
  
  allmer_ep_quant <- calculate_allmer_ep_lasso(y1, y2, lambda1, k, n, 
                                               gene_order_1, SBj_1, c1, c2, 
                                               lambda_1, TBj_1_2)
  
  L0 = allmer_ep_quant$likelihood
  g = allmer_ep_quant$gradient
  
  ##### GA #####
  learning_size = l_s
  
  # This value sets a bound of parameter precision
  tolerance <- t
  
  go_next <- T
  
  total_l = c(L0)
  total_g <- c(g)
  
  lastL_decrease <- F
  
  while(go_next == T){
    
    # Propose next kappa
    k1 = g*learning_size + k
    #initialize change_step for each iteration
    change_step <- F
    
    
    ## calculation for new log likelihood
    allmer_ep_quant <- calculate_allmer_ep_lasso(y1, y2, lambda1, k1, n, 
                                                 gene_order_1, SBj_1, c1, c2, 
                                                 lambda_1, TBj_1_2)
    
    L = allmer_ep_quant$likelihood
    
    print("Proposal Likelihood:")
    print(L)
    
    
    ## compare old likelihood and new likelihood
    if (lastL_decrease){
      while(L < L0){
        print("Decrease learning_size")
        change_step <- T
        learning_size = learning_size/2
        print("learning_size:")
        print(learning_size)
        
        # Propose next kappa
        k1 = g*learning_size + k
        
        allmer_ep_quant <- calculate_allmer_ep_lasso(y1, y2, lambda1, k1, n, 
                                                     gene_order_1, SBj_1, c1, c2, 
                                                     lambda_1, TBj_1_2)
        
        print("done: quantity")
        
        L = allmer_ep_quant$likelihood
        
        print("Likelihood increment:")
        print(L-L0)
        
      }
      lastL_decrease = F
    }
    if (!lastL_decrease & L < L0) {
      lastL_decrease = T
    }
    
    if (L > L0){
      lastL_decrease = F
    }
    
    if((L-L0) < tolerance & (L-L0) > 0 ){
      print("Stop!")
      go_next <- F
    }
    
    
    k = k1
    L0 = L
    
    g = allmer_ep_quant$gradient
    
    #record log likelihood
    total_l = c(total_l, L)
    total_g = c(total_g, g)
    
  }
  
  cur_k = tibble::tibble(kappa = k)
  
  return(cur_k)

}

