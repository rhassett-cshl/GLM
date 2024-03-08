## function to run glm of epigenomic model, l_s:learning_size, t:tolerance
do_glm <- function(gb, l_s, t){
  
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
