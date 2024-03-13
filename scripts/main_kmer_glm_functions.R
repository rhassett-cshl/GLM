##### This script is to store the main functions that is needed by GLM #######
##### under the kmer setting #################################################

## calculate lambda, SBj, gene_order, TBj
calculate_onceCompute <- function(gb, Yji){
  gene_rc <- gb %>% 
    dplyr::mutate(ensembl_gene_id = forcats::as_factor(ensembl_gene_id)) %>% ## maintain order of gene
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarize(score = sum(score))
  gene_length <- gb %>% 
    dplyr::mutate(ensembl_gene_id = forcats::as_factor(ensembl_gene_id)) %>% ## maintain order of gene
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarize(bin_num = dplyr::n())
  
  lambda <- sum(gene_rc$score)/sum(gene_length$bin_num)
  
  SBj <- gene_rc
  
  gene_order <- gb$ensembl_gene_id %>%
    match(., unique(.))

  TBj <- (Yji *gb$score) %>%
    Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum')

  return(list(lambda = lambda,
              SBj = gene_rc,
              gene_order = gene_order,
              TBj = TBj))
}

## calculate normalized items, standardized Yji' = c1 * Yji - c2
calculate_norm_item <- function(Yji){
  # colSums(Yji) * nrow(Yji) - colSums(Yji)^2 are negative for A/T/G/C features
  c1 = nrow(Yji) / sqrt(abs(colSums(Yji) * nrow(Yji) - colSums(Yji)^2))
  c2 = colSums(Yji) / sqrt(abs(colSums(Yji) * nrow(Yji) - colSums(Yji)^2))
  
  ## !!! when one kmer is not existing in current features, colSums(Yji) = 0, 
  ## which makes c1 and c2 to be infinity. So, when colSums(Yji) = 0, c1 = c2 = 0
  c1[is.infinite(c1)] = 1
  c2[is.na(c2)] = 1
  
  return(list(c1 = c1,
              c2 = c2))
}

## calculate lambda, SBj, gene_order, TBj(normalized)
calculate_onceCompute_norm <- function(gb, Yji, c1, c2){
  gene_rc <- gb %>% 
    dplyr::mutate(ensembl_gene_id = forcats::as_factor(ensembl_gene_id)) %>% ## maintain order of gene
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarize(score = sum(score))
  gene_length <- gb %>% 
    dplyr::mutate(ensembl_gene_id = forcats::as_factor(ensembl_gene_id)) %>% ## maintain order of gene
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarize(bin_num = dplyr::n())
  
  lambda <- sum(gene_rc$score)/sum(gene_length$bin_num)
  
  SBj <- gene_rc
  
  gene_order <- gb$ensembl_gene_id %>%
    match(., unique(.))
  
  # have the normalized TBj
  TBj_ori <- (Yji *gb$score) %>%
    Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum')

  TBj_norm = TBj_ori %*% as(diag(c1), "sparseMatrix") -
    as(SBj$score %*% t(c2), "sparseMatrix")
  
  return(list(lambda = lambda,
              SBj = gene_rc,
              gene_order = gene_order,
              TBj = TBj_norm))
}


# calculate e(-k.yji)
calculate_expNdot <- function(k, Yji){
  power <- Yji %*% k 
  #head(power)
  expNdot <- exp(-1 * power)
  #head(expNdot)
  return(expNdot)
}

# calculate normalized e(-k.yji)
calculate_expNdot_norm <- function(k, Yji, c1, c2){
  
  expNdot <- exp(sum(c2 * k)) * exp(-1 * Yji %*% (c1 * k))
  
  return(expNdot)
}


# calculate UBj
calculate_UBj <- function(expNdot, gene_order){
  UBj <- expNdot %>%
    Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum')
  
  #head(UBj)
  
  return(UBj)
}

# calculate alphaj
calculate_alphaj <- function(lambda, SBj, UBj){
  alphaj <- SBj$score / (lambda * UBj)
  #head(alphaj)
  
  return(alphaj)
}

# calculate VBj 
calculate_VBj <- function(expNdot, Yji, gene_order){
  VBj <- (Yji * as.vector(expNdot)) %>%
    Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum')
  
  #head(VBj)
  
  return(VBj)
}

# calculate normalized VBj 
calculate_VBj_norm <- function(expNdot, Yji, gene_order, c1, c2, UBj){
  
  item1 <- Yji %*% as(diag(c1), "sparseMatrix") * as.vector(expNdot)
  
  item2 <- UBj * t(replicate(nrow(UBj), c2))
  
  item3 <- item1 %>%
    Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum')

  VBj <- item3 - item2
  
  return(VBj)
}

# calculate simplified likelihood
# calculate_likelihood <- function(SBj, k, TBj, UBj){
#   item1 <- (-1)*SBj$score*log(UBj)
#   #head(item1)
#   
#   item2 <- TBj %*% k 
#   #head(TBj)
#   #head(item2)
#   
#   likelihood <- sum(item1-item2)
#   
#   return(likelihood)
# }

# calculate penalized likelihood, USE FULL LIKELIHOOD EQUATION
calculate_likelihood <- function(SBj, k, TBj, UBj, lambda1, lambda2, n){
  item1 <- (-1)*SBj$score*log(UBj)
  #item1 <- SBj$score*(log(SBj$score) - log(UBj))
  item2 <- TBj %*% k
  
  penalty <- n*(lambda1 * (lambda2 * sum(abs(k)) + (1-lambda2)/2 * sum(k^2))) # takes number of bins, n
  
  likelihood <- sum(item1-item2-SBj$score)
  p_likelihood <- likelihood - penalty
  
  #print(likelihood/(penalty + 0.01)) # prevent when penalty could be 0
  
  return(p_likelihood)
}


# calculate gradient 
calculate_gradient <- function(lambda, alphaj, VBj, TBj, lambda1, lambda2, k, n){
  #head(VBj)
  #head(alphaj)
  item1 <- as.vector(lambda * alphaj) * VBj
  #head(item1)
  #head(TBj) 
  
  ## change gradient due to penalized likelihood
  penalty_g <- function(nth_k, lambda1, lambda2, n){
    if(nth_k == 0){
      p_g = 0 # treat derivative(|k1|, k1=0) = 0
    }else{
      p_g = nth_k/abs(nth_k) * lambda1 * lambda2 + lambda1 * (1 - lambda2) * nth_k
      p_g = p_g * n # take in the number of bins
    }
    return(p_g)
  }
  
  p_gradient <- sapply(k, penalty_g, lambda1, lambda2, n) 
  #head(p_gradient) 
  
  gradient <- colSums(item1 - TBj) - p_gradient
  #head(gradient)
  
  return(gradient)
}

## use simplified likelihood to avoid the case when SBj$score == 0
calculate_lasso_likelihood <- function(SBj, k, TBj, UBj, lambda1, n){
  # item1 <- SBj$score*(log(SBj$score) - log(UBj))
  item1 <- (-1)*SBj$score*log(UBj)
  
  item2 <- TBj %*% k
  
  penalty <- n * lambda1 * sum(abs(k)) # takes number of bins, n
  
  likelihood <- sum(item1-item2-SBj$score)
  p_likelihood <- likelihood - penalty
  
  #print(likelihood/(penalty + 0.01)) # prevent when penalty could be 0
  
  return(p_likelihood)
}

# calculate gradient 
calculate_lasso_gradient <- function(lambda, alphaj, VBj, TBj, lambda1, k, n){
  #head(VBj)
  #head(alphaj)
  item1 <- as.vector(lambda * alphaj) * VBj
  #head(item1)
  #head(TBj) 
  
  ## change gradient due to penalized likelihood
  penalty_g <- function(nth_k, lambda1, n){
    if(nth_k == 0){
      p_g = 0 # treat derivative(|k1|, k1=0) = 0
    }else{
      p_g = nth_k/abs(nth_k) * lambda1
      p_g = p_g * n # take in the number of bins
    }
    return(p_g)
  }
  
  p_gradient <- sapply(k, penalty_g, lambda1, n) 
  #head(p_gradient) 
  
  gradient <- colSums(item1 - TBj) - p_gradient
  #head(gradient)
  
  return(gradient)
}

## calculate the lasso likelihood and gradient of allemr-ep combined matrix 
calculate_allmer_ep_lasso <- function(y1, y2, lambda1, k, n, 
                                      gene_order_1, SBj_1, c1, c2, lambda_1,
                                      TBj_1_2){
  
  # separate k into k1 and k2 for allmer and ep, respectively 
  k1 = k[1:ncol(y1)]
  k2 = k[(ncol(y1)+1):length(k)]
  
  # y1 based calculation
  expNdot1 <- calculate_expNdot_norm(k1, y1, c1, c2)
  
  # UBj_1 = calculate_UBj(expNdot1, gene_order_1)
  # alphaj_1 = calculate_alphaj(lambda_1, SBj_1, UBj_1)
  
  
  # y2 based calculation
  expNdot2 <- calculate_expNdot(k2, y2)
  
  # combined expNdot
  expNdot = expNdot1 * expNdot2
  
  ## remove useless object and release memory
  rm(expNdot1)
  rm(expNdot2)
  gc()
  
  ### test : print
  # print("done: expNdot")
  
  ## y1 y2 combined calculation ##
  UBj_1_2 = calculate_UBj(expNdot, gene_order_1)
  alphaj_1_2 = calculate_alphaj(lambda_1, SBj_1, UBj_1_2)
  
  ## special notice For VBj1 and VBj2, in calculating seperate VBj, 
  ## the expNdot should be expNdot1 * expNdot2
  
  # test: calculate the epigenomic features first 
  VBj_2 = calculate_VBj(expNdot, y2, gene_order_1)

  ### test : print
  # print("done: VBj_2")
  
  VBj_1 = calculate_VBj_norm(expNdot, y1, gene_order_1, c1, c2, UBj_1_2)
  
  ### test : print
  # print("done: VBj_1")
  
  # VBj_2 = calculate_VBj(expNdot, y2, gene_order_1)
  # 
  # ### test : print
  # print("done: VBj_2")
  
  VBj_1_2 = Matrix::cbind2(VBj_1, VBj_2)
  
  ### test : print
  # print("done: VBj_1_2")
  
  
  ## remove useless object and release memory
  rm(VBj_1)
  rm(VBj_2)
  gc()
  
  ## total likelihood and gradient for the allmer-ep combined matrix 
  L0_1_2  = calculate_lasso_likelihood(SBj_1, k, TBj_1_2, UBj_1_2, lambda1, n)
  g_1_2 = calculate_lasso_gradient(lambda_1, alphaj_1_2, VBj_1_2, TBj_1_2, lambda1, k, n)
  
  
  ## return all the quantities the GLM needs
  return(list(likelihood = L0_1_2,
              gradient = g_1_2))
}