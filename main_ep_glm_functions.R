##### This script is to store the main functions that is needed by GLM #######
##### under the epigenomics setting ##########################################

## calculate lambda, SBj,TBj
calculate_onceCompute <- function(gb){
  # calculate lambda: gb is binned into windows, so length l should be the number 
  # of windows per gb
  gene_rc <- gb %>% 
    dplyr::mutate(ensembl_gene_id = as_factor(ensembl_gene_id)) %>% ## maintain order of gene
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarize(score = sum(score))
  gene_length <- gb %>% 
    dplyr::mutate(ensembl_gene_id = as_factor(ensembl_gene_id)) %>% ## maintain order of gene
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarize(bin_num = dplyr::n())
  lambda <- sum(gene_rc$score)/sum(gene_length$bin_num)

  #Yji contains gene_id, xji and features
  Yji <- gb %>% 
    dplyr::select(5:last_col())
  
  # keep the gene order from gb
  gene_order = gene_rc$ensembl_gene_id
  
  #calculation of TBj
  TBj <- Yji %>% 
    dplyr::mutate(across(c(3: last_col()), ~ .x*score)) %>% 
    dplyr::select(-score) %>% 
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::summarise(across(.cols = everything(), sum)) %>% 
    dplyr::arrange(match(ensembl_gene_id, gene_order)) ### KEEP GENE ORDER!!!
  
  return(list(lambda = lambda,
              SBj = gene_rc,
              TBj = TBj,
              gene_order = gene_order))
}

# calculate e(-k.yji)
calculate_expNdot <- function(k, Yji){
  power <- Yji %>%
    dplyr::select(3:last_col()) %>% 
    as.matrix(.) %*% k %>% 
    as.vector()
  
  expNdot <- Yji %>%
    dplyr::select(ensembl_gene_id) %>% 
    dplyr::mutate(exp_power = exp(-1 * power))
  
  return(expNdot)
}

# calculate UBj
calculate_UBj <- function(expNdot, gene_order){
  UBj <- expNdot %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarise(UBj = sum(exp_power)) %>% 
    dplyr::arrange(match(ensembl_gene_id, gene_order)) ### KEEP GENE ORDER!!!
  
  return(UBj)
}

# calculate alphaj
calculate_alphaj <- function(lambda, SBj, UBj){
  alphaj <- SBj %>% 
    dplyr::inner_join(UBj, by = 'ensembl_gene_id') %>% 
    dplyr::mutate(alpha = score / (lambda * UBj)) %>% 
    dplyr::select(-score, -UBj)
  return(alphaj)
}

# calculate VBj 
calculate_VBj <- function(expNdot, Yji, gene_order){
  VBj <- Yji %>% 
    dplyr::select(-score) %>% 
    dplyr::mutate(across(where(is.numeric), ~.x*expNdot$exp_power)) %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::summarise(across(.cols = everything(), sum)) %>% 
    dplyr::arrange(match(ensembl_gene_id, gene_order)) ### KEEP GENE ORDER!!!
  
  return(VBj)
}

# calculate simplified likelihood
calculate_likelihood <- function(SBj, k, TBj, UBj){
  item1 <- (-1)*SBj$score*log(UBj$UBj)
  
  item2 <- TBj %>% 
    dplyr::select(2:last_col()) %>% 
    as.matrix(.) %*% k %>% 
    as.vector()
  
  likelihood <- sum(item1-item2)
  
  return(likelihood)
}

# calculate gradient 
calculate_gradient <- function(lambda, alphaj, VBj, TBj){
  item1 <- VBj %>%
    dplyr::mutate(across(where(is.numeric), ~.x * lambda * alphaj$alpha)) %>%
    dplyr::select(-ensembl_gene_id)
  
  TBj_number <- TBj %>% dplyr::select(-ensembl_gene_id)
  
  gradient <- colSums(item1 - TBj_number)
  
  return(gradient)
}
