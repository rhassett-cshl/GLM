############## This script is to do the GA of epigenomic features ############
library(tidyverse)
library(dplyr)
library(GenomicRanges)

root_dir = '/grid/siepel/home_norepl/liliu/projects/unified_model'

# path of sampled gb in bins with read counts and all epigenomic features
gb_ft_in = paste0(root_dir,
                  '/data/PROseq-RNA-K562-dukler-1_samp_epft_norm.Rdata')
  
  
# read in 
gb <- readRDS(gb_ft_in)

# use the main function 
source(paste0(root_dir, '/glm/main_ep_glm_functions.R'))


#Yji contains gene_id, xji and features
Yji <- gb %>% 
  dplyr::select(5:last_col())

# calculation of once computeted variables: lambda & SBj & gene_order & TBj
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



###################### original GA ####################################
learning_size = 1e-7 #previously 0.0001

increase_cut <- 1 #previously 0.01

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

  if(change_step == T & (L-L0)<increase_cut){
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

## see final g and k
print(g)
print(k)

## info that is needed to be saved
result = list(total_l = total_l,
              total_g = total_g,
              total_k = total_k)

## save 
result_out = paste0(root_dir, 
                    '/data/ep_result/PROseq-RNA-K562-dukler-1_samp_grocap2kb.RData')
saveRDS(result, result_out)
