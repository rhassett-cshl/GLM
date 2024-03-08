######### generate a single file for the input of simulated EP model   #########
######### Output files: either the total training set of the 10 simulations ####
######### Or just one round of simulation, for easiest use             #########
######### The single input file for script: estimate_kappa_simulation.R ########
library(tidyverse)

## version of introducing Gaussian noise 
version = 'v2'


############################## input #######################################
input_dir = paste0('C:/Users/ling/Desktop/simulation/', version)

# path of the read counts of the simulator
sim_rc_in = paste0(input_dir, '/dataset/sim_rc_train.Rdata')

## path of the features for the simulator 
sim_gb_in = paste0(input_dir, '/dataset/sim_gb_train.Rdata')

############################## output #######################################
output_dir =  'D:/pub_glm'
all_gb_out = paste0(output_dir, 
                    '/data/simEP_gbrcGaussian_trainAll.Rdata')
one_gb_out = paste0(output_dir, 
                    '/data/simEP_gbrcGaussian_trainOne.Rdata')

####################### start of : merge gb and rc ###############################
## read in 
sim_rc = readRDS(sim_rc_in)
sim_gb = readRDS(sim_gb_in)

## merge gb and rc, so it will have columns: gb_score_epft
all_gb = sim_gb %>%
  dplyr::select(-zeta) %>%
  dplyr::mutate(score = sim_rc$`1`) %>% # just select one Poisson sampled reads
  dplyr::relocate(score, .after = ensembl_gene_id) %>%
  dplyr::mutate(ensembl_gene_id = as.character(ensembl_gene_id))

## select one round of simulation 
i = 1 # define nth round of simulation 
sim_time = 10 # give total number of simulation 
sim_batch = length(unique(sim_rc$ensembl_gene_id)) / sim_time
gene_batch_id = ((i - 1) * sim_batch + 1) : (i * sim_batch) %>% 
  as.character()

one_gb = all_gb %>% 
  dplyr::filter(ensembl_gene_id %in% gene_batch_id)

####################### end of : merge gb and rc ###############################

## save the merged file 
saveRDS(all_gb, all_gb_out)
saveRDS(one_gb, one_gb_out)

