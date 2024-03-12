######## This script is to generate the single input of kappa in either ########
######## epigenomic model or combined model, in script predict_zeta.R   ########
library(tidyverse)

############################## INPUT PATH ##############################
## main path
root_dir = '/grid/siepel/home_norepl/liliu/projects/unified_model'

## the path of compare cell analysis
comp_dir = paste0(root_dir, '/compare_cell/shared_sites')

## path of the ep+allmer features
all_epAllmerK_df_in = paste0(comp_dir, '/allCell_ori_epAllmerKappaTb.Rdata')


############################## OUTPUT PATH ##############################
out_dir = '/grid/siepel/home_norepl/liliu/projects/pub_glm/data'
epKappa_fourCell_out = paste0(out_dir, '/epKappa_fourCell.csv')
epAllmerKappa_fourCell_out = paste0(out_dir, '/epAllmerKappa_fourCell.csv')

########### prepare epigenomic coefficients for all cell lines ###########
## K562&MCF7: internal TSS removed; CD14+&HELA: internal TSS not removed
## path of the shared gene (sg) analysis
sg_dir = paste0(root_dir, '/compare_cell/shared_genes')
## path of the shared sites (sg) analysis
ss_dir = paste0(root_dir, '/compare_cell/shared_sites')

## path of kappa of the shared gene analysis
sg_oriK_df_in = paste0(sg_dir, '/capCut_oriEpKappaTb.Rdata')
## path of kappa of the shared sites analysis
ss_oriK_df_in = paste0(ss_dir, '/allCell_oriEpKappaTb.Rdata')

## read in two sets of kappa
sg_oriK_df = readRDS(sg_oriK_df_in)
ss_oriK_df = readRDS(ss_oriK_df_in)


## merge two sets of kappa and change name
all_oriK_df = ss_oriK_df %>%
  dplyr::filter(type == 'hela_original' | type == 'cd14_original') %>%
  dplyr::bind_rows(sg_oriK_df)

## make file changes only for .csv format 
all_oriK_csv = all_oriK_df %>% 
  dplyr::select(-sd_coef) %>% 
  dplyr::rename(kappa = mean_coef)

## save .csv file for better illustration 
write.table(all_oriK_csv, epKappa_fourCell_out, sep = ",", row.names = F, col.names = T)



######### ep allmer combined model kappa ############
all_epAllmerK_df = readRDS(all_epAllmerK_df_in)

## make file changes only for .csv format 
all_epAllmerK_csv = all_epAllmerK_df %>% 
  dplyr::rename(kappa = coef)

## save .csv file for better illustration 
write.table(all_epAllmerK_csv, epAllmerKappa_fourCell_out, sep = ",", row.names = F, col.names = T)



## feature





