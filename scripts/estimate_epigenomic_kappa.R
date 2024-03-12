#!/usr/bin/env Rscript

#### load argparse ####
library(argparse)

#### load packages ####
suppressPackageStartupMessages({
  library(tidyverse)
})

#### parse command line arguments ####

parser <- ArgumentParser(prog = "./estimate_epigenomic_kappa.R",
                         description = "Estimate coefficients kappa of epigenomic model")

parser$add_argument("-i", type = "character", required = TRUE, 
                    metavar = "InputFile",
                    help = "input training set of gene coordination, read counts and covariates")

parser$add_argument("-c", type = "character", default = "simulation", 
                    metavar = "Category",
                    help = "category of input (name of cell line or simulation) [default \"%(default)s\"]")

parser$add_argument("-l_s", type = "double", default = 1e-7,
                    metavar = "LearningSize",
                    help = "learning size for gradient ascent [default %(default)s]")

parser$add_argument("-t", type = "double", default = 1e-1,
                    metavar = "Tolerance",
                    help = "tolerance for gradient ascent [default %(default)s]")

parser$add_argument("-o", type = "character", default=".", 
                    metavar = "OutputDir", 
                    help = "directory for saving the output of kappa")

args <- parser$parse_args()

## get individual parameters
gb_in = args$i
cate = args$c
l_s = args$l_s
t = args$t
out_dir = args$o

## read in gb
gb = readRDS(gb_in)

# #### test #####
# gb = readRDS('D:/pub_glm/data/k562_samp_epft_norm_train_1.Rdata')
# cate = 'k562'
# l_s = 1e-7
# t = 1000
# out_dir = 'D:/pub_glm/out'

# gb = readRDS('D:/pub_glm/data/simEP_gbrcGaussian_trainOne.Rdata')
# cate = "simulation"
# l_s = 1e-6
# t = 1e-2
# out_dir = 'D:/pub_glm/out'

#### set output path based on category ####
# name pattern
if(cate == "simulation"){
  out_pattern = "sim_epigenomic_kappa"
  
  #### keep all histone marks in simulation analysis
  gb = gb
}else if(cate != "simulation"){
  out_pattern = paste0(cate, "_epigenomic_kappa")
  
  #### remove the tss-specific promoter in real data analysis
  tss_hs <- c('h3k27ac', 'h3k4me2', 'h3k4me3', 'h3k9ac')
  gb = gb %>%
    dplyr::select(-dplyr::any_of(tss_hs))
}

# path of log file and .csv output
log_out = file.path(out_dir, paste0(out_pattern, '.log'))
kappa_out = file.path(out_dir, paste0(out_pattern, '.csv'))


#### select the calculation functions based on epigenomic models #####
source("main_ep_glm_functions.R")
source("main_ga_functions.R")


######### call function and run #########
## print on-screen short messages 
print("finish reading in dataframe")
print("start fitting ...")
print("fitting process is recorded in .log file")

## record the gradient ascent steps into a log file
sink(log_out, append = FALSE, split = FALSE) # suppress on-screen output

## call gradient ascent function 
cur_k = do_ep_glm(gb, l_s = l_s, t = t)
sink() ## stop and exit sink(), revert output back to the console

## print on-screen short messages 
print("finish fitting")

#### changes to the the final output (.csv) file
kappa_tb = cur_k %>% ## transpose the matrix so it's more readable in .csv format 
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "feature", values_to = "kappa") %>% 
  dplyr::mutate(feature = dplyr::case_when( ## change col names
    startsWith(feature, 'h') ~ gsub("([hk])", '\\U\\1\\U\\2', feature, perl=TRUE),
    feature == 'ctcf' ~ 'CTCF',
    feature == 'sj5' ~ '5\' spl',
    feature == 'sj3' ~ '3\' spl',
    feature == 'dms' ~ 'stem-loop',
    feature == 'wgbs' ~ 'DNAm',
    feature == 'rpts' ~ 'low complx'
  )) %>% 
  dplyr::arrange(kappa)

## save kappa values into a .csv file 
write.table(kappa_tb, kappa_out, sep = ",", row.names = F, col.names = T)
