#!/usr/bin/env Rscript

#### load argparse ####
library(argparse)

#### parse command line arguments ####

parser <- ArgumentParser(prog = "./estimate_kmer_kappa_simulation.R",
                         description = "Estimate coefficients kappa on simulated data of k-mer model")

parser$add_argument("-i", type = "character", required = TRUE, 
                    metavar = "InputFile",
                    help = "input training set of simulated data by SimPol")

parser$add_argument("-k_m", type = "character", required = TRUE, 
                    metavar = "Inputkmer",
                    help = "input training set of kmer sparse matrix")

parser$add_argument("-k_t", type = "character", required = TRUE, 
                    metavar = "InputkmerTypes",
                    help = "types for all kmers (AATTA, AATTT...)")

parser$add_argument("-m", type = "character", required = TRUE,
                    metavar = "Model", 
                    help = "select a model (epigenomic, kmer)")

parser$add_argument("-l_s", type = "double", required = TRUE,
                    metavar = "LearningSize",
                    help = "learning size for gradient ascent")

parser$add_argument("-t", type = "double", required = TRUE,
                    metavar = "Tolerance", 
                    help = "tolerance for gradient ascent")

parser$add_argument("-p", type = "double", required = TRUE,
                    metavar = "penalty",
                    help = "hyperparameter for L1 penalty, log10(nu)")

parser$add_argument("-o", type = "character", required = FALSE, 
                    metavar = "OutputDir", default=".",
                    help = "directory for saving the output of kappa")

args <- parser$parse_args()

## get individual parameters
gb_in = args$i
kmer_in = args$k_m
kmer_type_in = args$k_t
model = args$m
l_s = args$l_s
t = args$t
p = 10^(args$p) ## don't forget to take 10^(p)
out_dir = args$o

## read in gb
gb = readRDS(gb_in)
kmer = readRDS(kmer_in)
kmer_type = readRDS(kmer_type_in)


# #### test #####
# gb = readRDS('D:/pub_glm/data/simKmer_gbrc_trainOne.Rdata')
# kmer = readRDS('D:/pub_glm/data/simKmer_kmerMT_trainOne.Rdata')
# kmer_type = readRDS('D:/pub_glm/data/allmer_types.Rdata')
# model = "kmer"
# l_s = 1e-4
# t = 10
# p = 10^-5.8
# out_dir = 'D:/pub_glm/out'

#### set output path ####
kappa_out = file.path(out_dir, paste0('sim_', model, '_kappa.csv'))

# produce all combinations of lambda1
all_grid = tidyr::expand_grid(lambda1 = p)

#### load packages ####
suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(Matrix.utils)
})

#### select the calculation functions based on models #####
if(model == "epigenomic"){
  source("main_ep_glm_functions.R")
  source("main_ep_ga_functions.R")
}else if(model == "kmer"){
  source("main_kmer_glm_functions.R")
  source("main_ep_ga_functions.R")
}


######### call function and run #########
## print on-screen short messages 
print("finish reading in dataframe")
print("start fitting ...")
print("fitting process is recorded in .log file")

# ## record the gradient ascent steps into a log file
log_out = file.path(out_dir, paste0('sim_', model, '_kappa.log'))
sink(log_out, append = FALSE, split = FALSE) # suppress on-screen output

## call gradient ascent function 
cur_k = do_kmer_glm(grid = all_grid, Yji = kmer, gb = gb, l_s = l_s, t = t)
sink() ## stop and exit sink(), revert output back to the console

## print on-screen short messages 
print("finish fitting")

#### changes to the the final output (.csv) file
## change col names
col_n <- kmer_type[1:4^5]

## transpose the matrix so it's more readable in .csv format 
kappa_tb = cur_k %>% 
  dplyr::mutate(feature = col_n) %>% 
  dplyr::relocate(feature, .before = kappa)

## save kappa values into a .csv file 
write.table(kappa_tb, kappa_out, sep = ",", row.names = F, col.names = T)
