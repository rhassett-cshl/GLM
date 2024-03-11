#!/usr/bin/env Rscript

#### load argparse ####
library(argparse)

#### parse command line arguments ####

parser <- ArgumentParser(prog = "./estimate_kappa_simulation.R",
                         description = "Estimate coefficients kappa on simulated data")

parser$add_argument("-i", type = "character", required = TRUE, 
                    metavar = "InputFile",
                    help = "input training set of simulated data by SimPol")

parser$add_argument("-m", type = "character", required = TRUE,
                    metavar = "Model", 
                    help = "select a model (epigenomic, kmer)")

parser$add_argument("-l_s", type = "double", required = TRUE,
                    metavar = "LearningSize",
                    help = "learning size for gradient ascent")

parser$add_argument("-t", type = "double", required = TRUE,
                    metavar = "Tolerance", 
                    help = "tolerance for gradient ascent")

parser$add_argument("-o", type = "character", required = FALSE, 
                    metavar = "outputDir", default=".",
                    help = "directory for saving the output of kappa")

args <- parser$parse_args()

## get individual parameters
gb_in = args$i
model = args$m
l_s = args$l_s
t = args$t
out_dir = args$o

## read in gb
gb = readRDS(gb_in)

# #### test #####
# gb = readRDS('D:/pub_glm/data/simEP_gbrcGaussian_trainOne.Rdata')
# model = "epigenomic"
# l_s = 1e-6
# t = 1e-2
# out_dir = 'D:/pub_glm/out'

#### set output path ####
kappa_out = file.path(out_dir, paste0('sim_', model, '_kappa.csv'))

#### load packages ####
suppressPackageStartupMessages({
  library(tidyverse)
})

#### select the calculation functions based on models #####
if(model == "epigenomic"){
  source("main_ep_glm_functions.R")
  source("main_ep_ga_functions.R")
}else if(model == "kmer"){
  source("main_kmer_glm_functions.R")
}

######### call function and run #########
## print on-screen short messages 
print("finish reading in dataframe")
print("start fitting ...")
print("fitting process is recorded in .log file")

## record the gradient ascent steps into a log file
log_out = file.path(out_dir, paste0('sim_', model, '_kappa.log'))
sink(log_out, append = FALSE, split = FALSE) # suppress on-screen output

## call gradient ascent function 
cur_k = do_glm(gb, l_s = l_s, t = t)
sink() ## stop and exit sink(), revert output back to the console

## print on-screen short messages 
print("finish fitting")

#### changes to the the final output (.csv) file
## change col names
col_n <- c('CTCF', 'H3K36me3', 'stem-loop', 'H3K4me2', 'H3K9me3', 'H4K20me1')
colnames(cur_k) <-  col_n

## transpose the matrix so it's more readable in .csv format 
kappa_tb = cur_k %>% 
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "feature", values_to = "kappa")

## save kappa values into a .csv file 
write.table(kappa_tb, kappa_out, sep = ",", row.names = F, col.names = T)
