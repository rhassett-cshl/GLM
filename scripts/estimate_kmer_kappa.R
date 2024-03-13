#### load argparse ####
library(argparse)

#### load packages ####
suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(Matrix.utils)
})

#### parse command line arguments ####

parser <- ArgumentParser(prog = "Rscript ./estimate_kmer_kappa.R",
                         description = "Estimate coefficients kappa of k-mer model")

parser$add_argument("-i", type = "character", required = TRUE, 
                    metavar = "InputFile",
                    help = "input training set of gene coordination, read counts")

parser$add_argument("-k_m", type = "character", required = TRUE, 
                    metavar = "Inputkmer",
                    help = "input training set of k-mer sparse matrix")

parser$add_argument("-k_t", type = "character", required = TRUE, 
                    metavar = "InputkmerTypes",
                    help = "types for all kmers (AATTA, AATTT...)")

parser$add_argument("-c", type = "character", default = "simulation", 
                    metavar = "Category",
                    help = "category of input (name of cell line or simulation) [default \"%(default)s\"]")

parser$add_argument("-l_s", type = "double",  default = 1e-7,
                    metavar = "LearningSize",
                    help = "learning size for gradient ascent [default %(default)s]")

parser$add_argument("-t", type = "double",  default = 1e-2,
                    metavar = "Tolerance", 
                    help = "tolerance for gradient ascent [default %(default)s]")

parser$add_argument("-p", type = "double", default = -6,
                    metavar = "penalty",
                    help = "hyperparameter for L1 penalty, log10(nu) [default %(default)s]")

parser$add_argument("-o", type = "character", required = FALSE, 
                    metavar = "OutputDir", default=".",
                    help = "directory for saving the output of kappa")

args <- parser$parse_args()

## get individual parameters
gb_in = args$i
kmer_in = args$k_m
kmer_type_in = args$k_t
cate = args$c
l_s = args$l_s
t = args$t
p = 10^(args$p) ## don't forget to take 10^(p)
out_dir = args$o

## read in gb
gb = readRDS(gb_in)
kmer = readRDS(kmer_in)
kmer_type = readRDS(kmer_type_in)


#### test #####
# gb = readRDS('D:/pub_glm/data/simKmer_gbrc_trainAll.Rdata')
# kmer = readRDS('D:/pub_glm/data/simKmer_kmerMT_trainAll.Rdata')
# kmer_type = readRDS('D:/pub_glm/data/allmer_types.Rdata')
# cate = "simulation"
# l_s = 1e-4
# t = 10
# p = 10^-5.8
# out_dir = 'D:/pub_glm/out'

# gb = readRDS('D:/pub_glm/data/k562_subset1_gbrc_train.RData')
# kmer = readRDS('D:/pub_glm/data/k562_subset1_allmerMT_train.RData')
# kmer_type = readRDS('D:/pub_glm/data/allmer_types.Rdata')
# cate = "k562"
# l_s = 1e-7
# t = 2000 # optimized 1e-1
# p = 10^-3.6 ##optimized 10^-3.6
# out_dir = 'D:/pub_glm/out'


#### set output path based on category ####
# name pattern and kmer used in the current analysis
if(cate == "simulation"){
  out_pattern = "sim_kmer_kappa"
  kmer_inUse = kmer_type[c(1 : 4^5)] # simulation uses only 5-mers (k = 5)
}else if(cate != "simulation"){
  out_pattern = paste0(cate, "_kmer_kappa")
  kmer_inUse = kmer_type # simulation uses all k-mers (k <= 5)
}

# path of log file and .csv output
log_out = file.path(out_dir, paste0(out_pattern, '.log'))
kappa_out = file.path(out_dir, paste0(out_pattern, '.csv'))

#### select the calculation functions based on models #####
source("main_kmer_glm_functions.R")
source("main_ga_functions.R")


######### call function and run #########
## print on-screen short messages 
print("finish reading in dataframe")
print("start fitting ...")
print("fitting process is recorded in .log file")

## record the gradient ascent steps into a log file
sink(log_out, append = FALSE, split = FALSE) # suppress on-screen output

# call gradient ascent function
if(cate == "simulation"){
  cur_k = do_sim_kmer_glm(lambda1 = p, Yji = kmer, gb = gb, l_s = l_s, t = t)
}else if(cate != "simulation"){
  cur_k = do_allmer_glm(lambda1 = p, Yji = kmer, gb = gb, l_s = l_s, t = t)
}
sink() ## stop and exit sink(), revert output back to the console

## print on-screen short messages 
print("finish fitting")

#### changes to the the final output (.csv) file
## change col names
col_n <- kmer_inUse

## transpose the matrix so it's more readable in .csv format 
kappa_tb = cur_k %>% 
  dplyr::mutate(feature = col_n) %>% 
  dplyr::relocate(feature, .before = kappa)

## save kappa values into a .csv file 
write.table(kappa_tb, kappa_out, sep = ",", row.names = F, col.names = T)
