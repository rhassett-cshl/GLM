#!/usr/bin/env Rscript

#### load argparse ####
library(argparse)

#### load packages ####
suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(Matrix.utils)
})

#### parse command line arguments ####

parser <- ArgumentParser(prog = "./predict_zeta.R",
                         description = "Predict nucleotide-specific elongation rate zeta")

parser$add_argument("-g", type = "character", required = TRUE, 
                    metavar = "GeneBody",
                    help = "input of gene coordinates in the predicted regions")

parser$add_argument("-c", type = "character", default = "k562", 
                    metavar = "Category",
                    help = "category of input (name of cell line) [default \"%(default)s\"]")

parser$add_argument("-r", type = "character", required = TRUE, 
                    metavar = "ReadCount",
                    help = "input of PRO-seq read counts in the predicted regions")

parser$add_argument("-m", type = "character", default = "epigenomic", 
                    metavar = "Model",
                    help = "select a model (epigenomic, combined) [default \"%(default)s\"]")

parser$add_argument("-e_m", type = "character", required = TRUE, 
                    metavar = "EpMatrix",
                    help = "epigenomic covariates in the predicted regions")

parser$add_argument("-k", type = "character", required = TRUE, 
                    metavar = "Kappa",
                    help = "estimated kappa values for epigenomic or combined models")

parser$add_argument("-k_m", type = "character", default = NULL, 
                    metavar = "KmerMatrix",
                    help = "input of k-mer sparse matrix in the predicted regions")

parser$add_argument("-o", type = "character", default=".", 
                    metavar = "OutputDir", 
                    help = "directory for saving the predicted nucleotide-specific zeta")

args <- parser$parse_args()

## get individual parameters
shared_gb_in = args$g
cell = args$c
rc_in = args$r
model = args$m
ep_in = args$e_m
kappa_in = args$k
kmer_in = args$k_m
out_dir = args$o

## read in
shared_gb = readRDS(shared_gb_in)
rc = readRDS(rc_in)
ep = readRDS(ep_in)
kappa = readr::read_csv(kappa_in, col_types = cols())

## read allmer matrix if it's a combined model
if(model == "combined"){
  kmer = readRDS(kmer_in)  
}


##### test ######
# shared_gb = readRDS('../data/shared_gb_twoGenes.Rdata')
# cell = 'k562'
# rc = readRDS('../data/k562_rc_twoGenes.Rdata')
# model = 'epigenomic'
# ep = readRDS('../data/k562_epft_norm_twoGenes.Rdata')
# kappa = readr::read_csv('../data/epKappa_fourCell.csv', col_types = cols())
# out_dir = 'D:/pub_glm/out'

##### test ######
# shared_gb = readRDS('../data/shared_gb_twoGenes.Rdata')
# cell = 'k562'
# rc = readRDS('../data/k562_rc_twoGenes.Rdata')
# model = 'combined'
# ep = readRDS('../data/k562_epft_norm_twoGenes.Rdata')
# kappa = readr::read_csv('../data/epAllmerKappa_fourCell.csv')
# kmer = readRDS('../data/allmerMT_twoGenes.Rdata', col_types = cols())
# out_dir = 'D:/pub_glm/out'

## path of output bigwig, two strands merged into one file 
bw_out = file.path(out_dir, paste0(cell, '_', model, '_predZeta.bw'))

#### introduce main functions for prediction #####
if(model == "epigenomic"){
  source("main_predict_functions.R")
}else if(model == "combined"){
  source("main_predict_functions.R")
  source("main_kmer_glm_functions.R") ## introduce kmer glm calculation functions for combined model
}



######### call function and run #########
## print see
print("finish reading in")
print("start predicting")

## call function to predict elongation rate based on the selected model
if(model == "epigenomic"){
  pred_zeta_grng = pred_ep_zeta(shared_gb = shared_gb, rc = rc, 
                                ep = ep, ep_kappa = kappa)
}else if(model == "combined"){
  pred_zeta_grng = pred_epAllmer_zeta(shared_gb = shared_gb, rc = rc, ep = ep, 
                                      allmer = kmer, epAllmer_kappa = kappa)
}

## print see
print("finish predicting elongation rates")

## save predicted elongation rates to bigwig files
## export bw
rtracklayer::export.bw(pred_zeta_grng, bw_out)

## print see
print("finish saving prediction to one bigwig file")