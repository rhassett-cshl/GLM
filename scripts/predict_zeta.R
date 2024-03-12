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

parser$add_argument("-m", type = "character", default = "epigenomic", 
                    metavar = "Category",
                    help = "select a model (epigenomic, combined) [default \"%(default)s\"]")

parser$add_argument("-k_m", type = "character", default = NULL, 
                    metavar = "KmerMatrix",
                    help = "input of k-mer sparse matrix in the predicted regions")

parser$add_argument("-k_k", type = "character", default = NULL, 
                    metavar = "KmerKappa",
                    help = "estimated kappa values for kmer features")

parser$add_argument("-e_m", type = "character", required = TRUE, 
                    metavar = "EpigenomicMatrix",
                    help = "estimated kappa values for epigenomic features")

parser$add_argument("-e_k", type = "character", required = TRUE, 
                    metavar = "EpigenomicKappa",
                    help = "estimated kappa values for epigenomic features")

parser$add_argument("-o", type = "character", required = FALSE, 
                    metavar = "OutputDir", default=".",
                    help = "directory for saving the predicted nucleotide-specific zeta")

args <- parser$parse_args()

## get individual parameters
gb_in = args$g
cate = args$c
model = args$m
kmer_in = args$k_m
kmer_kappa_in = args$k_k
ep_in = args$e_m
ep_kappa_in = args$e_k
out_dir = args$o

## read in
gb = readRDS(gb_in)
kmer = readRDS(kmer_in)
kmer_kappa = readRDS(kmer_kappa_in)
ep = readRDS(ep_in)
ep_kappa = readRDS(ep_kappa_in)


##### test ######




