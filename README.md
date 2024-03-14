---
TITLE: "GLM"
AUTHORS: "Lingjie Liu; Yixin Zhao; Rebecca Hasset; Adam Siepel"
---

# Probabilistic and machine-learning methods for predicting local rates of transcription elongation from nascent RNA sequencing data

## Overview
In this study, we revisit the questions about the determinants of local elongation rates through gene bodies using a fundamentally different statistical modeling approach.  Our method is based on a recently developed "unified model" for nascent RNA sequencing (NRS) data, which describes both the kinetics of Pol II movement on the DNA template and the generation of NRS read counts (see [Siepel 2022](https://www.biorxiv.org/content/10.1101/2021.01.12.426408v1) and [Zhao et al 2023](https://academic.oup.com/nar/article/51/21/e106/7331020)). 

We adapt this model to allow for a continuously variable elongation rate along the genome, using a generalized linear model to capture the relationship between the local rate and nearby genomic features (see [Liu et al 2023](https://www.biorxiv.org/content/10.1101/2023.12.21.572932v1)). By accounting for differences across genes in initiation rate, we are able to efficiently pool information across genes, and extract high-resolution information about relative local elongation rates from steady-state data. 

Here, we provide scripts of the implementation, and demonstrate how quantify both previously known and novel correlates of elongation rate through estimation of coefficients $\kappa$. We then use our models to predict nucleotide-specific elongation rates $\zeta_i$ genome-wide and make our predictions available in a UCSC Genome Browser track.

<p align="center">
  <img src="figures/concept.png" alt="glm" width="1000"/>
</p>

<p align = "left">
	Fig. 1 Conceptual illustration of kinetic model for Pol II movement along DNA template in gene body. At nucleotide site $i$, relative elongation rate $\zeta_i$ is an exponentiated linear function of features $\vec Y_i$ and coefficients $\kappa$. 
	Promoter-proximal pausing and termination are ignored here. 
	Graphical model representation showing unobserved continuous-time Markov chain ($Z_i$) and observed NRS read counts ($X_i$).
	Conceptual illustration showing that differences in average gene-body read depth are explained by the scaled initiation rate $\chi$, while relative read depth is explained by the generalized linear model for local elongation rate $\zeta_i$. 
	Read count $X_i$ is assumed to be Poisson distributed with mean $\frac{\chi}{\zeta_i}$.
</p>



## Dependencies

GLM is implemented in the statistical programming language [R](https://www.r-project.org/), and depends on a couple of packages. One of the easiest ways to install them is via [conda](https://docs.conda.io/en/latest/).

```
conda env create --file pub_glm_env.yaml
```

Once installed, you can activate the environment then run the examples within it,

```
conda activate pub_glm
```

Data for the examples below can be downloaded from [Galaxy](https://usegalaxy.org/u/lingjie_liu/h/glm-test-data), and should be assumed to be stored within the data directory. 



## Examples

### Estimate coefficients $\kappa$ based on epigenomic model

```
usage: Rscript ./estimate_epigenomic_kappa.R [-h] -i InputFile [-c Category]
                                             [-l_s LearningSize]
                                             [-t Tolerance] [-o OutputDir]

Estimate coefficients kappa of epigenomic model

options:
  -h, --help         show this help message and exit
  -i InputFile       input training set of gene coordination, read counts and
                     covariates
  -c Category        category of input (name of cell line or simulation)
                     [default "simulation"]
  -l_s LearningSize  learning size for gradient ascent [default 1e-07]
  -t Tolerance       tolerance for gradient ascent [default 0.1]
  -o OutputDir       directory for saving the output of kappa
```

**Simulation data** <br>
We first tested our modeling approach on data simulated with [SimPol](https://github.com/CshlSiepelLab/SimPol), which tracks the movement of individual RNA polymerases along the DNA templates in thousands of cells under user-defined rates, and then samples NRS read counts in proportion to the steady-state polymerase density. For this study, we extended SimPol to consider synthetic epigenomic correlates, which we tiled along each synthetic DNA template using a block-sampling approach based on real data for CTCF transcription binding sites, four different histone marks, and RNA stem-loops. The gene coordinates of the simulated gene bodies, the synthetic NRS read counts generated by SimPol, and the simulated feature covariates are merged together. 80% of this data is used as the training data for estimating $\kappa$, and consequently serves as the input for `estimate_epigenomic_kappa.R`. <br>

The input file `-i` for the simulation analysis looks like below; `seqnames`, `start`, `end`, `strand`, `ensembl_gene_id` are gene coordinates; `score` is synthetic NRS read counts;  the remaining columns are feature covariates. Notably, each feature is standardized to have a mean of zero and a standard deviation of one.
```
# A tibble: 6 × 12
  seqnames start   end strand ensembl_gene_id score ctcf[,1] h3k36me3[,1] dms[,1] h3k4me2[,1] h3k9me3[,1] h4k20me1[,1]
  <chr>    <int> <int> <chr>  <chr>           <int>    <dbl>        <dbl>   <dbl>       <dbl>       <dbl>        <dbl>
1 1            1     1 +      1                   0  -0.0837       -0.341  -0.308      -0.160      -0.538         1.54
2 1            2     2 +      1                   1  -0.0837       -0.340  -0.308      -0.160      -0.539         1.56
3 1            3     3 +      1                   0  -0.0837       -0.338  -0.308      -0.161      -0.539         1.57
4 1            4     4 +      1                   0  -0.0837       -0.337  -0.308      -0.161      -0.539         1.58
5 1            5     5 +      1                   0  -0.0837       -0.335  -0.308      -0.161      -0.540         1.60
6 1            6     6 +      1                   0  -0.0837       -0.334  -0.308      -0.161      -0.540         1.61
```

The choice of `-l_s` for learning rate and `-t` for tolerance is essential to ensure both the sufficiency (precision of $\kappa$) and efficiency (running time) of likelihood convergence in the gradient ascent algorithm. So, the optimal optimial learning rate and tolerance are specified in the command line below for estimating coefficients $\kappa$ in the simulation.
```
Rscript ./estimate_epigenomic_kappa.R -i ../data/simEP_gbrcGaussian_trainOne.Rdata -c simulation -l_s 1e-6 -t 1e-2 
```

While fitting, the proposed log likelihood for each iteration is recorded in the `sim_epigenomic_kappa.log` file.  
After the completion of fitting, the estimated coefficients $\kappa$ are saved in the `sim_epigenomic_kappa.csv` file, with two columns as follows. Both the .log file and the .csv file can be found in the output directory.

```
"feature","kappa"
"H3K36me3",-0.260447787120152
"stem-loop",-0.0553034654846287
"H4K20me1",-0.0467166979577607
"CTCF",-0.0181587606643657
"H3K4me2",0.0341186759406782
"H3K9me3",0.0412922564113079
```

**Real PRO-seq data** <br>
Having established that our model works well with simulated data, we applied it to real PRO-seq data from K562 cells, for which abundant epigenomic data are available. Currently, the epigenomic model of K562 includes twelve features. Therefore, the input format is similar to the input simulation data shown above but with more feature columns. All features are standardized and subjected to smoothing filters designed to capture broader effects that extend to neighboring nucleotide positions centered on the features. In each replicate, we sampled 2000 genes, with 80% of the data serving as the training data, which serves as the input. One replicate is utilized in the command line and executed as an example (fitting can take hours due to the size of the real data).

```
Rscript ./estimate_epigenomic_kappa.R -i ../data/k562_samp_epft_norm_train_1.Rdata -c k562 -l_s 1e-7 -t 1e-1
```


### Estimate coefficients $\kappa$ based on k-mer model ($k \in \{1, 2, 3, 4, 5\}$)

```
usage: Rscript ./estimate_kmer_kappa.R [-h] -i InputFile -k_m Inputkmer -k_t
                                       InputkmerTypes [-c Category]
                                       [-l_s LearningSize] [-t Tolerance]
                                       [-p penalty] [-o OutputDir]

Estimate coefficients kappa of k-mer model

options:
  -h, --help           show this help message and exit
  -i InputFile         input training set of gene coordination, read counts
  -k_m Inputkmer       input training set of k-mer sparse matrix
  -k_t InputkmerTypes  types for all kmers (AATTA, AATTT...)
  -c Category          category of input (name of cell line or simulation)
                       [default "simulation"]
  -l_s LearningSize    learning size for gradient ascent [default 1e-07]
  -t Tolerance         tolerance for gradient ascent [default 0.01]
  -p penalty           hyperparameter for L1 penalty, log10(nu) [default -6]
  -o OutputDir         directory for saving the output of kappa
```

1. Work with simulated data.
```
Rscript ./estimate_kmer_kappa.R -i ../data/simKmer_gbrc_trainAll.Rdata -k_m ../data/simKmer_kmerMT_trainAll.Rdata -k_t ../data/allmer_types.RData -c simulation -l_s 1e-4 -t 1e-2 -p -5.8
```

2. Work with real data.
```
Rscript ./estimate_kmer_kappa.R -i ../data/k562_subset1_gbrc_train.RData -k_m ../data/k562_subset1_allmerMT_train.RData -k_t ../data/allmer_types.RData -c k562 -l_s 1e-7 -t 1e-1 -p -3.6
```


### Predicted nucleotide-specific elongation rates $\zeta_i$

```
usage: Rscript ./predict_zeta.R [-h] -g GeneBody [-c Category] -r ReadCount
                                [-m Model] -e_m EpMatrix -k Kappa
                                [-k_m KmerMatrix] [-o OutputDir]

Predict nucleotide-specific elongation rate zeta

options:
  -h, --help       show this help message and exit
  -g GeneBody      input of gene coordinates in the predicted regions
  -c Category      category of input (name of cell line) [default "k562"]
  -r ReadCount     input of PRO-seq read counts in the predicted regions
  -m Model         select a model (epigenomic, combined) [default
                   "epigenomic"]
  -e_m EpMatrix    epigenomic covariates in the predicted regions
  -k Kappa         estimated kappa values for epigenomic or combined models
  -k_m KmerMatrix  input of k-mer sparse matrix in the predicted regions
  -o OutputDir     directory for saving the predicted nucleotide-specific zeta
```

1. Predict using the epigenomic model.
```
Rscript ./predict_zeta.R -g ../data/shared_gb_twoGenes.Rdata -c k562 -r ../data/k562_rc_twoGenes.Rdata -m epigenomic -e_m ../data/k562_epft_norm_twoGenes.Rdata -k ../data/epKappa_fourCell.csv
```

2.  Predict using the combined model (epigenomic and k-mer features).
```
Rscript ./predict_zeta.R -g ../data/shared_gb_twoGenes.Rdata -c k562 -r ../data/k562_rc_twoGenes.Rdata -m combined -e_m ../data/k562_epft_norm_twoGenes.Rdata -k ../data/epAllmerKappa_fourCell.csv -k_m ../data/allmerMT_twoGenes.Rdata
```



## Citation
Liu, L., Zhao, Y. & Siepel, DNA-sequence and epigenomic determinants of local rates of transcription elongation. 2023. Preprint at [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.12.21.572932v1)

Zhao, Y., Liu, L., Hassett, R. & Siepel, A. Model-based characterization of the equilibrium dynamics of transcription initiation and promoter-proximal pausing in human cells. Nucleic Acids Res 51, 2023

Siepel, A. A unified probabilistic modeling framework for eukaryotic transcription based on nascent RNA sequencing data. 2021. Preprint at [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.01.12.426408v1)
