## Examples

### Estimate rates based on simulated data

```
usage: ./estimate_kappa_simulation.R [-h] -i InputFile -m Model -l_s
                                     LearningSize -t Tolerance [-o OutputDir]

Estimate coefficients kappa on simulated data

options:
  -h, --help         show this help message and exit
  -i InputFile       input training set of simulated data by SimPol
  -m Model           select a model (epigenomic, kmer)
  -l_s LearningSize  learning size for gradient ascent
  -t Tolerance       tolerance for gradient ascent
  -o OutputDir       directory for saving the output of kappa
```

The input data is produced by [SimPol](https://github.com/CshlSiepelLab/SimPol), a simulator we developed for simulating the dynamics of RNA Polymerase (RNAP) on DNA template. 
A .csv file is generated by SimPol for each simulated transcription unit (TU). 100 TUs are simulated and 80% serves as training data.

All data used in these exmaples can be accessed [here](https://usegalaxy.org/u/lingjie_liu/h/glm-test-data). 

```
./estimate_kappa_simulation.R -i ../data/simEP_gbrcGaussian_trainOne.Rdata -m epigenomic -l_s 1e-6 -t 1e-2 -o ../out
```
