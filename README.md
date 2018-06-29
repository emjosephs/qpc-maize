# qpc
scripts for detecting polygenic adaptation in maize using Q<sub>pc</sub>.


## Getting started
all processed data is in the **data** folder, generally in .rda files

**qpctools** contains what's needed for an R package with functions needed to do the analysis. 

This package is installed by running the following command: 
```
Rscript install_qpctools.R
```
You must have the R packages devtools installed already

Other required packages are *mass*, *dplyr*, *qvalue*, and *viridis*. I might make a nice script for installing these at some point.

## Q<sub>pc</sub>
**Qpc-gwaspanel.Rmd** has scripts for running Q<sub>pc</sub> in a panel of 240 maize lines with associated phenotypes.

**Qpc-ames.Rmd** has scripts for running Q<sub>pc</sub> on polygenic scores to test for polygenic adaptation in the Ames panel.

**Qpc-euro.Rmd** has scripts for running Q<sub>pc</sub> on polygenic scores to test for polygenic adaptation in European landraces.

## Simulations
**Simulations-traitqpc.Rmd** has code for running simulations on trait Q<sub>pc</sub> in the GWAS panel. 

**Simulations-polygenicqpc.Rmd** has code for running simulations on polygenic Q<sub>pc</sub> in the ames panel and European landraces. 



## Etc
Thanks to Hilary Parker for the guide on making a personal R package: https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
