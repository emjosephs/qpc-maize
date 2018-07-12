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

Other required packages are *mass*, *dplyr*, *qvalue*, [*LaCroixColoR*](https://github.com/johannesbjork/LaCroixColoR) and *viridis*

## Q<sub>pc</sub>
**Qpc-gwaspanel.Rmd** has scripts for running Q<sub>pc</sub> in a panel of 240 maize lines with associated phenotypes.

**Qpc-ames.Rmd** has scripts for running Q<sub>pc</sub> on polygenic scores to test for polygenic adaptation in the Ames panel.

**Qpc-euro.Rmd** has scripts for running Q<sub>pc</sub> on polygenic scores to test for polygenic adaptation in European landraces.

## Simulations
**Simulations-traitqpc.Rmd** has code for running simulations on trait Q<sub>pc</sub> in the GWAS panel. 

**Simulations-polygenicqpc.Rmd** has code for running simulations on polygenic Q<sub>pc</sub> in the ames panel and European landraces. 

I have not included the files needed to run the simulations because they are quite large, but I can provide them upon request or put them up on figshare.

## Data
Here is a description of the files in data

**263-gwas-results** has results from GWAS that correspond to Q<sub>pc</sub> analysis on polygenic scores in European landraces. These files only have significant hits, the full association mapping results are on [figshare](https://figshare.com/articles/263_GWAS/6807473)

**281-gwas-results** has results from GWAS that correspond to Q<sub>pc</sub> analysis on polygenic scores in the Ames panel. The entire output file from GEMMA is available on [figshare](https://figshare.com/articles/281_GWAS/6807461).

The kinship matrix used in Q<sub>pc</sub> on traits in the GWAS Panel is in **All\_240E.nomaf.nomissing.K**.

The kinship matrix used in Q<sub>pc</sub> on polygenic scores in the Ames panel is in **ames.281E.K.rda** and the matrix for the Ames panel only is in **amesOnly.E.K.rda**

The kinship matrix used in Q<sub>pc</sub> on on polygenic scores in the European landraces is in **euro.282.E.rda** and the landrace-only matrix (after eigendecomposition) is in **euro-only-eigen.rda**.

**240.names**, and **merged263Landraces.names** have the names of individuals in the kinship matrices **All\_240E.nomaf.nomissing.K** and **euro.282.E.rda** respectively. **blup.names** ahs the names of the traits in the order they are indexed. 

**FlintGarciaTableS1.csv** has data for the GWAS panel from [Flint-Garcia et al. 2005](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-313X.2005.02591.x) and **eurolandraceinfo.csv** has info about the European landraces.

## Figures
Here are scripts to generate figures used in the preprint


## Etc
Thanks to Hilary Parker for the guide on making a personal R package: https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
