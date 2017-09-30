### DeltaNeTS version 1.0

The R version of DeltaNeTS package (__v.1.0__) have been successfully tested on __R (>= 3.4.1)__  platforms. Please refer to the DeltaNeTS manuscript for more detailed information about the algorithm. Any questions regarding DeltaNeTS usage can be addressed to heeju.noh@chem.ethz.ch or to rudi.gunawan@chem.ethz.ch.


#### Installation instruction:

To install `deltanets` directly from github repository, `devtools` is required. 

1. Install and load `devtools` package.
2. Install the package, using `devtools::install_github("CABSEL/DeltaNeTS/DeltaNeTS_R/deltanets_v1.0_R")`. Your package is inatalled in R library directory.
3. Load the package, using `library(deltanets)`.


#### Example data in deltanets package

deltanets package includes gene expression data generated from GeneNetWeaver [1] software.:


`lfc`: time-series expression (log2FC) profiles for 1000 genes from 114 signle-gene knockout (KO) experiments. For each experiment, gene expressions at 6 time points (0, 20, 40, 60, 80, 100) are included. Please refere to GNW1 data in DeltaNeTS manuscript.

`glist`: The list of gene symbols corresponding to the rows in the log2FC data

`tobject`: The table of sample descriptions including time points (if in time-series) and group indices (same index for the same experiment)

`grn`: Transcription factor (TF)-gene interactions of the 10000-gene size network used for generating gene expression data.

#### Preparation for DeltaNeTS inputs

DeltaNeTS requires log2FC data and slope matrix (if data are time-series). 

The slope matrix can be calculated using `generateSlope` function.

```{r warning=FALSE,eval=FALSE,echo=TRUE}
tp <- tobject$Time ## a vector of time points of the samples in the matrix lfc
group <- tobject$Group ## a vector of indices of the grouped samples
slope <- generateSlope(lfc = lfc, tp = tp, group = group)
```

#### Implementing DeltaNeTS
DeltaNeTS infers the network model A and perturbation matrix P, using lasso regression with 10-fold cross validation (default).

```{r warning=FALSE,eval=FALSE,echo=TRUE}
result <- deltanets(lfc = lfc, slope = slope)
result$P ## perturbation matrix
result$A ## inferred network
```

Parallel computation is also available in deltanets. The following example shows parallel computation (`par`=`TRUE`) using 4 cores.  
```{r warning=FALSE,eval=FALSE,echo=TRUE}
result <- deltanets(lfc = lfc, slope = slope, par = TRUE, numCores = 4)
```

#### Ranking the genes based on perturbation scores
The greater magnitude of the perturbation scores by deltanets implies the higher perturbation caused in the experiment. The rank of genes or ranked gene lists can be obtained, using rankp function in deltanets.
* In the example data set, tobject includes target gene information for each sample, and can be compared with the rank predictions.

```{r warning=FALSE,eval=FALSE,echo=TRUE}
rankMat = rankp(result,glist)
head(rankMat$rankOfGenes[1:10,1:12],n=12)
```
#### __REFERENCES__:
[1]	Schaffter,T.,Marbach,D., and Floreano,D. (2011) GeneNetWeaver: In silico benchmark generation and performance profiling of network inference methods. Bioinformatics, 27(16), 2263–2270.
