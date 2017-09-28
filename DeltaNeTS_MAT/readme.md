
### DeltaNeTS version 1.0
The MATLAB subroutines in the DeltaNeTS package (__v.1.0__) have been successfully tested on __MATLAB® 2014b to 2016b__  platforms. Please refer to the [DeltaNeTS manuscript](http://www.sciencedirect.com/science/article/pii/S2405896316328154) for more detailed information about the algorithm. Any questions regarding DeltaNeTS usage can be addressed to heeju.noh@chem.ethz.ch or to rudi.gunawan@chem.ethz.ch.

#
> #### Installation instruction: 

1.	Unzip the package [___DeltaNeTS_v1.0_MAT.zip___](https://github.com/CABSEL/DeltaNeTS/blob/master/DeltaNeTS_MAT/DeltaNeTS_v1.0_MAT.zip) (1.94 MB) to a preferred folder.
2.	Download the MATLAB version of [GLMNET package](http://web.stanford.edu/~hastie/glmnet_matlab/download.html), and unzip the GLMNET package under a new subfolder in DeltaNeTS.  
3.	Set the current working directory to DeltaNeTS in MATLAB. 
4.	Add the path for GLMNET package.

#
 >  #### The DeltaNeTS package includes the following:

<br />

__1. example_data__

A subfolder in DeltaNeTS package, containing example data for testing DeltaNeTS

* ___log2FC_GNW1_114KO-EXP_6tpooints.txt___: log2FC of gene expression data (`lfc`) generated from GeneNetWeaver__[1]__ software. `lfc` consists of time-series expression profiles for 1000 genes from 114 signle-gene knockout (KO) experiments. For each experiment, gene expressions at 6 time points (0, 20, 40, 60, 80, 100) are included. Please refere to GNW1 data in DeltaNeTS manuscript.
* ___GList.txt___: The list of gene symbols corresponding to the rows in the log2FC data
* ___table_of_samples.txt___: The table of sample descriptions including time points (if in time-series) and group indices (same index for the same KO experiment)

<br /><br />
__2. findiff.m__

This function implements a 2nd order accurate finite difference for calculating slopes (time-derivatives) of the log2FC data using three time points. The function is used in _generateSlope_ function below. 

<br /><br />
__3. 	generateSlope.m__

The function for calculating slope matrix from log2FC data. If more than two time points are available for a given drug/compound treatment, then a 2nd order accurate finite difference approximation is used for calculating the slopes. If only two time points are available, then a linear slope between the two time points is used.

```
slope = generateSlope( lfc, tp, group )
```

REQUIRES: **findiff.m**  
<br />
INPUT ARGUMENTS:

* `lfc`:	The matrix of log2FC data. Each row represents a gene and each column represents a sample.
* `tp`:	A vector of time points of the samples in the matrix `lfc`. The length of the vector should be the same as the number of samples (i.e. the number of columns in the matrix `lfc`).
* `group`:	A vector of indices indicating the set of samples from a particular drug/compound treatment. The (time-series) samples from the same drug treatment experiment should have the same unique index. The length of the vector should be the same as the number of samples.
<br />
OUTPUT ARGUMENTS:

`slope`: the slope matrix

<br /><br />
__4.	run_deltanets_example.m__

An example script of running deltanets with example data (GNW1).

<br /><br />
__5.	deltanets.m__

The main function for DeltaNeTS for generating the gene perturbation impacts for each sample. 

```
[P, A] = deltanets( lfc, slope, kfold, par, numCores )
```
<br />
INPUT ARGUMENTS:

* `lfc`:	The matrix of log2FC data. Each row represents a gene and each column represents a sample.
* `slope`:	The slope matrix from log2FC data. This matrix can be obtained using the function _generateSlope_. If the data are not time-series, set slope to an empty matrix (i.e. `slope`=[]). 
* `kfold`:	The number of folds used in the k-fold cross validation.
* `par`:	A Boolean variable _TRUE_ or _FALSE_ indicating whether to use parallel computing. The default is _FALSE_ (no parallel computation).
* `numCores`: The number of CPU cores to be used for parallel computing. This parameter is considered only if par is _TRUE_. The default is 4.
<br />
OUTPUT ARGUMENTS:

* `P`: The matrix of perturbation impacts for each gene. Each row corresponds to a gene following the same order as the one in the log2FC data, while each column corresponds to samples listed in tobject.
* `A`: The _n_ x _n_ matrix of the inferred GRN. The (_i_, _j_)-th element of the matrix corresponds to the coefficient of regulaotry impact from gene _j_ to gene _i_ in the GRN. The rows and columns of the matrix correspond to genes following the same order as in GList.

<br /><br />
#### __REFERENCES__:
[1]	Schaffter,T.,Marbach,D., and Floreano,D. (2011) GeneNetWeaver: In silico benchmark generation and performance profiling of network inference methods. Bioinformatics, 27(16), 2263–2270.

