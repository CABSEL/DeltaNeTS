

# DeltaNeTS
<img style = "float: right;" src = "https://github.com/CABSEL/DeltaNeTS/blob/master/image_deltanet.png" width="200" height="200" align="right"> 

DeltaNets is an extension of our previous method DeltaNet, which involves a single-step inference of gene regulatory network and gene targets from (time-series) transcriptional profiles. If steady-state gene expression data are applied, it is same as running DeltaNet. 
DeltaNeTS generates a perturbation score for each gene in every sample. For a given sample, the magnitude of a gene perturbation scores reflects the confidence that this gene is directly perturbed (for example, by drug or chemical compounds), while the sign reflects the nature of the perturbation (positive for gene induction, negative for gene repression).

Please refer to [DeltaNeTS](http://www.sciencedirect.com/science/article/pii/S2405896316328154) manuscript for more detailed information.


### Prerequisites:
For MATLAB distribution,
* MATLAB® R2014b or later
* [GLMNET package](http://web.stanford.edu/~hastie/glmnet_matlab/)


### Last Update
* MATLAB package: **version 1.0** (28.09.2017)


### License
Redistribution and use in source and binary forms, with or without modification, are permitted provided agreeing to the *Simplified BSD Style License* (see [more](http://opensource.org/licenses/bsd-license.php)).

[License](copypath) (RTF, 2 KB)


### Download and installation
Please refer to [README for MATLAB](copypath)

### Acknowledgements
This work has been supported by the ETH Zurich Research Grant.
