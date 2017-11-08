# RefBool
RefBool - a reference-based algorithm for booleanizing gene expression data

## Requirements
RefBool was tested with Matlab 2015b and requires the "Statistics and Machine Learning Toolbox" for distribution fitting. No other toolboxes or third-party functions are needed to run RefBool.

## Usage
The core functionality of RefBool is the booleanization or discretization into three states of gene expression data and has been tested so far for human RNA-Seq and microarray data. We describe in the following the workflow to obtain discretized expression values for query genes.

### 1. Determine Threshold Distributions
RefBool requires the existence of a reference-distribution for assessing whether an expression measurement is associated to an active or inactive state of the corresponding gene. These distributions can be obtained by calling **DetermineThresholdDistributions.m** with the following inputs:
1. An m by n matrix (not cell array!) containing m genes (rows) with n samples (columns) each. This matrix contains the RNA-Seq measurements based on which the threshold distributions are determined. The datasets in the *DataAccess.txt* could be used for human RNA-Seq and microarray data.
2. The number of bootstrap samples to be drawn to approximate the parameter distributions of the probability distribution. We empirically determined that 1000 bootstrap samples were enough to obtain close-to-normal parameter distributions but a lower number may be acceptable here.
3. The precision of the thresholds. Since RefBool is based on samples of the real, unknown gene expression distribution, the thresholds can only be determined up to a certain precision. A value of *10^(-x)* indicates a precision up to the *x*-th decimal place. Here, we recommend to use 10^(-3) to avoid spurious results.
4. A goodness of fit criterion based on which the best fitting distribution is determined (as character). One of
	- **NLogL** - Negative log likelihood (not recommended)
	- **BIC** - Bayesian information criterion
	- **AIC** - Akaike information criterion
	- **AICc** - AIC with correction for finite samples sizes
5. A boolean value indicating whether to choose the maximum (TRUE) or minimum (FALSE) observed gene expression value if no extreme point is detected for the lower threshold. This option might be removed in a future update, since it can be shown that the extreme point is left of the lowest observed value. Thus, **minimum value (FALSE) is strongly recommended!**

RefBool then returns a (m by 1) cell array of lower and upper threshold distributions for further analysis. Importantly, the function does not know of gene IDs and the user has to keep track of them! It is, however, guaranteed that the threshold distributions in row i of the output correspond to gene i in the input. **Note: This step has to be performed only once for each gene and experimental method!**

### 2. Calculate p- and q- values for query expression values
Based on the previously obtained threshold distributions, RefBool derives p- and q-values for a query gene expression (pattern) by calling the function **CalculatePValues.m** with the following inputs:
- **1**: A cell array of previously obtained threshold distributions.
- **2**: A matrix of query expression values containing m genes (rows) in n samples (columns). Importantly, the user has to ensure that the query genes are in the same order as the threshold distributions and both have the same length. In particular, query expression values in row *i* and threshold distributions in row *i* have to belong to the same gene.

q-values are calculated at the end of the function based on Benjamini-Hochberg correction. Even though it is possible to call this function with a single query gene expression value, the user should derive the whole query expression matrix and corresponding threshold distributions and invoke the function with them to ensure that multiple testing correction takes effect.

### 3. Booleanize query expression values
Finally, booleanized or discretized query expression values are obtained by calling the function **Booleanize_Bypval.m** with parameters:
- **1-4**: The q-value matrices for lower, upper and intermediate thresholds obtained in the last step.
- **5**: A significance cutoff for lower and upper thresholds (e.g. 0.05)
- **6**: A significance cutoff for the intermediate region (e.g. 0.1)

The output is a single matrix object populated with 0s, 1s and 0.5s corresponding to inactive, active and intermediate states. Note: The discretization for some query expression values might be *NaN*. This happens if the query values are close to the mean of one of the threshold distributions. However, by relaxing the intermediate region cutoff to *1*, all these values will be assignet to the intermediate expression state.

## Runtime
The most time consuming part of RefBool is the computation of the threshold distributions. For the results described in the paper, the derivation for one gene took around 30 to 45 seconds, depending on the best fitting distribution. Therefore, if one wants to discretize genome-wide gene expression measurements, the use of a computing cluster may be necessary. However, this step needs to be performed only once! The runtime of the other functions is much faster and takes only several seconds to minutes.

## Availability of testdata
- Threshold distributions on a genome-wide scale for human RNA-Seq in TPM is available at:
https://webdav-r3lab.uni.lu/public/cbg/RefBool_ReferenceDistributions/
- A wroking test code is provided in the *example* folder
- See further information in the *DataAccess.txt* file 
- We are currently computing threshold distributions for microarray (probeset IDs) expression and will soon provide links for downloading them. 

## Citation
Sascha Jung, Andras Hartmann, Antonio del Sol; RefBool: a reference-based algorithm for discretizing gene expression data. Bioinformatics 2017 Jul 1;33(13):1953-1962. DOI: [10.1093/bioinformatics/btx111](https://doi.org/10.1093/bioinformatics/btx111).
