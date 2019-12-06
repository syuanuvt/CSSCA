
# ClusterSSCA: A Brief Tutorial (Current Version 0.6.0)
The package includes most essential functions to carry out the Cluster-wise Sparse Simultaneous Component Analysis (CSSCA). To learn more about the technical details of CSSCA, or to cite the methods in your publication, please kindly refer to the following citation. 

## citation
```
S.Yuan, K.De Roover, M.Dufner, J.J.A.Denissen, K. Van Deun (in press). Revealing subgroups that differ in common and distinctive variation in multi-block data: Clusterwise Sparse Simultaneous Component Analysis. 
```

## corresponding address
Issues on the packages could be filed as pull requests (recommended) or sent to s.yuan@uvt.nl. Specific uestions about the code and the method could be addressed at the same email address. 

## Installation
To install the Package, make sure that the R package "devtool" and "remotes" have been installed. Then, use the following code to install the Package ClusterSSCA. 
```{r eval = FALSE}
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
remotes::install_github("syuanuvt/CSSCA")
```

## Data simulation 

### Functions

**CSSCASimulation**: Data simulation following CSSCA models, as described in Yuan et al.(in press), with the options to specify a full set of parameters, including: (1) the number of clusters k, (2) the size of each cluster (allowed to have unbalanced cluster size), (3) the number of common and distinctive components, (4) the number of data blocks L, (5) the number of variables in each data block, (6) the sparsity level of the loading matrices S, (7) the noise level of the data e and (8) the congruence of loading matrices (note that becuase the additional sparseness will be imposed on the loading matrices, the empirical congruence, which can be obtained from a post-hoc manner after the data has been generated), and (8) the proportion of the total variance that is accounted for by the mean level cluster differences m. 

Please note that the calculation of m is after substracting e of the noise from the original data. To give a concrete example, if m=.9 while e=.3, then the total variance could be decomposed to three parts: 30% of the variance pertains to the error structure, 63% to the structural cluster differences and 7% to the mean-level cluster differences. 

### Example
```{r eval = FALSE}
n_cluster <- 3 # number of clusters
mem_cluster <- c(50,50,50) # number of observations in each cluster (of length n_cluster)
n_block <- 2  # number of blocks in the multi-block data
n_com <- 2  # number of common components
n_distinct <- c(1,1) # number of distinctive components in each data block (of length n_block)
n_var <- c(20,10)  # number of variables in each data block (of length n_var)
p_sparse <- 0.5  # level of sparsity
p_noise <- 0.3  # level of noise 
p_combase <- 0.5 # congruence level of the loading matrices 
p_fixzero <- 1 # we assume the positions of zero loadings are identical across different clusters 
mean  <- 0.1 # the proportion of the mean structure in total cluster differences 

CSSCASimulation(n_cluster, mem_cluster, n_block, n_com, n_distinct, n_var, p_sparse,
 p_noise, p_combase, p_fixzero, "both", mean_v)
####

```

## Work flow of analyses on real datasets

(1) Check whether the assumptions of CSSCA () are what you expect and whether the required values of the parameters are available.

As should be done for any other analyses, we strongly advise the users to first examine whether CSSCA  analyses could address the research questions properly and whether the assumptions underlying the CSSCA model are at least reasonable for the application therein. Most importantly, CSSCA assumes that *the same set of subgroups differ in both mean structure and co-variance structure*. Therefore, the methods would *NOT* be suitable if the user only expect mean differences among clusters (in which cases we would suggest the alternative methods (i.e. K-means) to avoid over-fitting) nor could CSSCA detects different clusters with regard to the mean structure and the convariance structure. Another important assumption is that *the same amount of common and distinctive components, which should be supplied by the user, pertain to all clusters*. Therefore, the current method would only be applicable in case the users have solid theoritical expectations towards the underlying structure of the components (we are currently developing new methods to relax such constraints, and will update the package when avilable). Last, as we have briefly discussed in the paper, *the performance of CSSCA deteriorates with increasing level of noise (unsatisfying results would be expected with 40% or higher amount of noise)*. Users are advised to take additional care if their data is prone to high level of noise.

(2) Define the parameters that are needed to run CSSCA: specify *the number of common components that underlie all data blocks and the number of disticntive components that underlie each data block* that underlie all clusters. 

(3) If the number of clusters k and the desired level of sparsity S (in case of CSSCA) are known or fixed by the user, FixedCSSCA could be directly used to estimate CSSCA results.

(4) If the number of clusters k and the desired level of sparsity S (in case of CSSCA) are *NOT* known or fixed by the user, a model selection procedure should first be applied to select the optimal value for each parameter. The model selection procedure is wrapped in the function ModelSelectionCSSCA. FixedCSSCA could then be used to estimate the CSSCA and CPCA solution, given the optimal values of the tuning paramter. Note that in the current version, an important drawback is that the *minimal and maximal values of k and S could not be selected*. As a result, *the current procedure prohibits the selection of one cluster and non-sparseness*, two particular yet important scenerios. We are currently testing alternative model selection scheme that would enable the selection of such extreme values. 

## CSSCA analysis on multi-block datasets (when the tuning parameters are known)

### Functions

**FixedCSSCA**: a wrapper that estimates the CSSCA results with fixed levels of the parameters. 

### Example

In the example, we utilize the same basic parameter setting as used in the simulation study.

```{r eval = FALSE}
resultsCSSCA <- FixedCSSCA(target_data, n_block, n_com, n_distinct, n_var, n_cluster,  p_sparse)
```

## CSSCA with model selection (when the tuning parameters are not known)

### Functions
  
**ComputationCSSCA**: a wrapper function that estimates the CSSCA results at each possible value of the grid

**ModelSelectionCPCA**: carry out the model selection procedure to automatically select the optimal number of clusters and level of sparsity.

Note that, at this stage, the function "ModelSelectionPCA" could only be applied when each point in the grid is equally spaced. We will soon add the additional feature that allows to analyze on any kind of grid.

### Examples

```{r eval = FALSE}
### first define the grid 
cluster_range <- 1:4
sparse_range <- c(0, 0.2, 0.4, 0.6)
## compute and store the original CSSCA estimates in the local drive 
Computation(target_data,  n_com, n_distinct, n_var, n_cluster,  p_sparse)
## select the values for the tuning parameters
select.opt <- ModelSelectionCSSCA(ncluster_range, psparse_range)
## use the values of the tuning paramters to re-estimate the fixed effects 
results_opt <- FixedCSSCA(target_data, n_block, n_com, n_distinct, n_var, select.opt[[1]],  select.opt[[2]])
```

## miscellaneous

**ListCongruence**: Compute the average Tucker's Congruence of corresponding components across components and clusters


## Depends

pracma, psych, iCluster, irlba, mclust, Rcpp, doParallel, RcppArmadillo

## License

GPL-2

## Contact

Shuai Yuan (s.yuan@uvt.nl)

## Scripts to reproduce the results of the manuscript "Revealing subgroups that differ in common and distinctive variation in multi-block data: Clusterwise Sparse Simultaneous Component Analysis"
**Application Folder**: process data for the demonstrative application of CSSCA on the linked personality measurements and non-verbal behavior data  
**Simulation Folder**: Script to create and estimate the two simulation studies (simulation script) and Script to reproduce the tables, figures and other reported numbers of the simulation studies (Simulation Results Processing)

## Planned adds-on

At the current stage, only CSSCA related functions that have been described in the article are fully ready to be used. CPCA related functions still need to be tested before completely usable.  

In the next update, we plan to add the following feature: summary tables and data visualizations;
alternative model selection schemes (e.g., cross-validation). Should you have any suggestions for future developments and (or) for implementations of CSSCA and CPCA, please kindly send your suggestions to Shuai Yuan at s.yuan@uvt.nl
