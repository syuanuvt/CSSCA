# CSSCA: a brief tutorial of application (current version 0.6.0; accompanying research article currently under review at Social Science Computer Review)
The package includes most essential functions to carry out the Cluster-wise Sparse Simultaneous Component Analysis and the Cluster-wise Principal Component Analysis (CPCA). 

## Running simulation with CSSCA model
**CSSCASimulation**: Data simulation following CSSCA models, as described in Yuan et al., (under review), with the options to specify a full set of parameters, including: (1) the number of clusters, (2) the size of each cluster (allowed to have unbalanced cluster size), (3) the number of common and distinctive components, (4) the number of data blocks, (5) the number of variables in each data block, (6) the sparsity level of the loading matrices, (7) the desirable noise level of the data and (8) the congruence of loading matrices. 

## CPCA: Clustering analysis on single dataset based on both mean structure and component structure  
**csca_cpp**:  For the users who want to conduct clustering analysis on single data block that accounts for both mean structure and component structure, we have implemented an efficient Cluster-wise Principal Component Analysis (CPCA algorithm). In essence, detect hidden subgroups that possess group-specific centroids and group-specific component loadings) on the single dataset, the ClusterSSCA package also implemented such analysis, which could also be viewed as a special case of CSSCA (i.e., without separating distinctive components and sparseness induced) Clusterwise PCA clusters observation in such a way that the observations belong to the same cluster would have the same loading matrices while the observations belong to different clusters would have different loading matrices. 

## CSSCA with fixed parameters
**FixedCSSCA**: a wrapper that estimates the CSSCA results based on fixed parameters (i.e. no model selection provided)

## CSSCA with model selection
**VariousCSSCA**: a wrapper that estimates the CSSCA results with varying levels of the parameters. The function will produce the resulting data for each combination of candidating values of the unknown parameters  
**ModelSelectionCSSCA**: carry out the model selection procedure to automatically select the optimal number of clusters and level of sparsity.

## Other useful functions
**Add**: In general, the cluster differences could be decomposed into mean-level differences and co-variance differences. The function helps the users to create, based on existing covariance data, certain proportion of mean-level cluster differences and noise  
**ListCongruence**: Compute the average Tucker's Congruence of corresponding components across components and clusters

## Scripts to reproduce the results of the manuscript "Revealing subgroups that differ in common and distinctive variation in multi-block data: Clusterwise Sparse Simultaneous Component Analysis"
**Application Folder**: process data for the demonstrative application of CSSCA on the linked personality measurements and non-verbal behavior data  
**Simulation Folder**: Script to create and estimate the two simulation studies (simulation script) and Script to reproduce the tables, figures and other reported numbers of the simulation studies (Simulation Results Processing)
