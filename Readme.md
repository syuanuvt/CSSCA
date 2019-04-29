# CSSCA
cluster-wise sparse simultaneous component analysis: inlcuding the non-sparse version (clusterwise PCA)

## Running simulation 
**simulation function_full_mean**: the full data simulation function, includes the settings of the sparse level, the number of clusters, the cluster assignment, the number of components, the proportion of mean struture to the covariance structure, the noise level, the number of blocks and the number of variables

## Clusterwise PCA: hybrid method that includes PCA and clustering. 
**csca_cpp**:  Clusterwise PCA clusters observation in such a way that the observations belong to the same cluster would have the same loading matrices while the observations belong to different clusters would have different loading matrices. The function is written based on rcpp

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
