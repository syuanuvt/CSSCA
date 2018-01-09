# CSSCA
cluster-wise sparse simultaneous component analysis: inlcuding the non-sparse version of the cluster-wise simultaneous component analysis as the special case

## Initial Set-up
**function package**: including the package that is required to loaded, and the functions that serve as the basic functions for the CSSCA analysis
**simulation function_full_mean**: the full data simulation function, includes the settings of the sparse level, the number of clusters, the cluster assignment, the number of components, the proportion of mean struture to the covariance structure, the noise level, the number of blocks and the number of variables

## CSCA: the non-sparse version of the cluster-wise sparse simultaneous component analysis
**SCA common:** (non-sparse version of) simultaneous component analysis for linked data. Fast SVD algorithm is used to compute the solutions, therefore the algorithm could not distinguish between common and specific components (only common components are obtained)
**final solution_no sparse:** Clusterwise-simultaneous component analysis for linked data. Given the number of components and the number of clusters, the function could compute the cluster assignments of the obervations and the cluster-specific loading matrices
**variable selection:** determine the best fitted number of components and number of clusters
**final simulation_non sparse version:** The try-out function for overall monitory simulation
