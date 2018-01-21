#### validate the simple version of the variable selection procedure (i.e. purely based on co-variance matrix)
#### validation of the complete version of the variable selction procedure would be done via the function (ValidateComplete)
#### the formal version of the methods could be seen in "final simulation"

# then estimate the results based on the selected variables
#results_csca <- csca(all_data, n_var, n_block, est_csca_component, est_csca_cluster, n_respondents, vip.iteration)
n_cluster_summary <- c(n_cluster_summary, est_csca_cluster)
n_component_summary <- c(n_component_summary, est_csca_component)
# then select the number of componets and clusters based on other mean-level clustering methods
#est_mean_cluster <- iCluster2(block_version_data, n_cluster)$clusters
# based on the number of cluster, determine the number of components