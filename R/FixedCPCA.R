#' FixedCPCA
#'
#' CPCA is a clustering method that could be viewed as a special case of CSSCA. In essence, CPCA works on 
#' single-block data, and extracts various unknown sub-groups that differ in both mean structure and component structure.
#' In other words, if the assumption - the underlying data generation mechanism suggests the existence of various clusters 
#' that differ in both mean structure and component structure - holds, CPCA is able to find back these clusters. AS discussed on 
#' the potential application of the CSSCA algorithm, we would like to suggest the researcher to use other distance-based clustering methods (e.g., k-means) or 
#' model-based clustering methods (e.g., GMM) if they deem an alternative data generation mechanism is more appropriate.
#' The current function implements CPCA with known values of paramters (i.e. the number of clusters). When such theoritical knowledge
#' is availble, the model selection is not necessary and the current function could be used to directly compute cluster-wise mean structure and
#' component structure
#'
#' @param all_data A matrix with concatenated data (the aggregation of the data blocks by rows (entries)). The CSSCA method will be performed on the data.
#' @param n_total An integer indicates the number of components assumed
#' @param n_cluster A positive integer indicates the number of clusters assumed. When ncluster = 1, CSSCA boils down to PCA analysis.
#' @param computation Either "easy", "medium" or "difficult" that indicates the desired complexity of the computation. "easy" mode in most cases
#' is already quite enough to produce accurate model estimates. However, when complex datasets are to be estimated, it is favorable to use the "difficult" mode so as to increase the
#' accuracy yet sacrifice the speed of the algorithm. The default computation mode is "easy"
#' @return a list of four elements. The first element is vector that indicates the resulting cluster partitions, the nth element refers
#' to the cluster assignment of the nth entry;
#' the second element is a numeric value that is the minimal loss function value;
#' the third element is a list that displays cluster-specific loading matrices;
#' the forth element is a list that displays cluster-specific score matrices;
#' @examples
#' the setting for simulation and calculation
#' n_cluster <- 3
#' mem_cluster <- c(30,30,30)# 50 entries in each cluster
#' n_obs <- sum(mem_cluster)
#' n_block <- 2
#' n_com <- 2
#' n_distinct <- c(1,1) #1 distinctive components in each block
#' n_var <- c(15,9)
#' p_sparse <- 0
#' p_noise <- 0.3
#' p_combase <- 0.5 # moderate similarity
#' p_fixzero <- 0.5 # moderate similarity
#' mean_v  <- 0.1 # co-variance structrue dominates
#' # extract the data from the simulation
#' (not run)  sim <- CSSCASimulation(n_cluster, mem_cluster, n_block, n_com, n_distinct, n_var, p_sparse,
#'  p_noise, p_combase, p_fixzero, "both", mean_v)
#'  target_data <- sim$concatnated_data
#'  # feed the data with original cluster assignment and estimate with the CSSCA method
#'  results <- FixedCPCA(target_data,  n_com, n_cluster)
#'
FixedCPCA <- function(all_data, n_total, n_cluster, computation = "easy"){
  
  ##############################################################################################################
  ########################### obtain additional information from the input #####################################
  ##############################################################################################################
  
  # additional information from the input
  n_observation <- nrow(all_data)
  n_var <- ncol(all_data)
  
  all_data <- base::scale(all_data)
  
  # the number of blocks is not important here
  n_block <- 2
  
  if(computation == "test"){
    csca_times <- 1
  }
  
  if(computation == "easy"){
    csca_times <- 30
  }
  if(computation == "medium"){
    csca_times <- 50
  }
  if(computation == "difficult"){
    csca_times <- 100
  }
  upper <- 1e9
  
  if(n_cluster != 1){
    global <- csca_cpp(all_data, n_var, n_block, n_total, n_cluster, csca_times)
  }
  
  if(n_cluster == 1){
    
    results <- sca_common_cpp(all_data, n_com)
    
    global <- list(cluster_mem = rep(1, n_observation), loadings = list(loading_result), scores = list(score_result)) 
  }
  
  return(global)
}