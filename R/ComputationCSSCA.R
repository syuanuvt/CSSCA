#' Calculate the CSSCA results in various conditions. Note that because of the large scale of the output, it is recommend to retain sufficient internal storage space for the current function.
#' To guaratee that the model selection algorithm could work properly, it is require to have at least 4 elements in the selection range
#' Note that the function would automatically do parallel computations use the multiple core of the computer
#'
#' @param all_data A matrix with concatenated data (the aggregation of the data blocks by rows (entries)). The CSSCA method will be performed on the data.
#' @param n_blcok A positive integer indicates the number of blocks (i.e. the number of data sources)
#' @param n_com A positive integer indicates the number of common components
#' @param n_distinct A vector of length nblock, with the ith element indicates the number of distinctive components assumed for the ith data block. It could also be an integer; in such cases, we assume all blocks have the same amount of distinctive components.
#' @param n_var A vector of length nblock, with the ith element indicates the number of variables assumed for the ith data block. It could also be an integer; in such cases, we assume all blocks have the same amount of variables.
#' @param ncluster_range a vector indicates the range of number of clusters that could be selected from.
#' All elements in the vector should be positive integers. Repeatation of the elements is not allowed in the vector.
#' @param psparse_range a vector indicates the range of the sparsiy level that could be selected from.
#' All elements in the vector should be within the range of [0,1]. Repeatation of the elements is not allowed in the vector.
#' @param computation The level of complexity the computation should bear (choose from "easy", "medium", "difficult").
#' Higher computational complexity would lead to longer computational time but could also always provide more accurate results.
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
#' p_sparse <- 0.5
#' p_noise <- 0.3
#' p_combase <- 0.5 # moderate similarity
#' p_fixzero <- 0.5 # moderate similarity
#' mean_v  <- 0.1 # co-variance structrue dominates
#' # the custimerized range for paramter selection
#' cluster_range <- 1:4
#' sparse_range <- c(0, 0.1, 0.3, 0.5)
#' simulate the data with the function CSSCASimulation
#' (not run)  CSSCASimulation(n_cluster, mem_cluster, n_block, n_com, n_distinct, n_var, p_sparse,
#'  p_noise, p_combase, p_fixzero, "both", mean_v)
#'  calculate the results of CSSCA in various conditions and save the results in the current working directory.
#'  Note that the function may take up very long time to finish
#' (not run) ComputationCSSCA(sim$concatnated_data, n_block, n_com, n_distinct, n_var, cluster_range, sparse_range, computation = "easy")

ComputationCSSCA <- function(all_data, n_block, n_com, n_distinct, n_var, ncluster_range, psparse_range, cutoff.prop = 1/6, n_replicate = 3, rate = 1/10, computation = "easy"){

   # number of candidates per parameter
   selection_table <- expand.grid(ncluster = ncluster_range, prob_sparsity = psparse_range)

   # parallel computation
   no_cores <- detectCores() - 1
   c1 <- makePSOCKcluster(no_cores)
   registerDoParallel(c1)
    foreach(i = 1:nrow(selection_table), .verbose = TRUE,
                      .packages = c( "psych", "irlba", "iCluster", "mclust",
                                                            "Rcpp",
                                                            "RcppArmadillo",
                                                            "foreach", "doParallel", "ClusterSSCA"), .combine=rbind) %dopar%{

              p_sparse <- selection_table$prob_sparsity[i]
              n_cluster <- selection_table$ncluster[i]

              results_fixed <- FixedCSSCA(all_data, n_block, n_com, n_distinct, n_var, n_cluster, p_sparse, cutoff.prop, n_replicate, rate, computation)
              save(results_fixed, file=paste0(i, ".RData"))
          }
   stopCluster(c1)
}


