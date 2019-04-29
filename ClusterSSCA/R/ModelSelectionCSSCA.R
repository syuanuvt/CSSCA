#' Conduct model selection for CSSCA. Note that the estimation results of all conditions (by the function
#' VariousCSSCA) should already exist in the current working directory
#'
#' @param ncluster_range a vector indicates the range of number of clusters that could be selected from.
#' All elements in the vector should be positive integers. Repeatation of the elements is not allowed in the vector.
#' @param psparse_range a vector indicates the range of the sparsiy level that could be selected from.
#' All elements in the vector should be within the range of [0,1]. Repeatation of the elements is not allowed in the vector.
#' @return The selected level of sparsity and number of clusters
#' @examples
#' (following the example used in demnstatig the simulation and calculation function)
#' n_cluster <- 3
#' mem_cluster <- c(50,50,50) # 50 entries in each cluster
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
#'
#' simulate the data with the function CSSCASimulation
#' (not run)  CSSCASimulation(n_cluster, mem_cluster, n_block, n_com, n_distinct, n_var, p_sparse,
#' p_noise, p_combase, p_fixzero, "both", mean_v)
#' compute the CSSCA results in various settings and save the data in local directory
#' (not run) VariousCSSCA(sim$concatnated_data, n_block, n_com, n_distinct, n_var, cluster_range, sparse_range, computation = "easy")
#' (not run) ModelSelectionCSSCA(cluster_range, sparse_range)

ModelSelectionCSSCA <- function(ncluster_range, psparse_range){

  compute_matrix <- expand.grid(ncluster_select = cluster_range, psparse_select = sparse_range)
  compute_matrix$loss <- rep(0, nrow(compute_matrix))

  # record all the loss function values from the estimation
  for (j in 1:nrow(compute_matrix)){
    load(paste0(j, ".RData"))
    compute_matrix$loss[j] <- results_fixed[[2]]
  }

  # copy the dataset
  a <- compute_matrix
  sr.result <- matrix(nrow = length(sparse_range), ncol = (length(cluster_range) - 2))
  # when the algorithm does not end normally, the result is set at NA
  a$loss[which(a$loss == upper)] <- NA
  ## first select the optimal number of clusters for each dataset
  for (i in 1:length(sparse_range)){
    filter.data <- a %>%
      dplyr::filter(psparse_select == sparse_range[i]) %>%
      dplyr::select(loss)
    filter.data <- as.matrix(filter.data)
    for (j in 1:ncol(sr.result)){
      if (!is.na(filter.data[j]) & !is.na(filter.data[j + 1]) & !is.na(filter.data[j+2]))
        sr.result[i, j] <- (filter.data[j] - filter.data[j+ 1]) / (filter.data[j+1] - filter.data[j+ 2])
    }
  }
  sum <- apply(sr.result, 2, sum, na.rm = TRUE)
  opt_cluster <- cluster_range[which(sum == max(sum)) + 1]

  ## then select the optimal number of sparsity for each dataset
  filter.data <- a %>%
    dplyr::filter(ncluster_select == opt_cluster) %>%
    dplyr::select(loss)
  filter.data <- as.matrix(filter.data)
  sr.vector <- rep(0, (length(sparse_range)-2))
  for (j in 1:(length(sparse_range)-2)){
    sr.vector[j] <- (filter.data[j+2] - filter.data[j+ 1]) / (filter.data[j+1] - filter.data[j])
  }
  opt_psparse <- sparse_range[which(sr.vector == max(sr.vector))+1]

  return (opt_cluster = opt_cluster, opt_psparse = opt_psparse)
}
