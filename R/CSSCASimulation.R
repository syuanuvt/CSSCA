#' Simulate the data according to the CSSCA model
#'
#' @param ncluster the number of clusters that should be simulated
#' @param memcluster A vector indicates the amount of entries in each cluster. The vector should be of length ncluste, with
#' the \code{nth} element indicates the amount of entries in the \code{nth} cluster. It could also be an integer; in such cases, we assume all clusters have the same amount of entries.
#' @param nblcok A positive integer indicates the number of blocks (i.e. the number of data sources)
#' @param ncom An integer indicates the number of common components
#' @param ndistinct A vector of length nblock, with the \code{ith} element indicates the number of distinctive components assumed for the \code{ith} data block. It could also be an integer; in such cases, we assume all blocks have the same amount of distinctive components.
#' @param nvar A vector of length nblock, with the \code{ith} element indicates the number of variables assumed for the \code{ith} data block. It could also be an integer; in such cases, we assume all blocks have the same amount of variables.
#' @param psparse A number within the range of [0,1] that indicates the sparsity level (i.e. the proportion of zero elements in the loading matrix)
#' @param p_noise A number within the range of [0,1] that indicates the percentage of noise structrue that should be added to the final data.
#' @param pcombase A number within the range of [0,1] that indicates the percentage of the "common"(i.e. identical) part
#' in the loading matrices of various clusters. The cluster-specific part would then be (1 - pcombase). It is one of the parameter that controls
#' for the similarities between loading matrices
#' @param pfixzero A number within the range of [0,1] that indicates the percentage of the zero loadings that share the same
#' positions over all clusters. It is one of the parameter that controls
#' for the similarities between loading matrices.
#' @param meancov Possible values: "mean' = only includes mean structure, "cov" = only includes covariance structure and "both" = includes both mean structure and co-variance structure
#' @param meanp A number within the range of [0,1] that indicates the proportion of mean structure
#' @return a list of six elements. The first element is a list that includes the generated final data per block;
#' the second element is the concatenated version of the final data (concatenate the block-version data into one single dataset);
#' the third element is the data that involves cluster difference only in co-variance structure (i.e. before adding mean structure and noise stucture)
#' the forth element is a list of cluster-specific score matrices
#' the fifth element is a list of cluster-specific loading matrices
#' the last element is a vector indicates the cluster assignment (the \code{nth} element of the vector indicates
#' the cluster assignment of the \code{nth} observation)
#' @examples
#'    n_cluster <- 3
#'    mem_cluster <- c(50,50,50) # 50 entries in each cluster
#'    n_block <- 2
#'    n_com <- 2
#'    n_distinct <- c(1,1) #1 distinctive components in each block
#'    n_var <- c(15,9)
#'    p_sparse <- 0.5
#'    p_noise <- 0.3
#'    p_combase <- 0.5 # moderate similarity
#'    p_fixzero <- 0.5 # moderate similarity
#'    mean_v  <- 0.1 # co-variance structrue dominates
#'  (not run)  CSSCASimulation(n_cluster, mem_cluster, n_block, n_com, n_distinct, n_var, p_sparse,
#'  p_noise, p_combase, p_fixzero, "both", mean_v)
#'
CSSCASimulation <- function(ncluster, memcluster, nblock, ncom, ndistinct, nvar,
                                      psparse = 0, pnoise = 0, pcombase = 0, pfixzero = 0, meancov, pmean) {

  # ncluster: the number of cluster
  # memcluster: number of participants in each cluster
  # ncom: number of common components
  # ndistinct: number of distinct components (a vector of length nblock, with the ith element indicating the number of distinctive components
  # in ith data block)
  # nblcok: number of blocks
  # nvar: number of variables (vector indicates the number of variables for in each data block)
  # psparse: percentage of sparsity (vector indicates the percentage for each block)
  # pnoise: percentage of noise
  # pcombase: percentage of the common part of the loading matrices for all clusters
  # pfixzero: percentage of the common "fixed zero" of the loading matrices for all clusters (note that it could not be too close to 1, otherwise no cluster difference could be presented)
  # meancov: whetherthe mean structure and the covariance struture are involved, "mean" = mean structure only, "cov" = covariance only, "both" = both
  # meanp: the proportion of mean struture (only when meancov = 'both')

  # specification of only mean struture
  if (meancov == "mean") {
    pcombase <- 1
    pfixzero <- 1
  }

  # irregular input
  if (length(memcluster) == 1){
    memcluster <- rep(memcluster, ncluster)
  }

  if (length(nvar) == 1){
    nvar <- rep(nvar, nblock)
  }

  if (length(ndistinct) ==1){
    ndistinct <- rep(ndistinct, nblock)
  }

  # the vector indicates the relationship between distinctive component and the number of data block
  distinct_index <- vector("numeric")
  sum_var <- 0
  ini_var <- 0
  for (i in 1:nblock){
    distinct_index <- c(distinct_index, rep(i, ndistinct[i]))
    ini_var <- ini_var + nvar[i]
    sum_var <- c(sum_var, ini_var)
  }


  # the aggregate level of information
  all_component <- sum(ncom, ndistinct)
  all_member <- sum(memcluster)
  all_var <- sum(nvar)

  # for the easy computation (fill in the 0th item)
  cluster_rep <- c(0, memcluster)
  block_rep <- c(0, nvar)

  # cluster assignment (becuase of the randomness, the assignment could be arbitrary, just from the first avilable
  # observation looping towards the last avilable observation)
  cluster_mem <- vector("numeric", length = all_member)
  for  (i in 1:ncluster){
    cluster_mem[(sum(cluster_rep[1:i]) + 1):sum(cluster_rep[1:i+1])] <- i
  }

  ## create the initial score matrix
  inidata <- matrix (rnorm(all_member * all_var), nrow = all_member, ncol = all_var)
  #initial_score <- matrix(nrow = all_member, ncol = all_component)
  true_score <- matrix(nrow = all_member, ncol = all_component)
  initial_p <- list()

  # the common base of the loading matrix
  common_p <- matrix(runif(all_var * all_component, min = -1, max = 1), nrow = all_var, ncol = all_component)
  scale_common_p <- sqrt(pcombase) * sweep(common_p, 2, sqrt(colSums(common_p ^ 2)), `/`)

  # generate the cluster-specific score matrix and loading matrix
  for  (i in 1:ncluster){
    #initial_score[(sum(cluster_rep[1:i]) + 1):sum(cluster_rep[1:i+1]),]
    score_cluster <- matrix(rnorm(memcluster[i] * all_component), nrow = memcluster[i], ncol = all_component)
    score.de <- gramSchmidt(score_cluster)
    true_score[(sum(cluster_rep[1:i]) + 1):sum(cluster_rep[1:i+1]),] <- MatrixCenter_cpp(score.de$Q, 1, 0) * sqrt(nrow(score_cluster))
    # the cluster-specific loading matrix
    specific_p <- matrix(runif(all_var * all_component, min = -1, max = 1), nrow = all_var, ncol = all_component)
    initial_p[[i]] <- scale_common_p + sqrt(1 - pcombase) * sweep(specific_p, 2, sqrt(colSums(specific_p ^ 2)), `/`)
  }

  # generate the common zero positions concerned about the distinctive component (position-based generation)
  distinct_zeros <- vector("numeric")
  num_var <- 0
  for (i in 1:sum(ndistinct)){
    distinct_zeros <- c(distinct_zeros, ((all_var * (ncom + i - 1)) + sum_var[distinct_index[i]] + 1): ((all_var * (ncom + i - 1)) + sum_var[(distinct_index[i] + 1)]))
  }
  # the positions that are not yet specified as zeros
  retain_zeros <- setdiff(1: (all_var * all_component), distinct_zeros)

  # generate the sparse component loadings
  # the number of zeros in component loadings
  num_zeros <- round(length(retain_zeros) * psparse)
  # generate the positions of common(fixed) zeros for all clusters
  ## first get the number of common zeros
  num_common_zeros <- round(num_zeros * pfixzero)
  while(TRUE){
    common_zeros <- sample(retain_zeros, num_common_zeros)
    test <- initial_p[[1]]
    test[common_zeros] <- 0
    test[distinct_zeros] <- 0
    if (sum(apply(test, 2, sum) == 0) == 0)
      break
  }
  ## the zero positions that are remained after the common zeros been imposed
  retain_zeros_cluster <- setdiff(retain_zeros, common_zeros)

  # list of component loadings for all clusters
  loading_per_cluster <- list()
  # list of component score for all clusters
  score_per_Cluster <- list()
  # list of the positions of 'zero' elements (e.g. if 1 is included, the first element of the loading matrix will be forced to 0)
  zero_number_per_cluster <- list()

  # not allow all the loadings of a specific component are zero
  forbiden_list <- lapply(seq(from = 0, by = all_var, length.out = ncom), function(x) x + 1:all_var)

  if (pcombase == 1 & pfixzero == 1){
    for (i in 1:ncluster){
      copy_p <- initial_p[[i]]
      copy_p[common_zeros] <- 0
      copy_p[distinct_zeros] <- 0
      # the loading matrix of ith cluster
      loading_per_cluster[[i]] <- copy_p
    }
  }

  #when not mean structure only
  if (pcombase != 1 | pfixzero != 1) {
    # generate the loading matrix
    for (i in 1:ncluster){
      while (TRUE){
        # random generate the positions of zero loadings
        loading_sample <- sample(retain_zeros_cluster, (num_zeros - num_common_zeros))
        # to make sure: (1) no components will have all zero loadings
        copy_p <- initial_p[[i]]
        copy_p[common_zeros] <- 0
        copy_p[distinct_zeros] <- 0
        copy_p[loading_sample] <- 0

        # no components should have all zero loadings
        if (sum(apply(copy_p, 2, sum) == 0) == 0)
          break
      }


      # the loading matrix of ith cluster
      loading_per_cluster[[i]] <- copy_p
      if (pmean != 1){

        # make the row-wise sum-of-square of the loading matrix equals to (1-pmean)
        loading_per_cluster[[i]] <- loading_per_cluster[[i]] * sqrt((1-pmean)*(1-pnoise) / mean(apply(loading_per_cluster[[i]] ^ 2, 1, sum)))
      }
    }
  }

  final_data <- matrix(nrow = all_member, ncol = all_var)
  # generate the data matrix (also adding some noise, and mean struture)
  for (i in 1:ncluster){
    final_cluster <- true_score[(sum(cluster_rep[1:i]) + 1):sum(cluster_rep[1:(i+1)]), ]%*% t(loading_per_cluster[[i]])
    final_data[(sum(cluster_rep[1:i]) + 1):sum(cluster_rep[1:(i+1)]), ] <- final_cluster
    score_per_Cluster[[i]] <- true_score[(sum(cluster_rep[1:i]) + 1):sum(cluster_rep[1:(i+1)]), ]
  }
  # add noise and mean structure
  final_data_pre <- Add(final_data, cluster_mem, ncluster, pmean, pnoise)[[1]]
  # make the total variance across variables equals to 1
  final_data <- t(sqrt(1 / apply(final_data_pre, 2, var)) * t(final_data_pre))

  # split the data matrix for each block
  final_data_block <- lapply(1: nblock, function(x) final_data[, (sum(block_rep[1:x]) + 1):sum(block_rep[1:x+1])])

  l <- list(block_data = final_data_block, concatnated_data = final_data, score = score_per_Cluster, loading = loading_per_cluster, cluster_assign = cluster_mem)
}
