# The simulation for cluster-wise sparse simultaneous component analysis
# Yuan Shuai 2017/5/2

CSSCASimulation.full.mean <- function(ncluster, memcluster, ncom, ndistinct, nblock, nvar,
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
  
   # example
   #ncluster <- 3
   #memcluster <- 50
   #common <- 2
   #distinct <- 1
   #nblock <- 2
   #nvar <- 10
   #pnoise <- 0.1
   #pcombase <- 0.5
   #pfixzero <- 0.5
   #meancov <- "both"
   #pmean <- 0.5
   #psparse <- 0.2
   
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
  true_score <- matrix(nrow = all_member, ncol = all_component)
  initial_p <- list()
  
  # the common base of the loading matrix
  common_p <- matrix(runif(all_var * all_component, min = -1, max = 1), nrow = all_var, ncol = all_component)
  scale_common_p <- sqrt(pcombase) * sweep(common_p, 2, sqrt(colSums(common_p ^ 2)), `/`)
  
  # generate the cluster-specific score matrix and loading matrix
  for  (i in 1:ncluster){
    true_score[(sum(cluster_rep[1:i]) + 1):sum(cluster_rep[1:i+1]),] <- matrix(rnorm(memcluster[i] * all_component), nrow = memcluster[i], ncol = all_component)
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
  num_common_zeros <- round(num_zeros * pfixzero)
  while(TRUE){
    common_zeros <- sample(retain_zeros, num_common_zeros)
    test <- initial_p[[i]]
    test[common_zeros] <- 0
    test[distinct_zeros] <- 0
    if (sum(apply(test, 1, sum) == 0) == 0)
    break
  }
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
        while (TRUE){
          # random generate the positions of zero loadings
          loading_sample <- sample(retain_zeros_cluster, (num_zeros - num_common_zeros))
          # to make sure: (1) the list of zero positions for each cluster is different (i.e. no two clusters have the complete same positions of zero)
          # and (2) no component will have all zero loadings
          if (sum(forbiden_list %in% list(loading_sample)) == 0 &&
              sum(zero_number_per_cluster %in% list(loading_sample)) == 0)
          break  
        }
        copy_p <- initial_p[[i]]
        copy_p[common_zeros] <- 0
        copy_p[distinct_zeros] <- 0
        copy_p[loading_sample] <- 0
        
        if (sum(apply(copy_p, 1, sum) == 0) == 0)
        break
      }
      
      
      # the loading matrix of ith cluster
      loading_per_cluster[[i]] <- copy_p 
      if (pmean != 1){
        # make the row-wise sum-of-square of the loading matrix equals to (1-pmean)
        loading_per_cluster[[i]] <- sqrt((1-pmean)*(1-pnoise) / apply(loading_per_cluster[[i]] ^ 2, 1, sum)) * loading_per_cluster[[i]]
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
    final_data <- Add(final_data, cluster_mem, ncluster, pmean, pnoise)
    # make the total variance across variables equals to 1
    final_data <- t(sqrt(1 / apply(final_data, 2, var)) * t(final_data))
  
  # split the data matrix for each block
  final_data_block <- lapply(1: nblock, function(x) final_data[, (sum(block_rep[1:x]) + 1):sum(block_rep[1:x+1])])
  
   l <- list(final_data_block, final_data, inidata, score_per_Cluster, loading_per_cluster, cluster_mem)
}  

