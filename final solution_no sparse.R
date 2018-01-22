# the initial version: known sparsity
## created by Yuan Shuai 2017/5/6
# update p after itration + update p after comparison

csca <- function(xc, nvar, nblock, all_components, ncluster, nrespondents, iteration){
  # xc: concatenated data
  # nvar: number of variables in each data block
  # nblock: number of blocks
  # all_components: number of all components
  # ncluster: number of clusters
  # memcluster: number of members in each cluster
  
  # set the upper bound of loss function
  upper <- 1e9
  # set the minimum converge error
  converge_stage1 <- 1e-2
  converge_stage2 <- 1e-5
    
  # set the maximum iteration time for one single random start
  MAXITER <- 100
  
  if (length(nvar) == 1){
    nvar <- rep(nvar, nblock)
  }
  
  # the aggregate level of information
  all_member <- nrespondents
  all_var <- sum(nvar)
  
  # set the upper bound of minimum loss 
  loss_min <- upper
  loss_all <- vector("numeric", length = iteration)
  
  # the starting partition
  start <- vector("numeric", length = ncluster)
  for (y in 1:ncluster){
    start[y] <- round(all_member / ncluster)
  }
  start[ncluster] <- all_member - (ncluster - 1) * round(all_member / ncluster)
  
for (v in 1:iteration){
  
  flag <- 0    # flag for breaking the loop
     
  loss_iteration <- vector("numeric")

  # step 1: initialize the cluster assignment
  random_start <- randomstart(start)
  cluster_yn <- random_start$mem
  cluster_mem <- random_start$ind
  
  # step 2: perform the SSCA for each cluster and calculate the initial loss function
  current_data <- list()
  current_loading <- list()
  current_score <- list()
  current_loss <- list()
  temp_score <- list()
  temp_loss <- list()
  temp_loading <- list()
  
  # the members in a certain cluster
  current_mem <- list()
  loss_tot <- 0
  
  # not random start (rational start from the based on the value of singular value decomposition)
   for (i in 1:ncluster){
    current_data[[i]] <- xc[which(cluster_mem == i), ]
    ssca_results <- Sca_common(current_data[[i]], nvar, nblock, all_components)
    loss_tot <- loss_tot + ssca_results$l
    current_loading[[i]] <- ssca_results$p
    current_score[[i]] <- ssca_results$t
    current_loss[[i]] <- ssca_results$l
    current_mem[[i]] <- which(cluster_mem == i)
  }
  
  # check the non-increase of the loss function
  cat("Update Loss: ",loss_tot, sep="\n")

  # iteration 
  # step 3: update cluster assignment for each observation
  conv <- 0
  ## note that loss_tot only update after one whole iteration, and loss_tot_ind update after every observation
  loss_tot_ind <- loss_tot
  iter <- 0
  
  while (conv == 0){
    iter <- iter + 1
    # the sign of whethe there ar emembership exchanges
    sign <- 0
    # check the fit of a particular observation with a particular cluster 
    for (j in 1:all_member){
      # set a impossible upper bound
      min_change <- upper
      for (k in 1:ncluster){
        # the member belongs to a certain cluster
        if (cluster_mem[j] == k){
          ## get the loss function without observation j
          # the data matrix without j
          new_data <- xc[setdiff(current_mem[[k]], j), ]
          
          # stop the loop when the new_data contains fewer than two rows
          if (nrow(new_data) < 3){
            flag <- 1  # flag for breaking the loop
            break      # for the loop k
          }
          
          # update both the T and P
          # compute the new score matrix (use rational start, only compute one time)
          new_ssca <- Sca_common(new_data, nvar, nblock, all_components)
          # update the current loading
          temp_loading[[k]] <- new_ssca$p
          # update the current score
          temp_score[[k]] <- new_ssca$t
          # update the current loss
          temp_loss[[k]] <- new_ssca$l
          
          # the loss function with observation j is the old loss function
          # we could thus compute the change in loss
          loss_change <- current_loss[[k]] - temp_loss[[k]]
        }
        
        # the member does not belong to a certain cluster
        if (cluster_mem[j] != k) {
          # the loss function with observation j
          # the data matrix with j
          new_data <- xc[c(current_mem[[k]], j), ]
          # compute the new score matrix (use rational start, only compute one time)
          new_ssca <- Sca_common(new_data, nvar, nblock, all_components)
          # update the current loading
          temp_loading[[k]] <- new_ssca$p
          # update the current score
          temp_score[[k]] <- new_ssca$t
          # update the current loss
          temp_loss[[k]] <- new_ssca$l
          
          # the loss function without observation j is the old loss function
          # we could thus compute the change in loss
          loss_change <- temp_loss[[k]] - current_loss[[k]]
        }
        
        # assign jth eobservation to the best fit cluster
        if (loss_change < min_change){
          min_change <- loss_change
          # update the membership status of this certain observation
          mem_update <- k
        }
      }
      
      # break the loop j
      if (flag == 1) break
        
      # update the cluster membership, loading matrix and the loss function (loss function only for indication) immediately 
      # (only when the member is in differenet cluster)
      if  (cluster_mem[j] == mem_update){
        # check the non-increase of the loss function
        cat("Update Loss: ",loss_tot_ind, sep="\n")
      }
      
      if (cluster_mem[j] != mem_update){
        sign <- 1
        # update the cluster membership
        c <- cluster_mem[j]
        current_mem[[c]] <- setdiff(current_mem[[c]], j)
        current_mem[[mem_update]] <- c(current_mem[[mem_update]], j)
        
        # the loss function of other clusters
        loss_others <- loss_tot_ind - current_loss[[c]] - current_loss[[mem_update]]
        
        # update the cluster which previously includes the current observation (exclude the current observation)
        # update the current loading
        current_loading[[c]] <- temp_loading[[c]]
        # update the current score
        current_score[[c]] <- temp_score[[c]]
        # update the current loss
        current_loss[[c]] <- temp_loss[[c]]
        
        # update the cluster which has the best fit with the current observation (include the current observation)
        # update the current loading
        current_loading[[mem_update]] <- temp_loading[[mem_update]]
        # updating the current score
        current_score[[mem_update]] <- temp_score[[mem_update]]
        # updating the current loss
        current_loss[[mem_update]] <- temp_loss[[mem_update]]
        
        # update the overall membership vector
        cluster_mem[j] <- mem_update
        
        # update the loss function
        loss_tot_ind <- loss_others + current_loss[[c]] + current_loss[[mem_update]]
        
        cat("Update Loss: ",loss_tot_ind, sep="\n")
      }
    }
    
    # break the while loop
    if (flag == 1) break
    
    # update the cluster membership
    # recompute the score matrix and loading matrix for every cluster
    old_loss_tot <- loss_tot
    loss_tot <- loss_tot_ind
    
    loss_tot_change <- old_loss_tot - loss_tot
    # check the non-increase of the loss function
    cat("Update Loss In this Iteration: ",loss_tot, sep="\n")
    
    # stop citerian
    if (sign == 0) {
      conv <- 1
    }
    if (iter == MAXITER){
      conv <- 1
    }
    
    loss_iteration <- c(loss_iteration, loss_tot)
  }
  # record all the loss function
  loss_all[v] <- loss_tot
  
  if (loss_tot < loss_min){
    loss_min <- loss_tot
    opt_loading <- current_loading
    opt_score <- current_score
    opt_mem <- current_mem
    opt_cluster_mem <- cluster_mem
  }
  
}

  # return
  final_score <- matrix(nrow = all_member, ncol = all_components)
  final_loading <- list()
  
  for (t in 1:ncluster){
    for (s in 1:length(opt_mem[[t]])){
      final_score[opt_mem[[t]][s], ] <- opt_score[[t]][s, ]
    }
    final_loading <- opt_loading[[t]]
  }
  
  results <- list(opt_cluster_mem, loss_all, loss_iteration)
  return (results)
}