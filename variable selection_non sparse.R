# performing variable selection in the non-sparse settings

varsel_nonsparse <- function(xc, nvar, nblock, max_components, max_cluster, nrespondents, iteration){
  # xc: concatenated data
  # nvar: number of variables in each data block
  # nblock: number of blocks
  # all_components: number of all components
  # ncluster: number of clusters
  
  ######################specification of the certain values########################
  #xc <- all_data
  #nvar <- n_var
  #nblock <- n_block
  #max_component <- 12
  #max_cluster <- 12
  #nrespondents <- n_res
  # first try out the quickest version
  #iteration <- 1
  
  #################################################################################
  
  vaf.table <- matrix(nrow = max_cluster, ncol = max_components)
  # the initiate values
  current.cluster <- max_cluster
  current.component <- max_component
  old.cluster <- 0
  old.component <- 0
  norm.data <- sum(all_data ^ 2)
  
  while (1){
    # first update the components, conditional on the number of clusters
    # the first three VAF will be calculated for sure (so the seletion will not be valid
    # if the number of clusters < 4 or the number of components < 4) 
    if (current.cluster == old.cluster) break
      for (i in 1:3){
        if (is.na(vaf.table[current.cluster, i])){
          min_loss <- min(csca(xc, nvar, nblock, i, current.cluster, nrespondents, iteration)[[2]])
          vaf.table[current.cluster, i] <- (norm.data - min_loss) / norm.data
        }
      }
      # initiate the scree ratio
      sr.current <- (vaf.table[current.cluster, 2] - vaf.table[current.cluster, 1]) / (vaf.table[current.cluster, 3] - vaf.table[current.cluster, 2]) 
      # find the best fitted number of components
      while (1){
        sr.old <- sr.current
        i <- i + 1
        if (is.na(vaf.table[current.cluster, i])){
          min_loss <- min(csca(xc, nvar, nblock, i, current.cluster, nrespondents, iteration)[[2]])
          vaf.table[current.cluster, i] <- (norm.data - min_loss) / norm.data
        }
        sr.current <- (vaf.table[current.cluster, (i-1)] - vaf.table[current.cluster, (i-2)]) / (vaf.table[current.cluster, i] - vaf.table[current.cluster, (i-1)])
        # stop the iteration 
        if (sr.current > sr.old)  {
          old.component <- current.component
          current.component <- (i - 1)
          break
        }
        if (i == max_component) {
          old.component <- current.component
          current.component <- (i - 1)
          break
        }
      }
      
    if (current.component == old.component) break
      for (i in 1:3){
        if (is.na(vaf.table[i, current.component])){
          min_loss <- min(csca(xc, nvar, nblock, current.component, i, nrespondents, iteration)[[2]])
          vaf.table[i, current.component] <- (norm.data - min_loss) / norm.data
        }
      }
      # initiate the scree ratio
      sr.current <- (vaf.table[2, current.component] - vaf.table[1, current.component]) / (vaf.table[3, current.component] - vaf.table[2, current.component]) 
      # find the best fitted number of components
      while (1){
        sr.old <- sr.current
        i <- i + 1
        if (is.na(vaf.table[i, current.component])){
          min_loss <- min(csca(xc, nvar, nblock, current.component, i, nrespondents, iteration)[[2]])
          vaf.table[i, current.component] <- (norm.data - min_loss) / norm.data
        }
        sr.current <- (vaf.table[(i-1), current.component] - vaf.table[(i-2), current.component]) / (vaf.table[i, current.component] - vaf.table[(i-1), current.component])
        # stop the iteration 
        if (sr.current > sr.old)  {
          old.cluster <- current.cluster
          current.cluster <- (i - 1)
          break
        }
        if (i == max_cluster) {
          old.cluster <- current.cluster
          current.cluster <- (i - 1)
          break
        }
      }
  }
  
  return(list(current.cluster, current.component, vaf.table))
}