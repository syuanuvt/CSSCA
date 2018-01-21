#### validate the simple version of the variable selection procedure (i.e. purely based on co-variance matrix)
#### validation of the complete version of the variable selction procedure would be done via the function (ValidateComplete)
#### the formal version of the methods could be seen in "final simulation"

############################# the parameters #########################
################### part1: the data #####################
# the idea is to test whether the algorithm could detect the correct number of clusters and components
# therefore, the number of clusters and number of components will be treated as variables
n_block <- 2
n_var <- c(20, 20)
n_respondents <- c(40, 40, 40)
p_noise <- 0.1
p_sparse <- 0
meanv_option <- c(0.1, 0.9)
# p_combase refers to the proportion of the identical non-zero part across clusters; 
# p_fixzero refers to the proportion of the zeros that stay in the same location
p_combase <- 0
p_fixzero <- 1
n_distinct <- 1

################# part2: the algorithm ##############
max_component <- 16
max_cluster <- 16
# the starting value which should be larger than any possible sum of loss values
upper <- 1e9

################ part3: the parameters for the iterations ############ 
condition <- 0
number.testing.clusters <- 3
number.testing.components <- 3 
# the number of simulations in one certain condition
number.time <- 1
# The number of iterations
iteration <- 1

#########################  data generation ###########################
for (e in 1:number.testing.clusters){
  n_cluster <- 5 * e
  mem_cluster <- rep(40, n_cluster)
  n_respondents <- sum(mem_cluster)
  clustermem <- vector()
  # assign the cluster membership
  for (i in 1:n_cluster){
    clustermem <- c(clustermem, rep(i, 40))
  }
  
  for (c in 1:number.testing.components){
    n_com <- 3 * c
    n_com1 <- 2 + n_com
    
    for (b in 1:2){
      meanv <- meanv_option[b]
      condition <- condition + 1
      
      # restore the data in the system
      used_data[[condition]] <- list()
      seeds[[condition]] <- vector("numeric", length = 25)
      setting[[condition]] <- list(num.cluster = n_cluster, meanv = meanv, num.component = n_com1)
      n_cluster_summary[[condition]] <- list()
      n_component_summary[[condition]] <- list()

      for (times in 1:number.time){
        # set the random seed
        seeds[[condition]][times] <- .Random.seed[times]
        set.seed(seeds[[condition]][times])
        toy <- CSSCASimulation.full.mean(n_cluster, mem_cluster, n_com, n_distinct, n_block, n_var, p_sparse, 
                                         p_noise, p_combase, p_fixzero, "both", meanv)
        block_version_data <- toy[[1]]
        all_data <- toy[[2]]
        true_score <- toy[[4]]
        true_loading <- toy[[5]]
        
        # record the data
        used_data[[condition]][[times]] <- list(data = toy, num.cluster = n_cluster, meanv = meanv, num.component = n_com1)

        
        # first select the number of components and clusters based on csca
        est_csca_cluster <- varsel_nonsparse_full(all_data, n_var, n_block, max_component, max_cluster, n_respondents, iteration)[[1]]
        est_csca_component <- varsel_nonsparse_full(all_data, n_var, n_block, max_component, max_cluster, n_respondents, iteration)[[2]]
        # then estimate the results based on the selected variables
        #results_csca <- csca(all_data, n_var, n_block, est_csca_component, est_csca_cluster, n_respondents, vip.iteration)
        n_cluster_summary <- c(n_cluster_summary, est_csca_cluster)
        n_component_summary <- c(n_component_summary, est_csca_component)
        # then select the number of componets and clusters based on other mean-level clustering methods
        #est_mean_cluster <- iCluster2(block_version_data, n_cluster)$clusters
        # based on the number of cluster, determine the number of components
        
        
      }
      
    }
  }
}



