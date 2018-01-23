#### validate the simple version of the variable selection procedure (i.e. purely based on co-variance matrix)
#### validation of the complete version of the variable selction procedure would be done via the function (ValidateComplete)
#### the formal version of the methods could be seen in "final simulation"

############################# the parameters #########################
################### part1: the data #####################
# the idea is to test whether the algorithm could detect the correct number of clusters and components
# therefore, the number of clusters and number of components will be treated as variables
n_block <- 2
n_var <- c(15, 15)
#n_respondents <- c(40, 40, 40)
p_noise <- 0.1
p_sparse <- 0
# in the simple version, the mean structure always constructs a minor part of the total structure
# meanv_option <- c(0.1, 0.9)
meanv <- 0.1

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
number.time <- 5
# The number of iterations
iteration <- 20

################ part4: the recording variables #######################

used_data <- list()
seeds<- list()
setting <- list()
# rows represent the selected number of clusters and columns represent the selected number of components
# create one table for each simulation
select.table <- list()
for (i in 1:number.time){
  select.table[[i]] <- matrix(nrow = number.testing.clusters, ncol = number.testing.components)
}


#########################  data generation ###########################
for (e in 1:number.testing.clusters){
  n_cluster <- 4 * e
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

      condition <- condition + 1
      
      # restore the data in the system
      used_data[[condition]] <- list()
      seeds[[condition]] <- vector("numeric", length = 25)
      setting[[condition]] <- list(num.cluster = n_cluster, meanv = meanv, num.component = n_com1)

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
        est_csca_cluster <- varsel_nonsparse(all_data, n_var, n_block, max_component, max_cluster, n_respondents, iteration)[[1]]
        est_csca_component <- varsel_nonsparse(all_data, n_var, n_block, max_component, max_cluster, n_respondents, iteration)[[2]]
        # record the results (cluster * 100 + component)
        select.table[[times]][e, c] <- est_csca_cluster * 100 + est_csca_component
        
      }
      
    }
  }



