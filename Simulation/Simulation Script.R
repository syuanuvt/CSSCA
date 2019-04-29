## Simulation script created for the paper CSSCA
## date of last edits : 31-03-2018
## Note that the simulation 1 has been splited into two parts
## further note that the simulation script does not utilizes two of the functions in package ClusterSSCA: FixedCSSCA and VariousCSSCA

######################################################################################
############## sim1 part 1: small-scale data (nvar = c(15,15) ########################
######################################################################################
setwd("~/sim1_data1")
## set a random seed
set.seed(912)
## fixed parameters
nblock <- 2
nvar <- c(15, 15)
ncom <- 2
distinct <- c(1, 1)

## factorial designs
n_cluster <- c(2, 4)
## 0 = unequal, 1 = equal
mem_equal <- c(0, 1)
p_sparse <- c(0.3, 0.5, 0.7)
p_noise <- c(0.1, 0.2, 0.3)
mean_level <- c(0.1, 0.5, 0.9)
cong_level <- c(0, 0.7)
## we design 10 replications in each condition
replication <- 10

# additional information from the input
sum_var <- sum(nvar)
n_total <- sum(ncom, distinct)

# running information (i.e. the "difficult" computation in ClusterSSCA)
partition_small_times <- 12 #the maximum number of replacement when using CSSCA in mean-dominant structure and (or) using iCluster in covariance-dominant structure (i.e. inconsistency of the method and the structure)
partition_big_times <- 25 #the maximum number of replacement when using CSSCA in mean-dominant structure and (or) using iCluster in covariance-dominant structure (i.e. consistency of the method and the structrue)
csca_times <- 25
upper <- 1e9

# compress all paarmeters in the design matrix
design_matrix <- expand.grid(n_cluster = n_cluster, mem_equal = mem_equal,
                             p_sparse = p_sparse, p_noise = p_noise, mean_level = mean_level,
                             cong_level = cong_level)
design_matrix_replication <- design_matrix[rep(1:nrow(design_matrix), times = replication), ]

# run the script with multi-core system
no_cores <- detectCores() - 1
c1 <- makePSOCKcluster(no_cores)
registerDoParallel(c1)
# to store the final results (including congruence between the simulation loading matrics, ARI of iCluster(which could already been obtained here),
# ARI of iCluster, and the congruence between obtained loading matrics and simulatd loading matrices)
results_sim1_data1 <- foreach(i = 1:nrow(design_matrix_replication),
                                       .packages = c("rARPACK", "psych", "irlba", "mclust",
                                                     "combinat", "GPArotation",
                                                     "tidyr", "ez", "dplyr",  "Rcpp",
                                                     "RcppArmadillo", "iCluster",
                                                     "foreach", "doParallel", "ClusterSSCA"), .combine=rbind) %dopar%{


                                                       # set the specific values of the parameters
                                                       ncluster <- design_matrix_replication$n_cluster[i]
                                                       psparse <- design_matrix_replication$p_sparse[i]
                                                       pnoise <- design_matrix_replication$p_noise[i]
                                                       mean <- design_matrix_replication$mean_level[i]
                                                       cong <- design_matrix_replication$cong_level[i]
                                                       memequal <- design_matrix_replication$mem_equal[i]

                                                       # the number of observations in each cluster
                                                       if (memequal == 1){
                                                         cluster_mem <- rep(50, ncluster)
                                                       }
                                                       if (memequal == 0){
                                                         cluster_mem <- c(rep(50, (ncluster - 1)), 30)
                                                       }


                                                       # data simulation
                                                       sim_data <- CSSCASimulation(ncluster, cluster_mem, nblock, ncom,
                                                                                   distinct, nvar, psparse, pnoise, cong, cong, "both", mean)

                                                       block_data <- sim_data[[1]]
                                                       all_data <- sim_data[[2]]
                                                       loadings <- sim_data[[4]]
                                                       cluster_assign <- sim_data[[5]]

                                                       # cluster analysis with iCluster and examine the clustering accuracy
                                                       icluster_results <- iCluster2(block_data, ncluster)$clusters
                                                       # r1 represents the similarities between the estimated partitons of iCluster and the true partitions
                                                       r1 <- adjustedRandIndex(cluster_assign, icluster_results)

                                                       # the average congruence between the loading matrices
                                                       n <- 0
                                                       sum <- 0
                                                       for (j in 1:(ncluster - 1)){
                                                         for (k in (j + 1):ncluster){
                                                           n <- n + 1
                                                           sum <- sum + TuckerCongruence(loadings[[j]], loadings[[k]])
                                                         }
                                                       }
                                                       # r2 represents the average congruence between cluster-specific true loading matrices
                                                       r2 <- sum / n


                                                       # computation with CSSCA (i.e. same procedures as described in FixedCSSCA)
                                                       n_observation <- nrow(all_data)

                                                       results1 <- csca_cpp(all_data, nvar, nblock, n_total, ncluster, csca_times)

                                                       # the two rational starts
                                                       results1_cssca <- cssca_quick_cpp(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, results1$cluster_mem, 1/6)
                                                       results2_cssca <- cssca_quick_cpp(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, icluster_results, 1/6)

                                                       # Compare the results of the two rational starts to determine the number of semi-rational starts
                                                       covariance_similarity <- adjustedRandIndex(results1$cluster_mem, results1_cssca$cluster_mem)

                                                       # compare between icluster and cssca
                                                       mean_similarity <- adjustedRandIndex(icluster_results, results2_cssca$cluster_mem)

                                                       ######### CSCA part
                                                       min_loss_our <- upper
                                                       if (results1_cssca$loss < min_loss_our) {
                                                         min_loss_our <- results1_cssca$loss
                                                         global1 <- results1_cssca
                                                       }

                                                       ######### iclust part
                                                       min_loss_iclust <- upper
                                                       if (results2_cssca$loss < min_loss_iclust) {
                                                         min_loss_iclust <- results2_cssca$loss
                                                         global2 <- results2_cssca
                                                       }

                                                       #### weight scheme
                                                       ## when covariance similarity < mean similarity, mean structure probably dominates the overall structure
                                                       if (covariance_similarity < mean_similarity){
                                                         mean <- 1
                                                         n_partition_csca <- partition_small_times
                                                         n_partition_iclust <- partition_big_times
                                                       }
                                                       ## when covariance similarity > mean similarity, m=covariance structure probably dominates the overall structure
                                                       if ((covariance_similarity > mean_similarity) | (covariance_similarity = mean_similarity)){
                                                         mean <- 0
                                                         n_partition_csca <- partition_big_times
                                                         n_partition_iclust <- partition_small_times
                                                       }

                                                       ### partition for the csca part
                                                       partition_results_csca <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, results1$cluster_mem, 1/6, n_partition_csca, 3, 1/10)
                                                       if (partition_results_csca$loss < min_loss_our){
                                                         min_loss_our <- partition_results_csca$loss
                                                         global1 <- partition_results_csca
                                                       }

                                                       #### partition for thr mclust part
                                                       partition_results_iclust <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, icluster_results, 1/6, n_partition_iclust, 3, 1/10)
                                                       if (partition_results_iclust$loss < min_loss_iclust){
                                                         min_loss_iclust <- partition_results_iclust$loss
                                                         global2 <- partition_results_iclust
                                                       }


                                                       if (min_loss_our < min_loss_iclust){
                                                         global <- global1
                                                       }  else {
                                                         global <- global2
                                                       }

                                                       ## summary of the results
                                                       est_loadings <- global$loadings
                                                       est_loss <- global$loss
                                                       est_cluster_mem <- global$cluster_mem

                                                       # r3 represents the average congruence between the estimated loadings and true loadings
                                                       r3 <- ListCongruence(loadings, est_loadings, ncluster)
                                                       # r4 represents the similarities between the estimated cluster partition of CSSCA and the true cluster partition
                                                       r4 <- adjustedRandIndex(est_cluster_mem, cluster_assign)
                                                       values <- c(r1, r2, mean(r3), r4, est_loss)

                                                       # save data
                                                       out <- list(all_data = all_data, sim_loadings = loadings, sim_assign = cluster_assign, est = global, values = values)

                                                       save(out, file=paste0(i, ".RData"))

                                                     }
stopCluster(c1)

######################################################################################
############## sim1 part 2: large-scale data nvar = c(15,150) ########################
######################################################################################
### note that except for the variable nvar, the values of other parameters stay unchanged comparing to the first part of simulation 1
setwd("~sim1_data2")
## set a random seed
set.seed(912)
## fixed parameters
nblock <- 2
nvar <- c(15, 150)
ncom <- 2
distinct <- c(1, 1)

## factorial designs
n_cluster <- c(2, 4)
## 0 = unequal, 1 = equal
mem_equal <- c(0, 1)
p_sparse <- c(0.3, 0.5, 0.7)
p_noise <- c(0.1, 0.2, 0.3)
mean_level <- c(0.1, 0.5, 0.9)
cong_level <- c(0, 0.7)
## we design 10 replications in each condition
replication <- 10

# additional information from the input
sum_var <- sum(nvar)
n_total <- sum(ncom, distinct)

# running information (i.e. the "difficult" computation in ClusterSSCA)
partition_small_times <- 12 #the maximum number of replacement when using CSSCA in mean-dominant structure and (or) using iCluster in covariance-dominant structure (i.e. inconsistency of the method and the structure)
partition_big_times <- 25 #the maximum number of replacement when using CSSCA in mean-dominant structure and (or) using iCluster in covariance-dominant structure (i.e. consistency of the method and the structrue)
csca_times <- 25
upper <- 1e9

# compress all paarmeters in the design matrix
design_matrix <- expand.grid(n_cluster = n_cluster, mem_equal = mem_equal,
                             p_sparse = p_sparse, p_noise = p_noise, mean_level = mean_level,
                             cong_level = cong_level)
design_matrix_replication <- design_matrix[rep(1:nrow(design_matrix), times = replication), ]

# run the script with multi-core system
no_cores <- detectCores() - 1
c1 <- makePSOCKcluster(no_cores)
registerDoParallel(c1)
# to store the final results (including congruence between the simulation loading matrics, ARI of iCluster(which could already been obtained here),
# ARI of iCluster, and the congruence between obtained loading matrics and simulatd loading matrices)
results_sim1_data1 <- foreach(i = 1:nrow(design_matrix_replication),
                              .packages = c("rARPACK", "psych", "irlba", "mclust",
                                            "combinat", "GPArotation",
                                            "tidyr", "ez", "dplyr",  "Rcpp",
                                            "RcppArmadillo", "iCluster",
                                            "foreach", "doParallel", "ClusterSSCA"), .combine=rbind) %dopar%{


                                              # set the specific values of the parameters
                                              ncluster <- design_matrix_replication$n_cluster[i]
                                              psparse <- design_matrix_replication$p_sparse[i]
                                              pnoise <- design_matrix_replication$p_noise[i]
                                              mean <- design_matrix_replication$mean_level[i]
                                              cong <- design_matrix_replication$cong_level[i]
                                              memequal <- design_matrix_replication$mem_equal[i]

                                              # the number of observations in each cluster
                                              if (memequal == 1){
                                                cluster_mem <- rep(50, ncluster)
                                              }
                                              if (memequal == 0){
                                                cluster_mem <- c(rep(50, (ncluster - 1)), 30)
                                              }


                                              # data simulation
                                              sim_data <- CSSCASimulation(ncluster, cluster_mem, nblock, ncom,
                                                                          distinct, nvar, psparse, pnoise, cong, cong, "both", mean)

                                              block_data <- sim_data[[1]]
                                              all_data <- sim_data[[2]]
                                              loadings <- sim_data[[4]]
                                              cluster_assign <- sim_data[[5]]

                                              # cluster analysis with iCluster and examine the clustering accuracy
                                              icluster_results <- iCluster2(block_data, ncluster)$clusters
                                              # r1 represents the similarities between the estimated partitons of iCluster and the true partitions
                                              r1 <- adjustedRandIndex(cluster_assign, icluster_results)

                                              # the average congruence between the loading matrices
                                              n <- 0
                                              sum <- 0
                                              for (j in 1:(ncluster - 1)){
                                                for (k in (j + 1):ncluster){
                                                  n <- n + 1
                                                  sum <- sum + TuckerCongruence(loadings[[j]], loadings[[k]])
                                                }
                                              }
                                              # r2 represents the average congruence between cluster-specific true loading matrices
                                              r2 <- sum / n


                                              # computation with CSSCA (i.e. same procedures as described in FixedCSSCA)
                                              n_observation <- nrow(all_data)

                                              results1 <- csca_cpp(all_data, nvar, nblock, n_total, ncluster, csca_times)

                                              # the two rational starts
                                              results1_cssca <- cssca_quick_cpp(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, results1$cluster_mem, 1/6)
                                              results2_cssca <- cssca_quick_cpp(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, icluster_results, 1/6)

                                              # Compare the results of the two rational starts to determine the number of semi-rational starts
                                              covariance_similarity <- adjustedRandIndex(results1$cluster_mem, results1_cssca$cluster_mem)

                                              # compare between icluster and cssca
                                              mean_similarity <- adjustedRandIndex(icluster_results, results2_cssca$cluster_mem)

                                              ######### CSCA part
                                              min_loss_our <- upper
                                              if (results1_cssca$loss < min_loss_our) {
                                                min_loss_our <- results1_cssca$loss
                                                global1 <- results1_cssca
                                              }

                                              ######### iclust part
                                              min_loss_iclust <- upper
                                              if (results2_cssca$loss < min_loss_iclust) {
                                                min_loss_iclust <- results2_cssca$loss
                                                global2 <- results2_cssca
                                              }

                                              #### weight scheme
                                              ## when covariance similarity < mean similarity, mean structure probably dominates the overall structure
                                              if (covariance_similarity < mean_similarity){
                                                mean <- 1
                                                n_partition_csca <- partition_small_times
                                                n_partition_iclust <- partition_big_times
                                              }
                                              ## when covariance similarity > mean similarity, m=covariance structure probably dominates the overall structure
                                              if ((covariance_similarity > mean_similarity) | (covariance_similarity = mean_similarity)){
                                                mean <- 0
                                                n_partition_csca <- partition_big_times
                                                n_partition_iclust <- partition_small_times
                                              }

                                              ### partition for the csca part
                                              partition_results_csca <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, results1$cluster_mem, 1/6, n_partition_csca, 3, 1/10)
                                              if (partition_results_csca$loss < min_loss_our){
                                                min_loss_our <- partition_results_csca$loss
                                                global1 <- partition_results_csca
                                              }

                                              #### partition for thr mclust part
                                              partition_results_iclust <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, icluster_results, 1/6, n_partition_iclust, 3, 1/10)
                                              if (partition_results_iclust$loss < min_loss_iclust){
                                                min_loss_iclust <- partition_results_iclust$loss
                                                global2 <- partition_results_iclust
                                              }


                                              if (min_loss_our < min_loss_iclust){
                                                global <- global1
                                              }  else {
                                                global <- global2
                                              }

                                              ## summary of the results
                                              est_loadings <- global$loadings
                                              est_loss <- global$loss
                                              est_cluster_mem <- global$cluster_mem

                                              # r3 represents the average congruence between the estimated loadings and true loadings
                                              r3 <- ListCongruence(loadings, est_loadings, ncluster)
                                              # r4 represents the similarities between the estimated cluster partition of CSSCA and the true cluster partition
                                              r4 <- adjustedRandIndex(est_cluster_mem, cluster_assign)
                                              values <- c(r1, r2, mean(r3), r4, est_loss)

                                              # save data
                                              out <- list(all_data = all_data, sim_loadings = loadings, sim_assign = cluster_assign, est = global, values = values)

                                              save(out, file=paste0(i, ".RData"))

                                            }
stopCluster(c1)


######################################################################################
############################ sim2: Model selection ###################################
######################################################################################
setwd("~/sim2")
set.seed(912)
## fixed parameters
nblock <- 2
ncom <- 2
distinct <- c(1, 1)

# factorial design (overall 48 conditions)
n_cluster <- c(2, 4)
p_sparse <- c(0.3, 0.7)
p_noise <- c(0.15, 0.3)
mean_level <- c(0.1, 0.5, 0.9)
cong_level <- c(0, 0.7)
dim <- c(0,1)

# the possible ranges of the values to be selected
cluster_range <- 1:7
sparse_range <- seq(0.2, 0.8, by = 0.1)

# additional information from the input
n_total <- sum(ncom, distinct)

# running information
partition_small_times <- 12 #the maximum number of replacement when using CSSCA in mean-dominant structure and (or) using iCluster in covariance-dominant structure (i.e. inconsistency of the method and the structure)
partition_big_times <- 25 #the maximum number of replacement when using CSSCA in mean-dominant structure and (or) using iCluster in covariance-dominant structure (i.e. consistency of the method and the structrue)
csca_times <- 25
upper <- 1e9

# specify the main directory, and the simulated datasets will be stored in the sub-directory
mainDir <- "~/sim2"
# matrix that indicates the specific parameter setting of each condition
design_matrix <- expand.grid(n_cluster = n_cluster, dim = dim,
                             p_sparse = p_sparse, p_noise = p_noise, mean_level = mean_level,
                             cong_level = cong_level)
results_design <- vector(length = nrow(design_matrix))

# simulate a total of 96 datasets that represent the combinations of different values of parameters
for (i in 1:nrow(design_matrix)){

  # creaete and check
  subDir <- paste0(i, "dataset")
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir, subDir))

  # first create the simulations
  # variables in specific settings
  ncluster <- design_matrix$n_cluster[i]
  psparse <- design_matrix$p_sparse[i]
  pnoise <- design_matrix$p_noise[i]
  mean <- design_matrix$mean_level[i]
  cong <- design_matrix$cong_level[i]
  dim.size <- design_matrix$dim[i]
  # small-size data
  if (dim.size == 0){
    nvar <- c(15, 15)
  }
  # large-size data
  if (dim.size == 1){
    nvar <- c(15, 50)
  }

  cluster_mem <- rep(40, ncluster)

  # data simulation
  sim_data <- CSSCASimulation(ncluster, cluster_mem, nblock, ncom,
                              distinct, nvar, psparse, pnoise, cong, cong, "both", mean)

  block_data <- sim_data[[1]]
  all_data <- sim_data[[2]]
  loadings <- sim_data[[4]]
  cluster_assign <- sim_data[[5]]

  # the average congruence between the loading matrices
  n <- 0
  sum <- 0
  for (j in 1:(ncluster - 1)){
    for (k in (j + 1):ncluster){
      n <- n + 1
      sum <- sum + TuckerCongruence(loadings[[j]], loadings[[k]])
    }
  }
  results_design[i] <- sum / n

  out <- list(all_data = all_data, block_data = block_data, sim_loadings = loadings, sim_assign = cluster_assign)
  save(out, file="sim.RData")
}

# each dataset will be analyzed in each combination of the number of clusters and the sparsity level
compute_matrix <- expand.grid(ncluster_select = cluster_range, psparse_select = sparse_range)
compute_matrix$ind <- 1:nrow(compute_matrix)
compute_matrix_replication <- cbind(compute_matrix, index = rep(1:nrow(design_matrix), each = nrow(compute_matrix)))

# data simulation with multi-core system
no_cores <- detectCores() - 1
c1 <- makePSOCKcluster(no_cores)
registerDoParallel(c1)
# to store the final results (including congruence between the simulation loading matrics, ARI of iCluster(which could already been obtained here),
# ARI of iCluster, and the congruence between obtained loading matrics and simulatd loading matrices)
foreach(i = 1:nrow(compute_matrix_replication),
                              .packages = c("rARPACK", "psych", "irlba", "mclust",
                                            "combinat", "GPArotation",
                                            "tidyr", "ez", "dplyr",  "Rcpp",
                                            "RcppArmadillo", "iCluster",
                                            "foreach", "doParallel", "ClusterSSCA")) %dopar%{

                                              # set working directory
                                              ## iter is the number of datasets
                                              ## iter.ind is the index of the estimation done to a specific dataset
                                              ## i is the overall index of the estimation
                                              iter <- compute_matrix_replication$index[i]
                                              iter.ind <- compute_matrix_replication$ind[i]
                                              subDir <- paste0(iter, "dataset")
                                              setwd(file.path(mainDir, subDir))

                                              # variables in specific settings
                                              ncluster <- compute_matrix_replication$ncluster_select[i]
                                              psparse <- compute_matrix_replication$psparse_select[i]
                                              pnoise <- design_matrix$p_noise[iter]
                                              mean <- design_matrix$mean_level[iter]
                                              cong <- design_matrix$cong_level[iter]
                                              dim.size <- design_matrix$dim[i]

                                              # small-size data
                                              if (dim.size == 0){
                                                nvar <- c(15, 15)
                                              }
                                              # large-size data
                                              if (dim.size == 1){
                                                nvar <- c(15, 50)
                                              }

                                              # load the data
                                              load(file="sim.RData")
                                              all_data <- out[[1]]
                                              block_data <- out[[2]]

                                              # additional information from the input
                                              sum_var <- sum(nvar)
                                              # create the structure of loading matrices
                                              distinct_index <- vector("numeric")
                                              all_var <- 0
                                              ini_var <- 0
                                              for (p in 1:nblock){
                                                distinct_index <- c(distinct_index, rep(p, distinct[p]))
                                                ini_var <- ini_var + nvar[p]
                                                all_var <- c(all_var, ini_var)
                                              }
                                              distinct_zeros <- vector("numeric")
                                              for (r.distinct in 1:sum(distinct)){
                                                distinct_zeros <- c(distinct_zeros, ((sum_var * (ncom + r.distinct - 1)) + all_var[distinct_index[r.distinct]] + 1): ((sum_var * (ncom + r.distinct - 1)) + all_var[(distinct_index[r.distinct] + 1)]))
                                              }

                                              # computation with CSSCA
                                                n_observation <- nrow(all_data)
                                                if(ncluster == 1){
                                                  ssca_results <- IntSparseSca_random_cpp(all_data, sum_var, nblock, n_total, distinct_zeros, psparse)
                                                  out <- list(loss = ssca_results$loss, est = ssca_results)
                                                  save(out, file=paste0(iter.ind, ".RData"))
                                                }
                                                if(ncluster != 1){

                                                  # two raational partition startes
                                                  icluster_results <- iCluster2(block_data, ncluster)$clusters
                                                  results1 <- csca_cpp(all_data, nvar, nblock, n_total, ncluster, csca_times)

                                                  # two results produced by inputting the rational starts
                                                  results1_cssca <- cssca_quick_cpp(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, results1$cluster_mem, 1/6)
                                                  results2_cssca <- cssca_quick_cpp(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, icluster_results, 1/6)

                                                  # Compare the results of the two rational starts to determine the number of semi-rational starts
                                                  covariance_similarity <- adjustedRandIndex(results1$cluster_mem, results1_cssca$cluster_mem)

                                                  # compare between icluster and cssca
                                                  mean_similarity <- adjustedRandIndex(icluster_results, results2_cssca$cluster_mem)

                                                  ## CSCA part
                                                  min_loss_our <- upper
                                                  if (results1_cssca$loss < min_loss_our) {
                                                    min_loss_our <- results1_cssca$loss
                                                    global1 <- results1_cssca
                                                  }

                                                  ## iclust part
                                                  min_loss_iclust <- upper
                                                  if (results2_cssca$loss < min_loss_iclust) {
                                                    min_loss_iclust <- results2_cssca$loss
                                                    global2 <- results2_cssca
                                                  }

                                                  ## weight scheme
                                                  ## when covariance similarity < mean similarity, mean structure probably dominates the overall structure
                                                  if (covariance_similarity < mean_similarity){
                                                    mean <- 1
                                                    n_partition_csca <- partition_small_times
                                                    n_partition_iclust <- partition_big_times
                                                  }
                                                  ## when covariance similarity > mean similarity, m=covariance structure probably dominates the overall structure
                                                  if ((covariance_similarity > mean_similarity) | (covariance_similarity = mean_similarity)){
                                                    mean <- 0
                                                    n_partition_csca <- partition_big_times
                                                    n_partition_iclust <- partition_small_times
                                                  }

                                                  ### partition for the csca part
                                                  partition_results_csca <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, results1$cluster_mem, 1/6, n_partition_csca, 3, 1/10)
                                                  if (partition_results_csca$loss < min_loss_our){
                                                    min_loss_our <- partition_results_csca$loss
                                                    global1 <- partition_results_csca
                                                  }

                                                  #### partition for thr mclust part
                                                  partition_results_iclust <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, icluster_results, 1/6, n_partition_iclust, 3, 1/10)
                                                  if (partition_results_iclust$loss < min_loss_iclust){
                                                    min_loss_iclust <- partition_results_iclust$loss
                                                    global2 <- partition_results_iclust
                                                  }


                                                  if (min_loss_our < min_loss_iclust){
                                                    global <- global1
                                                  }  else {
                                                    global <- global2
                                                  }

                                                  ## summary of the results
                                                  est_loadings <- global$loadings
                                                  est_loss <- global$loss
                                                  est_cluster_mem <- global$cluster_mem

                                                  # save data
                                                  out <- list(loss = est_loss, est = global, icluster_results = icluster_results)

                                                  save(out, file=paste0(iter.ind, ".RData"))
                                                }
                                              }


stopCluster(c1)

