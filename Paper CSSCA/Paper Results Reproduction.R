## The script to reporduce the results in the manuscript
library("rARPACK")
library("psych")
library("irlba")
library("mclust")
library("combinat")
library("GPArotation")
library("tidyr")
library("dplyr")
library("mclust")
library("ggplot2")
library("Rcpp")
library("RcppArmadillo")
library("foreach")
library("doParallel")
library("ClusterSSCA")
library("apaTables")
library(MBESS)
library(reshape2)
library(jtools)
library("plotrix")

#########################################################################################
################################# Simulation Studies ####################################
##################################### Sim 1 #############################################
######## important note: Due to the hardware avilability, we have run the sim1 on two
######## separate PC, devided by the two levels of the factor clustersize
######## follows is the code when level = 1 (large-scale data); the exact same code
######## is used for another half of the analysis with the only difference that
######## number of observations in each cluster equals 50 or 30
#########################################################################################
setwd("~/sim1_data1")
## set a random seed
set.seed(912)
## fixed parameters
nblock <- 2
ncom <- 2
distinct <- c(1, 1)

## factorial designs
n_var <- c(0, 1)
n_cluster <- c(2, 4)
## 0 = unequal, 1 = equal
mem_equal <- c(0, 1)
p_sparse <- c(0.3, 0.5, 0.7)
p_noise <- c(0.1, 0.2, 0.3)
mean_level <- c(0.1, 0.5, 0.9)
cong_level <- c(0, 0.7)
## we design 10 replications in each condition
replication <- 40

# additional information from the input
n_total <- 4

# running information (i.e. the "difficult" computation in ClusterSSCA)
partition_times <- 20
csca_times <- 25
upper <- 1e9

# compress all paarmeters in the design matrix
design_matrix <- expand.grid(n_var = n_var, n_cluster = n_cluster, mem_equal = mem_equal,
                             p_sparse = p_sparse, p_noise = p_noise, mean_level = mean_level,
                             cong_level = cong_level)
design_matrix_replication <- design_matrix[rep(1:nrow(design_matrix), times = replication), ]

# run the script with multi-core system
no_cores <- 20
c1 <- makePSOCKcluster(no_cores)
registerDoParallel(c1)
## from iteration 1576, we start recording the time
# to store the final results (including congruence between the simulation loading matrics, ARI of iCluster(which could already been obtained here),
# ARI of iCluster, and the congruence between obtained loading matrics and simulatd loading matrices)
results_sim1_data1 <- foreach(i = 1:nrow(design_matrix_replication),
                              .packages = c( "psych",  "mclust",
                                             "combinat",
                                             "tidyr",  "dplyr",  "Rcpp",
                                             "RcppArmadillo", "iCluster",
                                             "foreach", "doParallel", "ClusterSSCA"), .combine=rbind) %dopar%{


                                               # set the specific values of the parameters
                                               largenvar <- design_matrix_replication$n_var[i]
                                               ncluster <- design_matrix_replication$n_cluster[i]
                                               psparse <- design_matrix_replication$p_sparse[i]
                                               pnoise <- design_matrix_replication$p_noise[i]
                                               mean <- design_matrix_replication$mean_level[i]
                                               cong <- design_matrix_replication$cong_level[i]
                                               memequal <- design_matrix_replication$mem_equal[i]

                                               # the number of variables
                                               # low dimensions
                                               if(largenvar == 0){
                                                 nvar <- c(15, 15)
                                               }
                                               # high dimensions
                                               if (largenvar == 1){
                                                 nvar <- c(15, 50)
                                               }

                                               # the number of observations in each cluster
                                               # in the cases that the cluster size is small, each cluster has 50 or 30 observations
                                               if (memequal == 1){
                                                 cluster_mem <- rep(100, ncluster)
                                               }
                                               if (memequal == 0){
                                                 cluster_mem <- c(rep(100, (ncluster - 1)), 60)
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


                                               start.time <- proc.time()
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
                                               min_loss <- upper
                                               if (results1_cssca$loss < min_loss) {
                                                 min_loss <- results1_cssca$loss
                                                 global <- results1_cssca
                                               }

                                               ######### iclust paat
                                               if (results2_cssca$loss < min_loss) {
                                                 min_loss <- results2_cssca$loss
                                                 global <- results2_cssca
                                               }

                                               #### weight scheme
                                               ## when covariance similarity < mean similarity, mean structure probably dominates the overall structure
                                               if (covariance_similarity < mean_similarity){
                                                 #### partition for thr iclust part
                                                 partition_results_iclust <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, icluster_results, 1/6, partition_times, 3, 1/10)
                                                 if (partition_results_iclust$loss < min_loss){
                                                   min_loss <- partition_results_iclust$loss
                                                   global <- partition_results_iclust
                                                 }
                                               }
                                               ## when covariance similarity > mean similarity, m=covariance structure probably dominates the overall structure
                                               if ((covariance_similarity > mean_similarity) | (covariance_similarity = mean_similarity)){
                                                 ### partition for the csca part
                                                 partition_results_csca <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, results1$cluster_mem, 1/6, partition_times, 3, 1/10)
                                                 if (partition_results_csca$loss < min_loss){
                                                   min_loss <- partition_results_csca$loss
                                                   global <- partition_results_csca
                                                 }
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

                                               dur.time <- proc.time() - start.time
                                               # save data
                                               out <- list(all_data = all_data, sim_loadings = loadings, sim_assign = cluster_assign, est = global, values = values, time = dur.time)

                                               save(out, file=paste0(i, ".RData"))

                                             }
stopCluster(c1)


###################################################################################################################

#setwd("~/sim1_data2")
nblock <- 2
n_var <- c(0,1)
ncom <- 2
distinct <- c(1, 1)
n_cluster <- c(2, 4)
mem_equal <- c(0, 1)   ## 0 = unequal, 1 = equal
p_sparse <- c(0.3, 0.5, 0.7)
p_noise <- c(0.1, 0.2, 0.3)
mean_level <- c(0.1, 0.5, 0.9)
cong_level <- c(0, 0.7)
cluster_size <- c(0,1)
replication <- 1:40

# sumamry table of paramters and results
design_matrix <- expand.grid(n_var = n_var, n_cluster = n_cluster, mem_equal = mem_equal,
                             p_sparse = p_sparse, p_noise = p_noise, mean_level = mean_level,
                             cong_level = cong_level, cluster_size = cluster_size)
design_matrix_replication <- expand.grid(n_var = n_var, n_cluster = n_cluster, mem_equal = mem_equal,
                                         p_sparse = p_sparse, p_noise = p_noise, mean_level = mean_level,
                                         cong_level = cong_level, replication = replication,
                                         cluster_size = cluster_size)
# 0 = small, 1 = large

design_matrix_full <- expand.grid(n_var = n_var, n_cluster = n_cluster, mem_equal = mem_equal,
                                  p_sparse = p_sparse, p_noise = p_noise, mean_level = mean_level,
                                  cong_level = cong_level, cluster_size = cluster_size)

# four additional columns to record the value of loss functions under iCLuster and CSSCA estimations, as well as the
# value of the average congruence between loading matrices and the average congruence of true loading matrices and estimated loading matrices
design_matrix_replication$icluster <- rep(0, nrow(design_matrix_replication))
design_matrix_replication$cssca <- rep(0, nrow(design_matrix_replication))
design_matrix_replication$congruence <- rep(0, nrow(design_matrix_replication))
design_matrix_replication$loading_cong <- rep(0, nrow(design_matrix_replication))
design_matrix_replication$dur_time <- rep(0, nrow(design_matrix_replication))

# if necessay, read in the two data
#setwd("~/sim1_data1")
for (i in 1:(nrow(design_matrix_replication) / 2)){
  load(paste0(i, ".RData"))
  values <- out$values
  # optimal loss function value of iCluster
  design_matrix_replication$icluster[i] <- values[1]
  # optimal loss fucntion value of CSSCA
  design_matrix_replication$cssca[i] <- values[4]
  # congruence between cluster-specific loading matrices
  design_matrix_replication$congruence[i] <- values[2]
  # the execution time
  design_matrix_replication$dur_time[i] <- out$time[1]
  # congruence between the true loading matrices and the estimated loading matrices
  true_loadings <- out$sim_loadings
  est_loadings <- out$est$loadings
  design_matrix_replication$loading_cong[i] <- mean(ListCongruence(true_loadings, est_loadings, as.numeric(as.character(design_matrix_replication$n_cluster[i]))))
}

#setwd("~/sim1_data2")
for (j in 1:((nrow(design_matrix_replication) / 2))){
  load(paste0(j, ".RData"))
  i <- j + 17280
  values <- out$values
  # optimal loss function value of iCluster
  design_matrix_replication$icluster[i] <- values[1]
  # optimal loss fucntion value of CSSCA
  design_matrix_replication$cssca[i] <- values[4]
  # congruence between cluster-specific loading matrices
  design_matrix_replication$congruence[i] <- values[2]
  # the execution time
  design_matrix_replication$dur_time[i] <- out$time[1]
  # congruence between the true loading matrices and the estimated loading matrices
  true_loadings <- out$sim_loadings
  est_loadings <- out$est$loadings
  design_matrix_replication$loading_cong[i] <- mean(ListCongruence(true_loadings, est_loadings, as.numeric(as.character(design_matrix_replication$n_cluster[i]))))
}
###############################################################################
##### Reproduce the results of the analysis (including figures)
###############################################################################
# create the factors for the variables
cols <- colnames(design_matrix)
design_matrix_replication[cols]<- sapply(design_matrix_replication[cols], function(x) as.factor(x))

## the actual congruence level as a function of the factor congruence
summary_congruence <- design_matrix_replication_factor %>%
  group_by(cong_level) %>%
  dplyr::summarise(congruence_mean = mean(congruence), congruence_sd = sd(congruence))

# reproduce the results of the first section (time demanding)
time <- design_matrix_replication_factor %>%
  dplyr::summarise(dur_time = mean(dur_time))

# test ARI over different factors
design_matrix_replication_factor %>%
  group_by(n_var) %>%
  dplyr::summarise(cssca.mean = mean(cssca), cssca.sd = sd(cssca))
design_matrix_replication_factor %>%
  group_by(n_cluster) %>%
  dplyr::summarise(cssca.mean = mean(cssca), cssca.sd = sd(cssca))
design_matrix_replication_factor %>%
  group_by(mem_equal) %>%
  dplyr::summarise(cssca.mean = mean(cssca), cssca.sd = sd(cssca))
design_matrix_replication_factor %>%
  group_by(p_sparse) %>%
  dplyr::summarise(cssca.mean = mean(cssca), cssca.sd = sd(cssca))
design_matrix_replication_factor %>%
  group_by(p_noise) %>%
  dplyr::summarise(cssca.mean = mean(cssca), cssca.sd = sd(cssca))
design_matrix_replication_factor %>%
  group_by(mean_level) %>%
  dplyr::summarise(cssca.mean = mean(cssca), cssca.sd = sd(cssca))
design_matrix_replication_factor %>%
  group_by(cong_level) %>%
  dplyr::summarise(cssca.mean = mean(cssca), cssca.sd = sd(cssca))
design_matrix_replication_factor %>%
  group_by(cluster_size) %>%
  dplyr::summarise(cssca.mean = mean(cssca), cssca.sd = sd(cssca))

# Reproduce Table 1: the summary results of ARI summarized by multiple variables
a <- design_matrix_replication %>%
  group_by(mean_level,p_noise) %>%
  dplyr::summarise(cssca.mean = mean(cssca), cssca.sd = sd(cssca))

# the average congruence between the original loading matrices and the estimated loading matrices
design_matrix_replication_factor %>%
  dplyr::summarise(cong.mean = mean(loading_cong), cong.sd = sd(loading_cong))

# Reproduce Table 2: the summary results of Tucker's Phi summarized by multiple variables
b <- design_matrix_replication %>%
  group_by(mean_level, p_noise) %>%
  dplyr::summarise(cong.mean = mean(loading_cong), cong.sd = sd(loading_cong))

# comparison between cssca and icluster (preparation for the plot)
compare_results_sim1_mean <- design_matrix_replication %>%
  group_by(mean_level) %>%
  dplyr::summarise(iCluster = mean(icluster), CSSCA = mean(cssca))

compare_results_sim1_se <- design_matrix_replication %>%
  group_by(mean_level) %>%
  dplyr::summarise(iCluster = std.error(icluster), CSSCA = std.error(cssca))

## plot the results
compare_long_mean <- melt(compare_results_sim1_mean, id.vars = "mean_level")
compare_long_se <- melt(compare_results_sim1_se, id.vars = "mean_level")
compare_long_mean <- cbind(compare_long_mean, compare_long_se$value)
colnames(compare_long_mean) <- c("b", "Methods", "average_ARI", "se")

cbPalette <- c("#999999", "#000000")
ggplot(aes(x = b, y = average_ARI, group = Methods, color = Methods), data = compare_long_mean) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = average_ARI - 1.96 * se, ymax = average_ARI + 1.96 * se), width = .1, size = .5, position = position_dodge(0.05)) +
  #scale_color_brewer(palette = "Paired") +
  labs(x = "Proportion of mean-level differences b", y = "Average ARI") +
  theme_apa(legend.pos = "right", legend.use.title = TRUE) +
  scale_colour_manual(values=cbPalette) +
  theme(text = element_text(size = 12))


  #########################################################################################
  ################################# Simulation 2###########################################
  #########################################################################################

####################################################################################
######################### simulation part ##########################################
####################################################################################

##############################################################################################
##############################################################################################
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
replication <- 1:4

# the possible ranges of the values to be selected
cluster_range <- 1:7
sparse_range <- seq(0.2, 0.8, by = 0.1)

# additional information from the input
n_total <- sum(ncom, distinct)

# running information
partition_times <- 20
csca_times <- 25
upper <- 1e9

# specify the main directory, and the simulated datasets will be stored in the sub-directory
mainDir <- "~/sim2"
# matrix that indicates the specific parameter setting of each condition
design_matrix <- expand.grid(n_cluster = n_cluster, dim = dim,
                             p_sparse = p_sparse, p_noise = p_noise, mean_level = mean_level,
                             cong_level = cong_level, replication = replication)
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

  cluster_mem <- rep(50, ncluster)

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
foreach(i = (nrow(compute_matrix_replication) / 2):nrow(compute_matrix_replication),
        .packages = c("rARPACK", "psych", "irlba", "mclust",
                      "combinat", "GPArotation",
                      "tidyr", "ez", "dplyr",  "Rcpp",
                      "RcppArmadillo", "iCluster",
                      "foreach", "doParallel", "ClusterSSCA")) %dopar%{

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


                          start.time <- proc.time()
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
                          min_loss <- upper
                          if (results1_cssca$loss < min_loss) {
                            min_loss <- results1_cssca$loss
                            global <- results1_cssca
                          }

                          ######### iclust paat
                          if (results2_cssca$loss < min_loss) {
                            min_loss <- results2_cssca$loss
                            global <- results2_cssca
                          }

                          #### weight scheme
                          ## when covariance similarity < mean similarity, mean structure probably dominates the overall structure
                          if (covariance_similarity < mean_similarity){
                            #### partition for thr iclust part
                            partition_results_iclust <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, icluster_results, 1/6, partition_times, 3, 1/10)
                            if (partition_results_iclust$loss < min_loss){
                              min_loss <- partition_results_iclust$loss
                              global <- partition_results_iclust
                            }
                          }
                          ## when covariance similarity > mean similarity, m=covariance structure probably dominates the overall structure
                          if ((covariance_similarity > mean_similarity) | (covariance_similarity = mean_similarity)){
                            ### partition for the csca part
                            partition_results_csca <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, results1$cluster_mem, 1/6, partition_times, 3, 1/10)
                            if (partition_results_csca$loss < min_loss){
                              min_loss <- partition_results_csca$loss
                              global <- partition_results_csca
                            }
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

####################################################################################
######################### Analysis of results of simulation 2 ##########################################
####################################################################################

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
  replication <- 1:6

  # the possible ranges of the values to be selected
  cluster_range <- 1:7
  sparse_range <- seq(0.2, 0.8, by = 0.1)

  # additional information from the input
  n_total <- sum(ncom, distinct)

  # running information
  partition_times <- 20
  csca_times <- 25
  upper <- 1e9

  # specify the main directory, and the simulated datasets will be stored in the sub-directory
  mainDir <- "~/sim2"
  # matrix that indicates the specific parameter setting of each condition
  design_matrix <- expand.grid(n_cluster = n_cluster, dim = dim,
                               p_sparse = p_sparse, p_noise = p_noise, mean_level = mean_level,
                               cong_level = cong_level, replication = replication)

  # the matrix that records th results of the estimation
  compute_matrix <- expand.grid(ncluster_select = cluster_range, psparse_select = sparse_range)
  compute_matrix$loss <- rep(0, nrow(compute_matrix))
  compute_matrix$vaf <- rep(0, nrow(compute_matrix))

  # we will store the results of each condition in a list
  result_matrix <- list()
  ## restore the loss function values in the new matrix
  for (i in 1:nrow(design_matrix)){
    # set the working directory at the specific sub-directory (the analysis on one specific dataset)
    result_matrix[[i]] <- compute_matrix
    subDir <- paste0(i, "dataset")
    setwd(file.path(mainDir, subDir))
    for (j in 1:nrow(compute_matrix)){
      load(paste0(j, ".RData"))
      result_matrix[[i]]$loss[j] <- out$loss
      result_matrix[[i]]$vaf[j] <- mean(sapply(out$est$loadings, function(x) sum(x^2)))
    }
  }
  ## the matrix store the model selection result for each dataset created
  result_matrix_summary <- design_matrix
  ## the selected number of cluster
  result_matrix_summary$cluster_select <- rep(0, nrow(design_matrix))
  ## the selected level of sparsity
  result_matrix_summary$sparse_select <- rep(0, nrow(design_matrix))

  for (k in 1:nrow(design_matrix)){
    dataset1 <- result_matrix[[k]]
    # copy the dataset
    a <- dataset1
    sr.result <- matrix(nrow = length(sparse_range), ncol = (length(cluster_range) - 2))
    # when the algorithm does not end normally, the result is set at NA
    a$loss[which(a$loss == upper)] <- NA
    ## first select the optimal number of clusters for each dataset,
    ## following the model selection procedures described in the paper
    ## with the comparison of the average scree ratios across all possible values of level of sparsity
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
    ## Conditional on the optimal value of the nuber of cluster,
    ## select the optimal number of sparsity for each dataset
    filter.data <- dataset1 %>%
      dplyr::filter(ncluster_select == opt_cluster) %>%
      dplyr::select(loss)
    filter.data <- as.matrix(filter.data)
    sr.vector <- rep(0, (length(sparse_range)-2))
    for (j in 1:(length(sparse_range)-2)){
      sr.vector[j] <- (filter.data[j+2] - filter.data[j+ 1]) / (filter.data[j+1] - filter.data[j])
    }
    if (all(diff(filter.data) >= 0)){
      opt_psparse <- sparse_range[which(sr.vector == max(sr.vector))+1]
    }
    # sr < 0 rarely happens
    if (sum(diff(filter.data) < 0) > 0){
      opt_psparse <- sparse_range[which(diff(filter.data) < 0)[1] + 1]
    }
    ## restore the model selection results for various datasets
    result_matrix_summary$cluster_select[k] <- opt_cluster
    result_matrix_summary$sparse_select[k] <- opt_psparse
  }

  ## create loogical values to indicate whether the model selection procedure selects the true model parameter
  result_matrix_summary$cluster_cor <- result_matrix_summary$n_cluster == result_matrix_summary$cluster_select
  result_matrix_summary$sparse_cor <- result_matrix_summary$p_sparse == result_matrix_summary$sparse_select

  ## analysis on the predictive powewrs of various factors
  both <- sum(result_matrix_summary$cluster_cor & result_matrix_summary$sparse_cor)
  both / 576
  only.cluster <- sum(result_matrix_summary$cluster_cor & !result_matrix_summary$sparse_cor)
  only.cluster
  only.cluster / 576
  only.sparse <- sum(!result_matrix_summary$cluster_cor & result_matrix_summary$sparse_cor)
  only.sparse
  only.sparse / 576

  # Check the performance of succesful selections of paramters as function of various parameters
  result_matrix_summary %>%
    group_by(mean_level) %>%
    summarize(a = sum(cluster_cor & sparse_cor))

##########################################################################
############## Create the dataset that is to be analyzed #################
##########################################################################
###########################################################################
### read-in the data
load("~/application.RData")
data.app <- data.app[,3:ncol(data.app)]
## set up the parameters
nblock <- 2 # number of blocks
nvar <- c(6, 18) # number of variables in each data block
ncom <- 2 # number of common components
distinct <- c(1,1) # number of distinct components
n_total <- ncom+sum(distinct)
sum_var <- sum(nvar)
p_sparse_select <- seq(0,0.9,0.1) # possible values of sparsity
n_cluster_select <- 1:8

# running information
partition_times <- 20 #the maximum number of replacement when using CSSCA in mean-dominant structure and (or) using iCluster in covariance-dominant structure (i.e. inconsistency of the method and the structure)
csca_times <- 25
upper <- 1e9

# create the structure of loading matrices
distinct_index <- vector("numeric")
all_var <- 0
# useless onwards
ini_var <- 0
# p: the number of block
for (p in 1:nblock){
  distinct_index <- c(distinct_index, rep(p, distinct[p]))
  ini_var <- ini_var + nvar[p]
  all_var <- c(all_var, ini_var)
}
distinct_zeros <- vector("numeric")
# r.distinct: the number of distinctive component
for (r.distinct in 1:sum(distinct)){
  distinct_zeros <- c(distinct_zeros, ((sum_var * (ncom + r.distinct - 1)) + all_var[distinct_index[r.distinct]] + 1): ((sum_var * (ncom + r.distinct - 1)) + all_var[(distinct_index[r.distinct] + 1)]))
}

design_matrix <-expand.grid(ncluster_select = n_cluster_select, psparse_select = p_sparse_select)
#design_matrix$ind <- 1:nrow(design_matrix)
setwd("C:/Users/u1275970/Documents/CSSCA/NewSimulation/application1")
all_data <- data.app
block_data <- list(all_data[,1:6], all_data[,7:24])

no_cores <- detectCores() - 1
c1 <- makePSOCKcluster(no_cores)
registerDoParallel(c1)
# to store the final results (including congruence between the simulation loading matrics, ARI of iCluster(which could already been obtained here),
# ARI of iCluster, and the congruence between obtained loading matrics and simulatd loading matrices)
foreach(i = 1:nrow(design_matrix),
        .packages = c("rARPACK", "psych", "irlba", "mclust",
                      "combinat", "GPArotation",
                      "tidyr", "ez", "dplyr",  "Rcpp",
                      "RcppArmadillo", "iCluster",
                      "foreach", "doParallel", "ClusterSSCA")) %dopar%{

                        # variables in specific settings
                        ncluster <- design_matrix$ncluster_select[i]
                        psparse <- design_matrix$psparse_select[i]


                        #if(!file.exists(paste0(i, ".RData"))){
                          # computation with CSSCA
                          n_observation <- nrow(all_data)
                          if(ncluster == 1){
                            ssca_results <- IntSparseSca_random_cpp(all_data, sum_var, nblock, n_total, distinct_zeros, psparse)
                            out <- list(loss = ssca_results$loss, est = ssca_results)
                            ### the output are t,p,l respectively
                            save(out, file=paste0(i, ".RData"))
                          }
                          if(ncluster != 1){

                            icluster_results <- iCluster2(block_data, ncluster)$clusters
                            results1 <- csca_cpp(all_data, nvar, nblock, n_total, ncluster, csca_times)

                            # the two rational starts
                            results1_cssca <- cssca_quick_cpp(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, results1$cluster_mem, 1/6)
                            results2_cssca <- cssca_quick_cpp(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, icluster_results, 1/6)

                            # Compare the results of the two rational starts to determine the number of semi-rational starts
                            covariance_similarity <- adjustedRandIndex(results1$cluster_mem, results1_cssca$cluster_mem)

                            # compare between icluster and cssca
                            mean_similarity <- adjustedRandIndex(icluster_results, results2_cssca$cluster_mem)

                            min_loss <- upper
                            ######### CSCA part
                            if (results1_cssca$loss < min_loss) {
                              min_loss <- results1_cssca$loss
                              global <- results1_cssca
                            }

                            ######### iclust paat
                            if (results2_cssca$loss < min_loss) {
                              min_loss <- results2_cssca$loss
                              global <- results2_cssca
                            }
                            #### weight scheme
                            ## when covariance similarity < mean similarity, mean structure probably dominates the overall structure
                            if (covariance_similarity < mean_similarity){
                              #### partition for thr iclust part
                              partition_results_iclust <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, icluster_results, 1/6, partition_times, 3, 1/10)
                              if (partition_results_iclust$loss < min_loss){
                                min_loss <- partition_results_iclust$loss
                                global <- partition_results_iclust
                              }
                            }
                            ## when covariance similarity > mean similarity, m=covariance structure probably dominates the overall structure
                            if ((covariance_similarity > mean_similarity) | (covariance_similarity = mean_similarity)){
                              ### partition for the csca part
                              partition_results_csca <- MainCSSCA(all_data, nvar, nblock, ncom, distinct, ncluster, nrow(all_data), psparse, results1$cluster_mem, 1/6, partition_times, 3, 1/10)
                              if (partition_results_csca$loss < min_loss){
                                min_loss <- partition_results_csca$loss
                                global <- partition_results_csca
                              }
                            }


                            ## summary of the results
                            est_loadings <- global$loadings
                            est_loss <- global$loss
                            est_cluster_mem <- global$cluster_mem

                            # save data
                            out <- list(loss = est_loss, est = global, icluster_results = icluster_results)

                            save(out, file=paste0(i, ".RData"))
                          }
                        }


                      #}
stopCluster(c1)


###########################################################################
#################### Analysis on the empirical results  ###################
###########################################################################
setwd("C:/Users/u1275970/Documents/CSSCA/NewSimulation/application1")
design_matrix$loss <- rep(0, nrow(design_matrix))
#head(self.insight.imputed)
for (i in 1:nrow(design_matrix)){
  n_cluster <- design_matrix$ncluster_select[i]
  load(paste0(i, ".RData"))
  if (n_cluster == 1){
    design_matrix$loss[i] <- out[[1]]
  }
  else{
    design_matrix$loss[i] <- out[[1]]
  }
}

## select the optimal number of cluster based on the average scree ratio
a <- design_matrix
sr.result <- matrix(nrow = length(p_sparse_select), ncol = (length(n_cluster_select) - 2))
#a$loss[which(a$loss == 0)] <- NA
for (i in 1:length(p_sparse_select)){
  filter.data <- a %>%
    dplyr::filter(psparse_select == p_sparse_select[i]) %>%
    dplyr::select(loss)
  filter.data <- as.matrix(filter.data)
  for (j in 1:ncol(sr.result)){
    if (!is.na(filter.data[j]) & !is.na(filter.data[j + 1]) & !is.na(filter.data[j+2]))
      sr.result[i, j] <- (filter.data[j] - filter.data[j+ 1]) / (filter.data[j+1] - filter.data[j+ 2])
  }
}
sum <- apply(sr.result, 2, sum, na.rm = TRUE)
opt_cluster <- n_cluster_select[which(sum == max(sum)) + 1]

## select the optimal level of sparsity, conditional on the optimal number of cluster
filter.data <- a %>%
  dplyr::filter(ncluster_select == opt_cluster) %>%
  dplyr::select(loss)
filter.data <- as.matrix(filter.data)
#filter.data$loss[which(filter.data$loss == 0)] <- 10000
sr.vector <- rep(0, (length(p_sparse_select)-2))
#### select the number of clusters based on non-sparse solution
for (j in 1:(length(p_sparse_select)-2)){
  sr.vector[j] <- (filter.data[j+2] - filter.data[j+ 1]) / (filter.data[j+1] - filter.data[j])
}

opt_psparse <- p_sparse_select[which(sr.vector == max(sr.vector))+1]

# plot for the optimal clusters
plot.data <- filter.data <- a %>%
  dplyr::filter(ncluster_select == opt_cluster) %>%
  dplyr::select(loss, psparse_select)
qplot(x = psparse_select, y = loss, data = plot.data, main = "Selection of level of sparsity") +
  #geom_text(data = plot_material1, mapping = aes(x = complexity, y = fit, label = Rnames)) +
  #geom_line(data = plot_material1, mapping = aes(x = complexity, y = fit), colour = "red") +
  geom_point(data = plot.data, mapping = aes(x = psparse_select, y = loss), colour = "blue")

## plot the heatmap
number <- which(design_matrix$ncluster_select == 3 & design_matrix$psparse_select == 0.4)
load(paste0(35, ".RData"))
#ssca_results
results <- out[[2]]
cluster.assignment <- results$cluster_mem
loss <- results$loss
loadings <- results$loadings
bk <- seq(-0.5,0.5,by=.1)
mycols <- colorRampPalette(colors = c("red", "white","blue"))(length(bk)-1)
row.lab <- colnames(self.insight.imputed)
column.lab <- c("com1", "com2","com3", "dis1", "dis2")
# the loadings matrix of each cluster

heatmap(abs(loadings[[1]]), Colv = NA, Rowv = NA, col = mycols, breaks = bk, scale = "none", labRow = row.lab, labCol = column.lab)
heatmap(abs(loadings[[2]]), Colv = NA, Rowv = NA, col = mycols, breaks = bk, scale = "none", labRow = row.lab, labCol = column.lab)
#heatmap(loadings[[3]], Colv = NA, Rowv = NA, col = mycols, breaks = bk, scale = "none", labRow = row.lab, labCol = column.lab)

