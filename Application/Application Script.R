## Analysis of application created for the paper CSSCA
## date of last edits : 31-03-2018
## the personality data used in the current paper could be request at s.yuan@uvt.nl
##############################################################################################

library("rARPACK")
library("psych")
library("irlba")
library("mclust")
library("combinat")
library("GPArotation")
library("tidyr")
library("ez")
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
library(ClusterSSCA)
library(mice)
library(psych)
library(factoextra)

##################################################################################
############################# Data Preprocessing #################################
##################################################################################
## read in the personality and non-verbal behavior data
load("app4.RData")
app4 <- read.table("selfinsight_app3.csv", header = TRUE, sep = ",")
col.names <- colnames(app4)
app4 <- app4[,1:35]
# missing data imputation
imp <- mice(app4, m = 20, maxit = 10, seed = 912)
app4.imputed <- mice::complete(imp)
# only retain the relevant variables
app4.imputed <- app4.imputed[, -c(9, 13, 17:18, 22, 24:25, 33:34)]
# re-scale the original data to have the sum-of-squares of each column(variable) equals to 1 and the column mean equals to 0
app4.imputed <- apply(app4.imputed, 2, function(x) x / sqrt(sum(x ^ 2)))
save(app4.imputed, file="app4.RData")

##################################################################################
############################# Data estimation  ###################################
##################################################################################
## set up the parameters
nblock <- 2 # number of blocks
nvar <- c(8, 18) # number of variables in each data block
ncom <- 2 # number of common components
distinct <- c(1,1) # number of distinct components
n_total <- sum(n_com, n_distinct)
sum_var <- sum(nvar)
p_sparse_select <- seq(0,0.9,0.1) # possible values of sparsity
n_cluster_select <- 1:8

# running information
partition_small_times <- 12 #the maximum number of replacement when using CSSCA in mean-dominant structure and (or) using iCluster in covariance-dominant structure (i.e. inconsistency of the method and the structure)
partition_big_times <- 25 #the maximum number of replacement when using CSSCA in mean-dominant structure and (or) using iCluster in covariance-dominant structure (i.e. consistency of the method and the structrue)
csca_times <- 25
upper <- 1e9

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

# the matrix represents the possible combinations of differnt values of the model parameters
design_matrix <-expand.grid(ncluster_select = n_cluster_select, psparse_select = p_sparse_select)
all_data <- app4.imputed
block_data <- list()
block_data[[1]] <- all_data[,1:8]
block_data[[2]] <- all_data[,9:26]

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

                          # computation with CSSCA
                          n_observation <- nrow(all_data)
                          if(ncluster == 1){
                            ssca_results <- IntSparseSca_random_cpp(all_data, sum_var, nblock, n_total, distinct_zeros, psparse)
                            out <- list(loss = ssca_results$loss, est = ssca_results)
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

                            # save data
                            out <- list(loss = est_loss, est = global, icluster_results = icluster_results)

                            save(out, file=paste0(i, ".RData"))
                          }
                        }

stopCluster(c1)

##################################################################################
############################# Estimation Summary  ################################
##################################################################################
design_matrix <-expand.grid(ncluster_select = n_cluster_select, psparse_select = p_sparse_select)
# extend the design matrix to include one row that restores the loss function obtained in each combination of the parameter
design_matrix$loss <- rep(0, nrow(design_matrix))
# obtain the loss function from the stored data and restore it into the matrix
for (i in 1:nrow(design_matrix)){
  n_cluster <- design_matrix$ncluster_select[i]
  load(paste0(i, ".RData"))
  design_matrix$loss[i] <- out[[1]]
}

### Select the optimal number of cluster
# a copy of the design matrix
a <- design_matrix
sr.result <- matrix(nrow = length(p_sparse_select), ncol = (length(n_cluster_select) - 2))
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

## select the optimal level of sparsity
filter.data <- a %>%
  dplyr::filter(ncluster_select == opt_cluster) %>%
  dplyr::select(loss)
filter.data <- as.matrix(filter.data)
sr.vector <- rep(0, (length(p_sparse_select)-2))
for (j in 1:(length(p_sparse_select)-2)){
  sr.vector[j] <- (filter.data[j+2] - filter.data[j+ 1]) / (filter.data[j+1] - filter.data[j])
}
opt_psparse <- p_sparse_select[which(sr.vector == max(sr.vector))+1]

## plot the heatmap of the loading matrix
number <- which(design_matrix$ncluster_select == opt_cluster & design_matrix$psparse_select == opt_psparse)
load(paste0(number, ".RData"))
# the full estimation results
results <- out[[2]]
cluster.assignment <- results$cluster_mem
loss <- results$loss
loadings <- results$loadings
bk <- seq(-1,1,by=.01)
mycols <- colorRampPalette(colors = c("red", "white","blue"))(length(bk)-1)
row.lab <- colnames(app4.imputed)
column.lab <- c("com1", "com2","dis1", "dis2")
# the loadings matrix of each cluster
heatmap(loadings[[1]], Colv = NA, Rowv = NA, col = mycols, breaks = bk, scale = "none", labRow = row.lab, labCol = column.lab)
heatmap(loadings[[2]], Colv = NA, Rowv = NA, col = mycols, breaks = bk, scale = "none", labRow = row.lab, labCol = column.lab)

# number of observations in each cluster
cluster.mem <- out[[2]]$cluster_mem
n.cluster1 <- sum(cluster.mem == 1)
n.cluster2 <- sum(cluster.mem == 2)

