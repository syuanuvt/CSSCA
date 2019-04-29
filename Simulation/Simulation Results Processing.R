
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
#########################################################################################
################################# Simulation 1###########################################
#########################################################################################

setwd("~/sim1_data1")
## set a random seed
set.seed(912)
nblock <- 2
nvar <- c(15, 15)
ncom <- 2
distinct <- c(1, 1)
n_cluster <- c(2, 4)
mem_equal <- c(0, 1)   ## 0 = unequal, 1 = equal
p_sparse <- c(0.3, 0.5, 0.7)
p_noise <- c(0.1, 0.2, 0.3)
mean_level <- c(0.1, 0.5, 0.9)
cong_level <- c(0, 0.7)
replication <- 10

# sumamry table of paramters and results
design_matrix <- expand.grid(n_cluster = n_cluster, mem_equal = mem_equal,
                             p_sparse = p_sparse, p_noise = p_noise, mean_level = mean_level,
                             cong_level = cong_level)
design_matrix_replication <- design_matrix[rep(1:nrow(design_matrix), times = replication),]

# four additional columns to record the value of loss functions under iCLuster and CSSCA estimations, as well as the
# value of the average congruence between loading matrices and the average congruence of true loading matrices and estimated loading matrices
design_matrix_replication$icluster <- rep(0, nrow(design_matrix_replication))
design_matrix_replication$cssca <- rep(0, nrow(design_matrix_replication))
design_matrix_replication$congruence <- rep(0, nrow(design_matrix_replication))
design_matrix_replication$loading_cong <- rep(0, nrow(design_matrix_replication))
for (i in 1:nrow(design_matrix_replication)){
  load(paste0(i, ".RData"))
  values <- out$values
  # optimal loss function value of iCluster
  design_matrix_replication$icluster[i] <- values[1]
  # optimal loss fucntion value of CSSCA
  design_matrix_replication$cssca[i] <- values[4]
  # congruence between cluster-specific loading matrices
  design_matrix_replication$congruence[i] <- values[2]
  # congruence between the true loading matrices and the estimated loading matrices
  true_loadings <- out$sim_loadings
  est_loadings <- out$est$loadings
  design_matrix_replication2$loading_cong[i] <- mean(ListCongruence(true_loadings, est_loadings, as.numeric(as.character(simulation1_results$n_cluster[i]))))
}

##  integrate the results of the second part of simulation 1
setwd("~/sim1_data2")
design_matrix_replication2 <- design_matrix_replication
for (i in 1:nrow(design_matrix_replication2)){
  load(paste0(i, ".RData"))
  values <- out$values
  design_matrix_replication2$icluster[i] <- values[1]
  design_matrix_replication2$cssca[i] <- values[4]
  design_matrix_replication2$congruence[i] <- values[2]
  true_loadings <- out$sim_loadings
  est_loadings <- out$est$loadings
  design_matrix_replication2$loading_cong[i] <- mean(ListCongruence(true_loadings, est_loadings, as.numeric(as.character(simulation1_results$n_cluster[i]))))
}

## aggregate the two datasets together
design_matrix_replication$dimension <- rep("low", nrow(design_matrix_replication))
design_matrix_replication2$dimension <- rep("high", nrow(design_matrix_replication2))
simulation1_results <- as.data.frame(rbind(design_matrix_replication, design_matrix_replication2))
# create the factors for the variables
cols <- colnames(design_matrix)
simulation1_results[cols] <- lapply(simulation1_results[cols], factor)

# Reproduce table 1: regression analysis of the clustering accuracy
simulation1_reg <- lm(cssca ~ n_cluster * mem_equal * p_sparse * p_noise * mean_level * cong_level * dimension, data = simulation1_results)
# summary of the regression table
summary(simulation1_reg)
# the effect sizes of the main and interaction effects
omega_sq(simulation1_reg)

# Reproduce Table 2: the summary results summarized by multiple variables
summary_results_sim1 <- simulation1_results %>%
  group_by(mean_level, dimension, p_noise) %>%
  dplyr::summarise(cssca.mean = mean(cssca), cssca.sd = sd(cssca))

# Reproduce Figure 3: comparison between cssca and icluster
compare_results_sim1_mean <- simulation1_results %>%
  group_by(mean_level) %>%
  dplyr::summarise(icluster = mean(icluster), cssca = mean(cssca))
compare_results_sim1_se <- simulation1_results %>%
  group_by(mean_level) %>%
  dplyr::summarise(icluster = std.error(icluster), cssca = std.error(cssca))
# prepare for plotting
compare_long_mean <- melt(compare_results_sim1_mean, id.vars = "mean_level")
compare_long_se <- melt(compare_results_sim1_se, id.vars = "mean_level")
compare_long_mean <- cbind(compare_long_mean, compare_long_se$value)
colnames(compare_long_mean) <- c("b", "Methods", "average_ARI", "se")
# actual plotting
ggplot(aes(x = b, y = average_ARI, group = Methods, color = Methods), data = compare_long_mean) +
  geom_point() +
  geom_line(aes(linetype=Methods)) +
  geom_errorbar(aes(ymin = average_ARI - se, ymax = average_ARI + se), width = .1, size = .5, position = position_dodge(0.05)) +
  labs(x = "Proportion of mean-level cluster differences", y = "Average ARI") +
  theme_apa() +
  theme(text = element_text(size = 12)) +
  scale_color_manual(values=c( "black", "grey"))

###### analysis on the congruence of the resulting loading matrices and the true loading matrices

# Reproduce table 1: regression analysis of the congruence between loading matrices
simulation1_reg_cong <- lm(loading_cong ~ n_cluster * mem_equal * p_sparse * p_noise * mean_level * cong_level * dimension, data = simulation1_results)
summary(simulation1_reg_cong)
omega_sq(simulation1_reg_cong)
# some additional results reported in the paper
results_cong <- simulation1_results %>%
  dplyr::summarise(cong = mean(loading_cong))
compare_results_sim1_cong <- simulation1_results %>%
  group_by(mean_level, p_noise) %>%
  dplyr::summarise(cong = mean(loading_cong), cong_sd = sd(loading_cong))


#########################################################################################
################################# Simulation 2###########################################
#########################################################################################
setwd("~/sim2")
set.seed(912)
## fixed parameters
nblock <- 2
ncom <- 2
distinct <- c(1, 1)
upper <- 1e9

# fatorial design (overall 48 conditions)
n_cluster <- c(2, 4)
p_sparse <- c(0.3, 0.7)
p_noise <- c(0.15, 0.3)
mean_level <- c(0.1, 0.5, 0.9)
cong_level <- c(0, 0.7)
dim <- c(0,1)

design_matrix <- expand.grid(n_cluster = n_cluster, dim = dim,
                             p_sparse = p_sparse, p_noise = p_noise, mean_level = mean_level,
                             cong_level = cong_level)

# the possible ranges of the values to be selected
cluster_range <- 1:7
sparse_range <- seq(0.2, 0.8, by = 0.1)

# the matrix that records th results of the estimation
compute_matrix <- expand.grid(ncluster_select = cluster_range, psparse_select = sparse_range)
compute_matrix$loss <- rep(0, nrow(compute_matrix))

mainDir <- "~/sim2_2"

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
  filter.data <- dataset1 %>%
    dplyr::filter(ncluster_select == opt_cluster) %>%
    dplyr::select(loss)
  filter.data <- as.matrix(filter.data)
  sr.vector <- rep(0, (length(sparse_range)-2))
  for (j in 1:(length(sparse_range)-2)){
    sr.vector[j] <- (filter.data[j+2] - filter.data[j+ 1]) / (filter.data[j+1] - filter.data[j])
  }
  opt_psparse <- sparse_range[which(sr.vector == max(sr.vector))+1]

  ## restore the model selection results for various datasets
  result_matrix_summary$cluster_select[k] <- opt_cluster
  result_matrix_summary$sparse_select[k] <- opt_psparse
}

## create loogical values to indicate whether the model selection procedure selects the true model parameter
result_matrix_summary$cluster_cor <- result_matrix_summary$n_cluster == result_matrix_summary$cluster_select
result_matrix_summary$sparse_cor <- result_matrix_summary$p_sparse == result_matrix_summary$sparse_select

## create various values that are reported in the paper (following are only a few examples)
# variable selection is correct for at least one of the parameter
sum(sim2$cluster_cor | sim2$sparse_cor)
# variable seelction is correct for both parameters
sum(sim2$cluster_cor & sim2$sparse_cor)
# the times that variable selection is succesful in various levels of noise
sim2 %>%
  group_by(p_noise) %>%
  summarize(a = sum(cluster_cor & sparse_cor))


