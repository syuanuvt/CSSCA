#install.packages("ez")
# initial set-up
library("rARPACK")
library("psych")
library("iCluster")
library("r.jive")
library("irlba")
library("mclust")
#install.packages("combinat")
library("combinat")
#install.packages("GPArotation")
library("GPArotation")
library("tidyr")
#install.packages("ez")
library("ez")
library("dplyr")

Add <- function(inidata, cluster_assign, cluster_num, prob, p_noise){
  # inidata: the data generated via the covariance structure
  # cluster_assign: the cluster assignment of the original data points
  # cluster_num: the number of clusters
  # prob: the proportion of the variance of the mean structure in total variance 
  # p_noise: pencentage of noise
  n_row <- nrow(inidata)
  n_col <- ncol(inidata)
  
  mean_data <- matrix(nrow = n_row, ncol = n_col)
  mean_cluster <- matrix(nrow = cluster_num, ncol = n_col)
  
  for (i in 1:cluster_num){
    mean_cluster[i, ] <- runif(n_col, -1, 1)
    mean_data[which(cluster_assign == i), ] <- rep(mean_cluster[i, ], each = sum(cluster_assign == i))
  }
  
  finaldata <- matrix(nrow = n_row, ncol = n_col)
  finaldata <- inidata + t(sqrt(prob * (1 - p_noise) / apply(mean_data, 2, var)) * t(mean_data))
  
  noise <- matrix(rnorm(n_row * n_col, 0, sqrt(p_noise)), n_row, n_col)
  finaldata <- finaldata + noise
  return (finaldata)
}

exchangemem <- function(cluster_assign, ncluster, number){
  n <- length(cluster_assign)
  change <- sample(1:n, number)
  for (i in 1:number){
    record <- cluster_assign[change[i]]
    change_pool <- setdiff(1:ncluster, record)
    change_cluster <- sample(change_pool, 1)
    cluster_assign[change[i]] <- change_cluster
  }
  return(cluster_assign)
}

paircongruence <- function(lista, listb, cluster_num){
  
  cong <- vector("numeric", length = cluster_num)
  for (i in 1:cluster_num){
    p <- 0
    cong_i <- vector("numeric", length = (cluster_num - 1))
    current <- lista[[i]]
    for (j in 1:cluster_num){
      p <- p + 1
      cong_i[p] <- tuckercongruence(current, listb[[j]])
    }
    cong[i] <- max(cong_i)
  }
  return (cong)
}

randomstart <- function(mem_cluster){
  #i: number of observations
  #ncluster: number of clusters
  #ncom: number of components
  #i <- 20
  #ncluster <- 4
  #set.seed(seed)
  
  i <- sum(mem_cluster)
  ncluster <- length(mem_cluster)
  
  mem.cluster <- matrix(0, nrow = i, ncol = ncluster)
  ind.cluster <- vector("numeric")
  
  # initialize
  for (t in 1:ncluster){
    ind.cluster <- c(ind.cluster, rep(t, mem_cluster[t]))
  }
  
  # randomlize
  ind.cluster <- sample(ind.cluster)
  
  # set to the membership cluster
  for (j in 1:i){
    mem.cluster[j, ind.cluster[j]] <- 1
  }
  
  return (list(mem = mem.cluster, ind = ind.cluster))
}

tuckercongruence <- function(matrix1, matrix2){
  ## note that two matrices should be in the same size
  # matrix1: reference matrix
  # matrix2: matrix to compare to the reference
  
  m <- nrow(matrix1)
  n <- ncol(matrix2)
  
  indic_perms <- permn(1:n)
  tuck <- vector("numeric")
  for (i in 1:length(indic_perms)){
    matrix2_perm <- matrix2[ ,indic_perms[[i]]]
    tuck[i] <- tr(abs(factor.congruence(matrix1, matrix2_perm))) / n
  }
  
  tuck_max <- max(tuck)
  return (tuck_max)
}

MatrixCenter <- function(matrix, center, scale){
  # matrix: the matrix that need to be centered and (or) scaled
  # center: does the matrix need to be centered (1 = Yes, 0 = No)
  # scaled: does the matrix need to be scaled (1 = Yes, 0 = No)
  
  variable_mean <- apply(matrix, 2, mean)
  variable_sd <- apply(matrix, 2, sd)
  n_observation <- nrow(matrix)
  
  if (center == 1){
    matrix <- matrix - rep(1, n_observation) %*% t(variable_mean)
  }
  if (scale == 1){
    matrix <- t(t(matrix) / variable_sd)
  }
  
  return(matrix)
}
