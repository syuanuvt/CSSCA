#' In the simulations, add the mean structure and the noise structure to the original co-variance stucture
#'
#' @param inidata The generated data that includes cluster differences in co-variance structure (i.e. loading matices).
#' The input should be a matrix with the number of rows should be equal to the number of entries \code{N}
#' while the number of columns should be equal to the total amount of variables \code{p} (over all blocks).
#' @param cluster_num  A non-negive interger indicates the total number of clusters that would be generated.
#' @param cluster_assign  A vector indicates the assignment of each entry to one of the cluster. The vector
#' should be of length \code{N}, and the nth element refers to the number of cluster that the nth entry belong to.
#' @param prob A number within the range of [0,1] that indicates the proportion of mean-level differences
#' in the total cluster difference; the proportion of total difference that is accounted by co-variance difference
#' would therefore be \code(1 - prob).
#' @param p_noise A number within the range of [0,1] that indicates the percentage of noise structrue that should be added to the final data.
#' @return a list of two elements. The first element is the generated final data and the second element is the noise structure
#' @examples
#' m <- matrix(1:6,2,3) #N = 2
#' n_cluster <- 2
#' cluster_assignment <- c(2,1)
#' p_mean <- 0.9 #90% mean structure
#' p_noise <- 0.1 # 10% noise
#' Add(m, n_cluster, cluster_assignment, p_mean, p_noise)
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
  return (list(final_data = finaldata, noise = noise))
}
