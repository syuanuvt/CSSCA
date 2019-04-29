#' Exchange the cluster assignment in the cluster partition vector. In the context of the CSSCA method, the function is
#' esecially useful to add some randomness into the output of preliminary clustering analysis (e.g. Mcluster, iCluster and the non-sparse version of the CSSCA method)
#'
#' @param cluster_assign the original vector of cluster assignment. The input vector should have \code{N} elements, and
#' the \code{nth} element reresents the assigned cluster of the \code{nth} entry
#' @param ncluster A non-negative integer represents the amount of clusters
#' @param number A non-negative integer indicates the change of cluster memberships should apply to how many entries
#' @return the new cluster assignment as indicatd by a new cluster partition vector with the length of \code{N} (the number of observations)
#' @examples
#'    #assume 6 entries and 3 clusters. The membrships of 3 entries should be changed
#'    cluser_ass <- c(1,1,2,2,3,3)
#'    n_cluster <- 3
#'    exe_number <- 3
#'  (not run) exchangemem(cluster_ass, n_cluster, exe_number)
#'
ExchangeMem <- function(cluster_assign, ncluster, number){
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
