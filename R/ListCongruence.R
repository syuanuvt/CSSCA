#' Calculate the Tucker Congruence of two lists of matrices
#'
#' @param lista the first list with matrix elements
#' @param listb the second list with matrix elements
#' @param cluster_num the number of cluster (which should also be the maximum length of the two lists)
#' @return a vector of length cluster_num. \code{nth} elemen of the vector indicates the optimal congruence between
#' the \code{nth} matrix in lista and the listb (maximum possible value)
#' @examples
#' # calculate the simulated (original) data and the (re-covered) data
#' (not run)  sim <- CSSCASimulation(n_cluster, mem_cluster, n_block, n_com, n_distinct, n_var, p_sparse,
#'  p_noise, p_combase, p_fixzero, "both", mean_v)
#' (not run)  results <- MainCSSCA(sim$concatnated_data, n_var, n_block, n_com, n_distinct, n_cluster, n_obs, p_sparse, sim$cluster_assign, n_replace = 5)
#' # calculate the average congruence between the simulated loading matrices and estimated loading matrices
#' ListCongruence(sim$loading, results$loading)
#'

ListCongruence <- function(lista, listb, cluster_num){

  cong <- vector("numeric", length = cluster_num)
  for (i in 1:cluster_num){
    if(is.null(lista[[i]])){
      next
    }
    p <- 0
    cong_i <- vector("numeric", length = (cluster_num - 1))
    current <- lista[[i]]
    for (j in 1:cluster_num){
      if(is.null(listb[[j]])){
        next
      }
      p <- p + 1
      cong_i[p] <- TuckerCongruence(current, listb[[j]])
    }
    cong[i] <- max(cong_i)
  }
  return (cong)
}
