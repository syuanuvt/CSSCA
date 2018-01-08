## the code for integrated SSCA 
# Yuan Shuai 2017/5/5

# version control
Sca_common <- function(data_con, nvar, nblock, common){
  
  # known_sparse: whether the sparsity of the data is known (0 = unknown, otherwise the input is the sparsity of the data)
  # data_con: concatenated data
  # nvar: the vector indivates the number of variables in each block
  # common: number of common components
  # initial: the rational start of the loading matrix
  # sparsity: if known, the sparsity of the component loadings
  # min_lasso: minimum value of penalty paramter (only for unknown sparsity)
  # max_lasso: maximum value of penalty paramter (only for unknown sparsity)
  
  upper <- 1e9
  
  if (length(nvar) == 1){
    nvar <- rep(nvar, nblock)
  }
  
  all_var <- ncol(data_con)
  all_mem <- nrow(data_con)
  
  # centering
  data_con <-  MatrixCenter(data_con, 1, 0)
  
  # find the global minimum (multiple start)
  result_svd <- svds(data_con, common)
  score_result <- result_svd$u
  singular_result <- result_svd$d
  loading_result1 <- result_svd$v
  loading_result <- t(t(loading_result1) * singular_result)
  
  residual <- data_con - score_result %*% t(loading_result)
  sum_residual <- sum(residual ^ 2)
  
  result <- list(t = score_result, p = loading_result, l = sum_residual) 
  return(result)
}
  