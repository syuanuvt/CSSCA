########################## the special initiation
#record the results
final_results <- list()
used_data <- list()
congruence <- list()
setting <- list()
global_min <- list()
ari_all <- list()
ari_icluster <- list()
ari_icluster_related <- list()
ari_our <- list()
factor_cong <- list()
score_cong <- list()
seeds<- list()
n_cluster_summary <- list()
n_component_summary <- list()
congruence[[1]] <- vector("numeric")
congruence[[2]] <- vector("numeric")


#########################################################

## the common settings
#n_cluster <- c(3,6,9)
n_block <- 2
n_var <- c(20, 20)

max_component <- 12
max_cluster <- 12
#n_res <- c(120, 240, 360)
condition <- 0
# THE PARAMETES to determine the number of iterations
iteration <- 1
vip_iteration <- 30
# iv a1: the noise level # set to 0.1
p_noise <- c(0.1, 0.3)
# iv a2: the sparsity level
p_sparse <- 0
# the proportion of mean structure variance
meanv <- c(0.1, 0.5, 0.9)
## iv b: the cluster membership
#mem_cluster <- c(40, 40, 40, 40, 40, 40)
#clustermem <- c(rep(1,40), rep(2,40), rep(3,40), rep(4,40), rep(5, 40), rep(6, 40))
#a_p_noise <- p_noise[1]

#mem_cluster <- list(c(40, 40, 40), c(72,24,24))
#clustermem <- list(c(rep(1, 40), rep(2, 40), rep(3,40)), c(rep(1,72), rep(2,24), rep(3,24)))
## iv d: the congruence level (first p_combase, second p_fixzero)
#cong <- list(c(0, 1), c(0.85, 1))
p_combase <- 0
p_fixzero <- 1

## iv e: the components
#n_com <- 2
#n_com1 <- 4 # the number of overall components
# fix the number of distinct components
n_distinct <- 1

## record the number
upper <- 1e9

# first generate the data
condition <- 0
for (e in 1:2){
  n_cluster <- 5 * e
  mem_cluster <- rep(40,n_cluster)
  n_respondents <- sum(mem_cluster)
  clustermem <- c()
  for (i in 1:n_cluster){
    clustermem <- c(clustermem, rep(i, 40))
  }
  for (c in 1:2){
    n_com <- 3 * c
    n_com1 <- 2 + n_com

   for (b in 1:2){
       meanv <- meanv[b]
       condition <- condition + 1
       #final_results[[condition]] <- list()
       used_data[[condition]] <- list()
       seeds[[condition]] <- vector("numeric", length = 25)
       setting[[condition]] <- list(num.cluster = n_cluster, meanv = meanv, num.component = n_com1)
       n_cluster_summary[[condition]] <- list()
       n_component_summary[[condition]] <- list()
       #global_min[[condition]] <- vector("numeric", length = 25)
       #ari_all[[condition]] <- vector("numeric", length = 25)
       #ari_icluster[[condition]] <- vector("numeric", length = 25)
       #ari_our[[condition]] <- vector("numeric", length = 25)
       #ari_icluster_related[[condition]] <- vector("numeric", length = 25)      
       #factor_cong[[condition]] <- list()
       
       for (times in 1){
         #factor_cong[[condition]][[times]] <- vector("numeric", length = n_cluster)
         # set the random seed
         seeds[[condition]][times] <- .Random.seed[times]
         set.seed(seeds[[condition]][times])
         toy <- CSSCASimulation.full.mean(n_cluster, mem_cluster, n_com, n_distinct, n_block, n_var, p_sparse, 
                                          p_noise[1], p_combase, p_fixzero, "both", meanv)
         block_version_data <- toy[[1]]
         all_data <- toy[[2]]
         true_score <- toy[[4]]
         true_loading <- toy[[5]]
         
         # record the data
         used_data[[condition]][[times]] <- list(data = toy, num.cluster = n_cluster, meanv = meanv, num.component = n_com1)
         # record the congruence
         #sum <- 0
         #n <- 0
         #for (i in 1:(n_cluster - 1)){
           #for (j in (i + 1):n_cluster){
             #n <- n + n_com + n_distinct * n_block
             #sum <- sum + tr(abs(factor.congruence(true_loading[[i]], true_loading[[j]])))
           #}
         #}
         #congruence[[d]] <- c(congruence[[d]], (sum / n))
         
         #min_loss <- upper
         ## start with non-sparse solution
         # first select the number of components and clusters based on csca
         est_csca_cluster <- varsel_nonsparse_full(all_data, n_var, n_block, max_component, max_cluster, n_respondents, iteration)[[1]]
         est_csca_component <- varsel_nonsparse_full(all_data, n_var, n_block, max_component, max_cluster, n_respondents, iteration)[[2]]
         # then estimate the results based on the selected variables
         #results_csca <- csca(all_data, n_var, n_block, est_csca_component, est_csca_cluster, n_respondents, vip.iteration)
         n_cluster_summary <- c(n_cluster_summary, est_csca_cluster)
         n_component_summary <- c(n_component_summary, est_csca_component)
         # then select the number of componets and clusters based on other mean-level clustering methods
         #est_mean_cluster <- iCluster2(block_version_data, n_cluster)$clusters
         # based on the number of cluster, determine the number of components
         
         
      }
    
    }
  }
}












        ## start with non-sparse solution
        varsel_nonsparse <- function(xc, nvar, nblock, max_components, max_cluster, nrespondents, iteration)
        
        results2 <- cssca_feed(all_data, n_var, n_block, n_com, n_distinct, n_cluster, n_res, p_sparse, results1[[1]])
        if (results2[[2]] < min_loss) {
          min_loss <- results2[[2]]
          global1 <- results2
        }
        for (k in 1:n_partition){
          new_cluster1 <- exchangemem(results2[[1]], n_cluster, (5*k))
          results3 <- cssca_feed(all_data, n_var, n_block, n_com, n_distinct, n_cluster, n_res, p_sparse, new_cluster1)
          if (results3[[2]] < min_loss) {
            min_loss <- results3[[2]]
            global1 <- results3
          }
        }
        
        ari_our[[condition]][times] <-  adjustedRandIndex(global1[[1]], clustermem)
        
        ## icluster solution
        min_loss <- upper
        results4 <- iCluster2(block_version_data, n_cluster)$clusters
        results5 <- cssca_feed(all_data, n_var, n_block, n_com, n_distinct, n_cluster, n_res, p_sparse, results4)
        if (results5[[2]] < min_loss) {
          min_loss <- results5[[2]]
          global2 <- results5
        }
        for (k in 1:n_partition){
          new_cluster2 <- exchangemem(results4, n_cluster, (5*k))
          results6 <- cssca_feed(all_data, n_var, n_block, n_com, n_distinct, n_cluster, n_res, p_sparse, new_cluster2)
          if (results6[[2]] < min_loss) {
            min_loss <- results6[[2]]
            global2 <- results6
          }
        }
        
        ari_icluster_related[[condition]][times] <- adjustedRandIndex(global2[[1]], clustermem)
        ari_icluster[[condition]][times] <- adjustedRandIndex(results4, clustermem)
        
        if (global1[[2]] < global2[[2]]){
          global <- global1
        }
        else {
          global <- global2
        }
        
        ari_all[[condition]][times] <- adjustedRandIndex(global[[1]], clustermem)
  
  
  
  
}

# first generate the data
condition <- 1
seeds[[condition]] <- .Random.seed[1]
set.seed(seeds[[condition]])
toy <- CSSCASimulation.full.mean(n_cluster, mem_cluster, n_com, n_distinct, n_block, n_var, p_sparse, 
                                 p_noise, p_combase, p_fixzero, "both", meanv)
block_version_data <- toy[[1]]
all_data <- toy[[2]]
true_score <- toy[[4]]
true_loading <- toy[[5]]

# test whether the data generation is appropriate
# sum(true_loading[[1]][1,] ^ 2) == 0.81

# try to solve the clustering problem with the non-sparse version of the algorithm
# try two versions: results1 refers to the version with total number of components; results2 refers to the version with 
# (total number of components - 1)
results1 <- csca(all_data, n_var, n_block, n_com1, n_cluster, n_res, 20)
min(results1[[2]])

results2 <- csca(all_data, n_var, n_block, (n_com1 - 1), n_cluster, n_res, 20)



#####################################

condition <- 0

for (d in 1:2){
  d_p_combase <- cong[[d]][1]
  d_p_fixzero <- cong[[d]][2]
  
        
          for (a in 1:2){
    a_p_noise <- p_noise[a]
    for (b in 1:3){
      b_meanv <- meanv[[b]]
      condition <- condition + 1
      final_results[[condition]] <- list()
      used_data[[condition]] <- list()
      seeds[[condition]] <- vector("numeric", length = 25)
      setting[[condition]] <- list(noise = a, meanv = b, cong = d)
      global_min[[condition]] <- vector("numeric", length = 25)
      ari_all[[condition]] <- vector("numeric", length = 25)
      ari_icluster[[condition]] <- vector("numeric", length = 25)
      ari_our[[condition]] <- vector("numeric", length = 25)
      ari_icluster_related[[condition]] <- vector("numeric", length = 25)      
      factor_cong[[condition]] <- list()
      
      # times in each condition
      for (times in 22:24){
        factor_cong[[condition]][[times]] <- vector("numeric", length = n_cluster)
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
        used_data[[condition]][[times]] <- list(data = toy, cong = d, noise = a, mean = c)
        sum <- 0
        n <- 0
        for (i in 1:(n_cluster - 1)){
          for (j in (i + 1):n_cluster){
            n <- n + n_com + n_distinct * n_block
            sum <- sum + tr(abs(factor.congruence(true_loading[[i]], true_loading[[j]])))
          }
        }
        congruence[[d]] <- c(congruence[[d]], (sum / n))
        
        min_loss <- upper
        
        ## start with non-sparse solution
        results1 <- csca(all_data, n_var, n_block, n_com1, n_cluster, n_res, 15)

          results2 <- cssca_feed(all_data, n_var, n_block, n_com, n_distinct, n_cluster, n_res, p_sparse, results1[[1]])
          if (results2[[2]] < min_loss) {
            min_loss <- results2[[2]]
            global1 <- results2
          }
          for (k in 1:n_partition){
            new_cluster1 <- exchangemem(results2[[1]], n_cluster, (5*k))
          results3 <- cssca_feed(all_data, n_var, n_block, n_com, n_distinct, n_cluster, n_res, p_sparse, new_cluster1)
          if (results3[[2]] < min_loss) {
            min_loss <- results3[[2]]
            global1 <- results3
          }
        }
        
        ari_our[[condition]][times] <-  adjustedRandIndex(global1[[1]], clustermem)
        
        ## icluster solution
        min_loss <- upper
        results4 <- iCluster2(block_version_data, n_cluster)$clusters
        results5 <- cssca_feed(all_data, n_var, n_block, n_com, n_distinct, n_cluster, n_res, p_sparse, results4)
        if (results5[[2]] < min_loss) {
          min_loss <- results5[[2]]
          global2 <- results5
        }
        for (k in 1:n_partition){
          new_cluster2 <- exchangemem(results4, n_cluster, (5*k))
          results6 <- cssca_feed(all_data, n_var, n_block, n_com, n_distinct, n_cluster, n_res, p_sparse, new_cluster2)
          if (results6[[2]] < min_loss) {
            min_loss <- results6[[2]]
            global2 <- results6
          }
        }
        
        ari_icluster_related[[condition]][times] <- adjustedRandIndex(global2[[1]], clustermem)
        ari_icluster[[condition]][times] <- adjustedRandIndex(results4, clustermem)
        
        if (global1[[2]] < global2[[2]]){
          global <- global1
        }
        else {
          global <- global2
        }
        
        ari_all[[condition]][times] <- adjustedRandIndex(global[[1]], clustermem)
        
        results7 <- cssca_feed(all_data, n_var, n_block, n_com, n_distinct, n_cluster, n_res, p_sparse, clustermem)
        if (global[[2]] < results7[[2]] | global[[2]] == results7[[2]] | adjustedRandIndex(global[[1]], results7[[1]]) == 1){
          global_min[[condition]][times] <- 1
        }
        else{
          global_min[[condition]][times] <- 0
        }
        
        ## compute the congruence of loadings
        factor_cong[[condition]][[times]] <- paircongruence(global[[3]], true_loading, n_cluster)
        
        final_results[[condition]][[times]] <- list(all = global[[1]], icluster = results4, true = results7[[1]], full = global)
        
      }
      
    }
  }
}
#save(ari_all, ari_icluster, ari_icluster_related, ari_our, final_results, used_data, congruence, setting, global_min, factor_cong, file = "simulation.RData")




