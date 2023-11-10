library(treeio)
library(ape)
library(dplyr)
library(phytools)
library(tibble)
library(progress)
library(tidyr)

## Compute error metrics for all the simulated trees 
## and save to file
foobar <- function(model_index, tree_height, subdir){
  l1 <- list()
  l2 <- list()
  #tree_heights <- c(25.0, 50.0, 75.0, 100.0, 125.0)
  
  n_replicates <- 100
  pb <- progress_bar$new(total = n_replicates)
  for (i in 1:n_replicates){
    prefix <- "simulated_trees/"
    postfix <- paste0("s1_h", tree_height, ".0_m", model_index, "_", i, ".tre")
    file1 <- paste0(prefix, postfix)
    file2 <- paste0(prefix, "results/", subdir, "/", postfix)
    if (!file.exists(file1) | !file.exists(file2)){
      l[[i]] <- NULL
    }else{
      tr1 <- read.beast.newick(file1)
      tr2 <- read.beast.newick(file2)
      
      true_netdiv <- tr1@data$r
      inferred_netdiv <- tr2@data$mean_netdiv
      true_lambda <- tr1@data$lambda
      inferred_lambda <- tr2@data$mean_lambda
      true_mu <- tr1@data$mu
      inferred_mu <- tr2@data$mean_mu
      l2[[i]] <- inferred_netdiv
      n_branches <- nrow(tr2@data)
      r_mse <- mean((true_netdiv - inferred_netdiv)^2)
      
      true_N <- sapply(tr1@data$N, sum)
      inferred_N <- tr2@data$nshift
      sim_is_shift <- true_N > 0
      
      criteria <- c("N_over_half", "bayes_factor_10","N_half_and_bayes_factor")
      ls5 <- list()
      for (m in 1:3){
        if (m == 1){
          inference_is_shift <- tr2@data$nshift > 0.5  ## N > 0.5 criterion
        }else if(m == 2){
          inference_is_shift <- tr2@data$shift_bf > 10 ## bayes factor criterion  
        }else{
          inference_is_shift <- (tr2@data$shift_bf > 10) &  (tr2@data$nshift > 0.5) ## bayes factor criterion   and N_half
        }
        
        
        #height <- max(node.depth.edgelength(tr1@phylo))
        ntaxa <- Ntip(tr1)
        
        N_mse <- mean((true_N - inferred_N)^2)
        
        pos <- sum(sim_is_shift)
        neg <- sum(!sim_is_shift)
        
        tp <- sum(sim_is_shift & inference_is_shift)
        tn <- sum((!sim_is_shift) & (!inference_is_shift))
        fn <- sum(sim_is_shift & (!inference_is_shift))
        fp <- sum((!sim_is_shift) & inference_is_shift)
        
        n_edges <- nrow(tr1@data)
        
        log_prop_error_lambda <- log(inferred_lambda) - log(true_lambda)
        log_prop_error_mu <- log(inferred_mu) - log(true_mu)
        
        ## geometric mean
        prop_error_geomean_lambda <- exp(mean(log_prop_error_lambda))
        ## arithmetic mean
        prop_error_mean_lambda <- mean(exp(log_prop_error_lambda))
        
        ## geometric mean
        prop_error_geomean_mu <- exp(mean(log_prop_error_mu))
        ## arithmetic mean
        prop_error_mean_mu <- mean(exp(log_prop_error_mu))
        
        df5 <- tibble(
          "model" = model_index,
          "inference" = subdir,
          "replicate" = i,
          "height" = tree_height,
          "ntaxa" = ntaxa,
          "r_mse" = r_mse,
          "N_mse" = N_mse,
          "N_true_sum" = sum(true_N),
          "true positive" = tp,
          "true negative" = tn,
          "false negative" = fn,
          "false positive" = fp,
          "false positive rate" = fp / neg,
          "true positive rate" = tp / pos,
          "false negative rate" = fn / pos,
          "true negative rate" = tn / neg,
          "prop_error_geomean_lambda" = prop_error_geomean_lambda,
          "prop_error_mean_lambda" = prop_error_mean_lambda,
          "prop_error_geomean_mu" = prop_error_geomean_mu,
          "prop_error_mean_mu" = prop_error_mean_mu,
          "criterion" = criteria[m]
        )
        ls5[[m]] <- df5
      }
      df1 <- bind_rows(ls5)
      
      l1[[i]] <- df1

      pb$tick()
    }
  }
    
  res <- bind_rows(l1)
  return(res)
}

subdirs <- c("true_all", "true_div", "unknown_rates")
tree_heights <- c(25.0, 50.0, 75.0, 100.0, 125.0)

q <- 1.0
l <- list()
for (i in seq_along(subdirs)){
  subdir <- subdirs[i]
  for (j in seq_along(tree_heights)){
    tree_height <- tree_heights[j]
    for (model_index in 1:4){
      print(paste0("subdir = ", subdir, ", model = ", model_index, ", height = ", tree_height))
      l[[q]] <- foobar(model_index, tree_height, subdir)
      q <- q + 1
    }
  }
}
df1 <- bind_rows(l)
#write.csv(df1, file = "output/branch_specific_estimation_error.csv")
