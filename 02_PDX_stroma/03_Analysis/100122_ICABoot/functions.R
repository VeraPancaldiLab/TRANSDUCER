library(JADE)
library(tidyverse)


### Performn JADE ICA to get a list
###  of S matrixes of a range of components
jade_range <- function(df, range.comp, MARGIN) {
  mats <- list()
  
  for (n.comp in range.comp) {
    jade_result <- JADE(df,
                        n.comp = n.comp
    )
    l_name <- paste(c("nc", n.comp), collapse = "")
    if (MARGIN == 1) {
      mats[[l_name]] <- jade_result[["A"]]
      suffix <- 1:ncol(jade_result[["A"]])
      colnames(mats[[l_name]]) <- paste("nc", suffix, sep = "")
    }
    
    if (MARGIN == 2) {
      mats[[l_name]] <- jade_result[["S"]]
    }
  }
  
  return(mats)
}

### Performn a bootstrap of a given df
### in rows (MARGIN =1) or columns
### (MARGIN = 2) keeping a speccified
### percentage of the data (perc)
bootstrap_df <- function(df,
                         MARGIN,
                         perc = 0.9) {
  if (MARGIN == 1) {
    # rows
    n <- round(nrow(df) * (1 - perc))
    
    out <- sample(rownames(df), n)
    df_boot <- df[!(rownames(df) %in% out), ]
    boot <- sample(rownames(df_boot), n)
    df_boot[out, ] <- df_boot[boot, ]
  }
  
  if (MARGIN == 2) {
    # cols
    n <- round(ncol(df) * (1 - perc))
    
    out <- sample(colnames(df), n)
    df_boot <- df[, !(colnames(df) %in% out)]
    boot <- sample(colnames(df_boot), n)
    df_boot[, out] <- df_boot[, boot]
  }
  return(df_boot)
}

### Wrap up function of bootstrap analysis.
### Perform a bootstrap of the given df in
### samples or genes and compute the
### median correlation of the components
### most correlating with each other
jade_choosencom <- function(df,
                            base_res = base_res,
                            range.comp = 2:12,
                            MARGIN = 1,
                            iterations = 1,
                            seed = 0) {
  set.seed(seed)
  df <- bootstrap_df(df, MARGIN = MARGIN)
  listof_correlations <- list(list())
  listof_correlations_id <- list()
  
  for (i in 1:iterations) {
    # loop through iterations
    res <- jade_range(df, range.comp, MARGIN)
    correlations_id <- list()
    
    for (nc in names(base_res)) {
      # loop through ICA runs
      base_ic <- base_res[[nc]]
      res_ic <- res[[nc]]
      cor_list <- as.list(rep(0, ncol(base_ic)))
      cor_id_list <- as.list(rep(0, ncol(base_ic)))
      names(cor_list) <- colnames(base_ic)
      names(cor_id_list) <- colnames(base_ic)
      
      for (b_ic in colnames(base_ic)) {
        # loop through boot ICA components
        for (r_ic in colnames(res_ic)) {
          # loop through base ICA components
          spearm_c <- cor(base_ic[, b_ic], res_ic[, r_ic])
          
          if (abs(spearm_c) > abs(cor_list[[b_ic]])) {
            # save the highest correlation between each
            # base ICA component and every bootstrap ICA
            cor_list[[b_ic]] <- spearm_c
            cor_id_list[[b_ic]] <- r_ic
          }
        }
      }
      # stopifnot(length(unique(unlist(cor_id_list))) == ncol(base_ic))
      listof_correlations[[nc]][[i]] <- cor_list
      # correlations_id[[nc]] <- cor_id_list
    }
    # listof_correlations_id[[i]] <- correlations_id
  }
  return(listof_correlations)
}

### Generate metrics for the representation
### of the bootstrap results
get_metrics <- function(bootstrap_results) {
  metrics <- data.frame()
  
  for (nc in names(bootstrap_results[-1])) {
    cor_values <- abs(unlist(bootstrap_results[[nc]]))
    metrics[nc, "components"] <- nc
    metrics[nc, "mean"] <- mean(cor_values)
    metrics[nc, "sd"] <- sd(cor_values)
    metrics[nc, "median"] <- median(cor_values)
  }
  return(metrics)
}