library(JADE)
library(tidyverse)
library(reshape2)

### Performn JADE ICA to get a list
###  of S matrixes of a range of components
jade_range <- function(df, range.comp, MARGIN) {
  mats <- list()

  for (n.comp in range.comp) {
    l_name <- paste(c("nc", n.comp), collapse = "")
    jade_result <- tryCatch(
      {
        JADE(df, n.comp = n.comp)
      },
      error = function(cond) {
        message(paste("JADE with ", n.comp, "did not converge, skipping iter"))
        return(NA)
      }
    )
    if (jade_result[1] %>% is.na()) {
      print("Skipping iteration")
      next
    }

    if (MARGIN == 1) {
      mats[[l_name]] <- jade_result[["A"]]
      suffix <- 1:ncol(jade_result[["A"]])
      colnames(mats[[l_name]]) <- paste("IC", suffix, sep = ".")
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
                            MARGIN = 1,
                            iterations = 1,
                            seed = 0) {
  range.comp <- as.numeric(gsub("nc", "", names(base_res)))
  set.seed(seed)
  df <- bootstrap_df(df, MARGIN = MARGIN)
  listof_correlations <- list()
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

  for (nc in names(bootstrap_results)) {
    cor_values <- abs(unlist(bootstrap_results[[nc]]))
    metrics[nc, "components"] <- nc
    metrics[nc, "mean"] <- mean(cor_values)
    metrics[nc, "sd"] <- sd(cor_values)
    metrics[nc, "median"] <- median(cor_values)
  }
  return(metrics)
}

boot_plots <- function(s_boot, g_boot, line_stat = "mean"){
  
  range.comp <- as.numeric(gsub("nc", "", names(s_boot)))
  
  #### Boxplot
  all_boot <- list(g_boot, s_boot)
  names(all_boot) <- c("genes", "samples")
  all_melt <- melt(all_boot)
  
  colnames(all_melt) <- c("correlation", "c", "iteration", "components", "bootstrap")
  all_melt$correlation <- abs(all_melt$corr)
  level_order <- factor(all_melt$components, level = names(s_boot))
  
  print(
    ggplot(all_melt, aes(x = level_order, y = correlation, fill = bootstrap)) +
    geom_boxplot(width = 0.5) +
    #scale_x_discrete(limits = paste("nc", range.comp, sep = "")) +
    labs(y = "absolute pearson correlation", x = "number of components") +
    coord_cartesian(ylim = c(0.8, 1)) +
    theme_bw()
  )
  
  
  #### line plot of a measure
  gene_metrics <- get_metrics(g_boot)
  sample_metrics <- get_metrics(s_boot)
  corrlim <- min(c(min(gene_metrics[, line_stat]),
                  min(sample_metrics[, line_stat])))
  
  plot(gene_metrics[, line_stat], ylim = c(corrlim, 1),
       type = "b", lty = 1, pch = 19, col = "red",
       xaxt = "n", xlab = "n of components", ylab = "Absolute pearson correlation"
  )

  lines(sample_metrics[, line_stat], type = "b", lty = 2, pch = 8, col = "blue")

  legend("bottomleft",
         legend = c("probe", "sample"),
         col = c("red", "blue"), lty = 1:2, cex = 0.8
  )

  title(paste(line_stat, "distribution"))
  axis(side = 1, at = 1:nrow(gene_metrics), labels = gene_metrics[, "components"])

}
