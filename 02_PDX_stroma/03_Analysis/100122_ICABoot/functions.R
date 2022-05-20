library(JADE)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(scico)
library(ggVennDiagram)
library(RCy3)
#library(ggcorrplot)

### Performn JADE ICA to get a list
###  of S or A matrices of a range of components
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
    if (jade_result[[1]] %>% is.na() %>% all()) {
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
### in rows (MARGIN = 1) or columns
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
### samples or genes and return the correlation
### of each bootstrap ic with the most 
### alike original ic
jade_choosencom <- function(df,
                            base_res = base_res,
                            MARGIN = 1,
                            iterations = 1,
                            seed = 0,
                            perc = 0.9) {
  range.comp <- as.numeric(gsub("nc", "", names(base_res)))
  set.seed(seed)
  listof_correlations <- list()
  listof_correlations_id <- list()

  for (i in 1:iterations) {
    # loop through iterations
    df <- bootstrap_df(df, MARGIN = MARGIN, perc = perc)
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
          pears_c <- cor(base_ic[, b_ic], res_ic[, r_ic])

          if (abs(pears_c) > abs(cor_list[[b_ic]])) {
            # save the highest correlation between each
            # base ICA component and every bootstrap ICA
            cor_list[[b_ic]] <- pears_c
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

### Generate metrics for the line
### representation of the bootstrap results
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

### Function to produce complete boxplot
### of the correlations distribution outputted
### by bootstrapping both samples and genes, and 
### a simple lineplot of the choosen stat ("mean", "median")
boot_plots <- function(s_boot, g_boot, line_stat = "mean", name = "analysis"){
  plot_title <- paste("02_Output/", name, "_bootstrapICA", sep = "")
  range.comp <- as.numeric(gsub("nc", "", names(s_boot)))
  
  #### Boxplot
  all_boot <- list(g_boot, s_boot)
  names(all_boot) <- c("genes", "samples")
  all_melt <- melt(all_boot)
  
  colnames(all_melt) <- c("correlation", "c", "iteration", "components", "bootstrap")
  all_melt$correlation <- abs(all_melt$corr)
  level_order <- factor(all_melt$components, level = names(s_boot))
  
  
  ggplot(all_melt, aes(x = level_order, y = correlation, fill = bootstrap)) +
    geom_boxplot(width = 0.5) +
    #scale_x_discrete(limits = paste("nc", range.comp, sep = "")) +
    labs(y = "absolute pearson correlation", x = "number of components") +
    #coord_cartesian(ylim = c(0.8, 1)) +
    theme_bw()
  ggsave(paste(plot_title, "boxplot.pdf", sep = "_"), height = 10, width = 12)
  
  
  
  #### line plot of a measure
  gene_metrics <- get_metrics(g_boot)
  sample_metrics <- get_metrics(s_boot)
  corrlim <- min(c(min(gene_metrics[, line_stat]),
                   min(sample_metrics[, line_stat])))
  
  pdf(paste(plot_title, "lineplot.pdf", sep = "_"))
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
  dev.off()
}


### Gets an A matrix result from JADE (annotated with rownames/colnames)
### and plots each component density with a rugplot. This function dos as 
### many plots as columns the annotation df has. carefull. Sample order 
### should be checked beforehand.

plot_sample_weights <- function(A_mat, annotations, cont_names, analysis_name){
  pdf(file=paste("02_Output/", analysis_name, ".pdf", sep=""))
  
  # Correlation plot
  corrplot <- annotations %>%
    dplyr::select(all_of(cont_names)) %>%
    bind_cols(A_mat) %>%
    formatted_cors(cor.stat = "spearman") %>%
    filter(measure1 %in% names(annotations), measure2 %in% names(A_mat)) %>% #not square corr
    mutate(r = abs(r)) %>% #abscorr
    ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations continuous vars",
         subtitle="Only significant correlation coefficients shown (95% I.C.)") +
    scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(0,1)) +
    geom_text() +
    theme_classic() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggpubr::rotate_x_text(angle = 90)
  
  print(corrplot)
  
  # Sample weights
  for (ann in colnames(annotations)){

    rug_aes <- annotations[[ann]]
    rug_name <- ann
    comps_plots <- lapply(colnames(A_mat), function(ic){
      p <- 
        ggplot(A_mat) +
        aes_string(ic) +
        geom_density() + 
        geom_rug(aes(color = rug_aes), length = unit(0.1, "npc")) +
        labs(color = ann) +
        theme_classic() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank())
      
      if(is.integer(rug_aes)) {
        p <- p  +
          scale_color_discrete()
        
      } else if(is.double(rug_aes)){
        p <- p  +
          scico::scale_color_scico(palette = "berlin")
      } 
      
      # if (!(ic %in% c("IC.1", paste("IC", 1+trunc(elected_ncomp/2), sep = ".")))) { # this is to control wich cells should have a ylabel
      #   p <- p +
      #     theme(axis.title.y = element_blank())
      # }
      p
      
    })
    ggarrange(plotlist = comps_plots, common.legend = T,  legend = "bottom") %>% annotate_figure(
      top = text_grob(ann, color = "black", face = "bold", size = 14), ) %>% print()
  }
  dev.off()
}

# Customizable corrplot (modified from https://www.khstats.com/blog/corr-plots/corr-plots/)
cors <- function(df, cor.stat) {
  M <- Hmisc::rcorr(as.matrix(df), type = cor.stat)
  Mdf <- map(M, ~data.frame(.x))
  return(Mdf)
}

formatted_cors <- function(df, cor.stat){
  cors(df, cor.stat) %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    dplyr::rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
}


dfs_corrplot <- function(df_x, df_y, abscorr = FALSE) {
  
  stopifnot(rownames(df_x)==rownames(df_y))
  joint_df <- df_y %>% bind_cols(df_x)
  
  df_corr <- rcorr(data.matrix(joint_df), type = "spearman")
  df_corr$r <- df_corr$r[colnames(df_x), colnames(df_y)]
  df_corr$P <- df_corr$P[colnames(df_x), colnames(df_y)]
  
  if (abscorr == TRUE){
    cmap = c("white", "white","red")
    scale_lims = c(0,1)
  } else {
    cmap = c("#7D0E29", "white", "#004376")
    scale_lims = c(-1,1)
  }
  
  ggcorrplot(df_corr$r, p.mat = df_corr$P) +
    scale_fill_gradient2(limit = scale_lims, low = cmap[1], mid = cmap[2], high =  cmap[3])
  
}


Plot_deconv <- function(deconv, complete_annotation, analysis_name){
  # Formating possible blanks in colnames
  deconv <- deconv %>%
    rename_with(~str_replace_all(.," ", "."))
  
  # Heatmap
  deconv %>% t() %>%
    pheatmap(scale = "row",
             annotation_col = complete_annotation) %>% as.grob() -> decon_heatmap
  
  # Corplot
  stopifnot(rownames(complete_annotation)==rownames(deconv))
  decon_corrplot <- deconv %>%
    bind_cols(complete_annotation) %>% 
    formatted_cors(cor.stat = "spearman") %>%
    filter(measure1 %in% names(complete_annotation), measure2 %in% names(deconv)) %>% #not square corr
    #mutate(r = abs(r)) %>% #abscorr
    ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation", title="",
         subtitle="Only significant correlation coefficients shown (95% I.C.)") +
    scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
    geom_text() +
    theme_classic() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggpubr::rotate_x_text(angle = 90)
  
  
  fig <- ggarrange(decon_corrplot + rremove("xylab"), decon_heatmap, widths = c(0.6,1),
                   labels = c("A", "B"),
                   ncol = 2, nrow = 1)
  
  pdf(paste("02_Output/", analysis_name, ".pdf", sep=""), nrow(deconv), ncol(deconv))
  print(annotate_figure(fig, top = text_grob(paste(analysis_name, "deconvolution results"), face = "bold", size = 20)))
  dev.off()
  
}

#' Filter a dataframe to keep just half the most variable genes
#'@description
#' `Get_half_mostvar` requires a data.frame with genes as rows,
#'  and its ids as rownames. returns a dataframe

Get_mostvar <- function(df, n){
  
  abdv <- apply(df, 1, function(x) {
    sum(
      abs(
        x - mean(x)
      )
    ) / length(x)
  })
  
  # to input either % or absolute n
  if ( n>0 && n<1){
    th_abdv.i <- nrow(df)*n
  } else{
    th_abdv.i <- n
  }
  
  abdv %>% sort(decreasing = T) %>% .[th_abdv.i] -> th_abdv
  df.f <- df[abdv >= th_abdv, , drop=F]
  return(df.f)
}


Plot_general_TFs <- function(tf_activity, analysis_name, n_mostvar, complete_annotation){
  
  tf_activity %>%
    pheatmap(main=paste( "TF activity", analysis_name),
             scale = "row", show_rownames = FALSE,
             annotation_col = complete_annotation) -> full_heatmap
  
  mostvar_TF <- Get_mostvar(tf_activity, n_mostvar)
  
  mostvar_TF %>%
    pheatmap(scale = "row", annotation_col = complete_annotation) %>% as.grob() -> mostvar_heatmap
  
  # Corplot
  stopifnot(rownames(complete_annotation)==colnames(mostvar_TF))
  mostvar_corrplot <- mostvar_TF %>%
    t() %>%
    as_tibble() %>%
    bind_cols(complete_annotation) %>%
    formatted_cors(cor.stat = "spearman") %>%
    filter(measure1 %in% names(complete_annotation), measure2 %in% rownames(mostvar_TF)) %>% #not square corr
    #mutate(r = abs(r)) %>% #abscorr
    ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation", title="",
         subtitle="Only significant correlation coefficients shown (95% I.C.)") +
    scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
    geom_text() +
    theme_classic() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggpubr::rotate_x_text(angle = 90)
  
  
  fig <- ggarrange(mostvar_corrplot, mostvar_heatmap, widths = c(0.6,1),
                   labels = c("A", "B"),
                   ncol = 2, nrow = 1)
  
  pdf(paste("02_Output/TF_analysis_", analysis_name, ".pdf", sep=""), nrow(mostvar_TF)/2, ncol(mostvar_TF)/2)
  print(full_heatmap)
  print(annotate_figure(fig, top = text_grob(paste(analysis_name, "vs the most variable TF activities"), face = "bold", size = 20)))
  dev.off()
}


#' Select the nTFs most absolutely correlated TFs with each factor 
#' and do a heatmap of each of them.
#'@description
#'
PlotBestCorr <- function(complete_annotation, tf_activity, nTFs, analysis_name = "analysis"){
  for (comp in colnames(complete_annotation)){
    
    rcorr(t(tf_activity), complete_annotation[,comp,drop=T])$r[,"y"] %>%
      abs() %>% sort(decreasing = T) %>%
      .[2:(nTFs+1)] %>% names() -> best_tfs
    
    pdf(paste("02_Output/", analysis_name, comp, ".pdf", sep=""))
    tf_activity %>% .[best_tfs,] %>%
      pheatmap(main=paste("Best",analysis_name, comp), scale = "row", annotation_col = complete_annotation)
    
    dev.off()
  }
}


#' Select the nTFs most absolutely correlated TFs with each factor 
#' and do a heatmap of each of them.
#'@description
#'
PlotGeneWeights <- function(S_mat, ensembl_toplot, n_genes, translate, complete_annotation, analysis_name = "analysis"){
  for (comp in colnames(S_mat)){
    S_mat %>% arrange(get(comp)) -> S_sort
    S_sort %>% head(n_genes) -> S_mneg
    S_sort %>% tail(n_genes) -> S_mpos
    S_mat %>% ggplot() +
      aes_string(y = comp) +
      geom_density(alpha=.5, fill="#AED3FA") +
      geom_hline(yintercept = 0, colour = "black") +
      geom_hline(yintercept = max(S_mneg[comp]), colour = "#7DB0DD", linetype = "dashed") + 
      geom_hline(yintercept = min(S_mpos[comp]), colour = "#EAAFBB", linetype = "dashed") + 
      theme_classic() -> densplot
    
    annot_row_ <- tibble(name = rownames(S_mpos), class = "most possitive")
    tibble(name = rownames(S_mneg), class = "most negative") %>% bind_rows(annot_row_) %>% column_to_rownames("name") -> annot_row__
    annot_row <- annot_row__
    rownames(annot_row__) %>% translate[.] %>% make.names(unique = TRUE)  -> rownames(annot_row)
    
    ensembl_toplot %>% dplyr::filter(EnsemblID %in% c(rownames(S_mpos), rownames(S_mneg))) %>%
      arrange(match(EnsemblID, c(rownames(S_mpos), rownames(S_mneg)))) -> ensembl_toplot_
    
    genes_toplot <- ensembl_toplot_
    ensembl_toplot_$EnsemblID %>% translate[.] %>%
      make.names(unique = TRUE) -> genes_toplot$Genenames
    
    genes_toplot %>% dplyr::select(!EnsemblID) %>%
      relocate(Genenames) %>% column_to_rownames("Genenames") %>% 
      pheatmap(scale = "row", annotation_row = annot_row, cluster_rows = F, 
               annotation_col = complete_annotation[c("PAMG", "ISRact", comp)]) %>% as.grob() -> heatmap
    
    fig <- ggarrange(densplot + rremove("xylab"), heatmap, heights = c(1.5, 10), widths = c(0.2,1),
                     labels = c("A", "B"),
                     ncol = 2, nrow = 1)
    pdf(paste("02_Output/", analysis_name, comp, ".pdf", sep=""))
    print(annotate_figure(fig, top = text_grob(paste(comp, "Gene weights"), face = "bold", size = 20)))
    dev.off()
  }
}

# Generate a subnetwork of the s
# based in the X% of most extreme gene weights of the valid_comps.

PlotNetwork <- function(nodes, edges, S_mat, valid_comp, main_name){
  thresholds <- apply(S_mat, 2, function(x) quantile(abs(x), 0.95))
  bestgenes <- lapply(names(thresholds), function(x)
    rownames(dplyr::filter(S_mat, abs(get(x)) > thresholds[x])))
  names(bestgenes) <- names(thresholds)
  ggVennDiagram(bestgenes, label = "count", label_alpha = 0) +
    scale_fill_gradient(low="grey", high = "red",  trans = "log2")

  is_bestgene <- lapply(names(bestgenes), function(x) nodes$id %in% bestgenes[[x]]) %>% 
    as_tibble(.name_repair= "universal" ) %>% 
    rename_all(~paste0("best_IC.",valid_comp)) %>%
    mutate(best_any = rowSums(.), 
           best_for = ifelse(best_any == 1, apply(., 1, function(x) colnames(.)[as.logical(x)][1]), ifelse(best_any == 0, 0, "many")))
  
  nodes_ <- bind_cols(nodes, is_bestgene)
  network <- createNetworkFromDataFrames(nodes_, edges, title = main_name)
  createColumnFilter(filter.name='best_IC_filter', column='best_any', 0, 'IS_NOT', network = network)
  createSubnetwork(subnetwork.name = paste0(main_name,'_5%'))
}


