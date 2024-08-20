# List the n most absolute correlated genes with ICA3
#' @description
#' this function calculate the chosen absolute correlation between genes and ICA3 component
#' with either full samples (best overall) or ISR samples (best subset) and then subset the top n
#' @param normdf data frame of normalized RNA-seq count data. Lines are samples and columns are genes
#' @param select_method "bestsubset" to select ISR samples before doing the correlation or
#' "bestoverall" to keep all the samples
#' @param n number of gene to keep
#' @param cormethod "pearson" or "spearman"
#' @param top_samples data frame of top sample info
#' @param sample_info data frame of every sample info

mostCorrGenes <- function(normdf, select_method, n, cormethod, top_samples, sample_info) {
  # Keep only top samples if selection method is bestsubset
  if (select_method == "bestsubset") {
    normdf <- dplyr::filter(normdf, sample %in% top_samples$sample) %>%
      arrange(sample)
    sample_df <- top_samples
  } else if (select_method == "bestoverall") {
    normdf <- arrange(normdf, sample)
    sample_df <- sample_info
  }
  
  # Calculate absolute correlation for each gene
  genes_ic3_cor <- normdf %>%
    dplyr::filter(!sample == "!") %>%
    column_to_rownames("sample") %>% # select continuous variables
    as.matrix() %>%
    cor(y = sample_df$ICA3, method = cormethod) %>%
    as.data.frame()
  
  # keep the 1000 most absolute correlated genes
  top_genes <- mutate(genes_ic3_cor, V1 = abs(V1)) %>%
    arrange(desc(V1)) %>%
    dplyr::slice(1:n) %>%
    rownames_to_column(var = "EnsemblID")
  
  
  return(top_genes)
}