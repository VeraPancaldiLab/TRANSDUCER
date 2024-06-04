# Plots a scatterplot indicating also Pearson and spearman correlation coefficients of two given columns of a tibble

correlation_plotter <- function(data, col1, col2, data_name){
  corr_spearman <- rcorr(data[[col1]], data[[col2]], type = "spearman")
  corr_pearson <- rcorr(data[[col1]], data[[col2]], type = "pearson")
  stats <- paste0("Spearman: R = ", round(corr_spearman$r["x","y"], 2), ", pval = ", round(corr_spearman$P["x","y"], 4),
                 "\nPearson: R = ", round(corr_pearson$r["x","y"], 2), ", pval = ", round(corr_pearson$P["x","y"], 4))
  
  ggplot(data) +
    aes_string(col1, col2) +
    geom_point(shape = 16, size = 2, show.legend = FALSE) +
    geom_smooth(method=lm) +
    theme_minimal() +
    labs(title = paste0("Comparison between ", col1, " and ", col2, " in ", data_name),
         subtitle = stats)
}
