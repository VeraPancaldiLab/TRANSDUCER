library(RColorBrewer)
#' Plot PCA components in relation to a given factor
#' @description
#' this function does a parplot of the desired n of components and color
#' them according to a given factor.
#'
#' prints a list of spearman correlations of the desired metadata factor
#' with the IC that best correlate with it, from the ICAs with the different number of components
#'
#' @param pca_toplot data frame containing the PCs and the factors
#' you want to use to correlate
#' @param feat name of the feature to color
#' @param ncomp number of PCs to plot
#' @param dotsize to adjust the dotsize manually
#'
#' TBD
#' continuous coloring
#'
plot_PCs <- function(pca_toplot, feat, ncomp, dotsize) {
  col_factor <- as.factor(pca_toplot[[feat]])
  col_n <- nlevels(col_factor)
  cols <- brewer.pal(col_n, "Spectral")
  cols <- colorRampPalette(cols)(col_n)
  pairs(pca_toplot[, paste("Dim", 1:ncomp, sep = ".")],
        col = cols[as.numeric(col_factor)],
        pch = 19,
        cex = dotsize,
        lower.panel = NULL,
        main = feat
  )
  par(xpd = TRUE)
  # x11()
  plot.new()
  legend(x = "center", fill = cols, legend = levels(col_factor), horiz = F, title = feat)
}