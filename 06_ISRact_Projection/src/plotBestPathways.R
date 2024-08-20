## Visualisation
### fGSEA plots
# Function to plot the top 10 low ES and high ES pathways regarding to the pvalue of the gsea
#' @description
#' this function plots table of enrichment graphs for the 10
#' pathways with positive enrichment score and lowest pvalue
#' and the 10 pathways with negative enrichment score and lowest pvalue
#'
#' @param gsea.res data frame of a gsea result returned by fgsea function
#' @param msigdbr_list list of genes name for each pathway from msigdb
#' @param ranked_genes named vector where names are genes name and values are pca values
plotBestPathways <- function(gsea.res, msigdbr_list, ranked_genes) {
  topPathwaysUp <- gsea.res[ES > 0][head(order(pval), n = 10), pathway]
  topPathwaysDown <- gsea.res[ES < 0][head(order(pval), n = 10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plot.new()
  plotGseaTable(msigdbr_list[topPathways], ranked_genes, gsea.res,
                gseaParam = 0.5
  )
}
