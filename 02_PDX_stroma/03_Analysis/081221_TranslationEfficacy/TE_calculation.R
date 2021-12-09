library(tidyverse)
library(Hmisc)
library(corrplot)
library(lmtest)
library(ggridges)
library(Biobase)
library(factoextra)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/081221_TranslationEfficacy/")

# Data loading
normHost_cyt <- read_tsv("../../00_Data/Processed_data/normHost_Cyt.tsv")
normHost_pol <- read_tsv("../../00_Data/Processed_data/normHost_Pol.tsv")
sample_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>% column_to_rownames("sample")

# Joint dataframe creation
sn_cyt <- paste(colnames(normHost_cyt[-1]), "cyt", sep = "_")
sn_pol <- paste(colnames(normHost_pol[-1]), "pol", sep = "_")

normHost_cyt.tmp <- normHost_cyt
normHost_pol.tmp <- normHost_pol

colnames(normHost_cyt.tmp) <- c("EnsemblID",sn_cyt)
colnames(normHost_pol.tmp) <- c("EnsemblID",sn_pol)

norm_host <- inner_join(normHost_cyt.tmp, normHost_pol.tmp, by= "EnsemblID")

# eset creation
## PhenoData
sn <- c(colnames(normHost_cyt[-1]), colnames(normHost_pol[-1]))

frac <- c(rep("cyt", length(sn_cyt)),
          rep("pol", length(sn_pol)))

rna_conc <- c(sample_info[colnames(normHost_cyt[-1]), "RNAconc_cyt"],
              sample_info[colnames(normHost_pol[-1]), "RNAconc_pol"])

diab_stat <- c(sample_info[colnames(normHost_cyt[-1]), "Diabetes"],
               sample_info[colnames(normHost_pol[-1]), "Diabetes"])
# diab_stat %>% replace_na("Uknown") -> diab_stat

tumor_ica <- bind_rows(sample_info[colnames(normHost_cyt[-1]),str_subset(colnames(sample_info), "ICA")],
                       sample_info[colnames(normHost_pol[-1]),str_subset(colnames(sample_info), "ICA")],)

rownames(tumor_ica) <- c(sn_cyt, sn_pol)
  
phenotype.tmp <- data.frame(sample = factor(sn),
                        fraction = factor(frac),
                        rna_conc = as.double(rna_conc),
                        diabetes = factor(diab_stat))

rownames(phenotype.tmp) <- colnames(norm_host)[-1]

phenotype <- merge(phenotype.tmp, tumor_ica, by = "row.names",) %>% column_to_rownames("Row.names")
phenotype <- phenotype[colnames(norm_host)[-1],]


# PCA
norm_host %>%
  column_to_rownames("EnsemblID") %>% data.matrix() %>%
  t() %>% prcomp(scale=T) -> norm_pca

# Var explained
var_explained <- norm_pca$sdev^2/sum(norm_pca$sdev^2)
fviz_eig(norm_pca, barfill = "lightgrey",
         barcolor = "black", title = "PDX cyt/pol PCA",
         subtitle = paste("90% variance reached at the",
                          which(cumsum(var_explained)>0.9)[1],
                          "th component"))


all(rownames(norm_pca$x) == rownames(phenotype)) %>% stopifnot()
pca_toplot <- cbind(norm_pca$x, phenotype, by="row.names")


#' Plot PCA components in relation to a given factor
#'@description
#' this function does a parplot of the desired n of components and color 
#' them according to a given factor.
#'
#'prints a list of spearman correlations of the desired metadata factor
#'with the IC that best correlate with it, from the ICAs with the different number of components 
#'
#'@param pca_toplot data frame containing the PCs and the factors 
#'you want to use to correlate
#'@param feat name of the feature to color
#'@param ncomp number of PCs to plot
#'@param dotsize to adjust the dotsize manually
#'
#'TBD 
#'continuous coloring
#'

plot_PCs <- function(pca_toplot, feat, ncomp, dotsize){
  col_factor <- as.factor(pca_toplot[[feat]])
  col_n <- nlevels(col_factor)
  cols <- brewer.pal(col_n, "Spectral")
  cols <- colorRampPalette(cols)(col_n)
  pairs(pca_toplot[, paste("PC",1:ncomp,sep ="")],
        col = cols[as.numeric(col_factor)],
        pch = 19,
        cex = dotsize,
        lower.panel = NULL)
  
  par(xpd = TRUE)
  x11()  
  plot.new()
  legend(x = "center",fill = cols, legend = levels(col_factor), horiz = F)
}

plot_PCs(pca_toplot, "fraction", 10, 0.5)
plot_PCs(pca_toplot, "diabetes", 10, 0.5)
plot_PCs(pca_toplot, "sample", 10, 0.5)
plot_PCs(pca_toplot, "rna_conc", 10, 0.5)
plot_PCs(pca_toplot, "ICA1", 10, 0.5)
plot_PCs(pca_toplot, "ICA2", 10, 0.5)
plot_PCs(pca_toplot, "ICA3", 10, 0.5)
plot_PCs(pca_toplot, "ICA4", 10, 0.5)
plot_PCs(pca_toplot, "ICA5", 10, 0.5)
plot_PCs(pca_toplot, "ICA6", 10, 0.5)

# Similarity through whole genome correlation
statistic = "pearson"
norm_host %>% column_to_rownames("EnsemblID") %>%
  data.matrix() %>% rcorr(.,type = statistic) -> corr
corrplot(corr = corr$r,
         p.mat = corr$P,
         is.corr = F, order = "hclust",
         type = "lower",
         main = paste(statistic, sep = ": "))

# Translation efficacy analysis