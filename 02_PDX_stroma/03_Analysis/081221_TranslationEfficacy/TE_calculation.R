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


# function to parplot the desired ammount of PCs related to a given factor.
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
## samples (only in first two)
pca_toplot %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2, color = sample)) + geom_point(size=4) +
  theme_bw(base_size=15) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")



phenoData <- new("AnnotatedDataFrame", data = phenotype)
norm_host %>% column_to_rownames("EnsemblID") %>% data.matrix() -> assayData

dset <- ExpressionSet(assayData = assayData, phenoData = phenoData)
