library(tidyverse)
library(anota2seq)
library(anota2seqUtils)
library(biomaRt)
library(sva)
library(factoextra)
library(edgeR)
library(RColorBrewer)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/06_ISRact_Projection/")
################################################################################
# Plot ncomp features
plot_PCs <- function(pca_toplot, feat, ncomp, dotsize){
  col_factor <- as.factor(pca_toplot[[feat]])
  col_n <- nlevels(col_factor)
  cols <- brewer.pal(col_n, "Spectral")
  cols <- colorRampPalette(cols)(col_n)
  pairs(pca_toplot[, paste("PC",1:ncomp,sep ="")],
        col = cols[as.numeric(col_factor)],
        pch = 19,
        cex = dotsize,
        lower.panel = NULL, 
        main = feat)
  
  par(xpd = TRUE)
  x11()  
  plot.new()
  legend(x = "center",fill = cols, legend = levels(col_factor), horiz = F, title = feat)
}
################################################################################

# Data loading
## Sample info and selection of data
include_med = FALSE

sample_info <- read_tsv("data/Sauyeun_PDX/sample_info.tsv")

top_samples <- dplyr::arrange(sample_info, ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact = ifelse(ICA3 < 0, "low", "high")) %>%
  dplyr::filter(sample != "PDAC009T")

sample_info <- dplyr::select(top_samples, sample, ISRact) %>%
  right_join(sample_info, by= "sample") %>% 
  dplyr::mutate(ISRact = if_else(is.na(ISRact), "medium", ISRact))

if (include_med == FALSE){
  sample_info <- column_to_rownames(top_samples, "sample")
} else {
  sample_info <- column_to_rownames(sample_info, "sample")
}

# load and filter expression data
rowTumor_cyt <- read_tsv("data/Sauyeun_PDX/rawTumor_Cyt.tsv") %>% 
  dplyr::select(EnsemblID, all_of(rownames(sample_info)))

rowTumor_pol <- read_tsv("data/Sauyeun_PDX/rawTumor_Pol.tsv") %>%
  dplyr::select(EnsemblID, all_of(rownames(sample_info)))

# Translate to geneID
## get gene ID with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)

annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = ensembl75)

translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

rowTumor_cyt_gn <- dplyr::mutate(rowTumor_cyt, geneID = translate[EnsemblID]) %>%
  dplyr::select(geneID, matches("PDAC")) %>% drop_na(geneID) %>%
  distinct(geneID, .keep_all = T)

rowTumor_pol_gn <- dplyr::mutate(rowTumor_pol, geneID = translate[EnsemblID]) %>% 
  dplyr::select(geneID, matches("PDAC")) %>% drop_na(geneID) %>%
  distinct(geneID, .keep_all = T)

# Prepare data for anota2seq
## assert order
stopifnot(all(names(rowTumor_cyt_gn) == names(rowTumor_pol_gn)))
stopifnot(all(names(rowTumor_cyt_gn)[-1] == rownames(sample_info)))

phenoVec <- sample_info$ISRact
batchVec <- rownames(sample_info)

## PCA
tmpPCA <- inner_join(rowTumor_cyt_gn, rowTumor_pol_gn, by = "geneID", suffix = c("cyt", "pol")) %>%
  column_to_rownames("geneID") %>% as.matrix()


tmpPCANoZero<-tmpPCA[!apply(tmpPCA,1,FUN = function(x)any(x==0)),]
tmpPCANoNorm<-voom(calcNormFactors(DGEList(tmpPCANoZero)))$E
tmpSD<-apply(tmpPCANoZero,1,sd)
sdQuantiles<-quantile(tmpSD)
tmpPCASDFlit<-tmpPCANoNorm[tmpSD>sdQuantiles[4],]
tmpPCATans<-t(tmpPCASDFlit)
pcaOut<-prcomp(tmpPCATans)

### scaplot
var_explained <- pcaOut$sdev^2/sum(pcaOut$sdev^2)
fviz_eig(pcaOut, barfill = "lightgrey",
         barcolor = "black", title = "Insulinoma PCA",
         subtitle = paste("90% variance reached at the",
                          which(cumsum(var_explained)>0.9)[1],
                          "th component"))

### plots
pca_toplot <- as_tibble(pcaOut$x, rownames = NA) %>% 
  rownames_to_column() %>% 
  mutate(sample_name = str_remove_all(rowname, "cyt|pol"),
         fraction = if_else(str_detect(rowname, "cyt"), "cyt", "pol")) %>%
  inner_join(rownames_to_column(sample_info, "sample_name"), by ="sample_name") %>%
  column_to_rownames("rowname")

ggplot(pca_toplot, aes(x = PC1, y = PC3, color = ISRact, shape = fraction)) +
  geom_point()

plot_PCs(pca_toplot, "sample_name", 4, 2)
plot_PCs(pca_toplot, "fraction", 4, 2)
plot_PCs(pca_toplot, "ISRact", 4, 2)
##

ads <- anota2seqDataSetFromMatrix(
  dataP = column_to_rownames(rowTumor_pol_gn, "geneID"),
  dataT = column_to_rownames(rowTumor_cyt_gn, "geneID"),
  phenoVec = phenoVec,
  batchVec = NULL,
  dataType = "RNAseq",
  filterZeroGenes = TRUE, 
  normalize = TRUE,
  transformation = "TMM-log2",
  varCutOff = NULL)

ads <- anota2seqResidOutlierTest(ads)

## Analysis
ads <- anota2seqAnalyze(ads,
                        correctionMethod = "BH", useProgBar = TRUE, fileStem = "results/Anota2seq/",
                        analysis = c("translation", "buffering", "translated mRNA", "total mRNA"))

ads<-anota2seqSelSigGenes(ads,maxPAdj = 0.15,
                          selDeltaPT = log(1.2),
                          selDeltaTP = log(1.2),
                          selDeltaP = 0,
                          selDeltaT = 0,
                          minSlopeTranslation = -1,
                          maxSlopeTranslation = 2,
                          minSlopeBuffering = -2,
                          maxSlopeBuffering =  1)
## Result Plots
par(mfrow = c(1, 2))
anota2seqPlotPvalues(ads, selContrast = 1, useRVM = TRUE, plotToFile = FALSE, contrastName = paste(unique(phenoVec), collapse = " vs. "))

ads <- anota2seqRegModes(ads, c(TRUE, TRUE))
par(mfrow = c(1, 1))
anota2seqPlotFC(ads, selContrast = 1, plotToFile = FALSE,  contrastName = paste(unique(phenoVec), collapse = " vs. "))
