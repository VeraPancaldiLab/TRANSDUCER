library(tidyverse)
library(ggpubr)
library(factoextra)
library(RColorBrewer)
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

setwd("/home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/00_Data/")

# Data loading
## Load counts (and name repair)
counts <- read_tsv("240323_Alignment/featureCountsOUTPUT.csv",
         comment = "#") %>% 
  dplyr::select(-all_of(c("Chr", "Start", "End", "Strand", "Length"))) %>% 
  dplyr::rename_with( ~ str_remove_all(.,"output_bams/|_Aligned.sortedByCoord.out.bam|__Aligned.sortedByCoord.out.bam"))

### Export raw counts  
write_tsv(counts, "rawcounts.tsv")

## load metadata
manip_info <- read_tsv("Data_RNA_sample_Jacobo.tsv")
multiqc <- read_tsv("240323_FASTQC/02_Output/qc_summary.tsv")  
STAR_align <- read_tsv("240323_Alignment/output_bams/STAR_info.tsv")
picard_metrics <- read_tsv("240323_PicardTools/02_Output/picard_rnaseqmetrics_assignment_plot.tsv")

sample_info <- left_join(manip_info, multiqc, by = "sample_name") %>%
  left_join(STAR_align, by = "sample_name") %>%
  left_join(picard_metrics, by = "sample_name")

## initial boxplot
raw_box <- log2(column_to_rownames(counts, "Geneid")+1)  %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "sample_name", values_to = "expression") %>%
  left_join(manip_info, "sample_name") %>%
  ggplot(aes(x = sample_name, y = expression, fill = Condition)) +
  geom_boxplot() +
  ggtitle("log2RAW") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Normalization (as Anota2seq TMM)
## Separate names by fraction, and put universal names
names(counts) %>%
  str_detect("_Input") %>%
  names(counts)[.] -> names_total

names(counts) %>%
  str_detect("_F8") %>%
  names(counts)[.] -> names_polysome

### 2 polysomes missing
not_polysome <- !str_remove(names_total, "_Input") %in% str_remove(names_polysome, "_F8")
print(paste0(names_total[not_polysome], " is missing the polysome fraction"))

## Remove 0 genes
filt_tmp <- column_to_rownames(counts, "Geneid") %>%
    .[!(apply(., 1, function(x) {
      any(x == 0)
    })), ] # from anota2seqRemoveZeroSamples(). remove any 0 gene


filt_box <- log2(filt_tmp) %>%
  as_tibble(rownames = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "sample_name", values_to = "expression") %>%
  left_join(manip_info, "sample_name") %>%
  ggplot(aes(x = sample_name, y = expression, fill = Condition)) +
  geom_boxplot() +
  ggtitle("log2Filt") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Normalize
norm_tmp <-
    limma::voom(edgeR::calcNormFactors(edgeR::DGEList(filt_tmp)))$E ## TMM-log2 from anota2seqNormalize()
  

norm_box <- as_tibble(norm_tmp, rownames = "Gene")  %>%
  pivot_longer(cols = -Gene, names_to = "sample_name", values_to = "expression") %>%
  left_join(manip_info, "sample_name") %>%
  ggplot(aes(x = sample_name, y = expression, fill = Condition)) +
  geom_boxplot() +
  ggtitle("log2TMM") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Exploratory Analysis
## comparison Boxplots
ggarrange(raw_box + rremove("x.text"),
          filt_box + rremove("x.text"),
          norm_box,
          ncol = 1, nrow = 3, heights = c(1,1,2), common.legend = TRUE, legend = "right")


## Exploratory PCA
norm_pca <- t(norm_tmp) %>%
  prcomp(scale=T) 

### scaplot
var_explained <- norm_pca$sdev^2/sum(norm_pca$sdev^2)
fviz_eig(norm_pca, barfill = "lightgrey",
         barcolor = "black", title = "Insulinoma PCA",
         subtitle = paste("90% variance reached at the",
                          which(cumsum(var_explained)>0.9)[1],
                          "th component"))

### plots

pca_toplot <- merge(norm_pca$x,
                    column_to_rownames(sample_info, "sample_name"),
                    by="row.names")

plot_PCs(pca_toplot, "Fraction", 3, 2)
plot_PCs(pca_toplot, "CAF", 3, 2)
plot_PCs(pca_toplot, "Batch", 3, 2)
plot_PCs(pca_toplot, "Condition", 10, 0.6)
plot_PCs(pca_toplot, "Experimentalist", 3, 2)
