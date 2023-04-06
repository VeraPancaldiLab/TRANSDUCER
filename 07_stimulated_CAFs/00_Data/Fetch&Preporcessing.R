library(tidyverse)
library(ggpubr)
library(factoextra)
library(RColorBrewer)
library(sva)
library(pheatmap)
################################################################################
################################PARAMETERS######################################
filter_samples = "17AC" # NULL | 17AC | 02136
correct_batch = T # Should correct for batch effect?
sample_sample_corrplot_annot = "picard_metrics" # manip_info | picard_metrics | tech_info | STAR_info
################################################################################
#################################FUNCTIONS######################################

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

# Customizable corrplot (modified from https://www.khstats.com/blog/corr-plots/corr-plots/)
cors <- function(df, cor.stat) {
  M <- Hmisc::rcorr(as.matrix(df), type = cor.stat)
  Mdf <- map(M, ~data.frame(.x))
  return(Mdf)
}

formatted_cors <- function(df, cor.stat){
  cors(df, cor.stat) %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, names_to = "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    dplyr::rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
}

################################################################################
###################################MAIN#########################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/00_Data/")

# Data loading
## Load counts (and name repair)
counts <- read_tsv("240323_Alignment/featureCountsOUTPUT.csv",
         comment = "#") %>% 
  dplyr::select(-all_of(c("Chr", "Start", "End", "Strand", "Length"))) %>% 
  dplyr::rename_with( ~ str_remove_all(.,"output_bams/|_Aligned.sortedByCoord.out.bam|__Aligned.sortedByCoord.out.bam"))

### Export unfiltered raw counts for further analyses
write_tsv(counts, "rawcounts.tsv")

### filter to a subset of samples
if (!is.null(filter_samples)){
  counts <- dplyr::select(counts, Geneid, names(counts)[str_detect(names(counts), pattern = filter_samples)])
}

## load metadata
manip_info <- read_tsv("Data_RNA_sample_Jacobo.tsv")
multiqc <- read_tsv("240323_FASTQC/02_Output/qc_summary.tsv")  
STAR_align <- read_tsv("240323_Alignment/output_bams/STAR_info.tsv")
picard_metrics <- read_tsv("240323_PicardTools/02_Output/picard_rnaseqmetrics_assignment_plot.tsv")

sample_info <- left_join(manip_info, multiqc, by = "sample_name") %>%
  left_join(STAR_align, by = "sample_name") %>%
  left_join(picard_metrics, by = "sample_name") %>%
  dplyr::filter(sample_name %in% names(counts)) %>%
  dplyr::mutate(sample_name = fct(sample_name, levels=names(counts)[-1])) %>%
  dplyr::arrange(sample_name)

all(names(counts)[-1] == sample_info$sample_name) %>% stopifnot()

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

## Batch effect correction
if (correct_batch == F){
  filt_tmp -> filt_tmp
  batch <- "non corrected"
  
} else if (correct_batch == T){
  
  filt_tmp <- ComBat_seq(as.matrix(filt_tmp),
                         batch = sample_info$Batch,
                         full_mod = TRUE,
                         covar_mod = cbind(as.numeric(factor(sample_info$Condition)), 
                                           as.numeric(factor(sample_info$Fraction)))) 
  
  batch <- "batch corrected"
}

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
plot_PCs(pca_toplot, "Condition", 5, 1.6)
plot_PCs(pca_toplot, "Condition", 3, 2)
plot_PCs(pca_toplot, "Experimentalist", 3, 2)

#### custom plot
highlight = c("Batch_A_17AC_FAKi_Input", "Batch_A_17AC_TGF_F8")
as_tibble(pca_toplot) %>% 
  mutate(highlight = fct(if_else(Row.names %in% highlight, Row.names, "other"))) %>%
  ggplot(aes(x=PC1, y=PC2, color = highlight)) +
  geom_point() +
  theme_pubr()

## Exploratory corrplots
### sample vs sample
if (sample_sample_corrplot_annot == "manip_info") {
  annot <- dplyr::select(sample_info, sample_name, Batch, CAF, Condition, Experimentalist, Fraction) %>%
    column_to_rownames("sample_name")
  } else if (sample_sample_corrplot_annot == "tech_info") {
  annot <- dplyr::select(sample_info, sample_name, ng_ul, R260_280, R260_230, RQN, RQN_Integragen, ng_ul_Integragen) %>%
    column_to_rownames("sample_name")
  } else if (sample_sample_corrplot_annot == "STAR_info") {
    annot <- dplyr::select(sample_info, sample_name, Number_of_input_reads, Uniquely_mapped_reads_perc, average_mapped_input_diff_, Mismatch_rate_per_base_perc, Deletion_rate_per_base, Insertion_rate_per_base, perc_of_reads_mapped_to_multiple_loci, perc_of_reads_mapped_to_too_many_loci) %>%
      column_to_rownames("sample_name")
  } else if (sample_sample_corrplot_annot == "picard_metrics") {
  annot <- dplyr::select(sample_info, sample_name, Coding, UTR, Intronic, Intergenic, Ribosomal, PF_not_aligned) %>%
    column_to_rownames("sample_name")
}

formatted_cors(norm_tmp, "pearson") %>%
  dplyr::select(measure1, measure2, r) %>%
  pivot_wider(id_cols = measure1, names_from = measure2, values_from = r) %>%
  column_to_rownames("measure1") %>%
  pheatmap(annotation_col = annot)

### PCA vs contmetadata
column_to_rownames(pca_toplot, "Row.names") %>%
  dplyr::select(where(is.numeric)) %>%
  formatted_cors("spearman") %>%
  dplyr::filter(!str_detect(measure1, "PC")) %>%
  dplyr::filter(str_detect(measure2, "PC.$")) %>% 
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

  