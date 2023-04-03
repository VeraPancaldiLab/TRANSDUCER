library(tidyverse)
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


## initial boxplot
raw_box <- log2(column_to_rownames(counts, "Geneid")+1)  %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample_name", values_to = "expression") %>%
  left_join(manip_info, "Sample_name") %>%
  ggplot(aes(x = Sample_name, y = expression, fill = Condition)) +
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
  pivot_longer(cols = -Gene, names_to = "Sample_name", values_to = "expression") %>%
  left_join(manip_info, "Sample_name") %>%
  ggplot(aes(x = Sample_name, y = expression, fill = Condition)) +
  geom_boxplot() +
  ggtitle("log2Filt") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Normalize
norm_tmp <-
    limma::voom(edgeR::calcNormFactors(edgeR::DGEList(filt_tmp)))$E ## TMM-log2 from anota2seqNormalize()
  

norm_box <- as_tibble(norm_tmp, rownames = "Gene")  %>%
  pivot_longer(cols = -Gene, names_to = "Sample_name", values_to = "expression") %>%
  left_join(manip_info, "Sample_name") %>%
  ggplot(aes(x = Sample_name, y = expression, fill = Condition)) +
  geom_boxplot() +
  ggtitle("log2TMM") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# comparison Boxplots
ggarrange(raw_box + rremove("x.text"),
          filt_box + rremove("x.text"),
          norm_box,
          ncol = 1, nrow = 3, heights = c(1,1,2), common.legend = TRUE, legend = "right")
