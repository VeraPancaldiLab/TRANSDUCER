#Test the ICA toolkit on CPTAC-PDAC cohort

library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(readxl)

setwd("~/Documents/02_TRANSDUCER/06_Human_Cohort_Projection/01_PDXTranslation_to_PDXTrascription/")
source(file = "../The-Molecular-Signature-Generator/R/functions.R")


#Load data
PDX_tumor_cyt <- read_tsv("~/Documents/02_TRANSDUCER/02_PDX_stroma/00_Data/Processed_data/normTumor_Cyt.tsv")


sample_info <- read_tsv("~/Documents/02_TRANSDUCER/02_PDX_stroma/00_Data/Processed_data/sample_info.tsv") %>%
  dplyr::select(!MarseilleID) %>%
  dplyr::mutate(Diabetes = as_factor(Diabetes)) %>%
  dplyr::rename(ISRact = ICA3) %>% 
  arrange(factor(sample, levels = names(PDX_tumor_cyt)[-1]))

# load annotation with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)

#listAttributes(ensembl75, page="feature_page")
annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id',
                                        'entrezgene',
                                        'chromosome_name'), mart = ensembl75)

gene_to_ensembl = deframe(annot_ensembl75[c( "external_gene_id", "ensembl_gene_id")])
ensembl_to_gene = setNames(names(gene_to_ensembl), gene_to_ensembl)
# filtering
## XY
annot_ensembl75 %>% dplyr::filter(!chromosome_name %in% c("X", "Y")) %>%
  pull(ensembl_gene_id) -> non_sex

PDX_tumor_cyt %>% dplyr::filter(EnsemblID %in% non_sex) -> PDX_tumor_cyt_

## mean centring gene wise
PDX_tumor_cyt_ %>% column_to_rownames("EnsemblID") %>%
  apply(1, function(x) x - mean(x)) %>% t() %>%
  data.frame() -> PDX_tumor_cyt__

## Inter quartile range (measure variability)
PDX_tumor_cyt__ %>% apply(1, IQR) -> iqrs
mostvar <- iqrs[iqrs > median(iqrs)]
PDX_tumor_cyt__ %>% as_tibble(rownames = "EnsemblID") %>%
  dplyr::filter(EnsemblID %in% names(mostvar)) -> PDX_tumor_cyt_icaready

# Check distribution variability and number of genes of the different steps


#Find Molecular signature
icas_list <- Range_ICA(PDX_tumor_cyt_icaready, "EnsemblID", 2:20)

Best_nc(icas_list,
        range.comp = 2:20,
        metadata = sample_info,
        metadata_id = "sample",
        vars = c("ISRact"),
        is.categorical = F)

#Exploration of the component space
best_ica <- Range_ICA(PDX_tumor_cyt_icaready, "EnsemblID", 14)

Sampleweights_distribution(ica = best_ica,
                           df = dplyr::select(sample_info, "sample", "ISRact", "PAMG", "Diabetes"),
                           df_id = "sample")

## sample info
sample_info_plots <- ICA_explorator(
  ica = best_ica,
  df = sample_info,
  df_cont = c("ISRact", "PAMG"),
  df_id = "sample")

sample_info_plots$cont

# technical info
multiqc <- read_tsv("~/Documents/02_TRANSDUCER/02_PDX_stroma/00_Data/Processed_data/multiQC_summary.tsv") %>%
  dplyr::filter(CITID %in% names(PDX_tumor_cyt)[-1], Fraction == "Cytosolic") %>%
  dplyr::select(-c(fastq_name, Fraction)) %>%
  dplyr::rename(sample = CITID)

sample_info_plots <- ICA_explorator(
  ica = best_ica,
  df = multiqc,
  df_cont = names(multiqc)[-1],
  df_id = "sample")

sample_info_plots$cont

#In depth dimension analysis
ic3_diabetes <- Sampleweights_indepth(ica = best_ica,
                                     interest_IC = "IC.3",
                                     df = sample_info,
                                     df_id = "sample",
                                     var = "Diabetes",
                                     test = "welch")

ic3_diabetes + 
  geom_boxplot(aes(fill = Diabetes), position=position_dodge(0.8), width=0.1) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 2, alpha = 0.5) +
  theme_classic()

ic3_ISRact <- Sampleweights_indepth(ica = best_ica,
                                 interest_IC = "IC.3",
                                 df = sample_info,
                                 df_id = "sample",
                                 var = "ISRact",
                                 test = "spearman")
ic3_ISRact + geom_point() + geom_smooth(method=lm) +
  theme_bw()

ic3_ISRact + geom_point(aes(color = Diabetes)) + geom_smooth(method=lm) +
  theme_bw()

#Plot gene weights
column_annotation <- dplyr::select(sample_info, c(sample, Diabetes, PAMG, ISRact)) %>%
  column_to_rownames("sample")

annot_colors <- list(class = c(`most possitive` = "red", `most negative` = "blue"),
                     PAMG = c("#FF7F00", "white", "#377DB8"),
                     ISRact = c("#FFFFCC", "#006837"),
                     Diabetes = c(`0` = "#F8766D" , `1` = "#00BFC4", `NA` = "#7F7F7F")) 

gw <- PlotGeneWeights(ica = best_ica,
                      interest_IC = "IC.3",
                      expression = PDX_tumor_cyt,
                      df_id = "EnsemblID",
                      n_genes  = 25,
                      column_annotation = column_annotation,
                      annotation_colors = annot_colors)

ggarrange(gw[["densplot"]] + rremove("xylab"), gw[["heatmap"]], heights = c(1.5, 10), widths = c(0.2,1),
          ncol = 2, nrow = 1)



#GSVA
library(msigdbr)
library(GSVA)

all_genesets <- msigdbr("Homo sapiens")
use_genesets <- dplyr::filter(all_genesets, gs_cat %in% c("H", "C2"))
msigdbr_list <- split(x = use_genesets$gene_symbol, f = use_genesets$gs_name)
msigdb_descriptions <- use_genesets[c("gs_name", "gs_description")] %>%
  unique() %>% column_to_rownames("gs_name")

gsvaRes <- gsva(data.matrix(best_ica$S), msigdbr_list, min.sz = 15)

# best for IC.14
gsvaTop <- as_tibble(gsvaRes, rownames = "gs_name") %>%  mutate(the_rank = rank(IC.14, ties.method = "random")) %>% 
  dplyr::filter(the_rank < 26 | the_rank > (nrow(gsvaRes)-26)) %>% 
  dplyr::mutate(gs_name = fct_reorder(gs_name, IC.14)) %>% mutate(sign = case_when(IC.14 < 0 ~ "neg", IC.14 > 0 ~ "pos")) %>%
  dplyr::select(gs_name, IC.14, sign)

ggplot(gsvaTop, aes(x = IC.14, y = gs_name, color = sign)) + 
  geom_point() +
  ggtitle("Best H and C2 gene sets IC.14") +
  theme_bw() +
  rremove("legend") +
  rremove("ylab")
