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

## Translate to gene names
best_ica_gn <- best_ica
as_tibble(best_ica_gn$S, rownames="EnsemblID") %>% 
  dplyr::mutate(GeneName = ensembl_to_gene[EnsemblID]) -> best_ica_gn_ 

summarise(best_ica_gn_, nan_count = sum(is.na(GeneName)),
          dup_count = sum(duplicated(GeneName))) #0 NAs, 21Dups

best_ica_gn$S <- dplyr::mutate(best_ica_gn_,  GeneName = if_else(GeneName %in% best_ica_gn_$GeneName[duplicated(best_ica_gn_$GeneName)],
                                                                 paste(GeneName, EnsemblID, sep = "_"),
                                                                 GeneName)) %>%
  dplyr::select(- EnsemblID) %>% column_to_rownames("GeneName")

PDX_tumor_cyt_gn <- dplyr::mutate(PDX_tumor_cyt,
                                 GeneName = if_else(ensembl_to_gene[EnsemblID] %in% ensembl_to_gene[PDX_tumor_cyt$EnsemblID][duplicated(ensembl_to_gene[PDX_tumor_cyt$EnsemblID])],
                                                     paste(ensembl_to_gene[EnsemblID], EnsemblID, sep = "_"),
                                                    ensembl_to_gene[EnsemblID])) %>%
  dplyr::select(-EnsemblID) %>% 
  relocate(GeneName)

gw <- PlotGeneWeights(ica = best_ica_gn,
                      interest_IC = "IC.3",
                      expression = PDX_tumor_cyt_gn,
                      df_id = "GeneName",
                      n_genes  = 25,
                      column_annotation = column_annotation,
                      annotation_colors = annot_colors)

ggarrange(gw[["densplot"]] + rremove("xylab"), gw[["heatmap"]], heights = c(1.5, 10), widths = c(0.2,1),
          ncol = 2, nrow = 1)



#GSVA
library(msigdbr)
library(GSVA)

all_genesets <- msigdbr("Homo sapiens")
#all_genesets %>% distinct(gs_subcat) #Check options
gset_db = "CP:REACTOME"
gset_remove_prefix = "REACTOME_"

use_genesets <- all_genesets %>% filter(gs_subcat %in% c(gset_db))
msigdbr_list <- split(x = use_genesets$ensembl_gene, f = use_genesets$gs_name)
msigdb_descriptions <- use_genesets[c("gs_name", "gs_description")] %>%
  unique() %>% column_to_rownames("gs_name")

gsvaRes <- gsva(data.matrix(best_ica$S), msigdbr_list, min.sz = 15)

# best for IC.4
gsvaRes[order(gsvaRes[,"IC.3"]),]

gsvaTop <- as_tibble(gsvaRes, rownames = "gene_set") %>% 
  mutate(the_rank = rank(-IC.3, ties.method = "random"),
         #gene_set = if_else(str_count(gene_set, "_") < 10, gene_set, signature_dict[gene_set]),
         gene_set = str_remove(gene_set, gset_remove_prefix),
         gene_set = fct_reorder(gene_set, the_rank,.desc = T)) %>%
  pivot_longer(cols = -c(gene_set, the_rank), names_to = "component", values_to = "ES") %>% 
  dplyr::filter(the_rank < 15 | the_rank > (nrow(gsvaRes)-15)) %>% 
  mutate(component = if_else(component == "IC.3", "IC.3", "Other")) %>% 
  dplyr::select(!c(the_rank))

gsva_IC3 <- ggplot(gsvaTop, aes(x = ES, y = gene_set)) + 
  geom_point(aes(alpha = if_else(component == "IC.3", 0.9, 0.3),
                 color = if_else(ES > 0, "blue", "red"))) +
  theme_bw() +
  labs(title = "IC.3", subtitle = paste0("Best ",gset_db, " gene sets")) +
  rremove("legend") +
  rremove("ylab")
