#Test the ICA toolkit on CPTAC-PDAC cohort

library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(readxl)

setwd("~/Desktop/Internship/human_cohort_analysis/Data/")
source(file = "~/Desktop/Internship/human_cohort_analysis/ICA_Toolkit-main/R/functions.R")


#Load data
cptac_raw <- read_delim("CPTAC-PDAC-raw-RSEM-expected-counts-gene-level.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE) %>%
  column_to_rownames("idx")%>%
  dplyr::select(contains("tumor")) %>% 
  rename_all(funs(str_replace_all(., "_tumor", ""))) %>%
  rename_all(~str_replace_all(.,"-", ".")) %>%
  rownames_to_column("gene")

#Import PC weight for human and pdx samples
human_PC1 <- read_csv("human_PC1.csv")

#Clinical data sheet for human cohort data
clinical_data <- read_excel("mmc1.xlsx", 
                            sheet = "Clinical_data") %>%
  mutate(follow_up_days = as.numeric(follow_up_days)) %>%
  mutate(status = ifelse(vital_status == "Deceased", 2, 1))

#Molecular phenotype data sheet for human cohort data
Molecular_phenotype_data <- read_excel("mmc1.xlsx", 
                                       sheet = "Molecular_phenotype_data") %>% 
  mutate_at(vars(immune_deconv:`necrosis_(%OF_TUMOR_WITH_NECROSIS)_histology_estimate`, KRAS_VAF), as.numeric) %>%
  mutate_at(vars(Bailey:Moffitt), as.factor) #change type to avoid errors

#load sample information for annotations
metadata <- dplyr::rename(human_PC1, case_id = sample) %>%
  inner_join(clinical_data, by = "case_id") %>%
  inner_join(Molecular_phenotype_data, by = "case_id") %>%
  mutate(Diabetes = ifelse(str_detect(medical_condition, "Diabetes"), "yes", "no")) %>%
  dplyr::rename(bcr_patient_barcode = case_id) %>%
  mutate(bcr_patient_barcode = str_replace_all(bcr_patient_barcode, "-", "."))

#Remove low tumor purity sample
cptac <- cptac_raw %>% 
  dplyr::select(gene, metadata$bcr_patient_barcode)

##1 Filter the sexual chromosomes genes to avoid getting a sex component
annot <- AnnotationDbi::select(org.Hs.eg.db, keys=cptac$gene, columns=c("ENSEMBL", "MAP"), keytype="SYMBOL")
to_keep <- dplyr::filter(annot, !grepl('X|Y', MAP)) %>%
  drop_na(MAP) %>%
  dplyr::select("SYMBOL") %>%
  distinct()

cptac_gn <- dplyr::filter(cptac, gene %in% to_keep$SYMBOL)

#2 Mean center genewise and selection of half most variable genes. Finally save an .RDS object as a checkpoint
cptac__ <- cptac_gn %>% pivot_longer(-gene, 'variable', 'value') %>%
  pivot_wider(id_cols = variable, names_from = gene) %>%
  mutate(across(where(is.numeric), ~ . - mean(.))) %>%
  pivot_longer(-variable, "gene",  "value") %>% 
  pivot_wider(id_cols = gene, names_from = variable)

abdv <- apply(cptac__ %>% dplyr::select(-gene), 1, function(x) { #! could be piped with  add_row(variable = "mean", !!! colMeans(.[-1])) when transposed
  sum(
    abs(
      x - mean(x)
    )
  ) / length(x)
})

median_abdv <- median(abdv)
cptac_p <- cptac__[abdv > median_abdv, ]

# plots with A1CF as example gene
pivot_longer(cptac_gn, -gene) %>% 
  dplyr::filter(gene == "A1CF") %>% 
  ggplot(aes(x = value)) + geom_density()

pivot_longer(cptac_p, -gene) %>% 
  dplyr::filter(gene == "A1CF") %>% 
  ggplot(aes(x = value)) + geom_density()

write_rds(cptac_p, file = "cptac_icaready.RDS")

#Find the most robust number of components
boot_res <- Boot_ICA(expression = cptac_p, df_id  =  "gene", range.comp = 2:20, iterations = 20, seed = 0)

ggplot(boot_res) + aes(x = nc, y = correlation, fill = bootstrap) +
  geom_boxplot(width = 0.5) +
  labs(y = "absolute pearson correlation", x = "number of components") +
  coord_cartesian(ylim = c(0.8, 1)) +
  theme_bw() +
  rotate_x_text(45)

#Find Molecular signature
icas_list <- Range_ICA(cptac_p, "gene", 2:20)

Best_nc(icas_list,
        range.comp = 2:20,
        metadata = metadata,
        metadata_id = "bcr_patient_barcode",
        vars = c("PC1"),
        is.categorical = F)

#Exploration of the component space
best_ica <- Range_ICA(cptac_p, "gene", 19)

annotations <- metadata %>%
  dplyr::select("bcr_patient_barcode", "PC1", "sex")

Sampleweights_distribution(ica = best_ica,
                           df = annotations,
                           df_id = "bcr_patient_barcode")

## clinical data
metadata_cont <- colnames(dplyr::select(metadata, where(is.numeric)))

metadata_disc <- c("tumor_stage_pathological", 'PC1status', "sex", "participant_country", 
                   "Bailey", "Collisson", "Moffitt", "Diabetes")

### find IC responsible
clinical_data_plots <- ICA_explorator(
  ica = best_ica,
  df = metadata,
  df_cont = metadata_cont,
  df_disc = metadata_disc,
  df_id = "bcr_patient_barcode")

clinical_data_plots$cont  

clinical_data_plots$disc

library(dorothea)
# Generate TF act 
### Select the regulon data
data(dorothea_hs_pancancer, package = "dorothea")
regulons <- dorothea_hs_pancancer %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

## VIPER
minsize = 5
ges.filter = FALSE

tf_act__<- dorothea::run_viper(column_to_rownames(cptac, "gene"), regulons,
                               options =  list(minsize = minsize, eset.filter = ges.filter, 
                                               cores = 1, verbose = FALSE, nes = TRUE))

# Format for analysis
tf_act_ <- as_tibble(tf_act__, rownames = "TF") %>%
  pivot_longer(cols = -TF, names_to = "bcr_patient_barcode") %>%
  pivot_wider(id_cols = c(bcr_patient_barcode), names_from = TF)

# get 50 most variable
tf_50_var <- summarise(tf_act_, across(where(is.double), var)) %>%
  pivot_longer(cols = everything(), names_to = "TF", values_to = "var") %>%
  dplyr::arrange(desc(var)) %>% dplyr::slice_head(n=50)

tf_act <- tf_act_ %>% dplyr::select(bcr_patient_barcode, all_of(tf_50_var$TF))
### find IC responsible
TF_plots <- ICA_explorator(
  ica = best_ica,
  df = tf_act,
  df_cont = tf_50_var$TF,
  df_id = "bcr_patient_barcode")

TF_plots$cont

#In depth dimension analysis
ic14_diabetes <- Sampleweights_indepth(ica = best_ica,
                                     interest_IC = "IC.14",
                                     df = metadata,
                                     df_id = "bcr_patient_barcode",
                                     var = "Diabetes",
                                     test = "ttest")

ic14_diabetes + 
  geom_violin(aes(fill = Diabetes), position=position_dodge(0.8), width=0.5) +
  geom_boxplot(position=position_dodge(0.8), width=0.1) +
  theme_classic()

ic14_PC1 <- Sampleweights_indepth(ica = best_ica,
                                 interest_IC = "IC.14",
                                 df = metadata,
                                 df_id = "bcr_patient_barcode",
                                 var = "PC1",
                                 test = "spearman")

ic14_PC1 + geom_point() + geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw()

#Plot gene weights
gw <- PlotGeneWeights(ica = best_ica,
                      interest_IC = "IC.14",
                      expression = cptac_gn,
                      df_id = "gene",
                      n_genes  = 25,
                      column_annotation = NA)

ggarrange(gw[["densplot"]] + rremove("xylab"), gw[["heatmap"]], heights = c(1.5, 10), widths = c(0.2,1),
          labels = c("A", "B"),
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
