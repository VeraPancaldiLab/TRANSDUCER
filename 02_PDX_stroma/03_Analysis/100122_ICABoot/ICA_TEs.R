#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
library(JADE)
library(Hmisc)
library(pheatmap)
library(msigdbr)
library(GSVA)
library(ggplotify)
library(ggpubr)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/")
source("functions.R")

# PARAMETERS
#-------------------------------------------------------------------------------
# Bootstrap
run_boot <- FALSE
range.comp <- 2:15 # when ncomp is =< df Warning: In sqrt(puiss[rangeW]) : NaNs produced
boot.iter <- 500
boot.perc <- 0.95

# Analysis
elected_ncomp <- 6 # 4 if looking at distribution, 6 for standar, like in tumour way
component_reorientation = TRUE
reorient <- c(1, 1, -1, 1, 1, -1)
#-------------------------------------------------------------------------------

# DATA LOADING/PROCESSING
#-------------------------------------------------------------------------------
# load TEs
TEs <- read_tsv("../081221_TranslationEfficacy/02_Output/TEs.tsv")
sample_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>%
  column_to_rownames("sample") %>% .[colnames(TEs)[-1],]

sample_info$Diabetes <- as_factor(sample_info$Diabetes)

# load annotation with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "mmusculus_gene_ensembl",
                        version = 75)

#listAttributes(ensembl75, page="feature_page")
annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                          'external_gene_id',
                          'entrezgene',
                          'mgi_id',
                          'chromosome_name'), mart = ensembl75)

# filtering
## XY
annot_ensembl75 %>% dplyr::filter(!chromosome_name %in% c("X", "Y")) %>%
  pull(ensembl_gene_id) -> non_sex

TEs %>% dplyr::filter(EnsemblID %in% non_sex) -> TEs_

## mean centring gene wise
TEs_ %>% column_to_rownames("EnsemblID") %>%
  apply(1, function(x) x - mean(x)) %>% t() %>%
  data.frame() -> TEs__

## Inter quartile range (measure variability)
TEs__ %>% apply(1, IQR) -> iqrs
mostvar <- iqrs[iqrs > median(iqrs)]
TEs__ %>% rownames_to_column("EnsemblID") %>%
  dplyr::filter(EnsemblID %in% names(mostvar)) %>%
  column_to_rownames("EnsemblID") -> TEs_icaready

# ICA Bootstrapping
if (run_boot == TRUE){
  ## Baseline
  TEs_icaready %>% jade_range(range.comp, MARGIN = 1) -> base_res_gene
  TEs_icaready %>% jade_range(range.comp, MARGIN = 2) -> base_res_sample

  ## Bootstrap
  gene_boot <- jade_choosencom(TEs_icaready, base_res_gene,
                               MARGIN = 1,
                               iterations = boot.iter,
                               seed = 0,
                               perc = boot.perc
  )
  
  sample_boot <- jade_choosencom(TEs_icaready, base_res_sample,
                                 MARGIN = 2,
                                 iterations = boot.iter,
                                 seed = 0,
                                 perc = boot.perc
  )
  
  boot_plots(s_boot = sample_boot, g_boot = gene_boot, name = "TEs")
}

# Most robust ICA analysis
jade_result <- JADE(TEs_icaready, n.comp = elected_ncomp)
colnames(jade_result[["A"]]) <- paste("IC", 1:elected_ncomp, sep = ".")
rownames(jade_result[["A"]]) <- names(jade_result$Xmu)
write_rds(jade_result, "02_Output/ICA_TEs.RDS") # never save a reoriented ICA analysis

## component reorientation
if (component_reorientation == TRUE){
  stopifnot("vector of reorientation should have the same length as number of components" = length(reorient)==elected_ncomp)
  jade_result[["A"]] <- t(t(jade_result[["A"]]) * reorient)
  jade_result[["S"]] <- t(t(jade_result[["S"]]) * reorient)
}
A_mat <- as.data.frame(jade_result[["A"]])
S_mat <- as.data.frame(jade_result[["S"]])
#-------------------------------------------------------------------------------

# EXPLORATORY ANALYSIS
#-------------------------------------------------------------------------------

# Sample weight analysis
## Metadata
annotations <- sample_info[-1]
stopifnot(rownames(A_mat) == rownames(annotations))
annotations %>% dplyr::select(!Diabetes) %>%
  names() -> cont_names

## plot
plot_sample_weights(A_mat, annotations, cont_names, "sampleweights_TEs")

## Export for further correlations
stopifnot(rownames(A_mat) == rownames(annotations))
bind_cols(A_mat, annotations[,c("ICA3", "PAMG")]) %>%
  dplyr::rename(ISRact = ICA3) -> complete_annotation

## QC and alignment info
multiqc <- read_tsv("../../00_Data/Processed_data/multiQC_summary.tsv") %>%
  dplyr::filter(CITID %in% names(TEs)[-1])


seq_info <- read_tsv("../../00_Data/Processed_data/sequencing_info.tsv") %>%
  dplyr::filter(CITID %in% names(TEs)[-1])

## Cyt read info
multiqc_cyt <- dplyr::filter(multiqc, Fraction == "Cytosolic") %>%
  dplyr::select(-c(fastq_name, Fraction))

multiqc_cyt_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(multiqc_cyt) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(multiqc_cyt), measure2 %in% names(A_mat)) %>% #not square corr
  mutate(r = abs(r)) %>% #abscorr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICATEs ~ multiqc",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(0,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ggpubr::rotate_x_text(angle = 90)


seq_info_cyt <- dplyr::filter(seq_info, Fraction == "Cytosolic") %>%
  mutate(input2maped_length_diff = Average_input_read_length - Average_mapped_length) %>%
  dplyr::select(-c(fastq_id, Fraction, Average_input_read_length, Average_mapped_length))

star_cyt_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(seq_info_cyt) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(seq_info_cyt), measure2 %in% names(A_mat)) %>% #not square corr
  mutate(r = abs(r)) %>% #abscorr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICATEs ~ STAR metadata",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(0,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ggpubr::rotate_x_text(angle = 90)

techi_plots_cyt <- ggarrange(multiqc_cyt_corr, star_cyt_corr, align = "h", common.legend = T, legend = "right")

## Polysome (F8-9) read info
multiqc_pol <- dplyr::filter(multiqc, Fraction == "Polysome") %>%
  dplyr::select(-c(fastq_name, Fraction))

multiqc_pol_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(multiqc_pol) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(multiqc_pol), measure2 %in% names(A_mat)) %>% #not square corr
  mutate(r = abs(r)) %>% #abscorr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICATEs ~ multiqc",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(0,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ggpubr::rotate_x_text(angle = 90)

seq_info_pol <- dplyr::filter(seq_info, Fraction == "Polysome") %>%
  mutate(input2maped_length_diff = Average_input_read_length - Average_mapped_length) %>%
  dplyr::select(-c(fastq_id, Fraction, Average_input_read_length, Average_mapped_length))

star_pol_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(seq_info_pol) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(seq_info_pol), measure2 %in% names(A_mat)) %>% #not square corr
  mutate(r = abs(r)) %>% #abscorr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICATEs ~ STAR metadata",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(0,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ggpubr::rotate_x_text(angle = 90)

techi_plots_pol <- ggarrange(multiqc_pol_corr, star_pol_corr, align = "h", common.legend = T, legend = "right")

## Export
techi_plots <- ggarrange(techi_plots_cyt, techi_plots_pol,
                         labels = c("Cytosolic", "Polysome"),
                         font.label = list(face = "bold", color = "black", size = 20),
                         vjust=-1,  nrow = 2, align = "v")

annotate_figure(techi_plots,
                top = text_grob("", size = 40))

ggsave("02_Output/techdata_corr_TEs.pdf", dpi = "print", width = 210, height = 300, units = "mm", scale = 1.25)

## RNAseq celltype deconvolution
mMCPcounter_res <-  read_tsv("../180122_Various/02_Output/mMCPcounter_results.tsv") %>% column_to_rownames("cell_types") %>% t() %>% as_tibble(rownames = "samples") %>% column_to_rownames("samples") %>% .[rownames(A_mat),]
ImmuCC_res <-  read_tsv("../180122_Various/02_Output/ImmuCC_results.tsv") %>% column_to_rownames("samples") %>% .[rownames(A_mat),]

Plot_deconv(ImmuCC_res, complete_annotation, "ImmuCC_TEs")
Plot_deconv(mMCPcounter_res, complete_annotation, "mMCPcounter_TEs")

## TF activity
tumour_tf <- read_tsv("../180122_Various/02_Output/TFact_tumor_PanCan.tsv") %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))
stroma_tf <- read_tsv("../180122_Various/02_Output/TFact_stroma_Gtex.tsv") %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))

Plot_general_TFs(tumour_tf, "TF_tumor_vs_TEs_ICA", 25, complete_annotation)
Plot_general_TFs(stroma_tf, "TF_stroma_vs_TEs_ICA", 25, complete_annotation)
PlotBestCorr(complete_annotation, tumour_tf, 10, "best_TF_tumor_vs_TEs")
PlotBestCorr(complete_annotation, stroma_tf, 10, "best_TF_stroma_vs_TEs")

# Gene weight analysis
## translate to gene names
translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

## plot
PlotGeneWeights(S_mat, TEs, 25, translate, complete_annotation, analysis_name = "gene_weights_TEs")

## GSVA of best genes
### gene set preparation
all_genesets <- msigdbr("Mus musculus")
all_genesets %>% filter(gs_subcat %in% c("CP:REACTOME")) -> use_genesets
msigdbr_list = split(x = use_genesets$ensembl_gene, f = use_genesets$gs_name)

signature_dict <- dplyr::select(all_genesets, c(gs_name, gs_id)) %>%
  distinct(gs_id, .keep_all = T) %>%  deframe()

gsvaRes <- gsva(data.matrix(S_mat), msigdbr_list, min.sz = 15)

### Plot
for (comp in colnames(S_mat)){
  gsvaRes[order(gsvaRes[,comp]),]
  
  gsvaTop <- as_tibble(gsvaRes, rownames = "gene_set") %>% 
    mutate(the_rank = rank(-get(comp), ties.method = "random"),
           gene_set = if_else(str_count(gene_set, "_") < 10, gene_set, signature_dict[gene_set]),
           gene_set = str_remove(gene_set, "REACTOME_"),
           gene_set = fct_reorder(gene_set, the_rank,.desc = T)) %>%
    pivot_longer(cols = -c(gene_set, the_rank), names_to = "component", values_to = "ES") %>% 
    dplyr::filter(the_rank < 25 | the_rank > (nrow(gsvaRes)-25)) %>% 
    mutate(component = if_else(component == comp, comp, "Other")) %>% 
    dplyr::select(!c(the_rank))
  
  ggplot(gsvaTop, aes(x = ES, y = gene_set)) + 
    geom_point(aes(alpha = if_else(component == comp, 0.9, 0.3),
                   color = if_else(ES > 0, "blue", "red"))) +
    theme_bw() +
    labs(title = comp, subtitle = "Best Reactome gene sets") +
    rremove("legend") +
    rremove("ylab")
  ggsave(paste0("02_Output/gsva_TEs", comp, ".pdf"))
}

concat_pdf <- c("TEs_bootstrapICA_lineplot.pdf",
                "TEs_bootstrapICA_boxplot.pdf",
                "techdata_corr_TEs.pdf",
                "sampleweights_TEs.pdf",
                "ImmuCC_TEs.pdf",
                "mMCPcounter_TEs.pdf",
                "TF_analysis_TF_stroma_vs_TEs_ICA.pdf",
                "TF_analysis_TF_tumor_vs_TEs_ICA.pdf")

for (comp in 1:elected_ncomp){
  c(concat_pdf, paste("best_TF_stroma_vs_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
  c(concat_pdf, paste("best_TF_tumor_vs_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
  c(concat_pdf, paste("gene_weights_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
  c(concat_pdf, paste("gsva_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
}

#system(paste(c("cd 02_Output/ \n pdfunite", concat_pdf, "ICA_TEs.pdf"), collapse = " "))
system(paste(c("cd 02_Output/ \n convert", concat_pdf, "ICA_TEs.pdf"), collapse = " "))

# Network analysis
## Protein Protein Intreraction
MID <- read_csv("01_Input/MID2022.csv")

trans_mgi_name <- deframe(annot_ensembl75[c("mgi_id", "external_gene_id")])
trans_ensembl_name <- deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])
trans_name_ensembl <- deframe(annot_ensembl75[c("external_gene_id", "ensembl_gene_id")])

nodes <- as_tibble(S_mat, rownames = "ensembl_id") %>%
  mutate(id = trans_ensembl_name[ensembl_id]) %>%
  relocate(ensembl_id, .after = "id")

edges <- mutate(MID,
                gene1 = trans_mgi_name[gene1],
                gene2 = trans_mgi_name[gene2]) %>%
  dplyr::filter(gene1 %in% nodes$id, 
                gene2 %in% nodes$id) %>%
  dplyr::rename(source = gene1,
                target = gene2,
                interaction = type) %>%
  dplyr::select(source, target, interaction)

PlotNetwork(nodes, edges, S_mat,valid_comp = 1:6, main_name = "MID_ICAcyt")

## POSTARS3 mice RBP-target
POSTARS3 <- read_tsv("01_Input/POSTARS3/POSTARS3_mapped.tsv") %>%
  mutate(RBP_name = str_to_title(RBP_name),
         RBP_id = trans_name_ensembl[RBP_name])

### calculate node features 
cyt_m <- read_tsv("../../00_Data/Processed_data/normHost_Cyt.tsv")

RBPs_corrs <-  dplyr::filter(cyt_m, EnsemblID %in% POSTARS3$RBP_id) %>%
  mutate(gene_symbol = translate[EnsemblID]) %>%
  dplyr::select(gene_symbol, rownames(A_mat)) %>% 
  pivot_longer(-gene_symbol, names_to = "sample") %>%
  pivot_wider(values_from = value, names_from = gene_symbol) %>%
  inner_join(rownames_to_column(A_mat, "sample"), by = "sample") %>%
  column_to_rownames("sample") %>%
  formatted_cors(cor.stat = "spearman") %>%
  filter(measure1 %in% POSTARS3$RBP_name,
         measure2 %in% names(A_mat)) %>%
  mutate(FDR = p.adjust(p = p, method = "BH"),
         sig_FDR = FDR < 0.05)

absent_RBPs <- na.omit(deframe(distinct(POSTARS3, RBP_id))[!deframe(distinct(POSTARS3, RBP_id)) %in% cyt_m$EnsemblID])

RBPs_to_net <- dplyr::select(RBPs_corrs, measure1, measure2, r, sig_FDR) %>%
  pivot_wider(id_cols = measure1,
              names_from = measure2,
              values_from = c(r, sig_FDR),
              names_sep = "_") %>%
  dplyr::rename(id = measure1) %>% 
  dplyr::add_row(id = absent_RBPs) # Fix this to include all the RBPs

### Build Cytoscape object -> Function
#### nodes
nodes <- as_tibble(S_mat, rownames = "ensembl_id") %>%
  mutate(id = trans_ensembl_name[ensembl_id],) %>% 
  full_join(dplyr::select(POSTARS3, gene_name, gene_biotype) %>%
               distinct() %>%
               rename(id = gene_name)) %>%
  left_join(RBPs_to_net) %>%
  relocate(id, gene_biotype, .after = "ensembl_id") %>%
  distinct(id, .keep_all = T)

#### Is the gene best for any component?
thresholds <- apply(S_mat, 2, function(x) quantile(abs(x), 0.95))
bestgenes <- lapply(names(thresholds), function(x)
  rownames(dplyr::filter(S_mat, abs(get(x)) > thresholds[x])))
names(bestgenes) <- names(thresholds)

is_bestgene <- lapply(names(bestgenes), function(x) nodes$ensembl_id %in% bestgenes[[x]]) %>% 
  as_tibble(.name_repair= "universal" ) %>% 
  rename_all(~paste0("best_IC.",1:6)) %>%
  mutate(best_any = rowSums(.), 
         best_for = ifelse(best_any == 1, apply(., 1, function(x) colnames(.)[as.logical(x)][1]), ifelse(best_any == 0, 0, "many")))

nodes <- bind_cols(nodes, is_bestgene)

edges <- dplyr::filter(POSTARS3,
                       RBP_name %in% nodes$id, 
                       gene_name %in% nodes$id) %>%
  dplyr::rename(source = RBP_name,
                target = gene_name) %>%
  dplyr::select(source, target, experimentplussoftware, sample_origin)

### Clean and account for duplicated edges
edges_clean <- dplyr::group_by(edges, source, target, .add = TRUE) %>%
  summarise(n_peaks = n()) %>%
  dplyr::left_join(dplyr::select(nodes, id, best_for) %>%
                     dplyr::rename(target=id))
  
### Calculate general network measures
#### general and component specific out degree
general_node_stats <- dplyr::group_by(edges_clean, source, best_for, .add=T) %>%
  dplyr::summarise(out_degree=n(),
            total_peaks=sum(n_peaks)) %>% 
  pivot_wider(id_cols = source,
              names_from = best_for,
              values_from = c(out_degree, total_peaks),
              names_sep = "_") %>%
  ungroup() %>%
  dplyr::mutate_all(~ifelse(is.na(.), 0, .)) %>%
  dplyr::mutate(out_degree_total = purrr::reduce(dplyr::select(., starts_with("out_degree")), `+`),
                total_peaks_absolute = purrr::reduce(dplyr::select(., starts_with("total_peaks")), `+`)) %>%
  dplyr::rename(id=source) %>%
  right_join(nodes)

### Build network
network <- createNetworkFromDataFrames(general_node_stats, edges_clean, title = "POSTARS3")
setNodeLabelMapping("id")

### Fisher exact test for RBP n of experimentally detected targets
#### test 1 RBP 1 component (check)
contingency.table <- dplyr::filter(general_node_stats, id == "Apc") %>% 
  dplyr::select(id, out_degree_total, out_degree_best_IC.6)

A = contingency.table$out_degree_best_IC.6 #targets of AGO1 that are IC6
B = contingency.table$out_degree_total - A #targets of AGO1 that are not IC.6
C = sum(general_node_stats$best_IC.6 == T) - A # non targets of AGO1 that are IC6
D = nrow(general_node_stats) - (A+B) - C #Non targets of AGO1 that are not IC.6

test_contingency = matrix(data = c(A,B,C,D), nrow = 2, ncol = 2)
fisher_test(test_contingency)

#### full implementation in function #TBD use apply to make faster
contingency.table <- dplyr::filter(general_node_stats, id %in% RBPs_to_net$id) %>% 
  dplyr::select(id, out_degree_total, matches("out_degree_best")) %>%
  column_to_rownames("id")

fisher_results = tibble(IC=character(), RBP = character(), n=numeric(), p=numeric(), OR=numeric())
for (IC in paste("IC", 1:6, sep = ".")){
  total_IC_genes = sum(general_node_stats[[paste0("best_", IC)]] == T)
  subset_ct <- dplyr::select(contingency.table, out_degree_total, matches(IC)) %>%
    as.matrix()
  for (RBP in rownames(subset_ct)){
    A = subset_ct[RBP,2] #targets of AGO1 that are IC6
    B = subset_ct[RBP,1] - A #targets of AGO1 that are not IC.6
    C = total_IC_genes - A # non targets of AGO1 that are IC6
    D = nrow(general_node_stats) - (A+B) - C #Non targets of AGO1 that are not IC.6
  
  test_contingency = matrix(data = c(A,B,C,D), nrow = 2, ncol = 2)
  fisher <- fisher.test(test_contingency, conf.int = T, conf.level = 0.95, alternative = "greater")
  fisher_results = add_row(fisher_results,
                           IC = IC,
                           RBP = RBP,
                           p = fisher$p.value,
                           OR = fisher$estimate)
  }
}

##### plot of Fisher results
ggplot(fisher_results, aes(y = RBP, x = OR, color = IC, size = log10(p))) +
  geom_point() +
  scale_size(trans = 'reverse') +
  theme_bw()
#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: Clinical and technical data
#-------------------------------------------------------------------------------
columns_multiqc = c("CITID", "Dups_R1", "GC_R1")
multiqc_sp <- dplyr::inner_join(dplyr::select(multiqc_cyt, all_of(columns_multiqc)),
                                dplyr::select(multiqc_pol, all_of(columns_multiqc)),
                                by="CITID", suffix = c("_cyt", "_pol"))

clinical_technical_corrplot <- complete_annotation %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(multiqc_sp) %>%
  dplyr::select(CITID, matches("IC."), ISRact, PAMG, matches("Dups_R1"), matches("GC_R1")) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% c("ISRact", "PAMG", "Dups_R1_cyt", "GC_R1_cyt", "Dups_R1_pol", "GC_R1_pol"),
         measure2 %in% names(A_mat)) %>% #not square corr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAcyt ~ multiqc",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
  ggpubr::rotate_x_text(angle = 90)

ggsave(file="02_Output/Figures/clinical_technical_corrplot_TEs.svg", plot=clinical_technical_corrplot, width=4, height=6)
sed -i "s/ textLength='[^']*'//" file.svg
#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: MCPCounter deconvolution
#-------------------------------------------------------------------------------
deconvolution_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(as_tibble(mMCPcounter_res, rownames = "CITID")) %>%
  dplyr::select(CITID, names(A_mat), names(mMCPcounter_res)) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(mMCPcounter_res),
         measure2 %in% names(A_mat))

clust <- dplyr::select(deconvolution_corr, measure1, measure2, r) %>% 
  pivot_wider(names_from=measure2, values_from = r) %>%
  column_to_rownames("measure1") %>%
  dist() %>%
  hclust()

deconvolution_corrplot <- ggplot(deconvolution_corr, aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAcyt ~ MCPCounter",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0),limits = clust$labels[clust$order]) +
  scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
  ggpubr::rotate_x_text(angle = 90)

ggsave(file="02_Output/Figures/deconvolution_corrplot_TEs.svg", plot=deconvolution_corrplot, width=8, height=6)
#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: Most correlated stromal TF activities
#-------------------------------------------------------------------------------
stroma_mostvar_TF <- Get_mostvar(stroma_tf, 20) %>% rownames_to_column() %>% pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)

stroma_tf_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(as_tibble(t(stroma_tf), rownames = "CITID")) %>%
  dplyr::select(CITID, names(A_mat), names(stroma_mostvar_TF)[-1]) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(stroma_mostvar_TF)[-1],
         measure2 %in% names(A_mat))

clust <- dplyr::select(stroma_tf_corr, measure1, measure2, r) %>% 
  pivot_wider(names_from=measure2, values_from = r) %>%
  column_to_rownames("measure1") %>%
  dist() %>%
  hclust()

stroma_tf_corrplot <- ggplot(stroma_tf_corr, aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAcyt ~ MCPCounter",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0),limits = clust$labels[clust$order]) +
  scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
  ggpubr::rotate_x_text(angle = 90)

ggsave(file="02_Output/Figures/stroma_tf_corrplot_TEs.svg", plot=stroma_tf_corrplot, width=10, height=6)

#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: all the correlation plots of above together
#-------------------------------------------------------------------------------
sampleweight_corrplots <- ggarrange(plotlist = list(clinical_technical_corrplot, deconvolution_corrplot, stroma_tf_corrplot), ncol = 3, common.legend = T, align = "h", widths = c(0.165, 0.4, 0.5), legend = "right")
ggsave(file="02_Output/Figures/sampleweight_corrplots_TEs.svg", plot=sampleweight_corrplots, width=20, height=6)

#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: IC.6 Gene weights bottom and top genes
#-------------------------------------------------------------------------------
# Get best genes
n_genes = 10
S_mat %>% arrange(desc(IC.6)) -> S_sort
S_sort %>% tail(10) -> IC6_mneg
S_sort %>% head(10) -> IC6_mpos

# Density plot
S_mat %>% ggplot() +
  aes_string(y = "IC.6") +
  geom_density(alpha=.5, fill="lightgrey") +
  scale_x_reverse() +
  scale_y_continuous(expand = c(0,0))+
  geom_hline(yintercept = 0, colour = "black") +
  annotate("rect",xmin = -Inf, xmax = Inf,   ymin = min(IC6_mneg["IC.6"]), ymax = max(IC6_mneg["IC.6"]),   fill = "blue", alpha = 0.5) +
  annotate("rect",xmin = -Inf, xmax = Inf,   ymin =  min(IC6_mpos["IC.6"]), ymax = max(IC6_mpos["IC.6"]),   fill = "red", alpha = 0.5) +
  theme_classic() -> densplot

# Heatmap
## Annotation
annot_row_ <- tibble(name = rownames(IC6_mpos), class = "most possitive")
tibble(name = rownames(IC6_mneg), class = "most negative") %>% bind_rows(annot_row_) %>% column_to_rownames("name") -> annot_row__
annot_row <- annot_row__
rownames(annot_row) <- rownames(annot_row__) %>% translate[.] %>% make.names(unique = TRUE)

annot_col <- complete_annotation[c("PAMG", "ISRact", "IC.6")]

annot_colors <- list(class = c(`most possitive` = "red", `most negative` = "blue"),
                     PAMG = c("#FF7F00", "white", "#377DB8"),
                     ISRact = c("#FFFFCC", "#006837"),
                     IC.6 = c("#FFFFCC", "#5b0066")) 

## Get gene names
ensembl_toplot_ <- TEs %>% dplyr::filter(EnsemblID %in% c(rownames(IC6_mpos), rownames(IC6_mneg))) %>%
  arrange(match(EnsemblID, c(rownames(IC6_mpos), rownames(IC6_mneg))))

genes_toplot <- ensembl_toplot_
genes_toplot$Genenames <- ensembl_toplot_$EnsemblID %>% translate[.] %>%
  make.names(unique = TRUE)

## Plot
heatmap <- genes_toplot %>% dplyr::select(!EnsemblID) %>%
  relocate(Genenames) %>% column_to_rownames("Genenames") %>% 
  pheatmap(scale = "row", color = colorRampPalette(c("#0C6291", "#FBFEF9", "#A63446"))(100),
           annotation_row = annot_row, annotation_col = annot_col, annotation_colors = annot_colors,
           cluster_rows = F, show_colnames = TRUE) %>%
  as.grob()

# Export joint figure
density_heatmap_IC6 <- ggarrange(densplot + rremove("xylab"), heatmap, heights = c(1.5, 10), widths = c(0.2,1),
                                 ncol = 2, nrow = 1)

ggsave(file="02_Output/Figures/density_heatmap_IC6_TEs.svg", plot=density_heatmap_IC6, width=10, height=6)

#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: IC.4 GSVA bottom and top genes
#-------------------------------------------------------------------------------
gsvaRes[order(gsvaRes[,"IC.6"]),]

gsvaTop <- as_tibble(gsvaRes, rownames = "gene_set") %>% 
  mutate(the_rank = rank(-IC.6, ties.method = "random"),
         #gene_set = if_else(str_count(gene_set, "_") < 10, gene_set, signature_dict[gene_set]),
         gene_set = str_remove(gene_set, "REACTOME_"),
         gene_set = fct_reorder(gene_set, the_rank,.desc = T)) %>%
  pivot_longer(cols = -c(gene_set, the_rank), names_to = "component", values_to = "ES") %>% 
  dplyr::filter(the_rank < 25 | the_rank > (nrow(gsvaRes)-25)) %>% 
  mutate(component = if_else(component == "IC.6", "IC.6", "Other")) %>% 
  dplyr::select(!c(the_rank))

gsva_IC6 <- ggplot(gsvaTop, aes(x = ES, y = gene_set)) + 
  geom_point(aes(alpha = if_else(component == "IC.6", 0.9, 0.3),
                 color = if_else(ES > 0, "blue", "red"))) +
  theme_bw() +
  labs(title = "IC.6", subtitle = "Best Reactome gene sets") +
  rremove("legend") +
  rremove("ylab")

ggsave(file="02_Output/Figures/gsva_IC6_TEs.svg", plot=gsva_IC6, width=10, height=6)

#-------------------------------------------------------------------------------
# THESIS PLOTS: Clinical and technical data + other components
#-------------------------------------------------------------------------------
columns_multiqc = c("CITID", "Dups_R1", "GC_R1")
multiqc_sp <- dplyr::inner_join(dplyr::select(multiqc_cyt, all_of(columns_multiqc)),
                                dplyr::select(multiqc_pol, all_of(columns_multiqc)),
                                by="CITID", suffix = c("_cyt", "_pol"))

columns_seq = c("Number_of_reads_mapped_to_too_many_loci", "Per_of_reads_mapped_to_multiple_loci", 
                "Deletion_rate_per_base", "input2maped_length_diff", "Uniquely_mapped_reads_per")
seq_info_sp <- dplyr::inner_join(seq_info_cyt,
                                seq_info_pol,
                                by="CITID", suffix = c("_cyt", "_pol"))


levels = c(paste0("ICA",1:6), "PAMG",
           "RNAconc_cyt", "Dups_R1_cyt", "GC_R1_cyt",  # general and rRNA abundance cyt
           "Number_of_reads_mapped_to_too_many_loci_cyt", "Per_of_reads_mapped_to_multiple_loci_cyt", # rRNA abundance
           "Deletion_rate_per_base_cyt", "input2maped_length_diff_cyt", "Uniquely_mapped_reads_per_cyt",
           "RNAconc_pol", "Dups_R1_pol", "GC_R1_pol", # same but pol
           "Number_of_reads_mapped_to_too_many_loci_pol", "Per_of_reads_mapped_to_multiple_loci_pol",
           "Deletion_rate_per_base_pol", "input2maped_length_diff_pol", "Uniquely_mapped_reads_per_pol") 

clinical_technical_full_corrplot <- sample_info %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(rownames_to_column(A_mat, "CITID")) %>% 
  inner_join(multiqc_sp) %>%
  inner_join(seq_info_sp) %>%
  dplyr::select(CITID, matches("IC."), all_of(levels)) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% levels,
         measure2 %in% names(A_mat)) %>% #not square corr
  dplyr::mutate(measure1 = factor(measure1,
                                  levels = levels)) %>% 
  mutate(r = abs(r)) %>% 
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Translatome ICA vs clinical and technical data",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="black", limits=c(0,1)) +
  geom_text(colour = "white") +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
  ggpubr::rotate_x_text(angle = 90)

ggsave(file="02_Output/Figures/clinical_technical_full_corrplot_TEs.svg", plot=clinical_technical_full_corrplot, width=12, height=6)
sed -i "s/ textLength='[^']*'//" file.svg
#-------------------------------------------------------------------------------
# THESIS PLOTS: Deconvolution plots
#-------------------------------------------------------------------------------
for (method in c("mMCPcounter", "ImmuCC")){
  
  if (method == "mMCPcounter"){
    deconv = mMCPcounter_res
  } else  if (method == "ImmuCC"){
    deconv = ImmuCC_res
  }
  
  deconvolution_corr <- A_mat %>%
    as_tibble(rownames = "CITID") %>%
    inner_join(as_tibble(deconv, rownames = "CITID")) %>%
    dplyr::select(CITID, names(A_mat), names(deconv)) %>%
    column_to_rownames("CITID") %>%
    formatted_cors("spearman") %>%
    filter(measure1 %in% names(deconv),
           measure2 %in% names(A_mat))
  
  clust <- dplyr::select(deconvolution_corr, measure1, measure2, r) %>% 
    pivot_wider(names_from=measure2, values_from = r) %>%
    column_to_rownames("measure1") %>%
    dist() %>%
    hclust()
  
  deconvolution_corrplot <- ggplot(deconvolution_corr, aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title=paste0("Correlations ICAcyt ~ ", method),
         subtitle="Only significant correlation coefficients shown (95% I.C.)") +
    scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
    geom_text() +
    theme_classic() +
    scale_x_discrete(expand=c(0,0),limits = clust$labels[clust$order]) +
    scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
    ggpubr::rotate_x_text(angle = 90)
  
  if (method == "mMCPcounter"){
    ggsave(file=paste0("02_Output/Figures/", method, "_deconv_TEs.svg"), plot=deconvolution_corrplot, width=10, height=6)
    
  } else  if (method == "ImmuCC"){
    ggsave(file=paste0("02_Output/Figures/", method, "_deconv_TEs.svg"), plot=deconvolution_corrplot, width=7, height=6)
  }
  
}


stroma_mostvar_TF <- Get_mostvar(stroma_tf, 15) %>% rownames_to_column() %>% pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)

stroma_tf_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(as_tibble(t(stroma_tf), rownames = "CITID")) %>%
  dplyr::select(CITID, names(A_mat), names(stroma_mostvar_TF)[-1]) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(stroma_mostvar_TF)[-1],
         measure2 %in% names(A_mat))

clust <- dplyr::select(stroma_tf_corr, measure1, measure2, r) %>% 
  pivot_wider(names_from=measure2, values_from = r) %>%
  column_to_rownames("measure1") %>%
  dist() %>%
  hclust()

stroma_tf_small_corrplot <- ggplot(stroma_tf_corr, aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAcyt ~ Stroma most variable TF activities",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0),limits = clust$labels[clust$order]) +
  scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
  ggpubr::rotate_x_text(angle = 90)

ggsave(file="02_Output/Figures/stroma_tf_small_corrplot_TEs.svg", plot=stroma_tf_small_corrplot, width=10, height=6)


tumour_mostvar_TF <- Get_mostvar(tumour_tf, 15) %>% rownames_to_column() %>% pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)

tumour_tf_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(as_tibble(t(tumour_tf), rownames = "CITID")) %>%
  dplyr::select(CITID, names(A_mat), names(tumour_mostvar_TF)[-1]) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(tumour_mostvar_TF)[-1],
         measure2 %in% names(A_mat))

clust <- dplyr::select(tumour_tf_corr, measure1, measure2, r) %>% 
  pivot_wider(names_from=measure2, values_from = r) %>%
  column_to_rownames("measure1") %>%
  dist() %>%
  hclust()

tumour_tf_small_corrplot <- ggplot(tumour_tf_corr, aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAcyt ~ Tumour most variable TF Activities",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0),limits = clust$labels[clust$order]) +
  scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
  ggpubr::rotate_x_text(angle = 90)

ggsave(file="02_Output/Figures/tumour_tf_small_corrplot_TEs.svg", plot=tumour_tf_small_corrplot, width=10, height=6)


#-------------------------------------------------------------------------------
# THESIS PLOTS: Most correlated TFs
#-------------------------------------------------------------------------------
comp = "IC.6"
nTFs = 25
for (fraction in c("tumour", "stroma")){
  if (fraction == "tumour"){
    tf_activity <- tumour_tf
  } else if (fraction == "stroma"){
    tf_activity <- stroma_tf
  }
  
correlations <- rcorr(t(tf_activity), complete_annotation[,comp,drop=T])
comp_p.adj <- p.adjust(correlations$P[,"y"], "BH")
best_tfs <-  correlations$r[,"y"] %>%
  abs() %>% sort(decreasing = T) %>%
  .[2:(nTFs+1)] %>% names()

annotation_TFs <- tibble(TFs = best_tfs, R = correlations$r[best_tfs,"y"], 
                         p.value = correlations$P[best_tfs,"y"],
                         p.adj = as.numeric(comp_p.adj[best_tfs] < 0.05)) %>% 
  column_to_rownames("TFs")

tf_activity %>% .[best_tfs,] %>%
  pheatmap(main=paste(fraction," TFs most correlated with ", comp),
           scale = "row",
           annotation_col = complete_annotation,
           annotation_row = annotation_TFs, 
           filename = paste0("02_Output/Figures/TFact_",fraction,"_IC.6.pdf"),
           width = 10,
           height = 8,
           annotation_colors = list(R = c("red", "white", "blue"),
                                    p.value = c("black", "white"),
                                    p.adj = c("white", "black"),
                                    PAMG = c("#FF7F00", "white", "#377DB8"),
                                    ISRact = c("#FFFFCC", "#006837"),
                                    IC.1 = c("white", "black"),
                                    IC.2 = c("white", "black"),
                                    IC.3 = c("white", "black"),
                                    IC.4 = c("white", "black"),
                                    IC.5 = c("white", "black"),
                                    IC.6 = c("white", "black")))
}
#-------------------------------------------------------------------------------
sed -i "s/ textLength='[^']*'//" file.svg
#-------------------------------------------------------------------------------
