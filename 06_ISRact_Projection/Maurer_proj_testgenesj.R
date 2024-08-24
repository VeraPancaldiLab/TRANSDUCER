# DGEA on Maurer cohort RNA-seq rawcount #! Need to adapt PC1 -> ISRactPCA

################################################################################
#' Import rawcount
#' follow the limma voom workflow for DGEA analysis
#' GSEA on ranking regarding DGEA results
################################################################################

library(readr)
library(readxl)
library(edgeR)
library(limma)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(EnhancedVolcano)
library(fgsea)
library(msigdbr)
library(RColorBrewer)
library(dorothea)
library(pheatmap)


## setwd and import local functions
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/plot_PCs.R")
source("src/plotBestPathways.R")

################################################################################
subset_celltype <- "epithelium" # epithelium | stroma | bulk
################################################################################
# Import data
## Expression
load("data/Maurer2019/Start Maurer EvsSvsB.RData")

## Metadata
### Import ISRactPCS projection of different subset  #! This comes from 02_TRANSDUCER/06_Enzos_Internship/human_cohort_analysis/Data
if (subset_celltype == "epithelium") {
  Maurer_PC1 <- read_csv("results/Maurer/Maurer_epithelium_PC1.csv")
} else if (subset_celltype == "stroma") {
  Maurer_PC1 <- read_csv("results/Maurer/Maurer_stroma_PC1.csv")
} else if (subset_celltype == "bulk") {
  Maurer_PC1 <- read_csv("results/Maurer/Maurer_bulk_PC1.csv")
}

### load sample information for annotations
sample_info <- rename(Maurer_PC1, Sample = sample) %>%
  inner_join(phenoEvsS, by = "Sample")

# DGEA following the Limma Voom pipeline from https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
################################################################################
## Select samples corresponding to the subset cell type
counts <- dplyr::select(raw_E_Vs_S, Maurer_PC1$sample)

### Visualize raw counts
random10 <- sort(sample.int(dim(sample_info)[1], 10))

data <- rownames_to_column(counts[, random10], "Prot") %>%
  pivot_longer(cols = 2:11, names_to = "case_id", values_to = "Expected_count") # pivot_longer data to make them easily plotable with ggplot 2

ggplot(data, mapping = aes(x = case_id, y = Expected_count, fill = case_id)) +
  geom_violin() +
  theme(legend.position = "none") +
  coord_flip() +
  yscale("log10")

# Follow the Limma Voom pipeline from https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

### Create DGEList object
d0 <- dplyr::select(counts, sample_info$Sample) %>%
  DGEList()

### Calculate normalization factors
d0 <- calcNormFactors(d0)
d0

### Filter low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop, ]
dim(d) # number of genes left

### Create a new variable “group”
group <- as.factor(Maurer_PC1$PC1status)

### Specify the model to be fitted
mm <- model.matrix(~ 0 + group)

### Filter with filter by expression
keep <- filterByExpr(d0, mm)
d <- d0[keep, ]
dim(d) # number of genes left

### Visualize after normalization and filtering
dnorm <- cpm(d)

data2 <- as.data.frame(dnorm[, random10]) %>%
  rownames_to_column("Prot") %>%
  pivot_longer(cols = 2:11, names_to = "case_id", values_to = "count_norm") # pivot_longer data to make them easily plotable with ggplot 2

ggplot(data2, mapping = aes(x = case_id, y = count_norm, fill = case_id)) +
  geom_violin() +
  theme(legend.position = "none") +
  coord_flip() +
  yscale("log10")

### Multidimensional scaling (MDS) plot
plotMDS(d, col = as.numeric(group))

### PCA plot regarding clinical and phenotype data

dnorm.t <- t(dnorm)
dnorm.pca <- PCA(dnorm.t, ncp = 10, graph = FALSE)
summary(dnorm.pca)

fviz_eig(dnorm.pca, addlabels = TRUE, ylim = c(0, 50))
eig.val <- get_eigenvalue(dnorm.pca)
eig.val # 63 dimensions to reach 90% of explained variance

# Join PCA results and more relevant clinical data
pca_fulldf <- dnorm.pca$ind$coord %>%
  as_tibble(rownames = "Sample")

for (feat in colnames(pca_fulldf)[-c(1:11)]) {
  plot_PCs(pca_fulldf, feat, 10, 1)
}

## Voom transformation and calculation of variance weights
y <- voom(d, mm, plot = T)

## Fitting linear models in limma
fit <- lmFit(y, mm)
head(coef(fit))

### define contrast matrix
contr <- makeContrasts(grouphigh_PC1 - grouplow_PC1, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr) # Estimate contrast for each gene

tmp <- eBayes(tmp) # Empirical Bayes smoothing of standard errors

## Visualize top differentiated genes with vulcano plot
top_diff_gene <- topTreat(tmp, coef = 1, n = Inf) %>%
  rownames_to_column("external_gene_name")
EnhancedVolcano(top_diff_gene,
  lab = top_diff_gene$external_gene_name, x = "logFC", y = "P.Value",
  pCutoffCol = "adj.P.Val", pCutoff = 0.05, title = "Low PC1 vs. high PC1",
  subtitle = "adj.P.Val < 0.05"
)


### Other visualisations of the vignettes
top.table <- topTable(tmp, sort.by = "P", n = Inf) %>%
  mutate(
    Expression = case_when(
      logFC >= log(2) & adj.P.Val <= 0.05 ~ "Up-regulated",
      logFC <= -log(2) & adj.P.Val <= 0.05 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )
  ) %>%
  rownames_to_column("Gene")

head(top.table, 20)

# How many DE genes are there?
length(which(top.table$adj.P.Val < 0.05))

dt <- decideTests(fit)
summary(dt)
vennDiagram(dt[, 1:2], circle.col = c("turquoise", "salmon"))

# GSEA on l2FC
## load gene sets
msigdbr_reactome_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
msigdbr_hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")

## extract ranked list
ranked_genes <- top_diff_gene %>%
  as_tibble() %>%
  arrange(desc(logFC))
ranked_genes <- ranked_genes[, c(1, 2)] %>% deframe()

msigdbr_list_reactome <- split(x = msigdbr_reactome_df$gene_symbol, f = msigdbr_reactome_df$gs_name)
msigdbr_list_hallmark <- split(x = msigdbr_hallmark_df$gene_symbol, f = msigdbr_hallmark_df$gs_name)

## analysis
gsea_res_reactome <- fgsea(pathways = msigdbr_list_reactome, stats = ranked_genes, minSize = 15)
gsea_res_hallmark <- fgsea(pathways = msigdbr_list_hallmark, stats = ranked_genes, minSize = 15)


## Plots
### Official viggnette plots
plotBestPathways(gsea_res_reactome, msigdbr_list_reactome, ranked_genes)
plotBestPathways(gsea_res_hallmark, msigdbr_list_hallmark, ranked_genes)

### custom plots thesis
#### Reactome
topPathwaysUp <- gsea_res_reactome[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown <- gsea_res_reactome[ES < 0][head(order(pval), n = 10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

maurerbarplot_DEGEA_reactome <- dplyr::filter(gsea_res_reactome, pathway %in% topPathways) %>%
  arrange(NES) %>%
  mutate(
    pathway = str_remove(pathway, "REACTOME_"),
    pathway = as_factor(pathway)
  ) %>%
  # dplyr::filter(p.adjust < p.adjust_th) %>%
  # slice_tail(n=20) %>%
  ggplot(aes(x = NES, y = pathway, fill = padj)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(
    colours = c("red", "white", "grey"), guide = guide_colourbar(reverse = TRUE),
    rescaler = ~ scales::rescale_mid(., mid = 0.05)
  ) +
  theme_bw() +
  labs(title = "REACTOME") +
  ylab("")

print(maurerbarplot_DEGEA_reactome)
# ggsave(maurerbarplot_DEGEA_reactome,
#        filename = paste0("results/Figures/maurer",subset_celltype ,"barplot_DEGEA_reactome.svg"),
#        width = 8,
#        height = 5)

#### Hallmark
topPathwaysUp <- gsea_res_hallmark[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown <- gsea_res_hallmark[ES < 0][head(order(pval), n = 10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

maurerbarplot_DEGEA_hallmark <- dplyr::filter(gsea_res_hallmark, pathway %in% topPathways) %>%
  arrange(NES) %>%
  mutate(
    pathway = as_factor(pathway),
    padj = ifelse(padj < 0.1, padj, NA)
  ) %>% # else the colour bar is screwed
  # dplyr::filter(p.adjust < p.adjust_th) %>%
  # slice_tail(n=20) %>%
  ggplot(aes(x = NES, y = pathway, fill = padj)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(
    colours = c("red", "white"), guide = guide_colourbar(reverse = TRUE),
    rescaler = ~ scales::rescale_mid(., mid = 0.05)
  ) +
  theme_bw() +
  labs(title = "HALLMARK") +
  ylab("")

print(maurerbarplot_DEGEA_hallmark)
# ggsave(maurerbarplot_DEGEA_hallmark,
#        filename = paste0("results/Figures/maurer",subset_celltype ,"barplot_DEGEA_hallmark.svg"),
#        width = 8,
#        height = 5)


# most correlated TF activities
## Generate TF act
### Select the regulon data
data(dorothea_hs_pancancer, package = "dorothea")
regulons <- dorothea_hs_pancancer %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))

### VIPER Analysis
################################################################################
minsize <- 5
ges.filter <- FALSE
################################################################################
tf_act__ <- dorothea::run_viper(counts, regulons,
  options = list(
    minsize = minsize, eset.filter = ges.filter,
    cores = 1, verbose = FALSE, nes = TRUE
  )
)

# Format for analysis
tf_act_ <- as_tibble(tf_act__, rownames = "TF") %>%
  pivot_longer(cols = -TF, names_to = "bcr_patient_barcode") %>%
  pivot_wider(id_cols = c(bcr_patient_barcode), names_from = TF)

### Visualisation
#### Most variable TFs
##### Calculate the n most variable TFs
tf_50_var <- summarise(tf_act_, across(where(is.double), var)) %>%
  pivot_longer(cols = everything(), names_to = "TF", values_to = "var") %>%
  dplyr::arrange(desc(var)) %>%
  dplyr::slice_head(n = 50)

# Visualize TF activity
# Create the matrix with only 50 top variable TF
tf_matrix_var <- as_tibble(tf_act__, rownames = "TF") %>%
  dplyr::select(TF, Maurer_PC1$sample) %>%
  dplyr::filter(TF %in% tf_50_var$TF) %>%
  column_to_rownames("TF") %>%
  as.matrix()

# Create the matrix with only 50 most correlated TF with PC1


# Calculate absolute correlation for each gene
tf_PC1_cor <- tf_act_ %>%
  rename(sample = bcr_patient_barcode) %>%
  inner_join(dplyr::select(Maurer_PC1, sample, PC1), by = "sample") %>% # make the two table in same order
  arrange(PC1) %>%
  dplyr::select(-PC1) %>%
  column_to_rownames("sample") %>% # select continuous variables
  as.matrix() %>%
  cor(y = Maurer_PC1$PC1, method = "spearman") %>%
  as.data.frame()

# keep the 50 most absolute correlated TF
top_tf <- mutate(tf_PC1_cor, V1 = abs(V1)) %>%
  arrange(desc(V1)) %>%
  slice(1:50) %>%
  rename(cor = V1) %>%
  rownames_to_column(var = "TF")

# Final matrix
tf_matrix_cor <- as_tibble(tf_act__, rownames = "TF") %>%
  dplyr::select(TF, Maurer_PC1$sample) %>%
  dplyr::filter(TF %in% top_tf$TF) %>%
  column_to_rownames("TF") %>%
  as.matrix()


# Get a df of the variables I want to annotate the columns of the heatmap
annotation_col <- dplyr::select(sample_info, Sample, PC1) %>%
  column_to_rownames(var = "Sample")

# Get a df of the variables I want to annotate the rows of the heatmap
annotation_row_var <- tf_50_var %>%
  column_to_rownames(var = "TF")

# Choose the color scale for each annotation variable
ann_colors <- list(
  PC1 = c("white", "firebrick"),
  PC1_weight = c("white", "darkorange"),
  KRAS_VAF = c("white", "#0B2356"),
  Moffitt = c(`BASAL-LIKE` = "#FF7F00", CLASSICAL = "#0000FC", `NA` = "grey"),
  Bailey = c(ADEX = "#713E23", immunogenic = "#EF2C18", `NA` = "darkgrey", `pancreatic progenitor` = "#0000FC", squamous = "#FF7F00"),
  Diabetes = c(no = "grey", yes = "black"),
  var = c("white", "darkgreen"),
  cor = brewer.pal(n = 11, name = "BrBG")
)

# Create breaks and limits for the continuous expression scale
breaksList <- seq(-3.5, 3.5, by = 0.1)

# Plot the final heatmap and pray so it works
# Top variable tf
pheatmap(tf_matrix_var,
  annotation_col = annotation_col,
  annotation_row = annotation_row_var,
  annotation_colors = ann_colors, scale = "row", show_colnames = FALSE,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
  breaks = breaksList,
  fontsize = 12,
  annotation_names_row = FALSE,
  main = paste("Activity of the 50 most variable transcription factor on the Maurer", subset_celltype, "dataset", sep = " ")
)

# Get a df of the variables I want to annotate the rows of the heatmap
annotation_row_cor <- tf_PC1_cor %>%
  rownames_to_column("TF") %>%
  inner_join(dplyr::select(top_tf, TF), by = "TF") %>%
  rename(cor = V1) %>%
  column_to_rownames(var = "TF")

# Top correlated tf with PC1
pheatmap(tf_matrix_cor,
  annotation_col = annotation_col,
  annotation_row = annotation_row_cor,
  annotation_colors = ann_colors, scale = "row", show_colnames = FALSE,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
  breaks = breaksList,
  fontsize = 12,
  annotation_names_row = FALSE,
  main = paste("Activity of the 50 most correlated transcription factor with PC1 on the Maurer", subset_celltype, "dataset", sep = " ")
)
