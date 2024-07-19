# DGEA on Shin et al. PDX with respect to ISRactPCA

################################################################################
#' Import rawcount
#' follow the limma voom workflow for DGEA analysis
#' GSEA on ranking regarding DGEA results
################################################################################

library(readr)
library(readxl)
library(edgeR)
library(biomaRt)
library(limma)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(EnhancedVolcano)
library(FactoMineR)
library(factoextra)
library(fgsea)
library(msigdbr)
library(RColorBrewer)
library(dorothea)
library(pheatmap)
library(viridis)

setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")

# Import data
rawTumor_Cyt <- read_delim("data/Sauyeun_PDX/rawTumor_Cyt.tsv",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

## load annotation with Biomart
ensembl75 <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = 75
) # listAttributes(ensembl75, page="feature_page")

annot_ensembl75 <- getBM(attributes = c(
  "ensembl_gene_id",
  "external_gene_id"
), mart = ensembl75)

## Add a Gene names column to rawTumor_Cyt
translate <- deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

rawTumor_Cyt$EnsemblID %>%
  translate[.] %>%
  make.names(unique = TRUE) -> rawTumor_Cyt$Gene

counts <- dplyr::select(rawTumor_Cyt, -EnsemblID) %>%
  column_to_rownames("Gene")

# Import Projection
pdx_PC1 <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS")
Sauyeun_PC1 <- read_tsv("results/Sauyeun_PDX/Shin_ISRactPCA.tsv")

# DGEA following the Limma Voom pipeline from https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
################################################################################
filter_intermediate_beforeDGEA <- F # T/F
################################################################################
## initial sample filtering
if (filter_intermediate_beforeDGEA) {
  sample_info <- dplyr::filter(Sauyeun_PC1, ISRact_bin != "intermediate_ISRact")
} else {
  sample_info <- Sauyeun_PC1
}

## Processing and object building
### Visualize expected counts by taking 10 random samples
random10 <- sort(sample.int(dim(sample_info)[1], 10))

data <- rownames_to_column(counts[, random10], "Prot") %>%
  pivot_longer(cols = 2:11, names_to = "sample", values_to = "count") # pivot_longer data to make them easily plotable with ggplot 2

ggplot(data, mapping = aes(x = sample, y = count, fill = sample)) +
  geom_violin() +
  theme(legend.position = "none") +
  coord_flip() +
  yscale("log10")

### Create DGEList object
d0 <- dplyr::select(counts, sample_info$sample) %>%
  DGEList()

### Calculate normalization factors
d0 <- calcNormFactors(d0, method = "TMM")
d0

### Filter low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop, ]
dim(d) # number of genes left

### Create a new variable “group”
group <- fct(sample_info$ISRact_bin, levels = c("low_ISRact", "intermediate_ISRact", "high_ISRact"))

### Specify the model to be fitted
mm <- model.matrix(~ 0 + group)

### Filter by gene expression
keep <- filterByExpr(d0, mm)
d <- d0[keep, ]
dim(d) # number of genes left

## Visualize before DGEA
dnorm <- cpm(d)

data2 <- as.data.frame(dnorm[, random10]) %>%
  rownames_to_column("Prot") %>%
  pivot_longer(cols = 2:11, names_to = "sample", values_to = "count_norm")

### Custom ggplot counts
ggplot(data2, mapping = aes(x = sample, y = count_norm, fill = sample)) +
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

#' Plot PCA components in relation to a given factor
#' @description
#' this function does a parplot of the desired n of components and color
#' them according to a given factor.
#'
#' prints a list of spearman correlations of the desired metadata factor
#' with the IC that best correlate with it, from the ICAs with the different number of components
#'
#' @param pca_toplot data frame containing the PCs and the factors
#' you want to use to correlate
#' @param feat name of the feature to color
#' @param ncomp number of PCs to plot
#' @param dotsize to adjust the dotsize manually
#'
#' TBD
#' continuous coloring
#'
plot_PCs <- function(pca_toplot, feat, ncomp, dotsize) {
  col_factor <- as.factor(pca_toplot[[feat]])
  col_n <- nlevels(col_factor)
  cols <- brewer.pal(col_n, "Spectral")
  cols <- colorRampPalette(cols)(col_n)
  pairs(pca_toplot[, paste("Dim", 1:ncomp, sep = ".")],
    col = cols[as.numeric(col_factor)],
    pch = 19,
    cex = dotsize,
    lower.panel = NULL,
    main = feat
  )
  par(xpd = TRUE)
  # x11()
  plot.new()
  legend(x = "center", fill = cols, legend = levels(col_factor), horiz = F, title = feat)
}

pca_fulldf <- dnorm.pca$ind$coord %>%
  as_tibble(rownames = "sample") %>%
  inner_join(sample_info[, c("sample", "Diabetes", "PAMG", "ISRact")])

for (feat in colnames(pca_fulldf)[-c(1:(dim(dnorm.pca$ind$coord)[2] + 1))]) {
  plot_PCs(pca_fulldf, feat, dim(dnorm.pca$ind$coord)[2], 1)
}

### Voom transformation and calculation of variance weights
y <- voom(d, mm, plot = T)

## Fitting linear models in limma
fit <- lmFit(y, mm) # lmFit fits a linear model using weighted least squares for each gene
head(coef(fit))

### define contrast matrix
contr <- makeContrasts(grouphigh_ISRact - grouplow_ISRact, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr) # Estimate contrast for each gene

tmp <- eBayes(tmp) # Empirical Bayes smoothing of standard errors

## Visualize top differentiated genes with vulcano plot
top_diff_gene <- topTreat(tmp, coef = 1, n = Inf) %>%
  rownames_to_column("external_gene_name")
EnhancedVolcano(top_diff_gene,
  lab = top_diff_gene$external_gene_name, x = "logFC", y = "P.Value",
  pCutoffCol = "adj.P.Val", pCutoff = 0.05, FCcutoff = 2, title = "Low ISRact vs. high ISRact",
  subtitle = "adj.P.Val < 0.05"
)



### other visualisations of teh viggnettes
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
vennDiagram(dt[, c(1, 3)], circle.col = c("turquoise", "salmon"))

# GSEA
## load gene sets
msigdbr_reactome_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  mutate(gs_name = str_remove(gs_name, "REACTOME_"))
msigdbr_hallmark_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  mutate(gs_name = str_remove(gs_name, "HALLMARK_"))


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

head(gsea_res_reactome[order(pval), ])
head(gsea_res_hallmark[order(pval), ])

## Visualisation
### fGSEA plots
# Function to plot the top 10 low ES and high ES pathways regarding to the pvalue of the gsea
#' @description
#' this function plots table of enrichment graphs for the 10
#' pathways with positive enrichment score and lowest pvalue
#' and the 10 pathways with negative enrichment score and lowest pvalue
#'
#' @param gsea.res data frame of a gsea result returned by fgsea function
#' @param msigdbr_list list of genes name for each pathway from msigdb
#' @param ranked_genes named vector where names are genes name and values are pca values
plotBestPathways <- function(gsea.res, msigdbr_list, ranked_genes) {
  topPathwaysUp <- gsea.res[ES > 0][head(order(pval), n = 10), pathway]
  topPathwaysDown <- gsea.res[ES < 0][head(order(pval), n = 10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plot.new()
  plotGseaTable(msigdbr_list[topPathways], ranked_genes, gsea.res,
    gseaParam = 0.5
  )
}

plotBestPathways(gsea_res_reactome, msigdbr_list_reactome, ranked_genes)
plotBestPathways(gsea_res_hallmark, msigdbr_list_hallmark, ranked_genes)

### customn plots thesis
#### Reactome
topPathwaysUp <- gsea_res_reactome[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown <- gsea_res_reactome[ES < 0][head(order(pval), n = 10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

shinbarplot_DEGEA_reactome <- dplyr::filter(gsea_res_reactome, pathway %in% topPathways) %>%
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

print(shinbarplot_DEGEA_reactome)
ggsave(shinbarplot_DEGEA_reactome,
  filename = "results/Figures/shinbarplot_DEGEA_reactome.svg",
  width = 8,
  height = 5
)

#### Hallmark
topPathwaysUp <- gsea_res_hallmark[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown <- gsea_res_hallmark[ES < 0][head(order(pval), n = 10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

shinbarplot_DEGEA_hallmark <- dplyr::filter(gsea_res_hallmark, pathway %in% topPathways) %>%
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
    colours = c("red", "white", "grey"), guide = guide_colourbar(reverse = TRUE),
    rescaler = ~ scales::rescale_mid(., mid = 0.05)
  ) +
  theme_bw() +
  labs(title = "HALLMARK") +
  ylab("")

print(shinbarplot_DEGEA_hallmark)
ggsave(shinbarplot_DEGEA_hallmark,
  filename = "results/Figures/shinbarplot_DEGEA_hallmark.svg",
  width = 8,
  height = 5
)


# most associated TF activities
## Generate TF act
### Select the regulon data
data(dorothea_hs_pancancer, package = "dorothea")
regulons <- dorothea_hs_pancancer %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))

### VIPER Analysis
minsize <- 5
ges.filter <- FALSE

tf_act__ <- dorothea::run_viper(counts, regulons,
  options = list(
    minsize = minsize, eset.filter = ges.filter,
    cores = 1, verbose = FALSE, nes = TRUE
  )
)

tf_act_ <- as_tibble(tf_act__, rownames = "TF") %>%
  pivot_longer(cols = -TF, names_to = "bcr_patient_barcode") %>%
  pivot_wider(id_cols = c(bcr_patient_barcode), names_from = TF)

### Visualisation
#### Most variable TFs
##### Calculate the n most variable TFs
n <- 50
tf_50_var <- summarise(tf_act_, across(where(is.double), var)) %>%
  pivot_longer(cols = everything(), names_to = "TF", values_to = "var") %>%
  dplyr::arrange(desc(var)) %>%
  dplyr::slice_head(n = n)

###### Define Annotations
annotation_row_var <- tf_50_var %>%
  column_to_rownames(var = "TF")

annotation_col <- dplyr::select(sample_info, sample, Diabetes, PAMG, ISRact, ISRact_bin, ISRactPCA) %>%
  column_to_rownames(var = "sample")

ann_colors <- list(
  ISRactPCA = c("seagreen", "white", "tomato3"),
  ISRact = c("#FFFFCC", "#006837"),
  ISRact_bin = c(`high_ISRact` = "brown", `intermediate_ISRact` = "grey", `low_ISRact` = "#006837"),
  PAMG = c("#FF7F00", "white", "#0000FC"),
  Diabetes = c(`no` = "grey", `yes` = "black", `NA` = "white"),
  var = c("white", "darkgreen"),
  cor = brewer.pal(n = 11, name = "BrBG")
)

# Create breaks and limits for the continuous expression scale
breaksList <- seq(-3.5, 3.5, by = 0.1)

##### Heatmap
as_tibble(tf_act__, rownames = "TF") %>%
  dplyr::select(TF, sample_info$sample) %>%
  dplyr::filter(TF %in% tf_50_var$TF) %>%
  column_to_rownames("TF") %>%
  as.matrix() %>%
  pheatmap(
    annotation_col = annotation_col,
    annotation_row = annotation_row_var,
    annotation_colors = ann_colors, scale = "row", show_colnames = FALSE,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
    breaks = breaksList,
    fontsize = 12,
    annotation_names_row = FALSE,
    main = paste0("Activity of the ", n, " most variable transcription factor on the Sauyeun PDX dataset")
  )

#### Most correlated with ISRact
##### calculate correlation TF ISRact
tf_ISRact_cor <- tf_act_ %>%
  dplyr::rename(sample = bcr_patient_barcode) %>%
  inner_join(dplyr::select(sample_info, sample, ISRactPCA, ISRact, PAMG), by = "sample") %>% # make the two table in same order
  arrange(sample) %>%
  column_to_rownames("sample") %>% # select continuous variables
  as.matrix() %>%
  cor(method = "spearman") %>%
  as_tibble(rownames = "TF") %>%
  dplyr::select(TF, ISRact)

##### filter by absolute correlation
n <- 30
top_tf_ISRact <- mutate(tf_ISRact_cor, ISRact = abs(ISRact)) %>%
  dplyr::filter(!str_detect(TF, "ISRact")) %>%
  arrange(desc(ISRact)) %>%
  dplyr::slice(1:n) %>%
  dplyr::rename(abscor = ISRact)

##### Define annotations
annotation_row_cor_ISRact <- dplyr::filter(tf_ISRact_cor, TF %in% top_tf_ISRact$TF) %>%
  dplyr::rename(cor = ISRact) %>%
  column_to_rownames(var = "TF")

##### Heatmap
as_tibble(tf_act__, rownames = "TF") %>%
  dplyr::select(TF, sample_info$sample) %>%
  dplyr::filter(TF %in% top_tf_ISRact$TF) %>%
  column_to_rownames("TF") %>%
  as.matrix() %>%
  pheatmap(
    annotation_col = annotation_col,
    annotation_row = annotation_row_cor_ISRact,
    annotation_colors = ann_colors, scale = "row", show_colnames = FALSE,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
    breaks = breaksList,
    fontsize = 12,
    annotation_names_row = FALSE,
    main = paste0("Activity of the ", n, " most correlated transcription factor  with ICA3 on the Shin dataset"),
    height = 6,
    width = 10,
    filename = "results/Figures/shin_ISRactPCA_TFact.pdf"
  )

#### Most correlated with ISRactPCA
##### calculate correlation TF ISRactPCA
tf_ISRactPCA_cor <- tf_act_ %>%
  dplyr::rename(sample = bcr_patient_barcode) %>%
  inner_join(dplyr::select(sample_info, sample, ISRactPCA, ISRact, PAMG), by = "sample") %>% # make the two table in same order
  arrange(sample) %>%
  column_to_rownames("sample") %>% # select continuous variables
  as.matrix() %>%
  cor(method = "spearman") %>%
  as.data.frame() %>%
  dplyr::select("ISRactPCA")

##### filter by absolute correlation
n <- 30
top_tf_ISRactPCA <- mutate(tf_ISRactPCA_cor, ISRactPCA = abs(ISRactPCA)) %>%
  arrange(desc(ISRactPCA)) %>%
  dplyr::slice(1:n) %>%
  dplyr::rename(abscor = ISRactPCA) %>%
  rownames_to_column(var = "TF")

##### Define annotations
annotation_row_cor_ISRactPCA <- rownames_to_column(tf_ISRactPCA_cor, "TF") %>%
  inner_join(dplyr::select(top_tf_ISRactPCA, TF), by = "TF") %>%
  rename(cor = ISRactPCA) %>%
  column_to_rownames(var = "TF")

##### Heatmap
as_tibble(tf_act__, rownames = "TF") %>%
  dplyr::select(TF, sample_info$sample) %>%
  dplyr::filter(TF %in% top_tf_ISRactPCA$TF) %>%
  column_to_rownames("TF") %>%
  as.matrix() %>%
  pheatmap(
    annotation_col = annotation_col,
    annotation_row = annotation_row_cor_ISRactPCA,
    annotation_colors = ann_colors, scale = "row", show_colnames = FALSE,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
    breaks = breaksList,
    fontsize = 12,
    annotation_names_row = FALSE,
    main = "Activity of the 50 most correlated transcription factor with ISRactPCA on the Sauyeun PDX dataset"
  )


# Proportion of common top correlated TF with PC1 and ICA3
length(Reduce(intersect, list(top_tf_ISRact$TF, top_tf_ISRactPCA$TF))) / n
