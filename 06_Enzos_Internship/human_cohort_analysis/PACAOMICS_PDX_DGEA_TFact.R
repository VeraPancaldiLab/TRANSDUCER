#DGEA on Remy Nicolle PDX RNA-seq rawcount (THIS IS NOT PACAOMICS THIS IS SAUYEUNS PDXs RAW COUNTS FROM REMY NICOLLE)

################################################################################
#'Import rawcount
#'follow the limma voom workflow for DGEA analysis
#'GSEA on ranking regarding DGEA results
################################################################################

library(readr)
library(readxl)
library(edgeR)
library(biomaRt)
library(Glimma)
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

setwd("~/Documents/02_TRANSDUCER/06_Enzos_Internship/human_cohort_analysis/Data/")

# Import data
## expression
load("Alexias_53PDX/PDX_HUMAN_RAW.RData")
PACAOMICs_90_raw <- as_tibble(x, rownames ="EnsemblID")

RN2017_raw <- read_delim("Human-Tumor_rawcount_Transcriptome.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

## Import PC weight 
PACAOMICs_90_PC1 <- read_csv("../02_Output/PACAOMICS_PC1_extended90.csv")

## load sample inf
full_sample_info <- read_delim("sample_info.tsv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE) %>%
  dplyr::inner_join(PACAOMICS_PC1, by = "sample") %>%
  mutate(Diabetes = replace(Diabetes, Diabetes == 1, "yes"))%>%
  mutate(Diabetes = replace(Diabetes, Diabetes == 0, "no"))%>%
  dplyr::rename(case_id = sample) %>%
  arrange(case_id)

#Select extreme ICA3 samples
top_sample_info <- full_sample_info %>%
  arrange(ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact = ifelse(ICA3 < 0, "low_ICA3", "high_ICA3")) %>%
  arrange(case_id)

################################################################################
data <- "extended90" # "extended90" | "RemyNicolle27"
if (data == "extended90") {
  raw_data = PACAOMICs_90_raw
  sample_info = PACAOMICs_90_PC1
} else if (data == "RemyNicolle27") {
  raw_data = RN2017_raw
  sample_info = NULL # change parameters in PCA_and_projection.R and produce it
}
################################################################################
# Translate EnsemblID to gene names
## load annotation with Biomart Version 75 for PDX data
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)#listAttributes(ensembl75, page="feature_page")

annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = ensembl75)

## Add a Gene names column to rawTumor_Cyt
translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

raw_data <- dplyr::mutate(raw_data, Gene = make.names(translate[EnsemblID], unique = TRUE)) %>% dplyr::select(-EnsemblID) %>%
  column_to_rownames("Gene")

#Visualize expected counts 
random10 <- sort(sample.int(dim(raw_data)[2], 10))

rownames_to_column(raw_data[, random10], "Prot") %>%
  pivot_longer(cols = 2:11, names_to = "case_id", values_to = "count")  %>% 
  ggplot(data, mapping = aes(x = case_id, y = count, fill = case_id)) +
  geom_violin() +
  theme(legend.position = "none") +
  coord_flip() +
  yscale("log10")

#Follow the Limma Voom pipeline from https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html 

#Create DGEList object
d0 <- DGEList(raw_data) # Do a select here to filter and exclude some samples
dim(d0) # original n of genes

#Calculate normalization factors
d0 <- calcNormFactors(d0)
d0

#Create a new variable “group”
group <- as.factor(sample_info$PC1status)

#Specify the model to be fitted
mm <- model.matrix(~0 + group)

#Filter with filter by expression 
keep <- filterByExpr(d0, mm)
d <- d0[keep,]
dim(d) # number of genes left

#Visualize after normalization and filtering
dnorm <- cpm(d)

as.data.frame(dnorm[, random10]) %>%
  rownames_to_column("Prot") %>%
  pivot_longer(cols = 2:11, names_to = "case_id", values_to = "count_norm") %>%  #pivot_longer data to make them easily plotable with ggplot 2
  ggplot(data2, mapping = aes(x = case_id, y = count_norm, fill = case_id))+
  geom_violin()+
  theme(legend.position = "none")+
  coord_flip()+
  yscale("log10")

#Multidimensional scaling (MDS) plot
plotMDS(d, col = as.numeric(group))

#PCA plot regarding clinical and phenotype data

dnorm.t <-  t(dnorm)
dnorm.pca <- PCA(dnorm.t, ncp = 10, graph = FALSE)
summary(dnorm.pca)

fviz_eig(dnorm.pca, addlabels = TRUE, ylim = c(0, 50))
eig.val <- get_eigenvalue(dnorm.pca)
eig.val #63 dimensions to reach 90% of explained variance

#' Plot PCA components in relation to a given factor
#'@description
#' this function does a parplot of the desired n of components and color
#' them according to a given factor.
#'
#'prints a list of spearman correlations of the desired metadata factor
#'with the IC that best correlate with it, from the ICAs with the different number of components
#'
#'@param pca_toplot data frame containing the PCs and the factors
#'you want to use to correlate
#'@param feat name of the feature to color
#'@param ncomp number of PCs to plot
#'@param dotsize to adjust the dotsize manually
#'
#'TBD
#'continuous coloring
#'
plot_PCs <- function(pca_toplot, feat, ncomp, dotsize){
  col_factor <- as.factor(pca_toplot[[feat]])
  col_n <- nlevels(col_factor)
  cols <- brewer.pal(col_n, "Spectral")
  cols <- colorRampPalette(cols)(col_n)
  pairs(pca_toplot[, paste("Dim",1:ncomp,sep =".")],
        col = cols[as.numeric(col_factor)],
        pch = 19,
        cex = dotsize,
        lower.panel = NULL,
        main = feat)
  par(xpd = TRUE)
  #x11()
  plot.new()
  legend(x = "center",fill = cols, legend = levels(col_factor), horiz = F, title = feat)
}

#Join PCA results and more relevant clinical data
pca_fulldf <- dnorm.pca$ind$coord %>%
  as_tibble(rownames="sample") %>%
  inner_join(sample_info[, c("sample", "PAMG", "ICA3")])

#plot pca dimensions with the good number of dimensions
for (feat in colnames(pca_fulldf)[-c(1:(dim(dnorm.pca$ind$coord)[2]+1))]){
  plot_PCs(pca_fulldf,feat,dim(dnorm.pca$ind$coord)[2],1)
  
}

#Voom transformation and calculation of variance weights
y <- voom(d, mm, plot = T)

#Fitting linear models in limma
#lmFit fits a linear model using weighted least squares for each gene
fit <- lmFit(y, mm)
head(coef(fit))

#Comparison between high and low PC1 
#Specify group to compare
contr <- makeContrasts(grouphigh_PC1 - grouplow_PC1, levels = colnames(coef(fit)))
contr

#Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#Empirical Bayes smoothing of standard errors
tmp <- eBayes(tmp)

#Visualize top differentiated genes with vulcano plot
top_diff_gene <- topTreat(tmp, coef=1, n=Inf) %>%
  rownames_to_column("external_gene_name")
EnhancedVolcano(top_diff_gene, lab = top_diff_gene$external_gene_name, x = "logFC", y = "P.Value",
                pCutoffCol = "adj.P.Val", pCutoff = 0.05, title = "Low PC1 vs. high PC1",
                subtitle = "adj.P.Val < 0.05")



#What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf) %>%
  mutate(
    Expression = case_when(logFC >= log(2) & adj.P.Val <= 0.05 ~ "Up-regulated",
                           logFC <= -log(2) & adj.P.Val <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")) %>%
  rownames_to_column("Gene")

head(top.table, 20) %>% View()

#How many DE genes are there?
length(which(top.table$adj.P.Val < 0.05))

dt <- decideTests(fit)
summary(dt)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))

#GSEA with reactome db on DGEA results
msigdbr_reactome_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
msigdbr_hallmark_df = msigdbr(species = "Homo sapiens", category = "H")


ranked_genes <- top_diff_gene %>%  as_tibble() %>% arrange(desc(logFC)) 
ranked_genes <- ranked_genes[,c(1,2)] %>% deframe()

msigdbr_list_reactome = split(x = msigdbr_reactome_df$gene_symbol, f = msigdbr_reactome_df$gs_name)
msigdbr_list_hallmark = split(x = msigdbr_hallmark_df$gene_symbol, f = msigdbr_hallmark_df$gs_name)

gsea_res_reactome <- fgsea(pathways = msigdbr_list_reactome, stats = ranked_genes, minSize = 15)
gsea_res_hallmark <- fgsea(pathways = msigdbr_list_hallmark, stats = ranked_genes, minSize = 15)


head(gsea_res_reactome[order(pval), ])
head(gsea_res_hallmark[order(pval), ])

#Function to plot the top 10 low ES and high ES pathways regarding to the pvalue of the gsea
#'@description
#' this function plots table of enrichment graphs for the 10 
#' pathways with positive enrichment score and lowest pvalue
#' and the 10 pathways with negative enrichment score and lowest pvalue
#'
#'@param gsea.res data frame of a gsea result returned by fgsea function
#'@param msigdbr_list list of genes name for each pathway from msigdb
#'@param ranked_genes named vector where names are genes name and values are pca values
plotBestPathways <- function(gsea.res, msigdbr_list, ranked_genes){
  topPathwaysUp <- gsea.res[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- gsea.res[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plot.new()
  plotGseaTable(msigdbr_list[topPathways], ranked_genes, gsea.res, 
                gseaParam=0.5)
  
  
}

plotBestPathways(gsea_res_reactome, msigdbr_list_reactome, ranked_genes)
plotBestPathways(gsea_res_hallmark, msigdbr_list_hallmark, ranked_genes)

#list of genes for each pathway
#Reactome top 10
gsea_res_reactome[ES > 0][head(order(pval), n=10), leadingEdge]
#Down 10
gsea_res_reactome[ES < 0][head(order(pval), n=10), leadingEdge]
#Hallmark top 10
gsea_res_hallmark[ES > 0][head(order(pval), n=10), leadingEdge]
#Down 10
gsea_res_hallmark[ES < 0][head(order(pval), n=10), leadingEdge]

#TF activity 

# Generate TF act 
### Select the regulon data
data(dorothea_hs_pancancer, package = "dorothea")
regulons <- dorothea_hs_pancancer %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

## VIPER
minsize = 5
ges.filter = FALSE

tf_act__ <- dorothea::run_viper(raw_data, regulons, # Is this good with raw data?
                                options =  list(minsize = minsize, eset.filter = ges.filter, 
                                                cores = 1, verbose = FALSE, nes = TRUE))

# Format for analysis
tf_act_ <- as_tibble(tf_act__, rownames = "TF") %>%
  pivot_longer(cols = -TF, names_to = "sample") %>%
  pivot_wider(id_cols = c(sample), names_from = TF)

# get 50 most variable
tf_50_var <- summarise(tf_act_, across(where(is.double), var)) %>%
  pivot_longer(cols = everything(), names_to = "TF", values_to = "var") %>%
  dplyr::arrange(desc(var)) %>% dplyr::slice_head(n=50)

#Visualize TF activity
#Create the matrix with only 50 top variable TF 
tf_matrix_var <-  as_tibble(tf_act__, rownames = "TF" ) %>%
  dplyr::select(TF, sample_info$sample) %>%
  dplyr::filter(TF %in% tf_50_var$TF) %>%
  column_to_rownames("TF") %>%
  as.matrix()

#Create the matrix with only 50 most correlated TF with ICA3

#Calculate absolute correlation for each gene
tf_ICA3_cor <- tf_act_ %>% 
  inner_join(dplyr::select(sample_info, sample), by = "sample")  %>% #make the two table in same order
  arrange(sample) %>%
  column_to_rownames("sample") %>% #select continuous variables
  as.matrix() %>%
  cor(y = sample_info$ICA3, method = "spearman", use = "complete.obs") %>%
  as.data.frame()

#keep the 50 most absolute correlated TF
top_tf_ICA3 <- mutate(tf_ICA3_cor, V1 = abs(V1)) %>%
  dplyr::arrange(desc(V1)) %>%
  dplyr::slice(1:50) %>%
  dplyr::rename(cor = V1) %>%
  rownames_to_column(var = "TF") 

#Final matrix
tf_matrix_cor_ICA3 <-  as_tibble(tf_act__, rownames = "TF" ) %>%
  dplyr::select(TF, sample_info$sample) %>%
  dplyr::filter(TF %in% top_tf_ICA3$TF) %>%
  column_to_rownames("TF") %>%
  as.matrix()


#Get a df of the variables I want to annotate the columns of the heatmap
annotation_col <- dplyr::select(sample_info, sample, PAMG, ICA3, PC1) %>%
  column_to_rownames(var = "sample")



#Get a df of the variables I want to annotate the rows of the heatmap
annotation_row_var <- tf_50_var %>% 
  column_to_rownames(var = "TF") 

#Choose the color scale for each annotation variable
ann_colors <- list(
  PC1 = c("#00BFC4", "white", "#F8766D"),
  ICA3 = c("white", "darkgoldenrod"),
  PAMG = c("#FF7F00", "white", "#0000FC"),
  Diabetes = c(`no` = "grey", `yes` = "black", `NA` = "white"),
  var = c("white", "darkgreen"),
  cor = brewer.pal(n = 11, name = 'BrBG')
)

#Create breaks and limits for the continuous expression scale
breaksList = seq(-3.5, 3.5, by = 0.1)

#Plot the final heatmap and pray so it works
#Top variable tf
pheatmap(tf_matrix_var, annotation_col = annotation_col,
         annotation_row = annotation_row_var,
         annotation_colors = ann_colors, scale = "row", show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         fontsize = 12,
         annotation_names_row = FALSE,
         main = "Activity of the 50 most variable transcription factor on the PaCaOMICS PDX dataset")

#Get a df of the variables I want to annotate the rows of the heatmap
annotation_row_cor_ICA3 <- tf_ICA3_cor %>% 
  rownames_to_column("TF") %>%
  inner_join(dplyr::select(top_tf_ICA3, TF), by = "TF") %>%
  dplyr::rename(cor = V1) %>%
  column_to_rownames(var = "TF") 

#Top correlated tf with ICA3
pheatmap(tf_matrix_cor_ICA3, annotation_col = annotation_col,
         annotation_row = annotation_row_cor_ICA3,
         annotation_colors = ann_colors, scale = "row", show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         fontsize = 12,
         annotation_names_row = FALSE,
         main = "Activity of the 50 most correlated transcription factor with ICA3 on the full PaCaOMICS PDX dataset")


#Calculate absolute correlation for each gene
tf_PC1_cor <- tf_act_ %>% 
  dplyr::rename(sample = sample) %>%
  inner_join(dplyr::select(sample_info, sample), by = "sample")  %>% #make the two table in same order
  dplyr::arrange(sample) %>%
  column_to_rownames("sample") %>% #select continuous variables
  as.matrix() %>%
  cor(y = sample_info$PC1, method = "spearman") %>%
  as.data.frame()

#keep the 50 most absolute correlated TF
top_tf_PC1 <- mutate(tf_PC1_cor, V1 = abs(V1)) %>%
  dplyr::arrange(desc(V1)) %>%
  dplyr::slice(1:50) %>%
  dplyr::rename(cor = V1) %>%
  rownames_to_column(var = "TF") 

#Final matrix
tf_matrix_cor_PC1 <-  as_tibble(tf_act__, rownames = "TF" ) %>%
  dplyr::select(TF, sample_info$sample) %>%
  dplyr::filter(TF %in% top_tf_PC1$TF) %>%
  column_to_rownames("TF") %>%
  as.matrix()

#Get a df of the variables I want to annotate the rows of the heatmap
annotation_row_cor_PC1 <- tf_PC1_cor %>% 
  rownames_to_column("TF") %>%
  inner_join(dplyr::select(top_tf_PC1, TF), by = "TF") %>%
  dplyr::rename(cor = V1) %>%
  column_to_rownames(var = "TF") 

#Top correlated tf with PC1
pheatmap(tf_matrix_cor_PC1, annotation_col = annotation_col,
         annotation_row = annotation_row_cor_PC1,
         annotation_colors = ann_colors, scale = "row", show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         fontsize = 12,
         annotation_names_row = FALSE,
         main = "Activity of the 50 most correlated transcription factor with PC1 on the full PaCaOMICS PDX dataset")


#Proportion of common top correlated TF with PC1 and ICA3
length(Reduce(intersect, list(top_tf_ICA3$TF,top_tf_PC1$TF)))/50

