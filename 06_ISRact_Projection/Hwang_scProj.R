library(Seurat)
library(tidyverse)
library(biomaRt)
################################################################################
# Load data and build Seurat Object
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")

expression_matrix <- ReadMtx(
  mtx = "data/Hwang_Nature_2022/GEO/matrix.mtx",
  features = "data/Hwang_Nature_2022/GEO/features.tsv",
  cells = "data/Hwang_Nature_2022/GEO/barcodes.tsv",
  skip.cell = 0
)

metadata <- read_tsv("data/Hwang_Nature_2022/GEO/metadata.tsv")  %>% column_to_rownames("...1") 
seurat_object <- CreateSeuratObject(counts = expression_matrix)
seurat_object <- AddMetaData(seurat_object, metadata = metadata)

## clean unused data
rm(expression_matrix)
rm(metadata)
gc()


# UMAP plot
## Load
UMAP_Coord <- read_tsv("data/Hwang_Nature_2022/GEO/obsm/X_umap.tsv", col_types = c("cnn")) %>%
  rename(UMAP_1 = "0",UMAP_2 = "1") %>% 
  column_to_rownames("...1") 


UMAP_coordinates_mat <- as(UMAP_Coord, "matrix")

seurat_object[['UMAP']] <- CreateDimReducObject(embeddings = UMAP_coordinates_mat, key = "UMAP_", global = T, assay = "RNA")

## Create my own UMAP
# ### Normalise (log2(TP10K+1)) DOUBTS!
# seurat_object <- NormalizeData(seurat_object,  normalization.method = "LogNormalize", scale.factor = 10000, verbose = T)
# 
# ### Select and scale high variable genes (2K) DOUBTS!
# seurat_object <- FindVariableFeatures(seurat_object, selection.method = "dispersion", nfeatures = 2000)
# 
# top10 <- head(VariableFeatures(seurat_object), 10)
# plot1 <- VariableFeaturePlot(seurat_object)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# 
# ### Scale data (not mentioned methods)
# all.genes <- rownames(seurat_object)
# seurat_object <- ScaleData(seurat_object, features = all.genes)
# 
# ### PCA and num dim exploration FAILED FRO MEMORY ISSUE
# seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), verbose = T)
# DimHeatmap(seurat_object, dims = 1:15, cells = 500, balanced = TRUE)
# ElbowPlot(seurat_object)
# 
# 
# seurat_object <- RunUMAP(seurat_object, dims = 1:10)

## Celltype plot
DimPlot(seurat_object, reduction = "UMAP", group.by ="new_celltypes")
DimPlot(seurat_object, reduction = "UMAP", group.by ="pid")
DimPlot(seurat_object, reduction = "UMAP", group.by ="Level 1 Annotation")
DimPlot(seurat_object, reduction = "UMAP", group.by ="Level 2 Annotation")
DimPlot(seurat_object, reduction = "UMAP", group.by ="Level 3 Annotation")

# ISRact Projection
## Load ISRAct signature to explore in the snRNAseq
pca_pdx <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS")
ISRact_contributions <- sort(pca_pdx$rotation[,"PC1"])
### translate to Gene ID
# Translate EnsemblID to gene names
## Version 75 for PDX data
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)#listAttributes(ensembl75, page="feature_page")

annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = ensembl75)

translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

ISRact_info <- tibble(EnsemblID = names(ISRact_contributions),
                      gene_name = translate[names(ISRact_contributions)],
                      PC1_score = ISRact_contributions,
                      class = if_else(ISRact_contributions > 0, "ISRac_high", "ISRac_low"))

ggplot(ISRact_info, aes(x = PC1_score)) +
  geom_density() +
  geom_rug(aes(color = class)) +
  theme_bw()
### extract genes

ISRact_high <- dplyr::filter(ISRact_info, class == "ISRac_high") %>%
  dplyr::select("gene_name") %>% 
  deframe()

ISRact_low <- dplyr::filter(ISRact_info, class == "ISRac_low") %>%
  dplyr::select("gene_name") %>% 
  deframe()

## subset scRNAseq data for matching bulk composition
seurat_object <- subset(x = seurat_object, 
                  subset = new_celltypes == "Epithelial-Malignant")

table(seurat_object$new_celltypes)
table(seurat_object$treatment_status)

### create new identity status (what you want to stratify based on)
identity <- "treatment_status" # treatment_status | 
Idents(seurat_object) <- identity
Idents(seurat_object) <- factor(Idents(seurat_object),levels = sort(levels(seurat_object)))

### Normalize, find variable genes, scale and center
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000) 
seurat_object <- ScaleData(seurat_object) 

### PCA + UMAP & tSNE
seurat_object <- RunPCA(seurat_object)
ElbowPlot(seurat_object, ndims = 50)

seurat_object <- RunTSNE(seurat_object, dims = 1:40)
seurat_object <- RunUMAP(seurat_object, dims = 1:40)

#### plot
library(cowplot)
library(patchwork)
library(ggplot2)

cols <- c("limegreen", #CRT
          "steelblue", #CRTl
          "mediumorchid4", #CRTln
          "yellow3", #CRTn
          "firebrick2", #CRTx
          "magenta", #GART
          "tan2", #RT
          "gray52") #Untreated

pca_plot <- DimPlot(seurat_object, reduction = "pca",pt.size = 0.1, label = T, cols = cols)
tsne_plot <- DimPlot(seurat_object, reduction = "tsne",pt.size = 0.1, label = T, cols = cols)
umap_plot <- DimPlot(seurat_object, reduction = "umap",pt.size = 0.1, label = T, cols = cols)

legend <- get_legend(umap_plot)
layout <- c (area(1, 1, 1, 2), #PCA
             area(1, 3, 1, 4), #tSNE
             area(1, 5, 1, 6), #UMAP
             area(1, 7)) #Legend

fig1 <- pca_plot + tsne_plot + umap_plot + legend + plot_layout(design = layout) & NoLegend()
fig1

##### plot alternative factors
factor = "batch"  # "response | pid
pca_plot <- DimPlot(seurat_object, reduction = "pca",pt.size = 0.1, label = T, group.by = factor)
tsne_plot <- DimPlot(seurat_object, reduction = "tsne",pt.size = 0.1, label = T,group.by = factor)
umap_plot <- DimPlot(seurat_object, reduction = "umap",pt.size = 0.1, label = T, group.by = factor)

legend <- get_legend(umap_plot)
layout <- c (area(1, 1, 1, 2), #PCA
             area(1, 3, 1, 4), #tSNE
             area(1, 5, 1, 6), #UMAP
             area(1, 7)) #Legend

fig1 <- pca_plot + tsne_plot + umap_plot + legend + plot_layout(design = layout) & NoLegend()
fig1

### Save and reload object for later
saveRDS(seurat_object, "data/Hwang_Nature_2022/Snapshot_Hwang_Processed.RDS")
seurat_object <- readRDS("data/Hwang_Nature_2022/Snapshot_Hwang_Processed.RDS")

## Perform the deconvolution of cell-type-specific signal in bulk-tissue RNA-seq top DEGs using single-cell RNA-seq reference
### get the gene intersection between both lists
sc_gene_names <- rownames(seurat_object)
ISRact_low <- ISRact_low[ISRact_low %in% sc_gene_names]
ISRact_high <- ISRact_high[ISRact_high %in% sc_gene_names]

### equal to gene sets of the same size
min_size = min(length(ISRact_low), length(ISRact_high))

ISRact_low <- head(ISRact_low, min_size)
ISRact_high <- rev(tail(ISRact_high, min_size))

dplyr::mutate(ISRact_info, new_class = if_else(gene_name %in% c(ISRact_low, ISRact_high), class, NA)) %>% 
  ggplot(aes(x = PC1_score)) +
  geom_density() +
  geom_rug(aes(color = new_class)) +
  theme_bw()

### select top_DEG for detection + filter and process scRNAseq
topDEGs_list = ISRact_high # ISRact_high | ISRact_low

seurat_object <- ScaleData(seurat_object, features = topDEGs_list)
seurat_object <- RunPCA(seurat_object, features = topDEGs_list)
topDEGs_list <- rownames(seurat_object@reductions[["pca"]]@feature.loadings)

### Generate dimensionality reduction Plot only with these genes
pca_plot <- DimPlot(seurat_object, reduction = "pca", pt.size = 0.1, label = F, cols = cols)
legend <- get_legend(pca_plot)

pca_plot <- pca_plot & NoLegend()
pca_plot
plot(legend)

### Pairwise correlations and hierarchical clustering
#### random matrix creation
random.matrix <- matrix(runif(500, min = -1, max = 1), nrow = 50)

#### sample quantiles corresponding to probabilities
quantile.range <- quantile(random.matrix, probs = seq(0, 1, 0.01))

#### plot details like palette break definition
palette.breaks <- seq(quantile.range["35%"], quantile.range["83%"], 0.06)
color.palette <- colorRampPalette(c("#0571b0","#f7f7f7","#ca0020"))(length(palette.breaks)-1)

library(gplots)

clustFunction <- function(x)
  
  hclust(as.dist(1-cor(t(as.matrix(x)),method = "pearson")), method = "average")

heatmapPearson <- function(correlations)
  
  heatmap.2(x = correlations,
            
            col = color.palette,
            
            breaks = palette.breaks,
            
            trace = "none", symm = T,
            
            hclustfun = clustFunction)

correlations_DEGs_log <- cor(method = "pearson",
                             
                             log2(t(as.matrix(FetchData(object = seurat_object, vars = c(topDEGs_list), layer = "data")))+1))

heatmapPearson(correlations_DEGs_log)
