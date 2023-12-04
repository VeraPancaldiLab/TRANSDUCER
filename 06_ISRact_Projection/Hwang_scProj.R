library(Seurat)
library(tidyverse)

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


