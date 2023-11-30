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

## Celltype plot
DimPlot(seurat_object, reduction = "UMAP", group.by ="new_celltypes")


