library(Seurat)
library(tidyverse)

################################################################################

# Load data
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
expression_matrix <- ReadMtx(
  mtx = "data/Hwang_Nature_2022/snRNAseq_naive/gene_sorted-naivedata_scp.mtx",
  features = "data/Hwang_Nature_2022/snRNAseq_naive/naivedata_scp.genes.csv",
  cells = "data/Hwang_Nature_2022/snRNAseq_naive/naivedata_scp.barcodes.csv",
  feature.column = 1,
  skip.cell = 1,
  strip.suffix = TRUE
)

expression_matrix@Dimnames[[2]] = str_remove(string = expression_matrix@Dimnames[[2]], pattern="^\\d+,") # remove rownum included in rownames

seurat_object <- CreateSeuratObject(counts = expression_matrix)

metadata <- read_tsv("data/Hwang_Nature_2022/snRNAseq_naive/combinenaivedata-reprocessed-clean-detailed-annotations.tsv") %>% 
   dplyr::filter(NAME != "TYPE") %>% column_to_rownames("NAME") 

seurat_object <- AddMetaData(seurat_object, metadata = metadata)

seurat_object$orig.ident <- seurat_object$pid

# QC and filtering? 
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by= "orig.ident", ncol = 3)


# UMAP plot
## Load
UMAP_Coord <- read_tsv("data/Hwang_Nature_2022/snRNAseq_naive/combinenaivedata-reprocessed-clean-detailed-UMAP.tsv", col_types = c("cnn")) %>% 
  dplyr::filter(NAME != "TYPE") %>%
  rename(UMAP_1 = "X",UMAP_2 = "Y") %>% 
  column_to_rownames("NAME") 


UMAP_coordinates_mat <- as(UMAP_Coord, "matrix")

seurat_object[['UMAP']] <- CreateDimReducObject(embeddings = UMAP_coordinates_mat, key = "UMAP_", global = T, assay = "RNA")

## Celltype plot
DimPlot(seurat_object, reduction = "UMAP", group.by ="cell_subsets" )

## Plot genes 
FeaturePlot(seurat_object, slot = "counts", features = c("ATF4", "IL18", "IFNG", "PHGDH", "SERPINB8", "CBS", "YTHDF3"))
