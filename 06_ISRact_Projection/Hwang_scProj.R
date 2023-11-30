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

## clean unused data
rm(expression_matrix)
rm(metadata)
gc()

# QC and filtering? 
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by= "orig.ident", ncol = 3)

seurat_object[["is_lowQ"]] <- if_else(seurat_object$pid %in% c("003_10x", "2276_10x", "2591_10x"),seurat_object$pid, "no")
# UMAP plot
## Load
UMAP_Coord <- read_tsv("data/Hwang_Nature_2022/snRNAseq_naive/combinenaivedata-reprocessed-clean-detailed-UMAP.tsv", col_types = c("cnn")) %>% 
  dplyr::filter(NAME != "TYPE") %>%
  rename(UMAP_1 = "X",UMAP_2 = "Y") %>% 
  column_to_rownames("NAME") 


UMAP_coordinates_mat <- as(UMAP_Coord, "matrix")

seurat_object[['UMAP']] <- CreateDimReducObject(embeddings = UMAP_coordinates_mat, key = "UMAP_", global = T, assay = "RNA")

## Celltype plot
DimPlot(seurat_object, reduction = "UMAP", group.by ="is_lowQ" )

## Plot genes 
FeaturePlot(seurat_object, slot = "counts", features = c("ATF4", "IL18", "IFNG", "PHGDH", "SERPINB8", "CBS", "YTHDF3", "ACTA2"))


# Signature projection
## Step 1: Load differentially expressed genes (DEGs)
bulk_sign <- Get DEGs lists, extracting gene symbols from DEGs tables.

## Step 2: Load and process scRNAseq (done above!)
## Step 3: Subset scRNA-seq reference dataset to match cell-type populations in the bulk-tissue RNA-seq
### Filter scRNAseq data to cell types contained in dataset
if (filter == TRUE){
seurat_object <- subset(x = seurat_object,
                            cells = colnames(seurat_object)[seurat_object$cell_subsets %in%
                                                          c("Tumor")])
}

table(seurat_object$cell_subsets)

### Curate labels
Idents(seurat_object) <- "cell_subsets"
Idents(seurat_object) <- factor(Idents(seurat_object),
                                levels = sort(levels(seurat_object)))
### Normalize and Scale
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- ScaleData(seurat_object)

### Plot data
#### PCA
seurat_object <- RunPCA(seurat_object)
ElbowPlot(seurat_object, ndims = 50)

#### UMAP and tSNE
seurat_object <- RunTSNE(seurat_object, dims = 1:20)
seurat_object <- RunUMAP(seurat_object, dims = 1:20)

pca_plot <- DimPlot(seurat_object, reduction = "pca", pt.size = 0.1, label = T)
tsne_plot <- DimPlot(seurat_object, reduction = "tsne", pt.size = 0.1, label = T)
umap_plot <- DimPlot(seurat_object, reduction = "umap", pt.size = 0.1, label = T)
legend <- get_legend(umap_plot)

## Step 4:  Perform the deconvolution of cell-type-specific signal in bulk-tissue RNA-seq top DEGs using single-cell RNA-seq reference




################################################################################ (i can see an orphan dot a
#Treated
expression_matrix <- ReadMtx(
  mtx = "data/Hwang_Nature_2022/snRNAseq_treated/gene_sorted-treated_data_scp.mtx",
  features = "data/Hwang_Nature_2022/snRNAseq_treated/treateddata_scp.genes.csv",
  cells = "data/Hwang_Nature_2022/snRNAseq_treated/treateddata_scp.barcodes.csv",
  feature.column = 1,
  skip.cell = 0,
  strip.suffix = TRUE
)

seurat_object <- CreateSeuratObject(counts = expression_matrix)

metadata <- read_tsv("data/Hwang_Nature_2022/snRNAseq_treated/combinetreateddata-reprocessed-final-annotations.tsv") %>% 
  dplyr::filter(NAME != "TYPE") %>% column_to_rownames("NAME") 

seurat_object <- AddMetaData(seurat_object, metadata = metadata)

seurat_object$orig.ident <- seurat_object$pid

## clean unused data
rm(expression_matrix)
rm(metadata)
gc()

# QC and filtering? 
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by= "orig.ident", ncol = 3)

# UMAP plot
## Load
UMAP_Coord <- read_tsv("data/Hwang_Nature_2022/snRNAseq_treated/combinetreateddata-reprocessed-final-UMAP.tsv", col_types = c("cnn")) %>% 
  dplyr::filter(NAME != "TYPE") %>%
  rename(UMAP_1 = "X",UMAP_2 = "Y") %>% 
  column_to_rownames("NAME") 


UMAP_coordinates_mat <- as(UMAP_Coord, "matrix")

seurat_object[['UMAP']] <- CreateDimReducObject(embeddings = UMAP_coordinates_mat, key = "UMAP_", global = T, assay = "RNA")

## Celltype plot
DimPlot(seurat_object, reduction = "UMAP", group.by ="cell_subsets")

## Plot genes 
FeaturePlot(seurat_object, slot = "counts", features = c("ATF4", "IL18", "IFNG", "PHGDH", "SERPINB8", "CBS", "YTHDF3", "ACTA2"))


# Matthieu 
library(GEOquery)
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
GSE202051 <- getGEO('GSE202051',GSEMatrix=TRUE)

seuratObject <- LoadH5Seurat("GSE202051_totaldata-final-toshare.h5ad")


library(anndata)
seuratObject <- read_h5ad("GSE202051_totaldata-final-toshare.h5ad", backed = NULL)

# this is anndata notseurat. Needs to be imported to seurat
GSE202051_seurat <- sceasy::convertFormat(obj = "GSE202051_totaldata-final-toshare.h5ad",
                      from = "anndata",
                      to = "seurat")

# try to getting Ann data to Seurat
##https://satijalab.org/seurat/archive/v2.4/conversion_vignette
library(reticulate)
ad <- import("anndata", convert = FALSE)
GSE202051_ad <- ad$read_h5ad("GSE202051_totaldata-final-toshare.h5ad")
GSE202051_seurat <- Convert(GSE202051_ad, to = "seurat" , assay = "RNA")


##https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html 
SeuratDisk::Convert("GSE202051_totaldata-final-toshare.h5ad", dest = "h5seurat", overwrite = TRUE)
GSE202051_seurat <- LoadH5Seurat("GSE202051_totaldata-final-toshare.h5seurat",assays = "RNA")

## youtube tutorial using python WORKS!!!
expression_matrix <- ReadMtx(
  mtx = "data/Hwang_Nature_2022/GEO/matrix.mtx",
  features = "data/Hwang_Nature_2022/GEO/features.tsv",
  cells = "data/Hwang_Nature_2022/GEO/barcodes.tsv",
  skip.cell = 0
)


seurat_object <- CreateSeuratObject(counts = expression_matrix)

metadata <- read_tsv("data/Hwang_Nature_2022/GEO/metadata.tsv")  %>% column_to_rownames("...1") 

seurat_object <- AddMetaData(seurat_object, metadata = metadata)

### clean unused data
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





