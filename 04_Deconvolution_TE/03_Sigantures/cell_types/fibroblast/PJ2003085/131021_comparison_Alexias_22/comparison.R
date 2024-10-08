library(tidyverse)
library(edgeR)
library(sva)
library(factoextra)
library(pheatmap)
library(org.Hs.eg.db)
library(ggvenn) # devtools::install_github("yanlinlin82/ggvenn")
library(RColorBrewer)
library(msigdbr)
library(GSVA)
################################################################################
################################PARAMETERS######################################
joint = F # include the last set of CAFs sequenced asideAn innervated, vascularized and immunocompetent human skin model to study cutaneous neuroimmune interactions
correct_batch = T # Should correct for batch effect?
signatures_magali = F # version sent by Magali 11/10/21
################################################################################
#################################FUNCTIONS######################################

#' Filter a dataframe to keep just half the most variable genes
#'@description
#' `Get_half_mostvar` requires a data.frame with genes as rows,
#'  and its ids as rownames. returns a dataframe

Get_mostvar <- function(df, n){

  abdv <- apply(df, 1, function(x) {
    sum(
      abs(
        x - mean(x)
      )
    ) / length(x)
  })
  
  # to input either % or absolute n
  if ( n>0 && n<1){
    th_abdv.i <- nrow(df)*n
  } else{
    th_abdv.i <- n
  }
    
  abdv %>% sort(decreasing = T) %>% .[th_abdv.i] -> th_abdv
  df.f <- df[abdv >= th_abdv, , drop=F]
  return(df.f)
}


#' Filter a dataframe to keep genes with at least a defined % of non 0 expression samples
#'@description
#' `Exclude_0s` requires a data.frame with genes as rows,
#'  and its ids as rownames. returns a dataframe

Exclude_0s <- function(df, threshold){
  
  n0s <- apply(df, 1, function(x) {
    sum(
        x > 0
    )/length(x)
  })
  

  df.f <- df[which(n0s > threshold), ]
  return(df.f)
}

ID_converter <- function(df, # dataset with ProbeIDs as rownames
                         annotation_table, # Bioconductor annotation table ie. AnnotationDbi::select(hgu219.db, probes, c("SYMBOL", "ENSEMBL", "GENENAME"))
                         old_IDs, # Current IDs ("SYMBOL", "ENSEMBL", "GENENAME")
                         new_IDs # Desired IDs ("SYMBOL", "ENSEMBL", "GENENAME")
)
  #TBI: choose function of aggregation
{
  final_df <- merge(df,annotation_table, by.x=0, by.y=old_IDs)
  
  # Stats of conversion
  non_agg <- nrow(final_df)
  non_agg_uniq <- length(unique(final_df[new_IDs]))
  non_agg_nas <- sum(is.na(final_df[new_IDs]))
  non_agg_nonas <- non_agg-non_agg_nas
  
  # Conversion
  final_df <- aggregate(final_df,
                        final_df[new_IDs],
                        FUN = mean) %>% 
    column_to_rownames(new_IDs) %>% 
    subset(select = colnames(df))
  
  agg <- nrow(final_df)
  
  print(paste(100*(non_agg-agg)/non_agg,
              "% of the originally merged df have been agregated/removed"))
  
  print(paste(100*(non_agg_nas)/non_agg,
              "% of the originally merged df have been removed due no", new_IDs))
  
  print(paste(100*(non_agg_nonas - agg)/non_agg,
              "% of the non NAs df have been aggregated"))
  
  return(final_df)
}
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/04_Deconvolution_TE/03_Sigantures/cell_types/fibroblast/PJ2003085/131021_comparison_Alexias_22/")

# Data loading
read_tsv("01_Input/CAF_rawcounts.tsv") %>%
  column_to_rownames("EnsemblID") -> acafs.raw

acafs.info <- read_tsv("01_Input/CAF_metadata.tsv") %>% 
  column_to_rownames("rawID")


P85cafs.raw <- read_tsv("01_Input/PJ2003085_rawcounts.tsv") %>% 
  column_to_rownames("EnsemblID")

P85cafs.info <- read_tsv("01_Input/PJ2003085_metadata.tsv")


## Assembly check
common_ensemblids <- intersect(rownames(acafs.raw), rownames(P85cafs.raw))
c(length(common_ensemblids), " out of ",
  min(c(nrow(acafs.raw), nrow(P85cafs.raw))),
  " (smaller dataset) Ensembl IDs are common") %>% paste(collapse = "") %>% print()

# Processing
## Common IDs and sample exclusion
acafs.raw <- acafs.raw[common_ensemblids, rownames(acafs.info)]
colnames(acafs.raw) <- acafs.info[,"Official_name"]

totalRNAcols <- P85cafs.info[P85cafs.info$fraction == "Total", "sample", drop = T]
P85cafs.raw <- P85cafs.raw[common_ensemblids,totalRNAcols]
P85cafs.info <- P85cafs.info[P85cafs.info$sample %in% totalRNAcols,]

## Filtering by 0s & TMM norm
Exclude_0s(acafs.raw, 0.5) %>% DGEList() %>% 
  calcNormFactors(method = "TMM") %>% cpm() %>% 
  as_tibble(rownames = "EnsemblID") -> acafs.tmm

Exclude_0s(P85cafs.raw, 0.5) %>% DGEList() %>% 
  calcNormFactors(method = "TMM") %>% cpm() %>% 
  as_tibble(rownames = "EnsemblID") -> P85cafs.tmm

## joint or not 
cafs.info <- acafs.info
rownames(cafs.info) <- NULL
cafs.info %>% column_to_rownames("Official_name") -> cafs.info # Official_name or Name

if (joint == TRUE){
  cafs.tmm <- merge(acafs.tmm, P85cafs.tmm, by = "EnsemblID")
  cafs.info[P85cafs.info$sample, "extraction_group"] = "poly"
  title_prefix <- "joint"
} else {
  cafs.tmm <- acafs.tmm
  title_prefix <- "just Alexia's"
}

## batch effect correction (should be done before filtering)
if (correct_batch == F){
  cafs.tmm %>% column_to_rownames("EnsemblID") -> cafs.choose
  batch <- "non corrected"
  
} else if (correct_batch == T){
  
  cafs.tmm %>% column_to_rownames("EnsemblID") %>% as.matrix() %>% 
    ComBat(batch = cafs.info$extraction_group, 
           par.prior=TRUE, 
           prior.plots=FALSE) -> cafs.tmm.nobatch
  
  batch <- "batch corrected"
  cafs.tmm.nobatch -> cafs.choose
}

# Analysis
title_res <- paste(title_prefix, batch, collapse = " ")
annot <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(cafs.choose), columns="SYMBOL", keytype="ENSEMBL")
cafs.choose.sym <- ID_converter(df = cafs.choose,annotation_table = annot,
                                old_IDs = "ENSEMBL", new_IDs = "SYMBOL")

## General Heatmap of the 21k most variable genes
full_colanot <- c("proliferation", "passes", "Inmortalized", "tissue_origin",
                  "sex", "age", "tumor_differentiation", "initial_sample",
                  "treatment_neoadjuvant", "OS")

cafs.choose.sym %>% Get_mostvar(n = 1000) %>% pheatmap(main=title_res, scale = "row",  treeheight_row = 0,  cellwidth=15, cellheight=0.6, filename = "02_Output/plots/1k_mostvar.png",
                                                     cluster_cols = T, cluster_rows = T, show_rownames = F, annotation_col = cafs.info[full_colanot])
dev.off()


cafs.choose.sym %>% Get_mostvar(n = 50) %>% pheatmap(main=title_res,scale = "row",  cellwidth=15, cellheight=10, filename = "02_Output/plots/50_mostvar.png",
                                                   cluster_cols = T, cluster_rows = T, show_rownames = T, annotation_col = cafs.info[full_colanot])
dev.off()

## Barplot of whatever gene
gene = "SHMT1"
cafs.gene.expression <- as_tibble(cafs.choose.sym, rownames = "gene") %>% 
  pivot_longer(cols = -gene ,names_to = "CL",values_to = "expression" ) %>%
  pivot_wider(id_cols = "CL", names_from = "gene", values_from = "expression")

dplyr::mutate(cafs.gene.expression, CL = fct(CL, levels = arrange(cafs.gene.expression, desc(get(gene)))$CL)) %>% ggplot() + aes_string(x = gene, y = "CL") + geom_bar(stat = "identity")

## GSVA 
### General
all_genesets <- msigdbr("human")
reactome <- dplyr::filter(all_genesets, gs_subcat == "CP:REACTOME") %>% dplyr::mutate(gs_name = str_remove(gs_name,"REACTOME_")) %>% split(x = .$gene_symbol, f = .$gs_name)
column_annot <- dplyr::select(cafs.info, passes, proliferation)

gsvaRes <- gsva(cafs.choose.sym %>% data.matrix(), reactome, min.sz = 15) 
gsvaRes %>% Get_mostvar(50) %>% pheatmap(main=title_res,  cellwidth=10, cellheight=10, annotation_col = cafs.info[c("proliferation", "passes", "Inmortalized", "tissue_origin")],
                                         filename = "02_Output/plots/Reactomemostvar.png")


### stromal subtypes
#### Elayada et al signatures
if (signatures_magali == F){
  
  stroma_signatures <- read_rds("01_Input/PDACstromaSignatures.rds")
  caf.all <- c("iCAFcult", "iCAFsc",
               "myCAFsc", "myCAFcult")
  
  caf.labs <- c("myCAFsc", "iCAFsc") # this ones are newest and most reliable
} else
  #Magali richards signatures
  if (signatures_magali == T){

    pdac_signatures <- read_rds("01_Input/CellTypes_Markers_By_level.rds")
    pdac_CAF_signatures <- pdac_signatures$lvl3$fibro_stell
    
    stroma_signatures <- pdac_CAF_signatures
    caf.labs <- names(pdac_CAF_signatures)
}

#### Espinet et al 2019 signature
IFNsign <- read_tsv("01_Input/IFNsign_Espinet.tsv", col_names = "IFNsign")
stroma_signatures$IFNsign <- IFNsign$IFNsign

#### Huocong et al 2022 apCAF signature (adaptation human by case ~)
apCAF_Huocong <- read_tsv("01_Input/apCAF_Huocong.tsv", col_names = "apCAF_Huocong")
stroma_signatures$apCAF_Huocong <- str_to_upper(apCAF_Huocong$apCAF_Huocong)

#### Luo Nature Communications 2022 (Not only CAFs)
Luo_info <- read_tsv("01_Input/Luo_NatComs_2022_info.tsv")
Luo_clusters <- read_tsv("01_Input/Luo_NatComs_2022.tsv") %>%
  rename_with(~Luo_info$name) %>% 
  dplyr::select(Luo_info[Luo_info$Fibroblast == TRUE,]$name)

stroma_signatures <- c(stroma_signatures, Luo_clusters)

### Venn of the signature composition
myCol <- brewer.pal(4, "Pastel2")

ggvenn::ggvenn(data = stroma_signatures[caf.labs], fill_color = myCol)

### Analysis perse
gsvaRes <- gsva(cafs.choose.sym %>% data.matrix(), stroma_signatures) 
gsvaRes[caf.labs,] -> gsvaRes.cafs

gsvaRes %>% pheatmap(main=title_res,
                           cluster_cols = T, cluster_rows = T, annotation_col = cafs.info[full_colanot])


### Signature enrichment
#svg("02_Output/CAFs_class.svg")
gsvaRes.cafs %>% 
  pheatmap(main=title_res, cluster_cols = T, cluster_rows = T, cellwidth = 10, cellheight = 40,
           annotation_col = cafs.info[c("proliferation", "passes", "Inmortalized", "tissue_origin")])
#dev.off()

gsvaRes.cafs %>% t() -> cafs.info[, caf.labs]


### Marker genes + signature enrichments
#### function 
####(by default will drop duplicated genes if any)
Get_geneset_genes <- function(data, gene_set, gene_set_name){
  
  gene_set_genes <- tibble(name = gene_set_name, value = gene_set[gene_set_name]) %>%
    unnest(c("value")) %>%
    filter(value %in% rownames(data)) %>%
    group_by(value) %>% 
    mutate(name = ifelse(base::duplicated(value,fromLast=T), "multiple", name)) %>%
    distinct(value, .keep_all = T,) %>%
    arrange(name) %>% 
    column_to_rownames("value")
}

#### Elayada et al 2019 (myCAF iCAF)
elayada2019_genes <- Get_geneset_genes(data = cafs.choose.sym,
                                 gene_set = stroma_signatures,
                                 gene_set_name = caf.labs)

cafs.choose.sym[rownames(elayada2019_genes),] %>% 
  pheatmap(main=title_res,  cellwidth=15, cellheight=1, filename = "02_Output/plots/PDAC_elayada2019_markergenes.png",
           cluster_cols = T, cluster_rows = F, scale = "row", show_rownames = F, annotation_col = as.data.frame(t(gsvaRes[caf.labs,])),
           annotation_row = elayada2019_genes)
dev.off()

#### Espinet et al 2019 (IFNsign)
espinet2019_genes <- Get_geneset_genes(data = cafs.choose.sym,
                                       gene_set = stroma_signatures,
                                       gene_set_name = "IFNsign")

cafs.choose.sym[rownames(espinet2019_genes),] %>% 
  pheatmap(main=title_res,  cellwidth=15, cellheight=10, filename = "02_Output/plots/PDAC_espinet2019_markergenes.png",
           cluster_cols = T, cluster_rows = F, scale = "row", show_rownames = T, annotation_col = as.data.frame(t(gsvaRes[c("IFNsign", caf.labs),])),
           annotation_row = espinet2019_genes)
dev.off()

#### Huocong et al 2022 apCAF
huocong2022_genes <- Get_geneset_genes(data = cafs.choose.sym,
                                       gene_set = stroma_signatures,
                                       gene_set_name = "apCAF_Huocong")

cafs.choose.sym[rownames(huocong2022_genes),] %>% 
  pheatmap(main=title_res,  cellwidth=15, cellheight=10, filename = "02_Output/plots/PDAC_huocong2022_markergenes.png",
           cluster_cols = T, cluster_rows = F, scale = "row", show_rownames = T, annotation_col = as.data.frame(t(gsvaRes[c("apCAF_Huocong", caf.labs),])),
           annotation_row = huocong2022_genes)
dev.off()

#### Luo Nature Communications 2022 Pn cancer CAFs
Luo2022_genes <- Get_geneset_genes(data = cafs.choose.sym,
                                       gene_set = stroma_signatures,
                                       gene_set_name = names(Luo_clusters))

cafs.choose.sym[rownames(Luo2022_genes),] %>% 
  pheatmap(main=title_res,  cellwidth=15, cellheight=10, filename = "02_Output/plots/PDAC_Luo2022_markergenes.png",
           cluster_cols = T, cluster_rows = F, show_rownames = T, annotation_col = as.data.frame(t(gsvaRes[c(names(Luo_clusters), caf.labs),])),
           annotation_row = Luo2022_genes)
dev.off()


## PCA analysis
cafs.choose %>% prcomp() -> pca_res

### plots
pca_res %>% fviz_eig(main = title_res)
all(rownames(pca_res$rotation) == rownames(cafs.info)) %>% stopifnot()
pca_toplot <- cbind(pca_res$rotation, cafs.info)

### extraction group
pca_toplot %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2, color = extraction_group)) + geom_point(size=4) +
  theme_bw(base_size=15) + 
  theme(legend.position="top") +
  labs(title = title_res) 

### sample names
pca_toplot %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2, label = rownames(.), color = extraction_group)) + geom_text(size=4) +
  theme_bw(base_size=15) + 
  theme(legend.position="top") +
  labs(title = title_res) 

### proliferation
pca_toplot %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2, color = proliferation)) + geom_point(size=4) +
  theme_bw(base_size=15) + 
  theme(legend.position="top") +
  labs(title = title_res) 
