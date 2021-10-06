library(tidyverse)
library(fgsea)
library(GSVA)
library(org.Hs.eg.db)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/04_Deconvolution_TE/03_Sigantures/cell_types/fibroblast/PJ2003085/090921_FastqProcessing")

# preRanked GSEA of all types vs. Mean
read_tsv("02_Output/mean.foldchanges.tsv") -> foldchanges
read_rds("02_Output/multiassay.RDS") -> mae # TMM counts and TE at all levels

## Ensembl to gene names & Entrez

foldchanges_translator <- AnnotationDbi::select(org.Hs.eg.db,
                      keys=foldchanges$Row.names,
                      column=c("ENTREZID", "SYMBOL"),
                      keytype="ENSEMBL",
                      multiVals="first")


foldchanges_genenames = foldchanges
foldchanges_translator[c("ENSEMBL", "SYMBOL")] %>% deframe() -> translate
foldchanges$Row.names %>% translate[.] -> foldchanges_genenames$Row.names

foldchanges_genenames %>% distinct(Row.names, .keep_all= TRUE) %>%
    filter(!is.na(Row.names)) %>%
    column_to_rownames("Row.names") -> foldchanges_genenames


foldchanges_entrez = foldchanges
foldchanges_translator[c("ENSEMBL", "ENTREZID")] %>% deframe() -> translate
foldchanges$Row.names %>% translate[.] -> foldchanges_entrez$Row.names

foldchanges_entrez %>% distinct(Row.names, .keep_all= TRUE) %>%
  filter(!is.na(Row.names)) %>%
  column_to_rownames("Row.names") -> foldchanges_entrez


## GSEA
### Fold changes
for (lvl in c("Total", "Polysome", "TE")){
for (cl in c("s17AAO2007", "s17T", "sDALJO", "sMAYCL")){
comparison <- paste(cl,lvl,sep=".")

#### Stromal subtype
stroma_signatures <- read_rds("01_Input/PDACstromaSignatures.rds")

rank <- foldchanges_genenames %>% arrange(desc(get(comparison))) %>%
  dplyr::select(comparison) %>% 
  filter_all(all_vars(!is.infinite(.)))


ranked_nl <- rank[,comparison, drop =T]
names(ranked_nl) <- rownames(rank)

fgseaRes <- fgsea(pathways = stroma_signatures,
                  stats = ranked_nl)
plot.new()
plotGseaTable(stroma_signatures[grep("CAF", names(stroma_signatures))], ranked_nl, fgseaRes, 
              gseaParam = 0.5)
title(ylab = comparison,   font.lab = 4, col.lab = "red")

#### Reactome pathways
rank <- foldchanges_entrez %>% arrange(desc(get(comparison))) %>%
  dplyr::select(comparison) %>% 
  filter_all(all_vars(!is.infinite(.)))


ranked_nl <- rank[,comparison, drop =T]
names(ranked_nl) <- rownames(rank)


reactome <- reactomePathways(names(ranked_nl))

fgseaRes <- fgsea(pathways = reactome,
                  stats = ranked_nl)

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plot.new()
plotGseaTable(reactome[topPathways], ranked_nl, fgseaRes, 
              gseaParam=0.5)
title(ylab = comparison,   font.lab = 4, col.lab = "red")
}
### GSVA (ssGSEA)
mae_translator <- AnnotationDbi::select(org.Hs.eg.db,
                                                keys=mae@metadata$genes,
                                                column=c("ENTREZID", "SYMBOL"),
                                                keytype="ENSEMBL",
                                                multiVals="first")

#### signatures
mae_lvl_genenames <- mae@ExperimentList[[lvl]]
mae_translator[c("ENSEMBL", "SYMBOL")] %>% deframe() -> translate
rownames(mae_lvl_genenames) %>% translate[.] -> rownames(mae_lvl_genenames)



gsva(mae_lvl_genenames, stroma_signatures) %>% pheatmap(main=paste(lvl, "vs. stromal subtypes"),
                                                        cluster_cols = F, cluster_rows = F)



#### reactome
mae_lvl_entrez <- mae@ExperimentList[[lvl]]
mae_translator[c("ENSEMBL", "ENTREZID")] %>% deframe() -> translate
rownames(mae_lvl_entrez) %>% translate[.] -> rownames(mae_lvl_entrez)

reactome_mae <- reactomePathways(rownames(mae_lvl_entrez))


gsvaRes <- gsva(mae_lvl_entrez, reactome_mae, min.sz = 15)
gsvaRes %>% apply( 1, function(x) {sum(abs(x - mean(x))) / length(x)}) %>%
  order(decreasing = T) -> bestpath

gsvaRes %>% pheatmap(main=paste(lvl, "vs. all pathways"), show_rownames = F )
gsvaRes[bestpath[1:50],] %>% pheatmap(main=paste(lvl, "vs. mostvar pathways"))
}
