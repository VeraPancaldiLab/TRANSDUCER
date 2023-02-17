library(tidyverse)
library(reshape2)
library(ggplot2)
#library(ggsignif)
#library(Hmisc)
#library(corrplot)
#library(GSVA)
#library(pheatmap)
library(ComplexHeatmap)
#library(circlize)
#library(RColorBrewer)
#library(edgeR)
#library(scales)
library(factoextra)
library(pdacmolgrad)
library(DESeq2) 
library(EnhancedVolcano)
library(enrichR)


setwd("/home/jacobo/Documents/02_TRANSDUCER/03_IC3Characterization/02_CCLE/")

# Data loading
metadata <- read_tsv("01_Input/CCLE_Pancreatic_meta.tsv")
metadata <- as.data.frame(metadata[0:40,])
plot_interesting <- metadata$CellLine # from 40, they are missing in this datset
interesting_cl <- paste(plot_interesting,'_PANCREAS', sep = "")

rownames(metadata) <- interesting_cl
metadata$CBSexp[is.na(metadata$CBSexp)] <- "other"
metadata$CBSexp[metadata$CBSexp != "other"] <- "neg"
metadata$PHGDHexp[is.na(metadata$PHGDHexp)] <- "other"
metadata$PHGDHexp[metadata$PHGDHexp != "other"] <- "neg"
metadata$PHGDHCBSexp <- interaction(metadata$PHGDHexp, metadata$CBSexp)

## CCLE data
### TMM
ccle_counts <- readRDS("02_Output/CCLE_symbols_counts.RDS") #gene symbols, read counts
ccle_counts <- ccle_counts[,interesting_cl] 

#### normalization
ccle_counts_dge <- DGEList(ccle_counts)
ccle_TMM_dge <- calcNormFactors(ccle_counts_dge,
                                     method = "TMM")

ccle_TMM <- cpm(ccle_TMM_dge)

### TPM
ccle_TPM <- readRDS("02_Output/CCLE_symbols_TPM.RDS") #gene symbols, TPM
ccle_TPM <- ccle_TPM[,interesting_cl[0:40]] # from 40, they are missing in this datset
ccle_TPM_dge <- DGEList(ccle_TPM)


### CHOOSE NORM
ccle <- ccle_TMM
ccle_dge0 <- ccle_TPM_dge
norm <- "TMM"

# ccle <- ccle_TPM
# ccle_dge <- ccle_TPM_dge
# norm <- "TPM"


# Distribution check
ccle_counts.m <- melt(as.matrix(ccle_counts+ 1), measure.vars = 1:ncol(ccle_counts))
ccle.m <- melt(as.matrix(ccle+ 0.01), measure.vars = 1:ncol(ccle))

## Raw
ggplot(ccle_counts.m, aes(x=value, y=Var2)) + 
  geom_boxplot() + 
  #scale_x_continuous(trans = 'log10') + 
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(title = "Normalization check", subtitle = "Raw counts") + 
  xlab("counts") + ylab("") +
  theme_classic()

## TMM
ggplot(ccle.m, aes(x=value, y=Var2)) + 
  geom_boxplot() + 
  #scale_x_continuous(trans = 'log10') + 
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(title = "Normalization check", subtitle = "Normalized") + 
  xlab(norm) + ylab("") +
  theme_classic()


# Check PHGDH vs CBS expresion to find same plot as Yvan
phgdhvscbs <- as.data.frame(t(ccle[c("PHGDH", "CBS", "IDO1"),]))
doublelow_s <- phgdhvscbs
ggplot(phgdhvscbs, aes(x=PHGDH, y=CBS)) + 
  geom_point() + 
  geom_point(data=phgdhvscbs[metadata$PHGDHexp =="PHGDH Neg" & metadata$CBS == "CBS Neg",], aes(x=PHGDH, y=CBS), colour="red", size=2) +
  scale_x_continuous(trans = 'log2') +
 scale_y_continuous(trans = 'log2') +
 ggtitle("log(PHGDH VS CBS)")

# PAMG #!Should do ICA  and correlation with S mat.
ccle.pamg <- projectMolGrad(newexp = ccle, geneSymbols = row.names(ccle))
metadata$PAMG <- ccle.pamg[["PDX"]]

# PCA and PHGDH CBS expression 
ccle.pca <- prcomp(ccle ,scale. = TRUE, center = TRUE) # scaled PCAx

fviz_eig(ccle.pca, barfill = "lightgrey",
         barcolor = "black", ggtheme = theme_void(), title = "PDAC Cell lines PCA")

ccle.pca.toplot <- as.data.frame(ccle.pca["rotation"])
metadata.toplot <- metadata[rownames(ccle.pca.toplot),]


## Metastatic/Primary tumour
plot.new()
col_factor <- as.factor(metadata.toplot$TumorType)

pairs(ccle.pca.toplot[, 1:7],
      col = as.numeric(col_factor),
      pch = 19,
      lower.panel = NULL)
par(xpd = TRUE)
legend("bottom", fill = as.numeric(unique(col_factor)), legend = levels(col_factor))



## PAMG
plot.new()
col_factor <- as.factor(metadata.toplot$PAMG)
colfunc <- colorRampPalette(c("chocolate1", "dodgerblue"))
### check scale
plot(col_factor,col=colfunc(length(col_factor)), main = "Check Color Scale of PAMG")


plot.new()
pairs(ccle.pca.toplot[, 1:7],
      col = colfunc(length(col_factor))[as.numeric(col_factor)],
      main = "PCA CCLE PAMG",
      pch = 19,
      lower.panel = NULL)



## PHGDH/CBS exp
### Cualitative
plot.new()
col_factor <- as.factor(metadata.toplot$PHGDHexp)
pch_factor <- as.factor(metadata.toplot$CBSexp)

pairs(ccle.pca.toplot[, 1:7],
      col = as.numeric(col_factor),
      pch = 17 + 2*as.numeric(pch_factor),
      lower.panel = NULL)
par(xpd = TRUE)
legend("bottomleft", pch = 17 + 2*as.numeric(pch_factor), col = as.numeric(col_factor), legend = levels(pch_factor))
legend("bottom", fill = as.numeric(rev(unique(col_factor))), legend = levels(col_factor))

### Specific Genes Cuantitative
gene = "PHGDH"
plot.new()
col_factor <- as.factor(ccle[gene,rownames(ccle.pca.toplot)])
colfunc <- colorRampPalette(c("black", "red"))
plot(col_factor,col=colfunc(length(col_factor)), main = paste("Check Color Scale of", gene))

pairs(ccle.pca.toplot[, 1:7],
      col = colfunc(length(col_factor))[as.numeric(col_factor)],
      main = paste(gene, "Cuantitative"),
      pch = 19,
      lower.panel = NULL)


# DGE with DEseq2
library(DESeq2)
library(EnhancedVolcano)
all(colnames(ccle_counts) == rownames(metadata)) # if True, the order is alright, an metadata can be used top extract sample factors

dds <- DESeqDataSetFromMatrix(countData = ccle_counts,
                              colData = metadata,
                              design= ~ PHGDHCBSexp)

dds$PHGDHCBSexp <- relevel(dds$PHGDHCBSexp, ref = "other.other")

## Low expressed genes removal (keep it small and fast)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## Analysis
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="PHGDHCBSexp_neg.neg_vs_other.other")
# or to shrink log fold changes association with condition:
resLFC <- lfcShrink(dds, coef="PHGDHCBSexp_neg.neg_vs_other.other")

## Results exploration
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered),
          file="02_Output/CBSPHGDH_negnegvsotherother_PDACCCLE.csv")

  
summary(res)
plotMA(res, ylim=c(-2,2))

EnhancedVolcano(res,
                lab = rownames(res),
                title = "PHGDH&CBS high vs low",
                subtitle = "10 as min total reads.",
                x = 'log2FoldChange',
                y = 'pvalue')

EnhancedVolcano(resLFC,
                lab = rownames(resLFC),
                title = "PHGDH&CBS high vs low",
                subtitle = "10 as min total reads. + shrinked log fold changes",
                x = 'log2FoldChange',
                y = 'pvalue')


### ENrichr (without LFC!)
pval.th <- 0.05
fold_change.th <- 1.5
databases <- c("GO_Biological_Process_2018","KEGG_2019_Human")


res_nonna <- na.omit(res)
degs <- res_nonna[(abs(res_nonna[, "log2FoldChange"]) > fold_change.th) &
                   (res_nonna[, "padj"] < pval.th), ]

#### Underexpressed
underexp <- degs[degs[, "log2FoldChange"] < 0, ]
print(nrow(underexp))
u_enriched <- enrichr(
  genes = rownames(underexp),
  databases = databases
)

View(u_enriched)


#### Overexpressed
overexp <- degs[degs[, "log2FoldChange"] > 0, ]
print(nrow(overexp))
o_enriched <- enrichr(
  genes = rownames(overexp),
  databases = databases
)
View(o_enriched)

### GSEA
#### prepare data. General GSEA
deseq_counts <- counts(dds, normalized=T) 
deseq_counts[is.na(deseq_counts)] <- 0
GENES <- rownames(deseq_counts)
DESCRIPTION <- rep(NaN, nrow(deseq_counts))
deseq_counts <- cbind(GENES, DESCRIPTION, deseq_counts)
write.table(deseq_counts, "02_Output/count_matrix_gsea.txt", 
            sep = "\t", row.names=FALSE,quote = FALSE)


lvls <- paste(c(length(metadata$PHGDHCBSexp), length(unique(metadata$PHGDHCBSexp)), 1), sep = "\t")

cat(lvls,file="02_Output/SerEnz_state.cls", sep="\t")
cat("\n",file="02_Output/SerEnz_state.cls", append = TRUE)
cat(paste(c('#', unique(as.vector(metadata$PHGDHCBSexp)))),file="02_Output/SerEnz_state.cls",
    sep="\t", append = TRUE)

cat("\n",file="02_Output/SerEnz_state.cls", append = TRUE)
cat(paste(metadata$PHGDHCBSexp, sep = "\t"),file="02_Output/SerEnz_state.cls",
    sep="\t",append = TRUE)

#### prepare data. for balance tries
#### 3 times 3 random other.other samples GSEA
##### subset data
oo <- which(metadata$PHGDHCBSexp == "other.other")
set.seed(1)
oo.g1 <- sample(oo, 3)
oo.g2 <- sample(oo, 3)
oo.g3 <- sample(oo, 3)

nn <- which(metadata$PHGDHCBSexp == "neg.neg")
metadata_balance <- metadata[c(oo.g1,oo.g2,oo.g3,nn),]
deseq_counts_balance <- deseq_counts[,c("GENES", "DESCRIPTION", rownames(metadata_balance))]
metadata_balance["PHGDHCBSexp"] <- factor(c(paste("other.other", rep(1:3, each=3), sep = ""), rep("neg.neg",3)))


##### format and export
write.table(deseq_counts_balance, "02_Output/balanced_count_matrix_gsea.txt", 
            sep = "\t", row.names=FALSE,quote = FALSE)


lvls_balance <- paste(c(length(metadata_balance$PHGDHCBSexp), length(unique(metadata_balance$PHGDHCBSexp)), 1), sep = "\t")

cat(lvls_balance,file="02_Output/balanced_SerEnz_state.cls", sep="\t")
cat("\n",file="02_Output/balanced_SerEnz_state.cls", append = TRUE)
cat(paste(c('#', unique(as.vector(metadata_balance$PHGDHCBSexp)))),file="02_Output/balanced_SerEnz_state.cls",
    sep="\t", append = TRUE)

cat("\n",file="02_Output/balanced_SerEnz_state.cls", append = TRUE)
cat(paste(metadata_balance$PHGDHCBSexp, sep = "\t"),file="02_Output/balanced_SerEnz_state.cls",
    sep="\t",append = TRUE)

##### Check subsets Serine expression distribution
cn_balanced <- colnames(deseq_counts_balance)[-c(1,2)]
dcb_plot <- t(deseq_counts_balance[c("PHGDH", "CBS", "IDO1"), cn_balanced])
dcb_plot <- as.data.frame(dcb_plot)
dcb_plot[c("PHGDH", "CBS", "IDO1")] <-  lapply(dcb_plot[c("PHGDH", "CBS", "IDO1")], FUN = as.numeric)
dcb_plot <- merge(dcb_plot, metadata_balance, by = "row.names")

ggplot(dcb_plot, aes(x=PHGDH, y=CBS, color = PHGDHCBSexp)) +
  geom_point(size=2, shape=16) +
  ggtitle(label = "CCLE PDAC Balanced", subtitle = "neg.neg vs other.other subsets (log)") +
  scale_y_continuous(trans = "log2") +
  scale_x_continuous(trans = "log2") +
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"))

#### Add Gene name and description column to GSEA result matrices
reference <- read_tsv("01_Input/GRCh38pt13_Abreviations.tsv")
colnames(reference) <- c("SYMBOL", "Description")
target_files <- list.files(path = "02_Output/GSEA_CBSPHGDH_negnegvsotherother_PDACCCLE/", pattern = "\\.tsv$")
target_files <- target_files[grep("[A-Z]", target_files, ignore.case = FALSE,)] # exclude files without upper case
 


for (f in target_files) {
  f <-  paste ("02_Output/GSEA_CBSPHGDH_negnegvsotherother_PDACCCLE/", f, sep = "")
  f_table <- read_tsv(f)
  OG_columns <- colnames(f_table)
  f_sym <- f_table$SYMBOL
  reference_s <- reference[reference$SYMBOL %in% f_sym,]
  
  f_table <- merge(f_table, reference_s)
  f_table$TITLE <- f_table$Description
  f_table <- f_table[OG_columns]
  write_tsv(x = f_table, file = f)
}
############################## TRASH ###########################################
# # DGE with edgeR
# # https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# ## First Try: 3 lowPHGDHCBS vs the rest
# ### the DGE object (already calcnornmfactors done)
# ccle_dge0
# dim(ccle_dge0)
# 
# ### Filter low-expressed genes
# cutoff <- 50
# drop <- which(apply(cpm(ccle_dge0), 1, max) < cutoff)
# ccle_dge <- ccle_dge0[-drop,] 
# dim(ccle_dge) # number of genes left
# 
# ### get factors (CBS state, PHGDH state)
# all(rownames(ccle_dge$samples) == rownames(metadata)) # if True, the order is alright, an metadata can be used top extract sample factors
# phgdh_st <- metadata$PHGDHexp
# cbs_st <- metadata$CBSexp
# combined_st <- interaction(phgdh_st,cbs_st)
# 
# plotMDS(ccle_dge, col = as.numeric(combined_st))
# 
# ###voom transform nad variance weight calculation
# mm <- model.matrix(~0 + combined_st) #specify model (all groups)
# y <- voom(ccle_dge, mm, plot = T)
