library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggsignif)
library(Hmisc)
library(corrplot)
library(GSVA)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(edgeR)
library(scales)
library(factoextra)
library(biomaRt)
'%ni%' <- Negate('%in%')
setwd("/home/jacobo/Documents/03_Sauyeun_paper/02_CCLE/")

interesting_cl <- c(
  "ASPC1_PANCREAS", "HUPT3_PANCREAS", "PANC1005_PANCREAS",
  "CAPAN2_PANCREAS", "PATU8902_PANCREAS", "SW1990_PANCREAS", "PATU8988S_PANCREAS", # this is an outlier
  "DANG_PANCREAS", "PANC0327_PANCREAS", "HPAC_PANCREAS",
  "HUPT4_PANCREAS", "PANC0203_PANCREAS", "MIAPACA2_PANCREAS",
  "PATU8988T_PANCREAS", "PANC1_PANCREAS"
) # No PL45_PANCREAS (absent)

plot_interesting <- gsub("_PANCREAS", "", interesting_cl)


# Data loading
## Metadata
metadata <- read_tsv("01_Input/metadata_Bcl_ST1.tsv")
colnames(metadata) <- c("cell_line", "Ser", "cancer_type",
                        "tissue", "gender", "ethnicity", "age")

metadata$cell_line <- str_replace_all(metadata$cell_line, "[^[:alnum:]]", "") %>%
  toupper()

metadata <- which(metadata$cell_line %in% plot_interesting) %>%
  metadata[., colnames(metadata)]

all(metadata$cell_line == plot_interesting[which(plot_interesting %in% 
                                                   metadata$cell_line)]) # THIS MUST BE TRUE SO YOU CAN REPLACE NAMES
metadata$cell_line <- interesting_cl[which(plot_interesting %in%
                                             metadata$cell_line)] 




## CCLE exp
### TPM
#### Ensembl (only for PCA basal/classical)
ccle_TPM_ens <- readRDS("02_Output/CCLE_ensembl_TPM.RDS")
ccle_TPM_ens <- ccle_TPM_ens[, interesting_cl]
ccle_TPM_ens.non0 <- ccle_TPM_ens[rowSums(ccle_TPM_ens > 0) != 0, ]

#### Gene Symbol
ccle_TPM_sym <- readRDS("02_Output/CCLE_symbols_TPM.RDS")
ccle_TPM_sym <- ccle_TPM_sym[, interesting_cl]
ccle_TPM_sym.non0 <- ccle_TPM_sym[rowSums(ccle_TPM_sym > 0) != 0, ]

### TMM
ccle_counts <- readRDS("02_Output/CCLE_symbols_counts.RDS") # gene symbols, read counts
ccle_counts <- ccle_counts[, interesting_cl]

#### normalization
ccle_counts_dge <- DGEList(ccle_counts)
ccle_TMM_dge <- calcNormFactors(ccle_counts_dge,
  method = "TMM"
)

ccle_TMM <- cpm(ccle_TMM_dge)


### CHOOSE NORM
ccle_sym <- ccle_TMM

# ccle_sym <- ccle_TPM_sym
ccle_ens <- ccle_TPM_ens


# Normalization check (in ens as is the rawest)
ccle_sym.non0 <- ccle_sym[rowSums(ccle_sym > 0) != 0, ]
ccle_checklog <- ccle_sym.non0 + 10**(-2)
ccle_checklog.m <- melt(ccle_checklog,
  measure.vars = 1:ncol(ccle_checklog),
  variable.name = "cell_line"
)

ggplot(ccle_checklog.m, aes(x = value, y = cell_line)) +
  geom_violin() +
  geom_boxplot(width = 0.05) +
  # scale_x_continuous(trans = 'log10') +
  scale_x_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  labs(title = "Normalization check", subtitle = "Non all-zero genes") +
  xlab("log10(TPM + 0.001)") +
  ylab("") +
  theme_classic()

boxplot(ccle_sym, xaxt = "n")
boxplot(ccle_sym.non0, xaxt = "n")
boxplot(ccle_checklog, xaxt = "n", log = "y")

# Characterization
## Serine type
type_serine <- c(rep("Dependent", 6), NA, rep("Independent", 8))
names(type_serine) <- interesting_cl
## Collison
type_collison <- c(NA, rep("QM", 2), "Classical", rep("QM", 5), "Classical", rep(NA, 2), rep("QM", 3))
names(type_collison) <- interesting_cl

## Enzymes
enzymes_ens <- c("ENSG00000092621", "ENSG00000160200", "ENSG00000135069") # PHGDH, CBS, PSAT1
enzymes_sym <- c("PHGDH", "CBS", "PSAT1")

ccle_sym.enz <- ccle_sym[enzymes_sym, ]
ccle_sym.enz <- scale(t(ccle_sym.enz), center = TRUE, scale = TRUE)
ccle_sym.enz <- t(ccle_sym.enz)

## ssGSEA
moffit <- read.csv("01_Input/moffit2015sign.csv", row.names = 1)
moffit <- moffit[!(rownames(moffit) %in% c("CTSL2", "LOC400573")), ]

moffit[moffit$Basal == 1, "Group"] <- "Basal"
moffit[moffit$Classical == 1, "Group"] <- "Classical"

gene_sets <- list(
  "Classical" = row.names(moffit[moffit[, "Classical"] == 1, ]),
  "Basal" = row.names(moffit[moffit[, "Basal"] == 1, ])
)

type_moffit <- gsva(as.matrix(ccle_sym.non0), gene_sets, method = "ssgsea")
type_moffit.plot <- type_moffit[, order(as.numeric(as.character(type_moffit["Basal", ])))]
Heatmap(type_moffit.plot,
  cluster_columns = F
)

## PAMG
library(pdacmolgrad)
# devtools::install_github("RemyNicolle/pdacmolgrad")
type_pamg <- projectMolGrad(newexp = ccle_sym, geneSymbols = row.names(ccle_sym))
type_pamg.log <- projectMolGrad(newexp = log2(ccle_sym + 0.01), geneSymbols = row.names(ccle_sym))

metadata["PAMG"] <- type_pamg["PDX"]

# Basal/Classical importance check
cl.pca <- prcomp(ccle_ens, scale. = TRUE, center = TRUE) # ccle_ens bcos its the most unprocessed
fviz_eig(cl.pca,
  barfill = "lightgrey",
  barcolor = "black", ggtheme = theme_void()
)

## moffit
cl.pca.col <- apply(type_moffit, 2, function(x) which.max(x))

plot(PC1 ~ PC2, data = cl.pca$rotation, col = c("dodgerblue", "chocolate1")[cl.pca.col])
legend(x = "topright", legend = rownames(type_moffit), col = c("dodgerblue", "chocolate1"), pch = 1)

### statistics
pca <- as.data.frame(cl.pca$rotation)
pca$Classical <- type_moffit["Classical", ]
pca$Basal <- type_moffit["Basal", ]

pairs(pca[, c("Basal", "Classical", paste("PC", 1:10, sep = ""))],
  panel = panel.smooth
)
pca.cor <- rcorr(as.matrix(pca), type = "pearson")
corrplot(pca.cor$r,
  order = "hclust", p.mat = pca.cor$P,
  sig.level = 0.05, method = "color", insig = "label_sig",
  title = "Pearson, IC 95%", mar = c(0, 0, 1, 0)
) # http://stackoverflow.com/a/14754408/54964



## collison
cl.pca.col <- type_collison
cl.pca.col[cl.pca.col == "Classical"] <- 1
cl.pca.col[cl.pca.col == "QM"] <- 2
cl.pca.col[is.na(cl.pca.col)] <- 3
cl.pca.col <- as.numeric(cl.pca.col)

plot(PC1 ~ PC2, data = cl.pca$rotation, col = c("dodgerblue", "chocolate1", "grey")[cl.pca.col])
legend(x = "topright", legend = c("Classical", "QM", "NA"), col = c("dodgerblue", "chocolate1", "grey"), pch = 1)

### statistics
pca <- as.data.frame(cl.pca$rotation)
pca$type_collison <- type_collison

pca.m <- melt(pca, measure.vars = paste("PC", 1:15, sep = ""))
ggplot(pca.m, aes(x = variable, y = value, fill = type_collison)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Classical" = "dodgerblue", "QM" = "chocolate1")) +
  theme_classic()

t.test(PC1 ~ type_collison, data = pca)
t.test(PC15 ~ type_collison, data = pca)

## PAMG    # This dont work! is unfisnished and need adaptation from the first plot()
type_pamg <- type_pamg

### statistics
pca <- as.data.frame(cl.pca$rotation)
pca$type_pamg <- type_pamg

# plot(PC1~PC2, data = cl.pca$rotation, col = c("dodgerblue", "chocolate1")[cl.pca.col])
# legend(x="topright", legend = rownames(type_moffit), col=c("dodgerblue", "chocolate1"), pch=1)
#
#
# pairs(pca[,c("type_pamg.PDX", paste("PC", 1:10, sep = ""))],
#       panel=panel.smooth)
# pca.cor <- rcorr(as.matrix(pca), type = "pearson")
# corrplot(pca.cor$r, order="hclust",p.mat = pca.cor$P,
#          sig.level = 0.05, method="color", insig = "label_sig",
#          title = "Pearson, IC 95%", mar=c(0,0,1,0))# http://stackoverflow.com/a/14754408/54964


# Cytokines
cytokines_IC3 <- read_tsv("01_Input/ica3Cytokines.tsv", col_names = TRUE)
ccle_sym.citokines <- ccle_sym[cytokines_IC3$gene, ]


# Heatmap
## Data of heatmap perse
### Moffit 50 genes
# ccle_sym.transclass <- ccle_sym[rownames(moffit),]
# ccle_sym.transclass <- drop_na(ccle_sym.transclass)

### Cytokines
ccle_sym.transclass <- ccle_sym.citokines


ccle_sym.transclass <- scale(t(ccle_sym.transclass), center = TRUE, scale = TRUE)
ccle_sym.transclass <- t(ccle_sym.transclass)



### Column annotation

col_c <- colorRamp2(c(min(type_moffit["Classical", ]), max(type_moffit["Classical", ])), c("white", "dodgerblue"))
col_b <- colorRamp2(c(min(type_moffit["Basal", ]), max(type_moffit["Basal", ])), c("white", "chocolate1"))
col_collison <- c("NA" = "white", "QM" = "orange", "Classical" = "deepskyblue1")
col_serine <- c("NA" = "white", "Dependent" = "firebrick", "Independent" = "seagreen4")
col_PHGDH <- colorRamp2(c(min(ccle_sym.enz["PHGDH", ]), max(ccle_sym.enz["PHGDH", ])), c("white", "grey"))
col_CBS <- colorRamp2(c(min(ccle_sym.enz["CBS", ]), max(ccle_sym.enz["CBS", ])), c("white", "grey"))
col_PSAT1 <- colorRamp2(c(min(ccle_sym.enz["PSAT1", ]), max(ccle_sym.enz["PSAT1", ])), c("white", "grey"))

checkorder.c <- list(
  colnames(ccle_sym.transclass),
  names(type_serine),
  colnames(ccle_sym.enz),
  colnames(type_moffit),
  names(type_collison)
)
print("All must be TRUE")
lapply(
  seq_along(checkorder.c),
  function(n) all(checkorder.c[[n]] == unlist(checkorder.c[-n]))
)

column_ha <- HeatmapAnnotation(
  Serine = type_serine,
  PHGDH = ccle_sym.enz["PHGDH", ],
  CBS = ccle_sym.enz["CBS", ],
  PSAT1 = ccle_sym.enz["PSAT1", ],
  Classical_enr = type_moffit["Classical", ],
  Basal_enr = type_moffit["Basal", ],
  Collison = type_collison,
  col = list(
    Classical_enr = col_c,
    Basal_enr = col_b,
    Collison = col_collison,
    Serine = col_serine,
    PHGDH = col_PHGDH,
    CBS = col_CBS,
    PSAT1 = col_PSAT1
  )
)

### Row annotation
col_moffit_g <- c("Basal" = "chocolate1", "Classical" = "dodgerblue")

checkorder.r <- list(
  rownames(ccle_sym.transclass),
  rownames(moffit)
)
print("All must be TRUE")
lapply(
  seq_along(checkorder.r),
  function(n) all(checkorder.r[[n]] == unlist(checkorder.r[-n]))
)

# row_ha <- rowAnnotation(Moffit_geneset = moffit$Group,
#                            col = list(Moffit_geneset = col_moffit_g))

# ! CAREFUL the annotation is based in position not in labels
colnames(ccle_sym.transclass) <- plot_interesting
Heatmap(ccle_sym.transclass,
  top_annotation = column_ha, # left_annotation = row_ha, #This one is only for when ccle_sym.transclass is Moffit
  cluster_columns = F
)


# ICA analysis
icas_list <- readRDS(file = "02_Output/Bcl_bestfitICs.RDS")
nc_range <- 2:15
corr_sign <- list()
for (nc in nc_range) {
  nc_str <- paste("nc", nc, sep ="")
  nc_ica <- icas_list$A[[nc_str]]
  
  rownames(nc_ica) <- icas_list$samples

  corr_values <- merge(nc_ica, metadata, by.x=0,by.y="cell_line")
  row.names(corr_values) <- corr_values[,"Row.names"]
  corr_values <- subset(corr_values, select = -c(Row.names))
  continuous_var <- c(colnames(nc_ica), "age", "PAMG")

  ## Continuous
  corr_cont <- corr_values[,continuous_var]
  res2 <- rcorr(x=as.matrix(corr_cont))
  corrplot(res2$r, type = "upper", order="original", method = "circle", 
           p.mat = res2$P, sig.level = 0.05, insig = "blank", is.corr = T,
           main = nc_str, mar=c(0,0,1,0))# http://stackoverflow.com/a/14754408/54964
  
  ## Discrete
  corr_disc <- corr_values
  #corr_disc[continuous_var %ni% colnames(nc_ica)] <- NULL
  corr_disc.m <- melt(corr_disc, measure.vars = paste("ic", 1:nc, sep = ""))
  
  
  base_gg <- ggplot(data = corr_disc.m, aes(x=variable, y=value))
  
  ### Serine dependency
  print(base_gg + geom_boxplot(aes(fill = Ser)) + 
    theme_classic() +
    scale_fill_manual(values=c("#F70000", "#1DAB24", "#8B9495")) +
    ggtitle(nc_str, subtitle = "Serine dependency") +
    ylab("weight") +
    xlab("component"))
  
  ttest_data <- corr_disc[!is.na(corr_disc$Ser),]
  ttest_pvals <- data.frame()
  for (col in 1:ncol(nc_ica)){
    ttest_pvals[colnames(nc_ica)[col],"pval"] <- t.test(ttest_data[,col] ~ ttest_data$Ser)$p.value
  }
  corr_sign[nc_str] <- ttest_pvals
  
  ### gender
  print(base_gg + geom_boxplot(aes(fill = gender)) + 
    theme_classic() +
    scale_fill_manual(values=c("#F0BBD1", "#38A1E7")) +
    ggtitle(nc_str, subtitle = "Gender") +
    ylab("weight") +
    xlab("component"))
  
  ### ethnicity
  print(base_gg + geom_boxplot(aes(fill = ethnicity)) + 
    theme_classic() +
    scale_fill_brewer(palette="Dark2") +
    ggtitle(nc_str, subtitle = "Ethnicity ") +
    ylab("weight") +
    xlab("component"))
  
  ### tissue
  print(base_gg + geom_boxplot(aes(fill = tissue)) + 
    theme_classic() +
    scale_fill_brewer(palette="Set1") +
    ggtitle(nc_str, subtitle = "Tissue ") +
    ylab("weight") +
    xlab("component"))
}
  
## check most significant component
for (n in names(corr_sign)){
  print(n)
  print(sapply(corr_sign[n], min))
}








# Statistics
## Enzymes ~ Serine dependency
assay.data <- t(ccle_sym[enzymes_sym, ])
assay.data <- bind_cols(assay.data, as.data.frame(type_serine), .id = "")
assay.data.m <- melt(assay.data)

t.test(PHGDH ~ type_serine, data = assay.data)
t.test(CBS ~ type_serine, data = assay.data)
t.test(PSAT1 ~ type_serine, data = assay.data)

ggplot(drop_na(assay.data.m), aes(x = variable, y = value, fill = type_serine)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("Dependent" = "firebrick", "Independent" = "seagreen4")) +
  theme_classic()

## Cytokines ~ Serine dependency
assay.data <- t(ccle_sym[cytokines_IC3$gene, ])
assay.data <- bind_cols(assay.data, as.data.frame(type_serine), .id = "")
assay.data.m <- melt(assay.data)

t.test(LTB ~ type_serine, data = assay.data) # try of one cytokine that seemed important


ggplot(drop_na(assay.data.m), aes(x = variable, y = value, fill = type_serine)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("Dependent" = "firebrick", "Independent" = "seagreen4")) +
  theme_classic()

# Chemo resistance
chemo <- read_csv("data/cell_celllines/primary-screen-replicate-collapsed-logfold-change.csv")
treat_meta <- read_csv("data/cell_celllines/primary-screen-replicate-collapsed-treatment-info.csv")
cellline_meta <- read_csv("data/cell_celllines/primary-screen-cell-line-info.csv")

chemo <- merge(chemo, cellline_meta[, c("row_name", "ccle_name")],
  by.x = "X1", by.y = "row_name"
)
chemo <- as.data.frame(chemo)

assay.chemo <- chemo[chemo$ccle_name %in% interesting_cl, ]
rownames(assay.chemo) <- assay.chemo$ccle_name
assay.chemo <- subset(assay.chemo, select = -c(ccle_name, X1))
assay.chemo <- t(assay.chemo)
row.nans <- rowSums(is.na(assay.chemo))
assay.chemo.nans <- assay.chemo[row.nans > 0, ]
assay.chemo.nonas <- assay.chemo[row.nans == 0, ]

## Most variable response drugs
abdv <- apply(assay.chemo.nonas, 1, function(x) {
  sum(
    abs(
      x - mean(x)
    )
  ) / length(x)
})

assay.chemo.top <- assay.chemo.nonas[order(abdv, decreasing = TRUE)[1:1000], ] # top K most variant drugs

# Heatmap
Heatmap(assay.chemo.top)
