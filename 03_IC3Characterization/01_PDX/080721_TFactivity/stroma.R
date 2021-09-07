# Package loading
library(tidyverse)
library(dorothea)
library(viper)
library(circlize)
library(ComplexHeatmap)
library(Hmisc)

setwd("/home/jacobo/Documents/03_Sauyeun_paper/01_PDX/080721_TFactivity/")


# Data loading
## Expression data
load("../Remy_processed_data/all_RNA/Stroma/Mcpmallrna.RData")

Mcpmallrna # This must be quantile normalized data, probably log2
Mcpmallrna <- Mcpmallrna[,sort(colnames(Mcpmallrna))]
boxplot(Mcpmallrna)

## to gene-IDs
Mcpmallrna_genenames <- Mcpmallrna
ensbl_gene <- read_tsv("01_Input/GRCm39_ensemblvsgenename.txt")
translate <- deframe(ensbl_gene)

Mcpmallrna %>% rownames(.) %>% translate[.] -> rownames(Mcpmallrna_genenames)
Mcpmallrna_genenames <- Mcpmallrna_genenames[!is.na(rownames(Mcpmallrna_genenames)),]

# IC3 weights (choose)
## continuous
ic3 <- read_tsv("../Remy_processed_data/samplesIC3_custom.csv") %>%
  as.data.frame(.) %>%
  column_to_rownames(., var = "CITID") %>%
  .[colnames(Mcpmallrna),, drop = FALSE]
analysis <- "continuous_IC3"

## Discrete
# ### 4vs4noPDAC024T
# ic3 <- read_tsv("../020821_IC3_Into_Two_Groups/02_Output/samplesIC3_4vs4noPDAC024T.tsv") %>%
#   as.data.frame(.) %>%
#   column_to_rownames(., var = "rowname") %>% 
#   .[colnames(Mcpmallrna),, drop = FALSE] %>%  drop_na("groups") #sorted by IC3
# 
# analysis <- "4vs4noPDAC024T_IC3"
# # 
# ### 5vs5
# ic3 <- read_tsv("../020821_IC3_Into_Two_Groups/02_Output/samplesIC3_5vs5.tsv") %>%
#   as.data.frame(.) %>%
#   column_to_rownames(., var = "rowname") %>%
#   .[colnames(Mcpmallrna),, drop = FALSE] %>% drop_na("groups") #sorted by IC3
# 
# analysis <- "5vs5"

### common between different kinds
diab <-   ic3[,c("Diabetes"), drop = FALSE]
ic3 <- ic3[,c("ICA3SampleWeight"), drop = FALSE]

Mcpmallrna_genenames <- Mcpmallrna_genenames[,rownames(ic3)]
Mcpmallrna <- Mcpmallrna[,rownames(ic3)]


# PAMG (from tumour.R)
read_tsv("02_Output/tumorPAMG.tsv") %>% column_to_rownames() -> type_pamg

# TF activity
## load Dorothea Regulons
### GTex
# data(dorothea_mm, package = "dorothea")
# regulons <- dorothea_mm %>%
#   dplyr::filter(confidence %in% c("A", "B","C"))
# vipertype = "Gtex"

### PANcancer
data(dorothea_mm_pancancer, package = "dorothea")
regulons <- dorothea_mm_pancancer %>%
  dplyr::filter(confidence %in% c("A", "B","C"))
vipertype = "PANcancer"

## VIPER
minsize = 5
ges.filter = FALSE

tf_activities_stat <- dorothea::run_viper(Mcpmallrna_genenames, regulons,
                                          options =  list(minsize = minsize, eset.filter = ges.filter, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

## Find 40 TF most correlated with IC3 weights
corr_ic3tf <- rcorr(x = t(tf_activities_stat), y = data.matrix(ic3), type = "spearman")
ic3_tfs <- corr_ic3tf$r[,"ICA3SampleWeight", drop = F] %>% .[!(row.names(.) %in% c("ICA3SampleWeight")),] %>%
  abs(.) %>% sort(., decreasing = TRUE)  %>% .[1:40] %>% names(.) # exclude first one as its IC3SampleWeight

### explore the list of correlations
ic3_tfs_explore.r <- corr_ic3tf$r[,"ICA3SampleWeight"]
ic3_tfs_explore.P <- corr_ic3tf$P[,"ICA3SampleWeight"]

print(c("If false, theres something wrong:", 
        all(colnames(tf_activities_stat) == rownames(ic3))) )

toplot <- cbind(t(tf_activities_stat), ic3)
stats <- paste("Spearman: R = ", round(ic3_tfs_explore.r["Pax5"], 2),
               ", pval = ", round(ic3_tfs_explore.P["Pax5"], 4), sep = "")

ggplot(toplot, aes(x=ICA3SampleWeight, y=Pax5, label = rownames(toplot))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Stroma", subtitle = stats)


## Heatmap
ic3_sorted <- ic3[order(ic3$ICA3SampleWeight),,drop = F]
pamg_sorted <- type_pamg[rownames(ic3_sorted), "PDX", drop = F]
diab_stat <- diab[rownames(ic3_sorted),,drop = F]
tf_heatmap <- tf_activities_stat[ic3_tfs,rownames(ic3_sorted)]

### mean center and scale expression data
tf_heatmap <- scale(t(tf_heatmap), center = TRUE, scale = TRUE)
tf_heatmap <- t(tf_heatmap)

### annotations
col_IC3SampleWeight <- colorRamp2(c(min(ic3_sorted), max(ic3_sorted)), c("white", "green"))
col_PAMG <- colorRamp2(c(min(pamg_sorted),median(pamg_sorted$PDX), max(pamg_sorted)), c("orange", "white", "dodgerblue"))
col_diab <- c("NA" = "grey", "0" = "white", "1" = "black")

checkorder.c <- list(
  colnames(tf_heatmap),
  rownames(ic3_sorted),
  rownames(pamg_sorted),
  rownames(col_diab)
)
print("All must be TRUE")
lapply(
  seq_along(checkorder.c),
  function(n) all(checkorder.c[[n]] == unlist(checkorder.c[-n]))
)

column_ha <- HeatmapAnnotation(
  IC3SampleWeight = ic3_sorted$ICA3SampleWeight,
  PAMG = pamg_sorted$PDX,
  Diabetes = diab_stat$Diabetes,
  col = list(
    IC3SampleWeight = col_IC3SampleWeight,
    PAMG = col_PAMG,
    Diabetes = col_diab
    )
)

Heatmap(tf_heatmap,
        top_annotation = column_ha,
        cluster_rows = T,
        cluster_columns = F, 
        column_title = paste("TF activities most correlated with IC3:", vipertype, analysis)
)


## Table
ic3_corrTable <- tibble("Genes" =  ic3_tfs,
                        "R" = corr_ic3tf$r[ic3_tfs, "ICA3SampleWeight"],
                        "pvalue" = corr_ic3tf$P[ic3_tfs, "ICA3SampleWeight"])

write_tsv(ic3_corrTable, paste("02_Output/IC3BestTFActivity_", vipertype, "_", analysis  ,".tsv", sep = ""))







# msviper
df_4viper <- Mcpmallrna_genenames %>% as_tibble(rownames = "rowname") %>%
  group_by(rowname) %>% summarise_all(mean) %>% column_to_rownames() %>% data.matrix()

##  Serine dependance discrete groups
refg <- "Ser_i"
expg <- "Ser_d"
groups <- ifelse(ic3$ICA3SampleWeight < 0, refg, expg)

### eset creation
phenotype <- data.frame(Ser_dependance = factor(groups)) # por eso en Ttest "condition"
rownames(phenotype) <- colnames(Mcpmallrna_genenames)
phenoData <- new("AnnotatedDataFrame", data = phenotype)

dset_viper <- ExpressionSet(assayData = df_4viper, phenoData = phenoData)
dset_viper$sampleID <- factor(colnames(df_4viper))


# get genes significantly associated with each group
signature <- rowTtest(dset_viper, "Ser_dependance", expg, refg)
statistics_signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]



# Generate the null model with bootstrapping (1000 iterations). This are sample permutations
nullmodel <- ttestNull(dset_viper, "Ser_dependance", expg, refg, per = 1000, repos = T, verbose = F)


# Run msviper using the statistics signature, the regulons converted from dorothea table, the null model the minSize of regulon and the ges.filter
mrs <- msviper(ges = statistics_signature, regulon = df2regulon(regulons), nullmodel = nullmodel, minsize = minsize, ges.filter = ges.filter, verbose = F)
summary(mrs)
plot(mrs)

# shadow analysis to penalyze pleyotropic interactions
mrs_s <- shadow(mrs, 10)
summary(mrs_s)
