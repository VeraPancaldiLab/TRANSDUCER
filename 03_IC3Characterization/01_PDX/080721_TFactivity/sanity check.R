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
load("../Remy_processed_data/all_RNA/Tumeur/Hcpmallrna.RData")

Hcpmallrna # This must be quantile normalized data, probably log2
Hcpmallrna <- Hcpmallrna[,sort(colnames(Hcpmallrna))]
boxplot(Hcpmallrna)

## to gene-IDs
Hcpmallrna_genenames <- Hcpmallrna
ensbl_gene <- read_tsv("01_Input/GRCh38.p13_ensemblvsgenename.txt")
translate <- deframe(ensbl_gene)

Hcpmallrna %>% rownames(.) %>% translate[.] -> rownames(Hcpmallrna_genenames)
Hcpmallrna_genenames <- Hcpmallrna_genenames[!is.na(rownames(Hcpmallrna_genenames)),]

# IC3 weights (choose)
## continuous
ic3 <- read_csv2("../Remy_processed_data/samplesIC3.csv") %>%
  as.data.frame(.) %>%
  column_to_rownames(., var = "CITID") %>%
  .[colnames(Hcpmallrna),"ICA3SampleWeight", drop = FALSE]
analysis <- "continuous_IC3"

### Sanity Check Basal/classical TFs (COMME)
library(pdacmolgrad)
# devtools::install_github("RemyNicolle/pdacmolgrad")
type_pamg <- projectMolGrad(newexp = Hcpmallrna_genenames,
                            geneSymbols = row.names(Hcpmallrna_genenames)) %>% arrange(PDX)
sanitycheck_samples <- c(rownames(type_pamg)[1:4],
                         rownames(type_pamg)[-1:-(nrow(type_pamg)-4)])

Hcpmallrna_genenames <- Hcpmallrna_genenames[,sanitycheck_samples]
ic3 <- ic3[sanitycheck_samples,, drop =F]
type_pamg <- type_pamg[sanitycheck_samples,]

# TF activity
## load Dorothea Regulons
### PANcancer
data(dorothea_hs_pancancer, package = "dorothea")
regulons <- dorothea_hs_pancancer %>%
  dplyr::filter(confidence %in% c("A", "B","C"))
vipertype = "PANcancer"

## VIPER
minsize = 5
ges.filter = FALSE

tf_activities_stat <- dorothea::run_viper(Hcpmallrna_genenames, regulons,
                                          options =  list(minsize = minsize, eset.filter = ges.filter, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

## Find 40 TF most correlated with IC3 weights
corr_ic3tf <- rcorr(x = t(tf_activities_stat), y = data.matrix(type_pamg[,"PDX",drop=F]), type = "spearman")
ic3_tfs <- corr_ic3tf$r[,"PDX"] %>%
  abs(.) %>% sort(., decreasing = TRUE) %>% .[2:40] %>% names(.) # exclude first one as its IC3SampleWeight


## Heatmap
ic3_sorted <- ic3[order(ic3$ICA3SampleWeight),,drop = F]
pamg_sorted <- type_pamg[rownames(ic3_sorted), "PDX", drop = F]
tf_heatmap <- tf_activities_stat[ic3_tfs,rownames(ic3_sorted)]

### mean center and scale expression data
# tf_heatmap <- scale(t(tf_heatmap), center = TRUE, scale = TRUE)
# tf_heatmap <- t(tf_heatmap)

### annotations
col_IC3SampleWeight <- colorRamp2(c(min(ic3_sorted), max(ic3_sorted)), c("white", "green"))
col_PAMG <- colorRamp2(c(min(pamg_sorted),median(pamg_sorted$PDX), max(pamg_sorted)), c("orange", "white", "dodgerblue"))

checkorder.c <- list(
  colnames(tf_heatmap),
  rownames(ic3_sorted),
  rownames(pamg_sorted)
)
print("All must be TRUE")
lapply(
  seq_along(checkorder.c),
  function(n) all(checkorder.c[[n]] == unlist(checkorder.c[-n]))
)

column_ha <- HeatmapAnnotation(
  IC3SampleWeight = ic3_sorted$ICA3SampleWeight,
  PAMG = pamg_sorted$PDX,
  col = list(
    IC3SampleWeight = col_IC3SampleWeight,
    PAMG = col_PAMG
    )
)

Heatmap(tf_heatmap,
        top_annotation = column_ha,
        cluster_rows = T,
        cluster_columns = T, 
        column_title = paste("TF activities most correlated with Basal/Classical (sanitycheck):", vipertype, analysis)
)



# msviper
df_4viper <- Hcpmallrna_genenames %>% as_tibble(rownames = "rowname") %>%
  group_by(rowname) %>% summarise_all(mean) %>% column_to_rownames() %>% data.matrix()

##  Basal/classical dependance discrete groups
refg <- "Basal"
expg <- "Classical"
groups <- ifelse(type_pamg$PDX < 0, refg, expg)

### eset creation
phenotype <- data.frame(Ser_dependance = factor(groups))
rownames(phenotype) <- colnames(Hcpmallrna_genenames)
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
