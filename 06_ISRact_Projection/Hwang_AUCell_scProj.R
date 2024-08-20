suppressPackageStartupMessages({
  library(AUCell)
  #library(Biobase)
  library(GSEABase)
  library(data.table)
  library(DT)
  library(NMF)
  library(plotly)
  library(GEOquery)
  library(Matrix)
  # library(doMC);library(doRNG) # Loaded by AUCell, to avoid messages. Not available in windows?
  library(Seurat)
  library(future)
  library(tidyverse)
  library(biomaRt)
  library(readxl)
  library(UpSetR)
})
################################################################################

setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")

# Load sparse expression matrix
expression_matrix <- ReadMtx(
  mtx = "data/Hwang_Nature_2022/GEO/matrix.mtx",
  features = "data/Hwang_Nature_2022/GEO/features.tsv",
  cells = "data/Hwang_Nature_2022/GEO/barcodes.tsv",
  skip.cell = 0
)
## sample to speed up the process
set.seed(333)
expression_matrix <- expression_matrix[sample(rownames(expression_matrix), 5000),]

##clean
gc()

# Build gene sets
## Hwang Signatures
signatures_sheet <- read_excel("data/Hwang_Nature_2022/41588_2022_1134_MOESM4_ESM.xlsx", sheet = 2, skip = 3) %>%
  dplyr::rename(Gene_Rank = `...1`)
malignant_state <- signatures_sheet[1:8]
malignant_lineage <- signatures_sheet[c(1,9:15)]
CAF_lineage <- signatures_sheet[c(1,16:19)]
Hwang_Signatures <- list(malignant_state = malignant_state,
                         malignant_lineage = malignant_lineage,
                         CAF_lineage = CAF_lineage)

### Parse the signatures automatically
HwangGeneSets <- c()
for (l in seq_along(Hwang_Signatures)){
  for(s in seq_along(Hwang_Signatures[[l]])){
    sign = Hwang_Signatures[[l]][s]
    if(names(sign) == "Gene_Rank") next # Skip the lines that contain the ranking
    #print(paste0(names(Hwang_Signatures[l]),"-", names(sign))) #name of gene set
    #print(Hwang_Signatures[[l]][[s]]) # content of gene set
    HwangGeneSets <- c(HwangGeneSets,
                       GeneSet(Hwang_Signatures[[l]][[s]],
                               setName= paste0(names(Hwang_Signatures[l]),
                                               "-", 
                                               names(sign),
                                               "(", nrow(sign),")")))
  }
}


## ISRact
pca_pdx <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS")
ISRact_contributions <- sort(pca_pdx$rotation[,"PC1"])
### translate to Gene ID
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)#listAttributes(ensembl75, page="feature_page")

annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = ensembl75)

translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

ISRact_info <- tibble(EnsemblID = names(ISRact_contributions),
                      gene_name = translate[names(ISRact_contributions)],
                      PC1_score = ISRact_contributions,
                      class = if_else(ISRact_contributions > 0, "ISRac_high", "ISRac_low"))

ggplot(ISRact_info, aes(x = PC1_score)) +
  geom_density() +
  geom_rug(aes(color = class)) +
  theme_bw()

### extract genes
ISRact_high <- dplyr::filter(ISRact_info, class == "ISRac_high") %>%
  dplyr::select("gene_name") %>% 
  deframe()

ISRact_low <- dplyr::filter(ISRact_info, class == "ISRac_low") %>%
  dplyr::select("gene_name") %>% 
  deframe()

### Build Gene set objects
ISRactGeneSets <- c(
  GeneSet(ISRact_high, setName=paste0("ISRact_high(", length(ISRact_high),")")),
  GeneSet(ISRact_low, setName=paste0("ISRact_low(", length(ISRact_low),")")))

geneSets <- GeneSetCollection(c(ISRactGeneSets, HwangGeneSets))

### Quick comparison of gene composition
upset(fromList(geneIds(geneSets)[grep("ISRact|malignant_state",names(geneSets))]), nsets = 16)
upset(fromList(geneIds(geneSets)[grep("ISRact|malignant_lineage",names(geneSets))]), nsets = 16)
upset(fromList(geneIds(geneSets)[grep("ISRact|CAF",names(geneSets))]), nsets = 16)

# AUCell
## Calculate ranking
cells_rankings <- AUCell_buildRankings(expression_matrix, plotStats=TRUE)

## calculate AUC score
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

## Results
### Threshold exploration
set.seed(333)
par(mfrow=c(3,5)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 


### Celltpe asignment exploration
cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)

assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
assignmentMat[,1:2]

#### Plot
set.seed(123)
miniAssigMat <- assignmentMat[,sample(1:ncol(assignmentMat),500)]
library(NMF)
aheatmap(miniAssigMat, scale="none", color="black", legend=FALSE)

## Export for further analysis in the context of the annotation
write_rds(cells_AUC, "results/scRNAseq_proj/Hwang_AUCell_results.RDS")
