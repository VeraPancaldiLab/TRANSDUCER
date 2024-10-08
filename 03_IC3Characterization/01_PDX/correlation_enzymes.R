library(tidyverse)
library(reshape2)
library(Hmisc)
library(corrplot)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/03_IC3Characterization/01_PDX/")
# Load data
load("Remy_processed_data/all_RNA/Tumeur/Hcpmallrna.RData") # the same as Alexia confirmed as the normalized one 
Hcpmallrna %>% as_tibble(rownames = "EnsemblID") %>%
  column_to_rownames("EnsemblID") -> allrna

load("TranslatomeDataForShiny.RData") # shiny app data (Tumour residuals included)
resi %>% as_tibble(rownames = "EnsemblID") %>%
  column_to_rownames("EnsemblID") -> TEs

read_tsv("Remy_processed_data/samplesIC3_custom.csv") %>%
  dplyr::select(CITID,ICA3SampleWeight) %>% dplyr::rename(SerDep = ICA3SampleWeight) %>%
  column_to_rownames("CITID") -> ic3

# select the names of the genes to check
enzymes <- c("ENSG00000092621", "ENSG00000160200", "ENSG00000135069") # PHGDH, CBS, PSAT1
enzymes_l <- c("PHGDH", "CBS", "PSAT1")
names(enzymes) <- enzymes_l

# Distribution check
boxplot(allrna, xaxt = "n")
boxplot(TEs, xaxt = "n")

# Linearity check
assay.data <- allrna[enzymes,]
row.names(assay.data) <- enzymes_l

assay.data <- merge(t(assay.data), ic3, by="row.names")
rownames(assay.data) <- assay.data$Row.names
assay.data <- subset(assay.data, select = -Row.names)

pairs(assay.data, panel=panel.smooth)

# Correlation test
res_s <- rcorr(as.matrix(assay.data), type = "spearman")
corrplot(res_s$r, order="hclust",p.mat = res_s$P,
               sig.level = 0.05, method="color", insig = "label_sig",
               title = "Spearman", mar=c(0,0,1,0))# http://stackoverflow.com/a/14754408/54964

res_p <- rcorr(as.matrix(assay.data), type = "pearson")
corrplot(res_p$r, order="hclust",p.mat = res_p$P,
         sig.level = 0.05, method="color", insig = "label_sig",
         title = "Pearson", mar=c(0,0,1,0))# http://stackoverflow.com/a/14754408/54964

# Binarization in ISR+ and ISR-
bin.enzyme = "CBS"
#-
ic3_order <- order(ic3[,"SerDep"])
ic3_plus <- tail(rownames(ic3)[ic3_order],5)
ic3_minus <- head(rownames(ic3)[ic3_order],5)


bin.data <- allrna[enzymes[bin.enzyme],c(ic3_plus, ic3_minus)]
bin.data["ISR", ic3_plus] <- "ISR+"
bin.data["ISR", ic3_minus] <- "ISR-"
bin.data <- as.data.frame(t(bin.data))
colnames(bin.data) <- c(bin.enzyme, "ISR")
bin.data[,bin.enzyme]=as.numeric(bin.data[,bin.enzyme])
bin.data[,"ISR"]=as.factor(bin.data[,"ISR"])
#-                                                  
ggplot(bin.data, aes(x=ISR, y=CBS)) +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.1) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 1)) +
  ylab(label = "log2( CBS cpm + 1 )") +
  xlab(label = "") +
  theme_classic()


t.test(CBS ~ ISR, data = bin.data)


# The PHGDH correlation outlier
## zoom in PHGDH to find PDAC024T as the outlier
scatter.smooth(assay.data$SerDep, assay.data$PHGDH)
text(assay.data$SerDep, assay.data$PHGDH, labels=rownames(assay.data), cex=0.9, font=2)

assay.no024T <- assay.data[row.names(assay.data)!=c("PDAC024T"),]
scatter.smooth(assay.no024T$SerDep, assay.no024T$PHGDH)
text(assay.no024T$SerDep, assay.no024T$PHGDH, labels=rownames(assay.no024T), cex=0.9, font=2)

## correlation
resno024T_s <- rcorr(as.matrix(assay.no024T), type = "spearman")
corrplot(resno024T_s$r, order="hclust",p.mat = resno024T_s$P,
         sig.level = 0.05, method="color", insig = "label_sig",
         title = "Spearman: no O24T", mar=c(0,0,1,0))# http://stackoverflow.com/a/14754408/54964

resno024T_p <- rcorr(as.matrix(assay.no024T), type = "pearson")
corrplot(resno024T_p$r, order="hclust",p.mat = resno024T_p$P,
         sig.level = 0.05, method="color", insig = "label_sig",
         title = "Pearson: no O24T", mar=c(0,0,1,0))# http://stackoverflow.com/a/14754408/54964

