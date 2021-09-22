library(tidyverse)
library(Hmisc)
library(corrplot)
library(reshape)
library(MultiAssayExperiment) #BiocManager::install("MultiAssayExperiment")
library(lmtest)
library(ggridges)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/04_Deconvolution_TE/03_Sigantures/cell_types/fibroblast/PJ2003085/090921_FastqProcessing")

#load data and metadata
read_tsv("02_Output/rawcounts.tsv") %>% 
  column_to_rownames("EnsemblID") -> counts

read_tsv("02_Output/filteredcounts.tsv") %>%
  column_to_rownames("EnsemblID") -> counts_fil

read_tsv("02_Output/tmmcounts.tsv") %>%
  column_to_rownames("EnsemblID") -> counts_tmm 

read_tsv("01_Input/metadata.tsv") -> sample_info
# boxplots of each
## original
boxplot(log2(1+counts), main="log(raw)")
boxplot(log2(1+counts_fil), main="log(filtered)")
boxplot(log2(1+counts_tmm), main="log(TMM)")

# PCA
counts_pca <- prcomp(counts)
counts_fil_pca <- prcomp(counts_fil)
counts_tmm_pca <- prcomp(counts_tmm)


# Analysis
counts_res <- counts_tmm
pca_res <- counts_tmm_pca
norm_res <- "TMM"
## PCA 

var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
all(rownames(pca_res$rotation) == sample_info$sample) %>% stopifnot()
pca_toplot <- cbind(pca_res$rotation, sample_info[-1])

### origin
pca_toplot %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2, color = origin)) + geom_point(size=4) +
  theme_bw(base_size=15) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top") +
  labs(title = norm_res) 

### fraction
pca_toplot %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2, color = fraction)) + geom_point(size=4) +
  theme_bw(base_size=15) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top") +
  labs(title = norm_res)

### Quality
pca_toplot %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2, color = IG_RIN)) + geom_point(size=4) +
  theme_bw(base_size=15) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top") +
  labs(title = norm_res)

pca_toplot %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2, color = Yield)) + geom_point(size=4) +
  theme_bw(base_size=15) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top") +
  labs(title = norm_res)

## Similarity through whole genome correlation
stat = "spearman"
corr <- rcorr(data.matrix(counts_res), type = stat)
corrplot(corr = corr$r,
         p.mat = corr$P,
         is.corr = F, order = "hclust",
         type = "lower",
         main = paste(norm_res, stat, sep = ": "))


## TE analysis (genewise LRM)
### General efficacy plots

counts_res %>% 
  as.data.frame %>%
ggplot(aes(x=`sMAYCL-tot`, y=`sMAYCL-Pol`)) + geom_point(size=4) +
  theme_bw() + 
  theme(legend.position="top") +
  labs(title = norm_res)



### Multiexpression set adaptation
counts_res[,sample_info[sample_info$fraction == "Polysome", "sample", drop=T]] %>% 
  data.matrix() -> polysome_res
colnames(polysome_res) <- sample_info[sample_info$fraction == "Polysome", "origin", drop=T]

counts_res[,sample_info[sample_info$fraction == "Total", "sample", drop=T]] %>% 
  data.matrix() -> total_res
colnames(total_res) <- sample_info[sample_info$fraction == "Total", "origin", drop=T]

sample_info %>% group_by(origin) %>% summarise() %>%
  column_to_rownames("origin") -> colData_res

rownames(colData_res) -> colData_res$origin


mae <- MultiAssayExperiment(experiments = list(Total = total_res,
                                        Polysome = polysome_res),
                     colData = colData_res,
                     metadata = data.frame(normalization = norm_res,
                                           genes = rownames(counts_res)))
### TE calculation
calculaTE <- function(x)
  {
  fit <- lm(mae@ExperimentList$Polysome[x,] ~ mae@ExperimentList$Total[x,])
  residuals <- fit$residuals
  homoscedasticity <- bptest(fit,studentize = TRUE) # Koenker–Bassett test (homoscedasticity is H0)
  normality <- shapiro.test(residuals)
  
  return(c(x, residuals, homoscedasticity$p.value, normality$p.value))
}

lapply(mae@metadata$genes, calculaTE) %>% 
   as.data.frame(row.names = c("EnsemblID", rownames(colData(mae)),
                               "Phomo", "Pnorm"))  %>% t() %>% as_tibble() -> TEs.unf

TEs.unf %>% column_to_rownames("EnsemblID") %>%
  mutate_all(function(x) as.numeric(as.character(x))) -> TEs.unf
   
### Filter out of non valid genes lm and add TEs to Multiassay.
TEs.unf %>% filter(Phomo > 0.05 & Pnorm > 0.05) %>%
  select(!c(Phomo, Pnorm)) %>% data.matrix() -> TEs # Filter based in the pvalue

print(paste(nrow(TEs), " out of ",
            nrow(TEs.unf), " can be kept for TE analysis."))

mae <- c(mae, TE = TEs)

### non negative TE transformation
min_TEs <- min(TEs)
TEs.nneg <- TEs + abs(min_TEs) + 1

all(sort(TEs.nneg, index.return = T)$ix == sort(TEs, index.return = T)$ix) # false because there's some equal TEs that mess out the sorting 

ggplot(as.data.frame(TEs), aes(x = s17T)) + 
  geom_density()

ggplot(as.data.frame(TEs.nneg), aes(x = s17T)) + 
  geom_density()

### mean sample comparisons
vs <- "mean"
genes <- rownames(experiments(mae)$TE)
means.mae <- data.frame(row.names = genes)
means.mae$Total <- apply(experiments(mae)$Total[genes,],1,mean)
means.mae$Polysome <- apply(experiments(mae)$Polysome[genes,],1,mean)
means.mae$TE <- apply(TEs.nneg[genes,],1,mean)

#### calculate FC
l2fc.Total <- apply(mae@ExperimentList$Total[genes,], 2,
                        function(x) x/means.mae$Total) 
l2fc.Polysome <- apply(mae@ExperimentList$Polysome[genes,], 2,
                           function(x) x/means.mae$Polysome)
l2fc.TE <- apply(TEs.nneg[genes,], 2,
                     function(x) x/means.mae$TE)

l2fc.Total = log2(l2fc.Total[rowMins(l2fc.Total)!=0,])
l2fc.Polysome = log2(l2fc.Polysome[rowMins(l2fc.Polysome)!=0,])
l2fc.TE = log2(l2fc.TE[rowMins(l2fc.TE)!=0,])

ggplot(melt(l2fc.Total), aes(x = value, y = X2)) + geom_density_ridges(rel_min_height = 0.00000000000000001) + 
  scale_x_continuous("l2FC") + 
  labs(title = paste("TotalRNA vs.", vs)) +
  xlim(-13,5) +
  theme_classic()

ggplot(melt(l2fc.Polysome), aes(x = value, y = X2)) + geom_density_ridges(rel_min_height = 0.00000000000000001) + 
  scale_x_continuous("l2FC") + 
  labs(title = paste("Polysome vs.", vs)) +
  xlim(-13,5) +
  theme_classic()

ggplot(melt(l2fc.TE), aes(x = value, y = X2)) + geom_density_ridges(rel_min_height = 0.00000000000000001) + 
  scale_x_continuous("l2FC") + 
  labs(title = paste("TEs vs.", vs)) +
  xlim(-13,5) +
  theme_classic()

#### Saving of a single df for indepth analysis
Foldchanges <- merge(l2fc.Total, l2fc.Polysome, by="row.names",
                     all =T, suffixes = c(".Total", ".Polysome")) %>% column_to_rownames("Row.names")
Foldchanges <- merge(Foldchanges, l2fc.TE, by="row.names",
                     all =T, suffixes = c("", ".TE")) %>% column_to_rownames("Row.names")
filename <- paste(paste("02_Output/", vs, sep =""),"foldchanges.tsv", sep = ".")
Foldchanges %>% write_tsv(filename)
