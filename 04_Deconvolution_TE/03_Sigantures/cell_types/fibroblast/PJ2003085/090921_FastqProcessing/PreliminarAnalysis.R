library(tidyverse)
library(Hmisc)
library(corrplot)
library(reshape)
library(MultiAssayExperiment) #BiocManager::install("MultiAssayExperiment")
library(lmtest)
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


### data reshape
# tmpd <- t(counts_res)
# 
# tmpd.pol <- tmpd[sample_info[sample_info$fraction == "Polysome", "sample", drop=T],]
# rownames(tmpd.pol) <- sample_info[sample_info$fraction == "Polysome", "origin", drop=T]
# colnames(tmpd.pol) <- paste(colnames(tmpd.pol), ".pol", sep = "")
# 
# tmpd.tot <- tmpd[sample_info[sample_info$fraction == "Total", "sample", drop=T],]
# rownames(tmpd.tot) <- sample_info[sample_info$fraction == "Total", "origin", drop=T]
# colnames(tmpd.tot) <- paste(colnames(tmpd.tot), ".tot", sep = "")
# 
# tmpd <- merge(tmpd.tot, tmpd.pol, by="row.names") %>% 
#   column_to_rownames("Row.names")
# 
# a <- "ENSG00000278384.tot"
# b <- "ENSG00000278384.pol"
# 
# fit <- lm(get(a)~get(b), data = tmpd)
# plot(ENSG00000278384.tot~ENSG00000278384.pol, data = tmpd)
# text(tmpd[,b,drop =T], tmpd[,a,drop =T], labels = rownames(tmpd),  pos =1 )
# abline(fit)
# residuals(fit)


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

TEs.unf %>% column_to_rownames("EnsemblID") %>% # Turn the columns into numeruc values
   
   
TEs.unf %>% filter() # Filter based in the pvalue