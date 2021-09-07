setwd("/home/jacobo/Documents/03_Sauyeun_paper/01_PDX/270521_Stroma_IFNG/")
library(preprocessCore)
library(ggrepel)
load("Rstudio/TranslatomeDataForShiny.RData")
load("../Remy_processed_data/all_RNA/Tumeur/Hcpmallrna.RData")
Hcpmallrna_CHECKED <-  readRDS("../Hcpmallrna.rds") # Remy sent it to me explicitly to check


# same dimensions
tot <- as.data.frame(tot)[,sort(colnames(tot))]
Hcpmallrna <- as.data.frame(Hcpmallrna)[,sort(colnames(Hcpmallrna))]
Hcpmallrna_CHECKED <- as.data.frame(Hcpmallrna_CHECKED)[,sort(colnames(Hcpmallrna_CHECKED))]




# IC3 vector
ic3 <- read_csv2("../Remy_processed_data/samplesIC3.csv") %>%
  as.data.frame(.) %>%
  column_to_rownames(., var = "CITID") %>%
  .[colnames(Hcpmallrna),]


ic3_CHECKED <-  as.data.frame(readRDS("../ICAtranslation.rds"))


##ic3 rownames
print("mine then checked (rownames)")
print(row.names(ic3[order(ic3$ICA3SampleWeight),]))
print(row.names(ic3_CHECKED[order(ic3_CHECKED$ICA3),]))
all(ic3[order(ic3$ICA3SampleWeight),"CITID"] == row.names(ic3_CHECKED[order(ic3_CHECKED$ICA3),]))

##ic3 values
print("mine then checked (values)")
print(ic3[order(ic3$ICA3SampleWeight),"ICA3SampleWeight"])
print(ic3_CHECKED[order(ic3_CHECKED$ICA3),"ICA3"])
all.equal(ic3[order(ic3$ICA3SampleWeight),"ICA3SampleWeight"], ic3_CHECKED[order(ic3_CHECKED$ICA3),"ICA3"])

# The IL18 problematic
gene_plots_Hcp <- cbind(t(Hcpmallrna), as.data.frame(ic3))
gene_plots_Hcp <- gene_plots_Hcp[,!duplicated(colnames(gene_plots_Hcp))] # needed to plot

gene_plots_Hcp_CHECKED <- cbind(t(Hcpmallrna_CHECKED), as.data.frame(ic3))
gene_plots_Hcp_CHECKED <- gene_plots_Hcp_CHECKED[,!duplicated(colnames(gene_plots_Hcp_CHECKED))]


gene_plots_tot <- cbind(t(tot), as.data.frame(ic3))
gene_plots_tot <- gene_plots_tot[,!duplicated(colnames(gene_plots_tot))] # needed to plot

## Hcpmallrna
ggplot(gene_plots_Hcp, aes(x=ICA3SampleWeight, y=ENSG00000150782, label = rownames(gene_plots_Hcp))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="HCP")

## Hcpmallrna_CHECKED
ggplot(gene_plots_Hcp_CHECKED, aes(x=ICA3SampleWeight, y=ENSG00000150782, label = rownames(gene_plots_Hcp_CHECKED))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="HCP_CHECKED")

## Shiny
ggplot(gene_plots_tot, aes(x=ICA3SampleWeight, y=ENSG00000150782, label = rownames(gene_plots_tot))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="shiny")

# Expression
# Choose to compare
d1 <- Hcpmallrna_CHECKED
d1_l <- "Hcpmallrna_CHECKED"
d2 <- tot
d2_l <- "shiny"

## Distribution
boxplot(d1, main = d1_l)
boxplot(d2, main = d2_l)


colnames(d2) <- paste(colnames(d2), d2_l, sep="_")
colnames(d1) <- paste(colnames(d1), d1_l, sep="_")

# Correlations
# ## filtering half most variable genes
# half_most_variable_genes = function(d){
#   abdv <- apply(d, 1, function(x) {
#     sum(
#       abs(
#         x - mean(x)
#       )
#     ) / length(x)
# })
# 
# 
# median_abdv <- median(abdv) # median is the half
# d <- d[abdv > median_abdv, ]
# return(d)
# }
# 
# d2_f <- data.matrix(half_most_variable_genes(d2))
# d1_f <- data.matrix(half_most_variable_genes(d1))

## Stablishing a Z-score
# d2_f <- scale(d2)
# d1_f <- scale(d1)

## Nothing
d2_f <- data.matrix(d2)
d1_f <- data.matrix(d1)


## correlation perse
correlation <- rcorr(d2_f, d1_f, type = "spearman")
correlation$r <- correlation$r[colnames(d1_f), colnames(d2_f)]
correlation$P <- correlation$P[colnames(d1_f), colnames(d2_f)]

corrplot(correlation$r, order="original", , is.corr = FALSE, type = "full", p.mat = correlation$P,
         sig.level = 0.05, method="color", insig = "label_sig")# http://stackoverflow.com/a/

