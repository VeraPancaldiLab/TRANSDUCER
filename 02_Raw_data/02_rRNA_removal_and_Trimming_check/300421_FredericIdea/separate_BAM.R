wd <- "/home/jacobo/Documents/02_TRANSDUCER/02_Raw_data/02_rRNA_removal_and_Trimming_check/300421_FredericIdea"
setwd(wd)
#install.packages("https://bioconductor.org/packages/3.5/bioc/src/contrib/Rsubread_1.26.1.tar.gz", repos=NULL, type="source") # Try with the version of 2017

library(tidyverse)
library(reshape2)
source("SMAP/R/SMAPcount.R")


# Separation of reads
counts = SMAPcount(SMAPBAM = "/home/jacobo/Documents/02_TRANSDUCER/02_Raw_data/02_rRNA_removal_and_Trimming_check/300421_FredericIdea/output_bams",
                   GTF = "combined_reference/combined.gtf")
countsTumor = counts$FCcountsTumor
countsHost = counts$FCcountsHost

colnames(countsTumor) = c("GC50", "GC60")
colnames(countsHost) = c("GC50", "GC60")

# keep only non0s
Tnon0 <-rowSums(countsTumor) != 0
Hnon0 <-rowSums(countsHost) != 0

countsTumor <- countsTumor[Tnon0,]
countsHost <- countsHost[Hnon0,]

# tumor dotplot
countsTumor.m <- melt(countsTumor)
colnames(countsTumor.m) = c("gene", "GC", "reads")


ggplot(data = countsTumor.m, aes(x=GC, y=reads)) +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.05) +
  scale_y_continuous(trans = 'log2') +
  theme_classic() +
  labs(title = "Tumour reads")


# host dotplot
countsHost.m <- melt(countsHost)
colnames(countsHost.m) = c("gene", "GC", "reads")

ggplot(data = countsHost.m, aes(x=GC, y=reads)) +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.05) +
  scale_y_continuous(trans = 'log2') +
  theme_classic() +
  labs(title = "Host reads")

# barplot combined
## absolute
total_counts = data.frame()
total_counts["tumor",c("GC50", "GC60")] = colSums(countsTumor)
total_counts["host",c("GC50", "GC60")] = colSums(countsHost)

total_counts.m <- melt(as.matrix(total_counts))
colnames(total_counts.m) <- c("origin", "GC", "reads")
ggplot(data = total_counts.m, aes(x = GC, y = reads, fill = origin)) +
  geom_bar(stat="identity", position="dodge") +
  theme_classic() +
  labs(title = "count comparison (abs)")

## relative
total_counts.p = sweep(total_counts, 2, colSums(total_counts), FUN = '/')

total_counts.pm <- melt(as.matrix(total_counts.p))
colnames(total_counts.pm) <- c("origin", "GC", "reads")
ggplot(data = total_counts.pm, aes(x = GC, y = reads, fill = origin)) +
  geom_bar(stat="identity", position="dodge") +
  scale_y_continuous(labels=scales::percent) +
  theme_classic() +
  labs(title = "count comparison (%)")
