library(tidyverse)
library(edgeR)
################################################################################
# Data loading
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/00_Data/")
load("geneCount_raw_28s.RData")
equivalences <-
  read_tsv("sample_names_equivalences.tsv") %>% dplyr::select("fastq_name", "CITID")

cleanTumour <- c(
  "PDAC014T",
  "PDAC018T",
  "PDAC029T",
  "PDAC013T",
  "PDAC020T",
  "PDAC021T",
  "PDAC003T",
  "PDAC017T",
  "PDAC031T",
  "PDAC009T",
  "PDAC012T",
  "PDAC032T",
  "PDAC028T",
  "PDAC030T",
  "PDAC006T",
  "PDAC027T",
  "PDAC025T",
  "PDAC011T",
  "PDAC016T",
  "PDAC007T",
  "PDAC024T",
  "PDAC008T",
  "PDAC026T",
  "PDAC019T",
  "PDAC015T",
  "PDAC001T",
  "PDAC022T"
) # 27

cleanHost <- c(
  "PDAC018T",
  "PDAC029T",
  "PDAC013T",
  "PDAC020T",
  "PDAC021T",
  "PDAC017T",
  "PDAC031T",
  "PDAC012T",
  "PDAC032T",
  "PDAC028T",
  "PDAC030T",
  "PDAC027T",
  "PDAC025T",
  "PDAC011T",
  "PDAC024T",
  "PDAC008T",
  "PDAC026T",
  "PDAC019T",
  "PDAC015T",
  "PDAC022T"
) # 20

# Separate reads by fraction and organism
## by species (ID)
rownames(counts) %>% str_detect("ENSG") -> Tumor_or_stroma
counts %>%
  as_tibble(rownames = "EnsemblID") %>%
  .[Tumor_or_stroma, ] -> countsTumor_mix
counts %>%
  as_tibble(rownames = "EnsemblID") %>%
  .[!Tumor_or_stroma, ] -> countsHost_mix

## by fraction, and put universal names
colnames(counts) %>%
  str_detect("_F8-9", negate = T) %>%
  colnames(counts)[.] -> names_total

colnames(counts) %>%
  str_detect("_F8-9") %>%
  colnames(counts)[.] -> names_polysome

all(str_remove(names_polysome, "_F8.9") == names_total) %>% stopifnot()

translate <-
  deframe(equivalences) # will use this next to adapt all names

sample_names <- tibble(
  dataT = names_total,
  dataP = names_polysome,
  standard = translate[names_total]
) # OK


## Split Tumour
PDX_count_spliter <- function(dat, conversion) {
  splitted <- list()
  for (fraction in c("dataT", "dataP")) {
    dat %>%
      dplyr::select(c("EnsemblID", conversion[[fraction]])) %>%
      rename_at(conversion[[fraction]], ~ conversion[["standard"]]) %>%
      column_to_rownames("EnsemblID") -> splitted[[fraction]]
  }
  return(splitted)
}

## Split Tumour
PDX_count_spliter(countsTumor_mix, sample_names) -> countsTumor
countsTumor$dataT %>%
  rownames_to_column("EnsemblID") %>%
  write_tsv("Processed_data/rawTumor_Cyt.tsv")
countsTumor$dataP %>%
  rownames_to_column("EnsemblID") %>%
  write_tsv("Processed_data/rawTumor_Pol.tsv")

## Split Stroma
PDX_count_spliter(countsHost_mix, sample_names) -> countsHost
countsHost$dataT %>%
  rownames_to_column("EnsemblID") %>%
  write_tsv("Processed_data/rawHost_Cyt.tsv")
countsHost$dataP %>%
  rownames_to_column("EnsemblID") %>%
  write_tsv("Processed_data/rawHost_Pol.tsv")

## Get ncounts to test latter (just host samples)
ncounts <- tibble(samples = cleanHost)
countsHost$dataT %>% dplyr::select(cleanHost) %>% colSums() -> ncounts$cyt_host_counts
countsHost$dataP %>% dplyr::select(cleanHost) %>% colSums() -> ncounts$pol_host_counts
countsTumor$dataT %>% dplyr::select(cleanHost) %>% colSums() -> ncounts$cyt_tumor_counts
countsTumor$dataP %>% dplyr::select(cleanHost) %>% colSums() -> ncounts$pol_tumor_counts
ncounts %>% write_tsv("Processed_data/ncounts.tsv")

## Remove 0 genes and normalize
Clean_and_Norm <- function(data, clean_samples) {
  tmpdataT <- data$dataT[clean_samples]
  tmpdataP <- data$dataP[clean_samples]

  colnames(tmpdataT) <- paste(clean_samples, "Cyt", sep = "_")
  colnames(tmpdataP) <- paste(clean_samples, "Pol", sep = "_")

  data_tmp <- cbind(tmpdataT, tmpdataP)

  filt_tmp <-
    data_tmp[!(apply(data_tmp, 1, function(x) {
      any(x == 0)
    })), ] # from anota2seqRemoveZeroSamples()
  norm_tmp <-
    limma::voom(edgeR::calcNormFactors(edgeR::DGEList(filt_tmp)))$E ## TMM-log2 from anota2seqNormalize()

  norm_t <- norm_tmp[, colnames(tmpdataT)]
  norm_p <- norm_tmp[, colnames(tmpdataP)]

  colnames(norm_t) <- gsub("_Cyt", "", colnames(norm_t))
  colnames(norm_p) <- gsub("_Pol", "", colnames(norm_p))

  return(list(dataT = norm_t, dataP = norm_p))
}


normTumour <- Clean_and_Norm(
  data = countsTumor,
  clean_samples = cleanTumour
)
normTumour$dataT %>%
  as_tibble(rownames = "EnsemblID") %>%
  write_tsv("Processed_data/normTumor_Cyt.tsv")
normTumour$dataP %>%
  as_tibble(rownames = "EnsemblID") %>%
  write_tsv("Processed_data/normTumor_Pol.tsv")


normHost <- Clean_and_Norm(
  data = countsHost,
  clean_samples = cleanHost
)
normHost$dataT %>%
  as_tibble(rownames = "EnsemblID") %>%
  write_tsv("Processed_data/normHost_Cyt.tsv")
normHost$dataP %>%
  as_tibble(rownames = "EnsemblID") %>%
  write_tsv("Processed_data/normHost_Pol.tsv")


# Boxplots
plot_title <- list("Cytosolic", "Polysome")
for (i in 1:2) {
  countsTumor[[i]] %>%
    as_tibble(rownames = "EnsemblID") %>%
    pivot_longer(-EnsemblID, names_to = "samples") %>%
    ggplot(aes(x = samples, y = value + 1)) +
    geom_violin(fill = "gray") +
    geom_boxplot(width = 0.1, fill = "white") +
    scale_y_continuous(trans = "log2") +
    labs(title = paste(plot_title[[i]], "Tumour Raw", sep = " ")) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    )) -> ct
  print(ct)

  countsHost[[i]] %>%
    as_tibble(rownames = "EnsemblID") %>%
    pivot_longer(-EnsemblID, names_to = "samples") %>%
    ggplot(aes(x = samples, y = value + 1)) +
    geom_violin(fill = "gray") +
    geom_boxplot(width = 0.1, fill = "white") +
    scale_y_continuous(trans = "log2") +
    labs(title = paste(plot_title[[i]], "Host Raw", sep = " ")) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    )) -> ch
  print(ch)

  normTumour[[i]] %>%
    as_tibble(rownames = "EnsemblID") %>%
    pivot_longer(-EnsemblID, names_to = "samples") %>%
    ggplot(aes(x = samples, y = value)) +
    geom_violin(fill = "gray") +
    geom_boxplot(width = 0.1, fill = "white") +
    labs(title = paste(plot_title[[i]], "Tumour TMM-log2", sep = " ")) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    )) -> nt
  print(nt)

  normHost[[i]] %>%
    as_tibble(rownames = "EnsemblID") %>%
    pivot_longer(-EnsemblID, names_to = "samples") %>%
    ggplot(aes(x = samples, y = value)) +
    geom_violin(fill = "gray") +
    geom_boxplot(width = 0.1, fill = "white") +
    labs(title = paste(plot_title[[i]], "Host TMM-log2", sep = " ")) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    )) -> nh
  print(nh)
}
