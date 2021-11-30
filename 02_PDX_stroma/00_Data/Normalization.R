library(tidyverse)
library(edgeR)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/00_Data/")
load("geneCount_raw_28s.RData")

# Separate reads by fraction and organism
## by species (ID)
rownames(counts) %>% str_detect("ENSG") -> Tumor_or_stroma
counts %>% as_tibble(rownames = "EnsemblID") %>% .[Tumor_or_stroma,] -> countsTumor
counts %>% as_tibble(rownames = "EnsemblID") %>% .[!Tumor_or_stroma,] -> countsHost

## by fraction, and put universal names
colnames(counts) %>% str_detect("_F8-9", negate = T) %>% 
  colnames(counts)[.] -> names_total

colnames(counts) %>% str_detect("_F8-9") %>%
  colnames(counts)[.] -> names_polysome

all(str_remove(names_polysome, "_F8.9") == names_total) %>% stopifnot()

equivalences <- read_tsv("sample_names_equivalences.tsv") %>% .[,c("fastq_name", "CITID")]
translate <- deframe(equivalences) # will use this next to adapt all names

sample_names = tibble(total = names_total,
                      polysome = names_polysome,
                      standard = translate[names_total]) #OK

## Split Tumour
countsTumor %>% dplyr::select(c("EnsemblID", sample_names$total)) %>%
  rename_at(sample_names$total, ~ sample_names$standard) -> countsTumor_tot
countsTumor_tot %>% write_tsv("Processed_data/countsTumor_total.tsv")

countsTumor %>% dplyr::select(c("EnsemblID", sample_names$polysome)) %>%
  rename_at(sample_names$polysome, ~ sample_names$standard) -> countsTumor_pol  # this renames with official sample names
countsTumor_pol %>% write_tsv("Processed_data/countsTumor_polysome.tsv") 

## Split Stroma
countsHost %>% dplyr::select(c("EnsemblID", sample_names$total)) %>%
  rename_at(sample_names$total, ~ sample_names$standard) -> countsHost_tot
countsHost_tot %>% write_tsv("Processed_data/countsHost_total.tsv")  

countsHost %>% dplyr::select(c("EnsemblID", sample_names$polysome)) %>%
  rename_at(sample_names$polysome, ~ sample_names$standard) -> countsHost_pol # this renames with official sample names
countsHost_pol %>% write_tsv("Processed_data/countsHost_polysome.tsv")


# Filter (no need, already filtered for non 0)
for (x in list(countsTumor_tot, countsTumor_pol,
           countsHost_tot, countsHost_pol)){
  print(min(rowSums(x!=0) > 0))
}


# TMM norm
countsTumor_tot %>% column_to_rownames("EnsemblID") %>%
  DGEList() %>% calcNormFactors(method = "TMM") %>%
  cpm() -> countsTumor_tot.tmm

countsTumor_pol %>% column_to_rownames("EnsemblID") %>%
  DGEList() %>% calcNormFactors(method = "TMM") %>%
  cpm() -> countsTumor_pol.tmm

countsHost_tot %>% column_to_rownames("EnsemblID") %>%
  DGEList() %>% calcNormFactors(method = "TMM") %>%
  cpm() -> countsHost_tot.tmm

countsHost_pol %>% column_to_rownames("EnsemblID") %>%
  DGEList() %>% calcNormFactors(method = "TMM") %>%
  cpm() -> countsHost_pol.tmm


# Boxplots
raw = list(countsTumor_tot, countsTumor_pol,
           countsHost_tot, countsHost_pol)

tmm = list(countsTumor_tot.tmm, countsTumor_pol.tmm,
                 countsHost_tot.tmm, countsHost_pol.tmm)

plot_title = list("Tumor_tot", "Tumor_pol",
                  "Host_tot", "Host_pol")
for (i in 1:4){
  
raw[[i]] %>% column_to_rownames("EnsemblID") -> raw_i
boxplot(raw_i + 1, log = "y",
          main = paste(plot_title[[i]], "raw", sep = " "), las = 2)
  
boxplot((tmm[[i]] + 1), log = "y",
          main = paste(plot_title[[i]], "tmm", sep = " "),  las = 2)
  }



