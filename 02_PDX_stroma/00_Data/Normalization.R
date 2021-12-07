library(tidyverse)
library(edgeR)
library(anota2seq)
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

# 0 gene exclusion & TMM normalize 
countsTumor %>% column_to_rownames("EnsemblID") %>%
  .[!(apply(., 1, function(x) any(x == 0))),] %>%    # from anota2seqRemoveZeroSamples()
  DGEList() %>%
  calcNormFactors(method = "TMM") %>%
  cpm() -> countsTumor.tmm

countsHost %>% column_to_rownames("EnsemblID") %>%
  .[!(apply(., 1, function(x) any(x == 0))),] %>%   # from anota2seqRemoveZeroSamples()
  DGEList() %>% 
  calcNormFactors(method = "TMM") %>%
  cpm() -> countsHost.tmm


## Split Tumour
PDX_count_spliter <- function(dat, conversion, fraction){
  dat %>% as_tibble(rownames = "EnsemblID") %>% dplyr::select(c("EnsemblID", conversion[[fraction]])) %>%
    rename_at(conversion[[fraction]], ~ conversion[["standard"]])
}

## Split Tumour
PDX_count_spliter(countsTumor.tmm, sample_names, "total") -> countsTumor.tmm_tot
PDX_count_spliter(countsTumor, sample_names, "total") -> countsTumor_tot
countsTumor.tmm_tot %>% write_tsv("Processed_data/countsTumor_tmm_tot.tsv")
countsTumor_tot %>% write_tsv("Processed_data/countsTumor_raw_tot.tsv")

PDX_count_spliter(countsTumor.tmm, sample_names, "polysome") -> countsTumor.tmm_pol
PDX_count_spliter(countsTumor, sample_names, "polysome") -> countsTumor_pol
countsTumor.tmm_pol %>% write_tsv("Processed_data/countsTumor_tmm_pol.tsv")
countsTumor_pol %>% write_tsv("Processed_data/countsTumor_raw_pol.tsv")

## Split Stroma
PDX_count_spliter(countsHost.tmm, sample_names, "total") -> countsHost.tmm_tot
PDX_count_spliter(countsHost, sample_names, "total") -> countsHost_tot
countsHost.tmm_tot %>% write_tsv("Processed_data/countsHost_tmm_tot.tsv")
countsHost_tot %>% write_tsv("Processed_data/countsHost_raw_tot.tsv")

PDX_count_spliter(countsHost.tmm, sample_names, "polysome") -> countsHost.tmm_pol
PDX_count_spliter(countsHost, sample_names, "polysome") -> countsHost_pol
countsHost.tmm_pol %>% write_tsv("Processed_data/countsHost_tmm_pol.tsv")
countsHost_pol %>% write_tsv("Processed_data/countsHost_raw_pol.tsv")

# Boxplots
raw = list(countsTumor_tot, countsTumor_pol,
           countsHost_tot, countsHost_pol)

tmm = list(countsTumor.tmm_tot, countsTumor.tmm_pol,
                 countsHost.tmm_tot, countsHost.tmm_pol)

plot_title = list("Tumor_tot", "Tumor_pol",
                  "Host_tot", "Host_pol")
for (i in 1:4){
  
raw[[i]] %>% column_to_rownames("EnsemblID") -> raw_i
tmm[[i]] %>%  column_to_rownames("EnsemblID") -> tmm_i
boxplot(raw_i + 1, log = "y",
          main = paste(plot_title[[i]], "raw", sep = " "), las = 2)
  
boxplot(tmm_i + 1, log = "y",
          main = paste(plot_title[[i]], "tmm", sep = " "),  las = 2)
  }



