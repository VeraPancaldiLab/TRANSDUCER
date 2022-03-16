#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
source("functions.R")
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/")
ICA_cyt <- read_rds("02_Output/ICA_cyt.RDS")
ICA_pol <- read_rds("02_Output/ICA_pol.RDS")
ICA_TEs <- read_rds("02_Output/ICA_TEs.RDS")

# Comparison of sample weights component wide
A_cyt <- as_tibble(ICA_cyt$A, rownames = "sample") %>% rename_with( ~paste0(.,".cyt"), -sample)
A_pol <- as_tibble(ICA_pol$A, rownames = "sample") %>% rename_with( ~paste0(.,".pol"), -sample)
A_TEs <- as_tibble(ICA_TEs$A, rownames = "sample") %>% rename_with( ~paste0(.,".TEs"), -sample)

Acomparison <- inner_join(A_cyt, A_pol, by = "sample") %>%
  inner_join(A_TEs, by = "sample") %>%  column_to_rownames("sample")

## cyt ~ pol
dplyr::select(Acomparison, ends_with(c("cyt", "pol"))) %>%
  formatted_cors(cor.stat = "pearson") %>%
  dplyr::filter(str_detect(measure1, ".cyt"), str_detect(measure2, ".pol")) %>%
  ggplot(aes(measure2, measure1, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title= "cyt vs pol", 
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  rotate_x_text(angle = 45)

## cyt ~ TEs
dplyr::select(Acomparison, ends_with(c("cyt", "TEs"))) %>%
    formatted_cors(cor.stat = "pearson") %>%
  dplyr::filter(str_detect(measure1, ".cyt"), str_detect(measure2, ".TEs")) %>%
  ggplot(aes(measure2, measure1, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title= "cyt vs TEs", 
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  rotate_x_text(angle = 45)

## pol ~ TEs
dplyr::select(Acomparison, ends_with(c("pol", "TEs"))) %>%
  formatted_cors(cor.stat = "pearson") %>%
  dplyr::filter(str_detect(measure1, ".pol"), str_detect(measure2, ".TEs")) %>%
  ggplot(aes(measure2, measure1, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title= "pol vs TEs", 
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  rotate_x_text(angle = 45)
