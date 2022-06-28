# Package loading
library(tidyverse)
library(ggplot2)
setwd("/home/jacobo/Documents/02_TRANSDUCER/03_IC3Characterization/01_PDX/020821_IC3_Into_Two_Groups/")

group_type <- "5vs5" #4vs4noPDAC024T | 5vs5

# MAIN
ic3 <- read_tsv("../Remy_processed_data/samplesIC3_custom.csv") %>% mutate(as.double(ICA3SampleWeight)) %>%
  as.data.frame(.) %>%
  column_to_rownames(., var = "CITID") %>%
  .[,c("ICA3SampleWeight", "Diabetes"), drop = FALSE] %>% arrange(ICA3SampleWeight)


col_binary <- rep(NA,nrow(ic3))

if (group_type == "4vs4noPDAC024T"){
  col_binary[1:4] <- "low_IC3"
  col_binary[(nrow(ic3)-4):(nrow(ic3)-1)] <- "high_IC3"
  
} else if (group_type == "5vs5") {
  col_binary[1:5] <- "low_IC3"
  col_binary[(nrow(ic3)-4):nrow(ic3)] <- "high_IC3"
}

ic3$groups <- col_binary


ggplot(ic3, aes(x = ICA3SampleWeight)) +
  geom_density() + geom_rug(aes(color = groups)) + scale_colour_manual(values = c("green", "red")) + theme_classic()



write_tsv(ic3 %>% rownames_to_column(), paste(c("02_Output/samplesIC3_",group_type,".tsv"), collapse = ""))

          