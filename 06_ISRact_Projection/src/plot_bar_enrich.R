library(tidyverse)
## Alternative plotting of GO resulting of enrichPathway()
## as default Cluisterprofiler fuction do not work
convert_to_numeric <- function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
}
plot_bar_enrich <- function(x, p.adjust_th, Title) {
  mutate(x, GeneRatioNum = sapply(GeneRatio, convert_to_numeric)) %>% 
    arrange(GeneRatioNum) %>% 
    mutate(Description = as_factor(Description)) %>% 
    dplyr::filter(p.adjust < p.adjust_th) %>%
    slice_tail(n=20) %>%
    ggplot(aes(x = GeneRatioNum, y = Description, fill = p.adjust)) +
    geom_bar(stat="identity") +
    scale_fill_gradient(low = "red", high = "blue", guide=guide_colourbar(reverse = TRUE))+
    theme_bw() +
    labs(title = Title) +
    xlab("Gene Ratio") +
    ylab("")
}

# Example
# TEs_5perc_up_enrich_PA@result %>%  
#plot_bar_enrich(0.05,"Reactome Upregulated in IC.6 (extreme gene weights)")
