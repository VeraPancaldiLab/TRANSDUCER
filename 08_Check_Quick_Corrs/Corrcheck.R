setwd("/home/jacobo/Documents/02_TRANSDUCER/08_Check_Quick_Corrs/")
library(tidyverse)
library(Hmisc)

# In this script you can compare any of the components obtained from Tumour translatome/with transcription of tumour/stroma
# REMEMBER RUNNING the following command in the directory to make text editable in InkScape:
## for i in ./*.svg; do ./fixsvgs.sh "$i"; done

# Functions
plot_gene_correlation <- function(data, x_col, y_col, plot_title = fraction_expression, output_path = "02_Output/") {
  
  # Calculate Spearman correlation
  gene_cor <- rcorr(data[[x_col]], data[[y_col]], type = "spearman")
  
  # Format stats for the subtitle
  stats <- paste("Spearman: R = ", round(gene_cor$r["x","y"], 2),
                 ", pval = ", round(gene_cor$P["x","y"], 4), sep = "")
  
  # Create plot
  p <- ggplot(data, aes_string(x = x_col, y = y_col)) +
    geom_point(colour = "red") +
    geom_smooth(method = "lm") +
    theme_bw() +
    labs(title = plot_title, subtitle = stats, x = x_col, y = y_col)
  
  # Save plot
  output_filename <- paste0(output_path, "corrplot_", x_col, "_vs_", y_col, "_scatter.svg")
  ggsave(output_filename, plot = p, width = 3, height = 3)
  
  return(p)  # Return the plot object if you want to view it directly
}

# PARAMETERS
###############################################################################
fraction_expression = "Stroma" # "Tumour" | "Stroma" 
gene = "Cxcl13"
component = "IC.4.cyt" # IC.6.TEs | IC.4.cyt | ISRact | Tumour_IC.1
###############################################################################

# load the component information
ISRact_info <- read_tsv("../06_ISRact_Projection/data/Sauyeun_PDX/sample_info.tsv") %>% 
  dplyr::select(sample, PAMG, matches("ICA")) %>% 
  dplyr::rename(ISRact = "ICA3") %>% 
  dplyr::rename_with(~ str_replace(., "^ICA(\\d+)$", "Tumour_IC.\\1"), starts_with("ICA"))

IC_info <- read_tsv("../02_PDX_stroma/03_Analysis/100122_ICABoot/02_Output/all_ICA_samplescores.tsv") %>%
  dplyr::select(-ICA1, -ISRact, -PAMG) # remove as the are only in 20 samples

component_information <- inner_join(ISRact_info, IC_info)

# Load selected gene expression
if (fraction_expression == "Tumour"){
  load("../02_PDX_stroma/00_Data/Remy_processed_data/all_RNA/Tumeur/Hcpmallrna.RData")
  gene_expression <- Hcpmallrna
  
  ensbl_gene <- read_tsv("../02_PDX_stroma/00_Data/GRCh38.p13_ensemblvsgenename.txt")
  translate <- deframe(ensbl_gene)
  
} else if (fraction_expression == "Stroma"){
  load("../02_PDX_stroma/00_Data/Remy_processed_data/all_RNA/Stroma/Mcpmallrna.RData")
  gene_expression <- Mcpmallrna # This must be quantile normalized data, probably log2
  
  ensbl_gene <- read_tsv("../02_PDX_stroma/00_Data/GRCm39_ensemblvsgenename.txt")
  translate <- deframe(ensbl_gene)
  
}

rownames(gene_expression) <- gene_expression %>% rownames(.) %>% translate[.]
gene_expression <- gene_expression[!is.na(rownames(gene_expression)),]

# Join and prepare for plotting
gene_plots <- as_tibble(gene_expression, rownames = "geneID") %>%
  pivot_longer(-geneID, names_to = "sample", values_to = "expression") %>%
  pivot_wider( values_from = "expression", names_from = "geneID", values_fn = mean) %>% 
  left_join(component_information, by = "sample")



# Plot
plot_gene_correlation(gene_plots, component, gene) 

