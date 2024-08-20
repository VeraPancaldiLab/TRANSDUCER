#Filter for human cohort dataset

################################################################################
#'Create a function to select sample regarding :
#'proportion of tumor cells
#'proportion of stroma 
#'proportion of both 
#'threshold for acinar tissue
#'and islets 
#'with the choice of estimation between deconvolution or histology
################################################################################



#Filter function
#'@description
#' this function is made to filter the human cohort samples regarding to clinical data and molecular phenotype data
#' it is possible to chose thresholds for proportion of tumor cells, acinar cells, islets, and tumor tissue (tumor + stroma)
#' the estimation method is either deconvolution or histology
#'@param df data frame of human cohort RNA-seq data from http://linkedomics.org/data_download/CPTAC-PDAC/mRNA_RSEM_UQ_log2_Tumor.cct
#'@param clinical_data clinical data sheet of the human cohort metadata
#'@param molecular_phenotype_df molecular phenotype sheet of the human cohort metadata
#'@param estimation_method deconvolution or histology method to estimate proportion of each cell type
#'@param neoplastic_min minimum proportion of tumor cells wanted
#'@param acinar_max maximum of acinar cells 
#'@param islet_max maximum of islet cells
#'@param tumor_tissue_min minimum of tumoral tissue (tumor cells + stromal cells)

dataFilter <- function(df, clinical_data, molecular_phenotype_df, estimation_method, neoplastic_min, acinar_max, islet_max, tumor_tissue_min){
  
  clinical_df <- dplyr::select(clinical_data, case_id, histology_diagnosis, vital_status, is_this_patient_lost_to_follow_up)
  
  non_PDAC_sample <- dplyr::filter(clinical_df, histology_diagnosis != "PDAC")
  
  if (estimation_method == "deconvolution"){
    
    estimation_df <- dplyr::select(molecular_phenotype_df, case_id, contains("deconv")) %>%
      mutate(tumor_tissue = epithelial_cancer_deconv + stromal_deconv)
    
    filtered_sample <- dplyr::filter(estimation_df, epithelial_cancer_deconv >= neoplastic_min & 
                                       mature_exocrine_endocrine_deconv <= acinar_max + islet_max &
                                       tumor_tissue >= tumor_tissue_min) 
    
  } else if (estimation_method == "histology"){
    estimation_df <- dplyr::select(molecular_phenotype_df, case_id, contains("histology")) %>%
      mutate(tumor_tissue = neoplastic_cellularity_histology_estimate + stromal_histology_estimate)
    
    filtered_sample <- dplyr::filter(estimation_df, neoplastic_cellularity_histology_estimate >= neoplastic_min & 
                                       acinar_histology_estimate <= acinar_max &
                                       islet_histology_estimate <= islet_max &
                                       tumor_tissue >= tumor_tissue_min)
  }
  
  filtered_df <- dplyr::select(df, Gene, filtered_sample$case_id, -non_PDAC_sample$case_id)
  
  return(filtered_df)
}

#Exemple of use
# library(readr)
# 
# #Import dataset
# mRNA_tumor <- read_delim("PDAC_LinkedOmics_Data/mRNA_RSEM_UQ_log2_Tumor.cct",
#                          delim = "\t", escape_double = FALSE,
#                          trim_ws = TRUE) %>%
#   rename(Gene = ...1)
# 
# mmc1 <- read_excel("mmc1.xlsx", 
#                    sheet = "Clinical_data")
# mmc2 <- read_excel("mmc1.xlsx",
#                    sheet = "Molecular_phenotype_data") %>%
#   mutate_at(vars(immune_deconv:`necrosis_(%OF_TUMOR_WITH_NECROSIS)_histology_estimate`, KRAS_VAF), as.numeric) %>%
#   mutate_at(vars(Bailey:Moffitt), as.factor)
# 
# ################################################################################
# estimation_method = "histology" #histology | deconvolution
# neoplastic_min = 0.2
# acinar_max = 0.05
# islet_max = 0.1
# tumor_tissue_min = 0.8
# ################################################################################

# mRNA_tumor_filtered <- dataFilter(mRNA_tumor, mmc1, mmc2, estimation_method, neoplastic_min, acinar_max, islet_max, tumor_tissue_min)
# 
## Same but manually
# estimation_deconv <- dplyr::select(mmc2, case_id, contains("deconv"))%>%
#   mutate(tumor_tissue = epithelial_cancer_deconv + stromal_deconv)
# estimation_histo <- dplyr::select(mmc2, case_id, contains("histology"))%>%
#   mutate(tumor_tissue = neoplastic_cellularity_histology_estimate + stromal_histology_estimate)
# 
# 
# filtered_sample <- dplyr::filter(estimation_deconv, epithelial_cancer_deconv >= neoplastic_min&
#                                    mature_exocrine_endocrine_deconv <= acinar_max + islet_max &
#                                    tumor_tissue >= tumor_tissue_min)
# 
# filtered_sample2 <- dplyr::filter(estimation_histo, neoplastic_cellularity_histology_estimate >= neoplastic_min &
#                                    acinar_histology_estimate <= acinar_max &
#                                    islet_histology_estimate <= islet_max &
#                                    tumor_tissue >= tumor_tissue_min)
# 
# clinical_df <- dplyr::select(mmc1, case_id, histology_diagnosis, vital_status, is_this_patient_lost_to_follow_up)
# 
# non_PDAC_sample <- dplyr::filter(clinical_df, histology_diagnosis != "PDAC")
# 
# filtered_df <- dplyr::select(mRNA_tumor, Gene, filtered_sample2$case_id, -non_PDAC_sample$case_id)
