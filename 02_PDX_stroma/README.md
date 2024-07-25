<!-- PROJECT LOGO -->

<h3 align="center">

TRANSDUCER - ISRact stroma Analysis

</h3>

<h2 align="left">

Summary

</h2>

This Section of the repository covers the analysis performed on the PDX stroma to explore the context supporting ISRact tumour cells. In order to clarify this stroma, and taking teh advantage of the separation of tumour and stroma fractions of reads coming from Patient derived xenograft (PDX) models.

For this, first the stroma was analysed with supervised dimensionality reduction tools, where based in prior knowledge, gene transcription is codified into cell type proportions or transcription factor activities. For cell type proportions mMCPcounter and ImmuCC signatures were utilized. For estimating transcription factor activities we made use of the Dorothea PAN Cancer and GTEx regulons through Viper. Gtex used for stromal TF activity calculation, and PanCancer for the tumour TF activity, also included to asses whether some stromal phenotype could be linked to a tumour transcriptional feature.

Then, making use of Independent component analysis (ICA), components related to ISRact were scanned and processed in transcription and translation, exploring the relation with technical biasses, celltype proportions, tumoral and stromal TF activities as well as otehr tailored analysis depending on wether transcription or translation data was used. #! COMPLETE??


In the next sections, starting with the transcriptional analysis ofteh stroma and folloeing with the translational one the figures of the thesis and the scripts that originat them are depicted. In occassions further work has been done, and is indicated.

<h2 align="left">

Transcription

</h2>

<h3 align="left">

FigR1: IL18-IFNG related signals are associated with ISRact stroma, with potential implications in
immune infiltration (pg 80)

</h3>

<img src="replacewithpath/FigR1/FigR1.png">


<h3 align="left">

FigR2: IC.4_cyt in PDX stroma is the transcriptional component most relevant to ISRAct tumour
phenotype (pg 87)

</h3>

<img src="replacewithpath/FigR2/FigR2.png">


<h3 align="left">

FigR3: IC.4_cyt captures the relations with immunity observed in the preliminary analysis based on
correlations with Immune cell type proportions and TF activities (pg 89)

</h3>

<img src="replacewithpath/FigR3/FigR3.png">


FigR4: Analysis of IC.4cyt gene contributions confirms the relation of ISRact with immunity (pg 93)

</h3>

<img src="replacewithpath/FigR4/FigR4.png">




<h2 align="left">

Translation

</h2>





FigR5: Translation efficacy calculation and the exclusion of PDAC001..95
FigR6: IC.6_TEs as the stromal translational component representing ISRact lacks relation with TF
activities or Deconvolved cell type proportions....99
FigR7: mRNA levels of RBPs are specifically associated with IC.6_TEs suggesting an involvement of
m6A methylation, mRNA degradation and splicing in ISRact stroma..101
FigR8: Gene contribution analysis reveals a diverse array of signals being captured by IC.6_TEs,
among which is App 103
FigR9: The genes selected as Up translated and Downtranslated from IC.6_TEs show enrichment for
a wide range of functions.. 105
FigR10: Anota2seqUtils analysis of Translation related signatures suggest some cells of the stroma
could also be undergoing an ISR activation...107
FigR11: Isolated Anota2seqUtils analysis of the different gene regions shows a relation with length
and a spurious association with DRACH motif presence..108
8

FigR12: Isolated Anota2seqUtils analysis of the different gene regions ATtRACT motifs still shows
lack of specificity and points to an involvement of mRNA degradation/optimal codon usage110
FigR13: Integrated analysis of IC.6_TEs highly contributing genes shows sparsity, but gives clues to
an array of different functions, highlighting a potential role of matrix disorganisation, and
immunosuppressive signals through downregulation of Lum and Ccl11 respectively....112
FigR14: Correlation analysis between tumour translation components, and stromal Transcription
and translation ones highlight the need for a bigger cohort....114
FigR15: Literature defined CAF subtypes do not seem associated with ISRact nor IC.4_cyt..116
FigR16: CAFs show different enrichment scores depending on whether total, efficiently translated
mRNA levels, or Translation Efficacies are utilised..117
FigR17: Stimulated CAF data is under a strong batch effect and individual stimuli analysis suggests a
problem in the experiment rendering the data inadequate for analysis120
FigR18: Survival analysis reveals inconsistencies for detecting a relation between ISRact and survival
based on PHGDH/CBS average expression....126
FigR19: PCA allows the separation of ISRact high and low PDX samples based on transcriptomics.
... 130
FigR20: ISRactPCA captures signals consistently in Shin et al. PDX and CPTAC cohort,
recapitulating signals related to gemcitabine resistance....134
FigR21: ISRactPCA high vs low samples proteomic comparison does not show differential expressed
proteins, while supporting glucuronidation involvement..137
FigR22: ISRactPCA shows a negative correlation with Lymphocyte estimated proportions in the
CPTAC human cohort.. 138
FigR23: CCLE cell line data projection onto ISRactPCA is not consistent with the in-vitro defined
markers 141
FigR24: ISRactPCA could capture characteristics related to Basal/Classical classification....145
FigR25: Analysis of ISRactPCA on snRNAseq reveals two clusters corresponding to the
ISRactPCA_high and ISRactPCA_low phenotypes..146
FigR26: Cluster 5 and 10 marker genes have a substantial overlap with ISRactPCA genes and show
enrichment in relevant terms such as innate immunity or ISR activation.148
FigR27: ISRactPCA score is not predictive of survival across the major PDAC cohorts150
SupFigR1: Detailed Over Representation Analysis (ORA) of IC.6_TEs up and down-translated
contributing genes as summarised in FigR9 B....154
