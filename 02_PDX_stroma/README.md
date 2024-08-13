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

################################################################################

<h2 align="left">

Transcription

</h2>

<h3 align="left">

FigR1: IL18-IFNG related signals are associated with ISRact stroma, with potential implications in
immune infiltration (pg 80)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR1/FigR1.png" style="background-color:white;" />


<h3 align="left">

FigR2: IC.4_cyt in PDX stroma is the transcriptional component most relevant to ISRAct tumour
phenotype (pg 87)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR2/FigR2.png" style="background-color:white;" />


<h3 align="left">

FigR3: IC.4_cyt captures the relations with immunity observed in the preliminary analysis based on
correlations with Immune cell type proportions and TF activities (pg 89)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR3/FigR3.png" style="background-color:white;" />


<h3 align="left">

FigR4: Analysis of IC.4cyt gene contributions confirms the relation of ISRact with immunity (pg 93)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR4/FigR4.png" style="background-color:white;" />


################################################################################
<h2 align="left">

Translation

</h2>

<h3 align="left">

FigR5: Translation efficacy calculation and the exclusion of PDAC001 (pg 95)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR5/FigR5.png" style="background-color:white;" />

<h3 align="left">

FigR6: IC.6_TEs as the stromal translational component representing ISRact lacks relation with TF activities or Deconvolved cell type proportions (pg 99)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR6/FigR6.png" style="background-color:white;" />

<h3 align="left">

FigR7: mRNA levels of RBPs are specifically associated with IC.6_TEs suggesting an involvement of m6A methylation, mRNA degradation and splicing in ISRact stroma (pg 101)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR7/FigR7.png" style="background-color:white;" />

<h3 align="left">

FigR8: Gene contribution analysis reveals a diverse array of signals being captured by IC.6_TEs, among which is App (pg 103)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR8/FigR8.png" style="background-color:white;" />

<h3 align="left">

FigR9: The genes selected as Up translated and Downtranslated from IC.6_TEs show enrichment for a wide range of functions (pg 105)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR9/FigR9.png" style="background-color:white;" />

<h3 align="left">

FigR10: Anota2seqUtils analysis of Translation related signatures suggest some cells of the stroma could also be undergoing an ISR activation (pg 107)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR10/FigR10.png" style="background-color:white;" />

<h3 align="left">

FigR11: Isolated Anota2seqUtils analysis of the different gene regions shows a relation with length and a spurious association with DRACH motif presence (pg 108)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR11/FigR11.png" style="background-color:white;" />

<h3 align="left">

FigR12: Isolated Anota2seqUtils analysis of the different gene regions ATtRACT motifs still shows lack of specificity and points to an involvement of mRNA degradation/optimal codon usage (pg 110)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR12/FigR12.png" style="background-color:white;" />

<h3 align="left">

FigR13: Integrated analysis of IC.6_TEs highly contributing genes shows sparsity, but gives clues to an array of different functions, highlighting a potential role of matrix disorganisation, and immunosuppressive signals through downregulation of Lum and Ccl11 respectively (pg 112)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR13/FigR13.png" style="background-color:white;" />

################################################################################
<h2 align="left">

Comparsion Transcription VS Translation 

</h2>


<h3 align="left">

FigR14: Correlation analysis between tumour translation components, and stromal Transcription and translation ones highlight the need for a bigger cohort (pg 114)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR14/FigR14.png" style="background-color:white;" />


################################################################################
<h2 align="left">

CAF involvement

</h2>

<h3 align="left">

FigR15: Literature defined CAF subtypes do not seem associated with ISRact nor IC.4_cyt (pg 116)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR15/FigR15.png" style="background-color:white;" />

<h3 align="left">

FigR16: CAFs show different enrichment scores depending on whether total, efficiently translated mRNA levels, or Translation Efficacies are utilised (pg 117)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR16/FigR16.png" style="background-color:white;" />

<h3 align="left">

FigR17: Stimulated CAF data is under a strong batch effect and individual stimuli analysis suggests a problem in the experiment rendering the data inadequate for analysis (pg 120)

</h3>

<img src="https://github.com/VeraPancaldiLab/TRANSDUCER/tree/main/FIGURES_THESIS/FigR17/FigR17.png" style="background-color:white;" />

TEST

<img src="https://github.com/j-solor/Drosophila-miRNA-PcG-circuits/blob/main/misc/Figure1.png" style="background-color:white;" />

