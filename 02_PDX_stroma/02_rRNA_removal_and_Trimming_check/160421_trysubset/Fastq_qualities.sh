#!/bin/bash

# Title: Basic Run of SortMeRNA
# Author details: Author: Jacobo Solorzano, Contact details: jacobo.solorzano@inserm.fr
# Script and data info: This script takes the output of deplete_rRNAs.sh and agrupates all
###################### in the output folder of the trimming, runing fastqc and mulriqc afterwards

# Imput data: 
# INPUT: FASTQs from deplete_rRNAs.sh
# OUTPUT: FASTQC individual qualities and multiqc report
# Creation 16/04/21
# Modified -

CURRDIR=$(pwd)
INPUTFASTQDIR="${CURRDIR}/01_Input"
FASTQDIR="${CURRDIR}/02_Output/BBDuk"

#get all the FASTQs together
mv "${CURRDIR}/02_Output/sortmerna/out/rRNAfree/*" "${FASTQDIR}"


# Run fastqc
fastqc "${INPUTFASTQDIR}/" -o "${FASTQDIR}/"
fastqc "${FASTQDIR}/*.fastq"

# Run multiQC
multiqc ${FASTQDIR}/