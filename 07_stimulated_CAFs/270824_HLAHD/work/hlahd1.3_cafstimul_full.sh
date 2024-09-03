#!/bin/bash

#SBATCH --mem=300G
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=jacobo.solorzano@inserm.fr
#SBATCH --mail-type=FAIL

module=/save/jsolorzano/scripts/hlahd/program_module

module load -f ${module}

# Assign arguments
SAMPLE_NAME=$1
OUTPUT_PATH="/work/jsolorzano/results/hlahd_cafstimul"

# Construct the FASTQ file names
FASTQ_R1="/save/jsolorzano/data/cafstimul/${SAMPLE_NAME}_R1_001.fastq"
FASTQ_R2="/save/jsolorzano/data/cafstimul/${SAMPLE_NAME}_R2_001.fastq"

# Paths to the GENE_SPLIT_FLIT and DICTIONARY
GENE_SPLIT_FLIT=/save/jsolorzano/software/hlahd.1.3.0.2/HLA_gene.split.3.32.0.txt
DICTIONARY=/save/jsolorzano/software/hlahd.1.3.0.2/dictionary

hlahd.sh -t 20 -m 30 ${FASTQ_R1} ${FASTQ_R2} ${GENE_SPLIT_FLIT} ${DICTIONARY} ${SAMPLE_NAME} ${OUTPUT_PATH}
