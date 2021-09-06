#!/bin/bash
#SBATCH -J kf1_rRNAfree_quality  #jobname
#SBATCH -o stdout/kf1_rRNAfreeq.out #output file name
#SBATCH -e stdout/kf1_rRNAfreeq.out #error file name
#SBATCH --mem=64G #memory reservation
#SBATCH --cpus-per-task=14 #ncpu on the same node
#SBATCH --mail-type=BEGIN,END,FAIL

#Purge any previous modules
module purge

#Load the application
module load module load bioinfo/sortmerna-4.0.0
module load bioinfo/bbmap_38.31
module load bioinfo/FastQC_v0.11.7
module load bioinfo/MultiQC-v1.7

#execute

./fastq_processing.sh _R1.fastq.gz _R2.fastq.gz kf_fastqs kf_fastqs
