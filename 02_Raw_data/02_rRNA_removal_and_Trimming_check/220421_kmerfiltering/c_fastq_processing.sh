#!/bin/bash
#SBATCH -J kf_quality  #jobname
#SBATCH -o stdout/kfq1.out #output file name
#SBATCH -e stdout/kfq1.out #error file name
#SBATCH --mem=64G #memory reservation
#SBATCH --cpus-per-task=14 #ncpu on the same node
#SBATCH --mail-type=BEGIN,END,FAIL

#Purge any previous modules
module purge

#Load the application
module load bioinfo/bbmap_38.31
module load bioinfo/FastQC_v0.11.7
module load bioinfo/MultiQC-v1.7

#execute

./fastq_processing.sh _R1.fastq.gz _R2.fastq.gz raw_fastqs kf_fastqs

