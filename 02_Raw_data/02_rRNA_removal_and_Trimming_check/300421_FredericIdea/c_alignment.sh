#!/bin/bash
#SBATCH -J checGC_SMAP  #jobname
#SBATCH -o stdout/checGC_SMAP.out #output file name
#SBATCH -e stdout/checGC_SMAP.err #error file name
#SBATCH --mem=120G #memory reservation
#SBATCH --cpus-per-task=15 #ncpu on the same node
#SBATCH --mail-type=BEGIN,END,FAIL

#Purge any previous modules
module purge

#Load the application
module load bioinfo/STAR-2.6.0c

#execute
./alignment.sh
