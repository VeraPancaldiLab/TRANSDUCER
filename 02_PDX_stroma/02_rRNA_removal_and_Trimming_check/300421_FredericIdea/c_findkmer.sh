#!/bin/bash
#SBATCH -J findkmer  #jobname
#SBATCH -o stdout/findkmer.out #output file name
#SBATCH -e stdout/findkmer.err #error file name
#SBATCH --mem=120G #memory reservation
#SBATCH --cpus-per-task=15 #ncpu on the same node
#SBATCH --mail-type=BEGIN,END,FAIL

#Purge any previous modules
module purge

#Load the application
module load bioinfo/bbmap_38.31

#execute
kmercountexact.sh in=PDAC221Tkf1_fwd_P60.fq out=stdout fastadump=f | sort -k2,2rn - > findkmer/P60.out
kmercountexact.sh in=PDAC221Tkf1_fwd_P50.fq out=stdout fastadump=f | sort -k2,2rn - > findkmer/P50.out


