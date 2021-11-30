#!/bin/bash
# Title: Basic Run of SortMeRNA + two separate types of Trimming
# Author details: Author: Jacobo Solorzano, Contact details: jacobo.solorzano@inserm.fr
# Script and data info: This script does the simplest run for indexing a rRNA db
###################### and depleting the rRNA reads
# Imput data: 
# INPUT: rRNA db and FASTQs (NOT .GZ) where rRNA will be depleated, and afterwards trimmed
# OUTPUT: rRNA free FASTQs (r_R1.fastq ...), rRNA free + force trimmed FASTQs (rft_R1.fastq)
#         and rRNA free + Kmer filtered FASTQs (artifacts phfix and adapters).
# Creation 09/04/21
# Modified 20/04/21 by JS - added imput options

CURRDIR=$(pwd)
fwd_id=$1 # ie _R1.fastq
rev_id=$2 # ie _R2.fastq
INDIR=${CURRDIR}/$3
OUTDIR=${CURRDIR}/$4

echo $INDIR
echo ${OUTDIR}

# pending
# fix for choosing each type!
# some file with parameters for each type
# names for 

for file in $(ls ${INDIR} | grep "${fwd_id}")
do
	name=$(echo "$file" | awk -F_R1 '{print $1}')
	fastq1=${INDIR}/"${name}${fwd_id}"
	fastq2=${INDIR}/"${name}${rev_id}"
	
	echo $fastq1
	echo $fastq2
#rRNA depletion
sortmerna --ref ${INDIR}/PDXrRNAs_ensembl75.fasta \
--reads ${fastq1} \
--reads ${fastq2} \
--workdir ${OUTDIR}/sortmerna/ \
--fastx rRNA_aligned \
--other ${OUTDIR}/sortmerna/out/rRNAfree/${name} \
--paired_in \
--out2 

# 
# rm -rf ${OUTDIR}/sortmerna/kvdb/* #this dir has to be removed each round
# 
# # format names to add r and identify this as rRNA free fastq
# fastq1_r=${OUTDIR}/sortmerna/out/rRNAfree/${name}r_R1.fastq
# fastq2_r=${OUTDIR}/sortmerna/out/rRNAfree/${name}r_R2.fastq
# 
# mv ${OUTDIR}/sortmerna/out/rRNAfree/${name}_fwd.fq ${fastq1_r}
# mv ${OUTDIR}/sortmerna/out/rRNAfree/${name}_rev.fq ${fastq2_r}
# 
# # Trimming by BBDuk
# ## force trimming 12bp at the begining and 13bp in the end, so I end upo with 50BP reads
# ## these files will have rft before _R*.fastq
# bbduk.sh in1=${fastq1_r} in2=${fastq2_r} \
# out1=${OUTDIR}/BBDuk/${name}rft_R1.fastq out2=${OUTDIR}/BBDuk/${name}rft_R2.fastq \
# ftl=12 ftr=62


## Kmer filter will remove the Ilumina spike ins/artifacts and adapters 
## of at least 31kmer match. These files will have rkf before _R*.fastq
## CAREFULL CHOSE THE FASTQS ACDORDING THE PRIGIN YOU WANT
#bbduk.sh -Xmx54g in=${fastq1} in2=${fastq2} out=${OUTDIR}/${name}kf1${fwd_id} out2=${OUTDIR}/${name}kf1${rev_id} ref=adapters,artifacts,phix kmask=N k=25 mink=8 hdist=1 hdist2=0 stats=${OUTDIR}/${name}kf1_stats.txt tbo
# remove fastqs for space
#rm ${fastq1}
#rm ${fastq2}
 
done


fastqc $(ls ${OUTDIR}/sortmerna/out/rRNAfree/*.fastq.gz)
multiqc ${OUTDIR}/sortmerna/out/rRNAfree/*

exit 0
