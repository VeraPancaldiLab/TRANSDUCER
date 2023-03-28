#!/bin/bash
#01/04/21
#jacobo.solorzano@inserm.fr
FAFILE="/mnt/SERVER-CRCT-STORAGE/CRCT06/5To/Jacobo/00_Other/reference_genomes/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"
GTFFILE="/mnt/SERVER-CRCT-STORAGE/CRCT06/5To/Jacobo/00_Other/reference_genomes/Homo_sapiens.GRCh38.105.chr.gtf"
FASTQDIR="/home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/00_Data/PJ2212479/RAWDATA/"
index=0
seqlength=100 # set to the length of base pair sequencing
threads=20
WDIR="/home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/00_Data/240323_Alignment"

# Index 
if [ ! -d $WDIR/StarIndexed ]
then
  mkdir $WDIR/StarIndexed
  
  STAR --runMode genomeGenerate \
  --genomeDir $WDIR/StarIndexed \
  --genomeFastaFiles $FAFILE \
  --sjdbGTFfile $GTFFILE \
  --sjdbOverhang $seqlength 
else
  echo "Loading current indexed genome at" $WDIR
fi

# Load genome into memory
STAR --genomeDir $WDIR/StarIndexed \
--genomeLoad LoadAndExit

# Mapping
for file in $(ls $FASTQDIR | grep "R1")
do
	name=$(echo "$file" | awk -F_R1 '{print $1}')
	fastq1="${name}"_R1.fastq.gz
	fastq2="${name}"_R2.fastq.gz

STAR --limitBAMsortRAM 10000000000 \
--genomeLoad LoadAndKeep \
--genomeDir $WDIR/StarIndexed \
--runThreadN $threads \
--readFilesIn <(gunzip -c $FASTQDIR/$fastq1) <(gunzip -c $FASTQDIR/$fastq2)  \
--outSAMtype BAM SortedByCoordinate  \
--outFileNamePrefix $WDIR/output_bams/${name}_  

#STAR --genomeLoad LoadAndKeep \ #02/04/21
#--genomeDir $WDIR/StarIndexed\
#--runThreadN 14 \
#--readFilesIn <(gunzip -c $WDIR/$fastq1) <(gunzip -c $WDIR/$fastq2) \
#--sjdbGTFfile $WDIR/combined.gtf \
#--outFilterType BySJout \ #Default will be "normal". "BySJout" is for 2-pass m$
#--outFilterMultimapNmax 100 \
#--outSAMtype BAM SortedByCoordinate  \
#--outFileNamePrefix output_bams/${name}_

#!!! CAREFULL THIS REMOVE INPUT FQ (IN case theres a lack of space)
# rm $INDIR/$fastq1
# rm $INDIR/$fastq2
done

# Unload genome of memory

wait 

STAR --genomeLoad Remove --genomeDir $WDIR/StarIndexed

exit
