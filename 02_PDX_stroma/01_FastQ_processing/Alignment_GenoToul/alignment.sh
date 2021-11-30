#!/bin/bash
#01/04/21
#jacobo.solorzano@inserm.fr

INDIR="/home/jsolorzano/work/PDX_data"

# Index 
#mkdir combined_reference/StarIndexed_combined

#STAR --runMode genomeGenerate \
#--genomeDir combined_reference/StarIndexed_combined \
#--genomeFastaFiles combined_reference/combined.fasta \
#--sjdbGTFfile combined_reference/combined.gtf \
#--sjdbOverhang 75   # set to the length of base pair sequencing 

# Load genome into memory
STAR --genomeDir combined_reference/StarIndexed_combined \
--genomeLoad LoadAndExit

# Mapping
for file in $(ls $INDIR | grep "R1")
do
	name=$(echo "$file" | awk -F_R1 '{print $1}')
	fastq1="${name}"_R1.fastq.gz
	fastq2="${name}"_R2.fastq.gz

STAR --limitBAMsortRAM 10000000000 \
--genomeLoad LoadAndKeep \
--genomeDir combined_reference/StarIndexed_combined \
--runThreadN 14 \
--readFilesIn <(gunzip -c $INDIR/$fastq1) <(gunzip -c $INDIR/$fastq2)  \
--outFilterMultimapNmax 100 \
--outSAMtype BAM SortedByCoordinate  \
--outFileNamePrefix output_bams/${name}_  

#STAR --genomeLoad LoadAndKeep \ #02/04/21
#--genomeDir combined_reference/StarIndexed_combined \
#--runThreadN 14 \
#--readFilesIn <(gunzip -c $INDIR/$fastq1) <(gunzip -c $INDIR/$fastq2) \
#--sjdbGTFfile combined_reference/combined.gtf \
#--outFilterType BySJout \ #Default will be "normal". "BySJout" is for 2-pass m$
#--outFilterMultimapNmax 100 \
#--outSAMtype BAM SortedByCoordinate  \
#--outFileNamePrefix output_bams/${name}_


rm $INDIR/$fastq1
rm $INDIR/$fastq2
done

# Unload genome of memory

wait 

STAR --genomeLoad Remove --genomeDir combined_reference/StarIndexed_combined

exit
