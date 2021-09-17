#!/bin/bash
#01/04/21
#jacobo.solorzano@inserm.fr

INDIR="/home/jacobo/Documents/02_TRANSDUCER/04_Deconvolution_TE/03_Sigantures/cell_types/fibroblast/PJ2003085/090921_FastqProcessing/01_Input/"
length=100
# Index 
#mkdir 01_Input/StarIndexed_GRCh38.104

#STAR --runMode genomeGenerate \
#--genomeDir 01_Input/StarIndexed_GRCh38.104 \
#--genomeFastaFiles 01_Input/GRCh38.104/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa \
#--sjdbGTFfile 01_Input/GRCh38.104/Homo_sapiens.GRCh38.104.gtf \
#--sjdbOverhang $length   # set to the length of base pair sequencing 

# Load genome into memory
STAR --genomeDir 01_Input/StarIndexed_GRCh38.104 \
--genomeLoad LoadAndExit

echo $(ls $INDIR)
# Mapping
for file in $(ls $INDIR | grep "R1")
do
	echo "$file"
	name=$(echo "$file" | awk -F_R1 '{print $1}')
	fastq1="${name}"_R1.fastq.gz
	fastq2="${name}"_R2.fastq.gz

STAR --limitBAMsortRAM 10000000000 \
--genomeLoad LoadAndKeep \
--genomeDir 01_Input/StarIndexed_GRCh38.104 \
--runThreadN 14 \
--readFilesIn <(gunzip -c $INDIR/$fastq1) <(gunzip -c $INDIR/$fastq2)  \
--outSAMtype BAM SortedByCoordinate  \
--outFileNamePrefix 02_Output/output_bams/s${name}_  

#STAR --genomeLoad LoadAndKeep \ #02/04/21
#--genomeDir combined_reference/StarIndexed_combined \
#--runThreadN 14 \
#--readFilesIn <(gunzip -c $INDIR/$fastq1) <(gunzip -c $INDIR/$fastq2) \
#--sjdbGTFfile combined_reference/combined.gtf \
#--outFilterType BySJout \ #Default will be "normal". "BySJout" is for 2-pass m$
#--outFilterMultimapNmax 100 \
#--outSAMtype BAM SortedByCoordinate  \
#--outFileNamePrefix output_bams/${name}_


#rm $INDIR/$fastq1
#rm $INDIR/$fastq2
done

# Unload genome of memory

wait 

STAR --genomeLoad Remove --genomeDir 01_Input/StarIndexed_GRCh38.104

exit
