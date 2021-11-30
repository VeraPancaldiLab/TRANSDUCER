#!/bin/bash
# Download Human
## Genome Fasta :
wget -O - ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz  | gunzip -c > Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa  

## Transcriptome GTF :
wget -O - ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz | gunzip -c >Homo_sapiens.GRCh37.75.gtf


# Download Mice
## Genome Fasta :  
wget -O -  ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/dna/Mus_musculus.GRCm38.75.dna_sm.primary_assembly.fa.gz | gunzip -c   >Mus_musculus.GRCm38.75.dna_sm.primary_assembly.fa

## Transcriptome GTF :
wget -O -   ftp://ftp.ensembl.org/pub/release-75/gtf/mus_musculus/Mus_musculus.GRCm38.75.gtf.gz | gunzip -c > Mus_musculus.GRCm38.75.gtf   


# Combine 
sh 01_Input/SMAP/SMAP_prepareReference.sh \
-t Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa \
-u Homo_sapiens.GRCh37.75.gtf \
-m Mus_musculus.GRCm38.75.dna_sm.primary_assembly.fa \
-s Mus_musculus.GRCm38.75.gtf  \
-o 02_Output/combined_reference


# Cleaning
rm Homo_sapiens.*
rm Mus_musculus.*

exit 0
