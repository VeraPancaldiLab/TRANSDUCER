#!/bin/bash
#20/05/21
#jacobo.solorzano@inserm.fr

# Load genome into memory
STAR --genomeDir combined_reference/StarIndexed_combined \
--genomeLoad LoadAndExit

# Mapping

STAR --limitBAMsortRAM 10000000000 \
--genomeLoad LoadAndKeep \
--genomeDir combined_reference/StarIndexed_combined \
--runThreadN 14 \
--readFilesIn PDAC221Tkf1_fwd_P50.fq  \
--outFilterMultimapNmax 100 \
--outSAMtype BAM SortedByCoordinate  \
--outFileNamePrefix 50GC.bam


STAR --limitBAMsortRAM 10000000000 \
--genomeLoad LoadAndKeep \
--genomeDir combined_reference/StarIndexed_combined \
--runThreadN 14 \
--readFilesIn  PDAC221Tkf1_fwd_P60.fq \
--outFilterMultimapNmax 100 \
--outSAMtype BAM SortedByCoordinate  \
--outFileNamePrefix 60GC.bam 


# Unload genome of memory

wait 

STAR --genomeLoad Remove --genomeDir combined_reference/StarIndexed_combined

exit
