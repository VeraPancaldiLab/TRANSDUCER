#!/bin/bash
# 11/04/22
#jacobo.solorzano@inserm.fr
cd ../240323_Alignment/output_bams/
find *.bam -exec java -jar ~/Downloads/picard.jar CollectRnaSeqMetrics -I {} -O /home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/00_Data/240323_PicardTools/02_Output/{}.RnaSeqMetrics --REF_FLAT /mnt/SERVER-CRCT-STORAGE/CRCT06/5To/Jacobo/00_Other/reference_genomes/refFlat.hg38.nochr.txt.gz -STRAND NONE --RIBOSOMAL_INTERVALS /home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/00_Data/240323_PicardTools/01_Input/hg38.rRNA.interval_list \;
