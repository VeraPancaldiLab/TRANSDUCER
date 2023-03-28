#!/bin/bash
outdir="/home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/00_Data/240323_FASTQC/02_Output/"
fastq_regexp="/home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/00_Data/PJ2212479/RAWDATA/"
threads=12

fastqc -t $threads $(ls $fastq_regexp) -o $outdir
multiqc $outdir* -o $outdir

