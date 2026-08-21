#!/bin/bash


ASM=$(readlink -f $1)
TRANSCRIPTS=$(readlink -f $2)
MAP_THREADS=$3
SORT_THREADS=$4
STRINGTIE_THREADS=$5


minimap2 -t $MAP_THREADS -ax splice:hq -uf $ASM $TRANSCRIPTS | samtools view -Sb - > aln.bam
samtools sort -@ $SORT_THREADS -o aln_sort.bam aln.bam
samtools index aln_sort.bam
stringtie/stringtie-2.2.1/./stringtie aln_sort.bam -p $STRINGTIE_THREADS -L -o out.stringtie.gtf
