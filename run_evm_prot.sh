#!/bin/bash

PASA_ASSEMBLIES_GFF3=$(readlink -f $1)
TRANSDECODER_GFF3=$(readlink -f $2)
TRANS_GTF=$(readlink -f $3)
AUGUSTUS_GFF=$(readlink -f $4)
EVM_PATH=$(readlink -f $5)
THREADS=$6
SEG_SIZE=$7
OLP_SIZE=$8
ASM=$(readlink -f ../reference/soft.fasta)
MINIPROT_OUT=$(readlink -f $9 )


EVM_HOME=$EVM_PATH

$EVM_HOME/EVidenceModeler \
        --sample_id EVM_out \
        --genome $ASM \
        --gene_predictions gene_preds.gff3 \
        --transcript_alignments trans_preds.gff3 \
       --protein_alignments prot_preds.gff3 \
        --segmentSize $SEG_SIZE \
        --overlapSize $OLP_SIZE --CPU $THREADS --weights weights.txt
