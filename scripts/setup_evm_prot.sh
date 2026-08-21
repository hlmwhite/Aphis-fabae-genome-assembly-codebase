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

PASA_ASM_NAME=$(cat $PASA_ASSEMBLIES_GFF3 | head -n 1 | cut -f 2)

$EVM_PATH/EvmUtils/misc/./cufflinks_gtf_to_alignment_gff3.pl $TRANS_GTF > cufflinks.EVM.gff3
$EVM_PATH/EvmUtils/misc/./augustus_GFF3_to_EVM_GFF3.pl $AUGUSTUS_GFF > augustus.EVM.gff3


printf "ABINITIO_PREDICTION\tAugustus\t1\n" >> weights.txt
printf "TRANSCRIPT\tCufflinks\t1\n" >> weights.txt
printf "TRANSCRIPT\t$PASA_ASM_NAME\t10\n" >> weights.txt
printf "PROTEIN\tminiprot_protAln\t1\n" >> weights.txt
printf "OTHER_PREDICTION\ttransdecoder\t5\n" >> weights.txt


#../scripts/./miniprot_gff_2_EVMgff3.pl $MINIPROT_OUT > prot_preds.gff3
../scripts/./miniprot_gff_2_EVMgff3.pl $MINIPROT_OUT | awk '$6 >= 80' > prot_preds.gff3

cat augustus.EVM.gff3 $TRANSDECODER_GFF3 > gene_preds.gff3

cat $PASA_ASSEMBLIES_GFF3 cufflinks.EVM.gff3 > trans_preds.gff3



#EVM_HOME=$EVM_PATH

#$EVM_HOME/EVidenceModeler \
#        --sample_id EVM_out \
#        --genome $ASM \
#        --gene_predictions gene_preds.gff3 \
#        --transcript_alignments trans_preds.gff3 \
#	--protein_alignments prot_preds.gff3 \
#        --segmentSize $SEG_SIZE \
#        --overlapSize $OLP_SIZE --CPU $THREADS --weights weights.txt
