#!/bin/bash

PASAHOME_PATH=$1
ALN_CONFIG=$(readlink -f $2)
ASM=$(readlink -f $3)
ORIG_GFF3=$(readlink -f $4)
CLEAN_TRANSCRIPTS=$(readlink -f $5)
PASA_DB_NAME=$6
THREADS=$7

echo "round1"

$PASAHOME_PATH/scripts/Load_Current_Gene_Annotations.dbi \
     -c $ALN_CONFIG -g $ASM \
     -P $ORIG_GFF3

ln -s $CLEAN_TRANSCRIPTS


cat $PASAHOME_PATH/pasa_conf/pasa.annotationCompare.Template.txt |\
     sed "s/<__DATABASE__>/\/tmp\/$PASA_DB_NAME/" > annotCompare.config

$PASAHOME_PATH/./Launch_PASA_pipeline.pl \
        -c annotCompare.config -A \
        -g $ASM \
        -t $CLEAN_TRANSCRIPTS --CPU $THREADS


#round1_gff=$(readlink -f *.gene_structures_post_PASA_updates.*.gff3)

mv *.gene_structures_post_PASA_updates.*.gff3 round1.gene_structures_post_PASA_updates.gff3

round1_gff=$(readlink -f round1.gene_structures_post_PASA_updates.gff3)

echo "round2"

$PASAHOME_PATH/scripts/Load_Current_Gene_Annotations.dbi \
     -c $ALN_CONFIG -g $ASM \
     -P $round1_gff


$PASAHOME_PATH/./Launch_PASA_pipeline.pl \
        -c annotCompare.config -A \
        -g $ASM \
        -t $CLEAN_TRANSCRIPTS --CPU $THREADS

mv *.gene_structures_post_PASA_updates.*.gff3 round2.gene_structures_post_PASA_updates.gff3

