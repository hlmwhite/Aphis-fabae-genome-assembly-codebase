#!/bin/bash

THREADS=$1
AUG_OUT_DIR=$2

parallel -j $THREADS --bar --no-notice "nice ./{}" < jobs.lst

ls -lh $AUG_OUT_DIR/*.err 

cat $AUG_OUT_DIR/*gff > cat_aug.gff3

join_aug_pred.pl < cat_aug.gff3 > augustus.gff3 
