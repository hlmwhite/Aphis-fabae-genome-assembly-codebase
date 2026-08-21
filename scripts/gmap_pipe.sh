#!/bin/bash

DB_NAME=$1
ASM=$2
TRANSCRIPTS=$3
THREADS=$4

gmap_build -D index/ -d $DB_NAME $ASM
gmap -D index -d $DB_NAME $TRANSCRIPTS --min-intronlength=30 --intronlength=500000 -t $THREADS -f 1 -n 0 --trim-end-exons=20 > gmap.psl
cat gmap.psl | sort -n -k 16,16 | sort -s -k 14,14 | perl -ne '@f=split; print if ($f[0]>=100)' | blat2hints.pl --source=PB --nomult --ep_cutoff=20 --in=/dev/stdin --out=hints.gff
