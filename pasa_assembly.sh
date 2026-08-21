#!/bin/bash

PASAHOME_PATH=$1
ASM=$(readlink -f $2)
TRANSCRIPTS=$(readlink -f $3)
MAX_INTRON_LENGTH=$4
STRINGENT_ALN_OLP=$5
#TRANS_GTF=$(readlink -f $6)
CPUS=$6
PASA_DB=$7



if grep -q '^@' $TRANSCRIPTS; then
	grep -A1 '^@' $TRANSCRIPTS | grep -v '^--$' | sed 's/^@/>/' | oneline_script.sh - | fold -w 80 | awk '{print $1}' > transcripts.fa
else
	cat $TRANSCRIPTS | oneline_script.sh - | fold -w 80 | awk '{print $1}' > transcripts.fa
fi

ln -s $ASM assembly.fa

ln -s ../stringtie/out.stringtie.gtf trans.gtf

$PASAHOME_PATH/bin/./seqclean transcripts.fa -c 16

cat $PASAHOME_PATH/pasa_conf/pasa.alignAssembly.Template.txt |\
     sed "s/<__DATABASE__>/\/tmp\/$PASA_DB/" > alignAssembly.config

$PASAHOME_PATH/./Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g assembly.fa -t transcripts.fa.clean \
    -T -u transcripts.fa --ALIGNERS minimap2,gmap --CPU $CPUS --TRANSDECODER --MAX_INTRON_LENGTH $MAX_INTRON_LENGTH \
    --stringent_alignment_overlap $STRINGENT_ALN_OLP --trans_gtf trans.gtf

$PASAHOME_PATH/scripts/./pasa_asmbls_to_training_set.dbi \
    --pasa_transcripts_fasta "$PASA_DB".assemblies.fasta --pasa_transcripts_gff3 "$PASA_DB".pasa_assemblies.gff3
