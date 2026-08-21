#!/bin/bash

PROTS_FILE=$(readlink -f $1)
ASM=$(readlink -f $2)
THREADS=$3
MAX_INTRON=$4

echo 'concatenate protein fastas....'
while read line
do
cat $line
done < <(cat $PROTS_FILE) > prots.faa

echo 'running metaeuk....'
metaeuk easy-predict $ASM prots.faa predsResults temp --threads $THREADS

cat predsResults.gff | awk '{print $NF}'  | tr ';' '\t' | sed 's/Target_ID=//' | cut -f 1 | sort | uniq > metaeuk.hits

#faSomeRecords <(zcat uniprot_sprot.fasta.gz | tr '|' '_' |  sed 's/sp_//' | tr '_' '\t' | awk '{print $1}') metaeuk.hits metaeuk.hits.faa
faSomeRecords prots.faa metaeuk.hits metaeuk.hits.faa

echo 'running miniprot...'
miniprot -G $MAX_INTRON -t$THREADS --gff $ASM metaeuk.hits.faa > miniprot.paf.gff

echo 'filtering miniprot on 0.8+ identity alignments (this may take a while...)'

######## - note, this is a very slow way to do this

TEMPFILE=temp.file
echo 0 > $TEMPFILE


while read line
do
echo $line | if grep -q mRNA -; then
	COUNTER=$[$(cat $TEMPFILE) + 1]
	echo $COUNTER > $TEMPFILE
fi
echo $line | if grep -q mRNA -; then
	#echo "TRUE" $COUNTER
	#printf "$line\t$counter\n"
	printf "$line\t" | paste - $TEMPFILE
else
	#echo "FALSE" $counter
	#printf "$line\t$counter\n"
	printf "$line\t" | paste - $TEMPFILE
fi

done < miniprot.paf.gff > numbered.paf.gff

awk '$3 == "mRNA"' numbered.paf.gff | tr ';' '\t' | sed 's/Identity=//' | awk '$11 >= 0.8' | awk '{print $NF}' > mRNA.to.keep

awk 'FNR==NR{a[$1];next} ($NF in a)' mRNA.to.keep numbered.paf.gff | cut -f 1-9 | grep -v '^#' > filtered.miniprot.gff


#cat miniprot.paf.gff | grep -v '##PAF' > miniprot.gff

