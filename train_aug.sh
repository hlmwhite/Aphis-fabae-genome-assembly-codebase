#!/bin/bash

CDS=$(readlink -f $1)
GFF3=$(readlink -f $2)
ASM=$(readlink -f $3)
THREADS=$4
AUG_CONFIG_PATH=$5
SPECIES_NAME=$6

grep complete $CDS | perl -pe 's/>(\S+).*/$1/' > complete.orfs

grep -F -f complete.orfs $GFF3 | grep -P "(\tCDS\t|\texon\t)" | perl -pe 's/cds\.//; s/\.exon\d+//;' > trainingSetComplete.temp.gff3

cat trainingSetComplete.temp.gff3 | perl -pe 's/\t\S*(asmbl_\d+).*/\t$1/' | sort -n -k 4 | sort -s -k 9 | sort -s -k 1,1 > trainingSetComplete.gff3

ln -s trainingSetComplete.gff3 bonafide.gtf

computeFlankingRegion.pl bonafide.gtf > flank_region.txt

FLANK_VAL=$(cat flank_region.txt | tail -n 1 | awk '{print $5}')


gff2gbSmallDNA.pl bonafide.gtf $ASM $FLANK_VAL bonafide.gb

cat bonafide.gb | perl -ne 'if(m/\/gene=\"(\S+)\"/){ print "\"".$1."\"\n";}' | sort -u > traingenes.lst

grep -f <(cat traingenes.lst | sed 's/"//g' ) -F bonafide.gtf > bonafide.f.gtf

mkdir split
cd split

split -n l/$THREADS ../bonafide.f.gtf

IFS=

for file in x*; do while read line; do asmID=$(echo $line | awk '{print $NF}'); echo $line | cut -f 1-8 | awk -v ID=$asmID '{print $0"\tgene_id "ID"; transcript_id "ID"; gene_name "ID";"}' ;  done < $file > "$file"_temp.gtf & done

wait

#######

for file in *_temp.gtf
do
cat $file
done > mod_bonafide.f.gtf

head mod_bonafide.f.gtf > test.file

if grep -q 'ID=' test.file; then
    exit
else
    sed -i.bak -e 's/gene_id /gene_id ID=/' -e 's/transcript_id /transcript_id ID=/' -e 's/gene_name /gene_name ID=/' mod_bonafide.f.gtf 
fi

cd ..

gffread split/mod_bonafide.f.gtf -g $ASM -y prots.faa

cd-hit -i prots.faa -o 0.8.prots.faa -c 0.8

grep ">" 0.8.prots.faa | perl -pe 's/>//' > nonred.lst

cat bonafide.gb | perl -ne '
if ( $_ =~ m/LOCUS\s+(\S+)\s/ ) {
$txLocus = $1;
} elsif ( $_ =~ m/\/gene=\"(\S+)\"/ ) {
$txInGb3{$1} = $txLocus
}
if( eof() ) {
foreach ( keys %txInGb3 ) {
print "$_\t$txInGb3{$_}\n";
}
}' > loci.lst

grep -f nonred.lst loci.lst | cut -f2 > nonred.loci.lst

filterGenesIn.pl nonred.loci.lst bonafide.gb > bonafide.f.gb

grep -c LOCUS bonafide.gb bonafide.f.gb

AUGUSTUS_CONFIG_PATH=$AUG_CONFIG_PATH

new_species.pl --species=$SPECIES_NAME

etraining --species=$SPECIES_NAME bonafide.f.gb &> bonafide.out

### see notes about stopCodonExcludedFromCDS true, should be ok here, but worth keeping note of

etraining --species=$SPECIES_NAME bonafide.f.gb 2>&1 | grep "in sequence" | perl -pe 's/.*n sequence (\S+):.*/$1/' | sort -u > bad.lst

filterGenes.pl bad.lst bonafide.f.gb > bonafide.F.gb

randomSplit.pl bonafide.F.gb 200
mv bonafide.F.gb.test test.gb
mv bonafide.F.gb.train train.gb

etraining --species=$SPECIES_NAME train.gb &> etrain.out

augustus --species=$SPECIES_NAME test.gb > test.out

AMBER=$(tail -6 etrain.out | head -3 | awk '{print $NF}' | head -n 1 | sed -e 's/(//' -e 's/)//')
OCHRE=$(tail -6 etrain.out | head -3 | awk '{print $NF}' | tail -n 2 | head -n 1 | sed -e 's/(//' -e 's/)//')
OPAL=$(tail -6 etrain.out | head -3 | awk '{print $NF}' | tail -n 1 | sed -e 's/(//' -e 's/)//')

sed -i.bak -e "/^\/Constant\/amberprob/ s/0.33/$AMBER/" \
   -e "/^\/Constant\/ochreprob/ s/0.33/$OCHRE/" \
   -e "/^\/Constant\/opalprob/ s/0.34/$OPAL/" $AUGUSTUS_CONFIG_PATH/species/$SPECIES_NAME/*parameters.cfg 
