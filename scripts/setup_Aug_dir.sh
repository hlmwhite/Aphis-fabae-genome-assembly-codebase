#!/bin/bash

ASM=$(readlink -f $1)
HINTS=$(readlink -f $2)
AUG_CONFIG_PATH=$(readlink -f $3)
#EXTCFG=$(readlink -f $3)
AUG_OUT_DIR=$4
SPECIES_NAME=$5
AUG_OLP=$6
AUG_CHNK=$7

# suggested overlap as --overlap=100000
# suggested chunk as --chunksize=1100000
# made this larger if required

# e.g = miniconda_12jan_2023/pasa/pasa/config/extrinsic/extrinsic.M.RM.PB.cfg

summarizeACGTcontent.pl $ASM > summary.out

mkdir split
splitMfasta.pl $ASM --outputpath=split

for f in split/*.split.*;
do
NAME=$(grep ">" $f)
mv $f split/${NAME#>}.fa
done


grep "bases" summary.out | awk -v hints=$HINTS '{print "split/"$3".fa\t"hints"\t1\t"$1}' > chr.lst


aug_dir=$AUG_OUT_DIR

mkdir $aug_dir

augCall="augustus --gff3=on --species=$SPECIES_NAME --alternatives-from-evidence=1 --UTR=off --extrinsicCfgFile=$AUG_CONFIG_PATH/extrinsic/extrinsic.M.RM.PB.cfg --softmasking=1"
myPrefix_="aug_split"

createAugustusJoblist.pl --sequences=chr.lst --wrap="#" --overlap=$AUG_OLP \
--chunksize=$AUG_CHNK --outputdir=$aug_dir/ --joblist=jobs.lst \
--jobprefix=$myPrefix_ --partitionHints --command "$augCall"
