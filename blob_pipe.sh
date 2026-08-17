#!/bin/bash



if [ "$#" -le 1  ]; then
        echo ''
	echo 'activate environment:'
	echo '		'
	echo '		 conda activate blobtools'
	echo ''
        echo 'usage:'
        echo '          blob_pipe.sh fasta read1 read2 ill mapping_threads blast_threads - for illumina reads'
        echo '          blob_pipe.sh fasta reads pb(or)ont(or)hifi mapping_threads blast_threads - for single ended/long PB or ONT reads'
        echo '          blob_pipe.sh fasta interleaved_reads mapping_threads blast_threads - for interleaved paired end reads'
        echo ''
        echo ''
        exit
fi

# setup 

cp ~/ncbi_dbs/taxdb.b* .

## checking for previous files

asmfile=$(readlink -f $1)

if [ "$#" == 6 ]; then
mapping_threads=$5
blast_threads=$6
elif [ "$#" == 5 ]; then
mapping_threads=$4
blast_threads=$5
elif [ "$#" == 4 ]; then
mapping_threads=$3
blast_threads=$4
fi

SAM_FILE=out_sort.bam

if [ -f "$SAM_FILE" ]; then
        echo 'alignment output found (out_sort.bam), will use this'
        sleep 3

fi

FILE=megablast.out

if [ -f "$FILE" ]; then
        echo 'blast output found, will use this'
        sleep 3

fi

if [ -f "$FILE" ] && [ -f "$SAM_FILE" ]; then
	echo 'evidence of previous run found. To re-do analysis, run "rm -r megablast.out out.sam blob/"'
	sleep 5
fi


## run megablast

if [ -f "$FILE" ]; then
        echo ''
else

blastn \
-task megablast \
-query $asmfile \
-db ~ncbi_dbs/nt/nt \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads $blast_threads \
-evalue 1e-25 \
-out megablast.out &

fi

## run read mapping

if [ -f "$SAM_FILE" ]; then
	echo ''
else

if [ "$#" == 6 ]; then

bowtie2-build $asmfile BTind

# 1. map reads

bowtie2 -p $mapping_threads -x BTind -1 $2 -2 $3 | samtools view -Sb - > out.bam 
samtools sort -@ 16 -m 4G -o out_sort.bam out.bam 
samtools index out_sort.bam

elif [ "$#" == 5 ]; then

minimap2 -I 64G -ax map-"$3" -t $mapping_threads $asmfile $2 | samtools view -Sb -  > out.bam 
samtools sort -@ 16 -m 4G -o out_sort.bam out.bam 
samtools index out_sort.bam

else

bwa index $asmfile

bwa mem -t $mapping_threads $asmfile -p $2 | samtools view -Sb - > out.bam   
samtools sort -@ 16 -m 4G -o out_sort.bam out.bam 
samtools index out_sort.bam
fi

fi


wait

## 3. run blobplot



BLOB_FILE=blob/blob.blobDB.json

if [ -f "$BLOB_FILE" ]; then
        echo 'blob plot found, skipping blobtools create '
        sleep 3

else

mkdir blob

cd blob

blobtools create \
-i $asmfile \
-b ../out_sort.bam \
-t ../megablast.out \
-o blob \
--names ~/blobtools1.1/blobtools/data/names.dmp \
--nodes ~/blobtools1.1/blobtools/data/nodes.dmp

mkdir blobplot

for RANK in species genus family order phylum superkingdom
do

blobtools plot \
-i blob.blobDB.json -r $RANK \
-o blobplot/

blobtools view -i blob.blobDB.json --hits --concoct -r $RANK

mv blob.blobDB.table.txt blob."$RANK".blobDB.table.txt
mv blob.blobDB.concoct_coverage_info.tsv blob."$RANK".blobDB.concoct_coverage_info.tsv

done


fi

echo ''
echo 'blob plot finished!'
echo ''
