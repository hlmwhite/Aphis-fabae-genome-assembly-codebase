#!/bin/bash

# usage= split_blob.sh tax_level blob_dir

# install the following:    mamba install -c bioconda ucsc-fasomerecords
# pos4 = mkbams    , then add bam as pos 5, will generate bams and fqs for each contig where reads have mapped

blob_outdir=$(readlink -f $2)

tax_level=$1

genome_fasta=$(readlink -f $3)

if [ "$keep_reads_opt" == "mkbams" ]; then
keep_reads_opt=$4

bam=$(readlink -f $5)
fi

blob_full_table="$blob_outdir"/blob."$tax_level".blobDB.table.txt

blob_sum="$blob_outdir"/blobplot/blob*"$tax_level"*txt

mkdir "$tax_level"_out

while read line
do
	mkdir "$tax_level"_out/$line
	awk -v taxID=$line '$6 == taxID' $blob_full_table > "$tax_level"_out/"$line"/"$line".table.txt
	faSomeRecords $genome_fasta <(awk -v taxID=$line '$6 == taxID' $blob_full_table | cut -f 1) "$tax_level"_out/"$line"/"$line".fasta
	if [ "$keep_reads_opt" == "mkbams" ]; then
		echo 'filtering bam file, this may take a while...'
		awk -v taxID=$line '$6 == taxID' $blob_full_table | cut -f 1 | tr "\n" " " | xargs samtools view -bh $bam > "$tax_level"_out/"$line"/"$line".bam
		samtools index "$tax_level"_out/"$line"/"$line".bam
		samtools bam2fq "$tax_level"_out/"$line"/"$line".bam | gzip - > "$tax_level"_out/"$line"/"$line".fastq.gz
	fi
done < <(cat $blob_sum | grep -v -e '^#' -e '^all' | cut -f 1)
