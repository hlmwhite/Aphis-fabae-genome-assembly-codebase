### genome, and mitochindrial assembly - the methods described for assembly and annotation are applied to both assemblies

hifiasm -t 128 --h1 afrep2-331237_R1_001.fastq.gz,331221_R1_001.fastq.gz --h2 \
    afrep2-331237_R2_001.fastq.gz,331221_R2_001.fastq.gz m64147e_211013_175415.hifi_reads.fastq.gz

awk '/^S/{print ">"$2;print $3}' hifiasm.asm.hic.hap1.p_ctg.gfa > hifiasm.asm.hic.hap1.p_ctg.fasta # alternate assembly
awk '/^S/{print ">"$2;print $3}' hifiasm.asm.hic.hap2.p_ctg.gfa > hifiasm.asm.hic.hap2.p_ctg.fasta # primary assembly


### BUSCO

for genome in hic_hap1.ctg.fasta hic_hap2.ctg.fasta
do
    outdir=${genome%.fasta}
    busco -m geno -l arthropoda -c 64 -i $genome -o "$outdir"_busco
done

# hic_hap1.ctg		C:98.6%[S:95.8%,D:2.8%],F:0.5%,M:0.9%,n:1013    alternate assembly
# hic_hap2.ctg		C:98.7%[S:95.7%,D:3.0%],F:0.5%,M:0.8%,n:1013    primary assembly


### blobtools for contamination check/symbiont presence

#1. generate taxonomy hits with blast

blastn \
-task megablast \
-query $asmfile \ # where $asmfile is the priamry or alternate assembly fasta
-db CGR/progs/ncbi_dbs/nt/nt \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads 64 \
-evalue 1e-25 \
-out megablast.out &

#2. generate coverage info with minimap2
minimap2 -I 64G -ax map-"$3" -t 32 $asmfile $2 | samtools view -Sb -  > out.bam 
samtools sort -@ 16 -m 4G -o out_sort.bam out.bam 
samtools index out_sort.bam

#3. run blobtools
mkdir blob

cd blob

blobtools create \
-i $asmfile \
-b ../out_sort.bam \
-t ../megablast.out \
-o blob \
--names blobtools1.1/blobtools/data/names.dmp \
--nodes blobtools1.1/blobtools/data/nodes.dmp

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


### YaHS Hi-C scaffolding

echo "### Step 0: Index reference" 
bwa index -a bwtsw assembly.fasta

echo "### Step 1.A: FASTQ to BAM (1st)"
mkdir raw_bams
bwa mem -t 64 assembly.fasta all.r1.fq.gz | samtools view -@ 16 -Sb - > raw_bams/1.bam &

echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t 64 assembly.fasta all.r2.fq.gz | samtools view -@ 16 -Sb - > raw_bams/2.bam
#
wait 

echo "### Step 2.A: Filter 5' end (1st)"
mkdir filt_bams
samtools view -h raw_bams/1.bam | perl HiC/tools/arima/mapping_pipeline/filter_five_end.pl | \
    samtools view -Sb - > filt_bams/filt_1.bam &

echo "### Step 2.B: Filter 5' end (2nd)"
samtools view -h raw_bams/2.bam | perl HiC/tools/arima/mapping_pipeline/filter_five_end.pl | \
    samtools view -Sb - > filt_bams/filt_2.bam
#
wait

echo "### Step 3A: Pair reads & mapping quality filter"
SAMTOOLS=$(which samtools)
perl HiC/tools/arima/mapping_pipeline/two_read_bam_combiner.pl filt_bams/filt_1.bam filt_bams/filt_2.bam $SAMTOOLS 10 | \
    samtools view -bS -t assembly.fasta.fai - | samtools sort -@ 16 -o temp.bam -

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar picard-2.8.2.jar AddOrReplaceReadGroups \
    INPUT=temp.bam OUTPUT=final.bam ID=hic_R LB=hic_R SM=A_fabae PL=ILLUMINA PU=none

# scaffolding with YaHS tool

samtools index final.bam

samtools faidx assembly.fasta

mkdir yahs
cd yahs

yahs ../assembly.fasta ../final.bam  


### fixing broken chromosome in alternate assembly based on synteny to the alternate haplotype chromsome from primary

# pull out scaffolds of interest from both assemblies
faSomeRecords ../../hap2/hicExplorer/asm.fasta <(echo 'scaffold_1') scaff1_hap2.fa

faSomeRecords ../../hap1/hicExplorer/asm.fasta <(printf "scaffold_3\nscaffold_5") scaff3_5_hap1.fa

# run ragtag to scaffold together scaff 3 and 5 from the alternate haplotype
ragtag.py scaffold --aligner miniconda3/envs/ragtag/bin/nucmer -o scaff3_5.hap1.scaff --nucmer-params='--maxmatch -l 100 -c 500 -t 16' scaff1_hap2.fa scaff3_5_hap1.fa

# add in the fixed scaffold to the alternate assembly
faSomeRecords -exclude ../../hap1/hicExplorer/asm.fasta <(printf "scaffold_3\nscaffold_5") no_scaf3_5.hap1.fasta 

cat scaff3_5.hap1.scaff/ragtag.scaffold.fasta no_scaf3_5.hap1.fasta |  awk '/^>/{print ">scaffold_" ++i ; next}{print}' - | \
    fold -w 80 - > fixed.hap1.fa



############################################


# mitochondrial assembly


python tools/MitoHiFi/findMitoReference.py --species "Aphis fabae" --email some@email.ac.uk --outfolder data/ --min_length 14000
 
cd mitohifi_run/
 
# assembly with mitohifi

python ../tools/MitoHiFi/mitohifi.py -r m64147e_211013_175415.hifi_reads.fastq.gz -f ../data/NC_039988.1.fasta -g ../data/NC_039988.1.gb -o 5 -t 32

# annotation with mitofinder

singularity exec mitofinder_v1.4.2.sif  \
  mitofinder -j AFmit -a ../final_mitogenome.fasta -r ../NC_039988.1.gb -o 5 -p 32



