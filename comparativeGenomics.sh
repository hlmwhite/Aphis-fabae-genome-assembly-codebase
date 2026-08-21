
#### Comparative genomics work


### phylogentic analysis

# protein faa files were taken from assembly accessions liste din supplementary table 2 ( see ProTree_v2.sh script)

ProTree_v2.sh 128 32

###########################

# synteny analysis with A. gossypii

# generate a file of busco protein sequences from BUSCo outputs for each genome
for ORG in Agoss Afabae
do

for FASTA in "$ORG"_chr_busco/run_arthropoda_odb10/busco_sequences/single_copy_busco_sequences/*faa 
do
    BASE=$(basename $FASTA)
    FA_NAME=${BASE%.faa}
    cat $FASTA | sed "/>/ s/$/_$FA_NAME/"
done > "$ORG".busco.faa

done

# blast comparison between the two protein BUSCO sets
makeblastdb -in Afabae.busco.faa -dbtype prot -out Afabae_proteins_db

blastp -query Agoss.busco.faa -db Afabae_proteins_db -outfmt 6 -evalue 1e-10 -num_threads 16 -out output.blast

cat output.blast | awk '$1 != $2' | awk '$3 >= 60' |  sort -k1,1 -k11,11g | awk '{
  if (seen[$1] < 5) {
        print $0;
        seen[$1]++;
    }
}' | sed -e 's/-/_/g' -e 's/:/_/g' | tr ' ' '\t' > top.blast

# generate gff format of BUSCO genes, and run synteny analysis with mscanx

cat Afabae_chr_busco/run_arthropoda_odb10/full_table.tsv Agoss_chr_busco/run_arthropoda_odb10/full_table.tsv | grep -i 'Complete' | awk '{print $3"\t"$3"_"$4"_"$5"_"$1"\t"$4"\t"$5}' > top.gff

/pub65/markw/CGR/projects/RnD/whatshap/BUSCO_hifiasm/nhap4_s0.15_homCov_84_solanales_busco_old/mcscanx/MCScanX/./MCScanX ./top


###########################


##### orthologue analysis between four aphid species

# use agat to keep longest isoform for each species. for example:

agat_sp_keep_longest_isoform.pl --gff ../annotate/pasa_update/round2.gene_structures_post_PASA_updates.gff3 -o longest_iso.update.gff3

miniconda_12jan_2023/pasa/pasa/opt/pasa-2.5.2/./misc_utilities/gff3_file_to_proteins.pl longest_iso.update.gff3 ../annotate/reference/soft.fasta prot > A_fabae.faa


cp Aphid_downloads/A_glycines/Aglycines_longest/longest.Aglycines.faa A_glycines.faa
cp Aphid_downloads/A_gossypii/longest_Agossypii/longest.Agossypii.faa A_gossypii.faa
cp Aphid_downloads/R_maidis/Rmaidis_longest/longest.Rmaidis.faa R_maidis.faa
 
# run orthofinder to identify orthologues
for file in *faa
do
	name=${file%.faa}
	sed -i.bak "s/>/>$name/" $file
done
 
orthofinder -f . -t 128
 
# generate orthogroups as list for R plots
for file in *faa ; do
	name=${file%.faa}
	echo $name > ID.file
	grep "$name" OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv | cut -f 1 | cat ID.file - | grep -v Orthogroup  | sed '1d' >  "$name".list
done
 
#Afabae.list  Aglycines.list  Agossypii.list  Rmaidis.list
 
# in R
library("dplyr")
library("ggvenn")
  
x2 <- list("A. glycines" = readLines('A_glycines.list'), "A. gossypii" = readLines('A_gossypii.list'), "R. maidis" = readLines('R_maidis.list'), "A. fabae" = readLines('A_fabae.list'))
 
ggvenn(x2, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE) 


###########################

#### detoxification gene counts

# generate interproscan results across the aphid species on longest isoforms per gene only
/pub65/markw/bin/my_interproscan/interproscan-5.32-71.0/./interproscan.sh -cpu 32 -i A_gossypii.faa -b A_gossypii.ips -dp --goterms &
/pub65/markw/bin/my_interproscan/interproscan-5.32-71.0/./interproscan.sh -cpu 32 -i A_glycines.faa -b A_glycines.ips -dp --goterms &
/pub65/markw/bin/my_interproscan/interproscan-5.32-71.0/./interproscan.sh -cpu 32 -i R_maidis.faa -b R_maidis.ips -dp --goterms &
/pub65/markw/bin/my_interproscan/interproscan-5.32-71.0/./interproscan.sh -cpu 32 -i longest.Afabae.faa -b A_fabae.ips -dp --goterms &
#
 
# generate table of counts for each detox gene of interest
for file in *tsv
do
	name=${file%.ips.tsv}
	echo $name
	grep -e 'IPR001128' $file | cut -f 1 | sort | uniq | wc -l
	grep -e 'IPR004045' -e 'IPR004046' $file | cut -f 1 | sort | uniq | wc -l
	grep -e 'IPR002018' $file | cut -f 1 | sort | uniq | wc -l
	grep -e 'IPR002213' $file | cut -f 1 | sort | uniq | wc -l
	grep -e 'IPR003439' $file | cut -f 1 | sort | uniq | wc -l
done | paste - - - - - - | cat <(printf "function\tcP450\tGSTs\tCbxyl-est\tUDPgs\tABCts\n") - | sed '/^$/d' | \
awk '{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {  
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' | tr ' ' '\t' 
 
# function        A_fabae A_glycines      A_gossypii      R_maidis
# cP450   56      58      56      47
# GSTs    11      11      9       9
# Cbxyl-est       29      25      21      18
# UDPgs   65      49      51      41
# ABCts   69      70      71      67


