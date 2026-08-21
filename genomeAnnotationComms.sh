
### annotation - commands are representative and were used for annotation of both haplotype assemblies
 
# repeat mask the assembly with repeat masker

BuildDatabase -name $db_name -engine ncbi $genome

RepeatModeler -engine ncbi -database $db_name -pa $threads

RepeatMasker -lib RM*/consensi.fa.classified -pa $threads $genome -nolow -no_is -xsmall

cat "$genome".out | awk '{print $5"\t"$6"\t"$7}' | tail -n +4 > "$db_name".bed
 
ln -s asm.fasta.masked soft.fasta


# make directories for annotation and prepare hq transcripts

mkdir stringtie pasa_asm train_augustus gmap run_augustus prot_mapping evm pasa_update reads TE_blast reference

cat CGR/projects/RnD/isoAno/Afabae_hap2/prots.list > prots.file
 

cat ../../hq_transcripts.fasta | tr '/' '_' > reads/fixed.transcripts.fa
cd ..
 
#### STRINGTIE - map transcripts to the reference in prep for pasa assembly
cd stringtie
../scripts/./stringtie.sh ../reference/soft.fasta ../reads/fixed.transcripts.fa 32 4 32

#### GMAP - map transcripts to the reference in prep for pasa assembly
cd ../gmap
 
../scripts/./gmap_pipe.sh gmap_fabae_Hap2 ../reference/soft.fasta ../reads/fixed.transcripts.fa 48

### MAP PROTEINS - map proteins to the reference in prep for evidence modeler
cd ../prot_mapping
 
../scripts/./map_proteins.sh ../prots.file ../reference/soft.fasta 32 250000 &


#### PASA ASSEMBLY - generate high quality transcript assemblies for training augustus
cd ../pasa_asm
 
../scripts/./pasa_assembly.sh miniconda_12jan_2023/pasa/pasa/opt/pasa-2.5.2 ../reference/soft.fasta ../reads/fixed.transcripts.fa 250000 30.0 32 Afabae_hap2_pasa
 
#### TRAIN AUGUSTUS - train a genome specific model for A. fabae based on pasa assemblies
cd ../train_augustus
 
../scripts/./train_aug.sh ../pasa_asm/*.assemblies.fasta.transdecoder.cds ../pasa_asm/*.assemblies.fasta.transdecoder.genome.gff3 ../reference/soft.fasta 100 miniconda_12jan_2023/pasa/pasa/config Afabae_hap2_16jul
  
### RUN AUGUSTUS - run augustus gene finder with traine genome model
cd ../run_augustus
 
../scripts/./setup_Aug_dir.sh ../reference/soft.fasta ../gmap/hints.gff miniconda_12jan_2023/pasa/pasa/config aug_out Afabae_hap2_16jul 100000 1100000
 
../scripts/./run_augustus_parallel.sh 96 aug_out  

### EVM - generate consensus gene models
cd ../evm
 
../scripts/./setup_evm_prot.sh ../pasa_asm/*.pasa_assemblies.gff3 ../pasa_asm/*.assemblies.fasta.transdecoder.genome.gff3 ../pasa_asm/trans.gtf ../run_augustus/augustus.gff3 EVM/EVidenceModeler-v2.0.0 48 100000 10000 ../prot_mapping/filtered.miniprot.gff
 
../scripts/./run_evm_prot.sh ../pasa_asm/*.pasa_assemblies.gff3 ../pasa_asm/*.assemblies.fasta.transdecoder.genome.gff3 ../pasa_asm/trans.gtf ../run_augustus/augustus.gff3 EVM/EVidenceModeler-v2.0.0 48 100000 10000 prot_mapping/filtered.miniprot.gff

### TE IDENTIFCATION - remove gene models that are mostly transposon elements
cd ../TE_blast
 
../scripts/./TE_ident.sh progs/transposon_psi/t_psi_ref
  
 
#### PASA UPDATE - update final gene models with UTRs
cd ../pasa_update
 
../scripts/./update_script.sh miniconda_12jan_2023/pasa/pasa/opt/pasa-2.5.2 ../pasa_asm/alignAssembly.config ../reference/soft.fasta ../TE_blast/noTE.EVM.gff3 ../pasa_asm/transcripts.fa.clean Afabae_hap2_pasa 12


## generate protein and CDS files for predictions

miniconda_12jan_2023/pasa/pasa/opt/pasa-2.5.2/./misc_utilities/gff3_file_to_proteins.pl round2.gene_structures_post_PASA_updates.gff3 ../reference/soft.fasta prot > rnd2.update.pep

miniconda_12jan_2023/pasa/pasa/opt/pasa-2.5.2/./misc_utilities/gff3_file_to_proteins.pl round2.gene_structures_post_PASA_updates.gff3 ../reference/soft.fasta CDS > rnd2.update.CDS


#gff_stats.py <(cat round2.gene_structures_post_PASA_updates.gff3 | grep -v -e '#' | sed '/^$/d')
63      18874   8120.9617463176855      177773  264.23724637599634      46974248

### functional annotation - 

### eggnog

 
export EGGNOG_DATA_DIR=progs/eggnog_dbs
emapper.py -i pep.fasta -o eggnog_out --cpu 48

less -S eggnog_out.emapper.annotations | grep -v '#' | cut -f 1 | sort | uniq | wc -l
# 21367  ## number of transcipts/peps annotated

# signalp

export PATH=bin/signalp/signalp-5.0b/bin:$PATH
  
signalp -stdout -org euk -format short -fasta pep.fasta > signalp.out

## interproscan

bin/my_interproscan/interproscan-5.32-71.0/./interproscan.sh -cpu 128 -i pep.faa -b A_fabae.ips -dp --goterms


