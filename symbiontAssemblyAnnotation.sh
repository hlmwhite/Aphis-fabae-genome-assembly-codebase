


########################################################

### Symbiont assemblies

## symbiont assemblies were extracted from the primary assembly from hifiasm

blob_pipe.sh hifiasm.asm.hic.p_ctg.fasta m64147e_211013_175415.hifi_reads.fastq.gz hifi 32 64

# split contigs on genus based on blobtools assignment

split_blob.sh genus blob hifiasm.asm.hic.p_ctg.fasta

## Buchnera contigs (main chromosome and plasmids) were manually filtered from above outputs

cd genus_out/Buchnera

faSomeRecords Buchnera.fasta <(printf "ptg000053c\nptg000063l\nptg000098l\n") Buch.asm.fasta

bakta_db download --output bakta --type full

bakta --db bakta/db --verbose --output results_run2/ --prefix BuchAf --threads 16 Buch.asm.fasta \
  --locus-tag AFBUCH --genus Buchnera --species aphidicola

# genome statistics:
# 	Genome size: 678,006 bp
# 	Contigs/replicons: 3
# 	GC: 24.2 %
# 	N50: 634,752
# 	N ratio: 0.0 %
# 	coding density: 89.7 %

# annotation summary:
# 	tRNAs: 32
# 	tmRNAs: 1
# 	rRNAs: 3
# 	ncRNAs: 3
# 	ncRNA regions: 0
# 	CRISPR arrays: 0
# 	CDSs: 622
# 		hypotheticals: 15
# 		pseudogenes: 4
# 		signal peptides: 0
# 	sORFs: 0
# 	gaps: 0
# 	oriCs/oriVs: 1
# 	oriTs: 0


# Rickettsiella genome

cd blob/genus_out/Rickettsiella

bakta --db bakta/db --verbose --output results/ --prefix asmRk --threads 16 Rickettsiella.fasta  --locus-tag AFRIKT

# genome statistics:
# 	Genome size: 1,622,170 bp
# 	Contigs/replicons: 3
# 	GC: 39.3 %
# 	N50: 842,071
# 	N ratio: 0.0 %
# 	coding density: 90.0 %

# annotation summary:
# 	tRNAs: 46
# 	tmRNAs: 1
# 	rRNAs: 9
# 	ncRNAs: 3
# 	ncRNA regions: 1
# 	CRISPR arrays: 0
# 	CDSs: 1444
# 		hypotheticals: 209
# 		pseudogenes: 2
# 		signal peptides: 0
# 	sORFs: 0
# 	gaps: 0
# 	oriCs/oriVs: 1
# 	oriTs: 0

# pantoea genome

cd blob/genus_out/Pantoea

bakta --db bakta/db --verbose --output results/ --prefix asmPt --threads 16 Pantoea.fasta --locus-tag AFPANT

# genome statistics:
# 	Genome size: 5,356,453 bp
# 	Contigs/replicons: 3
# 	GC: 53.9 %
# 	N50: 4,195,533
# 	N ratio: 0.0 %
# 	coding density: 88.9 %

# annotation summary:
# 	tRNAs: 81
# 	tmRNAs: 1
# 	rRNAs: 22
# 	ncRNAs: 47
# 	ncRNA regions: 47
# 	CRISPR arrays: 0
# 	CDSs: 4943
# 		hypotheticals: 206
# 		pseudogenes: 21
# 		signal peptides: 0
# 	sORFs: 11
# 	gaps: 0
# 	oriCs/oriVs: 2
# 	oriTs: 0