#!/bin/bash

# ProTree is a simple automoated pipeline for generating phylogenetic trees based on protein data sets.
#       - It should suit quick and simple phylogenetic inferences without much input, however if there are any modifications to the commands required,
#       - you may do so within this script at their respective command.

# the program is run within a conda environment, with most tools installed through conda. Some other scripts are installed through github:

# conda install orthofinder
# conda install -c bioconda gblocks
# conda install -c bioconda mafft
# faSomeRecords  -  https://github.com/santiagosnchez/faSomeRecords
# fasta_to_phylip.py  -  https://github.com/audy/bioinformatics-hacks/blob/master/bin/fasta-to-phylip
# ModelTest - https://github.com/ddarriba/modeltest
# conda install -c genomedk raxml-ng 
# conda install -c bioconda fasttree
# conda install -c bioconda rename

#  using fasttree removes the use of protest3 and raxml, is useful for quick trees, or for focusing on 4 or less organisms
#  this requires the flag 'ft' at the end, but leave blank to use prottest and raxml as a default

if [ "$#" == 0 ]; then
        echo ''
        echo 'usage: ProTree.sh <threads to use> <raxml threads> <option - "ft" >'
        echo '  - run ProTree in current working directory with protein fastas to use.'
        echo '  - all protein fastas must end in ".faa"!'
        echo ''
        echo 'install:'
        echo '  conda install orthofinder'
        echo '  conda install -c bioconda gblocks'
        echo '  conda install -c bioconda mafft'
        echo '  faSomeRecords  -  https://github.com/santiagosnchez/faSomeRecords'
        echo '  fasta_to_phylip.py  -  https://github.com/audy/bioinformatics-hacks/blob/master/bin/fasta-to-phylip'
        echo '  ModelTest - https://github.com/ddarriba/modeltest'
        echo '  conda install -c genomedk raxml-ng '
        echo '  conda install -c bioconda fasttree'
        echo '  conda install -c bioconda rename'

        exit
fi


threads=$1
tree=$3
rax_threads=$2

echo ''
echo 'running ProTree_pipe...'
echo ''
echo '   WARNING'
echo ''
echo 'ensure command is the following format: (10 seconds to cancel)'
echo ''
echo 'ProTree_pipe.sh <threads to use> <rax_threads> <option - "ft">'
echo ''
echo 'note - raxml threads should probably not exceed 20 threads, or poorer results may be obtained...'
echo ''

sleep 5
echo '(5 seconds to cancel)'
sleep 5

# 1 prepping files

if [ -f "work_dir/prep.OK" ]; then
        echo 'fastas prepped, skipping to orthofinder...'
        sleep 2
cd work_dir

else
        echo 'renaming fastas...'
        sleep 2
        mkdir work_dir

        for file in *faa ; do
                name=${file%.faa}
                awk -v var="$name" '/^>/{print ">"var "_" ++i "_"; next}{print}' "$file" > work_dir/"$name"_mod.faa
        done

        echo 'generating header files...'
        sleep 2
        cd work_dir

        for file in *faa ; do 
                org=${file%_mod.faa}
                grep $org $file > $org.headers
        done

        for file in *headers; do if [ -s $file ]; then echo $file present 
                else 
                echo "$file is empty or not present"
                touch prep_NOT.OK
        fi
        done

        if [ -s prep_NOT.OK ]; then
                echo 'missing headers file, there is an issue somewhere...   exiting'
                exit
        else
                touch prep.OK
        fi

fi

# 2 run orthofinder to infer single copy orthogroups

if [ -f "orthofinder.OK" ]; then
        echo 'found orthofinder results, skipping to alignment'
        sleep 2
else

        echo 'running orthofinder...'
        sleep 2
        orthofinder -f . -t $threads -og 2>&1 | tee output_OF.log

        single_copy_count_OF=$(grep 'There were' output_OF.log | cut -d" " -f 47)

        if [ $single_copy_count_OF -gt 20 ]; then
                touch orthofinder.OK
        elif [ $single_copy_count_OF -lt 20 ] && [ $single_copy_count_OF -gt 1 ]; then
                echo ''
                echo 'only '$single_copy_count_OF' single copy othologues found. consider reducing number of samples/species'
                echo ''
                echo 'continuing for now...'
                touch orthofinder.OK
                echo ''
                sleep 5
        elif [ $single_copy_count_OF -eq 0 ]; then
                echo "NOT_OK" > orthofinder_NOT.OK
        fi

        if [ -s orthofinder_NOT.OK ]; then
                echo ''
                echo 'not OK, there is an issue somewhere, can not find any shared single copy orthologues...   exiting'
                echo ''
                exit
        fi 

fi

cd OrthoFinder/Results_*/

# 3 running alignment - mafft

if [ -d "single_OGs" ]; then
        echo 'found single orthogroups directory...'
        sleep 2
else
        mkdir single_OGs
fi
cd single_OGs/

if [ -f "mafft.OK" ]; then
        echo 'found mafft alignments, skipping to gblocks...'
        sleep 2
else

cp ../Single_Copy_Orthologue_Sequences/*fa .

echo 'running mafft...'
sleep 2

for file in OG*fa;
do
        base=${file%.fa}
        printf "mafft --auto $file > $base.mafft.out 2> $base.mafft.err\n" >> mafft.comms
        printf "mv $base.mafft.out $base.fasta\n" >> rename_mafft.comms
done

parallel < mafft.comms 
parallel < rename_mafft.comms


for file in *fasta; do if [ -s $file ]; then 
        echo $file 'not empty' 
else 
        echo "$file is empty"
        touch mafft_NOT.OK
fi
done

if [ -s mafft_NOT.OK ]; then
        echo 'not OK, there is an empty alignment that needs investigating, there is an issue somewhere...   exiting'
        exit
else
        touch mafft.OK
fi

fi

# 4 running gblocks trimming

if [ -f "gblocks.OK" ]; then
        echo 'found gblocks trimmed alignments, skipping to alignment format conversion...'
        ortho_count=$(ll -d *.fas | wc -l)
        sleep 2

else

        echo 'running Gblocks...'
        sleep 2
        for file in *.fasta ;
        do
                Gblocks $file -t=p -p=y
        done

        for file in *fasta-gb ;
        do
        base=${file%.fasta-gb}
        printf "mv $file $base.fas\n" >> rename_gblocks.comms
        done

        parallel < rename_gblocks.comms



        echo 'performing illegal character check...'
        for file in OG*fas
        do
                #illegal_char_check=$(grep -e 'J' -e 'B' -e 'Z' $file)
                illegal_char_check=$(grep -v '>' $file | grep -e 'J' -e 'B' -e 'Z' -e 'X')
                if [ -z "$illegal_char_check" ]
                then
                        echo "no illegal characters found in $file"
                else
                        echo $file
                        echo $illegal_char_check
                        echo "removing file"
                        #sleep 2
                        rm $file
        fi
        done

                echo 'checking for empty lines in gblocks output...'
                echo '	(normally means something wrong with alignmnet)'
        for file in OG*fas
        do
                empties=$(grep -cvP '\S' $file)
                if [ "$empties" -eq "0" ]; then
                        echo "no empty lines found in $file"
                else
                        echo $file
                        echo 'empty fasta entries found, removing file'
                        sleep 2
                        rm $file
                fi
        done

        ortho_count=$(ll -d *.fas | wc -l)

        if [ $ortho_count -gt 20 ]; then
                touch gblocks.OK
        elif [ $ortho_count -lt 20 ] && [ $ortho_count -gt 1 ]; then
                echo ''
                echo 'only '$ortho_count' single copy othologues after gblocks found. consider reducing number of samples/species'
                echo ''
                echo 'continuing for now...'
                touch gblocks.OK
                echo ''
                sleep 5
        elif [ $ortho_count -eq 0 ]; then
                echo "NOT_OK" > gblocks_NOT.OK
        fi

        if [ -s gblocks_NOT.OK ]; then
                echo ''
                echo 'not OK, there is an issue somewhere, can not find any shared single copy orthologues after gblocks...   exiting'
                echo ''
                exit
        fi 

fi

# 5 concatenating alignments and converting alignment file to phyllip format

if [ -f "../../alignment.fasta" ]; then
        echo 'found concatenated fasta alignmnt, converting to phyllip...'
        cd ../../

else

        echo 'concatenating alignments...'
        sleep 2
        mkdir gblocks_fas
        cd gblocks_fas

        cat ../*fas > all.fas

        cp ../../../../*headers .

        for file in *headers ; do 
                name=${file%.headers}
                faSomeRecords.py --fasta all.fas --list $file --outfile $name.faa
        done 

        for file in *.faa ; do
                name=${file%.faa}
                grep -v '>' $file > 1.file
                printf ">$name\n" > fas_head
                cat fas_head 1.file > "$name"_all.fa
        done

        cat *_all.fa > ../../../alignment.fasta

        cd ../../../

fi

# 6 converting alignment fasta to phylip format

if [ -f alignment.phy ]; then
        echo 'found alignment.phy, skipping to model estimation...'
        sleep 2

else
        echo 'running fasta_to_phylip...'

        fasta_to_phylip.py --input-fasta alignment.fasta --output-phy alignment.phy

        if [ -s alignment.fasta ]; then
                echo "alignment.fasta not empty"
        else
                echo "alignment.fasta NOT empty, there is an issue...  exiting"
                echo ''
                exit
        fi

        if [ -s alignment.phy ]; then
                echo "alignment.fasta not empty"
        else
                echo "alignment.phy NOT empty, there is an issue...  exiting"
                echo ''
                exit
        fi

fi

# 7a generating phylogenetic tree with fastree

if [ $tree == "ft" ]
then
        echo 'running fasttree...'

        mkdir fasttree_out
        cd fasttree_out

        model="N/A"
        fasttree ../alignment.fasta > fasttree.nwk

        tree_file="fasttree.nwk"
else 

# 7b generating model with modeltest

if [ -f modeltest.out.out ]; then
        model=$(grep BIC modeltest.out.log | tail -n 1 | awk '{print $2}' )
else
        export PATH=/CGR/progs/ModelTest:$PATH

        modeltest-ng-static -i alignment.phy -d aa -p $threads -o modeltest.out -T raxml

        model=$(grep BIC modeltest.out.log | tail -n 1 | awk '{print $2}' )
fi
# 8 generating tree with raxml-ng

if [ -z "$model" ]; then
        echo ''
        echo 'BIC model is empty! something has gone wrong. check the .fasta.phy alignnments, as well as the modeltest output'
        echo 'note, this could also be due to organisms with the same name, gblocks truncates the name to first 10 chars, so this fastatophy will moan about this'
        echo ''
        exit
fi

[[ -d raxml-ng-out ]] || mkdir raxml-ng-out

cd raxml-ng-out

mkdir bs_100
cd bs_100

raxml-ng --all --msa ../../alignment.phy --model $model --prefix alignment_raxml_100 --threads $rax_threads --tree pars{10} --bs-trees 100

cd ..

raxml-ng --all --msa ../alignment.phy --model $model --prefix alignment_raxml --threads $rax_threads --tree pars{10} --bs-trees 1000

tree_file="alignment_raxml.raxml.support"

fi

echo '######################'

echo '  '

echo 'ProTree complete, the final tree with support from bootstraps (alignment_raxml.raxml.support, or fasttree.nwk) is in newick format and can now be viewed with your favourite tree viewer'

echo ' ' >> ../../../model.file

echo 'alignment based on: ' $ortho_count '  Single copy orthologues' >> ../../../model.file
echo '' >> ../../../model.file
echo 'BIC model used = '$model >> ../../../model.file

echo ' ' >> ../../../model.file

cat ../../../model.file

echo 'Final tree with support:' 

cp $tree_file ../../../

cat $tree_file

cd ../../../

tar -zcf 'ProTree_out.tar.gz' work_dir

echo ''

echo 'use " rm -r work_dir " to delete the working directory'

echo ' '
