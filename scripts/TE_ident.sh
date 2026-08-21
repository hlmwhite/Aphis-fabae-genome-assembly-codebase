#!/bin/bash

T_PSI_REF_PATH=$(readlink -f $1)

diamond blastp --sensitive --query ../evm/EVM_out.EVM.pep --out repeats.out -d $T_PSI_REF_PATH --evalue 1e-10 --max-target-seqs 1 --outfmt 6 --header

grep -v -f <(cat repeats.out | grep -v '#'  | cut -f 1 | sed -e 's/evm.model.//' -e 's/$/;/') ../evm/EVM_out.EVM.gff3 | grep -v -f <(cat repeats.out | grep -v '#'  | cut -f 1 | sed -e 's/evm.model.//' -e 's/$/\.exon/') - > noTE.EVM.gff3


