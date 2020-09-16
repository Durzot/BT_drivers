#!/bin/bash

echo "working directory:" $PWD

# parameters 
path_mut='../MutSigCV/gpunit/example/data/LUSC.mutations.maf'
path_cov='../MutSigCV/gpunit/example/data/LUSC.coverage.txt'
path_gen='../MutSigCV/gpunit/example/data/genes.covariates.txt'
path_out='./results_LUSC/LUSC'
path_dict='./references/mutation_type_dictionary_file.txt';
path_ref='./references/chr_files_hg19';

matlab -nodisplay -nodesktop -r "run_MutSigCV('$path_mut', '$path_cov', '$path_gen', '$path_out', '$path_dict', '$path_ref')"
