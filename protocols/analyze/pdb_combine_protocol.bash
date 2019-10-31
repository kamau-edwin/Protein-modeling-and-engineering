#!/bin/bash
#$ -S /bin/bash
#$ -N combine
#$ -cwd
#$ -j y
#$ -o Combine.out
combine_silent.linuxgccrelease -in:file:s ../*.pdb -out:file:silent Combined.sile -out:file:silent_struct_type binary -database $ROSETTA_DATABASE >& Combine.log
exit 0;
