#!/bin/bash
#$ -S /bin/bash
#$ -N clus_pdb
#$ -cwd
#$ -j y
#$ -o cluster.out
cluster.linuxgccrelease -in:file:s ../*.pdb  -in:file:fullatom -cluster:gdtmm -cluster:radius -1 -cluster:population_weight 0.0 -cluster:sort_groups_by_energy -out:pdb -database $ROSETTA_DATABASE >& cluster.log
exit 0;
