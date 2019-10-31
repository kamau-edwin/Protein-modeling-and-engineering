#!/bin/bash
#$ -S /bin/bash
#$ -N pg9_fixbb
#$ -cwd
#$ -j y
#$ -pe openmpi 30-31
#$ -o log.out
module load openmpi/gcc/64/1.4.5
#mkdir -p fixbb 

mpirun -np $NSLOTS fixbb.mpi.linuxgccrelease -database $ROSETTA_DATABASE -s 3u4e.pdb  -in:file:fullatom -out:file:fullatom -out:pdb -out:path:pdb fixbb -out:suffix _MinPack -min_pack -ex1 -ex2 -extrachi_cutoff 0 -nstruct 150  -flip_HNQ -no_optH false -packing:repack_only >fix.log
cd fixbb
cluster.linuxgccrelease -s *.pdb -in:file:fullatom --cluster:radius 1.0 -cluster:population_weight 0.0 -cluster:sort_groups_by_energy -no_optH -out:pdb -database $ROSETTA_DATABASE -mute core.conformation.Conformation -mute core.scoring.etable -mute basic.io.database >& cluster.log
