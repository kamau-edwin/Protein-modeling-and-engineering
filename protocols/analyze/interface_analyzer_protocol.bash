#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -M etk243@nyumc.org
#$ -j y 
#$ -o inter.log

module load openmpi/gcc/64/1.10.1

# flags below used for cartesian and ddg calculation using new score_function ref2015==beta_nov15

echo "started run with $NSLOTS processor"
echo " "

#launch in relax directory or move to it 
project=`echo $PWD |awk -F"/" '{print$8}'`

echo $project

working_dir=`echo $PWD`

while read line
do
  mut=`echo $line|awk -F" " '{print$2$3$4}'`
  echo "Begin interface_analyzer analysis for  $mut mutation ...."
  cd $mut
  InterfaceAnalyzer.mpi.linuxgccrelease -s *.pdb -interface HL_P -tracer_data_print false -compute_packstat true -packstat::oversample 5 -out:file:score_only -out:file:scorefile $project-interface.sc -pack_input -pack_separated -ex1 -ex2 -use_input_sc -inter_group_neighbors_cutoff 6.0 -score:weights interface -overwrite >interface.log 
  awk '{print$NF","$6","$9","$8","$10","$7}' *sc|sed "s/^[M].*0001/$mut/"|sed 's/^W.*0001/WT/'|sed '1d' >$project-$mut-ana.sc
  cd $working_dir
done < input.txt

cd $working_dir/ala_scan

while read line
do
  mut=`echo $line|awk -F" " '{print$2$3$4}'`
  echo "Begin interface_analyzer analysis for  $mut mutation ...."
  cd $mut
  InterfaceAnalyzer.mpi.linuxgccrelease -s *.pdb -interface HL_P -tracer_data_print false -compute_packstat true -packstat::oversample 5 -out:file:score_only -out:file:scorefile $project-interface.sc -pack_input -pack_separated -ex1 -ex2 -use_input_sc -inter_group_neighbors_cutoff 6.0 -score:weights interface -overwrite >interface.log 
  awk '{print$NF","$6","$9","$8","$10","$7}' *sc|sed "s/^[M].*0001/$mut/"|sed 's/^W.*0001/WT/'|sed '1d' >$project-$mut-alascan.sc
  echo "Finished analyzing the interface metrics for $mut mutant"
 # echo ''
  cd $working_dir/ala_scan
done < *.input

cd $working_dir

ddg_analysis.py ddg_bind $project
