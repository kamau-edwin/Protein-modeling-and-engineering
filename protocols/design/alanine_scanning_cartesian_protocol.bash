#!/bin/bash
#$ -S /bin/bash
#$ -N ala_scan
#$ -cwd
#$ -M etk243@nyumc.org
#$ -j y 
#$ -o ala_scan.log

module load openmpi/gcc/64/1.10.1

# flags below used for cartesian and ddg calculation using new score_function ref2015==beta_nov15

echo "started run with $NSLOTS processor"
echo " "

project=`echo $PWD |awk -F"/" '{print$NF}'`

echo $project

mkdir -p analysis
working_dir=`echo $PWD`

output=`echo $working_dir/analysis`
echo " "
echo $output
while read line
do
  mut=`echo $line|awk -F" " '{print$2$3$4}'`
  pose=`echo $line|awk -F" " '{print$3}'`
 
  echo "Begin interface_analyzer analysis for  $mut mutation ...."
  make_mut_file.py $line
  
  mkdir -p $mut 
  
  cartesian_ddg.mpi.linuxgccrelease -s *relaxed.pdb -ddg:mut_file $mut.file  -ddg:iterations 3 -optimization:default_max_cycles 200 -bbnbr 1 -relax:min_type lbfgs_armijo_nonmonotone -fa_max_dis 9.0 -beta_cart >$project-cartesian_ddg.out 

  mv MUT_$pose*pdb $mut
  mv WT_bj*.pdb $mut  
  mv $mut.ddg $output
  
  cd $mut
  #echo $mut
  residue_energy_breakdown.mpi.linuxgccrelease -s *pdb -out:file:silent $output/$mut-energy_breakdown.out -beta_cart > $working_dir/residue_pair_energies.log
  echo "Finished generating mutation and analysis data for $mut mutant"
  echo ''
  cd $working_dir
done < input.txt

cd $working_dir 

echo "Plotting and summarizing  project $project analysis"
echo " "

ddg_analysis.py `echo "$FUN"`  $project-`echo "$NAME"`

