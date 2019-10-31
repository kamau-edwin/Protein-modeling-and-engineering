#!/bin/bash
#$ -S /bin/bash
#$ -N ddg_cart
#$ -pe openmpi 4
#$ -cwd
#$ -M etk243@nyumc.org
#$ -j y 
#$ -o cartesian.log

module load openmpi/gcc/64/1.10.1

# flags below used for cartesian and ddg calculation using new score_function ref2015==beta_nov15

echo "started run with $NSLOTS processor"
echo " "
project=`echo $PWD |awk -F"/" '{print$NF}'`
start_dir=`echo $PWD`

mv input.txt relax 
cd relax

mkdir -p $project-analysis

working_dir=`echo $PWD`

output=`echo $working_dir/$project-analysis`
echo " "
echo $output
lowest_model=`grep pose *pdb |awk '{print$1,$NF}' |sort -nk2|awk '{print$1}'|sed 's/\:pose//'|head -n1`

echo $lowest_model

tar -cf $project-relax.tgz *pdb --remove-files 

tar -xf $project-relax.tgz $lowest_model 

while read line
do
  mut=`echo $line|awk -F" " '{print$2$3$4}'`
  pose=`echo $line|awk -F" " '{print$3}'`
 
  echo "Began processing $mut mutation ...."
  
  make_mut_file.py $line
  
 # echo " making $mut directory"

  mkdir -p $mut 
  
  mpirun -np $NSLOTS cartesian_ddg.mpi.linuxgccrelease -s $lowest_model -ddg:mut_file $mut.file  -ddg:iterations 3 -optimization:default_max_cycles 200 -bbnbr 1 -relax:min_type lbfgs_armijo_nonmonotone -fa_max_dis 9.0 -beta_cart >$project-cartesian_ddg.out 

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

echo " "

