#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe openmpi 11
#$ -M etk243@nyumc.org
#$ -j y 
#$ -o relax_cart.log
module load openmpi/gcc/64/1.10.1
#module load python/2.7.3

#clean_pdb.py *.pdb HLP 


#project=`echo $PWD |awk -F"/" '{print$NF}'`

#ls *fasta |grep -v $project.fasta|xargs rm

#for i in *HLP.pdb
#do 
 #  mv $i `echo $i|sed 's/HLP/renumbered/'`
#done

#mutation_input_file.py $project.fasta 

mkdir -p relax

mpirun -np $NSLOTS relax.mpi.linuxgccrelease -s *.pdb -ignore_unrecognized_res -relax:constrain_relax_to_start_coords -relax:ramp_constraints false -relax:cartesian -relax:min_type lbfgs_armijo_nonmonotone -ex1 -ex2 -use_input_sc -nstruct 10  -out:path:pdb relax -beta_cart >relax.log
#mpirun -np $NSLOTS relax.mpi.linuxgccrelease -s *.pdb -alternate_3_letter_codes pdb_sugar -include_sugars -glycam_pdb_format -write_pdb_link_records -beta_cart -ignore_unrecognized_res -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains -relax:ramp_constraints false -relax:cartesian -relax:min_type lbfgs_armijo_nonmonotone -ex1 -ex2 -use_input_sc -nstruct 100  -out:path:pdb relax  >relax.log

#cd relax

#grep pose *_renumbered*pdb |awk '{print$1,$NF}'|sort -nk2|head -n50 > pdb_list
#tar -zcf $project-relax_cartesian.tgz *_renumbered*pdb --remove-files 
#while read line
#do
 #   tar $project-relax_cartesian.tgz `echo $line`
#done < pdb_list  
#qsub $SCRIPTS/cartesian.bash
