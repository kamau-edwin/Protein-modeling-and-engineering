#!/usr/bin/python
########################################################################
# File generator for comparative modeling protocol.py
#
# Last modified: 9.08.2015 by Edwin Kamau
#
# Usage: python Template.py project template  target_sequence template_sequence
########################################################################
# import all the necessary functions required to make functional calls
import string, sys, re, os
from subprocess import call
import random

if len(sys.argv)!=9:
    print 'Usage: python LoopKICwfragments.py <project> <Template> <target_sequence> <template_sequence> <start> <end> <nstruct> <native--None or path to native>' 
    sys.exit(0)
#global variables 
pwd = os.getcwd()
pdir = sys.argv[1]  # name of the job obtained from item 2 of usage(see above)
PDB = sys.argv[2]
start = int(sys.argv[5])
end = int(sys.argv[6])
nstruct=int(sys.argv[7])
native=sys.argv[8]

for i in sys.argv[3:4]: # item 3 of usage (see heading above)
    try: 
        target_sequence=i
    except ValueError:
        print 'alignment <target_sequence>'
        sys.exit(0)
for i in sys.argv[4:5]: # item 4 of usage (see heading above)
    try: 
        template_sequence=i
    except ValueError:
        print 'alignment <template_sequence>'
        sys.exit(0)
    else:
        print 'Aligning with clustalo' 
        f = open('./temp.fasta', 'w') 
        f.write('>Target\n%s\n>Template1\n%s' % (target_sequence, template_sequence))
        f.close()

#Run the alignment
return_code = call("clustalo -i temp.fasta -o temp.aln --outfmt=clu --force", shell=True)  # maybe True works better?
#return_code = call("clustalw2 temp.fasta",shell=True)

#Parse the alignment file
file = open("./temp.aln", 'r')
full_target_sequence = ''
full_template_sequence = ''

while 1:
    line = file.readline() # while loop that reads temp.aln line by line 
    if not line:
        break # stop while loop if no line is found
    if line.find('Target')==0: # find object in position one
        reline = re.sub(' +',' ',line) # substitute any space with single space
        reline_split = reline.split( ) #split reline on the space
        full_target_sequence += reline_split[1] #append position 2 of the reline_split to file full target ... 
    if line.find('Template1')==0:
        reline = re.sub(' +',' ',line)
        reline_split = reline.split( )
        full_template_sequence += reline_split[1]        

#create project directory std output/error and 
if not os.path.exists(pdir): # project directory path
       os.makedirs(pdir)
if not os.path.exists(pdir + '/' + pdir + '.out'): #loopmodeling output folder
       os.makedirs(pdir + '/'+ pdir + '.out')
if not os.path.exists(pdir + '/frags'): # Rosetta fragment file output folder
       os.makedirs(pdir + '/frags')
if not os.path.exists(pdir + '/' + pdir + '.out/Combined'): # Rosetta output analysis folder
       os.makedirs(pdir + '/'+ pdir +'.out/Combined')

f = open(pdir + '/'+  pdir + '.ali', 'w') # will create an alignment file used in modelling 
f.write('>%s\n%s\n>%s\n%s\n>pymol\n%s' % (pdir, full_target_sequence, PDB, full_target_sequence,full_template_sequence))
f.close()

# fasta file of the target sequence passed to the flag file 
f = open(pdir + '/' + pdir + '.fasta', 'w')
f.write('>'+ pdir + '\n%s\n\n' % (target_sequence))
f.close()

#print pdir + ' ' + full_target_sequence   
print PDB + '.pdb\n' + full_template_sequence
print ''
print pdir + '\n' + full_target_sequence      
print ' '
f = open(pdir + '/' + pdir + '.loops','w')
#for m in re.finditer('-+',str(full_target_sequence)):
#	f.write('LOOP %i %i %i  \n' % (m.start()-1,m.end(),m.end()-1))
for m in re.finditer('-+',str(full_template_sequence)):
        f.write('LOOP %i %i %i 0 1\n' % (m.start(),m.end()+1,m.end()))
	print 'LOOP %i %i %i 0 1' % (m.start(),m.end()+1,m.end())	
return_code = call("rm -f temp.*", shell=True)

# Create flag files
print 'Creating flag files...'
f = open(pdir + '/' + pdir + '.flags', 'w') # write below to pdir.flags in jobname folder
f.write('-s %s.pdb\n' %pdir)
f.write('-in:file:fullatom\n')
#f.write('-loops:frag_sizes 9 3 1\n')
#f.write('-loops:frag_files frags/' + pdir[:5] + '.200.9mers frags/' + pdir[:5] + '.200.3mers none\n')
f.write('-loops:loop_file ' + pdir + '.loops\n')
#f.write('-loops:remodel perturb_kic_with_fragments\n')
f.write('-loops:remodel perturb_kic\n')
#f.write('-loops:refine refine_kic_with_fragments\n')
f.write('-loops:refine refine_kic\n')
f.write('-out:nstruct %i\n' % nstruct)
f.write('-out:file:silent '+ pdir + '.out/' + pdir + '.silent.${SGE_TASK_ID}\n')
f.write('-out:file:fullatom\n')
f.write('-out:file:silent_struct_type binary\n')
f.write('-run:seed_offset ${SGE_TASK_ID}\n')
f.write('-kic_bump_overlap_factor 0.36\n')
f.write('-legacy_kic false\n')
f.write('-kic_min_after_repack true\n')
f.write('-corrections:score:use_bicubic_interpolation false\n')
#NGK flags
f.write('-cycles_outer 5\n')
f.write('-loops:kic_rama2b\n')
f.write('-loops:kic_omega_sampling\n')
f.write('-allow_omega_move\n')
f.write('-loops:ramp_fa_rep\n')
f.write('-loops:ramp_rama\n')
f.write('-ex1\n')
f.write('-ex2\n')
f.write('-extrachi_cutoff 0\n')
f.write('-overwrite\n')
f.write('-mute core.io.database\n')
f.write('-mute basic.io.database\n')
f.write('-mute protocols.looprelax.FragmentPerturber\n') 
f.write('-mute protocols.loops.loop_mover.refine.LoopMover_Refine_KIC\n') 
f.write('-mute core.fragments.ConstantLengthFragSet\n')
f.write('-mute core.conformation.conformation\n')
f.write('-mute core.scoring.etable\n')
f.write('-mute core.chemical.ResidueTypeSet\n')
f.write('-mute protocols.loop_modeling.loggers.ProgressBar\n')
f.close()
#Create a shell script bash file to generate fragments. 
f = open(pdir + '/frags/frags.bash', 'w')
f.write('#!/bin/bash\n') 
f.write('#$ -S /bin/bash\n') #use this shell to interpret script
f.write('#$ -N F_%s\n' % pdir) #name of the job
f.write('#$ -cwd\n') # run script in current working directory
f.write('#$ -M etk243@nyumc.org\n') #send status updates
f.write('#$ -j y \n') #output SGE input and output into the same file
f.write('#$ -o %s/frags/log.out\n' % pdir) # output file
f.write('module load openmpi/gcc/64/1.4.5\n')
f.write('\n') 
f.write('cd %s/frags\n' % pdir)
f.write('#make_fragments.pl -id %s -frag_sizes 9,3 ../%s.fasta\n' %(pdir,pdir))
f.write("#ls *.* |egrep -v 'bash|out|%s.200.*mers|%s.psipred_ss2' |xargs rm\n" %(pdir[:5],pdir[:5]))
f.write('cd ..\n')
f.write('thread_pdb_from_alignment.py --template=%s --target=%s --chain=G --align_format=fasta %s.ali $INPUT/%s.pdb %s.pdb\n' %(PDB,pdir,pdir,PDB,pdir))
f.write('qsub seed.bash\n')
f.close()
# create random number for seeding
f = open(pdir + '/seed', 'w')
for i in range(0,10):
   f.write('%s\n'%random.randrange(start,end))
f.close()
# create random seed generator
f = open(pdir + '/seed.bash', 'w')
f.write('#!/bin/bash\n') 
f.write('#$ -S /bin/bash\n') #use this shell to interpret script
f.write('#$ -cwd\n') # run script in current working directory
f.write('#$ -M etk243@nyumc.org\n') #send status updates
f.write('#$ -j y \n') #output SGE input and output into the same file
f.write('#$ -o log.out\n') # output file
f.write('\n') 
f.write('while read line\n')
f.write('do\n')
f.write("  seed=`echo $line`\n")
f.write('  sed -i "7s/t[.].*/t.$seed/" *.flags\n')
#f.write('  sed -i "9s/t[.].*/t.$seed/" *.flags\n')
f.write('  sed -i "10s/ .*/ $seed/" *.flags\n')
#f.write('  sed -i "24s/ .*/ $seed/" *.flags\n')
f.write('  qsub loopmodeling.bash\n')
f.write('  sleep 30\n')
f.write('done < seed\n')
f.close()
call('chmod +x '+ pdir + '/seed.bash',shell=True)
# loopmodeling launch
f = open(pdir + '/loopmodeling.bash', 'w')
f.write('#!/bin/bash\n') 
f.write('#$ -S /bin/bash\n') #use this shell to interpret script
#f.write('#$ -t %i-%i\n'%(start,end))
f.write('#$ -N L_%s\n' % pdir) #name of the job
f.write('#$ -cwd\n') # run script in current working directory
f.write('#$ -pe openmpi 21\n')
f.write('#$ -l mem_free=5.2G\n') #NGK use only
f.write('#$ -M etk243@nyumc.org\n') #send status updates
f.write('#$ -j y \n') #output SGE input and output into the same file
f.write('#$ -o log.out\n') # output file
f.write('module load openmpi/gcc/64/1.4.5\n')
f.write('\n') 
f.write('mpirun -np $NSLOTS loopmodel.mpi.linuxgccrelease @%s.flags -database $ROSETTA_DATABASE &>loop.log\n'% pdir)
f.close()

#Create a clustering and score file 
f=open(pdir + '/' + pdir + '.out/Combined/post_processing.bash', 'w')
f.write('#!/bin/bash\n')
f.write('#$ -S /bin/bash\n')
f.write('#$ -N P_%s\n'% pdir)
f.write('#$ -cwd\n')
f.write('#$ -j y\n')
f.write('#$ -o '+ pdir + '/'+ pdir +'.out/Combined/processing.out\n')
#f.write('#$ -o %s.out/Combined/processing.out\n'% pdir)
f.write('module load openmpi/gcc/64/1.4.5\n')
f.write('cd  %s/%s.out/Combined\n' % (pdir,pdir))
#f.write("ls * |egrep -v '*.silent|post*'|xargs rm -rf\n")
f.write('combine_silent.linuxgccrelease -in:file:silent ../*.silent.* -out:file:silent Combined.sile -out:file:silent_struct_type binary -silent_read_through_errors -database $ROSETTA_DATABASE -mute core.conformation.Conformation -mute core.scoring.etable -mute basic.io.database -mute core.io.silent.SilentFileData >&Combine.log\n')
f.write('Clusteranscore.py Combined.sile per.silent -1 0.5\n')
#f.write('rm -f ../*.silent.*\n')
f.write('mkdir -p clust\n')
f.write('cd clust\n')
f.write('cluster.linuxgccrelease -in:file:silent ../*.silent -in:file:silent_struct_type binary -in:file:fullatom -cluster:sort_groups_by_energy -out:pdb -database $ROSETTA_DATABASE -cluster:gdtmm -cluster:radius -1 -mute core.conformation.Conformation -mute protocols.cluster -mute basic.io.database -mute core.scoring.etable>cluster.log\n')
f.write('mpirun -np 10 score_jd2.mpi.linuxgccrelease -in:file:s *.pdb -out:pdb -database $ROSETTA_DATABASE >score.log\n')
f.write('for f in *.*_0001.pdb;do mv -f  $f "`echo $f |sed -r  s/.{9}$/.pdb/`";done\n')
f.write('\n')
f.write('pymol -cqd "run $SCRIPTS/analysis.py;analysis pdbs=*.pdb,native=%s,align_file=../../../*.ali"\n'% native)
f.write('tar -zcf %s_NGKF.tgz *.pdb --remove-files\n'% pdir)
f.write('tar -zcf ../combined_silent.tgz ../*.sile --remove-files\n')
f.write('rm -f *.pdb \n')
#f.write('while read line;do tar -zvxf *_NGKF.tgz $line;done<v1_bysecstruct.txt\n')
#f.write('while read line;do tar -zvxf *_NGKF.tgz $line;done<v2_bysecstruct.txt\n')
f.close()
call('chmod +x '+ pdir + '/' + pdir + '.out'+'/Combined/post_processing.bash',shell=True)

#Create the Launcher script to run sge
f = open('job.launch', 'a')
f.write('qsub %s/frags/frags.bash\n' %pdir)
f.close()
call('chmod +x job.launch', shell=True)

# Create analysis launch file:
f = open('analysis.launch', 'a')
f.write('qsub %s/%s.out/Combined/post_processing.bash\n' % (pdir,pdir))
f.close()
call('chmod +x analysis.launch', shell=True)

print ''
print 'NEXT STEPS: Run the following command to launch the parallel jobs'
print '        ./job.launch'
print ''
