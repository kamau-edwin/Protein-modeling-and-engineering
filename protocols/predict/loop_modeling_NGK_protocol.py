#!/usr/bin/python
########################################################################
# File generator for comparative modeling protocol.py
#
# Written by: Edwin Kamau and Brett Spurrier 
# Last modified: 5.28.2015 by Edwin Kamau
#
# Usage: python Template.py project template  target_sequence template_sequence
########################################################################
# import all the necessary functions required to make functional calls
import string, sys, re, os
from subprocess import call
from random import randint


if len(sys.argv)!=7:
    print 'Usage: python loonkic.py <project> <Template> <target_sequence> <template_sequence> <start> <nstruct>' 
    sys.exit(0)
#global variables 
pwd = os.getcwd()
pdir = sys.argv[1]  # name of the job obtained from item 2 of usage(see above)
PDB = sys.argv[2]
start = int(sys.argv[5])
#end = int(sys.argv[6])
nstruct=int(sys.argv[6])

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
#return_code = call("clustalo -i temp.fasta -o temp.aln --outfmt=clu --force", shell=True)  # maybe True works better?
return_code = call("clustalw2 temp.fasta",shell=True)

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
if not os.path.exists(pdir): # project directory
       os.makedirs(pdir)
if not os.path.exists(pdir + '/' + pdir + '.out'): # Rosetta output file
       os.makedirs(pdir + '/'+ pdir + '.out')
#if not os.path.exists(pdir + '/frags'): # Rosetta output file
#       os.makedirs(pdir + '/frags')
if not os.path.exists(pdir + '/' + pdir + '.out/Combined'): # Rosetta output file
       os.makedirs(pdir + '/'+ pdir +'.out/Combined')

f = open(pdir + '/'+  pdir + '.ali', 'w') # will create an alignment file used in modelling 
f.write('>%s\n%s\n>%s\n%s' % (pdir, full_target_sequence, PDB, full_template_sequence))
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
f.close()
#f = open(pdir + '/' + pdir + '.out/Combined/resfile','w')
#f.write('start\n')
#for index,res in enumerate(full_template_sequence):
#	if res !='-':
#		f.write('%s G\n' % (index+1))
#f.close()
#Cleanup
return_code = call("rm -f temp.*", shell=True)

# Create flag files
print 'Creating flag files...'
f = open(pdir + '/' + pdir + '.flags', 'w') # write below to pdir.flags in jobname folder
f.write('-s %s.pdb\n' % pdir)
f.write('-in:file:native ../starting_files/start_input/%s_MinPack.pdb\n' % PDB) 
f.write('-in:file:fullatom\n')
f.write('-loops:loop_file ' + pdir + '.loops\n')
f.write('-loops:remodel perturb_kic\n')
f.write('-loops:refine refine_kic\n')
f.write('-out:nstruct %i\n' % nstruct)
f.write('-out:file:silent '+ pdir + '.out/' + pdir + '.silent.%s \n' % start)
f.write('-out:file:fullatom\n')
f.write('-out:file:silent_struct_type binary\n')
f.write('-kic_bump_overlap_factor 0.36\n')
f.write('-legacy_kic false\n')
f.write('-kic_min_after_repack true\n')
f.write('-corrections:score:use_bicubic_interpolation false\n')
#NGK flags
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
f.write('-mute protocols.looprelax.FragmentPerturber\n') 
f.write('-mute core.fragments.ConstantLengthFragSet\n')
f.write('-mute protocols.loop_modeling.loggers.ProgressBar\n')
f.write('-run:seed_offset %s\n' %start)
f.close()
#Create a shell script bash file. Bash file called by qsub 
f = open(pdir + '/' + pdir + '.bash', 'w')
f.write('#!/bin/bash\n') 
f.write('#$ -S /bin/bash\n') #use this shell to interpret script
#f.write('#$ -t %i-%i\n'%(start,end))
f.write('#$ -N loop_%s\n' % start) #name of the job
f.write('#$ -cwd\n') # run script in current working directory
f.write('#$ -pe openmpi 50-51\n')
f.write('#$ -M etk243@nyumc.org\n') #send status updates
f.write('#$ -j y \n') #output SGE input and output into the same file
f.write('#$ -o %s/error.out\n' % pdir) # output file
f.write('#$ -l mem_free=6G\n')
f.write('module load openmpi/gcc/64/1.4.5\n')
f.write('\n') 
f.write('cd %s\n' % pdir)
f.write('thread_pdb_from_alignment.py --template=%s --target=%s --chain=G --align_format=fasta %s.ali $INPUT/%s_MinPack.pdb %s.pdb\n' %(PDB,pdir,pdir,PDB,pdir))
f.write('\n') 
f.write('mpirun -np $NSLOTS loopmodel.mpi.linuxgccrelease @%s.flags -database $ROSETTA_DATABASE &>loop.log\n'% pdir)
f.close()

#Create a clustering and score file 
f=open(pdir + '/' + pdir + '.out/Combined/post_processing.bash', 'w')
f.write('#!/bin/bash\n')
f.write('#$ -S /bin/bash\n')
f.write('#$ -N post_%s\n'% pdir)
f.write('#$ -cwd\n')
f.write('#$ -j y\n')
f.write('#$ -o '+ pdir + '/'+ pdir +'.out/Combined/processing.out\n')
#f.write('#$ -o %s.out/Combined/processing.out\n'% pdir)
f.write('module load openmpi/gcc/64/1.4.5\n')
f.write('cd  %s/%s.out/Combined\n' % (pdir,pdir))
f.write("ls * |egrep -v '*.silent|post*'|xargs rm -rf\n")
f.write('combine_silent.linuxgccrelease -in:file:silent ../*.silent.* -out:file:silent Combined.sile -out:file:silent_struct_type binary -silent_read_through_errors -database $ROSETTA_DATABASE -mute core.conformation.Conformation -mute core.scoring.etable -mute basic.io.database -mute core.io.silent.SilentFileData >&Combine.log\n')
f.write('Clusteranscore.py Combined.sile per.silent -1 1.0\n')
f.write('rm -f ../*.silent.*\n')
f.write('mkdir -p clust\n')
f.write('cd clust\n')
f.write('rm -rf *\n')
#f.write('cluster.linuxgccrelease -in:file:silent ../*.silent -in:file:fullatom --cluster:radius 1.0 -cluster:population_weight 0.0 -cluster:sort_groups_by_energy -no_optH -out:pdb -database $ROSETTA_DATABASE -mute core.conformation.Conformation -mute core.scoring.etable -mute basic.io.database >& cluster.log\n')
f.write('extract_pdbs.linuxgccrelease -in:file:silent ../*.silent -database $ROSETTA_DATABASE >extract.log')
f.write('mpirun -np 10 score_jd2.mpi.linuxgccrelease -in:file:s *.pdb -out:pdb -database $ROSETTA_DATABASE >score.log\n')
f.write('for f in *_*_0001.pdb;do mv -f  $f "`echo $f |sed -r  s/.{9}$/.pdb/`";done\n')
f.write('\n')
f.write('pymol -cqd "run $SCRIPTS/align_allfiles.py;align_allfiles `ls $INPUT/%s_MinPack.pdb,cycles=0,outfile=%s_loopmodel"\n'%(PDB,pdir))
#f.write('pymol -cqd "run $SCRIPTS/align_allfiles.py;align_allfiles `ls $WORK/%s.MinPack.pdb,cycles=0,outfile=%s;align_allfiles ../../fixbb/%s.pdb,cycles=0,res_sel=i. %s-%s,outfile=V2_loop"\n'%(PDB,pdir,beg_res,end_res))
f.write('tar -zcxf %s_NGK.tgz *.pdb --remove-files\n'% pdir)
f.write('tar -zcxf ../combined_silent.tgz ../*.sile --remove_files\n')
f.close()
call('chmod +x '+ pdir + '/' + pdir + '.out'+'/Combined/post_processing.bash',shell=True)

#Create the Launcher script to run sge
f = open('job.launch', 'w')
f.write('qsub %s/%s.bash\n' % (pdir,pdir))
f.close()
call('chmod +x job.launch', shell=True)

# Create analysis launch file:
f = open('analysis.launch', 'a')
f.write('qsub %s/%s.out/Combined/post_processing.bash\n' % (pdir,pdir))
f.close()
call('chmod +x analysis.launch', shell=True)

# Download fragment files from robetta
#f =open('frag_picker.launch','a')
#f.write('qsub '+ pdir + '/' + 'frags/frags.bash\n')
#f.close()
#call('chmod +x frag_picker.launch',shell=True)

#f = open(pdir + '/frags/frags.bash','w')
#f.write('#!/bin/bash\n')
#f.write('#$ -S /bin/bash\n')
#f.write('#$ -N frag_%s\n' % pdir)
#f.write('#$ -cwd\n')
#f.write('#$ -j y\n')
#f.write('#$ -o %s/frags/log.out\n' % pdir)
#f.write('cd ' + pdir + '/frags\n')
#f.write('#ls * |egrep -v "frags.bash"|xargs rm\n')
#f.write('#module load openmpi/gcc/64/1.4.5\n')
#f.write('cd ..\n')
#f.write('rm -f *.pdb \n')
#f.write("ls *.* |egrep -v '*.launch|%s'|xargs rm\n" % pdir)
#f.write('\n')
#f.write('#thread_pdb_from_alignment.py --template=%s --target=%s --chain=G --align_format=fasta %s.ali $INPUT/%s_MinPack.pdb %s.pdb\n' %(PDB,pdir,pdir,PDB,pdir))
#f.write('\n')
#f.write('./job.launch\n')
#f.close()
print ''
print 'NEXT STEPS: Run the following command to launch the parallel jobs'
print '        ./frag_picker.launch'
print ''
