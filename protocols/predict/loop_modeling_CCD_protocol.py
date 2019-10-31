#!/usr/bin/python
########################################################################
# File generator for comparative modeling protocol.py
#
# Written by: Edwin Kamau and Brett Spurrier 
# Last modified: 12.30.2014 by Edwin Kamau
#
# Usage: python Template.py project template  target_sequence template_sequence
########################################################################
# import all the necessary functions required to make functional calls
import string, sys, re, os
from subprocess import call
from random import randint


if len(sys.argv)!=8:
    print 'Usage: python homology.py <project> <Template> <target_sequence> <template_sequence> <start> <end> <nstruct>' 
    sys.exit(0)
#global variables 
pwd = os.getcwd()
pdir = sys.argv[1]  # name of the job obtained from item 2 of usage(see above)
PDB = sys.argv[2]
start = int(sys.argv[5])
end = int(sys.argv[6])
nstruct=int(sys.argv[7])

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
        print 'Aligning with clustalw2' 
        f = open('./temp.fasta', 'w') 
        f.write('>Target\n%s\n>Template1\n%s' % (target_sequence, template_sequence))
        f.close()

#Run the alignment
return_code = call("clustalo -i temp.fasta -o temp.aln --outfmt=clu", shell=True)  # maybe True works better?
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
if not os.path.exists(pdir): # project directory
       os.makedirs(pdir)
if not os.path.exists(pdir + '/' + pdir + '.out'): # Rosetta output file
       os.makedirs(pdir + '/'+ pdir + '.out')
if not os.path.exists(pdir + '/frags'): # Rosetta output file
       os.makedirs(pdir + '/frags')
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
	#f.write('LOOP %i %i 0 0 0\n' % (m.start()-1,m.end()))
for m in re.finditer('-+',str(full_template_sequence)):
	f.write('LOOP %i %i %i 0 1\n' % (m.start(),m.end()+1,m.end()))
	print 'LOOP %i %i %i 0 1' % (m.start(),m.end()+1,m.end())	
f.close()
#Cleanup
#return_code = call("rm -f temp.*", shell=True)

# Create flag files
print 'Creating flag files...'
f = open(pdir + '/' + pdir + '.flags', 'w') # write below to pdir.flags in jobname folder
f.write('-s %s.pdb\n' %pdir)
f.write('-in:file:fullatom\n')
f.write('#-in:file:native ../starting_files/start_input/%s_MinPack.pdb\n' % PDB)       
f.write('-loops:loop_file ' + pdir + '.loops\n')
f.write('-loops:frag_files frags/' + pdir + '.200.9mers frags/' + pdir + '.200.3mers none\n')
f.write('-loops:frag_sizes 9 3 1\n')
f.write('-loops:remodel quick_ccd\n')
f.write('-loops:refine refine_ccd\n')i
f.write('-loops:outer_cycles 5\n')
f.write('-out:nstruct %i\n' % nstruct)
f.write('-out:file:silent '+ pdir + '.out/' + pdir + '.silent.${SGE_TASK_ID} \n')
f.write('-out:file:fullatom\n')
f.write('-out:file:silent_struct_type binary\n')
f.write('-ex1\n')
f.write('-ex2\n')
f.write('-extrachi_cutoff 0\n')
f.write('-overwrite\n')
f.write('-mute core.io.database\n')
f.write('-mute protocols.looprelax.FragmentPerturber\n') 
f.write('-mute core.fragments.ConstantLengthFragSet\n')
f.write('-mute protocols.loop_modeling.loggers.ProgressBar\n')
f.write('-run:seed_offset ${SGE_TASK_ID}\n')
f.close()
#Create a shell script bash file. Bash file called by qsub 
f = open(pdir + '/' + pdir + '.bash', 'w')
f.write('#!/bin/bash\n') 
f.write('#$ -S /bin/bash\n') #use this shell to interpret script
f.write('#$ -t %i-%i\n'%(start,end))
f.write('#$ -N loop_%s\n' % pdir) #name of the job
f.write('#$ -cwd\n') # run script in current working directory
f.write('#$ -M etk243@nyumc.org\n') #send status updates
f.write('#$ -j y \n') #output input and output into the same file
f.write('#$ -l mem_free=6G\n') #output input and output into the same file
#f.write('#$ -o ' + os.getcwd() + '/' + pdir +'/error.out\n') # output file
f.write('#$ -o error.out\n') # output file
f.write('module load openmpi/gcc/64/1.4.5\n')
f.write('loopmodel.linuxgccrelease @%s.flags -database $ROSETTA_DATABASE &>loop.log\n'% pdir)
f.write('rm -f temp\n')
f.write('exit 0;\n')
f.close()

#Create a clustering and score file 
f=open(pdir + '/' + pdir + '.out/Combined/post_processing.bash', 'w')
f.write('#!/bin/bash\n')
f.write('#$ -S /bin/bash\n')
f.write('#$ -N post_%s\n'% pdir)
f.write('#$ -cwd\n')
f.write('#$ -j y\n')
f.write('#$ -o '+ pdir + '/'+ pdir +'.out/Combined/processing.out\n')
f.write('module load openmpi/gcc/64/1.4.5\n')
f.write('cd '+ pdir + '/' + pdir +'.out/Combined\n')
f.write('combine_silent.linuxgccrelease -in:file:silent ../*.silent.* -out:file:silent Combined.sile -out:file:silent_struct_type binary -silent_read_through_errors -database $ROSETTA_DATABASE -mute core.conformation.Conformation -mute core.scoring.etable -mute basic.io.database -mute core.io.silent.SilentFileData >&Combine.log\n')
f.write('Clusteranscore.py Combined.sile per.silent -1 0.5\n')
f.write('rm -f ../*.silent.*\n')
f.write('mkdir -p clust\n')
f.write('cd clust\n')
f.write('cluster.linuxgccrelease -in:file:silent ../*.silent -in:file:fullatom -cluster:gdtmm -cluster:radius 1.0 -cluster:population_weight 0 -cluster:sort_groups_by_energy -no_optH -out:pdb -database $ROSETTA_DATABASE -mute core.conformation.Conformation -mute core.scoring.etable -mute basic.io.database >& cluster.log\n')
f.write('score_jd2.linuxgccrelease -in:file:s *.pdb -no_scorefile -out:pdb -database $ROSETTA_DATABASE >score.log\n')
f.write('for f in *_0001.pdb;do mv -f  $f "`echo $f |sed s/_0001//`";done\n')
f.write('score_vs_rmsd_full.py -o %s -g %s -m -n ../*pdb  *pdb\n' % (pdir,pdir))
f.write('rm -f ../*pdb\n') 
#f.write('i=0\n')
#f.write('mkdir -p best20\n')
#f.write('while read line\n')
#f.write('do\n')
#f.write('    if [ "$i" -lt 20 ]; then\n')
#f.write('        mv $line best20\n')
#f.write('	((i++))\n')
#f.write('    fi\n')
#f.write('done < byscore\n')
f.write('exit 0;\n')
f.close()
call('chmod +x '+ pdir + '/' + pdir + '.out'+'/Combined/post_processing.bash',shell=True)

#Create the Launcher script to run sge
f = open(pdir + '/job.launch', 'w')
#f.write('qsub ' + pdir + '/' + pdir + '.bash\n')
f.write('qsub ' + pdir + '.bash\n')
f.close()
call('chmod +x '+ pdir + '/job.launch', shell=True)

# Create analysis launch file:
f = open( 'analysis.launch', 'a')
f.write('qsub '+ pdir + '/' + pdir + '.out/Combined/post_processing.bash\n')
f.close()
call("chmod +x analysis.launch", shell=True)

# Download fragment files from robetta
f =open('frag_picker.launch','a')
f.write('qsub '+ pdir + '/' + 'frags/frags.bash\n')
f.close()
call('chmod +x frag_picker.launch',shell=True)

f = open(pdir + '/frags/frags.bash','w')
f.write('#!/bin/bash\n')
f.write('#$ -S /bin/bash\n')
f.write('#$ -N frag_%s\n' % pdir)
f.write('#$ -cwd\n')
f.write('#$ -j y\n')
f.write('#$ -o %s/log.out\n' % pdir)
f.write('cd ' + pdir + '/frags\n')
f.write('module load openmpi/gcc/64/1.4.5\n')
f.write('make_fragments.pl -id %s -frag_sizes 9,3 ../%s.fasta\n' %(pdir,pdir))
f.write('ls *.* |egrep -v  "*mers|*ss2|*.bash" |xargs rm\n') 
f.write('cd ..\n')
f.write('rm -f *.pdb\n')
f.write('ls *.* |egrep -v "job.launch|%s" |xargs rm\n' % pdir)
f.write('\n')
f.write('fixbb.linuxgccrelease -database $ROSETTA_DATABASE -s $INPUT/../%s.pdb  -in:file:fullatom -out:file:fullatom -out:pdb -out:path:pdb $INPUT -out:suffix _MinPack -min_pack -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -flip_HNQ -no_optH false -packing:repack_only >fix.log\n'% PDB)
f.write('mv $INPUT/%s_MinPack_* $INPUT/%s_MinPack.pdb\n'% (PDB,PDB))) 
f.write('thread_pdb_from_alignment.py --template=%s --target=%s --chain=G --align_format=fasta %s.ali $INPUT/%s_MinPack.pdb %s.pdb\n' %(PDB,pdir,pdir,PDB,pdir))
f.write('\n')
f.write('./job.launch\n')
f.write('exit;\n')
print ''
print 'NEXT STEPS: Run the following command to launch the parallel jobs'
print '        ./frag_picker.launch'
print ''
