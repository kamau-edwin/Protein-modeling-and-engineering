#!/usr/bin/python
########################################################################
# File generator for premut protocol.py
#
# Written by: Edwin Kamau 
# Last modified: 8.8.2015 by Edwin Kamau
#
########################################################################
# import all the necessary functions required to make functional calls
import sys, re, os
from subprocess import call

if len(sys.argv)!=4:
   print 'Usage: seqrank_design.py <project> design seq'
   sys.exit(0)

#global variables 
pdir = sys.argv[1]  # name of the job obtained from item 2 of usage(see above)
pdir=pdir.lower()
seqs=sys.argv[3]
seqs=seqs.split('_')
design=sys.argv[2]
# make directories 
if not os.path.exists(pdir+'/input/seqrank'): # project directory
       os.makedirs(pdir+'/input/seqrank')
if not os.path.exists(pdir+'/input/seqrank/br'): # project directory
       os.makedirs(pdir+'/input/seqrank/br')
if not os.path.exists(pdir+'/input/seqrank/fix'): # project directory
       os.makedirs(pdir+'/input/seqrank/fix')
#Create a shell script bash file. Bash file called by qsub 
if design =='fix':
	f = open(pdir+'/input/seqrank/fix/seqrank_design.bash', 'w')
	f.write('#!/bin/bash\n') 
	f.write('#$ -S /bin/bash\n')
	f.write('#$ -N des_%s\n' % pdir[:3])
	f.write('#$ -cwd\n')
	f.write('#$ -j y\n')
	f.write('#$ -o processing.out\n')
	f.write('\n')
	f.write('cd %s/input/seqrank/fix\n'% pdir)
	f.write('\n')
	#i=0
	#while i<len(seqs):	
	 #  f.write('seqtol_resfile.py `ls /ifs/home/etk243/HIV/Rosetta/input/PG9/input_design/*fix*%s*pdb` 6.0 %s_%s PIKAA H:152 H:153 H:154 H:155 H:156 H:157 H:159 152:%s 153:%s 154:%s 155:%s 156:%s 157:%s 159:%s\n'%(pdir,pdir,seqs[i],seqs[i][0],seqs[i][1],seqs[i][2],seqs[i][3],seqs[i][4],seqs[i][5],seqs[i][6]))
	  # i +=1
	f.write('\n')
	f.write('qsub fix.bash\n')
	f.close()
	f=open(pdir + '/input/seqrank/fix/fix.bash','w')
	f.write('#!/bin/bash\n')
	f.write('#$ -S /bin/bash\n')
	f.write('#$ -N des_%s\n' % pdir[:3])
	f.write('#$ -cwd\n')
	f.write('#$ -j y\n')
	f.write('#$ -o post_processing.out\n')
	f.write('\n')
	#f.write('for res in *.resfile\n')
	#f.write('do\n')
	#f.write("    prefix=`echo $res |sed 's/%s//'|sed 's/_seqtol.resfile//'`\n" % pdir)	
	#f.write('    fixbb.linuxgccrelease -database $ROSETTA_DATABASE -s $WORK/PG9/input_design/*fix*%s*pdb -in:file:fullatom -out:file:fullatom -nstruct 1 -resfile $res -suffix $prefix -min_pack -ex1 -ex2 -extrachi_cutoff 0 >seqrank_design.log\n' % pdir)
	#f.write('done\n')
	#f.write('\n')
	f.write('pymol -cqd "run $SCRIPTS/seqtol.py"\n')
	f.close()
	call('chmod +x '+ pdir + '/input/seqrank/fix/seqrank_design.bash', shell=True)
	call('chmod +x '+ pdir + '/input/seqrank/fix/fix.bash', shell=True)

elif design =='br':
	f = open(pdir+'/input/seqrank/br/seqrank_design.bash', 'w')
	f.write('#!/bin/bash\n') 
	f.write('#$ -S /bin/bash\n')
	f.write('#$ -N des_%s\n' % pdir[:3])
	f.write('#$ -cwd\n')
	f.write('#$ -j y\n')
	f.write('#$ -o post_processing.out\n')
	f.write('\n')	
	f.write('cd %s/input/seqrank/br\n'% pdir)
	f.write('\n')
#	i=0
#	while i<len(seqs):	
#	   f.write('seqtol_resfile.py `ls /ifs/home/etk243/HIV/Rosetta/input/PG9/input_design/*br*%s*pdb` 6.0 %s_%s PIKAA H:152 H:153 H:154 H:155 H:156 H:157 H:159 152:%s 153:%s 154:%s 155:%s 156:%s 157:%s 159:%s\n'%(pdir,pdir,seqs[i],seqs[i][0],seqs[i][1],seqs[i][2],seqs[i][3],seqs[i][4],seqs[i][5],seqs[i][6]))
#	   i +=1
	f.write('qsub br.bash\n')
	f.close()
	f=open(pdir + '/input/seqrank/br/br.bash','w')
	f.write('#!/bin/bash\n')
	f.write('#$ -S /bin/bash\n')
	f.write('#$ -N des_%s\n' % pdir[:3])
	f.write('#$ -cwd\n')
	f.write('#$ -j y\n')
	f.write('#$ -o processing.out\n')
	f.write('\n')
	f.write('for res in *.resfile\n')
	f.write('do\n')
	f.write("    prefix=`echo $res |sed 's/%s//'|sed 's/_seqtol.resfile//'`\n" % pdir)	
	f.write('    fixbb.linuxgccrelease -database $ROSETTA_DATABASE -s $WORK/PG9/input_design/*br*%s*pdb -in:file:fullatom -out:file:fullatom -nstruct 1 -resfile $res -suffix $prefix -min_pack -ex1 -ex2 -extrachi_cutoff 0 >seqrank_design.log\n' % pdir)
	f.write('done\n')
	f.write('\n')
	f.write('pymol -cqd "run $SCRIPTS/seqtol.py"\n')
	f.close()
	call('chmod +x '+ pdir + '/input/seqrank/br/seqrank_design.bash', shell=True)
	call('chmod +x '+ pdir + '/input/seqrank/br/br.bash', shell=True)

#Create the Launcher script to run sge
f=open('seqrank.launch', 'a')
if design=='fix':
    f.write('qsub %s/input/seqrank/fix/seqrank_design.bash\n' % pdir)
elif design=='br':
    f.write('qsub %s/input/seqrank/br/seqrank_design.bash\n'% pdir)
f.close()
call('chmod +x seqrank.launch', shell=True)

print 'generated startup files for project %s' % pdir
print 'NEXT STEPS: Run the following command to launch the parallel jobs'
print '   ./design.launch'
print ''
