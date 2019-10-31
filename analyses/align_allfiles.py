#! /usr/bin/env python
# original Written by Jules Jacobsen (jacobsen@ebi.ac.uk). Feel free to do whatever you like with this code.
# modified by Edwin kamau on 6.28.15
# Does alignment and grabs score terms from clustered pdbs after design or loop modeling
# outputs a csv file with score terms and rmsd that can then be used for score vs rmsd plots
# extensively modified by Robert L. Campbell (rlc1@queensu.ca)
#from __future__ import print_function
from pymol import cmd
import glob
import pandas as pd
import re
import string
from os import popen,remove
import matplotlib as mpl
mpl.use('AGG')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tarfile 

def scorevsrmsd(input_path,name,out='_score_vs_rmsd_plot'):
    df=pd.read_csv(input_path)
    col=['score','bb','all','ca']
    ext=['REU','rmsd']
    x=df[df[col[0]]==min(df[col[0]])].index.tolist()[0]
    n=len(df)    
    title=name.split('_')
    xlower=min(df[col[1]])-0.009
    xupper=max(df[col[1]])+0.01
    ylower=min(df[col[0]])-1.0
    yupper=max(df[col[0]])+1.0
    fig,ax1=plt.subplots()
    fig2,ax2=plt.subplots()
    ax1.scatter(df[col[1]],df[col[0]],lw=0,color='black')
    ax2.scatter(df[col[2]],df[col[0]],lw=0,color='black')
    ax1.scatter(df.loc[x,col[1]],min(df[col[0]]),color='red')
    ax2.scatter(df.loc[x,col[2]],min(df[col[0]]),color='red')
    ax1.set_title('Backbone',fontsize=18,fontweight='bold')
    ax2.set_title('All Atoms',fontsize=18,fontweight='bold')
    ax1.set_xlabel(ext[1],fontsize=16,fontstyle='oblique',fontweight='bold')
    ax2.set_xlabel(ext[1],fontsize=16,fontstyle='oblique',fontweight='bold')
    ax1.set_ylabel(col[0],fontsize=16,fontstyle='oblique',fontweight='bold')
    ax2.set_ylabel(col[0],fontsize=16,fontstyle='oblique',fontweight='bold')
    ax1.text(0.98, 0.03,'n=%s'%n ,
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax1.transAxes,
        color='grey', fontsize=12)
    ax2.text(0.98, 0.03,'n=%s'%n ,
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax1.transAxes,
        color='grey', fontsize=12)
    #ax1.axis([xlower,xupper,ylower,yupper])
    #ax2.axis([xlower,xupper,ylower,yupper])
    fig.savefig(name+out+'_bb',dpi=300,bbox_inches='tight')
    fig2.savefig(name+out+'_all',dpi=300,bbox_inches='tight')
# original align_allfiles source of all pymol arguments passed in bash
def align_allfiles(target='',files='*.pdb',res_sel='i. 1-',outfile='',name='',cutoff=2.0, cycles=5,method='align'):
  """
  Aligns all models in a list of files to one target using either align or super method. 
  Outputs: Rosetta energy terms,c-alpha, backbone, and all atom rmsd for each model 
  sorted by c-alpha rmsd

  usage:
    align_allfiles target=target name,files=name.pdb,res_sel='i. 1-,outfile='file name'
    cutoff=2,cycles=5 method=align
        where target specifies the model id you want to align all others against,
        and  cutoff and cycles are options passed to the align command.
        You can specify the files to load and align using a wildcard. You should specify same
	number of residues for the target and other models. if outfile is not given output is
        stddount

    	By default all residues in the structure are aligned,cutoff is 2,number of cycles is 5,        and align method is used.

    Example:
      pymol -cdq "run/$SCRIPTS/align_allfiles;align_allfiles <target name>,[none default 
      arguments eg res_sel= i. 100-120,...]"

  """
  cutoff = int(cutoff)
  cycles = int(cycles)
  cmd.load(target) # change 1
  #print outfile, cutoff,cycles,res_sel 
  #Define rmsd and residue selection for mobile and target
  all='all & %s' % res_sel
  ca='n. ca & %s' % res_sel
  bb='n. n+c+ca+o and %s' % res_sel
#obtain files from folder and remove native
  #tar=tarfile.open('*.tgz','r:tgz')
  file_list=glob.glob(files)
  #file_list.remove(target)  
  #file_list.sort()
  for file in file_list:
    if file==target:
	file_list.remove(file)
    else:
    	file_list.sort()
  #print file_list
 # print len(file_list)
  extension = re.compile( '(^.*[\/]|\.(pdb|ent|brk))' )
  object_list = []
  target=extension.sub('',target) # change 3
  # Define rmsd storage variables
  rmsd_list=[]
  for i in range(len(file_list)):
    obj_name1 = extension.sub('',file_list[i])
    object_list.append(extension.sub('',file_list[i]))
  #  energy=popen('grep "total_energy"' file_list[i])
   # print energy 
    cmd.load(file_list[i],obj_name1)
    if method == 'align':
        rms_ca = cmd.align('%s & %s'%(object_list[i],ca),'%s & %s'%(target,ca),cutoff=cutoff,cycles=cycles)
        rms_bb = cmd.align('%s & %s'%(object_list[i],bb),'%s & %s'%(target,bb),cutoff=cutoff,cycles=cycles)
        rms_all = cmd.align('%s & %s'%(object_list[i],all),'%s & %s'%(target,all),cutoff=cutoff,cycles=cycles)
        #print '%s,%6.2f,%6.2f,%6.2f' %(object_list[i],rms_bb[0],rms_ca[0],rms_all[0])
    elif method == 'super':
        rms_ca = cmd.super('%s & %s'%(object_list[i],ca),'%s & %s'%(target,ca),cutoff=cutoff,cycles=cycles)
        rms_bb = cmd.super('%s & %s'%(object_list[i],bb),'%s & %s'%(target,bb),cutoff=cutoff,cycles=cycles)
        rms_all = cmd.super('%s & %s'%(object_list[i],all),'%s & %s'%(target,all),cutoff=cutoff,cycles=cycles)
        #print '%s,%6.2f,%6.2f,%6.2f' %(object_list[i],rms_ca[0],rms_bb[0],rms_all[0])
    elif method == 'cealign':
        rmsdict = cmd.cealign('%s & %s' % (target,target_selection),'%s & %s' % (object_list[i],mobile_selection))
        rms = [rmsdict['RMSD'],rmsdict['alignment_length'],1,0,0]
    else:
        print "only 'align', 'super' and 'cealign' are accepted as methods"
        sys.exit(-1)
#rms = cmd.align('%s & %s'%(object_list[i],mobile_selection),'%s & %s'%(target,target_selection),cutoff=cutoff,cycles=cycles,object=objectname)
#else:
#rms = cmd.align('%s & %s'%(object_list[i],mobile_selection),'%s & %s'%(target,target_selection),cutoff=cutoff,cycles=cycles)

    #rmsd[object_list[i]] = (rms_ca[0],rms_ca[1],rms_bb[0],rms_bb[1],rms_all[0],rms_all[1])
    
    rmsd_list.append((object_list[i],rms_ca[0],rms_bb[0],rms_all[0]))
    cmd.delete(obj_name1)
  #print rmsd
  rmsd_list.sort(lambda x,y: cmp(x[2],y[2]))#compare ca rms you can assign other rms indexes to sort  i.e index 2 for bb and index 3 for all atom
# loop over dictionary and print out matrix of final rms values
  outfp=outfile
  table=False 
  out=open(outfp +'_score_vs_rmsd.csv','w')
  out.write('model,score,ca,bb,all,fa_atr,fa_rep,fa_sol,fa_elec\n')
  for r in rmsd_list:
      for decoy in file_list:
        model=open(decoy,'r')
        for line in model:
          line_split=line.split()
          if line_split[0]=='#BEGIN_POSE_ENERGIES_TABLE':
             table=True
	     continue 
          if table and line_split[0]=='pose':
            score=float(line_split[-1])
            fa_atr=float(line_split[1])
            fa_rep=float(line_split[2])
            fa_elec=float(line_split[5])
            fa_sol=float(line_split[3])
            if r[0]==decoy.replace('.pdb',''):
             out.write('%s,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f\n' %(decoy,score,r[1],r[2],r[3],fa_atr,fa_rep,fa_sol,fa_elec))
             print '%s,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f\n' %(decoy,score,r[1],r[2],r[3],fa_atr,fa_rep,fa_sol,fa_elec)
  out.close()
#generate ensemble figure
  protocol=outfp.split('_')[1]
  cmd.bg_color(color='white')
  cmd.hide("all")
  cmd.show("cartoon")
  cmd.set("antialias",1)
  cmd.set("ray_trace_mode",0)
  cmd.set("depth_cue",0)
  cmd.set("ray_trace_fog",0)
  cmd.set("stick_radius", 0.1)
  cmd.set("cartoon_side_chain_helper",1)
  cmd.set("cartoon_flat_sheets",1)
#  cmd.set("cartoon_transparency",0.8)
  cmd.set("ray_shadow",0)
  if protocol =='loopmodel':  
  	cmd.rotate("x",70)
  	cmd.zoom("resi 143-163")
  	cmd.rotate("z",20)
  	cmd.rotate("y",45)
  	cmd.move("y",8)
  	cmd.move("x",3)
  	cmd.save(outfp+'ensemble.pse') 
  	cmd.delete("all")
  elif protocol=='fixed' or protocol =='backrub':
  	util.cbab("chain h")
  	util.cbac("chain l")
  	util.cbam("chain g")
  	cmd.turn("y",-140)
  	cmd.turn("x",-10)
  	cmd.turn("z",10)
  	cmd.move("y",-5)
  	cmd.move("z",30)
  	cmd.png(outfp + '_ensemble.png',2400,2400,dpi=300,ray=1)
  	cmd.save(outfp+'ensemble.pse') 
  	cmd.delete("all")
#generate score vs rmsd plot 
  scorevsrmsd(outfp+'_score_vs_rmsd.csv',name=outfp)
cmd.extend('align_allfiles',align_allfiles)
