#! /local/apps/python/2.7.3 
import operator
import re
import glob
from pymol import cmd
from subprocess import call
import time
import sys
import os
import matplotlib as mpl
mpl.use('AGG')
import matplotlib.pyplot as plt
sys.path.append('/ifs/home/etk243/HIV/scripts/rosetta/tools/protein_tools/scripts')
from amino_acids import longer_names
if sys.version_info < (2, 7):
    raise Exception("You must use python2.7 to run this script")
def loop_model_analysis(pdbs=None,align_file=None,native=None,name='n. CA',cycles=5,cutoff=2.0):
  """
  Calculates backbone and all atom rmsd for the modeled loop region,obtains model score and modeled loop structural characteristics
  SCRIPT methodology:
  inputs pdbs in the current directory,reads each pdb file filtered by the cluster it belongs to and extracts the score.
  Creates a sorted by score cluster dictionary.The lowest scoring decoy from each cluster is then used as the target for backbone and all atom rmsd calculations based on the modeled loop region only
  Generates an output file containing decoys with helical, strand or mixed conformation, with their score, rmsd, and modeled loop secondary structure 
  Calculates statistics for each cluster and overall loop modeling  
  Defaults:
	pymol alignment cutoff and cycles. 
  Needed
    At least one selection 
    Alignment file path
    Decoys
  Example
    pymol -cqd "run path to analysis file;analysis pdbs= *.pdb,cutoff,cycles,..."
    """    
  help='Usage:pymol -cqd "run analysis.py;analysis pdbs= path_to_pdbs or tar zipped pdb file,align_file=path_to_alignment_file, name=n. ca, cutoff=2.0, cycles=5"\nNote: Default alignment cutoff and cycles can be changed if necessary'
  start=float(time.time())
  id= os.getcwd().split('/')[7]
  #id= os.getcwd().split('/')[5]
  ID=id.upper()
  print 'Native %s' % native
  print ''
  if native=='None':
     native=None  
  if not pdbs:
    print 'Missing pdb files or zipped pdb file.See Usage'
    print ''
    print help
    cmd.quit()
  if pdbs.split('.')[1]=='pdb' and os.path.exists('c.0.0.pdb'):
     files=glob.glob(pdbs)
  elif pdbs.split('.')[1]=='pdb' and not os.path.exists('c.0.0.pdb'):
     print 'tarred file exist.Using it instead of given pdbs argument'
     pdbs='*.tgz'
     print 'untarring file...'
     call('tar -zxf %s' %pdbs,shell=True)
     files=glob.glob('*.pdb')
  elif pdbs.split('.')[1]=='tgz' and not os.path.exists('c.0.0.pdb'):
     print 'untarring file ...'
     call('tar -zxf %s' %pdbs,shell=True)
     files=glob.glob('*.pdb')
  else:
     files=glob.glob('*.pdb')
  if not align_file:
    print 'Target-template alignment file missing.See Usage'
    print ''
    print help
    cmd.quit()	
  ali_file=glob.glob(align_file) # alignment file with gaps
  extension = re.compile( '(^.*[\/]|\.(pdb|ent|brk))' )
    # dictionary to hold decoy and score
  clust0_dic={};clust1_dic={};clust2_dic={};clust3_dic={};clust4_dic={}
    # sort by energy variables
  cluster0="";cluster1="";Cluster2="";cluster3="";cluster4="";energy=[]
    # alignment parameters
  #bb='n. n+ca+c+o'
  bb='n. N+CA+C+O'
  all='all'
    # create model sel1 and sel2 files from alignment file
  beg=[] # Modeled loop region beginning index
  end=[] # Modeled loop region end index
    # open the file and extract beginning and end of gaps and append
  with open(ali_file[0],'r') as f:
    for line in f:
      if '-' in line:
        for m in re.finditer('-+',str(line)):
          beg.append(m.start()+1)
          end.append(m.end())
  if len(beg)==2 and len(end)==2: # check if more than one loop was modeled
      sel1='i. %s-%s' %(beg[0],end[0])
      sel2='i. %s-%s' %(beg[1],end[1])
  elif len(beg)==1 and len(end)==1: # check if one loop was modeled and whether V1 or V2
    if beg[0]<32:
      sel1='i. %s-%s' %(beg[0],end[0])
      sel2=None
    elif beg[0]>32:
      sel2='i. %s-%s' %(beg[0],end[0])
      sel1=None
    # obtain modeled residues from modeling selection 
  print 'Selection 1: %s' %sel1
  print 'Selection 2: %s\n' %sel2
  model=files[0]
  cmd.load(model)
  model=extension.sub('',model)
  v1_res={'res':[]}
  v2_res={'res':[]}
  if sel1: # pick residues in the modeled loop region
    cmd.iterate('%s & %s & %s'%(model,sel1,name),'res.append((%s[resn]))' % longer_names,space=v1_res)
    #cmd.iterate('%s & %s & %s'%(model,sel1,name),'res.append((%s[resn]))' % one_letter,space=v1_res)
  if sel2:
    cmd.iterate('%s & %s & %s'%(model,sel2,name),'res.append((%s[resn]))' % longer_names,space=v2_res)
    #cmd.iterate('%s & %s & %s'%(model,sel2,name),'res.append((%s[resn]))' % one_letter,space=v2_res)
  v1_loop=''.join(v1_res['res'])
  v2_loop=''.join(v2_res['res'])
  v1v2=v1_loop+v2_loop
# extract score for each decoy create a dictionary and sort it
  cluster_names=[]
  for pdb in files:
    clus=int(pdb.split('.')[1])
    cluster_names.append(clus)
    with open(pdb,'r') as f:
      for line in f:
        if line.startswith('pose'):
          score=float('%0.2f'%float(line.split()[-1]))
          energy.append(int(score))
          if pdb[2]=='0':
            clust0_dic[pdb]= score
            cluster0=sorted(clust0_dic.items(),key=operator.itemgetter(1))
          elif pdb[2]=='1':
            clust1_dic[pdb]=score
            cluster1=sorted(clust1_dic.items(),key=operator.itemgetter(1))
          elif pdb[2]=='2':
              clust2_dic[pdb]= score
              Cluster2=sorted(clust2_dic.items(),key=operator.itemgetter(1))
          elif pdb[2]=='3':
            clust3_dic[pdb]=score
            cluster3=sorted(clust3_dic.items(),key=operator.itemgetter(1))
          elif pdb[2]=='4':
            clust4_dic[pdb]= score
            cluster4=sorted(clust4_dic.items(),key=operator.itemgetter(1))          
	  else:
            print 'No file with decoys found or clusters exceed set limit'
            print help
  main=[]
  if cluster0:
   for rec in cluster0:
     main.append(rec)
  if cluster1:
   for rec in cluster1:
     main.append(rec)
  if Cluster2:
     for rec in Cluster2:
       main.append(rec)
  if cluster3:
     for rec in cluster3:
       main.append(rec)
  if cluster4:
     for rec in cluster4:
       main.append(rec)
  main.sort(lambda x,y:cmp(x[1],y[1]))
  low_score=min(energy)-1
  high_score=max(energy)+1
  cluster_names=sorted(list(set(cluster_names)))
  print 'Lowest scoring model %s score %s ' %(main[0][0],main[0][1])
  print 'Finished sorting and storing scores and decoys for all clusters\n'
  print 'Running alignment...'
#COLOR PALETTE FOR V1 PLOTS
  cb_palette = [(65, 68, 81), (31, 119, 180),(153, 86, 136), (123, 102, 210), (188, 189, 34), (23, 190, 207)]
  for i in range(len(cb_palette)):
    r,g,b=cb_palette[i]
    cb_palette[i]=(r/255.,g/255.,b/255.)
#COLOR PALETTE FOR V2 PLOTS
  cb1_palette = [(89, 89, 89), (0, 107, 164),(44, 160, 44), (148, 103, 189), (143, 135, 130), (23, 190, 207)]
  for i in range(len(cb1_palette)):
    r,g,b=cb1_palette[i]
    cb1_palette[i]=(r/255.,g/255.,b/255.)
  ss_palette = [(204,102,119),(68,170,153),(136,204,238)]
  for i in range(len(ss_palette)):
    r,g,b=ss_palette[i]
    ss_palette[i]=(r/255.,g/255.,b/255.)
  loaded=0
  def score_vs_rmsd():
    """
    Runs alignment based on lowest scoring model from each cluster and all cluster and outputs score vs rmsd figure
    Requires all structures and at least one selection 
    """
    fig,ax=plt.subplots()
    fig2,a=plt.subplots()
    fig3,bx=plt.subplots()
    g,bx1=plt.subplots()
    g2,bx2=plt.subplots()
    fig4,b=plt.subplots()
    for model,score in main:
      if loaded==len(main):
        if native:
         target=native
         cmd.load(target)
         target=extension.sub('',target.split('/')[-1])
         obj_name=extension.sub('',model)
        else:
         target=extension.sub('',main[0][0])
         obj_name=extension.sub('',model)
      elif loaded !=len(main):
        if native:
         target=native
         cmd.load(target)
         target=extension.sub('',target.split('/')[-1])
         obj_name=extension.sub('',model)
         cmd.load(model,obj_name)
        else:
         target=main[0][0]
         cmd.load(target)
         target=extension.sub('',target)
         obj_name=extension.sub('',model)
         cmd.load(model,obj_name)
      if sel1:
        rms_bb=cmd.align('%s & %s & %s'%(obj_name,bb,sel1),'%s & %s & %s'%(target,bb,sel1),cutoff=cutoff,cycles=cycles)
        rms_all = cmd.align('%s & %s & %s'%(obj_name,all,sel1),'%s & %s & %s'%(target,all,sel1),cutoff=cutoff,cycles=cycles)
        if cluster0 and obj_name[2]=='0':
          ax.scatter('%0.2f'%(rms_bb[0]),score,marker="o",color=cb_palette[1])
          a.scatter('%0.2f'%(rms_all[0]),score,marker="o",color=cb_palette[1])
          if obj_name==extension.sub('',cluster0[0][0]):
             bx1.scatter('%0.2f'%(rms_bb[0]),score,marker="o",color=cb_palette[1],label='0')
        elif cluster1 and obj_name[2]=='1':
          ax.scatter('%0.2f'%(rms_bb[0]),score,marker="h",color=cb_palette[2])
          a.scatter('%0.2f'%(rms_all[0]),score,marker="h",color=cb_palette[2])
          if obj_name==extension.sub('',cluster1[0][0]):
              bx1.scatter('%0.2f'%(rms_bb[0]),score,marker="h",color=cb_palette[2],label='1')
        elif Cluster2 and obj_name[2]=='2':
          ax.scatter('%0.2f'%(rms_bb[0]),score,marker="d",color=cb_palette[3])
          a.scatter('%0.2f'%(rms_all[0]),score,marker="d",color=cb_palette[3])
          if obj_name==extension.sub('',Cluster2[0][0]):
            bx1.scatter('%0.2f'%(rms_bb[0]),score,marker="d",color=cb_palette[3],label='2')
        elif cluster3 and obj_name[2]=='3':
          ax.scatter('%0.2f'%(rms_bb[0]),score,marker="p",color=cb_palette[4])
          a.scatter('%0.2f'%(rms_all[0]),score,marker="p",color=cb_palette[4])
          if obj_name==extension.sub('',cluster3[0][0]):
           bx1.scatter('%0.2f'%(rms_bb[0]),score,marker="p",color=cb_palette[4],label='3')
        elif cluster4 and obj_name[2]=='4':
          ax.scatter('%0.2f'%(rms_bb[0]),score,marker="x",color=cb_palette[5])
          a.scatter('%0.2f'%(rms_all[0]),score,marker="x",color=cb_palette[5])
          if obj_name==extension.sub('',cluster4[0][0]):
            bx1.scatter('%0.2f'%(rms_bb[0]),score,marker="x",color=cb_palette[5],label='4')
        if obj_name==extension.sub('',main[0][0]):
          ax.scatter('%0.2f'%(rms_bb[0]),score,marker="*",color='red',s=25)
          a.scatter('%0.2f'%(rms_all[0]),score,marker="*",color='red',s=25)
        #ax.set_title('ALL CLUSTERS',fontweight='medium',fontsize=8,loc='left')
        #a.set_title('ALL CLUSTERS',fontweight='medium',fontsize=8,loc='left')
        ax.tick_params(axis='both',top='off',right='off',labelsize=8)
        a.tick_params(axis='both',top='off',right='off',labelsize=8)
        #ax.set_xlim([-0.2,5])
        #a.set_xlim([-0.2,5])
        ax.set_ylim([low_score,high_score])
        a.set_ylim([low_score,high_score])
        ax.text(0.99, 0.01,'N=%s'%len(main),
           verticalalignment='bottom', horizontalalignment='right',
          transform=ax.transAxes,
          color='grey', fontsize=8)
        a.text(0.99, 0.01,'N=%s'%len(main),
           verticalalignment='bottom', horizontalalignment='right',
          transform=a.transAxes,
          color='grey', fontsize=8)
      if sel2:
        rms_bb=cmd.align('%s & %s &%s'%(obj_name,bb,sel2),'%s & %s & %s'%(target,bb,sel2),cutoff=cutoff,cycles=cycles)
        rms_all = cmd.align('%s & %s & %s'%(obj_name,all,sel2),'%s & %s & %s'%(target,all,sel2),cutoff=cutoff,cycles=cycles)
        if cluster0 and obj_name[2]=='0':
          bx.scatter('%0.2f'%(rms_bb[0]),score,marker="o",color=cb1_palette[1])
          b.scatter('%0.2f'%(rms_all[0]),score,marker="o",color=cb1_palette[1])
          if obj_name==extension.sub('',cluster0[0][0]):
             bx2.scatter('%0.2f'%(rms_bb[0]),score,marker="o",color=cb1_palette[1],label='0')
        elif cluster1 and obj_name[2]=='1':
          bx.scatter('%0.2f'%(rms_bb[0]),score,marker="h",color=cb1_palette[2])
          b.scatter('%0.2f'%(rms_all[0]),score,marker="h",color=cb1_palette[2])
          if obj_name==extension.sub('',cluster1[0][0]):
           bx2.scatter('%0.2f'%(rms_bb[0]),score,marker="h",color=cb1_palette[2],label='1')
        elif Cluster2 and obj_name[2]=='2':
          bx.scatter('%0.2f'%(rms_bb[0]),score,marker="d",color=cb1_palette[3])
          b.scatter('%0.2f'%(rms_all[0]),score,marker="d",color=cb1_palette[3])
          if obj_name==extension.sub('',Cluster2[0][0]):
            bx2.scatter('%0.2f'%(rms_bb[0]),score,marker="d",color=cb1_palette[3],label='2')
        elif cluster3 and obj_name[2]=='3':
          bx.scatter('%0.2f'%(rms_bb[0]),score,marker="p",color=cb1_palette[4])
          b.scatter('%0.2f'%(rms_all[0]),score,marker="p",color=cb1_palette[4])
          if obj_name==extension.sub('',cluster3[0][0]):
            bx2.scatter('%0.2f'%(rms_bb[0]),score,marker="p",color=cb1_palette[4],label='3')
        elif cluster4 and obj_name[2]=='4':
          bx.scatter('%0.2f'%(rms_bb[0]),score,marker="x",color=cb1_palette[5])
          b.scatter('%0.2f'%(rms_all[0]),score,marker="x",color=cb1_palette[5])
          if obj_name==extension.sub('',cluster4[0][0]):
            bx2.scatter('%0.2f'%(rms_bb[0]),score,marker="x",color=cb1_palette[5],label='4')
        if obj_name==extension.sub('',main[0][0]):
          bx.scatter('%0.2f'%(rms_bb[0]),score,marker="*",color='red',s=30)
          b.scatter('%0.2f'%(rms_all[0]),score,marker="*",color='red',s=30)
        #bx.set_title('Backbone',fontweight='medium',fontsize=8,loc='left')
        #b.set_title('All Atoms',fontweight='medium',fontsize=8,loc='left')
        bx.tick_params(axis='both',top='off',right='off',labelsize=8)
        b.tick_params(axis='both',top='off',right='off',labelsize=8)
    #      bx.set_xlim([-0.2,5])
    #      b.set_xlim([-0.2,5])
        bx.set_ylim([low_score,high_score])
        b.set_ylim([low_score,high_score])
        bx.text(0.99, 0.01,'N=%s'%len(main),
             verticalalignment='bottom', horizontalalignment='right',
            transform=bx.transAxes,
            color='grey', fontsize=8)
        b.text(0.99, 0.01,'N=%s'%len(main),
             verticalalignment='bottom', horizontalalignment='right',
            transform=b.transAxes,
            color='grey', fontsize=8)
    #bel2=list(set(label2))
  # MAKE plots
    if sel1:
      handle1,label1=bx1.get_legend_handles_labels()
      handle1,label1=zip(*(sorted(zip(handle1,label1),key=operator.itemgetter(1))))
      fig.text(0.5,0.050,'Backbone RMSD',ha='center',va='center',fontweight='semibold')
      fig.legend(handle1,label1,loc=1,scatterpoints=1,fontsize=12,title='cluster',bbox_to_anchor=(0.86,0.88))
      #fig.legend(loc='best',scatterpoints=1,fontsize=8)
      fig.text(0.06, 0.5, 'SCORE', ha='center', va='center', rotation='vertical',fontweight='semibold')
      fig.text(0.5, 0.920, ID+':V1', ha='center', va='center',fontweight='bold')
      fig.savefig(id+'_v1_backbone',dpi=300,bbox_inches='tight')
      fig2.text(0.5,0.050,'All Atom RMSD',ha='center',va='center',fontweight='semibold')
      fig2.legend(handle1,label1,loc=1,scatterpoints=1,fontsize=12,title='cluster',bbox_to_anchor=(0.86,0.88))
      fig2.text(0.08, 0.5, 'SCORE', ha='center', va='center', rotation='vertical',fontweight='semibold')
      fig2.text(0.5, 0.920, ID+':V1', ha='center', va='center',fontweight='bold')
      fig2.savefig(id+'_v1_all_atom',dpi=300,bbox_inches='tight')
    if sel2:
      handle2,label2=bx2.get_legend_handles_labels()
      handle2,label2=zip(*(sorted(zip(handle2,label2),key=operator.itemgetter(1))))
      fig3.text(0.5,0.05,'Backbone RMSD',ha='center',va='center',fontweight='semibold')
      fig3.legend(handle2,label2,loc=1,scatterpoints=1,fontsize=12,title='cluster',bbox_to_anchor=(0.86,0.88))
      fig3.text(0.06, 0.5, 'SCORE', ha='center', va='center', rotation='vertical',fontweight='semibold')
      fig3.text(0.5, 0.920, ID+':V2', ha='center', va='center',fontweight='bold')
      fig3.savefig(id+'_v2_backbone',dpi=300,bbox_inches='tight')
      fig4.text(0.5,0.05,'All Atom RMSD',ha='center',va='center',fontweight='semibold')
      fig4.legend(handle2,label2,loc=1,scatterpoints=1,fontsize=12,title='cluster',bbox_to_anchor=(0.86,0.88))
      #fig4.legend(loc='best',scatterpoints=1,fontsize=8)
      fig4.text(0.06, 0.5, 'SCORE', ha='center', va='center', rotation='vertical',fontweight='semibold')
      fig4.text(0.5, 0.920, ID+':V2', ha='center', va='center',fontweight='bold')
      fig4.savefig(id+'_v2_all_atom',dpi=300,bbox_inches='tight')
    if not sel1 and not sel2:
      print 'No selection passed to the function'
      cmd.quit()
# RUN SECONDARY STRUCTURE ASSIGNMENT
  cluster_names=sorted(list(set(cluster_names)))
  if os.path.exists('cluster0.csv'):
    print 'Found secondary structure assigned file'
    print 'Running score vs rmsd plots'
    start_align=time.time()
    score_vs_rmsd()
    end_align=(time.time()-start_align)
    end=time.time()-start    
    print 'Timings'
    print '  Total time %0.2fs' % float(end)
    print '  Align and plot %0.2fs' % float(end_align)
    cmd.quit()
  else:
    print 'No secondary structure assigned file found'
    print 'loading pdbs ...'
    for model,score in main:
	obj_name=extension.sub('',model)
	cmd.load(model,obj_name)
        loaded +=1
  print ' %s models loaded' % loaded
  print 'Assigning Secondary structure...'
  target=native
  cmd.load(target)
  target=extension.sub('',target.split('/')[-1])
  fig,ax=plt.subplots()
  fig2,bx=plt.subplots()
  #g,bx1=plt.subplots()
  g2,bx2=plt.subplots()
  total_loop_v1 =0;total_loop_v2 =0
  total_helix_v1 =0;total_helix_v2 =0
  total_strand_v1 =0;total_strand_v2 =0
  total_hs_v1 =0;total_hs_v2 =0
  if sel1 or sel2:
    resid=beg[0]+118
  if sel1 and sel2:
    resid=beg[0]+120
    resid2=beg[1]+116
  temp=sys.stdout 
  processed_all=[]
  if sel1:
    print 'Processing based on selection 1'
  if sel2:
    print 'Processing based on selection 2'
  if sel1 and sel2:
    print 'Processing based on selection 1 and 2'
  loop_v1_c0=0;loop_v1_c1=0;loop_v1_c2=0;loop_v1_c3=0;loop_v1_c4=0
  helix_v1_c0=0;helix_v1_c1=0;helix_v1_c2=0;helix_v1_c3=0;helix_v1_c4=0
  strand_v1_c0=0;strand_v1_c1=0;strand_v1_c2=0;strand_v1_c3=0;strand_v1_c4=0
  hs_v1_c0=0;hs_v1_c1=0;hs_v1_c2=0;hs_v1_c3=0;hs_v1_c4=0
  res_helix_v1_c0={};res_helix_v1_c1={};res_helix_v1_c2={};res_helix_v1_c3={};res_helix_v1_c4={}
  res_strand_v1_c0={};res_strand_v1_c1={};res_strand_v1_c2={};res_strand_v1_c3={};res_strand_v1_c4={}
  res_hs_v1_c0={};res_hs_v1_c1={};res_hs_v1_c2={};res_hs_v1_c3={};res_hs_v1_c4={}
  loop_v2_c0=0;loop_v2_c1=0;loop_v2_c2=0;loop_v2_c3=0;loop_v2_c4=0
  helix_v2_c0=0;helix_v2_c1=0;helix_v2_c2=0;helix_v2_c3=0;helix_v2_c4=0
  strand_v2_c0=0;strand_v2_c1=0;strand_v2_c2=0;strand_v2_c3=0;strand_v2_c4=0  
  hs_v2_c0=0;hs_v2_c1=0;hs_v2_c2=0;hs_v2_c3=0;hs_v2_c4=0
  res_helix_v2_c0={};res_helix_v2_c1={};res_helix_v2_c2={};res_helix_v2_c3={};res_helix_v2_c4={}
  res_strand_v2_c0={};res_strand_v2_c1={};res_strand_v2_c2={};res_strand_v2_c3={};res_strand_v2_c4={}
  res_hs_v2_c0={};res_hs_v2_c1={};res_hs_v2_c2={};res_hs_v2_c3={};res_hs_v2_c4={}
  v2_5H_c0=[];v2_2H_c0=[];v2_3H_c0=[];v2_4H_c0=[]
  v2_5S_c0=[];v2_2S_c0=[];v2_3S_c0=[];v2_4S_c0=[]
  v2_5HS_c0=[];v2_2HS_c0=[];v2_3HS_c0=[];v2_4HS_c0=[]
  v1_5H_c0=[];v1_2H_c0=[];v1_3H_c0=[];v1_4H_c0=[]
  v1_5S_c0=[];v1_2S_c0=[];v1_3S_c0=[];v1_4S_c0=[]
  v1_5HS_c0=[];v1_2HS_c0=[];v1_3HS_c0=[];v1_4HS_c0=[]
  if cluster0:
   log=open('cluster0.csv','w')
   print>>log,'Decoy,loop,sec_struct,score,bb,all,ss'
   start_c0=time.time()
   for model,score in cluster0:
     obj_name=extension.sub('',model)
     cluster={'sec':[]}
     cluster2={'sec':[]}
     if sel1:
      cmd.iterate('%s & %s & %s'%(obj_name,sel1,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        processed_all.append(model)
        helix_v1_c0 +=1
        total_helix_v1 +=1
        print>>log,model,',','v1',',','H',',',score,',',cluster
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)==2:
          v1_2H_c0.append(model)
        if len(ss)==3:
          v1_3H_c0.append(model)
        if len(ss)==4:
          v1_4H_c0.append(model)
        if len(ss)>=5:
          v1_5H_c0.append(model)
        for i in ss:
          res_helix_v1_c0[i+resid] =v1_loop[i]
      if 'S' in cluster and 'H' not in cluster:
        processed_all.append(model)
        strand_v1_c0 +=1
        total_strand_v1 +=1
        print>>log,model,',','v1',',','S',',',score,',',cluster
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)==2:
          v1_2S_c0.append(model)
        if len(ss)==3:
          v1_3S_c0.append(model)
        if len(ss)==4:
          v1_4S_c0.append(model)
        if len(ss)>=5:
          v1_5S_c0.append(model)
        for i in ss:
          res_strand_v1_c0[i+resid]=v1_loop[i]
      if 'H' in cluster and 'S' in cluster:
        processed_all.append(model)
        hs_v1_c0 +=1
        total_hs_v1 +=1
        print>>log,model,',','v1',',','HS',',',score,',',cluster
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)==2:
            v1_2HS_c0.append(model)
        if len(ss)==3:
            v1_3HS_c0.append(model)
        if len(ss)==4:
            v1_4HS_c0.append(model)
        if len(ss)>=5:
            v1_5HS_c0.append(model)
        for i in ss:
          res_hs_v1_c0[i+resid]=v1_loop[i]
      if 'H' not in cluster and 'S' not in cluster:
        processed_all.append((model,score))
        loop_v1_c0 +=1
        total_loop_v1 +=1
     if sel2:
      cmd.iterate('%s & %s & %s'%(obj_name,sel2,name),'sec.append(ss)',space=cluster2)
      cluster2=''.join(cluster2['sec'])
      rms_bb=cmd.align('%s & %s & %s'%(obj_name,bb,sel2),'%s & %s & %s'%(target,bb,sel2),cutoff=cutoff,cycles=cycles)
      if 'H' in cluster2 and 'S' not in cluster2:
        processed_all.append(model)
        helix_v2_c0 +=1
        total_helix_v2 +=1
        print>>log,model,',','v2',',','H',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\alpha$',s=70,color=ss_palette[1],label='$Helix$')
        ss=[m.start() for m in re.finditer('H',str(cluster2))]
        if len(ss)==2:
          v2_2H_c0.append(model)
        if len(ss)==3:
          v2_3H_c0.append(model)
        if len(ss)==4:
          v2_4H_c0.append(model)
        if len(ss)>=5:
          v2_5H_c0.append(model)
        for i in ss:
          res_helix_v2_c0[i+resid] =v2_loop[i]
      if 'S' in cluster2 and 'H' not in cluster2:
        processed_all.append(model)
        strand_v2_c0 +=1
        total_strand_v2 +=1
        print>>log,model,',','v2',',','S',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\beta$',s=70,color=ss_palette[0],label='$Strand$')
        ss=[m.start() for m in re.finditer('S',str(cluster2))]
        if len(ss)>=5:
            v2_5S_c0.append(model)
        if len(ss)==2:
            v2_2S_c0.append(model)
        if len(ss)==3:
            v2_3S_c0.append(model)
        if len(ss)==4:
            v2_4S_c0.append(model)
        for i in ss:
          res_strand_v2_c0[i+resid]=v2_loop[i]
      if 'H' in cluster2 and 'S' in cluster2:
        processed_all.append(model)
        hs_v2_c0 +=1
        total_hs_v2 +=1
        print>>log,model,',','v2',',','HS',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\alpha\beta$',s=70,color=ss_palette[2],label='$Helix-Strand$')
        ss=[m.start() for m in re.finditer('H|S',str(cluster2))]
        if len(ss)>=5:
            v2_5HS_c0.append(model)
        if len(ss)==2:
            v2_2HS_c0.append(model)
        if len(ss)==3:
            v2_3HS_c0.append(model)
        if len(ss)==4:
            v2_4HS_c0.append(model)
        for i in ss:
          res_hs_v2_c0[i+resid]=v2_loop[i]
      if 'H' not in cluster2 and 'S' not in cluster2:
        print>>log,model,',','v2',',','L',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker='$L$',s=70,color='b',label='$Loop$')
        processed_all.append(model)
        loop_v2_c0 +=1
        total_loop_v2 +=1
   log.close()
   time_c0=time.time()-start_c0
# CLUSTER 1 Sel1 only
  v1_5H_c1=[];v1_2H_c1=[];v1_3H_c1=[];v1_4H_c1=[]
  v1_5S_c1=[];v1_2S_c1=[];v1_3S_c1=[];v1_4S_c1=[]
  v1_5HS_c1=[];v1_2HS_c1=[];v1_3HS_c1=[];v1_4HS_c1=[]
  v2_5H_c1=[];v2_2H_c1=[];v2_3H_c1=[];v2_4H_c1=[]
  v2_5S_c1=[];v2_2S_c1=[];v2_3S_c1=[];v2_4S_c1=[]
  v2_5HS_c1=[];v2_2HS_c1=[];v2_3HS_c1=[];v2_4HS_c1=[]
  if cluster1:
   start_c1=time.time()
   log1=open('cluster1.csv','w')
   print>>log1,'Decoy,loop,sec_struct,score,bb,all,ss'
   for model,score in cluster1:
    obj_name=extension.sub('',model)
    cluster={'sec':[]}
    cluster2={'sec':[]}
    if sel1:
      cmd.iterate('%s & %s & %s'%(obj_name,sel1,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        processed_all.append(model)
        helix_v1_c1 +=1
        total_helix_v1 +=1
        print>>log1,model,',','v1',',','H',',',score,',',cluster
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)>=5:
          v1_5H_c1.append(model)
        if len(ss)==2:
          v1_2H_c1.append(model)
        if len(ss)==3:
          v1_3H_c1.append(model)
        if len(ss)==4:
          v1_4H_c1.append(model)
        for i in ss:
          res_helix_v1_c1[i +resid]=v1_loop[i]
      if 'S' in cluster and 'H' not in cluster:
        processed_all.append(model)
        strand_v1_c1 +=1
        total_strand_v1 +=1
        print>>log1,model,',','v1',',','S',',',score,',',cluster
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)>=5:
          v1_5S_c1.append(model)
        if len(ss)==2:
          v1_2S_c1.append(model)
        if len(ss)==3:
          v1_3S_c1.append(model)
        if len(ss)==4:
          v1_4S_c1.append(model)
        for i in ss:
          res_strand_v1_c1[i+resid]=v1_loop[i]
      if 'H' in cluster and 'S' in cluster:
        processed_all.append(model)
        hs_v1_c1 +=1
        total_hs_v1 +=1
        print>>log1,model,',','v1',',','HS',',',score,',',cluster
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)>=5:
          v1_5HS_c1.append(model)
        if len(ss)==2:
          v1_2HS_c1.append(model)
        if len(ss)==3:
          v1_3HS_c1.append(model)
        if len(ss)==4:
          v1_4HS_c1.append(model)
        for i in ss:
          res_hs_v1_c1[i+resid]=v1_loop[i]
      if 'H' not in cluster and 'S' not in cluster:
        processed_all.append(model)
        loop_v1_c1 +=1
        total_loop_v1 +=1
    if sel2:
      cmd.iterate('%s & %s & %s'%(obj_name,sel2,name),'sec.append(ss)',space=cluster2)
      cluster2=''.join(cluster2['sec'])
      rms_bb=cmd.align('%s & %s &%s'%(obj_name,bb,sel2),'%s & %s & %s'%(target,bb,sel2),cutoff=cutoff,cycles=cycles)
      if 'H' in cluster2 and 'S' not in cluster2:
        processed_all.append(model)
        helix_v2_c1 +=1
        total_helix_v2 +=1
        print>>log1, model,',','v2',',','H',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\alpha$',s=70,color=ss_palette[1],label='$Helix$')
        ss=[m.start() for m in re.finditer('H',str(cluster2))]
        if len(ss)>=5:
          v2_5H_c1.append(model)
        if len(ss)==2:
          v2_2H_c1.append(model)
        if len(ss)==3:
          v2_3H_c1.append(model)
        if len(ss)==4:
          v2_4H_c1.append(model)
        for i in ss:
          res_helix_v2_c1[i +resid]=v2_loop[i]
      if 'S' in cluster2 and 'H' not in cluster2:
        print>>log1,model,',','v2',',','S',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\beta$',s=70,color=ss_palette[0],label='$Strand$')
        processed_all.append(model)
        strand_v2_c1 +=1
        total_strand_v2 +=1
        ss=[m.start() for m in re.finditer('S',str(cluster2))]
        if len(ss)>=1:
          v2_5S_c1.append(model)
        if len(ss)==2:
          v2_2S_c1.append(model)
        if len(ss)==3:
          v2_3S_c1.append(model)
        if len(ss)==4:
          v2_4S_c1.append(model)
        for i in ss:
          res_strand_v2_c1[i+resid]=v2_loop[i]
      if 'H' in cluster2 and 'S' in cluster2:
        print>>log1,model,',','v2',',','HS',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\alpha\beta$',color=ss_palette[2],label='$Helix-Strand$')
        processed_all.append(model)
        hs_v2_c1 +=1
        total_hs_v2 +=1
        ss=[m.start() for m in re.finditer('H|S',str(cluster2))]
        if len(ss)>=5:
          v2_5HS_c1.append(model)
        if len(ss)==2:
          v2_2HS_c1.append(model)
        if len(ss)==3:
          v2_3HS_c1.append(model)
        if len(ss)==4:
          v2_4HS_c1.append(model)
        for i in ss:
          res_hs_v2_c1[i+resid]=v2_loop[i]
      if 'H' not in cluster2 and 'S' not in cluster2:
        print>>log1,model,',','v2',',','L',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker="$L$",color='b',label='$Loop$')
        processed_all.append(model)
        loop_v2_c1 +=1
        total_loop_v2 +=1
   log1.close()
   time_c1=time.time()-start_c1
# CLUSTER 2 Sel1 only
  v1_5H_c2=[];v1_2H_c2=[];v1_3H_c2=[];v1_4H_c2=[]
  v1_5S_c2=[];v1_2S_c2=[];v1_3S_c2=[];v1_4S_c2=[]
  v1_5HS_c2=[];v1_2HS_c2=[];v1_3HS_c2=[];v1_4HS_c2=[]
  v2_5H_c2=[];v2_2H_c2=[];v2_3H_c2=[];v2_4H_c2=[]
  v2_5S_c2=[];v2_2S_c2=[];v2_3S_c2=[];v2_4S_c2=[]
  v2_5HS_c2=[];v2_2HS_c2=[];v2_3HS_c2=[];v2_4HS_c2=[]
  if Cluster2:
   start_c2=time.time()
   log2=open('cluster2.csv','w')
   print>>log2,'Decoy,loop,sec_struct,score,bb,all,ss'
   for model,score in Cluster2:
    obj_name=extension.sub('',model)
    cluster={'sec':[]}
    cluster2={'sec':[]}
    if sel1:
      cmd.iterate('%s & %s & %s'%(obj_name,sel1,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        processed_all.append(model)
        helix_v1_c2 +=1
        total_helix_v1 +=1
        print>>log2,model,',','v1',',','H',',',score,',',cluster
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)>=1:
          v1_5H_c2.append(model)
        if len(ss)==2:
          v1_2H_c2.append(model)
        if len(ss)==3:
          v1_3H_c2.append(model)
        if len(ss)==4:
          v1_4H_c2.append(model)
        for i in ss:
          res_helix_v1_c2[i+resid]=v1_loop[i]
      if 'S' in cluster and 'H' not in cluster:
        processed_all.append(model)
        strand_v1_c2 +=1
        total_strand_v1 +=1
        print>>log2,model,',','v1',',','S',',',score,',',cluster
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)>=5:
          v1_5S_c2.append(model)
        if len(ss)==2:
          v1_2S_c2.append(model)
        if len(ss)==3:
          v1_3S_c2.append(model)
        if len(ss)==4: 
          v1_4S_c2.append(model)     
        for i in ss:
          res_strand_v1_c2[i+resid]=v1_loop[i]
      if 'H' in cluster and 'S' in cluster:
        processed_all.append(model)
        hs_v1_c2 +=1
        total_hs_v1 +=1
        print>>log2,model,',','v1',',','HS',',',score,',',cluster
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)>=5:
          v1_5HS_c2.append(model)
        if len(ss)==2:
          v1_2HS_c2.append(model)
        if len(ss)==3:
          v1_3HS_c2.append(model)
        if len(ss)==4:
          v1_4HS_c2.append(model)
        for i in ss:
          res_hs_v1_c2[i+resid]=v1_loop[i]
      if 'H' not in cluster and 'S' not in cluster:
        processed_all.append(model)
        loop_v1_c2 +=1
        total_loop_v1 +=1
    if sel2:
      cmd.iterate('%s & %s & %s'%(obj_name,sel2,name),'sec.append(ss)',space=cluster2)
      cluster2=''.join(cluster2['sec'])
      rms_bb=cmd.align('%s & %s &%s'%(obj_name,bb,sel2),'%s & %s & %s'%(target,bb,sel2),cutoff=cutoff,cycles=cycles)
      if 'H' in cluster2 and 'S' not in cluster2:
        processed_all.append(model)
        helix_v2_c2 +=1
        total_helix_v2 +=1
        print>>log2,model,',','v2',',','H',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\alpha$',s=70,color=ss_palette[1],label='$Helix$')
        ss=[m.start() for m in re.finditer('H',str(cluster2))]
        if len(ss)>=5:
          v2_5H_c2.append(model)
        if len(ss)==2:
          v2_2H_c2.append(model)
        if len(ss)==3:
          v2_3H_c2.append(model)
        if len(ss)==4:
          v2_4H_c2.append(model)
        for i in ss:
          res_helix_v2_c2[i+resid]=v2_loop[i]
      if 'S' in cluster2 and 'H' not in cluster2:
        processed_all.append(model)
        strand_v2_c2 +=1
        total_strand_v2 +=1
        print>>log2,model,',','v2',',','S',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\beta$',s=70,color=ss_palette[0],label='$Strand$')
        ss=[m.start() for m in re.finditer('S',str(cluster2))]
        if len(ss)>=5:
          v2_5S_c2.append(model)
        if len(ss)==2:
          v2_2S_c2.append(model)
        if len(ss)==3:
          v2_3S_c2.append(model)
        if len(ss)==4: 
          v2_4S_c2.append(model)     
        for i in ss:
          res_strand_v2_c2[i+resid]=v2_loop[i]
      if 'H' in cluster2 and 'S' in cluster2:
        processed_all.append(model)
        hs_v2_c2 +=1
        total_hs_v2 +=1
        print>>log2,model,',','v2',',','HS',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\alpha\beta$',color=ss_palette[2],label='$Helix-Strand$')
        ss=[m.start() for m in re.finditer('H|S',str(cluster2))]
        if len(ss)>=5:
          v2_5HS_c2.append(model)
        if len(ss)==2:
          v2_2HS_c2.append(model)
        if len(ss)==3:
          v2_3HS_c2.append(model)
        if len(ss)==4:
          v2_4HS_c2.append(model)
        for i in ss:
          res_hs_v2_c2[i+resid]=v2_loop[i]
      if 'H' not in cluster2 and 'S' not in cluster2:
        processed_all.append(model)
        loop_v2_c2 +=1
        total_loop_v2 +=1
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker="$L$",color='b',label='$Loop$')
        print>>log2,model,',','v2',',','L',',',score,',','%0.2f'%rms_bb[0],',',cluster2
   time_c2=time.time()-start_c2
   log2.close()
# CLUSTER 3 Sel1 only
  v1_5H_c3=[];v1_2H_c3=[];v1_3H_c3=[];v1_4H_c3=[]
  v1_5S_c3=[];v1_2S_c3=[];v1_3S_c3=[];v1_4S_c3=[]
  v1_5HS_c3=[];v1_2HS_c3=[];v1_3HS_c3=[];v1_4HS_c3=[]
  v2_5H_c3=[];v2_2H_c3=[];v2_3H_c3=[];v2_4H_c3=[]
  v2_5S_c3=[];v2_2S_c3=[];v2_3S_c3=[];v2_4S_c3=[]
  v2_5HS_c3=[];v2_2HS_c3=[];v2_3HS_c3=[];v2_4HS_c3=[]
  if cluster3:
   start_c3=time.time()
   log3=open('cluster3.csv','w')
   print>>log3,'Decoy,loop,sec_struct,score,bb,all,ss'
   for model, score in cluster3:
    obj_name=extension.sub('',model)
    cluster={'sec':[]}
    cluster2={'sec':[]}
    if sel1:
      cmd.iterate('%s & %s & %s'%(obj_name,sel1,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        processed_all.append(model)
        helix_v1_c3 +=1
        total_helix_v1 +=1
        print>>log3,model,',','v1',',','H',',',score,',',cluster
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)>=5:
          v1_5H_c3.append(model)
        if len(ss)==2:
          v1_2H_c3.append(model)
        if len(ss)==3:
          v1_3H_c3.append(model)
        if len(ss)==4:
          v1_4H_c3.append(model)
        for i in ss:
          res_helix_v1_c3[i+resid]=v1_loop[i]
      if 'S' in cluster and 'H' not in cluster:
        processed_all.append(model)
        strand_v1_c3 +=1
        total_strand_v1 +=1
        print>>log3,model,',','v1',',','S',',',score,',',cluster
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)>=5:
          v1_5S_c3.append(model)
        if len(ss)==2:
          v1_2S_c3.append(model)
        if len(ss)==3:
          v1_3S_c3.append(model)
        if len(ss)==4:
          v1_4S_c3.append(model)
        for i in ss:
          res_strand_v1_c3[i+resid]=v1_loop[i]
      if 'H' in cluster and 'S' in cluster:
        processed_all.append(model)
        hs_v1_c3 +=1
        total_hs_v1 +=1
        print>>log3,model,',','v1',',','HS',',',score,',',cluster
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)>=5:
          v1_5HS_c3.append(model)
        if len(ss)==2:
          v1_2HS_c3.append(model)
        if len(ss)==3:
          v1_3HS_c3.append(model)
        if len(ss)==4:
          v1_4HS_c3.append(model)
        for i in ss:
          res_hs_v1_c3[i+resid]=v1_loop[i]
      if 'H' not in cluster and 'S' not in cluster:
        processed_all.append((model,score))
        loop_v1_c3 +=1
        total_loop_v1 +=1
    if sel2:
      cmd.iterate('%s & %s & %s'%(obj_name,sel2,name),'sec.append(ss)',space=cluster2)
      cluster2=''.join(cluster2['sec'])
      rms_bb=cmd.align('%s & %s &%s'%(obj_name,bb,sel2),'%s & %s & %s'%(target,bb,sel2),cutoff=cutoff,cycles=cycles)
      if 'H' in cluster2 and 'S' not in cluster2:
        processed_all.append(model)
        helix_v2_c3 +=1
        total_helix_v2 +=1
        print>>log3,model,',','v2',',','H',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\alpha$',color=ss_palette[1],label='$Helix$')
        ss=[m.start() for m in re.finditer('H',str(cluster2))]
        if len(ss)>=5:
          v2_5H_c3.append(model)
        if len(ss)==2:
          v2_2H_c3.append(model)
        if len(ss)==3:
          v2_3H_c3.append(model)
        if len(ss)==4:
          v2_4H_c3.append(model)
        for i in ss:
          res_helix_v2_c3[i+resid]=v2_loop[i]
      if 'S' in cluster2 and 'H' not in cluster2:
        processed_all.append(model)
        strand_v2_c3 +=1
        total_strand_v2 +=1
        print>>log3,model,',','v2',',','S',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\beta$',color=ss_palette[0],label='$Strand$')
        ss=[m.start() for m in re.finditer('S',str(cluster2))]
        if len(ss)>=5:
          v2_5S_c3.append(model)
        if len(ss)==2:
          v2_2S_c3.append(model)
        if len(ss)==3:
          v2_3S_c3.append(model)
        if len(ss)==4:
          v2_4S_c3.append(model)
        for i in ss:
          res_strand_v2_c3[i+resid]=v2_loop[i]
      if 'H' in cluster2 and 'S' in cluster2:
        processed_all.append(model)
        hs_v2_c3 +=1
        total_hs_v2 +=1
        print>>log3,model,',','v2',',','HS',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\alpha\beta$',color=ss_palette[2],label='$Helix-Strand$')
        ss=[m.start() for m in re.finditer('H|S',str(cluster2))]
        if len(ss)>=5:
          v2_5HS_c3.append(model)
        if len(ss)==2:
          v2_2HS_c3.append(model)
        if len(ss)==3:
          v2_3HS_c3.append(model)
        if len(ss)==4:
          v2_4HS_c3.append(model)
        for i in ss:
          res_hs_v2_c3[i+resid]=v2_loop[i]
      if 'H' not in cluster2 and 'S' not in cluster2:
        processed_all.append(model)
        loop_v2_c3 +=1
        total_loop_v2 +=1
        print>>log3,model,',','v2',',','L',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,s=70,marker="$L$",color='b',label='$Loop$')
   time_c3=time.time()-start_c3
   log3.close()
# CLUSTER  4  Sel1 only
  v1_5H_c4=[];v1_2H_c4=[];v1_3H_c4=[];v1_4H_c4=[]
  v1_5S_c4=[];v1_2S_c4=[];v1_3S_c4=[];v1_4S_c4=[]
  v1_5HS_c4=[];v1_2HS_c4=[];v1_3HS_c4=[];v1_4HS_c4=[]
  v2_5H_c4=[];v2_2H_c4=[];v2_3H_c4=[];v2_4H_c4=[]
  v2_5S_c4=[];v2_2S_c4=[];v2_3S_c4=[];v2_4S_c4=[]
  v2_5HS_c4=[];v2_2HS_c4=[];v2_3HS_c4=[];v2_4HS_c4=[]
  if cluster4:
   start_c4=time.time()
   log4=open('cluster4.csv','w')
   print>>log4,'Decoy,loop,sec_struct,score,bb,all,ss'
   for model,score in cluster4:
    obj_name=extension.sub('',model)
    cluster={'sec':[]}
    cluster2={'sec':[]}
    if sel1:
      cmd.iterate('%s & %s & %s'%(obj_name,sel1,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        processed_all.append(model)
        helix_v1_c4 +=1
        total_helix_v1 +=1
        print>>log4,model,',','v1',',','H',',',score,',',cluster
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)>=5:
          v1_5H_c4.append(model)
        if len(ss)==2:
          v1_2H_c4.append(model)
        if len(ss)==3:
          v1_3H_c4.append(model)
        if len(ss)==4:
          v1_4H_c4.append(model)
        for i in ss:
          res_helix_v1_c4[i+resid]=v1_loop[i]
      if 'S' in cluster and 'H' not in cluster:
        processed_all.append(model)
        strand_v1_c4 +=1
        total_strand_v1 +=1
        print>>log4,model,',','v1',',','S',',',score,',',cluster
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)>=5:
          v1_5S_c4.append(model)
        if len(ss)==2:
          v1_2S_c4.append(model)
        if len(ss)==3:
          v1_3S_c4.append(model)
        if len(ss)==4:
          v1_4S_c4.append(model)
        for i in ss:
          res_strand_v1_c4[i+resid]=v1_loop[i]
      if 'H' in cluster and 'S' in cluster:
        processed_all.append(model)
        hs_v1_c4 +=1
        total_hs_v1 +=1
        print>>log4,model,',','v1',',','HS',',',score,',',cluster
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)>=5:
          v1_5HS_c4.append(model)
        if len(ss)==2:
          v1_2HS_c4.append(model)
        if len(ss)==3:
          v1_3HS_c4.append(model)
        if len(ss)==4:
          v1_4HS_c4.append(model)
        for i in ss:
          res_hs_v1_c4[i+resid]=v1_loop[i]
      if 'H' not in cluster and 'S' not in cluster:
        processed_all.append(model)
        loop_v1_c4 +=1
        total_loop_v1 +=1
    if sel2:
      cmd.iterate('%s & %s & %s'%(obj_name,sel2,name),'sec.append(ss)',space=cluster2)
      cluster2=''.join(cluster2['sec'])
      rms_bb=cmd.align('%s & %s &%s'%(obj_name,bb,sel2),'%s & %s & %s'%(target,bb,sel2),cutoff=cutoff,cycles=cycles)
      if 'H' in cluster2 and 'S' not in cluster2:
        processed_all.append(model)
        helix_v2_c4 +=1
        total_helix_v2 +=1
        print>>log4,model,',','v2',',','H',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker="x",color='g',label='4:H')
        bx.scatter('%0.2f'%(rms_bb[0]),score,s=70,marker=r'$\alpha$',color=ss_palette[1],label='$Helix$')
        ss=[m.start() for m in re.finditer('H',str(cluster2))]
        if len(ss)>=5:
          v2_5H_c4.append(model)
        if len(ss)==2:
          v2_2H_c4.append(model)
        if len(ss)==3:
          v2_3H_c4.append(model)
        if len(ss)==4:
          v2_4H_c4.append(model)
        for i in ss:
          res_helix_v2_c4[i+resid]=v2_loop[i]
      if 'S' in cluster2 and 'H' not in cluster2:
        processed_all.append(model)
        strand_v2_c4 +=1
        total_strand_v2 +=1
        print>>log4,model,',','v2',',','S',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\beta$',s=70,color=ss_palette[0],label=r'$Strand$')
        ss=[m.start() for m in re.finditer('S',str(cluster2))]
        if len(ss)>=5:
          v2_5S_c4.append(model)
        if len(ss)==2:
          v2_2S_c4.append(model)
        if len(ss)==3:
          v2_3S_c4.append(model)
        if len(ss)==4:
          v2_4S_c4.append(model)
        for i in ss:
          res_strand_v2_c4[i+resid]=v2_loop[i]
      if 'H' in cluster2 and 'S' in cluster2:
        processed_all.append(model)
        hs_v2_c4 +=1
        total_hs_v2 +=1
        print>>log4,model,',','v2',',','HS',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker=r'$\alpha\beta$',s=70,color=ss_palette[2],label='$Helix-Strand$')
        ss=[m.start() for m in re.finditer('H|S',str(cluster2))]
        if len(ss)>=5:
          v2_5HS_c4.append(model)
        if len(ss)==2:
          v2_2HS_c4.append(model)
        if len(ss)==3:
          v2_3HS_c4.append(model)
        if len(ss)==4:
          v2_4HS_c4.append(model)
        for i in ss:
          res_hs_v2_c4[i+resid]=v2_loop[i]
      if 'H' not in cluster2 and 'S' not in cluster2:
        processed_all.append(model)
        loop_v2_c4 +=1
        total_loop_v2 +=1
        print>>log4,model,',','v2',',','L',',',score,',','%0.2f'%rms_bb[0],',',cluster2
        bx.scatter('%0.2f'%(rms_bb[0]),score,marker="$L$",s=70,color='b',label='$Loop$')
   time_c4=time.time()-start_c4
   log4.close()
   bx.set_ylim([low_score,high_score])
   bx.text(0.99, 0.01,'N=%s'%len(main),
        verticalalignment='bottom', horizontalalignment='right',
        transform=bx.transAxes,
        color='grey', fontsize=8)
  if sel2:
   handle1,label1=bx.get_legend_handles_labels()
   newleg=dict()
   for h,l in zip(handle1,label1):
     if l not in newleg.keys():
       newleg[l]=h
   handles=[]
   labels=[]
   for l in newleg.keys():
     handles.append(newleg[l])
     labels.append(l)
   print handles,labels
   bx.set_ylim([low_score,high_score])
   bx.text(0.99, 0.01,'N=%s'%len(main),
      verticalalignment='bottom', horizontalalignment='right',
      transform=bx.transAxes,
      color='grey', fontsize=8)
   fig2.text(0.5,0.040,'Backbone RMSD',ha='center',va='center',fontweight='semibold')
   fig2.legend(handles,labels,loc=8,ncol=len(labels),scatterpoints=1,fontsize=12,title='Sec.Struc',bbox_to_anchor=(0.5,0.084),frameon=False)
   fig2.text(0.04, 0.5, 'SCORE', ha='center', va='center', rotation='vertical',fontweight='semibold')
   #fig2.text(0.5, 0.920, ID+':V2', ha='center', va='center',fontweight='bold')
   fig2.savefig(id+'_v2_backbone',dpi=300,bbox_inches='tight')
#GATHERING modeling statistics 
  log=open(id+'_modeling_stat.csv','w')
  if sel1:
    print>>log,'Cluster,Loop,Total,Loop(Per),Helix(Per),Strand(Per),HelixRes,Strand_Res,HelixStrandRes'
    total_v1=loop_v1_c0+helix_v1_c0+strand_v1_c0+hs_v1_c0
    if total_v1>0:
     print>>log,'%s,V1,%s,%s(%0.1f),%s(%0.1f),%s(%0.1f),%s,%s,%s'%(cluster_names[0],total_v1,loop_v1_c0,float(loop_v1_c0*100.0/total_v1),helix_v1_c0,float(helix_v1_c0*100.0/total_v1),strand_v1_c0,float(strand_v1_c0*100.0/total_v1),res_helix_v1_c0,res_strand_v1_c0,res_hs_v1_c0)
  #cluster 1 stats
    total_v1=loop_v1_c1+helix_v1_c1+strand_v1_c1+hs_v1_c1
    if total_v1>0:
      print>>log,'%s,V1,%s,%s(%0.1f),%s(%0.1f),%s(%0.1f),%s,%s,%s'%(cluster_names[1],total_v1,loop_v1_c1,float(loop_v1_c1*100.0/total_v1),helix_v1_c1,float(helix_v1_c1*100.0/total_v1),strand_v1_c1,float(strand_v1_c1*100.0/total_v1),res_helix_v1_c1,res_strand_v1_c1,res_hs_v1_c1)
  #cluster 2 stats
    total_v1=loop_v1_c2+helix_v1_c2+strand_v1_c2+hs_v1_c2
    if total_v1>0:
      print>>log,'%s,V1,%s,%s(%0.1f),%s(%0.1f),%s(%0.1f),%s,%s,%s'%(cluster_names[2],total_v1,loop_v1_c2,float(loop_v1_c2*100.0/total_v1),helix_v1_c2,float(helix_v1_c2*100.0/total_v1),strand_v1_c2,float(strand_v1_c2*100.0/total_v1),res_helix_v1_c2,res_strand_v1_c2,res_hs_v1_c2)
  #Cluster 3 stats
    total_v1=loop_v1_c3+helix_v1_c3+strand_v1_c3+hs_v1_c3
    if total_v1>0:
      print>>log,'%s,V1,%s,%s(%0.1f),%s(%0.1f),%s(%0.1f),%s,%s,%s'%(cluster_names[3],total_v1,loop_v1_c3,float(loop_v1_c3*100.0/total_v1),helix_v1_c3,float(helix_v1_c3*100.0/total_v1),strand_v1_c3,float(strand_v1_c3*100.0/total_v1),res_helix_v1_c3,res_strand_v1_c3,res_hs_v1_c3)
  #Cluster 4 stats
    total_v1=loop_v1_c4+helix_v1_c4+strand_v1_c4+hs_v1_c4
    if total_v1>0:
      print>>log,'%s,V1,%s,%s(%0.1f),%s(%0.1f),%s(%0.1f),%s,%s,%s'%(cluster_names[4],total_v1,loop_v1_c4,float(loop_v1_c4*100.0/total_v1),helix_v1_c4,float(helix_v1_c4*100.0/total_v1),strand_v1_c4,float(strand_v1_c4*100.0/total_v1),res_helix_v1_c4,res_strand_v1_c4,res_hs_v1_c4)
  if sel2:
    print>>log,'Cluster,Loop,Total,Loop(Per),Helix(Per),Strand(Per),HelixRes,Strand_Res,HelixStrandRes'
    total_v2=loop_v2_c0+helix_v2_c0+strand_v2_c0+hs_v2_c0
    if total_v2>0:
       print>>log,'%s,V2,%s,%s(%0.1f),%s(%0.1f),%s(%0.1f),%s,%s,%s'%(cluster_names[0],total_v2,loop_v2_c0,float(loop_v2_c0*100.0/total_v2),helix_v2_c0,float(helix_v2_c0*100.0/total_v2),strand_v2_c0,float(strand_v2_c0*100.0/total_v2),res_helix_v2_c0,res_strand_v2_c0,res_hs_v2_c0)
  #cluster 1 stats
    total_v2=loop_v2_c1+helix_v2_c1+strand_v2_c1+hs_v2_c1
    if total_v2>0:
      print>>log,'%s,V2,%s,%s(%0.1f),%s(%0.1f),%s(%0.1f),%s,%s,%s'%(cluster_names[1],total_v2,loop_v2_c1,float(loop_v2_c1*100.0/total_v2),helix_v2_c1,float(helix_v2_c1*100.0/total_v2),strand_v2_c1,float(strand_v2_c1*100.0/total_v2),res_helix_v2_c1,res_strand_v2_c1,res_hs_v2_c1)
  #cluster 2 stats
    total_v2=loop_v2_c2+helix_v2_c2+strand_v2_c2+hs_v2_c2
    if total_v2>0:
      print>>log,'%s,V2,%s,%s(%0.1f),%s(%0.1f),%s(%0.1f),%s,%s,%s'%(cluster_names[2],total_v2,loop_v2_c2,float(loop_v2_c2*100.0/total_v2),helix_v2_c2,float(helix_v2_c2*100.0/total_v2),strand_v2_c2,float(strand_v2_c2*100.0/total_v2),res_helix_v2_c2,res_strand_v2_c2,res_hs_v2_c2)
  #Cluster 3 stats
    total_v2=loop_v2_c3+helix_v2_c3+strand_v2_c3+hs_v2_c3
    if total_v2>0:
      print>>log,'%s,V2,%s,%s(%0.1f),%s(%0.1f),%s(%0.1f),%s,%s,%s'%(cluster_names[3],total_v2,loop_v2_c3,float(loop_v2_c3*100.0/total_v2),helix_v2_c3,float(helix_v2_c3*100.0/total_v2),strand_v2_c3,float(strand_v2_c3*100.0/total_v2),res_helix_v2_c3,res_strand_v2_c3,res_hs_v2_c3)
  #Cluster 4 stats
    total_v2=loop_v2_c4+helix_v2_c4+strand_v2_c4+hs_v2_c4
    if total_v2>0:
      print>>log,'%s,V2,%s,%s(%0.1f),%s(%0.1f),%s(%01f),%s,%s,%s'%(cluster_names[4],total_v2,loop_v2_c4,float(loop_v2_c4*100.0/total_v2),helix_v2_c4,float(helix_v2_c4*100.0/total_v2),strand_v2_c4,float(strand_v2_c4*100.0/total_v2),res_helix_v2_c4,res_strand_v2_c4,res_hs_v2_c4)
  log.close()
  print 'Done processing all clusters'
# get pdb with 
  v1_pdbs=[]
  #log=open(id+'_v1_bysecstruct.txt','w')
  input_list=[v1_5H_c0,v1_2H_c0,v1_3H_c0,v1_4H_c0,v1_5S_c0,v1_2S_c0,v1_3S_c0,v1_4S_c0,v1_5HS_c0,v1_2HS_c0,v1_3HS_c0,v1_4HS_c0,v1_5H_c1,v1_2H_c1,v1_3H_c1,v1_4H_c1,v1_5S_c1,v1_2S_c1,v1_3S_c1,v1_4S_c1,v1_5HS_c1,v1_2HS_c1,v1_3HS_c1,v1_4HS_c1,v1_5H_c2,v1_2H_c2,v1_3H_c2,v1_4H_c2,v1_5S_c2,v1_2S_c2,v1_3S_c2,v1_4S_c2,v1_5HS_c2,v1_2HS_c2,v1_3HS_c2,v1_4HS_c2,v1_5H_c3,v1_2H_c3,v1_3H_c3,v1_4H_c3,v1_5S_c3,v1_2S_c3,v1_3S_c3,v1_4S_c3,v1_5HS_c3,v1_2HS_c3,v1_3HS_c3,v1_4HS_c3,v1_5H_c4,v1_2H_c4,v1_3H_c4,v1_4H_c4,v1_5S_c4,v1_2S_c4,v1_3S_c4,v1_4S_c4,v1_5HS_c4,v1_2HS_c4,v1_3HS_c4,v1_4HS_c4]
  for l in input_list:
    if not l:continue
    #print>>log,l[0]
    v1_pdbs.append(l[0])
  log.close()
  v2_pdbs=[]
  #log=open(id+'_v2_bysecstruct.txt','w')
  input_list2=[v2_5H_c0,v2_2H_c0,v2_3H_c0,v2_4H_c0,v2_5S_c0,v2_2S_c0,v2_3S_c0,v2_4S_c0,v2_5HS_c0,v2_2HS_c0,v2_3HS_c0,v2_4HS_c0,v2_5H_c1,v2_2H_c1,v2_3H_c1,v2_4H_c1,v2_5S_c1,v2_2S_c1,v2_3S_c1,v2_4S_c1,v2_5HS_c1,v2_2HS_c1,v2_3HS_c1,v2_4HS_c1,v2_5H_c2,v2_2H_c2,v2_3H_c2,v2_4H_c2,v2_5S_c2,v2_2S_c2,v2_3S_c2,v2_4S_c2,v2_5HS_c2,v2_2HS_c2,v2_3HS_c2,v2_4HS_c2,v2_5H_c3,v2_2H_c3,v2_3H_c3,v2_4H_c3,v2_5S_c3,v2_2S_c3,v2_3S_c3,v2_4S_c3,v2_5HS_c3,v2_2HS_c3,v2_3HS_c3,v2_4HS_c3,v2_5H_c4,v2_2H_c4,v2_3H_c4,v2_4H_c4,v2_5S_c4,v2_2S_c4,v2_3S_c4,v2_4S_c4,v2_5HS_c4,v2_2HS_c4,v2_3HS_c4,v2_4HS_c4]
  for l in input_list2:
    if not l:continue
    #print>>log,l[0]
    v2_pdbs.append(l[0])
  log.close()
# SELECTION ONE AND TWO EXISTS
  v1v2_pdbs=[]
  start_v1v2=time.time()
  if sel1 and sel2 and len(processed_all)!=loaded:
    helix_v1v2=0;helixstrand_v1v2=0;strand_v1v2=0;res_hs_v1v2={};res_helix_v1v2={};res_strand_v1v2={}
    log=open('cluster_v1v2.csv','w')
    print>>log,'decoy,loop,sec_struct,score,v1_ss,v2_ss'
    start_v1v2=time.time()
    for model,score in main:
      if model not in processed_all:
        v1v2_pdbs.append(model)
        obj_name=extension.sub('',model)
        cluster={'sec':[]}
        cluster2={'sec':[]}
        cmd.iterate('%s & %s & %s'%(obj_name,sel1,name),'sec.append(ss)',space=cluster)
        cluster=''.join(cluster['sec'])
        cmd.iterate('%s & %s & %s'%(obj_name,sel2,name),'sec.append(ss)',space=cluster2)
        cluster2=''.join(cluster2['sec'])
        combo=cluster+cluster2
        if 'H' in cluster and 'S' not in cluster and 'H' in cluster2 and 'S' not in cluster2:#Helix v1 and v2 no strand
          helix_v1v2 +=1
          total_helix_v1 +=1
          total_helix_v2 +=1
          print>>log,model,',','V1V2',',','H',',',score,',',cluster,',',cluster2
          ss=[m.start() for m in re.finditer('H',str(combo))]
          for i in ss:
            if i<len(v1_loop):
              res_helix_v1v2[i+resid]=v1_loop[i]
            if i>=len(v1_loop):
              res_helix_v1v2[i+resid2]=v1v2[i]
  #STRAND CONFORMATION
        if 'S' in cluster and 'H' not in cluster and 'S' in cluster2 and 'H' not in cluster2:#Strand v1 and v2 
          strand_v1v2 +=1
          total_strand_v1 +=1
          total_strand_v2 +=1
          print>>log,model,',','V1V2',',','S',',',score,',',cluster,',',cluster2
          ss=[m.start() for m in re.finditer('H',str(combo))]
          for i in ss:
            if i <len(v1_loop):
              res_strand_v1v2[i+resid]=v1_loop[i]
            if i >=len(v1_loop):
              res_strand_v1v2[i+resid2]=v1v2[i]      
      #MIXED CONFORMATION
        if 'H' in cluster and 'S' not in cluster and 'S' in cluster2 and 'H' not in cluster2:# Helix v1 and strand V2
          helixstrand_v1v2 +=1
          total_strand_v2 +=1
          total_helix_v1 +=1
          print>>log,model,',','V1V2',',','HS',',',score,',',cluster,',',cluster2
          ss=[m.start() for m in re.finditer('H',str(cluster))]
          ss2=[m.start() for m in re.finditer('S',str(cluster2))]
          for i in ss:
            res_helix_v1v2[i+resid]=v1_loop[i]
          for i in ss2:
            res_strand_v1v2[i+resid2]=v2_loop[i]
        if 'S' in cluster and 'H' not in cluster and 'H' in cluster2 and 'S' not in cluster2:#Strand v1 helix v2
          helixstrand_v1v2 +=1
          total_strand_v1 +=1
          total_helix_v2 +=1 
          print>>log, model,',','V1V2',',','SH',',',score,',',cluster,',',cluster2
          ss=[m.start() for m in re.finditer('S',str(cluster))]
          ss2=[m.start() for m in re.finditer('H',str(cluster2))]
          for i in ss:
            res_strand_v1v2[i+resid]=v1_loop[i]
          for i in ss2:
            res_helix_v1v2[i+resid2]=v2_loop[i]
        if 'S' in cluster and 'H' in cluster and 'S' in cluster2 and 'H' in cluster2:#Helix and strand v1 and helix and strand v2
          helixstrand_v1v2 +=1
          total_strand_v1 +=1
          total_strand_v2 +=1
          total_helix_v1 +=0
          total_helix_v2 +=0
          print>>log,model,',','V1V2',',','HS&HS',',',score,',',cluster,',',cluster2
          ss=[m.start() for m in re.finditer('H|S',str(combo))]
          for i in ss:
            if i <len(v1_loop):
              res_hs_v1v2[i+resid]=v1_loop[i]
            if i >=len(v1_loop):
              res_hs_v1v2[i+resid2]=v1v2[i]   
        if 'H' in cluster and 'S' not in cluster and 'H' in cluster2 and 'S' in cluster2:# Helix v1 and Helix and strand V2
          helixstrand_v1v2 +=1
          total_strand_v2 +=1
          total_helix_v1 +=1
          total_helix_v2 +=1
          print>>log,model,',','V1V2',',','H&HS',',',score,',',cluster,',',cluster2
          ss=[m.start() for m in re.finditer('H',str(cluster))]
          ss2=[m.start() for m in re.finditer('H|S',str(cluster2))]
          for i in ss:
            res_helix_v1v2[i+resid]=v1_loop[i]
          for i in ss2:
            res_hs_v1v2[i+resid2]=v2_loop[i]
        if 'S' in cluster and 'H' not in cluster and 'H' in cluster2 and 'S' in cluster2:# Strand v1 and Helix and strand V2
          helixstrand_v1v2 +=1
          total_strand_v2 +=1
          total_strand_v1 +=1
          total_helix_v2 +=1
          print>>log,model,',','V1V2',',','S&HS',',',score,',',cluster,',',cluster2
          ss=[m.start() for m in re.finditer('S',str(cluster))]
          ss2=[m.start() for m in re.finditer('H|S',str(cluster2))]
          for i in ss:
            res_strand_v1v2[i+resid]=v1_loop[i]
          for i in ss2:
            res_hs_v1v2[i+resid2]=v2_loop[i]
        if 'H' in cluster and 'S' in cluster and 'H' in cluster2 and 'S' not in cluster2:# Helix and strand V1 helix v2
          helixstrand_v1v2 +=1
          total_strand_v1 +=1
          total_helix_v1 +=1
          total_helix_v2 +=1
          print>>log,model,',','V1V2',',','HS&H',',',score,',',cluster,',',cluster2
          ss=[m.start() for m in re.finditer('H|S',str(cluster))]
          ss2=[m.start() for m in re.finditer('H',str(cluster2))]
          for i in ss:
            res_hs_v1v2[i+resid]=v1_loop[i]
          for i in ss2:
            res_helix_v1v2[i+resid2]=v2_loop[i]
        if 'S' in cluster and 'H' in cluster and 'H' not in cluster2 and 'S' in cluster2:# Helix  and strand v1 and strand V2
          helixstrand_v1v2 +=1
          total_strand_v2 +=1
          total_strand_v1 +=1
          total_helix_v1 +=1
          print>>log,model,',','V1V2',',','HS&S',',',score,',',cluster,',',cluster2
          ss=[m.start() for m in re.finditer('H|S',str(cluster))]
          ss2=[m.start() for m in re.finditer('S',str(cluster2))]
          for i in ss:
            res_hs_v1v2[i+resid]=v1_loop[i]
          for i in ss2:
            res_strand_v1v2[i+resid2]=v2_loop[i] 
  time_v1v2=time.time()-start_v1v2
  log.close()
# PERFORM CLUSTER STATISTICS
  if len(v1v2_pdbs)>0:
    log=open(id+'_modeling_stat.csv','a')
    total_v1v2=helix_v1v2+strand_v1v2+helixstrand_v1v2
    print>>log,'All,v1v2,%s,%s(%0.1f),%s(%0.1f),%s(%0.1f),%s,%s,%s'%(total_v1v2,helix_v1v2,float(helix_v1v2*100.0/total_v1v2),strand_v1v2,float(strand_v1v2*100.0/total_v1v2),helixstrand_v1v2,float(helixstrand_v1v2*100.0/total_v1v2),res_helix_v1v2,res_strand_v1v2,res_hs_v1v2)
    log.close()
    print 'Done processing v1v2'
    # store global modeling statistics i
  total_V1=total_loop_v1+total_helix_v1+total_strand_v1+total_hs_v1
  total_V2=total_loop_v2+total_helix_v2+total_strand_v2+total_hs_v2
  log=open(id+'_modeling_stat.csv','a')
  print>>log,'Cluster,Loop,Total,Total:Loops(Per),Total:Helix(Per),Total:Strand(Per),Total:HS(Per)'
  if not sel2:
    per_loop_v1='%0.1f'%(total_loop_v1*100.0/total_V1)
    per_helix_v1='%0.1f'%(total_helix_v1*100.0/total_V1)
    per_strand_v1='%0.1f'%(total_strand_v1*100.0/total_V1)
    per_hs_v1='%0.1f'%(total_hs_v1*100.0/total_V1)
    print>>log,'All,V1,%s,%s(%s),%s(%s),%s(%s),%s(%s)' %(total_V1,total_loop_v1,per_loop_v1,total_helix_v1,per_helix_v1,total_strand_v1,per_strand_v1,total_hs_v1,per_hs_v1)
  if not sel1:
    per_loop_v2='%0.1f'%(total_loop_v2*100.0/total_V2)
    per_helix_v2='%0.1f'%(total_helix_v2*100.0/total_V2)
    per_strand_v2='%0.1f'%(total_strand_v2*100.0/total_V2)
    per_hs_v2='%0.1f'%(total_hs_v2*100.0/total_V2)
    print>>log,'All,V2,%s,%s(%s),%s(%s),%s(%s),%s(%s)' %(total_V2,total_loop_v2,per_loop_v2,total_helix_v2,per_helix_v2,total_strand_v2,per_strand_v2,total_hs_v2,per_hs_v2)
  if sel1 and sel2:
    per_loop_v1='%0.1f'%(total_loop_v1*100.0/total_V1)
    per_helix_v1='%0.1f'%(total_helix_v1*100.0/total_V1)
    per_strand_v1='%0.1f'%(total_strand_v1*100.0/total_V1)
    per_hs_v1='%0.1f'%(total_hs_v1*100.0/total_V1)
    per_loop_v2='%0.1f'%(total_loop_v2*100.0/total_V2)
    per_helix_v2='%0.1f'%(total_helix_v2*100.0/total_V2)
    per_strand_v2='%0.1f'%(total_strand_v2*100.0/total_V2)
    per_hs_v2='%0.1f'%(total_hs_v2*100.0/total_V2)
    print>>log,'All,V1,%s,%s(%s),%s(%s),%s(%s),%s(%s)' %(total_V1,total_loop_v1,per_loop_v1,total_helix_v1,per_helix_v1,total_strand_v1,per_strand_v1,total_hs_v1,per_hs_v1)
    print>>log,'All,V2,%s,%s(%s),%s(%s),%s(%s),%s(%s)' %(total_V2,total_loop_v2,per_loop_v2,total_helix_v2,per_helix_v2,total_strand_v2,per_strand_v2,total_hs_v2,per_hs_v2)
  log.close()
  print 'Running Score vs RMSD plots...'
  start_align=time.time()
  #score_vs_rmsd()
  end_align=time.time()-start_align
  def image_getter():
    cmd.reinitialize()
    cmd.bg_color('white')
    cmd.set('depth_cue',0)
    cmd.set("antialias",1)
    cmd.set("ray_trace_mode",0)
    cmd.set("depth_cue",0)
    cmd.set("ray_trace_fog",0)
    cmd.set("stick_radius", 0.1)
    cmd.set('ray_shadows',0)
    cmd.set('label_color','black')
    #cmd.set("cartoon_side_chain_helper",1)
    #cmd.set("cartoon_flat_sheets",1)
    if len(v1_pdbs)>0:
      nterm,cterm=sel1[3:].split('-')
      if not os.path.exists('v1'):
        os.makedirs('v1')
      for pdb in v1_pdbs:
        target=v1_pdbs[0]
        cmd.load(target)
        target=extension.sub('',target)
        obj=extension.sub('',pdb)
        cmd.load(pdb,obj)
        cmd.hide('lines')
        cmd.show('cartoon','%s & %s' %(obj,sel1))
        cmd.align('%s & %s & %s'%(obj,bb,sel1),'%s & %s & %s'%(target,bb,sel1))
        cmd.zoom('%s'%sel1)
      cmd.save('v1/'+id+'_v1_pdbs.pse')
      cmd.delete('all')
    if len(v2_pdbs)>0:
      nterm,cterm=sel2[3:].split('-')
      if not os.path.exists('v2'):
        os.makedirs('v2')
      for pdb in v2_pdbs:
        target=v2_pdbs[0]
        target=extension.sub('',target)
        obj=extension.sub('',pdb)
        cmd.load(pdb,obj)
        cmd.hide('lines')
        cmd.show('cartoon','%s & %s' %(obj,sel2))
        cmd.align('%s & %s & %s'%(obj,bb,sel2),'%s & %s & %s'%(target,bb,sel2))
        cmd.zoom('%s'%sel2)
      cmd.save('v2/'+id+'_v2_pdbs.pse')
      cmd.delete('all')
    if len(v1v2_pdbs)>0:
      for pdb in v1v2_pdbs:
        target=v1v2_pdbs[0]
        target=extension.sub('',target)
        obj=extension.sub('',pdb)
        cmd.load(pdb,obj)
        cmd.hide('lines')
        cmd.show('cartoon','%s & %s' %(obj,sel2))
        cmd.align('%s & %s & %s'%(obj,bb,sel2),'%s & %s & %s'%(target,bb,sel2))
        cmd.zoom('%s'%sel1)
        #cmd.png('v2/'+ id+ '_' + obj+ '_v1v2',3600,2400,dpi=300,ray=1)
      cmd.save('v2/'+id+'_v1v2_pdbs.pse')
      cmd.delete('all')
  #image_getter()
  end = float(time.time())
  Total_time='%0.2f'%(float(end-start))
  print 'Timings'
  print 'Total time: %ss' % Total_time
  print '  Plotting: %0.2fs' % float(end_align)
  if cluster0:
   print '  Cluster 0: %0.2fs' % float(time_c0)
  if cluster1:
   print '  Cluster 1: %0.2fs' % float(time_c1)
  if Cluster2:
   print '  Cluster 2: %0.2fs' % float(time_c2)
  if cluster3:
   print '  Cluster 3: %0.2fs' % float(time_c3)
  if cluster4:
   print '  Cluster 4: %0.2fs' % float(time_c4)
  if sel1 and sel2:
   print '  Cluster v1v2: %0.2fs' % float(time_v1v2)

cmd.extend('analysis',loop_model_analysis)
