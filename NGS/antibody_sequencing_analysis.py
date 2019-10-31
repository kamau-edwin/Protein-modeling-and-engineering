#!/usr/bin/python
########################################################################
# antibody sequencing analysis, pairing setup igblast imgt 
#
# Written by: Edwin Kamau 
# Last modified: 6.28.2018 by Edwin Kamau
#
########################################################################
# import all the necessary functions required to make functional calls

import pandas as pd
from collections import Counter
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys
import os
import re
import subprocess
import platform
import glob 

if len(sys.argv)!=6 :
    print ('Usage: abtool.py  <path>  <output_name> <project> <species> <common_name> <receptor_type <igblast_version>\npath: NGS sequencing output file in JSON format eg from 10x \noutput_name:prefix that will be used to name the output defaults to output\nproject:name of the project eg. HIV,...\nspecies: Species of interest eg. Homo sapiens\ncommon_name:common name of the species of interest e.g humans\nreceptor: IG or TR\nigblast release version one want to use set up at '1,8.0' \n' )
    sys.exit(0)

path=sys.argv[1]
outname=sys.argv[2]
project=sys.argv[3]
working_dir=os.getcwd()
set1='imgt'
sp=sys.argv[4] # Homo sapiens
cn=sys.argv[5] # common name humans
receptor=sys.argv[6] # 'IG/TR'
release=sys.argv[7]
if not release:
    release='1.8.0' # version of igblast 


def abtool(path,outname='output'):
    '''
    function requires file to be analyzed. Assumes a particular annotation format thus may need to be modified if the format changes 
    if outname is not given it will use output as the prefix
    outputs: DNA and protein fasta sequences as single or multiple paired
             single and multiple paired chains as csv format  
    '''
    try:
        if path:
            # load file
            db=pd.read_json(path,orient='records')
            # REMOVE ALL NON-PAIRED chains
            df=db.copy()
            df=df[df.duplicated('clonotype',keep=False)]
            # get chain type, gene name and gene start and end from annotation
            for rec in df.index:
                sequence=str(df.loc[rec]['sequence']) 
                df.loc[rec,'chain']=df.loc[rec,'annotations'][0]['feature']['chain']# 
                gene_id=list(set([d['feature']['gene_name'] for d in df.loc[rec,'annotations']]))# grab gene names 
                position=[(s['contig_match_start'],s['contig_match_end']) for s in df.loc[rec,'annotations']]
                combine={} # Combine gene names and start and end and filter 5' UTR 
                for ids in range(len(gene_id)):
                    if gene_id[ids] not in combine:
                        combine[gene_id[ids]]=position[ids]
                for key,value in combine.items():
                    if 'V' in key:
                        df.loc[rec,'IGHV']=key
                        df.loc[rec,'v_start']=value[0]
                        df.loc[rec,'sequence']=sequence[value[0]:]# remove 5' UTR
                        df.loc[rec,'v_end']=value[1]
                    if 'D' in key and len(key)>4:
                        df.loc[rec,'IGHD']=key
                        df.loc[rec,'d_start']=value[0]
                        df.loc[rec,'d_end']=value[1]
                    if 'J' in key:
                        df.loc[rec,'IGHJ']=key
                        df.loc[rec,'j_start']=value[0]
                        df.loc[rec,'j_end']=value[1]
                    if len(key)==5 and key[3]!='J' or len(key)==4:
                        df.loc[rec,'type']=key
                        df.loc[rec,'t_start']=value[0]
                        df.loc[rec,'t_end']=value[1]
            columns=['clonotype','contig_name','chain','IGHV','IGHD','IGHJ','type','cdr3','v_start','v_end','d_start','d_end','j_start','j_end','t_start','t_end','sequence','aa_sequence']# get useful columns
            # set index to clonotype makes it easier to perform join operations
            df=df[columns]
            df=df.set_index('clonotype') # change index to 'clonotype'
            # separate single paired and multiple paired abs 
            single_pair=[key for key,value in Counter(df.index).items() if value==2]
            multiple_pair=[key for key,value in Counter(df.index).items() if value>2]
            # PAIRING SINGLE and MUltiple
            single_paired=df[df.index.isin(single_pair) & df.chain.isin(['IGH'])].join(df[df.index.isin(single_pair) &\
                  df.chain.isin(['IGL','IGK'])],how='outer', rsuffix='_lc').dropna()
            multiple_paired=df[df.index.isin(multiple_pair) & df.chain.isin(['IGH'])].join(df[df.index.isin(multiple_pair) \
                 & df.chain.isin(['IGK','IGL'])],how='outer',rsuffix='_lc').join(df[df.index.isin(multiple_pair) &\
                 df.chain.isin(['IGK','IGL'])],how='outer',rsuffix='_lc2').dropna()
            # remove duplicated rows
            multiple_paired=multiple_paired[multiple_paired.v_gene_lc !=multiple_paired.v_gene_lc2].drop_duplicates('cdr3')
            # Assign columns to output 
            create_fasta_files(df,single_paired,multiple_paired,project,outname)
            # output to csv format 
            multiple_paired[multiple_col].to_csv(outname + '_multiple_paired.csv')
            single_paired[single_col].to_csv(outname + '_single_paired.csv')
    except Exception as e:
        print(e)
        sys.exit()

def create_fasta_files(dataframe,single,multiple,project,outname='output'):
    '''
    Generates DNA fasta sequences from created dataframes 
    Assumes dataframe is indexed using contig_name and the antibody chains have been paired as single or multiple
    outputs one fasta files
    '''
    try: 
        j=1
        df=dataframe.copy()
        df=df.set_index('contig_name') # 
        f=open(output_name +'_dna_fasta.txt','w')
        #f2=open(output_name +'_protein_fasta.txt','w')
        for index in df.index:
            gene=df.loc[index,'IGHV']
            end=int(df.loc[index,'j_end'])
            dna=str(Seq(df.loc[index,'sequence'][:end]))#ENDS at end of CDR h3
            #prot=df.loc[index,'aa_sequence']
            if index in list(single.contig_name.values):
                #f2.write('>sp|%s00%s|%s|%s|HC\n%s\n' %(project,j,index,gene,prot))
                f.write('>sp|%s00%s|%s|%s|HC\n%s\n' %(project,j,index,gene,dna))
                j+=1
            if index in list(single.contig_name_lc.values):
                #f2.write('>sp|%s00%s|%s|%s|LC\n%s\n' %(project,j,index,gene,prot))
                f.write('>sp|%s00%s|%s|%s|LC\n%s\n' %(project,j,index,gene,dna))
                j+=1
            if index in list(multiple.contig_name.values):
                #f2.write('>mp|%s00%s|%s|%s|HC\n%s\n' %(project,j,index,gene,prot))
                f.write('>mp|%s00%s|%s|%s|HC\n%s\n' %(project,j,index,gene,dna))
                j+=1
            if index in list(multiple.contig_name_lc.values):
                #f2.write('>mp|%s00%s|%s|%s|LC\n%s\n' %(project,j,index,gene,prot))
                f.write('>mp|%s00%s|%s|%s|LC\n%s\n' %(project,j,index,gene,dna))
                j+=1
            if index in list(multiple.contig_name_lc2.values):
                #f2.write('>mp|%s00%s|%s|%s|LC2\n%s\n' %(project,j,index,gene,prot))
                f.write('>mp|%s00%s|%s|%s|LC2\n%s\n' %(project,j,index,gene,dna))
                j+=1
        f.close()
        #f2.close()
    except Exception as e:
        print(e)

def runigblast(igblastcmd,germpath)
    '''
    generate igblast file 
    '''
    file_root='imgt'+ '_' + cn + '_' + receptor
    vfile = germpath + '/' + file_root + 'V'
    dfile = germpath + '/' + file_root + 'D'
    jfile = germpath + '/' + file_root + 'J'
    seqfile=glob.glob(working_dir + '/' + outname'*dna.fasta')
    logfile=outname + 'imgt_mapping.log'

    cmd = '%s -germline_db_V %s -germline_db_D %s -germline_db_J %s -organism %s -domain_system %s  -query %s -auxiliary_data %s/%s_gl.aux -outfmt 3 >%s' % (igblastcmd, vfile, dfile, jfile, cn, set1, seqfile, optional,cn, logfile)
    print 'Running command %s' % cmd
    subprocess.call(cmd, shell=True)

    print('Finished mapping sequences to imgt')

def imgt(igblastpath,igblastcmd):
    '''
    downloads imgt germline sequences from imgt and makes database base on specified organism
    '''
    # make germline folder path
    try:
        germpath=igblastpath + '/%s' % set1

        def make_db():
            recs = list(SeqIO.parse('%s/imgt_germlines.fasta' % germpath, 'fasta'))
            fastafiles = []
            outrecs = {}
            for rec in recs:
                d = rec.description.split('|')
                if len(d) > 4 and receptor in d[1] and sp in d[2] and d[4] in ['V-REGION', 'D-REGION', 'J-REGION']:
                    seg = d[4][0]
                    rec.description = d[1]
                    rec.id = d[1]
                    key = d[1][:2] + seg
                    if key not in outrecs:
                        outrecs[key] = []
                    outrecs[key].append(rec)
            for fn, r in outrecs.items():
                fastafile = '%s/imgt_'% germpath + cn + '_' + fn + '.fasta'
                SeqIO.write(r, fastafile, 'fasta')
                fastafiles.append(fastafile)
            for fn in fastafiles:
                print('Processing germline file %s' % fn)
                cmd = '%s/makeblastdb -parse_seqids -dbtype nucl -in %s -out %s' % (igblastpath, fn, fn.replace('.fasta', ''))
                subprocess.call(cmd, shell=True)   
            print('Germline file processing complete.')
            #subprocess.call('touch %s/complete' % path, shell=True)

        if os.path.isfile(germpath + '/%s_germlines.fasta' %(set1)):
            if os.path.isfile(germpath + '/%s_%s_%sV.nhr' %(set1,cn,receptor)):
                print('imgt_germline file previously downloaded and processed for species : %s' % cn)
            if not os.path.isfile(germpath + '/%s_%s_%sV.nhr' %(set1,cn,receptor)):
                print('making %s %s database for' % (cn,set1))
                make_db()

        if not os.path.exists(germpath):
            os.makedirs(germpath)
            print ('Installing germline files in  %s.' % germpath)

            subprocess.call("wget -P %s -O %s/imgt_germlines.fasta http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP" % (germpath,germpath), shell=True)
            make_db()
    except exception as e:
        print(e)
    abtool(path,outname)
    runigblast(igblastcmd,germpath)

def igblast():
    '''
    sets up igblast executables in home director and the required files internal_data, optional files and imgt database in $HOME/igblast
    '''
    # dependent files generator
    def dependent(igblastpath,folder):
        '''
        creates dependent files for igblast 
        '''
        path=igblastpath +'/' + folder
        if not os.path.exists(path):
            os.makedirs(path)
            print ('Igblast executables installed but no %s dependency files\nInstalling the dependency' % folder)
            subprocess.call("wget -P %s -nH --cut-dirs=5 -r ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/%s/ " %(path,folder),shell=True)
        if 'internal_data' in path:
            print('linking %s to %s as %s' % (path,working_dir,folder))
            subprocess.call('ln -s %s %s/%s' % (path,working_dir,folder),shell=True)
            print ('Finished installing %s files' % folder)

    home=os.environ['HOME']

    igblastpath=home + '/ncbi-igblast-' + release + '/bin/' # path of the executables

    igblastcmd=igblastpath + '/igblastn' # executable file

    if os.path.exists(igblastpath):
        if os.path.isfile(igblastcmd): # is executable installed
            if os.path.exists(igblastpath + '/internal_data'): # is dependent file installed
                if os.path.exists(igblastpath + '/optional_file'): # is optional file installed
                    print('igblast executables and dependent files previously installed')
                if not os.path.exists(igblastpath + '/optional_file'):
                    dependent(igblastpath,'optional_file')
            if not os.path.exists(igblastpath + '/internal_data'): # are dependent file installed
                dependent(igblastpath,'internal_data')

    if not os.path.exists(igblastpath): # no igblast executables
        print('Installing igblast executables and dependent files in the %s directory' % home)

        url= 'ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/%s/ncbi-igblast-%s-x64-%s.tar.gz' %(release,release,system)

        extension=url.split('/')[-1]   

        subprocess.call("wget -P %s %s ; tar -xf %s/%s -C %s" %(home,url,home,extension,home),shell=True)
        for depend in ['internal_data','optional_file']:
            dependent(igblastpath,depend)

    imgt(igblastpath,igblastcmd) # pass to imgt to make db a

def setup():
# check platform 
    try:
        if platform.system()=='Darwin':
            system='macosx'
        if platform.system()=='Linux':
            system='linux'
        brew=ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
        # check and install dependencies
        for prog in ['xcode-select','brew','python3',]:
            if subprocess.call('which %s' % prog,shell=True)==0:continue
                print('%s already installed'% prog)
                if prog == 'brew':
                    print('Updating %s please wait\n' % prog)
                    subprocess.call('%s update' % prog,shell=True)
            else:
                if prog == 'xcode-select':
                    print('installing %s\will require permission to install\nplease respond to dialog boxes that will popup\n Password: The password used to login into the machine otherwise you should get the machine owner/institution and rerun the script to install %s' %(prog,prog))
                    subprocess.call('%s --install',shell=True)
                    print('%s installed' % prog )
                if prog =='brew':
                    print('installing %s' % prog)
                    subprocess.call('%s' % brew, shell=True)
                if prog.startswith('python'):
                    print('installing %s' % prog)
                    if subprocess.call('which brew',shell=true)==0:
                        subprocess.call('brew install %s' % prog, shell=True) # should install pip3 as well
                    else:
                        print('brew did not install properly\nrun the following in terminal to install brew\n%s\n' % brew)
                    # install python dependent packages
                    if subprocess.call('which pip3',shell=True)==0:
                        subprocess.call('pip3 install scipy numpy pandas biopython',shell=True)

    except exception as e:
        print(e)
    print('Your system has the required modules to run this program\nPlease update the modules regularly')
    igblast() # make igblast executables and prepare dependency file 

def main()
    """
    finds files in the folder that has the file extension 
    """
    try:
        for root,folder,files in os.walk(os.environ['HOME']):
            for file in files:
                if 'igblastn' in file and 'ncbi' in root and release in root:
                    igblastpath=root
                    igblastcmd=os.path.join(root,file)
                if cn +'_gl.aux' in file and 'optional_file' in root:
                    optional=root
                if cn + '_' + receptor + 'V.fasta' in file and set1 in root:
                    germpath=root
                if cn + '_V.phr' in file and 'internal_data' in root:
                    internal=root
        if igblastcmd and optional and internal and germpath: # previous download good to go
            abtool(path,outname)
            runigblast(igblastcmd,germpath)
        else:
            setup()
    except exception as e:
        print(e)

main()



# BELOW IS HOW TO IMPREMENT NOT JSON FILE 
#sequence=str(df.loc[rec]['sequence']) # grab all clonal dna sequence
#aa=str(df.loc[rec]['aa_sequence'])# grab each clonal prot sequence
#codons=[m.start() for m in re.finditer('ATG',sequence)]# find all start codons
#for codon in codons:
  #  prot_seq=Seq(sequence[codon:],IUPAC.ambiguous_dna).translate()#translate dna
  #  if str(prot_seq)==aa: # find codon that match the given prot_sequence
   #     df.loc[rec,'sequence']=sequence[codon:] # update the sequence by removing 5' UTR sequence.

   #df.drop(['annotations','barcode','cdr3_seq','filtered','frame','high_confidence','info','is_cell'
   #'primer_annotations','productive','quals','read_count','start_codon_pos',
 #'stop_codon_pos','umi_count'],axis=1,inplace=True)
