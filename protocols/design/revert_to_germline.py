#!/usr/bin/python
########################################################################
# Makes input file for flex and fixed ddg
# created 6.15.18 by Edwin Kamau 
########################################################################
# import all the necessary functions required to make functional calls
import  sys, os
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_one
import glob
import os 
import numpy as np 

if len(sys.argv)<3:
     print ('requires 3 arguments %s given\nUsage : revert_to_germline.py output_name fastafile\nname : name of the resfile to be created\nfastafile : fasta formated file that contains the heavy and light chain germline amino acid sequences. heavy chain first and light chain second'  % (len(sys.argv[1:]))) 
     sys.exit(0)

output_name=sys.argv[1]
fastafile=sys.argv[2]

def mutate(name,vh,vl):
    '''
    creates a resfile that can be utilized to convert mature antibody into germline. Assumes antibody crystal structure exists in the current directory and is chothia numbered  
    '''
    parser=PDBParser()
    structure=parser.get_structure('X',glob.glob('*pdb')[0])
    f=open(name+'_2_germ.resfile','w')
    f.write('NATAA\nstart\n\n')
    for model in structure:
        for chain in model:
            if chain.id=='H':
                for index,res in enumerate(chain):
                    if index <len(vh):
                        if three_to_one(res.resname)!=vh[index]:
                            f.write('%s %s PIKAA %s\n' %(''.join(str(res.id[1])+res.id[-1]),chain.id,vh[index]))
                            print('%s %s PIKAA %s' %(''.join(str(res.id[1])+res.id[-1]),chain.id,vh[index]))
            if chain.id=='L':
                for index,res in enumerate(chain):
                    if index <len(vl):
                        if three_to_one(res.resname)!=vl[index]:
                            f.write('%s %s PIKAA %s\n' %(''.join(str(res.id[1])+res.id[-1]),chain.id,vl[index]))
                            print('%s %s PIKAA %s' %(''.join(str(res.id[1])+res.id[-1]),chain.id,vl[index]))
    f.close()

def parse_fasta_file(output_name,fastafile):
    ''' Parse FASTA file and return list of sequences it contains as a list that can be fed to mutate protocol above
    '''
    try:
        sequences = []    
        with open(fastafile) as f:
            for line in f:
                if line.startswith('>'): sequences.append('')
                else: sequences[-1] += line[:-1]  # removing end-line
        
        if sequences:
           mutate(output_name,sequences[0],sequences[1])
    except Exception as e:
        print(e)
        sys.exit(0)

parse_fasta_file(output_name,fastafile)

print('Done making  input file')


