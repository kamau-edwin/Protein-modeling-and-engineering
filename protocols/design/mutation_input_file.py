#!/usr/bin/python
########################################################################
# Makes input file for ddg cartesian protocol 
########################################################################
# import all the necessary functions required to make functional calls
import  sys, os


if len(sys.argv)<2:
     print 'requires 2 arguments %s given' % (len(sys.argv[1:]))
     print 'Usage : mutation_input_file.py  fasta file'
     sys.exit(0)

file_name=sys.argv[1]

def parse_fasta_file(file_name):
    ''' Parse FASTA file and return list of sequences it contain
    assumes all chains are merged into a single fasta file. useful for creating
    an input file for cartesian ddg application.
    fasta format should be pdb sequence at position one and mutation sequence at
    position two. Sequences need to be of equal length or indexerror will occur 
    or need to handle different length with an exception 
    '''
    try:
        sequences = []    
        with file(file_name) as f:
            for line in f:
                if line.startswith('>'): sequences.append('')
                else: sequences[-1] += line[:-1]  # removing end-line
        
        if sequences:
            f=open('input.txt','w')
            for pos in range(len(sequences[0])):
                if sequences[0][pos]!=sequences[1][pos]:
                    #print ('1 {0} {1} {2}'.format(sequences[0][pos],pos+1,sequences[1][pos]))
                    print ('1 %s %s %s' %(sequences[0][pos],pos+1,sequences[1][pos]))
                    #f.write('1 {0} {1} {2}\n'.format(sequences[0][pos],pos+1,sequences[1][pos]))
                    f.write('1 %s %s %s\n'%(sequences[0][pos],pos+1,sequences[1][pos]))
            f.close()
        else:
            print('No fasta file with sequences provided')
            sys.exit(0)
        
    except IndexError:
        print('sequences have different length')
        os.remove('input.txt')
        sys.exit(0)

parse_fasta_file(file_name)

print('Done making  input file')


