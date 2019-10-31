# Protein-modeling and engineering

Python and Bash scripts utilized in conducting research in on HIV antibody prediction, engineering and HIV sequences during my thesis. 
The project involved modeling and engineering antibody structure computationallyusing ROSETTA macromodeling software

# BRANCHES
## protocols
Contains Python and Bash scripts used to generate folders and files used to run ROSETTA executables for structure prediction, preparation, design and analyses in SGE and SLURM environment 
## analyses
Contains code that that clusters,aligns, and perform secondary structure analysis of the generated protein models, analyze the effect of mutation locally and globally. Each code has data visualization plots 

## NGS--next generation sequencing
antibody_pairing_analysis.py generates DNA and Protein fasta files of paired antibodies from next generation sequencing project. Initial setup for Linux and MacOS for the required python modules included as well as igblast executables,imgt germlines and optional files

## design
Contains code that aid in generating input files necessary to perform design or mutation on the protein structure of interest 
