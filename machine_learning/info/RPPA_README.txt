~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DREAM7_DrugSensitivity1_RPPA.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
Reverse protein lysate array (RPPA) data is a method to quantify
protein abundance.  Cells were exponentially growing prior to
harvesting, and not subject to any particular treatment (i.e., these
are baseline data).

Note that for each protein, there can be multiple isoforms (eg, native
and phosphorylated forms).  Also, an antibody designed for a protein
can pull down paralogs.  These cases have been encoded as follows:

- a phosphorylated protein is designated as the protein name followed
  by a "p."  Additionally, if the phosphorylated residue is know, this
  is listed.  For example, the phosphorylated version of MYC is
  represented as MYC-p. Residue 73 is phosphorylated on the protein
  JUN and is represented as JUN-p73.
- Instances where an antibody pulls down multiple protein paralogs is
  represented by a comma-separated list of these proteins.  For
  example, the antibody for AKT pulls down AKT1, AKT2, and AKT3, which
  is represented as "AKT1,AKT2,AKT3"


The second column, titled "FullyValidated," indicates whether the
antibody used in the assay was fully validated.  Every antibody on the
RPPA platform undergoes rigorous validation to ensure that there is no
cross-reactivity, and that the antibodies behave in a linear manner
following the dilution series.  Data derived from antibodies that were
not fully validated should be interpreted with caution.  The most
robust dataset would include only those data from fully validated
antibodies.


-- RPPA data generated by Gordon Mills Lab at MD Anderson in 2011
-- Created by Laura Heiser on 1/1/11.
-- Edited by Jim Costello for the DREAM7 Drug Sensitivity challenge on 6/12/12.
