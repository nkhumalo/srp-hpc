# 1ai0: R6 Insulin hexamer

Insulin has been found in 3 different conformations: T6 --> T3R3 --> R6.
To obtain insulin in the R6 conformation there needs to be Phenol
present as well as a high concentration of Chloride ions. It is also
well known that the presence of Zn(2+) stabelizes the hexamer, as
this has been used to create slow acting insulin in the past.

In this directory we start from the 1ai0 entry in the Protein
Data Bank. The Phenol groups have been removed, the system solvated
at physiological pH of 7.4, and Na(+) counter ions added to 
neutralize the charge. Subsequently, a dynamics simulation has
been run to find a conformation change.

Files:

* 1ai0.pdb

  * Original PDB file.

* 1ai0-fixed.pdb

  * Removed Phenol groups.

* 1ai0-fixed-renum.pdb

  * Used `pdb4amber -i 1ai0-fixed.pdb -o 1ai0-fixed-renum.pdb` to fix
    residue numbers.

* 1ai0-fixed-renum-nw.pdb

  * Converted Amber atom labels to NWChem atom labels.

