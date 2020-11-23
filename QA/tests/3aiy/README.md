# 3AIY: Human Insulin Hexamer

This directory contains data for a simulation on Human Insulin Hexamer.
The idea is to separate the 6 proteins by a little bit and see if they
recombine in a molecular dynamics simulation. The name of this 
directory corresponds to the PDB identifier of the protein complex.
As the structure was determined by NMR it conveniently contains all 
atoms.

The structure was processed according to following steps:

* 3aiy.pdb - the original structure the Protein Data Bank

* 3aiy-fixed.pdb - removed the phenol groups present in the structure
  (simply deleting all "HETATM" lines, and the corresponding "CONNCT" lines)

* 3aiy-fixed-scaled.pdb - scaled the atom coordinates by a factor of 2
  to separate the insulin proteins (use the pdb_scale program from 
  nwchem/contrib/pdb_tools).

* 3aiy-fixed-scaled-nw.pdb - translate the PDB to the NWChem dialect using
  
  * pdb_amber2nwchem -i 3aiy-fixed-scaled.pdb -o 3aiy-fixed-scaled-nw.pdb

