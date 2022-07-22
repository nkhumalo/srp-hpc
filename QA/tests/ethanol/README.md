# Ethanol test case

This is a simple test case of a single ethanol molecule in water.
The PDB file contains the entire set of atoms included in the simulation
including both solute and solvent waters.

## Analysis configuration file

For the Online Adaptive Sampler (OAS) it is convenient to have a 
configuration file that identifies relevant groups of atoms. As the
molecular system for the MD simulation is defined in a PDB file it seems
sensible to include that information as annotation. So the configuration
file is simply a copy of the PDB file where each line is prepended with one
of the following:

- "#" a comment (PDB files may contain a lot of information that does specify 
  atom positions
- "0" an unimportant atom
- _n_ a relevant atom belonging to group _n_ of relevant atoms

For example in the `ethanol.oas` configuration file in this directory the
atoms of the ethanol molecule are identified as belonging to group 1.

As another example one could imagine protein-ligand docking and one might
identify the protein atoms that delineate the binding pocket as belonging to
group 1, and the ligand atoms as belonging to group 2.

Even more complicated scenarios are possible, such as a molecule that binds
in a pocket that is formed at the interface where two proteins interact.
The pocket delineating atoms of the first protein could labeled as belonging
to group 1, those of the second protein as group 2, and the ligand atoms 
as group 3.
