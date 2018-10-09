Tools to manipulate PDB files
=============================

General PDB files need tools to apply transformations to them. As this need is
quite common there a number of tools out there already and where applicable
you should use them. Examples of such tools are:

* convpdb.pl
* pdb4amber
* charmmlipid2amber.py
* pdbfixer

Unfortunately some of these tools try to be too helpful. For example `convpdb.pl`
will not only scale the coordinates in a PDB file if you ask for it, but it will
also drop all TER records, rename HOH to TIP3, drop the CRYST1 record, and drop
the chemical symbol at the end of the line. This might leave you with a useless
result. 

The tools in this directory do exactly what it says on the tin, nothing less
and nothing more.
