#!/bin/bash
../../../bin/LINUX64/nwchem 1nt5_prep1.nw    > 1nt5_prep1.out
../../../bin/LINUX64/nwchem 1nt5_eq.nw       > 1nt5_eq1.out
../../../bin/LINUX64/nwchem 1nt5_md.nw       > 1nt5_md1.out
../../../bin/LINUX64/nwchem 1nt5_md_ana.nw   > 1nt5_md_ana1.out
mv 1nt5_md.pdb 1nt5_md1.pdb
mv 1nt5_md.crd 1nt5_md1.crd
#
../../../bin/LINUX64/nwchem 1nt5_prep2.nw    > 1nt5_prep2.out
../../../bin/LINUX64/nwchem 1nt5_eq.nw       > 1nt5_eq2.out
../../../bin/LINUX64/nwchem 1nt5_md.nw       > 1nt5_md2.out
../../../bin/LINUX64/nwchem 1nt5_md_ana.nw   > 1nt5_md_ana2.out
mv 1nt5_md.pdb 1nt5_md2.pdb
mv 1nt5_md.crd 1nt5_md2.crd
#
../../../bin/LINUX64/nwchem 1nt5_prep3.nw    > 1nt5_prep3.out
../../../bin/LINUX64/nwchem 1nt5_eq.nw       > 1nt5_eq3.out
../../../bin/LINUX64/nwchem 1nt5_md.nw       > 1nt5_md3.out
../../../bin/LINUX64/nwchem 1nt5_md_ana.nw   > 1nt5_md_ana3.out
mv 1nt5_md.pdb 1nt5_md3.pdb
mv 1nt5_md.crd 1nt5_md3.crd
#
../../../bin/LINUX64/nwchem 1nt5_prep4.nw    > 1nt5_prep4.out
../../../bin/LINUX64/nwchem 1nt5_eq.nw       > 1nt5_eq4.out
../../../bin/LINUX64/nwchem 1nt5_md.nw       > 1nt5_md4.out
../../../bin/LINUX64/nwchem 1nt5_md_ana.nw   > 1nt5_md_ana4.out
mv 1nt5_md.pdb 1nt5_md4.pdb
mv 1nt5_md.crd 1nt5_md4.crd
#
../../../bin/LINUX64/nwchem 1nt5_prep5.nw    > 1nt5_prep5.out
../../../bin/LINUX64/nwchem 1nt5_eq.nw       > 1nt5_eq5.out
../../../bin/LINUX64/nwchem 1nt5_md.nw       > 1nt5_md5.out
../../../bin/LINUX64/nwchem 1nt5_md_ana.nw   > 1nt5_md_ana5.out
mv 1nt5_md.pdb 1nt5_md5.pdb
mv 1nt5_md.crd 1nt5_md5.crd
