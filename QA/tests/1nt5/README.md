# (1NT5)[https://www.rcsb.org/structure/1NT5]

This protein forms a transmembrane Sodium channel. It consists of two proteins that bind 
through Hydrogen bonds to make the channel. In this directory two sets of calculations
are given. One set running the dynamics in solution, an another running the dynamics in
vacuo.

Different simulations are set up by initially separating the two proteins by either
1, 2, 3, 4, or 5 nm. Subsequently the dynamics is run to see if the proteins 
recombine.

The calculations in solution are prepared by inputs labeled: `1nt5_prep?.nw`.
The calculations in vacuo are prepared by inputs labeled: `1nt5_prep1?.nw`.
The digit represented by "?" is replaced by the separation distance.

In solution only the structures with the proteins separated by 1 or 2 nm recombine.
In vacuo all structure recombine, and do so quickly.

## Results

Upon executing:

- run.sh
- run_vacuo.sh

the resulting trajectories and initial structure will be labeled:

- in solution:

  - 1nt5_md1.crd  1nt5_md1.pdb
  - 1nt5_md2.crd  1nt5_md2.pdb
  - 1nt5_md3.crd  1nt5_md3.pdb
  - 1nt5_md4.crd  1nt5_md4.pdb
  - 1nt5_md5.crd  1nt5_md5.pdb

- in vacuo:

  - 1nt5_md11.crd 1nt5_md11.pdb
  - 1nt5_md12.crd 1nt5_md12.pdb
  - 1nt5_md13.crd 1nt5_md13.pdb
  - 1nt5_md14.crd 1nt5_md14.pdb
  - 1nt5_md15.crd 1nt5_md15.pdb
