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
