# 1H9T Protein bound to DNA

This structure is from the protein data bank, see [1H9T](https://www.rcsb.org/structure/1H9T). 
The system consists of a protein that binds to DNA as part of regulating the transcription of
genes. 

In this case we simulate the protein binding to the DNA. The process in this simulation happens
in vacuum, and we drive this process by giving the DNA an initial speed. None of this is
physically realistic but this way you get a result within a few thousand time steps.

The simulation executes by running the `run_restrain.sh` script. This script was written to
run on the Erato machine at BNL and relies on its queueing system. Hence there are a few
sub-scripts:
- `run_restrain.sh`
  - `1h9t_min_restrain.sh`
  - `1h9t_eq_restrain.sh`
  - `1h9t_md_restrain.sh`

In more detail the `run_restrain.sh` script does:
1. Run the prepare step
   - Inputs `start_md_neutral.pdb`
   - Outputs `1h9t-restrain.top` (the topology file)
   - Outputs `1h9t-restrain_md.rst` (the restart file)
2. Submit the minimization calculation
   - Inputs `1h9t-restrain.top`
   - Inputs `1h9t-restrain_md.rst`
   - Outputs `1h9t-restrain_md.qrs` (a new restart file)
3. Copy `1h9t-restrain_md.qrs` over `1h9t-restrain_md.rst` so that the minimization results 
   can be used downstream
4. Submit the equilibration calculation
   - Inputs `1h9t-restrain.top`
   - Inputs `1h9t-restrain_md.rst`
   - Outputs `1h9t-restrain_md.trj` (a trajectory file of what the equilibration does)
5. Analyze the equilibration results
   - Inputs `1h9t-restrain_md.rst`
   - Inputs `1h9t-restrain_md.trj`
   - Outputs `1h9t-restrain_md.pdb` (the molecular structure in PDB format)
   - Outputs `1h9t_eq_restrain.crd` (the trajectory in Amber with periodic coordinates format)
6. Submit the molecular dynamics simulation
   - Inputs `1h9t-restrain.top`
   - Inputs `1h9t-restrain_md.rst`
   - Outputs `1h9t-restrain_md.trj` (a trajectory file of what the molecular dynamics does)
7. Analyze the molecular dynamics results
   - Inputs `1h9t-restrain_md.rst`
   - Inputs `1h9t-restrain_md.trj`
   - Outputs `1h9t-restrain_md.pdb` (the molecular structure in PDB format)
   - Outputs `1h9t_md_restrain.crd` (the trajectory in Amber with periodic coordinates format)
8. Tar the results
   - Inputs `1h9t_md.pdb`
   - Inputs `1h9t_md_restrain.crd`
   - Outputs `1h9t_md_restrain.tgz`

The final tar-file contains everything you need to visualize the process using VMD.
