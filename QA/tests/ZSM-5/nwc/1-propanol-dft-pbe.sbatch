#!/bin/bash
#SBATCH --ntasks=18
#SBATCH --nodes=1
#SBATCH --output=%x.o%j
#SBATCH --time=00:30:00
#SBATCH --account=csigeneral
srun ../../../../bin/LINUX64/nwchem 1-propanol-dft-pbe.nw

