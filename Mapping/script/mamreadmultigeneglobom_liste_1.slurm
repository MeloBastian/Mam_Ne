#!/bin/bash
#SBATCH --time=80:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=20gb
#SBATCH -o readmultigeneglobom_subset_1.out
#SBATCH -e readmultigeneglobom_subset_1.err

mpirun -np 16 /beegfs/data/mbastian/bayescode/data/readmultigeneglobom -dsomss -x 75 5 -1 mam_multigeneglobom_liste_1
