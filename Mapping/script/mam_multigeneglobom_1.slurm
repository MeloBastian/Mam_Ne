#!/bin/bash
#SBATCH --time=160:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=20gb
#SBATCH -o multigeneglobom_subset_1.out
#SBATCH -e multigeneglobom_subset_1.err

mpirun -np 16 /beegfs/data/mbastian/bayescode/data/multigeneglobom -t 1007forphylo.conc.treefile_rooted -d mam_subset_1_1000genes.ali -x 1 200 mam_multigeneglobom_liste_1
