# **/Mapping** (total size = 2.8G)

After filtering the alignments, we end up with a set of 6000 genes.

[1] We aggregate these single-gene aligments in 6 groups of 1000 multi-alignement files (.ali files, a special format used by Bayescode, tabulating all alignments in a single file with heaser) and 6x1000 concatened files (.conc files, in phylip format).

[2] One of the 6 concatenation is used to compute a phylogeny using Iqtree2 and a GTR+4 model (iqtree2 -s {multigene.conc}  -m GTR+G4 -nt 8 -alrt 1000 -B 1000 --boot-trees).

[3] The multialignement and the tree are used for the mapping of substitutions.
   The mapping is realised using Bayescode Version: 0d767fc (https://github.com/bayesiancook/bayescode.git)
   6 independent mapping were realised corresponding to the 6 set of 1000 genes. The output are then merged in a unique file.

"X" refers to the genes subset id

- ### */script*
 
For the step 3, there is one script per gene list, they are all similar, only the script for the first list is provided.

	- **/ali_and_conc_afterphylterandbayescode.py** (input: {144_splist},{1000_genelist}; output: {multialignement_name.ali}, {multialignement_name.conc})
  	Write a multialignement and a concatenat of 1000 genes. (step 1)

	- **/mam_multigeneglobom_liste_1** (input: {tree}, {alignement}; output: {out_name})
   	Script for multigenemapping, addapted to run on a cluster (step 3)

	- **/mamreadmultigeneglobom_liste_1.slurm** (input: {out_mapping_name})
   	Treatment of the mappings (step 3).

- ### */data*

	- **/list6002genes**: final 6001 gene list

	- **6001genes.tar.gz**: the 6001 genes alignements

	- **1007forphylo.conc.treefile_rooted**: Phylogeny (output of step 2)

	- **/mam_merge_genedsomsuffstat**: Suffstat from the 6 * 1000 genes mapping, merged in a unique file for FastCoevol (output of step 3)
