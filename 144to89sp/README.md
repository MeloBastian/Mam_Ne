# **/144to89species** (total size = 1.7G)

In this part of the analysis, the 144 species tree was reduced to remove closely related species. This was done using the cuttree program implemented in the bayescode suite : /bayescode/data/cuttree <chronogram> <cutoff> <out>.
We use 0.1 as a cutoff.
When 2 or more species are flagged as too close, one of them is choosen based on its coverage and on the avaibility of information about life history traits.
This new specie subset countain 89 species.
Stochastic mapping of substitution events is then realised using the same procedure as for the 144 species set (see /Mapping).
As there is less species, the subset of genes that can be simultaneously processed for the mapping run can be larger, so the genes are clustererd in sets of  4 * 1500 genes rather than 6 * 1000 genes.

"X" correspond to the subset id

- ### **/data**	

	- **/90splist_10ma_**:  89 species list
	
	- **/1500genes_subset4.conc.treefile_rooted**: rooted treefile from a Iqtree2 run using 1500 genes.
	
	- **/1500genes_subsetX.ali**: 4 multialignement of 1500 genes. 
	
	- **/1500genes_subsetX.conc**: 4 concatenation of 1500 genes.
		
	- **/merge_mam_multigeneglobom.genedsomsuffstat**: Sufficient statistics for the substitution mapping on the 4 sets of 1500 genes, merged in a unique file for fastcoevol.
